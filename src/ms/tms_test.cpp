/**
 * Tests for the Orbit structures TmsIndex relies on (FL and phi built from
 * run heads and lengths, integrated data columns, serialization) and for
 * tms_query against naive matching statistics and ms_query.
 */

#include "tms_test.hpp"
#include "tms_index.hpp"
#include "tms_query.hpp"
#include "tms_smem.hpp"
#include "ms_rlbwt.hpp"
#include "tsv.hpp"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace tms_test {

TextBwt make_text_bwt(const std::string& s) {
    TextBwt t;
    t.str = s + "$";
    const size_t n = t.str.size();
    t.text.resize(n);
    for (size_t i = 0; i < n; ++i) {
        const char ch = t.str[i];
        t.text[i] = ch == '$' ? orbit::TERMINATOR : ch == '%' ? orbit::SEPARATOR : static_cast<uchar>(ch);
    }
    t.sa.resize(n);
    std::iota(t.sa.begin(), t.sa.end(), 0);
    const auto& x = t.text;
    std::sort(t.sa.begin(), t.sa.end(), [&](ulint a, ulint b) {
        while (x[a] == x[b]) { ++a; ++b; }  // the unique terminator stops this
        return x[a] < x[b];
    });
    t.isa.resize(n);
    for (size_t j = 0; j < n; ++j) t.isa[t.sa[j]] = j;
    t.lcp.assign(n, 0);
    for (size_t j = 1; j < n; ++j) {
        ulint a = t.sa[j - 1], b = t.sa[j], l = 0;
        while (x[a + l] == x[b + l]) ++l;
        t.lcp[j] = l;
    }
    for (size_t j = 0; j < n; ++j) {
        const uchar c = x[t.sa[j] == 0 ? n - 1 : t.sa[j] - 1];
        if (j == 0 || c != t.heads.back()) {
            t.heads.push_back(c);
            t.lens.push_back(0);
            t.lcps_per_run.emplace_back();
        }
        ++t.lens.back();
        t.lcps_per_run.back().push_back(t.lcp[j]);
    }
    return t;
}

std::vector<ulint> naive_ms(const std::string& T, const std::string& P) {
    std::vector<ulint> ms(P.size(), 0);
    for (size_t i = 0; i < P.size(); ++i) {
        size_t lo = 0;
        // MS[i] >= MS[i-1] - 1, so start the search there.
        if (i > 0 && ms[i - 1] > 1) lo = ms[i - 1] - 1;
        size_t L = lo;
        while (i + L < P.size() && T.find(P.substr(i, L + 1)) != std::string::npos) ++L;
        ms[i] = L;
    }
    return ms;
}

namespace {

/**
 * FL read forward from row 0 visits rows ISA[0], ISA[1], ... with F
 * characters T[0], T[1], ..., and LF read from row 0 visits rows ISA[n - 2],
 * ISA[n - 3], ... with BWT characters.  Together these say that FL inverts
 * LF on every row and that FL's character is the F character.
 */
void check_fl_and_lf(const TextBwt& t, const orbit::split_params& sp, const char* name) {
    const ulint n = t.text.size();
    orbit::rlbwt::fl_permutation<orbit::empty_data_columns, false, true> fl(t.heads, t.lens, sp);
    orbit::rlbwt::lf_move<true> lf(t.heads, t.lens, sp);
    auto q = fl.first();
    assert(q.idx == 0);
    for (ulint i = 0; i < n; ++i) {
        const ulint p = (n - 1 + i) % n;  // the suffix at row q
        assert(q.idx == t.isa[p] && "FL visits the next text position");
        assert(fl.get_character(q) == t.text[p] && "FL character is the F character");
        q = fl.FL(q);
    }
    assert(q.idx == 0);
    auto l = lf.first();
    for (ulint i = 0; i < n; ++i) {
        const ulint p = n - 1 - i;
        assert(l.idx == t.isa[p] && "LF visits the previous text position");
        assert(lf.get_character(l) == t.text[p == 0 ? n - 1 : p - 1]);
        l = lf.LF(l);
    }
    (void)name;
}

enum class PhiCols { PLCP, COUNT };

/**
 * phi from rlbwt_to_phi maps SA[j] to SA[j - 1] (SA[0] to SA[n - 1]), and a
 * PLCP column holding the value at each interval's start gives PLCP anywhere
 * by subtracting the offset.  The integrated column survives serialization.
 */
void check_phi(const TextBwt& t, const orbit::split_params& sp) {
    const ulint n = t.text.size();
    std::vector<ulint> plcp(n);
    for (ulint j = 0; j < n; ++j) plcp[t.sa[j]] = t.lcp[j];
    auto enc = orbit::rlbwt::rlbwt_to_phi(t.heads, t.lens, sp);
    std::vector<orbit::columns_tuple<PhiCols>> cols(enc.intervals());
    ulint start = 0;
    for (ulint i = 0; i < enc.intervals(); ++i) {
        cols[i][0] = plcp[start];
        start += enc.get_length(i);
    }
    assert(start == n);
    using Phi = orbit::rlbwt::phi_permutation<PhiCols, true>;
    Phi phi(enc, cols);
    std::stringstream ss;
    phi.serialize(ss);
    Phi loaded;
    loaded.load(ss);
    for (Phi* p : {&phi, &loaded}) {
        auto pos = p->first();
        for (ulint x = 0; x < n; ++x) {
            assert(pos.idx == x);
            const ulint j = t.isa[x];
            const ulint want = t.sa[j == 0 ? n - 1 : j - 1];
            assert(p->phi(pos).idx == want && "phi(SA[j]) = SA[j - 1]");
            assert(p->template get<PhiCols::PLCP>(pos) - pos.offset == plcp[x] && "PLCP from the sample");
            if (pos.offset + 1 < p->get_length(pos.interval)) {
                ++pos.offset;
                ++pos.idx;
            } else if (x + 1 < n) {
                pos = p->down(pos);
            }
        }
    }
}

/** Texts from the families that found bugs in TeraMS, plus haplotype sets. */
std::vector<std::string> fuzz_texts(std::mt19937& rng) {
    std::vector<std::string> texts;
    auto rnd = [&](size_t len, const char* alpha, size_t k) {
        std::string s(len, 'A');
        for (auto& c : s) c = alpha[rng() % k];
        return s;
    };
    auto mutate = [&](std::string s, double rate) {
        std::uniform_real_distribution<double> u(0, 1);
        for (auto& c : s)
            if (u(rng) < rate) c = "ACGT"[rng() % 4];
        return s;
    };
    for (size_t len : {1, 2, 3, 4, 5}) texts.push_back(rnd(len, "ACGT", 4));
    texts.push_back("A");
    texts.push_back("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    texts.push_back(std::string(500, 'C'));
    texts.push_back(std::string(200, 'A') + std::string(200, 'C') + std::string(200, 'A'));
    for (int k = 0; k < 6; ++k) {
        std::string unit = rnd(1 + rng() % 7, "ACGT", 4), s;
        while (s.size() < 1200) s += unit;
        texts.push_back(k % 2 ? mutate(s, 0.002) : s);
    }
    for (int k = 0; k < 4; ++k) {
        std::string s = rnd(800, "ACGT", 4);
        s.replace(100 + rng() % 200, 50, std::string(50, 'N'));
        s.replace(500 + rng() % 200, 5, std::string(5, 'N'));
        texts.push_back(s);
    }
    for (int k = 0; k < 6; ++k) {
        const std::string base = rnd(300 + rng() % 700, "ACGT", 4);
        const int h = 2 + rng() % 6;
        std::string s = base;
        for (int j = 1; j < h; ++j) s += "%" + (k % 3 == 0 ? base : mutate(base, 0.01));
        texts.push_back(s);
    }
    for (int k = 0; k < 4; ++k) texts.push_back(rnd(2000, "ACGT", 4));
    return texts;
}

/** Patterns drawn from T with mutations, some with N or absent bytes. */
std::vector<std::string> fuzz_patterns(const std::string& T, std::mt19937& rng, int count) {
    std::vector<std::string> pats;
    std::uniform_real_distribution<double> u(0, 1);
    const double rates[] = {0.0, 0.01, 0.05, 0.3, 1.0};
    for (int k = 0; k < count; ++k) {
        const size_t len = (k % 40 == 0) ? 0 : (k % 40 == 1) ? 1 : rng() % 150;
        const size_t tlen = T.size() - 1;  // without the terminator
        std::string P;
        if (len <= tlen) P = T.substr(rng() % (tlen - len + 1), len);
        else P = T.substr(0, tlen);
        for (auto& c : P)
            if (c == '%' || u(rng) < rates[k % 5]) c = "ACGT"[rng() % 4];
        if (k % 7 == 0 && !P.empty()) P[rng() % P.size()] = 'N';
        if (k % 11 == 0 && !P.empty()) P[rng() % P.size()] = 'X';
        pats.push_back(std::move(P));
    }
    return pats;
}

/** An ms index over the same RLBWT, from full per-run LCPs. */
MSIndexSpillLCP<false> ms_index_for(const TextBwt& t) {
    auto [run_data, spill, max_top, max_sub, skinny, jumbo] = build_spill_data(t.lcps_per_run);
    return MSIndexSpillLCP<false>(t.heads, t.lens, run_data, std::move(spill), max_top, max_sub);
}

std::vector<TmsBuildOptions> split_variants() {
    TmsBuildOptions a;
    TmsBuildOptions b;
    b.fl_split = orbit::NO_SPLITTING;
    TmsBuildOptions c;
    c.lf_split = orbit::split_params{};
    TmsBuildOptions d;
    d.lf_split = orbit::split_params(std::nullopt, 2);
    d.fl_split = orbit::split_params(std::nullopt, 2);
    return {a, b, c, d};
}

/** An RLBWT with each run's top LCP and the text it came from. */
struct Input { std::vector<uchar> heads; std::vector<ulint> lens, tops; std::string text; };

/** The fuzz texts and, if its LF is a single cycle, minishred. */
std::vector<Input> batch_inputs(const std::string& data_dir, std::mt19937& rng) {
    std::vector<Input> inputs;
    for (const auto& s : fuzz_texts(rng)) {
        TextBwt t = make_text_bwt(s);
        Input in{t.heads, t.lens, {}, t.str};
        for (const auto& l : t.lcps_per_run) in.tops.push_back(l[0]);
        inputs.push_back(std::move(in));
    }
    {
        std::vector<uchar> heads;
        std::vector<ulint> lens;
        std::vector<std::vector<ulint>> lcps;
        if (tsv::load_tsv(data_dir + "/minishred1_20_002_lcp.tsv", heads, lens, lcps)) {
            Input in{heads, lens, {}, {}};
            for (const auto& l : lcps) in.tops.push_back(l[0]);
            try {
                TmsIndex probe(heads, lens, TmsBuildOptions{}, &in.tops);
                auto pos = probe.first();
                // Row 0 holds text position n - 1, so the s-th row LF visits
                // from it has BWT character T[n - 2 - s].
                const ulint n = probe.domain();
                std::string T(n, ' ');
                for (ulint st = 0; st < n; ++st) {
                    const uchar c = probe.get_character(pos.interval);
                    T[(2 * n - 2 - st) % n] = c == orbit::TERMINATOR ? '$' : c == orbit::SEPARATOR ? '%' : static_cast<char>(c);
                    pos = probe.LF_step(pos);
                }
                in.text = T;
                inputs.push_back(std::move(in));
            } catch (const std::exception& e) {
                std::cout << "  minishred skipped: " << e.what() << std::endl;
            }
        }
    }
    return inputs;
}

}  // namespace

bool test_orbit_structures() {
    std::cout << "Testing Orbit FL and phi built from run heads and lengths" << std::endl;
    std::mt19937 rng(5);
    size_t count = 0;
    for (const auto& s : fuzz_texts(rng)) {
        TextBwt t = make_text_bwt(s);
        for (const auto& sp : {orbit::NO_SPLITTING, orbit::split_params{}, orbit::split_params(std::nullopt, 2)}) {
            check_fl_and_lf(t, sp, "");
            check_phi(t, sp);
        }
        ++count;
    }
    std::cout << "  " << count << " texts PASSED" << std::endl;
    return true;
}

bool test_tms_vs_naive() {
    std::cout << "Testing tms_query against naive matching statistics and ms_query" << std::endl;
    std::mt19937 rng(17);
    size_t checked = 0;
    for (const auto& s : fuzz_texts(rng)) {
        TextBwt t = make_text_bwt(s);
        auto ms_idx = ms_index_for(t);
        std::vector<TmsIndex> idxs;
        for (const auto& o : split_variants()) idxs.emplace_back(t.heads, t.lens, o);
        for (const auto& P : fuzz_patterns(t.str, rng, 60)) {
            const auto want = naive_ms(t.str, P);
            if (ms_query(ms_idx, P) != want) {
                std::cout << "  ms_query disagrees with naive (not a tms failure) on text of length "
                          << s.size() << std::endl;
            }
            for (auto& idx : idxs) {
                if (tms_query(idx, P) != want) {
                    std::cout << "  FAILED: text=" << (s.size() < 80 ? s : s.substr(0, 80) + "...")
                              << " pattern=" << P << std::endl;
                    assert(false && "tms_query must match naive MS");
                    return false;
                }
                checked += P.size();
            }
        }
    }
    std::cout << "  " << checked << " values PASSED" << std::endl;
    return true;
}

bool test_tms_vs_ms_minishred(const std::string& data_dir) {
    std::cout << "Testing tms_query against ms_query on minishred" << std::endl;
    std::vector<uchar> heads;
    std::vector<ulint> lens;
    std::vector<std::vector<ulint>> lcps;
    if (!tsv::load_tsv(data_dir + "/minishred1_20_002_lcp.tsv", heads, lens, lcps)) {
        std::cout << "  DID NOT RUN" << std::endl;
        return false;
    }
    auto [run_data, spill, max_top, max_sub, skinny, jumbo] = build_spill_data(lcps);
    MSIndexSpillLCP<false> ms_idx(heads, lens, run_data, std::move(spill), max_top, max_sub);
    // The text, read back with LF from the first row.
    std::string T;
    {
        auto pos = ms_idx.first();
        for (ulint i = 0; i < ms_idx.domain(); ++i) {
            T += static_cast<char>(ms_idx.get_character(pos.interval));
            pos = ms_idx.LF(pos);
        }
        std::reverse(T.begin(), T.end());
        for (auto& c : T) if (static_cast<uchar>(c) == orbit::SEPARATOR) c = '%';
    }
    std::mt19937 rng(23);
    auto pats = fuzz_patterns(T, rng, 400);
    size_t checked = 0;
    for (const auto& o : split_variants()) {
        TmsIndex idx(heads, lens, o);
        std::stringstream ss;
        idx.serialize(ss);
        TmsIndex loaded;
        loaded.load(ss);
        for (const auto& P : pats) {
            const auto want = ms_query(ms_idx, P);
            assert(tms_query(idx, P) == want && "tms_query must match ms_query");
            assert(tms_query(loaded, P) == want && "a loaded index must answer identically");
            checked += P.size();
        }
    }
    std::cout << "  " << checked << " values PASSED" << std::endl;
    return true;
}

/**
 * TmsIndex::PackedAccess fits the index and agrees with the index's own
 * column reads for every LF row (length, character, the LF step from every
 * offset, psi_at_tail, resolving a position past a row's end) and for every
 * FL row (psi_step from every offset, and from one past the row's end).
 * Returns the number of rows checked.
 */
size_t check_packed_access(TmsIndex& idx) {
    const TmsIndex::PackedAccess pa = idx.packed_access();
    assert(pa.fits() && "the test indexes must fit packed access");
    std::array<uchar, 256> lf_code, fl_code;
    for (size_t c = 0; c < 256; ++c) {
        lf_code[c] = idx.lf_code(static_cast<uchar>(c));
        fl_code[c] = idx.fl_code(static_cast<uchar>(c));
        assert((lf_code[c] != TmsIndex::no_code) == idx.occurs(static_cast<uchar>(c)));
    }
    size_t rows = 0;
    for (ulint i = 0; i < idx.move_runs(); ++i, ++rows) {
        const auto r = pa.row(i);
        assert(pa.length(r) == idx.get_length(i));
        assert(pa.code(r) == lf_code[idx.get_character(i)]);
        const TmsIndex::FLPos q = pa.psi_at_tail(r), want_q = idx.psi_at_tail(i);
        assert(q.interval == want_q.interval && q.offset == want_q.offset);
        for (ulint o = 0; o < idx.get_length(i); ++o) {
            TmsIndex::LFPos p{};
            p.interval = i;
            p.offset = o;
            const TmsIndex::LFPos got = pa.start_LF(r, p), want = idx.start_LF(p);
            assert(got.interval == want.interval && got.offset == want.offset);
        }
        if (i + 1 < idx.move_runs()) {
            TmsIndex::LFPos p{};
            p.interval = i;
            p.offset = idx.get_length(i);
            const TmsIndex::LFPos want = idx.finish_LF(p);
            const auto rr = pa.resolve_LF(p);
            assert(p.interval == want.interval && p.offset == want.offset && rr == pa.row(want.interval));
        }
    }
    for (ulint i = 0; i < idx.psi_intervals(); ++i, ++rows) {
        const ulint len = idx.fl().get_length(i);
        for (ulint o = 0; o <= len; ++o) {
            if (o == len && i + 1 == idx.psi_intervals()) break;
            TmsIndex::FLPos q{};
            q.interval = i;
            q.offset = o;
            uchar c = 0, want_c = 0;
            const TmsIndex::FLPos got = pa.psi_step(q, c), want = idx.psi_step(q, want_c);
            assert(got.interval == want.interval && got.offset == want.offset && c == fl_code[want_c]);
        }
    }
    return rows;
}

/**
 * tms_query_batch in every mode, with and without positions, gives the
 * lengths tms_query gives, for several numbers of patterns in flight.  The
 * positions are the same in every mode and for every k, and each one is an
 * occurrence: T[pos .. pos + len) = P[i .. i + len), with TMS_NO_POS exactly
 * where the length is 0.  Runs on the fuzz texts and, if its LF is a single
 * cycle, on minishred.
 */
bool test_tms_batch(const std::string& data_dir) {
    std::cout << "Testing tms_query_batch in every mode against tms_query" << std::endl;
    std::mt19937 rng(41);
    const auto inputs = batch_inputs(data_dir, rng);
    const TmsMode modes[] = {TmsMode::PSI, TmsMode::PHI, TmsMode::PHISKIP, TmsMode::DUAL};
    size_t checked = 0, positions = 0, access_rows = 0;
    for (const auto& in : inputs) {
        const std::string& T = in.text;
        auto pats = fuzz_patterns(T, rng, 200);
        auto variants = split_variants();
        TmsBuildOptions nophi;
        nophi.phi_split = orbit::NO_SPLITTING;
        variants.push_back(nophi);
        for (const auto& o : variants) {
            TmsIndex idx(in.heads, in.lens, o, &in.tops);
            std::stringstream ss;
            idx.serialize(ss);
            TmsIndex loaded;
            loaded.load(ss);
            access_rows += check_packed_access(idx);
            std::vector<std::vector<ulint>> want;
            for (const auto& P : pats) want.push_back(tms_query(idx, P));
            std::vector<std::vector<ulint>> first_pos;
            for (TmsMode mode : modes) {
                for (size_t k : {1, 3, 32, 1000}) {
                    for (int access = 0; access < 4; ++access) {
                        // Packed and column access, without and with positions.
                        const bool packed = access < 2, with_pos = access % 2 == 1;
                        std::vector<std::vector<ulint>> got, pos;
                        tms_query_batch(k == 3 ? loaded : idx, pats, k, got, mode, with_pos ? &pos : nullptr, packed);
                        if (got != want) {
                            std::cout << "  FAILED lengths: mode " << int(mode) << ", k " << k << ", packed " << packed
                                      << std::endl;
                            assert(false && "tms_query_batch lengths must match tms_query");
                            return false;
                        }
                        ++checked;
                        if (!with_pos) continue;
                        if (first_pos.empty()) {
                            first_pos = pos;
                            for (size_t j = 0; j < pats.size(); ++j)
                                for (size_t i = 0; i < pats[j].size(); ++i) {
                                    const ulint l = got[j][i], p = pos[j][i];
                                    const bool ok = (l == 0) ? (p == TMS_NO_POS)
                                                             : (p + l <= T.size() && T.compare(p, l, pats[j], i, l) == 0);
                                    if (!ok) {
                                        std::cout << "  FAILED position: pattern " << pats[j] << " at " << i << " len " << l
                                                  << " pos " << p << std::endl;
                                        assert(false && "each position must be an occurrence");
                                        return false;
                                    }
                                    ++positions;
                                }
                        } else if (pos != first_pos) {
                            std::cout << "  FAILED: positions differ in mode " << int(mode) << ", k " << k << std::endl;
                            assert(false && "positions must not depend on mode or k");
                            return false;
                        }
                    }
                }
            }
        }
    }
    std::cout << "  " << access_rows << " rows of packed access, " << checked << " batches and " << positions
              << " positions PASSED" << std::endl;
    return true;
}

namespace {

// The most intervals a walker step fast-forwarded over in check_walker.
ulint g_max_fast_forward = 0;

/**
 * A TmsIndex walker over a permutation with the given interval starts, next
 * and lcp: locate finds the same interval as a plain binary search over the
 * starts for every text position, and from there step gives the text
 * position and LCP of next, twice in a row.
 */
template <typename Walker, typename StartOf, typename Next, typename Lcp>
void check_walker(const Walker& wk, ulint n, StartOf start_of, ulint intervals, Next next, Lcp lcp) {
    auto plain = [&](ulint x) {
        ulint lo = 0, hi = intervals;
        while (hi - lo > 1) {
            const ulint mid = lo + (hi - lo) / 2;
            (start_of(mid) <= x ? lo : hi) = mid;
        }
        return lo;
    };
    for (ulint x = 0; x < n; ++x) {
        const auto p = wk.locate(x);
        assert(p.interval == plain(x) && p.offset == x - start_of(p.interval) && p.idx == x);
        assert(wk.lcp(p) == lcp(p));
        auto q = wk.start_step(p);
        auto a = p;
        for (int s = 0; s < 2; ++s) {
            const ulint from = q.interval;
            ulint l;
            const ulint y = wk.step(q, l);
            a = next(a);
            assert(y == a.idx && l == lcp(a) && "a walker step agrees with phi or phi_inv");
            g_max_fast_forward = std::max(g_max_fast_forward, a.interval - from);
        }
    }
}

/**
 * TmsIndex's phi_inv maps SA[j] to SA[j + 1] (SA[n - 1] to SA[0]) and its
 * PLCPB is the LCP with the row below (0 for row n - 1); phi_at and
 * phi_inv_at locate every text position; both survive serialization.
 * Unsplit, phi_inv has the intervals Orbit's rlbwt_to_phi_inv gives.
 */
void check_phi_inv(const TextBwt& t, const TmsBuildOptions& o) {
    const ulint n = t.text.size();
    std::vector<ulint> tops;
    for (const auto& l : t.lcps_per_run) tops.push_back(l[0]);
    TmsIndex idx(t.heads, t.lens, o, &tops);
    assert(idx.has_phi_inv());
    std::stringstream ss;
    idx.serialize(ss);
    TmsIndex loaded;
    loaded.load(ss);
    assert(loaded.has_phi_inv());
    if (o.phi_split == orbit::NO_SPLITTING)
        assert(idx.phi_inv_intervals() == orbit::rlbwt::rlbwt_to_phi_inv(t.heads, t.lens, orbit::NO_SPLITTING).intervals());
    for (TmsIndex* x : {&idx, &loaded}) {
        for (ulint p = 0; p < n; ++p) {
            const ulint j = t.isa[p];
            auto a = x->phi_at(p);
            assert(a.idx == p);
            assert(x->plcp(a) == t.lcp[j] && "PLCP at a located phi point");
            assert(x->phi(a).idx == t.sa[j == 0 ? n - 1 : j - 1]);
            auto b = x->phi_inv_at(p);
            assert(b.idx == p);
            assert(x->plcpb(b) == (j + 1 < n ? t.lcp[j + 1] : 0) && "PLCPB is the LCP with the row below");
            assert(x->phi_inv(b).idx == t.sa[j + 1 < n ? j + 1 : 0] && "phi_inv(SA[j]) = SA[j + 1]");
        }
        // The walkers, with start tables of several densities and none.
        for (int shift : {-1, 0, 1, 3, 20}) {
            x->build_start_tables(shift);
            check_walker(x->phi_walker(), n, [&](ulint i) { return x->phi_start(i); }, x->phi_intervals(),
                         [&](TmsIndex::PhiPos p) { return x->phi(p); }, [&](TmsIndex::PhiPos p) { return x->plcp(p); });
            check_walker(x->phi_inv_walker(), n, [&](ulint i) { return x->phi_inv_start(i); }, x->phi_inv_intervals(),
                         [&](TmsIndex::PhiInvPos p) { return x->phi_inv(p); },
                         [&](TmsIndex::PhiInvPos p) { return x->plcpb(p); });
        }
    }
    // Loaded without phi_inv, phi has no start table.
    std::stringstream again(ss.str());
    TmsIndex phi_only;
    phi_only.load(again, TmsParts{true, true, false});
    assert(!phi_only.has_phi_inv() && phi_only.start_table_bytes() == 0);
    check_walker(phi_only.phi_walker(), n, [&](ulint i) { return phi_only.phi_start(i); }, phi_only.phi_intervals(),
                 [&](TmsIndex::PhiPos p) { return phi_only.phi(p); }, [&](TmsIndex::PhiPos p) { return phi_only.plcp(p); });
}

/**
 * SMEMs of P against T by brute force: every substring that occurs and
 * extends in neither direction (a maximal exact match), minus those inside
 * another.  Returns (start, length) pairs sorted by start.
 */
std::vector<std::pair<size_t, size_t>> naive_smems(const std::string& T, const std::string& P) {
    const size_t m = P.size();
    auto occurs = [&](size_t a, size_t b) { return T.find(P.substr(a, b - a)) != std::string::npos; };
    std::vector<std::pair<size_t, size_t>> mems;
    for (size_t a = 0; a < m; ++a)
        for (size_t b = a + 1; b <= m; ++b) {
            if (!occurs(a, b)) break;
            if ((b == m || !occurs(a, b + 1)) && (a == 0 || !occurs(a - 1, b))) mems.push_back({a, b - a});
        }
    std::vector<std::pair<size_t, size_t>> out;
    for (const auto& x : mems) {
        bool inside = false;
        for (const auto& y : mems)
            if (y != x && y.first <= x.first && x.first + x.second <= y.first + y.second) inside = true;
        if (!inside) out.push_back(x);
    }
    return out;
}

}  // namespace

bool test_phi_inv() {
    std::cout << "Testing phi_inv and PLCPB in TmsIndex" << std::endl;
    std::mt19937 rng(7);
    size_t count = 0;
    for (const auto& s : fuzz_texts(rng)) {
        TextBwt t = make_text_bwt(s);
        for (const auto& sp : {orbit::NO_SPLITTING, orbit::split_params{}, orbit::split_params(std::nullopt, 2)}) {
            TmsBuildOptions o;
            o.phi_split = sp;
            check_phi_inv(t, o);
        }
        ++count;
    }
    // Some walker steps fast-forward over several intervals.
    assert(g_max_fast_forward >= 4);
    std::cout << "  " << count << " texts PASSED (walker steps fast-forward up to " << g_max_fast_forward
              << " intervals)" << std::endl;
    return true;
}

bool test_smem_detection() {
    std::cout << "Testing SMEM detection from matching statistics against brute force" << std::endl;
    std::mt19937 rng(29);
    size_t checked = 0;
    for (const auto& s : fuzz_texts(rng)) {
        const std::string T = s + "$";
        auto pats = fuzz_patterns(T, rng, 40);
        for (auto& P : pats) {
            if (P.size() > 60) P.resize(60);
            const auto ms = naive_ms(T, P);
            const auto want = naive_smems(T, P);
            for (ulint min_len : {0, 1, 3, 10}) {
                std::vector<std::pair<size_t, size_t>> got, want_min;
                for (size_t i : tms_smems(ms, min_len)) got.push_back({i, ms[i]});
                for (const auto& x : want)
                    if (x.second > min_len) want_min.push_back(x);
                if (got != want_min) {
                    std::cout << "  FAILED: pattern " << P << " min_len " << min_len << std::endl;
                    assert(false && "tms_smems must match brute-force SMEMs");
                    return false;
                }
                checked += got.size();
            }
        }
    }
    std::cout << "  " << checked << " SMEMs PASSED" << std::endl;
    return true;
}

/**
 * tms_report_smems on the fuzz texts and minishred: for SMEM_ALL, each
 * SMEM's positions are exactly the start positions of P[i .. i + L) in T,
 * without repeats, with count c; on the fuzz texts they are also in
 * increasing row order and consecutive.  For SMEM_ONE, the position is one
 * of them.  Both are the same in every mode and for every k, and the SMEMs
 * are those tms_smems gives.  SMEM_ALL is also the same for every number of
 * walks in flight in tms_report_smems_batch (0, which lists one SMEM at a
 * time; 1; a few; and more than there are walks), and with max_listed N it
 * lists the first N positions of the full list and keeps the count.
 */
bool test_smem_report(const std::string& data_dir) {
    std::cout << "Testing SMEM reports against every occurrence" << std::endl;
    std::mt19937 rng(43);
    const auto inputs = batch_inputs(data_dir, rng);
    const TmsMode modes[] = {TmsMode::PSI, TmsMode::PHI, TmsMode::PHISKIP, TmsMode::DUAL};
    size_t smems = 0, positions = 0, capped = 0;
    // SMEMs with no rows above the toehold's, none below, and some on both sides.
    size_t none_up = 0, none_down = 0, both = 0;
    for (const auto& in : inputs) {
        const std::string& T = in.text;
        // Row order is checked where the suffix array is at hand.
        std::vector<ulint> isa;
        if (T.size() < 5000) {
            TextBwt t = make_text_bwt(T.substr(0, T.size() - 1));
            isa = t.isa;
        }
        auto pats = fuzz_patterns(T, rng, in.text.size() > 5000 ? 60 : 100);
        for (const auto& o : split_variants()) {
            TmsIndex idx(in.heads, in.lens, o, &in.tops);
            for (ulint min_len : {0, 5}) {
                std::vector<std::vector<TmsSmemHits>> first[2];
                for (TmsMode mode : modes) {
                    for (size_t k : {1, 32}) {
                        std::vector<std::vector<ulint>> len, pos;
                        tms_query_batch(idx, pats, k, len, mode, &pos);
                        for (int r = 0; r < 2; ++r) {
                            const TmsReport report = r ? TmsReport::SMEM_ALL : TmsReport::SMEM_ONE;
                            std::vector<TmsSmemHits> hits;
                            tms_report_smems_batch(idx, len, pos, min_len, report, 0, hits);
                            if (report == TmsReport::SMEM_ALL) {
                                for (size_t sk : {size_t(1), size_t(2), size_t(7), size_t(32), size_t(1) << 20}) {
                                    std::vector<TmsSmemHits> other;
                                    tms_report_smems_batch(idx, len, pos, min_len, report, sk, other);
                                    for (size_t j = 0; j < pats.size(); ++j) {
                                        std::string a, b;
                                        tms_format_smems(hits[j], report, a);
                                        tms_format_smems(other[j], report, b);
                                        if (a != b) {
                                            std::cout << "  FAILED: batched SMEM listing differs with " << sk << " walks" << std::endl;
                                            assert(false && "SMEM positions must not depend on the walks in flight");
                                            return false;
                                        }
                                    }
                                }
                                for (ulint cap : {0, 1, 2, 5})
                                    for (size_t sk : {size_t(0), size_t(3), size_t(32)}) {
                                        std::vector<TmsSmemHits> other;
                                        tms_report_smems_batch(idx, len, pos, min_len, report, sk, other, cap);
                                        for (size_t j = 0; j < pats.size(); ++j) {
                                            const auto& full = hits[j];
                                            const auto& part = other[j];
                                            assert(part.smems.size() == full.smems.size());
                                            size_t at_full = 0, at_part = 0;
                                            for (size_t x = 0; x < full.smems.size(); ++x) {
                                                const auto& f = full.smems[x];
                                                const auto& g = part.smems[x];
                                                const ulint want = std::min(f.count, cap);
                                                bool ok = g.start == f.start && g.len == f.len && g.count == f.count &&
                                                          g.listed == want;
                                                for (ulint y = 0; ok && y < want; ++y)
                                                    ok = part.pos[at_part + y] == full.pos[at_full + y];
                                                if (!ok) {
                                                    std::cout << "  FAILED: max_listed " << cap << " with " << sk << " walks" << std::endl;
                                                    assert(false && "a capped listing must be the first positions in row order");
                                                    return false;
                                                }
                                                capped += f.count > cap;
                                                at_full += f.listed;
                                                at_part += g.listed;
                                            }
                                            assert(at_part == part.pos.size());
                                        }
                                    }
                            }
                            if (!first[r].empty()) {
                                std::vector<std::string> a, b;
                                for (size_t j = 0; j < pats.size(); ++j) {
                                    a.emplace_back();
                                    b.emplace_back();
                                    tms_format_smems(hits[j], report, a.back());
                                    tms_format_smems(first[r][0][j], report, b.back());
                                }
                                if (a != b) {
                                    std::cout << "  FAILED: SMEMs differ in mode " << int(mode) << ", k " << k << std::endl;
                                    assert(false && "SMEM reports must not depend on mode or k");
                                    return false;
                                }
                                continue;
                            }
                            first[r].push_back(hits);
                            for (size_t j = 0; j < pats.size(); ++j) {
                                const auto starts = tms_smems(len[j], min_len);
                                assert(starts.size() == hits[j].smems.size());
                                size_t at = 0;
                                for (size_t x = 0; x < starts.size(); ++x) {
                                    const auto& h = hits[j].smems[x];
                                    assert(h.start == starts[x] && h.len == len[j][h.start]);
                                    const std::string sub = pats[j].substr(h.start, h.len);
                                    std::vector<ulint> all;
                                    for (size_t f = T.find(sub); f != std::string::npos; f = T.find(sub, f + 1)) all.push_back(f);
                                    assert(h.listed == h.count);
                                    std::vector<ulint> got(hits[j].pos.begin() + at, hits[j].pos.begin() + at + h.count);
                                    at += h.count;
                                    if (report == TmsReport::SMEM_ALL) {
                                        const bool up = got.front() != pos[j][h.start];
                                        const bool down = got.back() != pos[j][h.start];
                                        none_up += !up;
                                        none_down += !down;
                                        both += up && down;
                                    }
                                    bool ok;
                                    if (report == TmsReport::SMEM_ONE) {
                                        ok = h.count == 1 && std::binary_search(all.begin(), all.end(), got[0]);
                                    } else {
                                        if (!isa.empty())
                                            for (size_t y = 1; y < got.size(); ++y)
                                                if (isa[got[y]] != isa[got[y - 1]] + 1) {
                                                    std::cout << "  FAILED: positions not in row order" << std::endl;
                                                    assert(false && "positions must be in row order");
                                                    return false;
                                                }
                                        std::vector<ulint> sorted = got;
                                        std::sort(sorted.begin(), sorted.end());
                                        ok = h.count == got.size() && sorted == all;
                                    }
                                    if (!ok) {
                                        std::cout << "  FAILED: pattern " << pats[j] << " SMEM " << h.start << ":" << h.len
                                                  << " has " << h.count << " positions, text has " << all.size() << std::endl;
                                        assert(false && "SMEM positions must be the occurrences");
                                        return false;
                                    }
                                    ++smems;
                                    positions += got.size();
                                }
                                assert(at == hits[j].pos.size());
                            }
                        }
                    }
                }
            }
        }
    }
    assert(none_up > 0 && none_down > 0 && both > 0 && capped > 0 && "the inputs must exercise every kind of walk");
    std::cout << "  " << smems << " SMEMs and " << positions << " positions (" << none_up << " with none above, "
              << none_down << " with none below, " << both << " with some on both sides) PASSED" << std::endl;
    return true;
}

namespace {

/** Serialize idx and load it back, with every part or only those in need. */
TmsIndex round_trip(TmsIndex& idx, const TmsParts* need = nullptr) {
    std::stringstream ss;
    idx.serialize(ss);
    TmsIndex loaded;
    if (need) loaded.load(ss, *need);
    else loaded.load(ss);
    return loaded;
}

/** One query's output, as tms-batch would print it: lengths, positions and SMEMs per pattern. */
std::vector<std::string> query_output(TmsIndex& idx, const std::vector<std::string>& pats, size_t k, TmsMode mode,
                                      bool positions, TmsReport report, bool packed) {
    std::vector<std::vector<ulint>> len, pos;
    const bool toehold = positions || report != TmsReport::MS;
    tms_query_batch(idx, pats, k, len, mode, toehold ? &pos : nullptr, packed);
    std::vector<TmsSmemHits> hits;
    if (report != TmsReport::MS) tms_report_smems_batch(idx, len, pos, 0, report, k, hits);
    std::vector<std::string> out(pats.size());
    for (size_t j = 0; j < pats.size(); ++j) {
        for (ulint x : len[j]) out[j] += std::to_string(x) + ' ';
        if (positions)
            for (ulint x : pos[j]) out[j] += (x == TMS_NO_POS ? std::string("-1") : std::to_string(x)) + ' ';
        if (report != TmsReport::MS) tms_format_smems(hits[j], report, out[j]);
    }
    return out;
}

/** Whether f() throws an E whose message contains what. */
template <typename E, typename F>
bool throws(F f, const char* what) {
    try {
        f();
    } catch (const E& e) {
        return std::string(e.what()).find(what) != std::string::npos;
    }
    return false;
}

}  // namespace

/**
 * Indexes built with the psi and phi layouts, and full indexes loaded with
 * only the parts a query needs, give the same output as the full index for
 * every query they can serve, with packed and column access and several
 * numbers of patterns in flight.  The LF columns a layout omits have width
 * 0, and the parts it omits are absent.  Queries that need a missing part are
 * refused, at load and at query time, and so are layouts given the wrong LCP
 * input.
 */
bool test_layouts(const std::string& data_dir) {
    std::cout << "Testing psi and phi layouts and partial loads against the full layout" << std::endl;
    std::mt19937 rng(47);
    const auto inputs = batch_inputs(data_dir, rng);
    size_t compared = 0;
    for (const auto& in : inputs) {
        auto pats = fuzz_patterns(in.text, rng, in.text.size() > 5000 ? 60 : 100);
        for (const auto& base : {split_variants()[0], split_variants()[3]}) {
            TmsBuildOptions psi_opts = base, phi_opts = base, phi_noinv_opts = base;
            psi_opts.layout = TmsLayout::PSI;
            phi_opts.layout = phi_noinv_opts.layout = TmsLayout::PHI;
            phi_noinv_opts.phi_inv = false;
            TmsIndex full(in.heads, in.lens, base, &in.tops);
            TmsIndex psi_built(in.heads, in.lens, psi_opts);
            TmsIndex phi_built(in.heads, in.lens, phi_opts, &in.tops);
            TmsIndex phi_noinv_built(in.heads, in.lens, phi_noinv_opts, &in.tops);
            TmsIndex psi = round_trip(psi_built), phi = round_trip(phi_built), phi_noinv = round_trip(phi_noinv_built);
            // Parts held, and LF columns of width 0 exactly where a part is missing.
            auto parts_are = [](const TmsIndex& x, bool a, bool b, bool c) {
                return x.has_psi() == a && x.has_phi() == b && x.has_phi_inv() == c;
            };
            assert(parts_are(full, true, true, true) && parts_are(psi, true, false, false));
            assert(parts_are(phi, false, true, true) && parts_are(phi_noinv, false, true, false));
            const auto& fw = full.lf().get_widths();
            const auto& sw = psi.lf().get_widths();
            const auto& hw = phi.lf().get_widths();
            constexpr size_t PSI_INT = TmsIndex::LF::template data_column<TmsLFCols::PSI_INT>();
            for (size_t c = 0; c < fw.size(); ++c) {
                const bool psi_col = c == PSI_INT || c == PSI_INT + 1, phi_col = c == PSI_INT + 2 || c == PSI_INT + 3;
                assert(fw[c] > 0);
                assert(sw[c] == (phi_col ? 0 : fw[c]) && hw[c] == (psi_col ? 0 : fw[c]));
            }
            assert(psi.packed_access().fits() && phi.packed_access().fits());

            // Full indexes loaded with only what psi MS, and what phi queries, read.
            const TmsParts psi_need = tms_query_parts(TmsMode::PSI, false, TmsReport::MS);
            const TmsParts phi_need = tms_query_parts(TmsMode::PHI, true, TmsReport::SMEM_ALL);
            TmsIndex full_psi = round_trip(full, &psi_need), full_phi = round_trip(full, &phi_need);
            assert(parts_are(full_psi, true, false, false) && parts_are(full_phi, false, true, true));

            for (size_t k : {1, 7, 32}) {
                for (bool packed : {true, false}) {
                    const auto want = query_output(full, pats, k, TmsMode::PSI, false, TmsReport::MS, packed);
                    for (TmsIndex* x : {&psi, &full_psi, &psi_built}) {
                        assert(query_output(*x, pats, k, TmsMode::PSI, false, TmsReport::MS, packed) == want);
                        ++compared;
                    }
                    for (bool positions : {false, true})
                        for (TmsReport report : {TmsReport::MS, TmsReport::SMEM_ONE, TmsReport::SMEM_ALL}) {
                            const auto want_phi = query_output(full, pats, k, TmsMode::PHI, positions, report, packed);
                            std::vector<TmsIndex*> xs = {&phi, &full_phi, &phi_built};
                            if (report != TmsReport::SMEM_ALL) xs.push_back(&phi_noinv);
                            for (TmsIndex* x : xs) {
                                if (query_output(*x, pats, k, TmsMode::PHI, positions, report, packed) != want_phi) {
                                    std::cout << "  FAILED: phi layout output differs, k " << k << std::endl;
                                    assert(false && "a phi index must answer as the full index does");
                                    return false;
                                }
                                ++compared;
                            }
                        }
                }
            }

            // Refusals.
            const TmsParts phiskip_need = tms_query_parts(TmsMode::PHISKIP, false, TmsReport::MS);
            const TmsParts smem_one_psi = tms_query_parts(TmsMode::PSI, false, TmsReport::SMEM_ONE);
            const TmsParts pos_psi = tms_query_parts(TmsMode::PSI, true, TmsReport::MS);
            assert(throws<std::runtime_error>([&] { round_trip(phi_built, &psi_need); }, "needs psi"));
            assert(throws<std::runtime_error>([&] { round_trip(phi_built, &phiskip_need); }, "needs psi"));
            assert(throws<std::runtime_error>([&] { round_trip(psi_built, &phi_need); }, "needs phi"));
            assert(throws<std::runtime_error>([&] { round_trip(psi_built, &smem_one_psi); }, "needs phi"));
            assert(throws<std::runtime_error>([&] { round_trip(psi_built, &pos_psi); }, "needs phi"));
            assert(throws<std::runtime_error>([&] { round_trip(phi_noinv_built, &phi_need); }, "needs phi_inv"));
            std::vector<std::vector<ulint>> len, pos;
            assert(throws<std::invalid_argument>([&] { tms_query_batch(phi, pats, 8, len, TmsMode::PSI); }, "psi"));
            assert(throws<std::invalid_argument>([&] { tms_query_batch(phi, pats, 8, len, TmsMode::DUAL); }, "psi"));
            assert(throws<std::invalid_argument>([&] { tms_query(phi, "ACGT"); }, "psi"));
            assert(throws<std::invalid_argument>([&] { tms_query_batch(psi, pats, 8, len, TmsMode::PHI); }, "phi"));
            assert(throws<std::invalid_argument>([&] { tms_query_batch(psi, pats, 8, len, TmsMode::PSI, &pos); }, "phi"));
            tms_query_batch(phi_noinv, pats, 8, len, TmsMode::PHI, &pos);
            std::vector<TmsSmemHits> hits;
            assert(throws<std::invalid_argument>(
                [&] { tms_report_smems_batch(phi_noinv, len, pos, 0, TmsReport::SMEM_ALL, 8, hits); }, "phi_inv"));
            assert(throws<std::invalid_argument>([&] { TmsIndex(in.heads, in.lens, psi_opts, &in.tops); }, "psi layout"));
            assert(throws<std::invalid_argument>([&] { TmsIndex(in.heads, in.lens, phi_opts); }, "phi layout"));
        }
    }
    std::cout << "  " << compared << " batches PASSED" << std::endl;
    return true;
}

bool run_all_tests(const std::string& data_dir) {
    bool all_ran = true;
    test_orbit_structures();
    std::cout << std::endl;
    test_tms_vs_naive();
    std::cout << std::endl;
    if (!test_tms_vs_ms_minishred(data_dir)) all_ran = false;
    std::cout << std::endl;
    test_tms_batch(data_dir);
    std::cout << std::endl;
    test_phi_inv();
    std::cout << std::endl;
    test_smem_detection();
    std::cout << std::endl;
    test_smem_report(data_dir);
    std::cout << std::endl;
    test_layouts(data_dir);
    std::cout << std::endl;
    std::cout << (all_ran ? "All tms_test checks PASSED" : "SOME TMS TESTS NOT RUN") << std::endl;
    return all_ran;
}

}  // namespace tms_test
