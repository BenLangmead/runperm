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
#include "ms_build.hpp"
#include "ms_io.hpp"
#include "orient.hpp"
#include "ftab_build.hpp"
#include "ftab_file.hpp"
#include <algorithm>
#include <array>
#include <set>
#include <cassert>
#include <cstring>
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

}  // namespace

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

namespace {

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

/** LF split in several ways, and phi split with balancing 2. */
std::vector<TmsBuildOptions> split_variants() {
    TmsBuildOptions a;
    TmsBuildOptions b;
    b.phi_split = orbit::split_params(std::nullopt, 2);
    TmsBuildOptions c;
    c.lf_split = orbit::split_params(orbit::DEFAULT_LENGTH_CAPPING, std::nullopt);
    TmsBuildOptions d;
    d.lf_split = orbit::split_params(std::nullopt, 2);
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
        if (ms_io::load_rlbwt_lcps(data_dir + "/minishred1_20_002", heads, lens, lcps)) {
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
    if (!ms_io::load_rlbwt_lcps(data_dir + "/minishred1_20_002", heads, lens, lcps)) {
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
 * TmsIndex::PackedAccess fits an index with psi and agrees with the index's
 * own column reads for every row: length, character, the LF step from every
 * offset, resolving a position past a row's end, psi_at_tail and
 * psi_at_head, and the psi step from every offset and from one past the
 * row's end.  psi inverts LF on every position, and its character is the
 * BWT character of the row psi reaches, the F character.  Returns the
 * number of rows checked.
 */
size_t check_packed_access(TmsIndex& idx) {
    assert(idx.has_psi());
    const TmsIndex::PackedAccess pa = idx.packed_access();
    assert(pa.fits() && "the test indexes must fit packed access");
    std::array<uchar, 256> lf_code, fl_code;
    for (size_t c = 0; c < 256; ++c) {
        lf_code[c] = idx.lf_code(static_cast<uchar>(c));
        fl_code[c] = idx.fl_code(static_cast<uchar>(c));
        assert((lf_code[c] != TmsIndex::no_code) == idx.occurs(static_cast<uchar>(c)));
        assert((fl_code[c] != TmsIndex::no_code) == idx.occurs(static_cast<uchar>(c)));
    }
    const ulint rows = idx.move_runs();
    assert(idx.psi_intervals() == rows);
    std::vector<ulint> start(rows + 1, 0);
    for (ulint i = 0; i < rows; ++i) start[i + 1] = start[i] + idx.get_length(i);
    auto absolute = [&](TmsIndex::LFPos p) { return start[p.interval] + p.offset; };
    for (ulint i = 0; i < rows; ++i) {
        const auto r = pa.row(i);
        const ulint len = idx.get_length(i);
        assert(pa.length(r) == len);
        assert(pa.code(r) == lf_code[idx.get_character(i)]);
        const TmsIndex::FLPos t = pa.psi_at_tail(i), want_t = idx.psi_at_tail(i);
        assert(t.interval == i && t.offset + 1 == len && want_t.interval == i && want_t.offset == t.offset);
        const TmsIndex::FLPos h = pa.psi_at_head(i), want_h = idx.psi_at_head(i);
        assert(h.interval == i && h.offset == 0 && want_h.interval == i && want_h.offset == 0);
        for (ulint o = 0; o <= len; ++o) {
            if (o == len && i + 1 == rows) break;
            TmsIndex::LFPos p{};
            p.interval = i;
            p.offset = o;
            if (o < len) {
                const TmsIndex::LFPos got = pa.start_LF(r, p), want = idx.start_LF(p);
                assert(got.interval == want.interval && got.offset == want.offset);
            } else {
                TmsIndex::LFPos q = p;
                const TmsIndex::LFPos want = idx.finish_LF(p);
                const auto rr = pa.resolve_LF(q);
                assert(q.interval == want.interval && q.offset == want.offset && rr == pa.row(want.interval));
            }
            uchar c = 0, want_c = 0;
            const TmsIndex::FLPos got = pa.psi_step(p, c), want = idx.psi_step(p, want_c);
            assert(got.interval == want.interval && got.offset == want.offset && c == fl_code[want_c]);
            // psi then LF returns to the position, and the F character is
            // the BWT character of the row psi reaches.
            const TmsIndex::LFPos y = idx.resolve_psi(got);
            assert(idx.get_character(y.interval) == want_c);
            assert(absolute(idx.LF_step(y)) == absolute(idx.finish_LF(p)));
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

/**
 * --orient on texts that hold both strands: each text is sequences and
 * their reverse complements, separated, and patterns are drawn from it and
 * from its reverse complement, with mutations, N and absent bytes.  The
 * matching statistics from matching every read's reverse complement (rc),
 * and from flipping reads by random base qualities (qual), converted back,
 * must equal brute force and the forward query, for the ms index and for
 * tms indexes in every mode.
 */
bool test_orient_both_strands() {
    std::cout << "Testing --orient rc and qual on texts holding both strands" << std::endl;
    std::mt19937 rng(73);
    const TmsMode modes[] = {TmsMode::PSI, TmsMode::PHI, TmsMode::PHISKIP, TmsMode::DUAL};
    size_t checked = 0, flipped = 0;
    for (int k = 0; k < 30; ++k) {
        std::string text;
        const int count = 1 + rng() % 4;
        for (int i = 0; i < count; ++i) {
            std::string w(1 + rng() % (k % 4 == 0 ? 8 : 300), 'A');
            for (auto& c : w) c = "ACGT"[rng() % (k % 5 == 0 ? 2 : 4)];
            if (k % 3 == 0 && w.size() > 10) w[rng() % w.size()] = 'N';
            std::string r = w;
            revcomp_in_place(r);
            text += (text.empty() ? "" : "%") + w + "%" + r;
        }
        TextBwt t = make_text_bwt(text);
        std::vector<ulint> tops;
        for (const auto& l : t.lcps_per_run) tops.push_back(l[0]);
        auto pats = fuzz_patterns(t.str, rng, 40);
        for (auto P : fuzz_patterns(t.str, rng, 40)) {
            revcomp_in_place(P);
            pats.push_back(P);
        }
        std::vector<std::string> quals;
        for (const auto& P : pats) {
            std::string q(P.size(), 'I');
            for (auto& c : q) c = static_cast<char>('!' + rng() % 41);
            quals.push_back(q);
        }
        std::vector<std::vector<ulint>> want;
        for (const auto& P : pats) want.push_back(naive_ms(t.str, P));

        // Runs query(seqs, out) on the reads as oriented by o and converts
        // the results back.
        auto oriented = [&](ReadOrient o, auto query) {
            std::vector<std::string> seqs = pats;
            std::vector<char> flips;
            flipped += orient_reads(seqs, &quals, o, flips);
            std::vector<std::vector<ulint>> out;
            query(seqs, out);
            std::vector<ulint> scratch;
            restore_orientation(flips, out, scratch);
            return out;
        };
        auto check = [&](const std::vector<std::vector<ulint>>& got, const char* what) {
            if (got != want) {
                std::cout << "  FAILED: " << what << ", text of length " << text.size() << std::endl;
                assert(false && "oriented matching statistics must equal brute force");
                return false;
            }
            checked += pats.size();
            return true;
        };
        auto ms_idx = ms_index_for(t);
        for (ReadOrient o : {ReadOrient::FWD, ReadOrient::RC, ReadOrient::QUAL}) {
            for (size_t kk : {1, 8}) {
                auto got = oriented(o, [&](const std::vector<std::string>& seqs, std::vector<std::vector<ulint>>& out) {
                    ms_query_batch(ms_idx, seqs, kk, out);
                });
                if (!check(got, "ms index")) return false;
            }
            auto one = oriented(o, [&](const std::vector<std::string>& seqs, std::vector<std::vector<ulint>>& out) {
                for (const auto& P : seqs) out.push_back(ms_query(ms_idx, P));
            });
            if (!check(one, "ms_query")) return false;
        }
        for (const auto& opts : split_variants()) {
            TmsIndex idx(t.heads, t.lens, opts, &tops);
            for (TmsMode mode : modes)
                for (ReadOrient o : {ReadOrient::FWD, ReadOrient::RC, ReadOrient::QUAL})
                    for (size_t kk : {1, 8}) {
                        auto got =
                            oriented(o, [&](const std::vector<std::string>& seqs, std::vector<std::vector<ulint>>& out) {
                                tms_query_batch(idx, seqs, kk, out, mode, nullptr);
                            });
                        if (!check(got, "tms index")) return false;
                    }
        }
    }
    std::cout << "  " << checked << " reads (" << flipped << " flipped) PASSED" << std::endl;
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
    phi_only.load(again, TmsParts{true, true, false, true});
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
 * The oracle modes: reposition outcomes recorded by every other mode, for
 * several numbers of patterns in flight, are the same, and replaying them
 * in the oracle and oracle-scan modes, with any number in flight, packed or
 * column access, with and without positions, through a file and on an index
 * loaded with LF alone, gives the lengths and positions of the other modes.
 * Replays that do not fit the query are refused, and oracle-scan refuses an
 * interval its scan does not find.
 */
bool test_tms_oracle(const std::string& data_dir) {
    std::cout << "Testing the oracle mode against the recording modes" << std::endl;
    std::mt19937 rng(53);
    const auto inputs = batch_inputs(data_dir, rng);
    const TmsMode modes[] = {TmsMode::PSI, TmsMode::PHI, TmsMode::PHISKIP, TmsMode::DUAL};
    size_t replays = 0, decisions = 0;
    for (const auto& in : inputs) {
        auto pats = fuzz_patterns(in.text, rng, in.text.size() > 5000 ? 80 : 150);
        for (const auto& o : {split_variants()[0], split_variants()[3]}) {
            TmsIndex idx(in.heads, in.lens, o, &in.tops);
            std::vector<std::vector<ulint>> want, want_pos;
            tms_query_batch(idx, pats, 8, want, TmsMode::PSI, &want_pos);
            std::vector<std::vector<RepositionDecision>> first;
            for (TmsMode mode : modes)
                for (size_t k : {1, 3, 32}) {
                    std::vector<std::vector<ulint>> len;
                    std::vector<std::vector<RepositionDecision>> rec;
                    tms_query_batch(idx, pats, k, len, mode, nullptr, true, nullptr, nullptr, &rec);
                    assert(len == want);
                    assert(rec.size() == pats.size());
                    if (first.empty()) first = rec;
                    for (size_t j = 0; j < rec.size(); ++j) {
                        if (rec[j].size() != first[j].size()) {
                            std::cout << "  FAILED: recorded outcomes differ in mode " << int(mode) << ", k " << k
                                      << std::endl;
                            assert(false && "recorded outcomes must not depend on mode or k");
                            return false;
                        }
                        for (size_t i = 0; i < rec[j].size(); ++i)
                            assert(rec[j][i].interval == first[j][i].interval && rec[j][i].len == first[j][i].len);
                    }
                }
            RepositionOracle oracle;
            oracle.move_runs = idx.move_runs();
            oracle.append(first);
            decisions += oracle.decisions.size();
            // Through the file format.
            std::stringstream file;
            oracle.write(file);
            RepositionOracle read;
            read.read(file);
            assert(read.move_runs == oracle.move_runs && read.offsets == oracle.offsets);
            assert(read.decisions.size() == oracle.decisions.size());
            const TmsParts lf_need = tms_query_parts(TmsMode::ORACLE, false, TmsReport::MS);
            const TmsParts pos_need = tms_query_parts(TmsMode::ORACLE, true, TmsReport::MS);
            assert(!lf_need.psi && !lf_need.phi && !pos_need.psi && pos_need.phi && pos_need.phi_starts);
            TmsIndex lf_only = round_trip(idx, &lf_need), with_pos = round_trip(idx, &pos_need);
            assert(!lf_only.has_psi() && !lf_only.has_phi());
            for (TmsMode replay_mode : {TmsMode::ORACLE, TmsMode::ORACLE_SCAN})
            for (size_t k : {1, 2, 3, 32, 1000})
                for (int access = 0; access < 4; ++access) {
                    const bool packed = access < 2, positions = access % 2 == 1;
                    // Replay a slice of the patterns too, as a block of reads does.
                    const size_t first_pat = k == 3 ? pats.size() / 3 : 0;
                    const std::vector<std::string> sub(pats.begin() + first_pat, pats.end());
                    const RepositionOracleView view = RepositionOracleView::of(read, first_pat, sub.size());
                    std::vector<std::vector<ulint>> len, pos;
                    tms_query_batch(positions ? with_pos : lf_only, sub, k, len, replay_mode,
                                    positions ? &pos : nullptr, packed, nullptr, &view);
                    bool ok = len.size() == sub.size();
                    for (size_t j = 0; ok && j < sub.size(); ++j) {
                        ok = len[j] == want[first_pat + j];
                        if (ok && positions) ok = pos[j] == want_pos[first_pat + j];
                    }
                    if (!ok) {
                        std::cout << "  FAILED replay: mode " << int(replay_mode) << ", k " << k << ", packed "
                                  << packed << ", positions " << positions << std::endl;
                        assert(false && "the oracle mode must give the other modes' output");
                        return false;
                    }
                    ++replays;
                }
            // Refusals: no oracle, an oracle for other intervals, and too few
            // or too many outcomes for a pattern.
            std::vector<std::vector<ulint>> len;
            assert(throws<std::invalid_argument>([&] { tms_query_batch(lf_only, pats, 8, len, TmsMode::ORACLE); },
                                                 "replay"));
            RepositionOracleView bad = RepositionOracleView::of(read, 0, pats.size());
            bad.move_runs += 1;
            assert(throws<std::invalid_argument>(
                [&] { tms_query_batch(lf_only, pats, 8, len, TmsMode::ORACLE, nullptr, true, nullptr, &bad); }, "LF intervals"));
            size_t with_rep = 0;
            while (with_rep < pats.size() && first[with_rep].empty()) ++with_rep;
            if (with_rep < pats.size()) {
                // A pattern's outcomes without its last one, and with its
                // last one twice.
                const std::vector<std::string> one{pats[with_rep]};
                RepositionOracle fewer, more;
                fewer.move_runs = more.move_runs = read.move_runs;
                fewer.append({std::vector<RepositionDecision>(first[with_rep].begin(), first[with_rep].end() - 1)});
                auto extra = first[with_rep];
                extra.push_back(extra.back());
                more.append({extra});
                const RepositionOracleView fv = RepositionOracleView::of(fewer, 0, 1), mv = RepositionOracleView::of(more, 0, 1);
                assert(throws<std::runtime_error>(
                    [&] { tms_query_batch(lf_only, one, 1, len, TmsMode::ORACLE, nullptr, true, nullptr, &fv); }, "fewer"));
                assert(throws<std::runtime_error>(
                    [&] { tms_query_batch(lf_only, one, 1, len, TmsMode::ORACLE, nullptr, true, nullptr, &mv); }, "more"));
                for (TmsMode m : {TmsMode::ORACLE, TmsMode::ORACLE_SCAN}) {
                    assert(throws<std::invalid_argument>(
                        [&] { tms_query_batch(lf_only, one, 1, len, m); }, "replay"));
                    assert(throws<std::runtime_error>(
                        [&] { tms_query_batch(lf_only, one, 1, len, m, nullptr, true, nullptr, &fv); }, "fewer"));
                }
                // oracle-scan refuses any interval but the one its scan
                // finds, the nearest candidate on the recorded side: here
                // the first outcome's interval moved one above and one below.
                for (int delta : {-1, 1}) {
                    auto moved = first[with_rep];
                    if (delta < 0 && moved[0].interval == 0) continue;
                    moved[0].interval += delta;
                    RepositionOracle bad_one;
                    bad_one.move_runs = read.move_runs;
                    bad_one.append({moved});
                    const RepositionOracleView bv = RepositionOracleView::of(bad_one, 0, 1);
                    assert(throws<std::runtime_error>(
                        [&] { tms_query_batch(lf_only, one, 1, len, TmsMode::ORACLE_SCAN, nullptr, true, nullptr, &bv); },
                        "oracle"));
                }
            }
            assert(throws<std::runtime_error>([&] { RepositionOracleView::of(read, 1, pats.size()); }, "fewer patterns"));
        }
    }
    std::cout << "  " << decisions << " recorded repositions and " << replays << " replays PASSED" << std::endl;
    return true;
}

/**
 * Indexes built with the psi and phi layouts, and full indexes loaded with
 * only the parts a query needs, give the same output as the full index for
 * every query they can serve, with packed and column access and several
 * numbers of patterns in flight.  The LF columns a layout omits have width
 * 0, and the parts it omits are absent.  Queries that need a missing part are
 * refused, at load and at query time, and so are layouts given the wrong LCP
 * input, and unsplit LF.
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
            // Parts held, and LF columns of width 0 exactly where a part is
            // missing, PHI_OFF, whose phi points all have offset 0, and the
            // LCP summary columns, built without --lcp-summary.
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
                const bool phi_off = c == PSI_INT + 3, summary_col = c == PSI_INT + 4 || c == PSI_INT + 5;
                assert((fw[c] > 0) == !(phi_off || summary_col));
                // The phi layout's LF has other rows, without the union
                // with the images, so only its zero widths are compared.
                assert(sw[c] == (phi_col ? 0 : fw[c]) && (hw[c] == 0) == (psi_col || phi_off || summary_col));
            }
            assert(psi.packed_access().fits() && phi.packed_access().fits());

            // Full indexes loaded with only what psi MS, and what phi queries, read.
            const TmsParts psi_need = tms_query_parts(TmsMode::PSI, false, TmsReport::MS);
            const TmsParts phi_need = tms_query_parts(TmsMode::PHI, true, TmsReport::SMEM_ALL);
            TmsIndex full_psi = round_trip(full, &psi_need), full_phi = round_trip(full, &phi_need);
            assert(parts_are(full_psi, true, false, false) && parts_are(full_phi, false, true, true));
            // Matching statistics alone load phi's lengths-based form, and
            // positions and SMEMs its starts-based form; a whole load has both.
            const TmsParts phi_ms_need = tms_query_parts(TmsMode::PHI, false, TmsReport::MS);
            assert(!phi_ms_need.phi_starts && phi_need.phi_starts);
            TmsIndex full_phi_len = round_trip(full, &phi_ms_need), phi_len = round_trip(phi_built, &phi_ms_need);
            for (const TmsIndex* x : {&full_phi_len, &phi_len})
                assert(parts_are(*x, false, true, false) && x->has_phi_lengths() && !x->has_phi_starts());
            assert(full_phi.has_phi_starts() && !full_phi.has_phi_lengths());
            assert(phi.has_phi_starts() && phi.has_phi_lengths() && full.has_phi_starts() && full.has_phi_lengths());
            assert(full_phi_len.phi_intervals() == full.phi_intervals());

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
                            if (!positions && report == TmsReport::MS) {
                                xs.push_back(&full_phi_len);
                                xs.push_back(&phi_len);
                            }
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
            assert(throws<std::invalid_argument>([&] { tms_query_batch(phi_len, pats, 8, len, TmsMode::PHI, &pos); },
                                                 "starts-based"));
            assert(throws<std::logic_error>([&] { round_trip(phi_len); }, "only some of its parts"));
            tms_query_batch(phi_noinv, pats, 8, len, TmsMode::PHI, &pos);
            std::vector<TmsSmemHits> hits;
            assert(throws<std::invalid_argument>(
                [&] { tms_report_smems_batch(phi_noinv, len, pos, 0, TmsReport::SMEM_ALL, 8, hits); }, "phi_inv"));
            assert(throws<std::invalid_argument>([&] { TmsIndex(in.heads, in.lens, psi_opts, &in.tops); }, "psi layout"));
            assert(throws<std::invalid_argument>([&] { TmsIndex(in.heads, in.lens, phi_opts); }, "phi layout"));
            for (TmsLayout layout : {TmsLayout::FULL, TmsLayout::PSI, TmsLayout::PHI}) {
                TmsBuildOptions unsplit;
                unsplit.layout = layout;
                unsplit.lf_split = orbit::NO_SPLITTING;
                const std::vector<ulint>* tops = layout == TmsLayout::PSI ? nullptr : &in.tops;
                assert(throws<std::invalid_argument>([&] { TmsIndex(in.heads, in.lens, unsplit, tops); }, "LF must be split"));
            }
        }
    }
    std::cout << "  " << compared << " batches PASSED" << std::endl;
    return true;
}

/**
 * The psi and full layouts, from LF split in several ways, give the same
 * output for every query as the full index from LF split with Orbit's
 * default: every mode, with and without positions, every report, packed and
 * column access, and 1 and 32 patterns in flight, also after serialization
 * and a partial load.  Their rows are the union of the split LF intervals'
 * starts and their images' starts, so there are as many as the phi
 * layout's split LF intervals (whose LF has no union) and at most twice
 * that.
 */
bool test_lf_rows(const std::string& data_dir) {
    std::cout << "Testing LF rows from LF split several ways" << std::endl;
    std::mt19937 rng(53);
    const auto inputs = batch_inputs(data_dir, rng);
    const TmsMode modes[] = {TmsMode::PSI, TmsMode::PHI, TmsMode::PHISKIP, TmsMode::DUAL};
    size_t compared = 0, access_rows = 0, rows = 0, split_rows = 0;
    for (const auto& in : inputs) {
        const bool big = in.text.size() > 5000;
        auto pats = fuzz_patterns(in.text, rng, big ? 60 : 100);
        TmsIndex ref(in.heads, in.lens, TmsBuildOptions{}, &in.tops);
        for (const auto& sp : {orbit::split_params(std::nullopt, 2), orbit::split_params(orbit::DEFAULT_LENGTH_CAPPING, std::nullopt),
                               orbit::split_params(2.0, 64), orbit::split_params(std::nullopt, 1024)}) {
            TmsBuildOptions o, psi_opts, phi_opts;
            o.lf_split = psi_opts.lf_split = phi_opts.lf_split = sp;
            psi_opts.layout = TmsLayout::PSI;
            phi_opts.layout = TmsLayout::PHI;
            TmsIndex full_built(in.heads, in.lens, o, &in.tops), psi_built(in.heads, in.lens, psi_opts);
            TmsIndex phi(in.heads, in.lens, phi_opts, &in.tops);
            TmsIndex full = round_trip(full_built), psi = round_trip(psi_built);
            const TmsParts psi_need = tms_query_parts(TmsMode::PSI, false, TmsReport::MS);
            TmsIndex full_psi = round_trip(full_built, &psi_need);
            // Every mode's matching statistics walk phi's lengths-based form
            // when it is present, so a load with only the starts-based form
            // checks that both forms agree.
            const TmsParts starts_need{true, true, false, true};
            TmsIndex full_starts = round_trip(full_built, &starts_need);
            assert(full_starts.has_phi_starts() && !full_starts.has_phi_lengths());
            assert(full.move_runs() == psi.move_runs() && full.move_runs() >= phi.move_runs() &&
                   full.move_runs() <= 2 * phi.move_runs());
            rows += full.move_runs();
            split_rows += phi.move_runs();
            access_rows += check_packed_access(full) + check_packed_access(psi_built);
            for (const auto& P : pats) assert(tms_query(full, P) == tms_query(ref, P) && tms_query(psi, P) == tms_query(ref, P));
            for (size_t k : {1, 32}) {
                for (bool packed : {true, false}) {
                    const auto want_psi = query_output(ref, pats, k, TmsMode::PSI, false, TmsReport::MS, packed);
                    for (TmsIndex* x : {&psi, &full_psi, &psi_built}) {
                        assert(query_output(*x, pats, k, TmsMode::PSI, false, TmsReport::MS, packed) == want_psi);
                        ++compared;
                    }
                    for (TmsMode mode : modes)
                        for (bool positions : {false, true})
                            for (TmsReport report : {TmsReport::MS, TmsReport::SMEM_ONE, TmsReport::SMEM_ALL}) {
                                if (big && report == TmsReport::SMEM_ONE) continue;
                                const auto want = query_output(ref, pats, k, mode, positions, report, packed);
                                const bool differs = query_output(full, pats, k, mode, positions, report, packed) != want ||
                                                     (!positions && report == TmsReport::MS &&
                                                      query_output(full_starts, pats, k, mode, false, report, packed) != want);
                                if (differs) {
                                    std::cout << "  FAILED: output depends on the LF split, mode " << int(mode) << ", k "
                                              << k << ", packed " << packed << std::endl;
                                    assert(false && "output must not depend on the LF split");
                                    return false;
                                }
                                ++compared;
                            }
                }
            }
        }
    }
    std::cout << "  " << access_rows << " rows of packed access, " << compared << " batches PASSED (" << rows
              << " rows for " << split_rows << " split LF intervals)" << std::endl;
    return true;
}

namespace {

/** A text of several strings, each ended by its own terminator. */
struct MultiText {
    std::string str;  // the strings, each followed by '$'
    std::vector<uchar> text;
    std::vector<ulint> sa, isa, lcp;
    std::vector<uchar> heads;  // maximal runs: adjacent terminators form one run
    std::vector<ulint> lens;
};

/**
 * The text S_0 $ S_1 $ ... S_{k-1} $ with distinct terminators, $_i < $_j for
 * i < j, all smaller than every other character, as in a ropebwt3 FMD and as
 * TeraLCP numbers its positions.  LCP never extends over a terminator.
 */
MultiText make_multi_text(const std::vector<std::string>& strings) {
    MultiText t;
    for (const auto& w : strings) t.str += w + "$";
    const size_t n = t.str.size();
    t.text.resize(n);
    for (size_t i = 0; i < n; ++i)
        t.text[i] = t.str[i] == '$' ? orbit::TERMINATOR : t.str[i] == '%' ? orbit::SEPARATOR : static_cast<uchar>(t.str[i]);
    const auto& x = t.text;
    t.sa.resize(n);
    std::iota(t.sa.begin(), t.sa.end(), 0);
    std::sort(t.sa.begin(), t.sa.end(), [&](ulint a, ulint b) {
        while (x[a] == x[b] && x[a] != orbit::TERMINATOR) { ++a; ++b; }
        // Two terminators compare by string, that is by position.
        if (x[a] == orbit::TERMINATOR && x[b] == orbit::TERMINATOR) return a < b;
        return x[a] < x[b];
    });
    t.isa.resize(n);
    for (size_t j = 0; j < n; ++j) t.isa[t.sa[j]] = j;
    t.lcp.assign(n, 0);
    for (size_t j = 1; j < n; ++j) {
        ulint a = t.sa[j - 1], b = t.sa[j], l = 0;
        while (x[a + l] == x[b + l] && x[a + l] != orbit::TERMINATOR) ++l;
        t.lcp[j] = l;
    }
    for (size_t j = 0; j < n; ++j) {
        const uchar c = x[t.sa[j] == 0 ? n - 1 : t.sa[j] - 1];
        if (j == 0 || c != t.heads.back()) {
            t.heads.push_back(c);
            t.lens.push_back(0);
        }
        ++t.lens.back();
    }
    return t;
}

/** The samples TeraLCP -ophi writes for a BWT with this suffix array and LCP array. */
TmsPhiSamples samples_for(const std::vector<uchar>& heads, const std::vector<ulint>& lens,
                          const std::vector<ulint>& sa, const std::vector<ulint>& lcp) {
    std::vector<std::array<ulint, 3>> recs;
    ulint row = 0;
    for (size_t r = 0; r < heads.size(); ++r) {
        if (heads[r] == orbit::TERMINATOR)
            for (ulint j = 0; j < lens[r]; ++j) recs.push_back({sa[row + j], sa[row + j], lcp[row + j]});
        else
            recs.push_back({sa[row], sa[row + lens[r] - 1], lcp[row]});
        row += lens[r];
    }
    ulint max_top = 0;
    for (const auto& x : recs) max_top = std::max(max_top, x[2]);
    TmsPhiSamples s(recs.size(), sa.size(), max_top);
    for (size_t i = 0; i < recs.size(); ++i) {
        s.head.set(i, recs[i][0]);
        s.tail.set(i, recs[i][1]);
        s.top.set(i, recs[i][2]);
    }
    return s;
}

/** samples without its last record. */
TmsPhiSamples drop_last(const TmsPhiSamples& samples, ulint n) {
    const size_t m = samples.size() - 1;
    ulint max_top = 0;
    for (size_t i = 0; i < m; ++i) max_top = std::max<ulint>(max_top, samples.top[i]);
    TmsPhiSamples s(m, n, max_top);
    for (size_t i = 0; i < m; ++i) {
        s.head.set(i, samples.head[i]);
        s.tail.set(i, samples.tail[i]);
        s.top.set(i, samples.top[i]);
    }
    return s;
}

/** phi, PLCP, phi_inv and PLCPB of idx at every text position, against the suffix and LCP arrays. */
void check_phi_against(TmsIndex& idx, const std::vector<ulint>& sa, const std::vector<ulint>& isa,
                       const std::vector<ulint>& lcp) {
    const ulint n = sa.size();
    for (ulint p = 0; p < n; ++p) {
        const ulint j = isa[p];
        auto a = idx.phi_at(p);
        assert(a.idx == p && idx.plcp(a) == lcp[j] && "PLCP at a located phi point");
        assert(idx.phi(a).idx == sa[j == 0 ? n - 1 : j - 1] && "phi(SA[j]) = SA[j - 1]");
        if (!idx.has_phi_inv()) continue;
        auto b = idx.phi_inv_at(p);
        assert(b.idx == p && idx.plcpb(b) == (j + 1 < n ? lcp[j + 1] : 0) && "PLCPB is the LCP with the row below");
        assert(idx.phi_inv(b).idx == sa[j + 1 < n ? j + 1 : 0] && "phi_inv(SA[j]) = SA[j + 1]");
    }
}

}  // namespace

namespace {
/**
 * Check that every counting method gives the walk's SMEM_ALL hits for the
 * statistics len and toeholds pos of pats on idx, for several caps and
 * numbers in flight and, for ADAPTIVE, budgets and length thresholds.  The
 * counter's row-start table samples every fourth interval, so that counts
 * of intervals more than four intervals apart read it.  Returns the
 * comparisons made.
 */
size_t check_counting(TmsIndex& idx, const std::vector<std::string>& pats, const std::vector<std::vector<ulint>>& len,
                      const std::vector<std::vector<ulint>>& pos, ulint min_len) {
    const TmsCounter ctr(idx, 2);
    std::vector<TmsCountOptions> opts;
    TmsCountOptions o;
    o.method = TmsCount::BACKWARD;
    opts.push_back(o);
    o.method = TmsCount::ADAPTIVE;
    o.budget_per_len = 0;
    for (ulint budget : {0, 1, 3, 20})
        for (ulint short_len : {0, 4}) {
            o.budget = budget;
            o.short_len = short_len;
            opts.push_back(o);
        }
    o.short_len = 0;
    for (double per_len : {0.25, 1.0}) {
        o.budget = 2;
        o.budget_per_len = per_len;
        opts.push_back(o);
    }
    const std::vector<std::vector<ulint>> none(pats.size());
    size_t compared = 0;
    for (ulint cap : {ulint(0), ulint(1), ulint(2), ulint(5), TMS_ALL_POSITIONS})
        for (size_t k : {size_t(0), size_t(1), size_t(7), size_t(32)}) {
            std::vector<TmsSmemHits> want, got;
            tms_report_smems_batch(idx, len, pos, min_len, TmsReport::SMEM_ALL, k, want, cap);
            for (const auto& co : opts) {
                const TmsCounting counting{co, &ctr, &pats};
                tms_report_smems_batch(idx, len, co.method == TmsCount::BACKWARD ? none : pos, min_len,
                                       TmsReport::SMEM_ALL, k, got, cap, &counting);
                for (size_t j = 0; j < pats.size(); ++j) {
                    std::string a, b;
                    tms_format_smems(want[j], TmsReport::SMEM_ALL, a);
                    tms_format_smems(got[j], TmsReport::SMEM_ALL, b);
                    if (a != b) {
                        std::cout << "  FAILED: counting method " << int(co.method) << " budget " << co.budget
                                  << " short_len " << co.short_len << " cap " << cap << " k " << k << " on "
                                  << pats[j] << ":\n    walk     " << a << "\n    counting " << b << std::endl;
                        assert(false && "every counting method must give the walk's counts and positions");
                        return compared;
                    }
                    ++compared;
                }
            }
        }
    return compared;
}
}  // namespace

/**
 * SeqIdTable's bucket lookup locates every text position as a search of all
 * strand starts does, giving the same contig, strand and offset and the same
 * label, and failing with the same message for a match that leaves its
 * strand, on random tables with long, short and one-base contigs, at
 * random positions and at every strand's first and last positions.  A table
 * extended after its lookup is built, and a table made with add alone, are
 * searched in full.
 */
bool test_seq_id_lookup() {
    std::cout << "Testing the sequence-ID table's bucket lookup against a full search" << std::endl;
    std::mt19937_64 rng(71);
    size_t checked = 0;
    for (int trial = 0; trial < 60; ++trial) {
        SeqIdTable plain;
        const size_t contigs = 1 + rng() % (trial % 4 == 0 ? 3 : 400);
        for (size_t k = 0; k < contigs; ++k) {
            const uint64_t len = rng() % 3 == 0 ? 1 + rng() % 4 : rng() % 5 == 0 ? 1 + rng() % 1000000 : 1 + rng() % 500;
            plain.add("c" + std::to_string(k) + std::string(rng() % 30, 'x'), len);
        }
        SeqIdTable fast = plain;
        fast.build_lookup();
        assert(plain.bucket.empty() && !fast.bucket.empty());
        auto same = [&](uint64_t p, uint64_t len) {
            std::string a, b;
            bool ea = false, eb = false;
            try { plain.append_location(a, p, len); } catch (const std::exception& e) { ea = true; a = e.what(); }
            try { fast.append_location(b, p, len); } catch (const std::exception& e) { eb = true; b = e.what(); }
            if (ea != eb || a != b) {
                std::cout << "  FAILED: position " << p << " length " << len << ": " << a << " | " << b << std::endl;
                assert(false && "the bucket lookup must locate as the full search does");
            }
            ++checked;
        };
        const uint64_t n = fast.text_length();
        for (int i = 0; i < 2000; ++i) same(rng() % n, 1 + rng() % 5);
        for (size_t j = 0; j < fast.starts.size(); ++j) {
            const uint64_t first = fast.starts[j], last = first + fast.lens[j / 2];  // last is the terminator
            for (uint64_t p : {first, first + 1, last - 1, last})
                if (p < n) same(p, 1);
        }
        // Extending the table drops its lookup until it is built again.
        fast.add("extra", 7);
        plain.add("extra", 7);
        assert(fast.bucket.empty());
        same(fast.starts.back() + 3, 2);
    }
    std::cout << "  " << checked << " positions PASSED" << std::endl;
    return true;
}

/**
 * Counting SMEM occurrences by backward search and adaptively gives the
 * same counts and listed positions as the walk, on the fuzz texts and
 * minishred with one terminator, in every split variant, and on string
 * collections with many terminators, whose phi and phi_inv come from
 * samples.  The inputs make some top and bottom ends move to other runs and
 * some counts read the row-start table, and the adaptive method both walk
 * and switch.
 */
bool test_smem_count(const std::string& data_dir) {
    std::cout << "Testing SMEM counting by backward search against the walk" << std::endl;
    std::mt19937 rng(47);
    tms_count_stats = TmsCountStats{};
    size_t compared = 0;
    for (const auto& in : batch_inputs(data_dir, rng)) {
        auto pats = fuzz_patterns(in.text, rng, in.text.size() > 5000 ? 40 : 60);
        for (const auto& o : split_variants()) {
            TmsIndex idx(in.heads, in.lens, o, &in.tops);
            std::vector<std::vector<ulint>> len, pos;
            tms_query_batch(idx, pats, 8, len, TmsMode::PSI, &pos);
            for (ulint min_len : {0, 5}) compared += check_counting(idx, pats, len, pos, min_len);
        }
    }
    for (int k = 0; k < 20; ++k) {
        std::vector<std::string> strings;
        const int count = 1 + rng() % 12;
        for (int i = 0; i < count; ++i) {
            if (!strings.empty() && rng() % 4 == 0) {
                strings.push_back(strings[rng() % strings.size()]);
                continue;
            }
            std::string w(1 + rng() % (k % 3 == 0 ? 6 : 80), 'A');
            for (auto& c : w) c = "ACGT"[rng() % (k % 5 == 0 ? 2 : 4)];
            strings.push_back(w);
        }
        MultiText t = make_multi_text(strings);
        std::vector<std::string> pats;
        for (std::string q : fuzz_patterns(t.str, rng, 30)) {
            for (auto& c : q)
                if (c == '$') c = 'A';
            pats.push_back(q);
        }
        for (TmsBuildOptions o : split_variants()) {
            TmsIndex idx(t.heads, t.lens, o, samples_for(t.heads, t.lens, t.sa, t.lcp));
            std::vector<std::vector<ulint>> len, pos;
            tms_query_batch(idx, pats, 8, len, TmsMode::PSI, &pos, true);
            compared += check_counting(idx, pats, len, pos, 0);
        }
    }
    const TmsCountStats& st = tms_count_stats;
    assert(st.top_moves > 0 && st.bot_moves > 0 && st.far_ends > 0 && st.switched > 0 && st.walked > 0 &&
           st.list_steps > 0 && "the inputs must exercise every part of counting");
    std::cout << "  " << compared << " read listings agree (" << st.searched << " searches, " << st.top_moves
              << " top and " << st.bot_moves << " bottom moves, " << st.far_ends << " counts from the table, "
              << st.switched << " adaptive switches) PASSED" << std::endl;
    return true;
}

/**
 * Indexes whose phi and phi_inv come from per-run suffix array and LCP
 * samples (as TeraLCP -ophi writes them) instead of a walk of LF.  With one
 * terminator they answer every query exactly as indexes built by the walk, in
 * the full and phi layouts and for each split variant, also after
 * serialization.  With many terminators, where rank-based LF does not follow
 * the text across a terminator and so the walk cannot be used, phi and
 * phi_inv are right at every text position, matching statistics equal naive
 * ones, every reported position is an occurrence, and SMEM_ALL lists every
 * occurrence.  Samples that do not fit the RLBWT, and samples with the psi
 * layout, are refused.
 */
namespace {
/** A stream buffer over a string that cannot seek, as a pipe. */
struct NoSeekBuf : std::streambuf {
    explicit NoSeekBuf(std::string& s) { setg(&s[0], &s[0], &s[0] + s.size()); }
};
}  // namespace

/**
 * The ORBITIDX header of a tms index: its sections agree with the index, a
 * load reads the sections its query needs and refuses a header whose column
 * codes disagree with the code, a load skips sections on a stream that cannot
 * seek, and a truncated file is refused.
 */
bool test_index_format() {
    std::cout << "Testing the ORBITIDX header of tms indexes" << std::endl;
    using namespace orbit_idx;
    TextBwt t = make_text_bwt("ACGTTGCAACGGATTACAGATTACAACGTACGTTTGCA");
    TmsIndex idx(t.heads, t.lens, TmsBuildOptions{}, samples_for(t.heads, t.lens, t.sa, t.lcp));
    idx.set_info("built_by=test\n");
    std::stringstream ss;
    idx.serialize(ss);
    const std::string file = ss.str();
    const Header h = decode_header(file.data());
    assert(h.tool == Tool::TMS && h.layout == Layout::FULL && h.n == idx.domain() && h.bwt_runs == idx.bwt_runs());
    assert(h.sections.size() == 5 && h.info == "built_by=test\n");
    const SectionEntry* lf = h.find(Section::LF);
    assert(lf && lf->rows == idx.lf().intervals() && lf->cols.size() == 10);
    assert(lf->cols[4].first == static_cast<uint8_t>(Col::PSI_INT) && lf->cols[4].second == idx.lf().get_widths()[4]);
    for (const auto& e : h.sections) assert(e.offset % SECTION_ALIGN == 0);
    assert(file.size() == h.sections.back().offset + h.sections.back().bytes);
    // What a query loads.
    const uint32_t lf_bit = section_bit(Section::LF);
    assert(TmsIndex::sections_for(TmsParts{}) == lf_bit);
    assert(TmsIndex::sections_for(TmsParts{true, false, false, false}) == (lf_bit | section_bit(Section::F_STARTS)));
    assert(TmsIndex::sections_for(TmsParts{false, true, false, false}) == (lf_bit | section_bit(Section::PHI_LEN)));
    assert(TmsIndex::sections_for(TmsParts{false, true, true, true}) ==
           (lf_bit | section_bit(Section::PHI_START) | section_bit(Section::PHI_INV)));
    assert(TmsIndex::sections_for(TmsParts{false, false, true, false}) == lf_bit);
    {
        std::stringstream in(file);
        TmsIndex x;
        x.load(in, TmsParts{false, true, true, true});
        const SectionEntry* st = h.find(Section::PHI_START);
        const SectionEntry* inv = h.find(Section::PHI_INV);
        assert(x.has_phi_starts() && !x.has_phi_lengths() && x.has_phi_inv());
        assert(x.start_table_bytes() == TmsIndex::start_table_bytes(st->domain, st->rows) +
                                            TmsIndex::start_table_bytes(inv->domain, inv->rows));
    }

    const TmsParts phi_need{false, true, false, false};
    // A column code other than the one this build reads.
    std::string bad_col = file;
    bad_col[FIXED_BYTES + 48 + 2 * 4] = static_cast<char>(Col::PLCP);
    {
        std::stringstream in(bad_col);
        TmsIndex x;
        assert(throws<std::runtime_error>([&] { x.load(in, phi_need); }, "has column plcp where this build reads psi_int"));
    }
    // Skipping sections on a stream that cannot seek reads past them.
    for (const TmsParts& need : {phi_need, TmsParts{false, true, true, true}, TmsParts{true, false, false, false}}) {
        std::string copy = file;
        NoSeekBuf buf(copy);
        std::istream in(&buf);
        TmsIndex piped, seeked;
        piped.load(in, need);
        std::stringstream again(file);
        seeked.load(again, need);
        std::ostringstream a, b;
        piped.describe(a);
        seeked.describe(b);
        assert(a.str() == b.str());
    }
    // A truncated file.
    {
        std::stringstream in(file.substr(0, file.size() - 100));
        TmsIndex x;
        assert(throws<std::runtime_error>([&] { x.load(in); }, "truncated"));
    }
    std::cout << "  " << h.sections.size() << " sections PASSED" << std::endl;
    return true;
}

bool test_phi_import() {
    std::cout << "Testing phi and phi_inv from imported samples" << std::endl;
    std::mt19937 rng(61);
    size_t compared = 0, multi_checked = 0;
    const TmsMode modes[] = {TmsMode::PSI, TmsMode::PHI, TmsMode::PHISKIP, TmsMode::DUAL};
    for (const auto& s : fuzz_texts(rng)) {
        TextBwt t = make_text_bwt(s);
        std::vector<ulint> tops;
        for (const auto& l : t.lcps_per_run) tops.push_back(l[0]);
        auto pats = fuzz_patterns(t.str, rng, 30);
        for (TmsBuildOptions o : split_variants()) {
            for (TmsLayout layout : {TmsLayout::FULL, TmsLayout::PHI}) {
                o.layout = layout;
                TmsIndex walked(t.heads, t.lens, o, &tops),
                    built(t.heads, t.lens, o, samples_for(t.heads, t.lens, t.sa, t.lcp));
                TmsIndex imported = round_trip(built);
                // Both routes store a phi interval at each row that starts a
                // run, 0 elsewhere, and no offsets: their LFs are the same.
                {
                    const auto &a = walked.lf(), &b = imported.lf();
                    constexpr size_t phi_off = TmsIndex::LF::template data_column<TmsLFCols::PHI_OFF>();
                    assert(a.get_widths() == b.get_widths() && a.get_widths()[phi_off] == 0);
                    assert(a.intervals() == b.intervals());
                    for (ulint k = 0; k < a.intervals(); ++k) {
                        assert(a.get_length(k) == b.get_length(k));
                        assert(a.template get<TmsLFCols::PHI_INT>(k) == b.template get<TmsLFCols::PHI_INT>(k));
                        assert(a.template get<TmsLFCols::PSI_INT>(k) == b.template get<TmsLFCols::PSI_INT>(k));
                        assert(a.template get<TmsLFCols::PSI_OFF>(k) == b.template get<TmsLFCols::PSI_OFF>(k));
                    }
                }
                assert(imported.phi_intervals() == walked.phi_intervals());
                assert(imported.phi_inv_intervals() == walked.phi_inv_intervals());
                check_phi_against(imported, t.sa, t.isa, t.lcp);
                for (TmsMode mode : modes) {
                    if (layout == TmsLayout::PHI && mode != TmsMode::PHI) continue;
                    for (bool positions : {false, true})
                        for (TmsReport report : {TmsReport::MS, TmsReport::SMEM_ONE, TmsReport::SMEM_ALL}) {
                            const auto want = query_output(walked, pats, 8, mode, positions, report, true);
                            if (query_output(imported, pats, 8, mode, positions, report, true) != want) {
                                std::cout << "  FAILED: imported phi answers differently, text of length " << s.size()
                                          << std::endl;
                                assert(false && "imported phi must answer as walked phi does");
                                return false;
                            }
                            ++compared;
                        }
                }
            }
        }
    }

    // Many terminators: string collections with repeats and duplicates.
    for (int k = 0; k < 40; ++k) {
        std::vector<std::string> strings;
        const int count = 1 + rng() % 12;
        for (int i = 0; i < count; ++i) {
            if (!strings.empty() && rng() % 4 == 0) {
                strings.push_back(strings[rng() % strings.size()]);
                continue;
            }
            std::string w(1 + rng() % (k % 3 == 0 ? 6 : 80), 'A');
            for (auto& c : w) c = "ACGT"[rng() % (k % 5 == 0 ? 2 : 4)];
            if (!strings.empty() && rng() % 3 == 0) {
                w = strings[rng() % strings.size()];
                if (!w.empty()) w[rng() % w.size()] = "ACGT"[rng() % 4];
            }
            strings.push_back(w);
        }
        MultiText t = make_multi_text(strings);
        std::vector<std::string> pats;
        for (const auto& P : fuzz_patterns(t.str, rng, 30)) {
            std::string q = P;
            for (auto& c : q)
                if (c == '$') c = 'A';
            pats.push_back(q);
        }
        for (TmsBuildOptions o : split_variants()) {
            for (TmsLayout layout : {TmsLayout::FULL, TmsLayout::PHI}) {
                o.layout = layout;
                TmsIndex built(t.heads, t.lens, o, samples_for(t.heads, t.lens, t.sa, t.lcp));
                TmsIndex idx = round_trip(built);
                check_phi_against(idx, t.sa, t.isa, t.lcp);
                for (TmsMode mode : modes) {
                    if (layout == TmsLayout::PHI && mode != TmsMode::PHI) continue;
                    std::vector<std::vector<ulint>> len, pos;
                    tms_query_batch(idx, pats, 8, len, mode, &pos, true);
                    // Matching statistics alone walk phi's lengths-based form.
                    const TmsParts ms_need = tms_query_parts(mode, false, TmsReport::MS);
                    TmsIndex lengths_only = round_trip(built, &ms_need);
                    assert(lengths_only.has_phi_lengths() == (mode != TmsMode::PSI) && !lengths_only.has_phi_starts());
                    std::vector<std::vector<ulint>> len_lengths;
                    tms_query_batch(lengths_only, pats, 8, len_lengths, mode, nullptr, true);
                    assert(len_lengths == len && "phi's two forms give the same matching statistics");
                    std::vector<TmsSmemHits> hits;
                    tms_report_smems_batch(idx, len, pos, 0, TmsReport::SMEM_ALL, 8, hits);
                    for (size_t j = 0; j < pats.size(); ++j) {
                        const auto want = naive_ms(t.str, pats[j]);
                        assert(len[j] == want && "matching statistics with many terminators");
                        for (size_t i = 0; i < want.size(); ++i)
                            if (want[i] > 0)
                                assert(t.str.compare(pos[j][i], want[i], pats[j], i, want[i]) == 0 &&
                                       "a reported position is an occurrence");
                        size_t at = 0;
                        for (const auto& h : hits[j].smems) {
                            const std::string m = pats[j].substr(h.start, h.len);
                            ulint occ = 0;
                            for (size_t p = t.str.find(m); p != std::string::npos; p = t.str.find(m, p + 1)) ++occ;
                            assert(h.count == occ && h.listed == occ && "SMEM_ALL counts every occurrence");
                            std::vector<ulint> listed(hits[j].pos.begin() + at, hits[j].pos.begin() + at + h.listed);
                            at += h.listed;
                            std::sort(listed.begin(), listed.end());
                            assert(std::adjacent_find(listed.begin(), listed.end()) == listed.end());
                            for (ulint p : listed) assert(t.str.compare(p, h.len, m) == 0);
                        }
                        multi_checked += want.size();
                    }
                }
            }
        }
        // Samples for another RLBWT, and samples with a psi layout, are refused.
        const ulint n = t.text.size();
        assert(throws<std::invalid_argument>(
            [&] { TmsIndex(t.heads, t.lens, TmsBuildOptions{}, drop_last(samples_for(t.heads, t.lens, t.sa, t.lcp), n)); },
            "records"));
        TmsBuildOptions psi_opts;
        psi_opts.layout = TmsLayout::PSI;
        assert(throws<std::invalid_argument>(
            [&] { TmsIndex(t.heads, t.lens, psi_opts, samples_for(t.heads, t.lens, t.sa, t.lcp)); }, "psi layout"));
    }
    std::cout << "  " << compared << " batches against walked phi, " << multi_checked
              << " values with many terminators PASSED" << std::endl;
    return true;
}

namespace {

/** The ORBITIDX header of an index, as its file would hold it. */
template <class Index>
orbit_idx::Header header_of(Index& idx) {
    std::stringstream ss;
    idx.serialize(ss);
    return ftab::index_header_of(ss);
}

/** A table written to its file format and read back, with its header in h. */
ftab::Ftab ftab_round_trip(const ftab::Ftab& tab, const orbit_idx::Header& ih, ftab::FileHeader& h) {
    std::stringstream ss;
    ftab::FileHeader w = ftab::make_header(ih, tab);
    w.present = tab.present();
    w.info = "test=1\n";
    ftab::write(ss, w, tab);
    h = ftab::read_header(ss);
    ftab::Ftab back = ftab::read_payload(ss, h, true);
    assert(back.table_words() == tab.table_words() && back.window_words() == tab.window_words());
    assert(h.present == w.present && h.info == w.info && ftab::link_mismatch(h, ih).empty());
    return back;
}

/** Number of distinct j-mers over ACGT in T, for each j from 0 to k. */
std::vector<uint64_t> brute_present(const std::string& T, unsigned k) {
    std::vector<uint64_t> out(k + 1, 0);
    out[0] = 1;
    for (unsigned j = 1; j <= k; ++j) {
        std::set<std::string> seen;
        for (size_t i = 0; i + j <= T.size(); ++i) {
            const std::string w = T.substr(i, j);
            if (w.find_first_not_of("ACGT") == std::string::npos) seen.insert(w);
        }
        out[j] = seen.size();
    }
    return out;
}

/**
 * Lengths and positions from tms_query_batch with table ft in every mode
 * idx serves, with and without positions, packed and column access and
 * several patterns in flight, against want and want_pos (positions without
 * a table).  Returns the batches compared.
 */
size_t check_tms_with_ftab(TmsIndex& idx, const std::vector<std::string>& pats, const ftab::Ftab& ft,
                           const std::vector<std::vector<ulint>>& want, const std::vector<std::vector<ulint>>* want_pos,
                           const std::string& what) {
    size_t checked = 0;
    for (TmsMode mode : {TmsMode::PSI, TmsMode::PHI, TmsMode::PHISKIP, TmsMode::DUAL}) {
        if ((mode != TmsMode::PHI && !idx.has_psi()) || (mode != TmsMode::PSI && !idx.has_phi())) continue;
        for (size_t k : {1, 3, 32})
            for (int access = 0; access < 4; ++access) {
                const bool packed = access < 2, with_pos = access % 2 == 1;
                if (with_pos && !want_pos) continue;
                std::vector<std::vector<ulint>> got, pos;
                tms_query_batch(idx, pats, k, got, mode, with_pos ? &pos : nullptr, packed, &ft);
                if (got != want || (with_pos && pos != *want_pos)) {
                    std::cout << "  FAILED tms (" << what << "): ftab k " << ft.k() << ", mode " << int(mode) << ", k " << k
                              << ", packed " << packed << ", positions " << with_pos << std::endl;
                    assert(false && "tms_query_batch with an ftab must match it without");
                    return checked;
                }
                ++checked;
            }
    }
    // Recording takes no table, which skips some of the repositions.
    if (idx.has_psi()) {
        std::vector<std::vector<ulint>> got;
        std::vector<std::vector<RepositionDecision>> rec;
        assert(throws<std::invalid_argument>(
            [&] { tms_query_batch(idx, pats, 8, got, TmsMode::PSI, nullptr, true, &ft, nullptr, &rec); }, "ftab"));
    }
    return checked;
}

}  // namespace

/**
 * Tables of several k, built on one thread and on several, for the tms
 * indexes (full, psi and phi layouts, several splittings) of each fuzz text
 * and minishred and the ms index of each fuzz text, each written to the
 * file format and read back.  The walk's counts of occurring j-mers at every
 * depth match brute force; tables built on different thread counts are the
 * same bits; and the batched queries with the table give exactly what they
 * give without it: lengths in every mode the layout serves, and positions
 * where it has phi, with packed and column access.  An index loaded with
 * phi's lengths-based form alone builds a table with a toehold and no window
 * positions, which serves the modes that walk phi without positions.
 * Patterns are the fuzz patterns plus exact substrings of each length up to
 * 14, so that table hits, reads of exactly k characters and shorter ones all
 * occur.
 */
bool test_ftab(const std::string& data_dir) {
    std::cout << "Testing ftab against queries without it" << std::endl;
    std::mt19937 rng(43);
    const auto inputs = batch_inputs(data_dir, rng);
    size_t checked = 0, hits = 0, tables = 0;
    for (const auto& in : inputs) {
        const std::string& T = in.text;
        auto pats = fuzz_patterns(T, rng, in.text.size() > 5000 ? 100 : 200);
        const size_t tlen = T.size() - 1;
        for (size_t len = 1; len <= 14; ++len)
            for (int r = 0; r < 4 && len <= tlen; ++r) pats.push_back(T.substr(rng() % (tlen - len + 1), len));
        std::vector<TmsIndex> idxs;
        for (const auto& o : split_variants()) idxs.emplace_back(in.heads, in.lens, o, &in.tops);
        for (size_t v : {0, 1}) {
            TmsBuildOptions psi_opts = split_variants()[v], phi_opts = split_variants()[v];
            psi_opts.layout = TmsLayout::PSI;
            phi_opts.layout = TmsLayout::PHI;
            TmsIndex psi_built(in.heads, in.lens, psi_opts), phi_built(in.heads, in.lens, phi_opts, &in.tops);
            idxs.push_back(round_trip(psi_built));
            idxs.push_back(round_trip(phi_built));
        }
        std::vector<std::vector<ulint>> want, want_pos;
        for (const auto& P : pats) want.push_back(tms_query(idxs[0], P));
        {
            std::vector<std::vector<ulint>> got;
            tms_query_batch(idxs[0], pats, 32, got, TmsMode::PSI, &want_pos);
            assert(got == want);
        }
        const std::vector<uint64_t> bf = brute_present(T, 5);
        for (auto& idx : idxs) {
            const orbit_idx::Header ih = header_of(idx);
            for (unsigned fk : {1u, 2u, 5u, 9u}) {
                const ftab::Built b1 = ftab::build_tms(idx, fk, 1), b4 = ftab::build_tms(idx, fk, 4);
                assert(b1.tab.table_words() == b4.tab.table_words() && b1.tab.window_words() == b4.tab.window_words());
                assert(b1.present == b4.present && b1.present[fk] == b1.tab.present());
                assert(b1.tab.has_toehold() == idx.has_phi() && b1.tab.has_window() == idx.has_phi());
                for (unsigned j = 0; j <= std::min(fk, 5u); ++j)
                    if (b1.present[j] != bf[j]) {
                        std::cout << "  FAILED tms: occurring " << j << "-mers differ from brute force" << std::endl;
                        assert(false && "the walk must count exactly the occurring j-mers");
                        return false;
                    }
                ftab::FileHeader h;
                const ftab::Ftab ft = ftab_round_trip(fk % 2 ? b4.tab : b1.tab, ih, h);
                ++tables;
                checked += check_tms_with_ftab(idx, pats, ft, want, idx.has_phi() ? &want_pos : nullptr, "built");
                for (const auto& P : pats) {
                    ftab::Entry e;
                    hits += ft.lookup(ft.key(P.data(), P.size()), e);
                }
                // A query that needs a toehold or window positions refuses a table without them.
                if (!idx.has_phi()) {
                    std::vector<std::vector<ulint>> got, pos;
                    bool refused = false;
                    try {
                        TmsIndex& full = idxs[0];
                        tms_query_batch(full, pats, 8, got, TmsMode::PHI, nullptr, true, &ft);
                    } catch (const std::invalid_argument&) {
                        refused = true;
                    }
                    assert(refused && "a mode that walks phi needs a table with a toehold");
                }
            }
            // phi's lengths-based form alone: a toehold without window positions.
            if (idx.has_phi()) {
                TmsParts need = idx.parts();
                need.phi_starts = false;
                need.psi = false;
                TmsIndex len_only = round_trip(idx, &need);
                const ftab::Built b = ftab::build_tms(len_only, 5, 3);
                assert(b.tab.has_toehold() && !b.tab.has_window());
                checked += check_tms_with_ftab(len_only, pats, b.tab, want, nullptr, "lengths form");
                bool refused = false;
                try {
                    std::vector<std::vector<ulint>> got, pos;
                    tms_query_batch(idx, pats, 8, got, TmsMode::PHI, &pos, true, &b.tab);
                } catch (const std::invalid_argument&) {
                    refused = true;
                }
                assert(refused && "positions need a table with window positions");
            }
        }
    }
    // The ms index over each fuzz text: table counts and batch results.
    std::mt19937 rng2(47);
    for (const auto& s : fuzz_texts(rng2)) {
        TextBwt t = make_text_bwt(s);
        auto idx = ms_index_for(t);
        const orbit_idx::Header ih = header_of(idx);
        auto pats = fuzz_patterns(t.str, rng2, 200);
        const size_t tlen = t.str.size() - 1;
        for (size_t len = 1; len <= 14; ++len)
            for (int r = 0; r < 4 && len <= tlen; ++r) pats.push_back(t.str.substr(rng2() % (tlen - len + 1), len));
        std::vector<std::vector<ulint>> want;
        for (const auto& P : pats) want.push_back(ms_query(idx, P));
        for (unsigned fk : {1u, 2u, 5u, 12u}) {
            const ftab::Built b = ftab::build_ms(idx, fk, fk % 2 ? 4 : 1);
            if (fk <= 5 && b.tab.present() != ftab_survey::brute_force(idx, fk)) {
                std::cout << "  FAILED ms: occurring " << fk << "-mers differ from brute force" << std::endl;
                assert(false && "ftab must hold exactly the occurring k-mers");
                return false;
            }
            assert(b.present[fk] == b.tab.present() && !b.tab.has_toehold());
            ftab::FileHeader h;
            const ftab::Ftab ft = ftab_round_trip(b.tab, ih, h);
            ++tables;
            for (size_t k : {1, 3, 32}) {
                std::vector<std::vector<ulint>> got, got_narrow;
                ms_query_batch(idx, pats, k, got, &ft);
                ms_query_batch_impl(idx, idx.column_access(), pats, k, got_narrow, &ft);
                if (got != want || got_narrow != want) {
                    std::cout << "  FAILED ms: ftab k " << fk << ", k " << k << std::endl;
                    assert(false && "ms_query_batch with an ftab must match ms_query");
                    return false;
                }
                ++checked;
            }
        }
    }
    std::cout << "  " << tables << " tables, " << checked << " batches, " << hits << " tms table hits PASSED" << std::endl;
    return true;
}

/**
 * An ftab's link to its index: it matches the index it was built from, and
 * one whose header differs only in its info text (as a rewrap's does), and
 * is refused, naming the field, by an index of the same text with other
 * split settings or layout and by an index of another text.
 */
bool test_ftab_link(const std::string& data_dir) {
    std::cout << "Testing ftab links to their indexes" << std::endl;
    std::mt19937 rng(53);
    const auto inputs = batch_inputs(data_dir, rng);
    size_t refused = 0;
    auto mismatch = [&](const ftab::FileHeader& h, const orbit_idx::Header& ih, const std::string& field) {
        const std::string bad = ftab::link_mismatch(h, ih);
        if (bad.empty() || bad.find(field) == std::string::npos) {
            std::cout << "  FAILED: expected a mismatch in " << field << ", got \"" << bad << "\"" << std::endl;
            assert(false && "a link must name the field that differs");
            return false;
        }
        ++refused;
        return true;
    };
    for (size_t a = 0; a + 1 < inputs.size(); ++a) {
        const Input& in = inputs[a];
        const Input& other = inputs[a + 1];
        if (in.text.size() < 50 || other.text.size() < 50) continue;
        auto variants = split_variants();
        TmsIndex x(in.heads, in.lens, variants[0], &in.tops), y(in.heads, in.lens, variants[3], &in.tops);
        TmsIndex z(other.heads, other.lens, variants[0], &other.tops);
        TmsBuildOptions phi_split = variants[1], psi_opts = variants[0];
        psi_opts.layout = TmsLayout::PSI;
        TmsIndex w(in.heads, in.lens, phi_split, &in.tops), p(in.heads, in.lens, psi_opts);
        const orbit_idx::Header hx = header_of(x);
        ftab::FileHeader h;
        ftab_round_trip(ftab::build_tms(x, 4, 2).tab, hx, h);
        // Its own index, and the same with other info text.
        orbit_idx::Header rewrapped = hx;
        rewrapped.info = "rewrapped_from=elsewhere\n";
        if (!ftab::link_mismatch(h, rewrapped).empty()) {
            std::cout << "  FAILED: a rewrapped header must match" << std::endl;
            assert(false);
            return false;
        }
        const orbit_idx::Header hy = header_of(y), hz = header_of(z), hw = header_of(w), hp = header_of(p);
        if (hy.find(orbit_idx::Section::LF)->crc != hx.find(orbit_idx::Section::LF)->crc && !mismatch(h, hy, "section LF"))
            return false;
        if (!mismatch(h, hw, "section ")) return false;
        if (!mismatch(h, hp, "layout")) return false;
        if (!mismatch(h, hz, other.text.size() != in.text.size() ? "n" : "")) return false;
        // ms: the same text with another LF split, and another text.
        if (in.text.find_first_not_of("ACGT$") == std::string::npos) {
            TextBwt t = make_text_bwt(in.text.substr(0, in.text.size() - 1));
            auto mi = ms_index_for(t);
            const orbit_idx::Header hm = header_of(mi);
            ftab::FileHeader hmf;
            ftab_round_trip(ftab::build_ms(mi, 3, 1).tab, hm, hmf);
            if (!mismatch(hmf, hx, "tool") || !mismatch(h, hm, "tool")) return false;
        }
    }
    std::cout << "  " << refused << " mismatches named PASSED" << std::endl;
    return true;
}

/**
 * The auto rule: the largest k whose table fits in the budget, then down
 * while occupancy is under the minimum, on constructed sizes and counts and
 * on the tms indexes of the fuzz texts.
 */
bool test_ftab_auto(const std::string& data_dir) {
    std::cout << "Testing ftab-build --auto's choice of k" << std::endl;
    ftab::Shape s;
    s.interval = 10;
    s.offset = 5;
    for (unsigned k = 1; k <= 12; ++k) {
        s.k = k;
        const double bytes = double(s.payload_bytes());
        assert(ftab::largest_fitting_k(s, bytes) == k);
        assert(ftab::largest_fitting_k(s, bytes - 1) == std::max(1u, k - 1));
    }
    s.k = 0;
    assert(ftab::largest_fitting_k(s, 0) == 1);
    auto present_with = [](std::vector<double> occ) {
        std::vector<uint64_t> p(occ.size());
        for (size_t j = 0; j < occ.size(); ++j) p[j] = uint64_t(occ[j] * double(uint64_t{1} << (2 * j)) + 0.5);
        return p;
    };
    assert(ftab::step_down(6, present_with({1, 1, 1, 0.9, 0.8, 0.79, 0.5}), 0.8) == 4);
    assert(ftab::step_down(6, present_with({1, 1, 1, 0.9, 0.85, 0.8, 0.8}), 0.8) == 6);
    assert(ftab::step_down(6, present_with({1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5}), 0.8) == 1);
    std::mt19937 rng(59);
    const auto inputs = batch_inputs(data_dir, rng);
    size_t cases = 0;
    for (const auto& in : inputs) {
        if (in.text.size() < 100) continue;
        TmsIndex idx(in.heads, in.lens, split_variants()[0], &in.tops);
        const ftab::Shape base = ftab::tms_shape(idx, 2);
        const std::vector<uint64_t> bf = brute_present(in.text, 7);
        for (unsigned target : {3u, 5u, 7u}) {
            ftab::Shape t = base;
            t.k = target;
            const uint64_t index_bytes = uint64_t(double(t.payload_bytes()) / 0.05) + 1;
            ftab::AutoChoice choice;
            auto build = [&](const ftab::Shape& sh) { return ftab::build_tms(idx, sh, 2); };
            const ftab::Built b = ftab::build_auto(base, index_bytes, ftab::AutoRule{}, build, choice);
            unsigned want = target;
            while (want > 1 && double(bf[want]) < 0.8 * double(uint64_t{1} << (2 * want))) --want;
            if (choice.k_fit != target || choice.k != want || b.tab.k() != want) {
                std::cout << "  FAILED auto: fit " << choice.k_fit << " (want " << target << "), k " << choice.k
                          << " (want " << want << ")" << std::endl;
                assert(false && "--auto must follow its rule");
                return false;
            }
            assert(b.tab.table_words() == ftab::build_tms(idx, want, 1).tab.table_words());
            ++cases;
        }
    }
    std::cout << "  " << cases << " indexes PASSED" << std::endl;
    return true;
}

/**
 * Exact LCP summaries of every LF row of an index built with heads, lens
 * and o: the LCP at the row's head and the least of its other rows' (the
 * greatest LCP of all for a row of one), from each run's full LCPs.
 */
TmsLcpSummary exact_summaries(const std::vector<uchar>& heads, const std::vector<ulint>& lens,
                              const std::vector<std::vector<ulint>>& lcps, const TmsBuildOptions& o) {
    const orbit::int_vector rows = TmsIndex::lf_row_lengths(heads, lens, o);
    std::vector<ulint> all;
    for (const auto& l : lcps) all.insert(all.end(), l.begin(), l.end());
    const ulint big = *std::max_element(all.begin(), all.end());
    const uchar w = std::max<uchar>(1, orbit::bit_width(big));
    TmsLcpSummary sum{orbit::int_vector(rows.size(), w), orbit::int_vector(rows.size(), w)};
    for (size_t k = 0, at = 0; k < rows.size(); at += rows[k++]) {
        sum.top.set(k, all[at]);
        ulint m = big;
        for (ulint j = 1; j < rows[k]; ++j) m = std::min(m, all[at + j]);
        sum.int_min.set(k, m);
    }
    return sum;
}

/**
 * tms_query_batch in PHISUM gives the lengths and positions PHI gives, for
 * several numbers of patterns in flight, with packed and column access,
 * after a serialize and load, in the full and phi layouts, and under
 * splittings that cut runs into many LF rows.  The summaries come from each
 * run's minima-only pairs, as from a minima file, and also, exactly, from
 * the full LCPs.  PHISUM refuses an index without summaries.  Runs on the
 * fuzz texts and, if its LF is a single cycle, on minishred.
 */
bool test_tms_phisum(const std::string& data_dir) {
    std::cout << "Testing tms_query_batch in phisum mode against phi" << std::endl;
    std::mt19937 rng(43);
    struct LcpInput { std::vector<uchar> heads; std::vector<ulint> lens; std::vector<std::vector<ulint>> lcps; std::string str; };
    std::vector<LcpInput> inputs;
    for (const auto& s : fuzz_texts(rng)) {
        TextBwt t = make_text_bwt(s);
        inputs.push_back({t.heads, t.lens, t.lcps_per_run, t.str});
    }
    {
        LcpInput in;
        if (ms_io::load_rlbwt_lcps(data_dir + "/minishred1_20_002", in.heads, in.lens, in.lcps)) {
            std::vector<ulint> tops;
            for (const auto& l : in.lcps) tops.push_back(l[0]);
            try {
                TmsIndex probe(in.heads, in.lens, TmsBuildOptions{}, &tops);
                // The text read back with LF, as batch_inputs reads it.
                auto pos = probe.first();
                const ulint n = probe.domain();
                in.str.assign(n, ' ');
                for (ulint st = 0; st < n; ++st) {
                    const uchar c = probe.get_character(pos.interval);
                    in.str[(2 * n - 2 - st) % n] =
                        c == orbit::TERMINATOR ? '$' : c == orbit::SEPARATOR ? '%' : static_cast<char>(c);
                    pos = probe.LF_step(pos);
                }
                inputs.push_back(std::move(in));
            } catch (const std::exception& e) {
                std::cout << "  minishred skipped: " << e.what() << std::endl;
            }
        }
    }
    auto variants = split_variants();
    for (double cap : {1.0, 0.25}) {
        TmsBuildOptions o;
        o.lf_split = orbit::split_params(cap, std::nullopt);
        variants.push_back(o);
    }
    size_t checked = 0, cont_rows = 0;
    bool refused = false;
    for (const auto& in : inputs) {
        std::vector<ulint> tops;
        for (const auto& l : in.lcps) tops.push_back(l[0]);
        const auto pats = fuzz_patterns(in.str, rng, 120);
        for (TmsBuildOptions o : variants) {
            for (TmsLayout layout : {TmsLayout::FULL, TmsLayout::PHI}) {
                o.layout = layout;
                o.lcp_summary = nullptr;
                TmsIndex plain(in.heads, in.lens, o, &tops);
                if (!refused) {
                    std::vector<std::vector<ulint>> len;
                    assert(throws<std::invalid_argument>([&] { tms_query_batch(plain, pats, 8, len, TmsMode::PHISUM); },
                                                         "--lcp-summary"));
                    refused = true;
                }
                std::vector<std::vector<ulint>> want, want_pos;
                tms_query_batch(plain, pats, 8, want, TmsMode::PHI, &want_pos);
                const orbit::int_vector rows = TmsIndex::lf_row_lengths(in.heads, in.lens, o);
                auto cut = ms_build::row_summaries(in.lens, rows, ms_build::lcp_vectors_source(in.lcps, true));
                const TmsLcpSummary from_minima{std::move(cut.first), std::move(cut.second)};
                const TmsLcpSummary exact = exact_summaries(in.heads, in.lens, in.lcps, o);
                for (const TmsLcpSummary* sum : {&from_minima, &exact}) {
                    o.lcp_summary = sum;
                    TmsIndex idx(in.heads, in.lens, o, &tops);
                    assert(idx.has_lcp_summary() && idx.move_runs() == plain.move_runs());
                    for (ulint k = 1; k < idx.move_runs(); ++k) cont_rows += idx.get_character(k) == idx.get_character(k - 1);
                    std::stringstream ss;
                    idx.serialize(ss);
                    TmsIndex loaded;
                    loaded.load(ss);
                    assert(loaded.has_lcp_summary());
                    for (size_t k : {1, 3, 32}) {
                        for (int access = 0; access < 4; ++access) {
                            const bool packed = access < 2, with_pos = access % 2 == 1;
                            std::vector<std::vector<ulint>> got, pos;
                            tms_query_batch(k == 3 ? loaded : idx, pats, k, got, TmsMode::PHISUM, with_pos ? &pos : nullptr,
                                            packed);
                            if (got != want || (with_pos && pos != want_pos)) {
                                std::cout << "  FAILED: " << (got != want ? "lengths" : "positions") << ", k " << k
                                          << ", packed " << packed << ", layout " << int(layout) << ", summaries "
                                          << (sum == &exact ? "exact" : "from minima") << std::endl;
                                assert(false && "phisum must give phi's lengths and positions");
                                return false;
                            }
                            ++checked;
                        }
                    }
                }
            }
        }
    }
    std::cout << "  " << checked << " batches over " << cont_rows << " rows continuing a run PASSED" << std::endl;
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
    test_tms_phisum(data_dir);
    std::cout << std::endl;
    test_ftab(data_dir);
    std::cout << std::endl;
    test_ftab_link(data_dir);
    std::cout << std::endl;
    test_ftab_auto(data_dir);
    std::cout << std::endl;
    test_orient_both_strands();
    std::cout << std::endl;
    test_phi_inv();
    std::cout << std::endl;
    test_smem_detection();
    std::cout << std::endl;
    test_smem_report(data_dir);
    std::cout << std::endl;
    test_smem_count(data_dir);
    std::cout << std::endl;
    test_seq_id_lookup();
    std::cout << std::endl;
    test_layouts(data_dir);
    std::cout << std::endl;
    test_tms_oracle(data_dir);
    std::cout << std::endl;
    test_lf_rows(data_dir);
    std::cout << std::endl;
    test_phi_import();
    std::cout << std::endl;
    test_index_format();
    std::cout << std::endl;
    std::cout << (all_ran ? "All tms_test checks PASSED" : "SOME TMS TESTS NOT RUN") << std::endl;
    return all_ran;
}

}  // namespace tms_test
