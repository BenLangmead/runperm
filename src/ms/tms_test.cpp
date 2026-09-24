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
    size_t checked = 0, positions = 0;
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
            std::vector<std::vector<ulint>> want;
            for (const auto& P : pats) want.push_back(tms_query(idx, P));
            std::vector<std::vector<ulint>> first_pos;
            for (TmsMode mode : modes) {
                for (size_t k : {1, 3, 32, 1000}) {
                    for (bool with_pos : {false, true}) {
                        std::vector<std::vector<ulint>> got, pos;
                        tms_query_batch(k == 3 ? loaded : idx, pats, k, got, mode, with_pos ? &pos : nullptr);
                        if (got != want) {
                            std::cout << "  FAILED lengths: mode " << int(mode) << ", k " << k << std::endl;
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
    std::cout << "  " << checked << " batches and " << positions << " positions PASSED" << std::endl;
    return true;
}

namespace {

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
    }
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
    std::cout << "  " << count << " texts PASSED" << std::endl;
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
 * are those tms_smems gives.
 */
bool test_smem_report(const std::string& data_dir) {
    std::cout << "Testing SMEM reports against every occurrence" << std::endl;
    std::mt19937 rng(43);
    const auto inputs = batch_inputs(data_dir, rng);
    const TmsMode modes[] = {TmsMode::PSI, TmsMode::PHI, TmsMode::PHISKIP, TmsMode::DUAL};
    size_t smems = 0, positions = 0;
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
                            std::vector<TmsSmemHits> hits(pats.size());
                            for (size_t j = 0; j < pats.size(); ++j)
                                tms_report_smems(idx, len[j], pos[j], min_len, report, hits[j]);
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
                                    std::vector<ulint> got(hits[j].pos.begin() + at, hits[j].pos.begin() + at + h.count);
                                    at += h.count;
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
    std::cout << "  " << smems << " SMEMs and " << positions << " positions PASSED" << std::endl;
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
    std::cout << (all_ran ? "All tms_test checks PASSED" : "SOME TMS TESTS NOT RUN") << std::endl;
    return all_ran;
}

}  // namespace tms_test
