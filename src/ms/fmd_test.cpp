/**
 * Tests for string collections in the FMD form (see fmd_test.hpp) and for
 * BWT lengths past 2^41.
 *
 * A collection's BWT is built here as an FM-index of a collection is: the
 * strings S_1 .. S_k each end in their own terminator, terminator i sorts
 * below terminator j for i < j and below every letter, and the BWT is that of
 * S_1 $_1 S_2 $_2 ... S_k $_k with the first string preceded by $_k.  All
 * terminators get the same head byte, so LF ranks them, and its permutation
 * generally has several cycles.  N sorts after T, as in ropebwt3.
 */

#include "fmd_test.hpp"
#include "ms_rlbwt.hpp"
#include "rlbwt_io.hpp"
#include "serialize.hpp"
#include "tms_index.hpp"
#include "tms_query.hpp"
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace fmd_test {

using uchar = orbit::uchar;
using ulint = orbit::ulint;

namespace {

/** An RLBWT as files hold it (raw head bytes) with per-row LCPs by run. */
struct Coll {
    std::vector<uchar> raw_heads;
    std::vector<ulint> lens;
    std::vector<std::vector<ulint>> lcps_per_run;
    ulint n = 0;
    size_t lf_cycles = 0;
};

/** The head byte fmd2rlbwt writes for a character of a string, or for a terminator. */
uchar head_byte(char c) { return c == '$' ? '$' : c == 'N' ? ms_alpha::TEXT_N : static_cast<uchar>(c); }

/**
 * The collection BWT of strs (over ACGTN) by sorting suffixes.  If n_matches,
 * N matches N when computing LCPs (as TeraLCP does); otherwise LCPs stop at
 * N.  Matching statistics are the same either way, since no read character
 * matches an N.
 */
Coll build_collection(const std::vector<std::string>& strs, bool n_matches) {
    const int k = static_cast<int>(strs.size());
    std::vector<int> t;      // ranks: terminator i is i, then A C G T N
    std::vector<char> ch;    // the character, '$' for terminators
    for (int i = 0; i < k; ++i) {
        for (char c : strs[i]) {
            const int r = c == 'A' ? 0 : c == 'C' ? 1 : c == 'G' ? 2 : c == 'T' ? 3 : 4;
            t.push_back(k + r);
            ch.push_back(c);
        }
        t.push_back(i);
        ch.push_back('$');
    }
    const size_t n = t.size();
    std::vector<size_t> sa(n);
    std::iota(sa.begin(), sa.end(), 0);
    // Distinct terminators end every comparison inside the text.
    std::sort(sa.begin(), sa.end(), [&](size_t a, size_t b) {
        while (t[a] == t[b]) { ++a; ++b; }
        return t[a] < t[b];
    });
    auto same = [&](size_t a, size_t b) {
        return t[a] == t[b] && (n_matches || ch[a] != 'N');
    };
    Coll c;
    c.n = n;
    std::vector<char> bwt(n);
    for (size_t j = 0; j < n; ++j) {
        bwt[j] = ch[sa[j] == 0 ? n - 1 : sa[j] - 1];
        ulint l = 0;
        if (j > 0) while (same(sa[j - 1] + l, sa[j] + l)) ++l;
        if (j == 0 || bwt[j] != bwt[j - 1]) {
            c.raw_heads.push_back(head_byte(bwt[j]));
            c.lens.push_back(0);
            c.lcps_per_run.emplace_back();
        }
        ++c.lens.back();
        c.lcps_per_run.back().push_back(l);
    }
    // Cycles of rank-based LF, with characters in byte order of their heads.
    std::vector<ulint> count(256, 0), C(256, 0), seen(256, 0);
    for (size_t j = 0; j < n; ++j) ++count[head_byte(bwt[j])];
    for (size_t b = 1; b < 256; ++b) C[b] = C[b - 1] + count[b - 1];
    std::vector<ulint> lf(n);
    for (size_t j = 0; j < n; ++j) {
        const uchar b = head_byte(bwt[j]);
        lf[j] = C[b] + seen[b]++;
    }
    std::vector<bool> done(n, false);
    for (size_t j = 0; j < n; ++j) {
        if (done[j]) continue;
        ++c.lf_cycles;
        for (size_t x = j; !done[x]; x = lf[x]) done[x] = true;
    }
    return c;
}

/** Matching statistics of P against the strings; N matches nothing. */
std::vector<ulint> naive_ms(const std::vector<std::string>& strs, const std::string& P) {
    std::vector<ulint> ms(P.size(), 0);
    for (size_t i = 0; i < P.size(); ++i)
        for (const auto& s : strs)
            for (size_t st = 0; st < s.size(); ++st) {
                size_t l = 0;
                while (i + l < P.size() && st + l < s.size() && P[i + l] == s[st + l] && P[i + l] != 'N') ++l;
                ms[i] = std::max<ulint>(ms[i], l);
            }
    return ms;
}

/** Load raw heads and lengths through files, as build-rlbwt and tms-build do. */
void load_through_files(const std::string& dir, const std::vector<uchar>& raw_heads, const std::vector<ulint>& lens,
                        std::vector<uchar>& heads, std::vector<ulint>& lens_out) {
    const std::string hp = dir + "/.fmd_test.heads", lp = dir + "/.fmd_test.len";
    assert(rlbwt_io::write_rlbwt(hp, lp, raw_heads, lens));
    std::string err;
    const bool ok = rlbwt_io::load_rlbwt(hp, lp, heads, lens_out, err);
    std::remove(hp.c_str());
    std::remove(lp.c_str());
    if (!ok) std::cout << "  load_rlbwt: " << err << std::endl;
    assert(ok && "load_rlbwt must accept one terminator per string and N as 0x7f");
    assert(lens_out == lens);
}

/** An ms index from loaded heads and lengths and retained (minima-only) LCP pairs, as build-rlbwt makes it. */
MSIndexSpillLCP<false> build_ms(std::vector<uchar> heads, std::vector<ulint> lens, std::vector<RunLcpPairs> pairs) {
    apply_lf_splitting(heads, lens, pairs, orbit::split_params{}, true);
    auto [run_data, spill, max_top, max_sub, skinny, jumbo] =
        build_spill_data_from_pairs(std::move(pairs), 0.9, true, false, 0, 6);
    return MSIndexSpillLCP<false>(heads, lens, run_data, std::move(spill), max_top, max_sub, 0, 6);
}

/** A psi-layout tms index, serialized and loaded back. */
TmsIndex build_tms_psi(const std::vector<uchar>& heads, const std::vector<ulint>& lens) {
    TmsBuildOptions o;
    o.layout = TmsLayout::PSI;
    TmsIndex built(heads, lens, o);
    std::stringstream ss;
    built.serialize(ss);
    TmsIndex loaded;
    loaded.load(ss);
    assert(loaded.domain() == built.domain() && loaded.move_runs() == built.move_runs());
    return loaded;
}

/** Check every way of computing matching statistics against want, for each pattern. */
void check_queries(MSIndexSpillLCP<false>& ms, TmsIndex& tms, const std::vector<std::string>& pats,
                   const std::vector<std::vector<ulint>>& want, const char* what) {
    for (size_t j = 0; j < pats.size(); ++j) {
        const auto a = ms_query(ms, pats[j]);
        const auto b = tms_query(tms, pats[j]);
        if (a != want[j] || b != want[j]) {
            std::cout << "  FAILED (" << what << ") on pattern " << pats[j] << std::endl;
            assert(false && "matching statistics must equal the naive ones");
        }
    }
    for (size_t k : {1, 3, 32}) {
        std::vector<std::vector<ulint>> got;
        ms_query_batch(ms, pats, k, got);
        assert(got == want && "ms_query_batch must equal the naive matching statistics");
        got.clear();
        tms_query_batch(tms, pats, k, got);
        assert(got == want && "tms_query_batch must equal the naive matching statistics");
    }
}

/**
 * Read every string back with LF from the terminator rows, as tms-text does
 * for a collection, and compare with the strings as a multiset.
 */
void check_strings(TmsIndex& tms, std::vector<std::string> strs) {
    ulint m = 0;
    for (ulint i = 0; i < tms.move_runs(); ++i)
        if (tms.get_character(i) == orbit::TERMINATOR) m += tms.get_length(i);
    assert(m == strs.size());
    std::vector<std::string> got;
    auto row = tms.first();
    for (ulint k = 0; k < m; ++k) {
        if (k > 0) {
            if (row.offset + 1 < tms.get_length(row.interval)) ++row.offset;
            else row = tms.down(row);  // the next interval's first row
        }
        std::string s;
        for (auto pos = row;;) {
            const uchar c = tms.get_character(pos.interval);
            if (c == orbit::TERMINATOR) break;
            s.push_back(c == ms_alpha::TEXT_N ? 'N' : static_cast<char>(c));
            pos = tms.LF_step(pos);
        }
        got.emplace_back(s.rbegin(), s.rend());
    }
    std::sort(got.begin(), got.end());
    std::sort(strs.begin(), strs.end());
    if (got != strs) {
        for (auto& g : got) std::cout << "  got  " << g << std::endl;
        for (auto& g : strs) std::cout << "  want " << g << std::endl;
    }
    assert(got == strs && "LF from the terminator rows must read back every string");
}

std::string rand_dna(std::mt19937& rng, size_t len, const char* alpha = "ACGT", size_t k = 4) {
    std::string s(len, 'A');
    for (auto& c : s) c = alpha[rng() % k];
    return s;
}

std::string revcomp(const std::string& s) {
    std::string r(s.rbegin(), s.rend());
    for (auto& c : r) c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : c == 'T' ? 'A' : c;
    return r;
}

/** Random collections: copies, prefixes, suffixes, N runs, reverse complements. */
std::vector<std::string> random_collection(std::mt19937& rng, int trial) {
    std::vector<std::string> strs;
    const int k = std::vector<int>{1, 2, 3, 5, 8, 20}[trial % 6];
    for (int i = 0; i < k; ++i) {
        const int kind = static_cast<int>(rng() % 8);
        std::string s;
        if (strs.empty() || kind < 3) s = rand_dna(rng, 1 + rng() % 200);
        else {
            const std::string& b = strs[rng() % strs.size()];
            if (kind == 3) s = b;
            else if (kind == 4) s = b.substr(0, 1 + rng() % b.size());
            else if (kind == 5) s = b.substr(rng() % b.size());
            else if (kind == 6) s = revcomp(b);
            else {
                s = b;
                for (auto& c : s) if (rng() % 50 == 0) c = "ACGT"[rng() % 4];
            }
        }
        if (rng() % 4 == 0 && s.size() > 4) {
            const size_t at = rng() % s.size(), len = 1 + rng() % std::min<size_t>(20, s.size() - at);
            s.replace(at, len, std::string(len, 'N'));
        }
        if (rng() % 15 == 0) s = std::string(1 + rng() % 5, 'N');
        if (rng() % 15 == 0) s = std::string(1, "ACGT"[rng() % 4]);
        strs.push_back(s);
    }
    if (trial % 3 == 0) {
        const size_t m = strs.size();
        for (size_t i = 0; i < m; ++i) strs.push_back(revcomp(strs[i]));
    }
    return strs;
}

/** Patterns from the strings (some running across a terminator), mutated, some with N. */
std::vector<std::string> patterns_for(std::mt19937& rng, const std::vector<std::string>& strs, int count) {
    std::string all;
    for (const auto& s : strs) all += s + "$";
    std::vector<std::string> pats;
    for (int k = 0; k < count; ++k) {
        const size_t len = k % 30 == 0 ? 0 : k % 30 == 1 ? 1 : rng() % 120;
        std::string P = len <= all.size() ? all.substr(rng() % (all.size() - len + 1), len) : all;
        const double rate = std::vector<double>{0.0, 0.01, 0.05, 0.3}[k % 4];
        for (auto& c : P)
            if (c == '$' || c == 'N' || (rng() % 1000) < rate * 1000) c = "ACGT"[rng() % 4];
        if (k % 5 == 0 && !P.empty()) P[rng() % P.size()] = 'N';
        if (k % 9 == 0) P = rand_dna(rng, rng() % 60);
        pats.push_back(P);
    }
    return pats;
}

void test_collections(const std::string& dir) {
    std::cout << "Testing string collections with one terminator per string and N after T" << std::endl;
    std::mt19937 rng(2026);
    size_t multi_cycle = 0, trials = 0, with_n = 0;
    for (int trial = 0; trial < 120; ++trial) {
        const auto strs = random_collection(rng, trial);
        const bool n_matches = trial % 2 == 0;
        const Coll c = build_collection(strs, n_matches);
        std::vector<uchar> heads;
        std::vector<ulint> lens;
        load_through_files(dir, c.raw_heads, c.lens, heads, lens);
        for (uchar h : heads) with_n += h == ms_alpha::TEXT_N;
        // Minima as a TeraLCP -ominima file would carry them.
        auto pairs = retained_lcp_pairs(c.lcps_per_run, true);
        const std::string mp = dir + "/.fmd_test.min";
        assert(rlbwt_io::write_minima(mp, pairs, c.n));
        std::vector<RunLcpPairs> read_back;
        std::string err;
        const bool ok = rlbwt_io::read_minima(mp, heads.size(), c.n, read_back, err);
        std::remove(mp.c_str());
        assert(ok && read_back == pairs);
        auto ms = build_ms(heads, lens, std::move(read_back));
        auto tms = build_tms_psi(heads, lens);
        check_strings(tms, strs);
        const auto pats = patterns_for(rng, strs, 60);
        std::vector<std::vector<ulint>> want;
        for (const auto& P : pats) want.push_back(naive_ms(strs, P));
        check_queries(ms, tms, pats, want, "collection");
        multi_cycle += c.lf_cycles > 1;
        ++trials;
    }
    assert(multi_cycle > trials / 4 && "the collections must include LF permutations with several cycles");
    assert(with_n > 0);
    std::cout << "  " << trials << " collections PASSED (" << multi_cycle << " with LF of several cycles)" << std::endl;
}

/**
 * The RLBWT of k copies of S, as heads, lengths and, per run, its rows' LCPs
 * as (rows, value) segments, the first of which starts with the run's top.
 * Each suffix x of S$ gives k adjacent rows, one per copy, in terminator
 * order: the first has LCP lcp(x, previous suffix) and the others |x| - 1
 * (x up to its terminator), and their BWT character is the one before x in S,
 * or a terminator for x = S$.  Built in O(|S|^2) for any k.
 */
struct Segment { ulint rows, value; };
void copies_rlbwt(const std::string& S, ulint k, std::vector<uchar>& raw_heads, std::vector<ulint>& lens,
                  std::vector<std::vector<Segment>>& segs) {
    const size_t L = S.size();
    std::vector<size_t> suf(L + 1);
    std::iota(suf.begin(), suf.end(), 0);
    // ropebwt3's order: the terminator, at position L, first, and N after T.
    auto rank = [](char c) { return c == 'A' ? 0 : c == 'C' ? 1 : c == 'G' ? 2 : c == 'T' ? 3 : 4; };
    std::sort(suf.begin(), suf.end(), [&](size_t a, size_t b) {
        while (a < L && b < L && S[a] == S[b]) { ++a; ++b; }
        if (a == L || b == L) return a == L && b != L;
        return rank(S[a]) < rank(S[b]);
    });
    raw_heads.clear();
    lens.clear();
    segs.clear();
    for (size_t x = 0; x <= L; ++x) {
        const size_t j = suf[x];
        ulint top = 0;
        if (x > 0) {
            const size_t p = suf[x - 1];
            while (p + top < L && j + top < L && S[p + top] == S[j + top]) ++top;
        }
        const uchar h = j == 0 ? '$' : head_byte(S[j - 1]);
        if (raw_heads.empty() || raw_heads.back() != h) {
            raw_heads.push_back(h);
            lens.push_back(0);
            segs.emplace_back();
        }
        lens.back() += k;
        segs.back().push_back({1, top});
        if (k > 1) segs.back().push_back({k - 1, static_cast<ulint>(L - j)});
    }
}

/** Retained minima-only pairs of runs given as segments (see lcp_keep_mask): prefix and suffix minima. */
std::vector<RunLcpPairs> pairs_from_segments(const std::vector<std::vector<Segment>>& segs) {
    std::vector<RunLcpPairs> out(segs.size());
    for (size_t i = 0; i < segs.size(); ++i) {
        const ulint top = segs[i][0].value;
        // Interior segments as (first offset, rows, value).
        std::vector<std::array<ulint, 3>> in;
        ulint off = 0;
        for (size_t s = 0; s < segs[i].size(); ++s) {
            ulint rows = segs[i][s].rows, start = off;
            if (s == 0) { ++start; --rows; }
            if (rows > 0) in.push_back({start, rows, segs[i][s].value});
            off += segs[i][s].rows;
        }
        std::vector<std::pair<ulint, ulint>> keep;
        ulint m = top;
        for (const auto& g : in)
            if (g[2] < m) { keep.emplace_back(g[0], g[2]); m = g[2]; }
        m = i + 1 < segs.size() ? segs[i + 1][0].value : std::numeric_limits<ulint>::max();
        for (size_t s = in.size(); s-- > 0;)
            if (in[s][2] < m) { keep.emplace_back(in[s][0] + in[s][1] - 1, in[s][2]); m = in[s][2]; }
        std::sort(keep.begin(), keep.end());
        keep.erase(std::unique(keep.begin(), keep.end()), keep.end());
        out[i].emplace_back(0, top);
        out[i].insert(out[i].end(), keep.begin(), keep.end());
    }
    return out;
}

void test_copies_generator() {
    std::cout << "Testing the k-copies RLBWT generator against suffix sorting" << std::endl;
    std::mt19937 rng(5);
    for (int trial = 0; trial < 60; ++trial) {
        const std::string S = trial % 3 ? rand_dna(rng, 1 + rng() % 12) : rand_dna(rng, 1 + rng() % 12, "ACGTN", 5);
        const ulint k = 1 + rng() % 6;
        const std::vector<std::string> strs(k, S);
        const Coll c = build_collection(strs, true);
        std::vector<uchar> h;
        std::vector<ulint> l;
        std::vector<std::vector<Segment>> segs;
        copies_rlbwt(S, k, h, l, segs);
        assert(h == c.raw_heads && l == c.lens);
        std::vector<std::vector<ulint>> lcps(segs.size());
        for (size_t i = 0; i < segs.size(); ++i)
            for (const auto& s : segs[i]) lcps[i].insert(lcps[i].end(), s.rows, s.value);
        assert(lcps == c.lcps_per_run && "the generator's LCPs must equal suffix sorting's");
        assert(pairs_from_segments(segs) == retained_lcp_pairs(c.lcps_per_run, true) &&
               "segment minima must equal retained_lcp_pairs");
    }
    std::cout << "  PASSED" << std::endl;
}

/**
 * E5: BWTs longer than 2^41 (and 2^42), from k copies of a short string with
 * k near 2^40, so runs are near the 5-byte length limit and offsets, LF
 * pointers and n need more than 40 bits.  Checks the 5-byte length files,
 * the ms build, serialize and load, the tms psi build, serialize and load,
 * and queries against naive matching statistics over S.
 */
void test_large_n(const std::string& dir) {
    std::cout << "Testing BWTs longer than 2^41 (k copies of short strings)" << std::endl;
    const ulint max5 = (ulint{1} << 40) - 1;
    struct Case { std::string S; ulint k; };
    const std::vector<Case> cases = {
        {"ACGT", max5},                         // runs of exactly 2^40 - 1, n = 5 (2^40 - 1)
        {"ACGTTGCAAC", (ulint{1} << 38) + 3},   // n = 11 (2^38 + 3), merged runs
        {"GATTACA", (ulint{1} << 39) - 17},     // n = 8 (2^39 - 17)
        {"CANNAT", (ulint{1} << 39) - 3},       // N in the text, as 0x7f
    };
    std::mt19937 rng(41);
    for (const auto& cs : cases) {
        std::vector<uchar> raw;
        std::vector<ulint> lens;
        std::vector<std::vector<Segment>> segs;
        copies_rlbwt(cs.S, cs.k, raw, lens, segs);
        ulint n = 0, maxlen = 0;
        for (ulint l : lens) { n += l; maxlen = std::max(maxlen, l); }
        assert(maxlen <= max5 && "the case's runs must fit 5-byte lengths");
        assert(n > (ulint{1} << 41));
        std::vector<uchar> heads;
        std::vector<ulint> lens2;
        load_through_files(dir, raw, lens, heads, lens2);
        // Minima through a TLMSM file, whose integers are 64-bit and ULEB128.
        auto pairs = pairs_from_segments(segs);
        const std::string mp = dir + "/.fmd_test_large.min";
        assert(rlbwt_io::write_minima(mp, pairs, n));
        std::vector<RunLcpPairs> read_back;
        std::string err;
        const bool ok = rlbwt_io::read_minima(mp, heads.size(), n, read_back, err);
        std::remove(mp.c_str());
        assert(ok && read_back == pairs);
        auto ms = build_ms(heads, lens2, std::move(read_back));
        assert(ms.domain() == n);
        const std::string ip = dir + "/.fmd_test_large.idx";
        assert(ms_serialize::write_index(ip, ms));
        auto loaded = ms_serialize::read_index(ip);
        std::remove(ip.c_str());
        assert(loaded && loaded->domain() == n && loaded->get_widths() == ms.get_widths());
        auto tms = build_tms_psi(heads, lens2);
        assert(tms.domain() == n);
        // Offsets reach k - 1, so the length and offset columns need over 37 bits.
        const auto& w = ms.get_widths();
        assert(w[0] >= 38 && w[2] >= 38 && "length and offset columns must hold the long runs");
        const std::vector<std::string> strs = {cs.S};
        std::vector<std::string> pats = {cs.S, cs.S + cs.S, cs.S.substr(1) + cs.S + "N" + cs.S, "N", "",
                                         std::string(30, cs.S[0]), "ACGTNACGTACGTTTTGCAACGATTACA"};
        for (int k = 0; k < 40; ++k) {
            std::string P;
            const size_t len = rng() % 50;
            for (size_t j = 0; j < len; ++j)
                P.push_back(rng() % 3 ? cs.S[(j + k) % cs.S.size()] : "ACGTN"[rng() % 5]);
            pats.push_back(P);
        }
        std::vector<std::vector<ulint>> want;
        for (const auto& P : pats) want.push_back(naive_ms(strs, P));
        check_queries(ms, tms, pats, want, "large n");
        check_queries(*loaded, tms, pats, want, "large n, loaded");
        std::cout << "  S=" << cs.S << " k=" << cs.k << " n=" << n << " runs=" << heads.size()
                  << " max_run=" << maxlen << " PASSED" << std::endl;
    }
}

}  // namespace

bool run_all_tests(const std::string& data_dir) {
    test_collections(data_dir);
    std::cout << std::endl;
    test_copies_generator();
    std::cout << std::endl;
    test_large_n(data_dir);
    return true;
}

}  // namespace fmd_test
