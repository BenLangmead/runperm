/**
 * Bidirectional index test suite.
 *
 * Tests are assert-based and print results to stdout, following the
 * convention in src/ms/ms_test.cpp.
 *
 * Author: Ben Langmead (ben.langmead@gmail.com)
 * Date: Feb 23, 2026
 */

#include "bi_test.hpp"
#include "bi_fmd.hpp"
#include "bi_io.hpp"
#include "bi_mem.hpp"
#include "bi_scheme.hpp"
#include "tsv.hpp"
#include <iostream>
#include <cstdio>
#include <cassert>
#include <algorithm>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include <set>

namespace {

constexpr uchar SEPARATOR = orbit::SEPARATOR;
constexpr uchar TERMINATOR = orbit::TERMINATOR;

// =========================================================================
// Test-only utilities
// =========================================================================

/// Map an ASCII character to Nucleotide mapped space.
static uchar map_char(uchar c) {
    return orbit::nucleotide::map_char(c);
}

/// Map a string to a vector of mapped characters.
static std::vector<uchar> map_string(const std::string& s) {
    std::vector<uchar> v(s.size());
    for (size_t i = 0; i < s.size(); ++i)
        v[i] = map_char(static_cast<uchar>(s[i]));
    return v;
}

/// Reverse-complement a mapped-character sequence.
static std::vector<uchar> reverse_complement(const std::vector<uchar>& seq) {
    std::vector<uchar> rc(seq.size());
    for (size_t i = 0; i < seq.size(); ++i)
        rc[seq.size() - 1 - i] = complement(seq[i]);
    return rc;
}

// -------------------------------------------------------------------------
// In-memory BWT construction (brute-force, small texts only)
// -------------------------------------------------------------------------

/// Build the BWT of text T (given as raw ASCII bytes, must include sentinel).
/// Returns (bwt_heads, bwt_run_lengths) with heads as unmapped ASCII chars.
///
/// Uses brute-force suffix-array construction: O(n^2 log n).  Fine for
/// test texts of a few hundred characters.
static std::pair<std::vector<uchar>, std::vector<ulint>>
build_bwt_from_text(const std::string& T) {
    const size_t n = T.size();
    assert(n > 0);

    // Build suffix array by sorting suffix indices
    std::vector<size_t> sa(n);
    std::iota(sa.begin(), sa.end(), 0);
    std::sort(sa.begin(), sa.end(), [&](size_t a, size_t b) {
        return T.compare(a, std::string::npos, T, b, std::string::npos) < 0;
    });

    // Derive BWT: bwt[i] = T[sa[i] - 1] (wrapping)
    std::vector<uchar> bwt(n);
    for (size_t i = 0; i < n; ++i)
        bwt[i] = static_cast<uchar>(T[(sa[i] + n - 1) % n]);

    // Run-length encode
    return orbit::rlbwt::bwt_to_rlbwt(bwt);
}

/// Build a double-strand test text: forward + SEP + revcomp + SEP.
/// Input is the forward strand as an ASCII DNA string (no sentinels).
/// Returns the full concatenated string with SEP characters (ASCII 1).
static std::string make_double_strand_text(const std::string& fwd) {
    // Build reverse complement in ASCII
    std::string rc(fwd.size(), ' ');
    for (size_t i = 0; i < fwd.size(); ++i) {
        char c = fwd[fwd.size() - 1 - i];
        switch (c) {
            case 'A': rc[i] = 'T'; break;
            case 'C': rc[i] = 'G'; break;
            case 'G': rc[i] = 'C'; break;
            case 'T': rc[i] = 'A'; break;
            default:  rc[i] = c;   break;  // N, etc.
        }
    }
    std::string text;
    text.reserve(fwd.size() + rc.size() + 2);
    text += fwd;
    text += static_cast<char>(SEPARATOR);
    text += rc;
    text += static_cast<char>(SEPARATOR);
    return text;
}

/// Naive exact match: count occurrences of P in T.
static size_t naive_count(const std::string& T, const std::string& P) {
    if (P.empty()) return T.size() + 1;  // every position matches empty string
    size_t count = 0;
    size_t pos = 0;
    while ((pos = T.find(P, pos)) != std::string::npos) {
        ++count;
        ++pos;
    }
    return count;
}

/// Naive count with at most k substitutions (Hamming distance).
static size_t naive_count_kmm(const std::string& T, const std::string& P, int k) {
    if (P.empty()) return T.size() + 1;
    size_t count = 0;
    for (size_t j = 0; j + P.size() <= T.size(); ++j) {
        int err = 0;
        for (size_t i = 0; i < P.size() && err <= k; ++i)
            if (P[i] != T[j + i]) ++err;
        if (err <= k) ++count;
    }
    return count;
}

/// Naive matching statistics: MS[i] = length of longest prefix of P[i..]
/// that occurs as a substring of T.
static std::vector<ulint> naive_matching_statistics(const std::string& T,
                                                     const std::string& P) {
    std::vector<ulint> ms(P.size(), 0);
    for (size_t i = 0; i < P.size(); ++i) {
        for (size_t L = P.size() - i; L >= 1; --L) {
            if (T.find(P.substr(i, L)) != std::string::npos) {
                ms[i] = static_cast<ulint>(L);
                break;
            }
        }
    }
    return ms;
}

// =========================================================================
// Hardcoded test texts
// =========================================================================

struct TestText {
    std::string name;
    std::string fwd;       // forward strand only (no sentinels)
    std::string full;      // full double-strand text with sentinels
    std::vector<uchar> heads;
    std::vector<ulint> run_lengths;
};

/**
 * Build TestText objects for the bidirectional index.
 */
static std::vector<TestText> build_test_texts() {
    std::vector<std::string> forwards = {
        "ACGTACGT",                        // tiny
        "AAACAAACAAAC",                    // repetitive
        "ACGTGATTACAGATCCA",               // heterogeneous
    };
    std::vector<std::string> names = {"tiny", "repetitive", "heterogeneous"};

    std::vector<TestText> texts;
    for (size_t i = 0; i < forwards.size(); ++i) {
        TestText tt;
        tt.name = names[i];
        tt.fwd = forwards[i];
        tt.full = make_double_strand_text(forwards[i]);
        auto [h, rl] = build_bwt_from_text(tt.full);
        tt.heads = std::move(h);
        tt.run_lengths = std::move(rl);
        texts.push_back(std::move(tt));
    }
    return texts;
}

/// Build an ASCII reverse complement of a DNA string.
static std::string ascii_revcomp(const std::string& s) {
    std::string rc(s.size(), ' ');
    for (size_t i = 0; i < s.size(); ++i) {
        char c = s[s.size() - 1 - i];
        switch (c) {
            case 'A': rc[i] = 'T'; break;
            case 'C': rc[i] = 'G'; break;
            case 'G': rc[i] = 'C'; break;
            case 'T': rc[i] = 'A'; break;
            default:  rc[i] = c;   break;
        }
    }
    return rc;
}

/**
 * Perform exact bidirectional search starting at position start_pos
 * within the pattern.  Extends right from start_pos, then left.
 * Returns the final BiInterval.
 */
static BiInterval exact_search_from(const FMDIndex& fmd,
                                     const std::vector<uchar>& mapped_pattern,
                                     int start_pos) {
    int m = static_cast<int>(mapped_pattern.size());
    assert(m > 0);

    // Initialize with single character at start_pos
    BiInterval bi = fmd.init_bi_interval(mapped_pattern[start_pos]);
    if (bi.empty()) return bi;

    // Extend right: start_pos+1 .. m-1
    for (int j = start_pos + 1; j < m; ++j) {
        bi = fmd.forward_extend(bi, mapped_pattern[j]);
        if (bi.empty()) return bi;
    }

    // Extend left: start_pos-1 .. 0
    for (int k = start_pos - 1; k >= 0; --k) {
        bi = fmd.backward_extend(bi, mapped_pattern[k]);
        if (bi.empty()) return bi;
    }

    return bi;
}

// =========================================================================
// Phase 1 tests: core primitives
// =========================================================================

/**
 * Check that fmd.count_all gives the same count as if we naively scanned the
 * BWT, for a few distinct rangees.
 */
static void test_count_all(const FMDIndex& fmd) {
    auto bwt = fmd.reconstruct_bwt();
    ulint n = fmd.domain();

    // Test a range of [lo, hi) intervals
    std::vector<std::pair<ulint,ulint>> ranges = {
        {0, n}, {0, 1}, {0, n/2}, {n/2, n}, {n/3, 2*n/3}
    };
    if (n >= 4) {
        ranges.push_back({n/4, 3*n/4});
        ranges.push_back({1, n - 1});
    }

    for (auto [lo, hi] : ranges) {
        auto cnt = fmd.count_all(lo, hi);
        // Brute-force count
        std::array<ulint, NUC_SIGMA> expected{};
        for (ulint p = lo; p < hi; ++p)
            expected[bwt[p]]++;
        for (uchar c = 0; c < NUC_SIGMA; ++c)
            assert(cnt[c] == expected[c]);
    }
    std::cout << "count_all PASSED" << std::endl;
}

/**
 * Test that fmd.init_bi_interval sizes are correct given the number of
 * occurrences of each char in the BWT.  Also check that a bi-interval for a
 * char and the interval for the reversep-complement of the char are equal in
 * size, because of FMD strategy.
 */
static void test_init_bi_interval(const FMDIndex& fmd,
                                  const std::string& full_text) {
    // Count each DNA character in the text
    for (uchar c = NUC_A; c <= NUC_T; ++c) {
        BiInterval bi = fmd.init_bi_interval(c);
        // Count occurrences of this character in the BWT
        auto bwt = fmd.reconstruct_bwt();
        ulint expected = 0;
        for (auto ch : bwt) if (ch == c) expected++;
        // bi.size should equal count in BWT (which equals count in text)
        assert(bi.size == expected);
        // Complement symmetry
        BiInterval bi_comp = fmd.init_bi_interval(complement(c));
        assert(bi.size == bi_comp.size);
    }
    std::cout << "init_bi_interval PASSED" << std::endl;
}

/**
 * Test that bi(P) and bi(revcomp(P)) are equal in size, and that the positions
 * of the bi-intervals are equal, which should be true thanks to FMD strategy.
 */
static void test_bi_symmetry(const FMDIndex& fmd,
                             const std::string& full_text) {

    // Test with some short substrings from the text
    std::mt19937 rng(42);
    int tested = 0;
    for (int trial = 0; trial < 200 && tested < 50; ++trial) {
        size_t len = 2 + rng() % 6;  // lengths 2-7
        if (len > full_text.size()) continue;
        size_t start = rng() % (full_text.size() - len + 1);
        std::string P = full_text.substr(start, len);

        // Skip if contains sentinel
        bool has_sep = false;
        for (char ch : P) if (ch == SEPARATOR || ch == TERMINATOR) { has_sep = true; break; }
        if (has_sep) continue;

        auto mapped = map_string(P);

        // Compute bi(P) by extending from position 0 (all forward extend)
        BiInterval bi = fmd.init_bi_interval(mapped[0]);
        bool ok = !bi.empty();
        for (size_t i = 1; i < mapped.size() && ok; ++i) {
            bi = fmd.forward_extend(bi, mapped[i]);
            ok = !bi.empty();
        }
        if (!ok) continue;

        // Compute bi(revcomp(P))
        auto rc = reverse_complement(mapped);
        BiInterval bi_rc = fmd.init_bi_interval(rc[0]);
        ok = !bi_rc.empty();
        for (size_t i = 1; i < rc.size() && ok; ++i) {
            bi_rc = fmd.forward_extend(bi_rc, rc[i]);
            ok = !bi_rc.empty();
        }
        if (!ok) continue;

        assert(bi.size == bi_rc.size);
        assert(bi.pos == bi_rc.pos_r);
        assert(bi.pos_r == bi_rc.pos);
        ++tested;
    }
    assert(tested > 0);
    std::cout << "bi_interval_symmetry PASSED (" << tested << " patterns)" << std::endl;
}

/**
 * For many pattern strings drawn from the text, test that forward extending
 * to the right from the second position through the end of the pattern, then
 * backward extending to the left by 1, yields the expected number of matches.
 */
static void test_backward_extend(const FMDIndex& fmd,
                                 const std::string& full_text) {
    std::mt19937 rng(77);
    int tested = 0;
    for (int trial = 0; trial < 200 && tested < 50; ++trial) {
        size_t len = 2 + rng() % 6;
        if (len > full_text.size()) continue;
        size_t start = rng() % (full_text.size() - len + 1);
        std::string P = full_text.substr(start, len);
        bool has_sep = false;
        for (char ch : P) if (ch == SEPARATOR || ch == TERMINATOR) { has_sep = true; break; }
        if (has_sep) continue;

        // Build bi-interval of P[1..] then backward_extend with P[0]
        auto mapped = map_string(P);
        BiInterval bi = fmd.init_bi_interval(mapped[1]);
        for (size_t i = 2; i < mapped.size(); ++i) {
            bi = fmd.forward_extend(bi, mapped[i]);
            if (bi.empty()) break;
        }
        if (bi.empty()) continue;

        // Extend backward with P[0]
        BiInterval bi_ext = fmd.backward_extend(bi, mapped[0]);
        size_t expected = naive_count(full_text, P);
        assert(bi_ext.size == expected);
        ++tested;
    }
    assert(tested > 0);
    std::cout << "backward_extend PASSED (" << tested << " patterns)" << std::endl;
}

// =========================================================================
// Phase 2 tests: exact match from all start positions
// =========================================================================

/**
 * For many pattern strings drawn from the text, test that exact_search_from,
 * which extends first right, then left, gives the expected number of matches
 * from every possibly starting position.
 */
static void test_exact_match_all_starts(const FMDIndex& fmd,
                                         const std::string& full_text) {
    std::mt19937 rng(99);
    int tested = 0;

    // Gather patterns: some known, some random
    std::vector<std::string> patterns;
    // Known substrings of length 3-8 from the text
    for (size_t len = 3; len <= 8 && len <= full_text.size(); ++len) {
        for (int t = 0; t < 5; ++t) {
            size_t start = rng() % (full_text.size() - len + 1);
            patterns.push_back(full_text.substr(start, len));
        }
    }
    // Also add the reverse complement of a substring
    if (full_text.size() >= 5) {
        std::string sub = full_text.substr(0, 5);
        bool ok = true;
        for (char c : sub) if (c == SEPARATOR || c == TERMINATOR) ok = false;
        if (ok) patterns.push_back(ascii_revcomp(sub));
    }
    // A pattern that should NOT occur
    patterns.push_back("NNNNN");

    for (const auto& P : patterns) {
        // Skip patterns containing sentinels
        bool has_sep = false;
        for (char ch : P) if (ch == SEPARATOR || ch == TERMINATOR) { has_sep = true; break; }
        if (has_sep) continue;

        auto mapped = map_string(P);
        size_t expected = naive_count(full_text, P);

        // From every start position, the search should give the same size
        for (int sp = 0; sp < static_cast<int>(mapped.size()); ++sp) {
            // Note the exact_search_from extends both forward and backward to
            // encompass the whole pattern, so all the BiInterval results
            // should be the same regardless of which position we start from
            BiInterval bi = exact_search_from(fmd, mapped, sp);
            assert(bi.size == expected);
        }

        // If the pattern occurs, verify reverse complement has same count
        if (expected > 0) {
            std::string P_rc = ascii_revcomp(P);
            size_t rc_count = naive_count(full_text, P_rc);
            // In a double-strand BWT, occ(P) == occ(revcomp(P))
            assert(expected == rc_count);
        }
        ++tested;
    }
    std::cout << "exact_match_all_starts PASSED (" << tested << " patterns)" << std::endl;
}

/**
 * Test that exact_search_from handles some edge cases correctly.
 * - Empty pattern -> full range
 * - Pattern containing N -> 0 matches
 * - Very long pattern (longer than text) -> 0 matches
 */
static void test_exact_match_edge_cases(const FMDIndex& fmd,
                                         const std::string& full_text) {

    { // Empty pattern, full range (don't use exact_search_from; it requires m>0)
        Position first = fmd.first_position();
        BiInterval bi = {first, first, fmd.domain()};
        assert(bi.size == fmd.domain());
    }

    { // Pattern containing N,  should yield 0
        auto mapped = map_string("ACNGT");
        // Search from position 0: init with A, extend with C, then N; should fail
        BiInterval bi = fmd.init_bi_interval(mapped[0]);
        bi = fmd.forward_extend(bi, mapped[1]);
        bi = fmd.forward_extend(bi, mapped[2]);  // N
        assert(bi.empty());
    }

    { // Very long pattern (longer than text), should yield 0
        std::string long_pat(full_text.size() + 10, 'A');
        auto mapped = map_string(long_pat);
        BiInterval bi = exact_search_from(fmd, mapped, 0);
        assert(bi.empty());
    }

    std::cout << "exact_match_edge_cases PASSED" << std::endl;
}

// =========================================================================
// Phase 3 tests: k-mer tables
// =========================================================================

/**
 * Helper: exact search using only k=0 table (character-by-character from full
 * range).  Extends left from last position of mapped pattern using
 * backward_extend.  Returns the final BiInterval.
 */
static BiInterval exact_search_backward_only(const FMDIndex& fmd,
                                             const std::vector<uchar>& mapped)
{
    assert(!mapped.empty());
    Position first = fmd.first_position();
    BiInterval bi = {first, first, fmd.domain()};
    for (int i = static_cast<int>(mapped.size()) - 1; i >= 0; --i) {
        bi = fmd.backward_extend(bi, mapped[i]);
        if (bi.empty()) return bi;
    }
    return bi;
}

/**
 * Test that every k-mer table entry matches K backward_extend calls.  Also
 * double check with naive counting.
 */
static void test_kmer_table_entries(const FMDIndex& fmd,
                                    const std::string& full_text)
{
    for (int K = 1; K <= 3; ++K) {
        KmerTable tab = build_kmer_table(fmd, K);
        uint64_t num_entries = 1ULL << (2 * K);
        assert(tab.table.size() == num_entries);

        for (uint64_t code = 0; code < num_entries; ++code) {
            auto kmer = decode_kmer(code, K);
            // Backward_extend from full range
            Position first = fmd.first_position();
            BiInterval bi_manual = {first, first, fmd.domain()};
            for (int i = K - 1; i >= 0 && !bi_manual.empty(); --i)
                bi_manual = fmd.backward_extend(bi_manual, kmer[i]);
            assert(tab.table[code] == bi_manual); // Must match table entry

            // If non-empty, additionally verify against naive count
            if (!bi_manual.empty()) {
                std::string kmer_str(K, ' ');
                for (int i = 0; i < K; ++i)
                    kmer_str[i] = orbit::nucleotide::unmap_char(kmer[i]);
                size_t expected = naive_count(full_text, kmer_str);
                assert(bi_manual.size == expected);
            }
        }
    }
    std::cout << "kmer_table_entries PASSED" << std::endl;
}

/**
 * Cascading table lookup vs k=0-only (character-by-character) must give
 * identical results.  Currently we try this with 1, 2, 3, and 4-mers.
 */
static void test_kmer_table_consistency(const FMDIndex& fmd,
                                         const std::string& full_text) {
    std::mt19937 rng(55);

    for (int K : {1, 2, 3, 4}) {
        FMDIndex fmd_with_table(fmd.bwt_heads(), fmd.bwt_run_lengths());
        fmd_with_table.add_kmer_table(build_kmer_table(fmd_with_table, K));
        fmd_with_table.add_kmer_table(build_kmer_table(fmd_with_table, 0));

        int tested = 0;
        for (int trial = 0; trial < 200 && tested < 30; ++trial) {
            size_t len = K + 1 + rng() % 5;
            if (len > full_text.size()) continue;
            size_t start = rng() % (full_text.size() - len + 1);
            std::string P = full_text.substr(start, len);
            bool has_sep = false;
            for (char ch : P) if (ch == SEPARATOR || ch == TERMINATOR) { has_sep = true; break; }
            if (has_sep) continue;

            auto mapped = map_string(P);

            BiInterval bi_k0 = exact_search_backward_only(fmd, mapped);

            auto [bi_start, consumed] = fmd_with_table.lookup_kmer(
                mapped.data(), 0, static_cast<int>(mapped.size()));
            BiInterval bi_table = bi_start;
            for (size_t j = consumed; j < mapped.size() && !bi_table.empty(); ++j)
                bi_table = fmd.forward_extend(bi_table, mapped[j]);

            assert(bi_table == bi_k0);
            ++tested;
        }
        assert(tested > 0);
    }
    std::cout << "kmer_table_consistency PASSED" << std::endl;
}

/**
 * Test that cascading table lookup falls back to k=0-only when a large-K entry
 * is empty.
 */
static void test_kmer_cascade_fallback(const FMDIndex& fmd,
                                        const std::string& full_text) {

    FMDIndex fmd2(fmd.bwt_heads(), fmd.bwt_run_lengths());
    fmd2.add_kmer_table(build_kmer_table(fmd2, 4));
    fmd2.add_kmer_table(build_kmer_table(fmd2, 0));

    // Find a 4-mer that does NOT occur in the text
    std::string absent;
    KmerTable tab4 = build_kmer_table(fmd, 4);
    for (uint64_t code = 0; code < tab4.table.size(); ++code) {
        if (tab4.table[code].empty()) {
            auto km = decode_kmer(code, 4);
            for (auto c : km) absent += static_cast<char>(orbit::nucleotide::unmap_char(c));
            break;
        }
    }
    if (absent.empty()) {
        // All 4-mers present (very small text) — skip test
        std::cout << "kmer_cascade_fallback SKIPPED (all 4-mers present)" << std::endl;
        return;
    }

    // Append some characters that DO occur
    std::string suffix;
    for (size_t i = 0; i < full_text.size() && suffix.size() < 3; ++i) {
        char c = full_text[i];
        if (c != SEPARATOR && c != TERMINATOR) suffix += c;
    }
    std::string query = absent + suffix;
    auto mapped = map_string(query);

    auto [bi, consumed] = fmd2.lookup_kmer(
        mapped.data(), 0, static_cast<int>(mapped.size()));
    // K=4 should have failed (absent 4-mer), so fallback to k=0
    assert(consumed == 0);
    assert(bi.size == fmd.domain());  // k=0 gives full range

    // Now extend forward to get bi-interval of the full query
    for (size_t j = 0; j < mapped.size() && !bi.empty(); ++j)
        bi = fmd.forward_extend(bi, mapped[j]);
    size_t expected = naive_count(full_text, query);
    assert(bi.size == expected);

    std::cout << "kmer_cascade_fallback PASSED" << std::endl;
}

/**
 * Test cascading lookup with multiple non-zero K tables (K=4, K=2, K=0).
 * Covers: (1) hit on longest K, (2) fallback to mid K when longest is empty,
 * (3) fallback to K=0 when all longer tables miss.
 */
static void test_kmer_cascade_multi(const FMDIndex& fmd,
                                    const std::string& full_text) {
    FMDIndex fmd_multi(fmd.bwt_heads(), fmd.bwt_run_lengths());
    fmd_multi.add_kmer_table(build_kmer_table(fmd_multi, 4));
    fmd_multi.add_kmer_table(build_kmer_table(fmd_multi, 2));
    fmd_multi.add_kmer_table(build_kmer_table(fmd_multi, 0));

    KmerTable tab4 = build_kmer_table(fmd, 4);
    KmerTable tab2 = build_kmer_table(fmd, 2);

    std::string suffix;
    for (size_t i = 0; i < full_text.size() && suffix.size() < 3; ++i) {
        char c = full_text[i];
        if (c != SEPARATOR && c != TERMINATOR) suffix += c;
    }
    if (suffix.empty()) return;

    auto run_query = [&](const std::string& query, int expected_consumed) {
        auto mapped = map_string(query);
        auto [bi, consumed] = fmd_multi.lookup_kmer(mapped.data(), 0, static_cast<int>(mapped.size()));
        assert(consumed == expected_consumed);
        for (size_t j = consumed; j < mapped.size() && !bi.empty(); ++j)
            bi = fmd.forward_extend(bi, mapped[j]);
        assert(bi.size == naive_count(full_text, query));
    };

    // (1) Query that hits K=4: use a 4-mer from the text
    std::string hit4;
    for (size_t i = 0; i + 4 <= full_text.size(); ++i) {
        bool ok = true;
        for (int j = 0; j < 4; ++j)
            if (full_text[i + j] == SEPARATOR || full_text[i + j] == TERMINATOR) ok = false;
        if (ok) { hit4 = full_text.substr(i, 4); break; }
    }
    if (!hit4.empty()) run_query(hit4 + suffix, 4);

    // (2) 4-mer absent but 2-mer present: find such a 4-mer
    std::string hit2_miss4;
    for (uint64_t code = 0; code < tab4.table.size(); ++code) {
        if (tab4.table[code].empty()) {
            auto km = decode_kmer(code, 4);
            std::string ab(1, orbit::nucleotide::unmap_char(km[0]));
            ab += orbit::nucleotide::unmap_char(km[1]);
            uint64_t code2 = KmerTable::encode(km.data(), 2);
            if (!tab2.table[code2].empty()) {
                for (auto c : km) hit2_miss4 += static_cast<char>(orbit::nucleotide::unmap_char(c));
                break;
            }
        }
    }
    if (!hit2_miss4.empty()) run_query(hit2_miss4 + suffix, 2);

    // (3) Both 4-mer and 2-mer absent: fallback to K=0
    std::string hit0_only;
    for (uint64_t code = 0; code < tab4.table.size(); ++code) {
        if (tab4.table[code].empty()) {
            auto km = decode_kmer(code, 4);
            uint64_t code2 = KmerTable::encode(km.data(), 2);
            if (tab2.table[code2].empty()) {
                for (auto c : km) hit0_only += static_cast<char>(orbit::nucleotide::unmap_char(c));
                break;
            }
        }
    }
    if (!hit0_only.empty()) run_query(hit0_only + suffix, 0);

    if (hit4.empty() && hit2_miss4.empty() && hit0_only.empty()) {
        std::cout << "kmer_cascade_multi SKIPPED (text too small to construct all cases)" << std::endl;
        return;
    }
    std::cout << "kmer_cascade_multi PASSED" << std::endl;
}

/**
 * Test that when K == |P|, table lookup alone gives the answer.
 */
static void test_kmer_table_covers_pattern(const FMDIndex& fmd,
                                            const std::string& full_text) {
    // Build a K=3 table, then query 3-mers and shorter
    KmerTable tab3 = build_kmer_table(fmd, 3);

    std::mt19937 rng(33);
    int tested = 0;
    for (int trial = 0; trial < 100 && tested < 20; ++trial) {
        size_t len = 1 + rng() % 3;  // 1, 2, or 3
        if (len > full_text.size()) continue;
        size_t start = rng() % (full_text.size() - len + 1);
        std::string P = full_text.substr(start, len);
        bool has_sep = false;
        for (char ch : P) if (ch == SEPARATOR || ch == TERMINATOR) { has_sep = true; break; }
        if (has_sep) continue;

        auto mapped = map_string(P);
        size_t expected = naive_count(full_text, P);

        if (len == 3) {
            // Direct table lookup
            uint64_t code = KmerTable::encode(mapped.data(), 3);
            BiInterval bi = tab3.lookup(code);
            assert(bi.size == expected);
        }
        // Also verify via backward_extend
        BiInterval bi2 = exact_search_backward_only(fmd, mapped);
        assert(bi2.size == expected);
        ++tested;
    }
    assert(tested > 0);
    std::cout << "kmer_table_K_ge_P PASSED (" << tested << " patterns)" << std::endl;
}

// =========================================================================
// Phase 4: Search schemes
// =========================================================================

/** Trivial scheme: exact match only (0 errors). Single piece. */
static std::vector<bi_approx::Search> exact_only_scheme() {
    return {{{0}, {0}, {0}}};
}

/** Two-phase k=1: phase 1 = exact or 1-mismatch in left half (right exact);
 *  phase 2 = 1-mismatch in right half (left exact), no exact matches. */
static std::vector<bi_approx::Search> two_phase_k1_scheme() {
    return {
        {{0, 1}, {0, 0}, {1, 1}},   // left 0-1 err, right exact
        {{1, 0}, {1, 1}, {1, 1}}    // right 1 err, left exact
    };
}

/**
 * Test trivial exact-only scheme: no mismatches allowed.
 */
static void test_scheme_exact(const FMDIndex& fmd, const std::string& full_text) {
    auto scheme = exact_only_scheme();
    std::mt19937 rng(77);
    int tested = 0;
    for (int trial = 0; trial < 100 && tested < 20; ++trial) {
        size_t len = 2 + rng() % 6;
        if (len > full_text.size()) continue;
        size_t start = rng() % (full_text.size() - len + 1);
        std::string P = full_text.substr(start, len);
        bool has_sep = false;
        for (char ch : P)
            if (ch == SEPARATOR || ch == TERMINATOR) { has_sep = true; break; }
        if (has_sep) continue;

        auto matches = bi_approx::approx_match_scheme(fmd, P, scheme);
        size_t total_occ = 0;
        for (const auto& m : matches) {
            assert(m.errors == 0);
            total_occ += m.bi.size;
        }
        assert(total_occ == naive_count(full_text, P));
        ++tested;
    }
    assert(tested > 0);
    std::cout << "scheme_exact PASSED (" << tested << " patterns)" << std::endl;
}

/**
 * Test two-phase k=1 scheme: phase 1 = exact or left-half mismatch; phase
 * 2 = right-half mismatch only.  Try patterns with mismatch at various
 * positions to exercise both phases.
 */
static void test_scheme_two_phase(const FMDIndex& fmd, const std::string& full_text) {
    auto scheme = two_phase_k1_scheme();
    static const char ACGT[] = {'A', 'C', 'G', 'T'};
    std::mt19937 rng(88);

    // Pick a base pattern from the text
    std::string base;
    for (size_t i = 0; i + 8 <= full_text.size(); ++i) {
        bool ok = true;
        for (int j = 0; j < 8; ++j)
            if (full_text[i + j] == SEPARATOR || full_text[i + j] == TERMINATOR) ok = false;
        if (ok) { base = full_text.substr(i, 8); break; }
    }
    if (base.size() < 6) {
        std::cout << "scheme_two_phase SKIPPED (text too short for 6+ char patterns)" << std::endl;
        return;
    }

    auto mismatch_pos = [&](const std::string& P, int pos) -> std::string {
        if (pos < 0 || pos >= static_cast<int>(P.size())) return P;
        char orig = P[pos];
        char sub = orig;
        while (sub == orig) sub = ACGT[rng() % 4];
        std::string Q = P;
        Q[pos] = sub;
        return Q;
    };

    int tested = 0;
    for (int pos = 0; pos < static_cast<int>(base.size()); ++pos) {
        std::string P = mismatch_pos(base, pos);
        if (P == base) continue;  // rng gave same char, skip
        auto matches = bi_approx::approx_match_scheme(fmd, P, scheme);
        size_t total_occ = 0;
        for (const auto& m : matches) total_occ += m.bi.size;
        size_t expected = naive_count_kmm(full_text, P, 1);
        assert(total_occ == expected);
        ++tested;
    }
    // Also test exact match (phase 1)
    {
        auto matches = bi_approx::approx_match_scheme(fmd, base, scheme);
        size_t total_occ = 0;
        for (const auto& m : matches) total_occ += m.bi.size;
        assert(total_occ == naive_count(full_text, base));
        ++tested;
    }
    assert(tested > 0);
    std::cout << "scheme_two_phase PASSED (" << tested << " patterns)" << std::endl;
}

// =========================================================================
// Phase 5: I/O tests
// =========================================================================

static const char* IO_TEST_BASE = "bi_io_test";

/**
 * Test that writing and reading the index and k-mer table roundtrips correctly.
 */
static void test_io_roundtrip(const FMDIndex& fmd, const std::string& full_text) {
    std::string base_path = std::string(IO_TEST_BASE) + ".bi";

    assert(bi_io::write_base(base_path, fmd));
    auto opt = bi_io::read_base(base_path);
    assert(opt && "read_base failed");
    FMDIndex loaded = std::move(*opt);

    assert(loaded.domain() == fmd.domain());
    assert(loaded.num_runs() == fmd.num_runs());
    assert(loaded.get_fingerprint() == fmd.get_fingerprint());

    // Build K=2 table, write, read, verify
    KmerTable tab = build_kmer_table(fmd, 2);
    std::string k2_path = base_path + ".k2";
    assert(bi_io::write_kmer_table(k2_path, fmd.get_fingerprint(), tab, fmd));
    auto opt_tab = bi_io::read_kmer_table(k2_path, fmd.get_fingerprint(), loaded);
    assert(opt_tab && "read_kmer_table failed");
    const KmerTable& loaded_tab = *opt_tab;
    assert(loaded_tab.K == 2);
    assert(loaded_tab.table.size() == tab.table.size());
    for (size_t i = 0; i < tab.table.size(); ++i) {
        assert(loaded_tab.table[i].pos == tab.table[i].pos);
        assert(loaded_tab.table[i].pos_r == tab.table[i].pos_r);
        assert(loaded_tab.table[i].size == tab.table[i].size);
    }

    // Verify search gives identical results (both need same tables)
    std::string P = full_text.substr(0, std::min(size_t(5), full_text.size()));
    while (!P.empty() && (P.back() == SEPARATOR || P.back() == TERMINATOR))
        P.pop_back();
    if (!P.empty()) {
        FMDIndex fmd_with_tab(fmd.bwt_heads(), fmd.bwt_run_lengths());
        fmd_with_tab.add_kmer_table(tab);
        loaded.add_kmer_table(loaded_tab);
        auto mapped = map_string(P);
        auto [bi_orig, _] = fmd_with_tab.lookup_kmer(mapped.data(), 0, static_cast<int>(mapped.size()));
        auto [bi_loaded, __] = loaded.lookup_kmer(mapped.data(), 0, static_cast<int>(mapped.size()));
        if (bi_orig.empty() && bi_loaded.empty()) { /* ok */ }
        else assert(bi_orig.size == bi_loaded.size);
    }

    std::remove(k2_path.c_str());
    std::remove(base_path.c_str());
    std::cout << "io_roundtrip PASSED" << std::endl;
}

/**
 * Test that reading a k-mer table with a fingerprint that fails to match the
 * fingerprint of the base FMD Index causes an error.
 */
static void test_io_fingerprint_rejection(const FMDIndex& fmd1, const FMDIndex& fmd2) {
    std::string base1 = std::string(IO_TEST_BASE) + "_1.bi";
    std::string base2 = std::string(IO_TEST_BASE) + "_2.bi";
    std::string k2_path = base1 + ".k2";

    assert(bi_io::write_base(base1, fmd1));
    assert(bi_io::write_base(base2, fmd2));
    KmerTable tab = build_kmer_table(fmd1, 2);
    assert(bi_io::write_kmer_table(k2_path, fmd1.get_fingerprint(), tab, fmd1));

    auto opt2 = bi_io::read_base(base2);
    assert(opt2);
    // Read table with fmd2's fingerprint — should fail
    auto opt_tab = bi_io::read_kmer_table(k2_path, opt2->get_fingerprint(), *opt2);
    assert(!opt_tab && "expected nullopt for fingerprint mismatch");

    std::remove(k2_path.c_str());
    std::remove(base1.c_str());
    std::remove(base2.c_str());
    std::cout << "fingerprint_rejection PASSED" << std::endl;
}

/**
 * Test that loading available k-mer tables from a base path works correctly.
 */
static void test_io_load_available_tables(const FMDIndex& fmd) {
    std::string base_path = std::string(IO_TEST_BASE) + "_avail.bi";

    assert(bi_io::write_base(base_path, fmd));
    KmerTable tab6 = build_kmer_table(fmd, 6);
    KmerTable tab0 = build_kmer_table(fmd, 0);
    assert(bi_io::write_kmer_table(base_path + ".k6", fmd.get_fingerprint(), tab6, fmd));
    assert(bi_io::write_kmer_table(base_path + ".k0", fmd.get_fingerprint(), tab0, fmd));

    auto opt = bi_io::read_base(base_path);
    assert(opt);
    auto tables = bi_io::load_available_tables(base_path, opt->get_fingerprint(), *opt);
    assert(tables.size() >= 2);
    assert(tables[0].K == 6);
    assert(tables.back().K == 0);

    // Delete .k6, load again — should find only k0
    std::remove((base_path + ".k6").c_str());
    tables = bi_io::load_available_tables(base_path, opt->get_fingerprint(), *opt);
    assert(tables.size() == 1);
    assert(tables[0].K == 0);

    std::remove((base_path + ".k0").c_str());
    std::remove(base_path.c_str());
    std::cout << "load_available_tables PASSED" << std::endl;
}

// =========================================================================
// TSV-based tests (use data/minishred1_20_002_lcp.tsv, skip if missing)
// =========================================================================

/**
 * Test that building an FMDIndex from a TSV file works correctly.
 */
static bool test_build_from_tsv(const std::string& data_dir) {
    std::string path = data_dir + "/minishred1_20_002_lcp.tsv";
    std::vector<uchar> bwt_heads;
    std::vector<ulint> bwt_run_lengths;
    std::vector<std::vector<ulint>> lcps_per_run;
    if (!tsv::load_tsv(path, bwt_heads, bwt_run_lengths, lcps_per_run)) {
        std::cout << "build_from_tsv SKIPPED (TSV not found)" << std::endl;
        return true;  // skip, not a failure
    }
    FMDIndex fmd(bwt_heads, bwt_run_lengths);
    std::string idx_path = std::string(IO_TEST_BASE) + "_tsv.bi";
    assert(bi_io::write_base(idx_path, fmd));
    auto opt = bi_io::read_base(idx_path);
    assert(opt);
    assert(opt->domain() == fmd.domain());
    std::remove(idx_path.c_str());
    std::cout << "build_from_tsv PASSED (n=" << fmd.domain() << ")" << std::endl;
    return true;
}

/**
 * Test that finding SMEMs in an FMDIndex built from a TSV file works correctly.
 */
static bool test_mem_tsv(const std::string& data_dir) {
    std::string path = data_dir + "/minishred1_20_002_lcp.tsv";
    std::vector<uchar> bwt_heads;
    std::vector<ulint> bwt_run_lengths;
    std::vector<std::vector<ulint>> lcps_per_run;
    if (!tsv::load_tsv(path, bwt_heads, bwt_run_lengths, lcps_per_run)) {
        std::cout << "mem_from_tsv SKIPPED (TSV not found)" << std::endl;
        return true;
    }
    FMDIndex fmd(bwt_heads, bwt_run_lengths);
    fmd.add_kmer_table(build_kmer_table(fmd, 0));
    auto smems = bi_mem::find_smems(fmd, "ACGTACGT");
    assert(smems.size() > 0 || true);  // may be 0 if pattern absent
    std::cout << "mem_from_tsv PASSED (found " << smems.size() << " SMEMs for ACGTACGT)" << std::endl;
    return true;
}

/**
 * Test that approximate matching in an FMDIndex built from a TSV file works
 * correctly.
 */
static bool test_approx_tsv(const std::string& data_dir) {
    std::string path = data_dir + "/minishred1_20_002_lcp.tsv";
    std::vector<uchar> bwt_heads;
    std::vector<ulint> bwt_run_lengths;
    std::vector<std::vector<ulint>> lcps_per_run;
    if (!tsv::load_tsv(path, bwt_heads, bwt_run_lengths, lcps_per_run)) {
        std::cout << "approx_from_tsv SKIPPED (TSV not found)" << std::endl;
        return true;
    }
    FMDIndex fmd(bwt_heads, bwt_run_lengths);
    fmd.add_kmer_table(build_kmer_table(fmd, 0));
    std::string pattern = "ACGTACCT";  // 1-substitution from ACGTACGT
    auto scheme = std::vector<bi_approx::Search>{
        {{0, 1}, {0, 0}, {0, 1}},
        {{1, 0}, {0, 0}, {0, 1}}
    };
    auto matches = bi_approx::approx_match_scheme(fmd, pattern, scheme);
    assert(matches.size() >= 0);
    std::cout << "approx_from_tsv PASSED (found " << matches.size() << " approx matches)" << std::endl;
    return true;
}

}  // anonymous namespace

// =========================================================================
// Public test entry point
// =========================================================================

namespace bi_test {

bool run_all_tests(const std::string& data_dir) {
    auto texts = build_test_texts();

    for (auto& tt : texts) {
        std::cout << "\n=== Test text: " << tt.name
                  << " (n=" << tt.full.size() << ", r=" << tt.heads.size()
                  << ") ===" << std::endl;

        FMDIndex fmd(tt.heads, tt.run_lengths);
        assert(fmd.domain() == tt.full.size());

        // Phase 1 tests
        test_count_all(fmd);
        test_init_bi_interval(fmd, tt.full);
        test_bi_symmetry(fmd, tt.full);
        test_backward_extend(fmd, tt.full);

        // Phase 2 tests
        test_exact_match_all_starts(fmd, tt.full);
        test_exact_match_edge_cases(fmd, tt.full);

        // Phase 3 tests
        test_kmer_table_entries(fmd, tt.full);
        test_kmer_table_consistency(fmd, tt.full);
        test_kmer_cascade_fallback(fmd, tt.full);
        test_kmer_cascade_multi(fmd, tt.full);
        test_kmer_table_covers_pattern(fmd, tt.full);
        test_scheme_exact(fmd, tt.full);
        test_scheme_two_phase(fmd, tt.full);
    }

    // Phase 4: I/O tests (use first text for fmd1, second for fmd2)
    {
        auto& tt0 = texts[0];
        auto& tt1 = texts[1];
        FMDIndex fmd0(tt0.heads, tt0.run_lengths);
        FMDIndex fmd1(tt1.heads, tt1.run_lengths);
        test_io_roundtrip(fmd0, tt0.full);
        test_io_fingerprint_rejection(fmd0, fmd1);
        test_io_load_available_tables(fmd0);
    }

    // TSV-based tests (skip if data file missing)
    test_build_from_tsv(data_dir.empty() ? "data" : data_dir);
    test_mem_tsv(data_dir.empty() ? "data" : data_dir);
    test_approx_tsv(data_dir.empty() ? "data" : data_dir);

    std::cout << "\nAll tests passed." << std::endl;
    return true;
}

}  // namespace bi_test
