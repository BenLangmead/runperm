/**
 * MS internal test suite functions.
 *
 * Author: Ben Langmead (ben.langmead@gmail.com)
 * Date: Feb 17, 2026
 */

#include "ms_test.hpp"
#include "ms_rlbwt.hpp"
#include "ms_io.hpp"
#include "tsv.hpp"
#include <iostream>
#include <cassert>
#include <cstdio>
#include <functional>
#include <random>
#include <string>
#include <utility>

namespace {

/**
 * Decode a ULEB128-encoded spillover vector (like those created by
 * build_spill_data) back into a vector of ulints. All sentinel/count/interior
 * values are preserved; the output is simply the list of all decoded integers
 * in the order in which they appear in the buffer.
 */
/** Reconstruct the original text by LF-walking from first() for n steps. */
static std::string reconstruct_text(MSIndexSpillLCP<false>& idx) {
    std::string t;
    t.reserve(static_cast<size_t>(idx.domain()));
    auto pos = idx.first();
    for (ulint i = 0; i < idx.domain(); ++i) {
        t += static_cast<char>(idx.get_character(pos.interval));
        pos = idx.LF(pos);
    }
    std::reverse(t.begin(), t.end());
    return t;
}

/**
 * Naive and slow -- O(|T| * |P|^2) -- matching statistics, for testing
 */
static std::vector<ulint> naive_matching_statistics(const std::string& T, const std::string& P) {
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

/**
 * Decode a ULEB128-encoded spillover vector back into a vector of ulints, for
 * easier testing.
 */
static std::vector<ulint>
decode_spillover_vector(const std::vector<uchar>& spill) {
    std::vector<ulint> result;
    size_t pos = 0;
    while (pos < spill.size()) {
        auto [value, new_pos] = decode_uleb128(spill, pos);
        result.push_back(value);
        pos = new_pos;
    }
    return result;
}

/**
 * Very simple test to see that we can store and retrieve TOP_LCP values.
 */
template <bool SP>
void test_top_lcp() {
    std::cout << "Testing TOP_LCP storage (StoreAbsolutePositions=" << SP << ")" << std::endl;
    const size_t r = 9;
    std::vector<uchar> bwt_heads =       {'T','C','G','A','T', 1 ,'A','T','A'};
    std::vector<ulint> bwt_run_lengths = { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };
    std::vector<std::array<ulint, 1>> run_data(r);
    // Make fake top-LCP information
    for (size_t i = 0; i < r; ++i) run_data[i][0] = i;
    MSIndexTopLCP<SP> index(bwt_heads, bwt_run_lengths, run_data);
    // Make sure we can read out the same fake top-LCP info
    for (ulint i = 0; i < index.intervals(); ++i)
        assert(index.template get<LCPRunCols::TOP_LCP>(i) == i);
    std::cout << "  PASSED" << std::endl;
}

/**
 * Test a case where we have spillover data and we still want to compute the
 * boundary and row-wise min LCPs correctly.
 */
template <bool SP>
void test_spill_lcp_queries() {
    std::cout << "Testing spill LCP (StoreAbsolutePositions=" << SP << ")" << std::endl;
    std::vector<uchar> bwt_heads =       {'T','C','G','A','T', 1 ,'A','T','A'};
    std::vector<ulint> bwt_run_lengths = { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };
    std::vector<std::vector<ulint>> lcps_per_run = {
        {0, 1, 2, 2, 1},      // [1] 0,0,-,-,-,-   min_sub=0  no spill  skinny
        {1, 0, 0},            // [0] 1,0,0,0       min_sub=1  no spill  jumbo
        {0, 1, 0},            // [1] 0,0,-,0       min_sub=0  no spill  skinny
        {1, 0, 1},            // [0] 1,0,0,-       min_sub=1  no spill  jumbo
        {0},                  // [0] 0,0           min_sub=0  no spill  skinny
        {0},                  // [0] 0,0           min_sub=0  no spill  skinny
        {0},                  // [1] 0,0           min_sub=0  no spill  skinny
        {1, 0, 0, 0},         // [0] 1,0,-,-       min_sub=1  no spill  jumbo
        {0, 1, 2, 1, 0, 0},   // [Inf] 0,0,-,-,-,- min_sub=0  no spill  skinny
    };
    auto [run_data, spill_vectors, max_top, max_sub, skinny_count, jumbo_count] = build_spill_data(lcps_per_run);
    std::vector<ulint> decoded_spillover = decode_spillover_vector(spill_vectors[0]);
    assert(skinny_count == 6);
    assert(jumbo_count == 3);
    assert((decoded_spillover == std::vector<ulint>{0, 1, 0, 2, 1, 0, 2, 0, 1, 2, 0, 1, 0, 1, 1, 0, 1, 0, 3, 1, 0, 2, 0, 3, 0}));
    MSIndexSpillLCP<SP> index(bwt_heads, bwt_run_lengths, run_data, std::move(spill_vectors), max_top, max_sub);
    // Make sure we can recovery the boundary LCPs correctly
    for (ulint i = 0; i < index.move_runs(); ++i) {
        assert(boundary_lcp(index, i) == lcps_per_run[i][0]);
        assert(row_min_lcp(index, i) == 0);
    }
    std::cout << "  PASSED" << std::endl;
}

/**
 * Test that the spillover vector is correct, that skinny/jumbo designations
 * are also correct.
 */
template <bool SP>
void test_spillover() {
    std::cout << "Testing spillover correctness (StoreAbsolutePositions=" << SP << ")" << std::endl;
    std::vector<uchar> bwt_heads =       {'T','C','G','A','T', 1 ,'A'};
    std::vector<ulint> bwt_run_lengths = { 4 , 3 , 2 , 3 , 2 , 2 , 3 };
    std::vector<std::vector<ulint>> lcps_per_run = {
        {2, 1, 3, 2}, // [1] 2,1,-,-   min_sub=1  skinny
        {1, 0, 1},    // [0] 1,0,-     min_sub=1  skinny
        {0, 2},       // [2] 0,0,-     min_sub=0  skinny
        {2, 1, 0},    // [1] 2,0,1,0   min_sub=2  skinny
        {3, 0},       // [0] 3,0,-     min_sub=3  jumbo
        {0, 1},       // [1] 0,0,-     min_sub=0  skinny
        {1, 0, 0},    // [inf] 1,0,-,- min_sub=1  skinny
    };
    auto [run_data, spill_vectors, max_top, max_sub, skinny_count, jumbo_count] = build_spill_data(lcps_per_run);
    assert(skinny_count == 6);
    assert(jumbo_count == 1);
    MSIndexSpillLCP<SP> index(bwt_heads, bwt_run_lengths, run_data, std::move(spill_vectors), max_top, max_sub);
    assert(boundary_lcp(index, 0) == 2);
    assert(row_min_lcp(index, 0) == 1);
    assert(range_min_to_top(index, 0, 2) == 1);
    assert(boundary_lcp(index, 1) == 1);
    assert(row_min_lcp(index, 1) == 0);
    std::cout << "  PASSED" << std::endl;
}

/**
 * Test ms_query consistency between lengths- and starts-based indexes.
 */
void test_ms_query_consistency() {
    std::cout << "Testing ms_query consistency (lengths vs starts)" << std::endl;
    std::vector<uchar> bwt_heads =       {'T','C','G','A','T', 1 ,'A'};
    std::vector<ulint> bwt_run_lengths = { 4 , 3 , 2 , 3 , 2 , 2 , 3 };
    std::vector<std::vector<ulint>> lcps_per_run = {
        {2, 1, 3, 2}, {1, 0, 1}, {0, 2}, {2, 1, 0}, {1, 0}, {0, 1}, {1, 0, 0},
    };
    [[maybe_unused]] auto [run_data, spill_vectors, max_top, max_sub, skinny_count, jumbo_count] = build_spill_data(lcps_per_run);
    [[maybe_unused]] auto [run_data2, spill_vectors2, max_top2, max_sub2, skinny2, jumbo2] = build_spill_data(lcps_per_run);
    MSIndexSpillLCP<false> idx_len(bwt_heads, bwt_run_lengths, run_data, std::move(spill_vectors), max_top, max_sub);
    MSIndexSpillLCP<true> idx_starts(bwt_heads, bwt_run_lengths, run_data2, std::move(spill_vectors2), max_top2, max_sub2);
    auto ms_len = ms_query(idx_len, "GAT");
    auto ms_starts = ms_query(idx_starts, "GAT");
    for (size_t i = 0; i < ms_len.size(); ++i) assert(ms_len[i] == ms_starts[i]);
    std::cout << "  PASSED" << std::endl;
}

/**
 * Test coalesce-spillover build from TSV.
 * TODO: This current just checks that there's no assert!  We should also check
 * that the coalescing is correct.
 */
bool test_coalesce_spillover_from_tsv(const std::string& data_dir) {
    std::cout << "Testing coalesce-spillover build from TSV" << std::endl;
    std::string path = data_dir + "/minishred1_20_002_lcp.tsv";
    std::vector<uchar> bwt_heads;
    std::vector<ulint> bwt_run_lengths;
    std::vector<std::vector<ulint>> lcps_per_run;
    if (!tsv::load_tsv(path, bwt_heads, bwt_run_lengths, lcps_per_run)) {
        std::cout << "  DID NOT RUN" << std::endl;
        return false;
    }
    auto [run_data, spill_vectors, max_top, max_sub, skinny_count, jumbo_count] =
        build_spill_data(lcps_per_run, 0.98, true /* coalesce_spillover */);
    MSIndexSpillLCP<false> idx(bwt_heads, bwt_run_lengths, run_data, std::move(spill_vectors), max_top, max_sub);
    auto ms = ms_query(idx, "GATTACA");
    assert(ms.size() == 7);
    std::cout << "  PASSED" << std::endl;
    return true;
}

void test_build_spill_data_coalesce() {
    std::cout << "Testing build_spill_data with coalesce_spillover" << std::endl;
    std::vector<std::vector<ulint>> lcps_per_run = {
        {0, 1, 2, 2, 1}, {1, 0, 0}, {0, 1, 0}, {1, 0, 1}, {0}, {0}, {0},
        {1, 0, 0, 0}, {0, 1, 2, 1, 0, 0},
    };
    auto [run_data, spill_vectors, max_top, max_sub, skinny_count, jumbo_count] =
        build_spill_data(lcps_per_run, 0.98, true);
    assert(skinny_count + jumbo_count == lcps_per_run.size());
    std::cout << "  PASSED" << std::endl;
}

void test_naive_matching_statistics_known() {
    std::cout << "Testing naive_matching_statistics (known answer)" << std::endl;
    const std::string T = "abracadabradad";
    const std::string P = "caradadabrd";
    const std::vector<ulint> expected = {2, 1, 5, 4, 3, 5, 4, 3, 2, 1, 1};
    auto ms = naive_matching_statistics(T, P);
    assert(ms.size() == expected.size());
    for (size_t i = 0; i < ms.size(); ++i)
        assert(ms[i] == expected[i]);
    std::cout << "  PASSED" << std::endl;
}

/**
 * Trace the reposition path (up and down) and print LCP values at each run.
 * Helps diagnose why base under-reports vs coal_split.
 */
static void trace_reposition_path(MSIndexSpillLCP<false>& idx, ulint interval, ulint offset,
                                  uchar c, std::ostream& out) {
    using Position = typename MSIndexSpillLCP<false>::Position;
    Position cur{interval, offset};
    const Position first_pos = idx.first();
    const Position last_pos = idx.last();
    const ulint max_walk = std::min(idx.domain(), static_cast<ulint>(50));  /* limit for output */

    out << "  Reposition from (int=" << interval << " off=" << offset << ") target='" << (char)c << "'" << std::endl;

    out << "  UP path:" << std::endl;
    ulint min_up = LCP_GAP;
    Position pos_up = cur;
    for (ulint k = 0; k < max_walk; ++k) {
        uchar ch = idx.get_character(pos_up.interval);
        ulint lcp_val = (k == 0) ? range_min(idx, pos_up.interval, pos_up.offset, true)
                                : row_min_lcp(idx, pos_up.interval);
        if (lcp_val != LCP_GAP) min_up = std::min(min_up, lcp_val);
        out << "    k=" << k << " int=" << pos_up.interval << " off=" << pos_up.offset
            << " c='" << (char)ch << "' lcp=" << (lcp_val == LCP_GAP ? ulint(-1) : lcp_val)
            << " min_up=" << (min_up == LCP_GAP ? ulint(-1) : min_up) << std::endl;
        if (ch == c) {
            ulint lcp_arrive = range_min(idx, pos_up.interval, pos_up.offset, false);
            out << "    -> arrived at target, lcp_arrive=" << lcp_arrive << std::endl;
            break;
        }
        pos_up = idx.up(pos_up);
        if (pos_up.interval == last_pos.interval && pos_up.offset == last_pos.offset) break;
    }

    out << "  DOWN path:" << std::endl;
    ulint min_down = LCP_GAP;
    Position pos_down = cur;
    for (ulint k = 0; k < max_walk; ++k) {
        uchar ch = idx.get_character(pos_down.interval);
        ulint lcp_val = (k == 0) ? range_min(idx, pos_down.interval, pos_down.offset, false)
                                : row_min_lcp(idx, pos_down.interval);
        if (lcp_val != LCP_GAP) min_down = std::min(min_down, lcp_val);
        out << "    k=" << k << " int=" << pos_down.interval << " off=" << pos_down.offset
            << " c='" << (char)ch << "' lcp=" << (lcp_val == LCP_GAP ? ulint(-1) : lcp_val)
            << " min_down=" << (min_down == LCP_GAP ? ulint(-1) : min_down);
        if (ch == c) {
            ulint lcp_arrive = range_min(idx, pos_down.interval, pos_down.offset, true);
            out << " -> arrived at target, lcp_arrive=" << lcp_arrive;
            if (lcp_arrive != LCP_GAP) min_down = std::min(min_down, lcp_arrive);
            break;
        }
        out << std::endl;
        Position next = idx.down(pos_down);
        if (next.interval == first_pos.interval && next.offset == first_pos.offset) break;
        ulint lcp_boundary = boundary_lcp(idx, next.interval);
        if (lcp_boundary != LCP_GAP) min_down = std::min(min_down, lcp_boundary);
        out << "      boundary at next int=" << next.interval << " = " << lcp_boundary << std::endl;
        pos_down = next;
    }

    out << "  Result: min_up=" << (min_up == LCP_GAP ? ulint(-1) : min_up)
        << " min_down=" << (min_down == LCP_GAP ? ulint(-1) : min_down) << std::endl;
}

/**
 * Step-by-step trace of ms_query for investigation. Prints each backward step,
 * distinguishing LF vs reposition, and when repositioning shows the cap (match_len).
 */
static void trace_ms_query_steps(MSIndexSpillLCP<false>& idx, const std::string& P, std::ostream& out,
                                 int trace_reposition_at = -1) {
    using Position = typename MSIndexSpillLCP<false>::Position;
    Position pos = idx.first();
    ulint match_len = 0;
    int step = 0;
    for (size_t i = P.size(); i > 0; --i) {
        const size_t out_idx = i - 1;
        uchar c = static_cast<uchar>(P[out_idx]);
        if (idx.get_character(pos.interval) == c) {
            pos = idx.LF(pos);
            match_len++;
            out << "  i=" << out_idx << " P[i]='" << (char)c << "': LF  (int=" << pos.interval
                << " off=" << pos.offset << ") match_len=" << match_len << std::endl;
            ++step;
            continue;
        }
        if (trace_reposition_at >= 0 && step == trace_reposition_at) {
            trace_reposition_path(idx, pos.interval, static_cast<ulint>(pos.offset), c, out);
        }
        auto opt = reposition_with_lcp(idx, pos.interval, static_cast<ulint>(pos.offset), c);
        if (!opt) {
            out << "  i=" << out_idx << " P[i]='" << (char)c << "': reposition FAILED" << std::endl;
            break;
        }
        ulint from_int = pos.interval, from_off = pos.offset;
        ulint cap = opt->second;
        match_len = std::min(match_len, cap) + 1;
        pos = opt->first;
        out << "  i=" << out_idx << " P[i]='" << (char)c << "': REPOSITION from (int=" << from_int
            << " off=" << from_off << ") cap=" << cap << " -> (int=" << pos.interval << " off="
            << pos.offset << ") match_len=" << match_len << std::endl;
        ++step;
    }
}

/**
 * Single-string example where ms_query and naive_matching_statistics differ
 * due to LCP compression: the index under-reports. Uses minishred1_20_002_lcp.tsv.
 */
bool test_ms_query_vs_naive_single_example(const std::string& data_dir) {
    std::cout << "Testing ms_query vs naive (single example from TSV)" << std::endl;
    std::string path = data_dir + "/minishred1_20_002_lcp.tsv";
    auto idx = ms_io::build_ms_index_spill_from_tsv<false>(path);
    if (!idx) {
        std::cout << "  DID NOT RUN" << std::endl;
        return false;
    }
    std::string T = reconstruct_text(*idx);
    const std::string P = "GATTGGCTCTCCCTACTCCT";
    assert(T.find(P) != std::string::npos);
    auto ms_idx = ms_query(*idx, P);
    auto ms_naive = naive_matching_statistics(T, P);
    assert(ms_idx == ms_naive && "index and naive MS must match");
    std::cout << "  PASSED" << std::endl;
    return true;
}

/**
 * Test ms_query vs naive matching statistics on a 100 randomly extracted
 * patterns that match exactly.
 */
 bool test_ms_query_vs_naive_tsv(const std::string& data_dir) {
    std::cout << "Testing ms_query vs naive (TSV, extract_errory_string, no reposition)" << std::endl;
    std::string path = data_dir + "/minishred1_20_002_lcp.tsv";
    auto idx = ms_io::build_ms_index_spill_from_tsv<false>(path);
    if (!idx) {
        std::cout << "  DID NOT RUN" << std::endl;
        return false;
    }
    std::string T = reconstruct_text(*idx);
    std::mt19937 rng(42);
    const int num_strings = 100;
    const ulint pattern_len = 100;
    for (int k = 0; k < num_strings; ++k) {
        std::string P = extract_errory_string(*idx, 0.0f, pattern_len, rng);
        assert(T.find(P) != std::string::npos && "extracted string must appear in text");
        auto ms_idx = ms_query(*idx, P);
        auto ms_naive = naive_matching_statistics(T, P);
        assert(ms_idx.size() == ms_naive.size());
        for (size_t i = 0; i < P.size(); ++i)
            assert(ms_idx[i] == ms_naive[i] && "index must match naive");
    }
    std::cout << "  PASSED" << std::endl;
    return true;
}

/**
 * Test ms_query vs naive matching statistics on a 100 randomly extracted
 * "errory" patterns, with 30% "error" rate.
 */
bool test_ms_query_vs_naive_tsv_reposition(const std::string& data_dir) {
    std::cout << "Testing ms_query vs naive (extract_errory_string, with reposition)" << std::endl;
    std::string path = data_dir + "/minishred1_20_002_lcp.tsv";
    auto idx = ms_io::build_ms_index_spill_from_tsv<false>(path);
    if (!idx) {
        std::cout << "  DID NOT RUN" << std::endl;
        return false;
    }
    std::string T = reconstruct_text(*idx);
    std::mt19937 rng(123);
    const int num_strings = 100;
    const ulint pattern_len = 100;
    int checked = 0;
    for (int k = 0; k < num_strings; ++k) {
        std::string P = extract_errory_string(*idx, 0.3f, pattern_len, rng);
        if (T.find(P) == std::string::npos) continue;  // skip: reposition can yield non-substrings
        auto ms_idx = ms_query(*idx, P);
        auto ms_naive = naive_matching_statistics(T, P);
        assert(ms_idx.size() == ms_naive.size());
        for (size_t i = 0; i < P.size(); ++i)
            assert(ms_idx[i] == ms_naive[i] && "Mismatch in MS vectors");
        ++checked;
    }
    std::cout << "  PASSED" << std::endl;
    return true;
}

/**
 * Compare two indexes by running ms_query on many extract_errory_string patterns
 * and asserting identical matching statistics. Uses extractor for pattern generation.
 */
template <typename MSCallable1, typename MSCallable2>
static int compare_index_pair(MSCallable1 ms_a, MSCallable2 ms_b,
                             MSIndexSpillLCP<false>& extractor,
                             const std::string& T, int num_patterns,
                             float reposition_rate, ulint pattern_len, std::mt19937& rng) {
    int checked = 0;
    for (int k = 0; k < num_patterns; ++k) {
        std::string P = extract_errory_string(extractor, reposition_rate, pattern_len, rng);
        if (P.empty()) continue;
        if (reposition_rate > 0.0f && T.find(P) == std::string::npos) continue;
        auto a = ms_a(P);
        auto b = ms_b(P);
        assert(a.size() == b.size() && "MS length mismatch");
        for (size_t i = 0; i < P.size(); ++i)
            assert(a[i] == b[i] && "Mismatch in MS vectors");
        ++checked;
    }
    return checked;
}

/** Compare two indexes; return (checked_count, first_failing_pattern). For discover tool only. */
template <typename MSCallable1, typename MSCallable2>
static std::pair<int, std::string> compare_find_first_mismatch(
    MSCallable1 ms_a, MSCallable2 ms_b,
    MSIndexSpillLCP<false>& extractor,
    const std::string& T, int num_patterns,
    float reposition_rate, ulint pattern_len, std::mt19937& rng) {
    for (int k = 0; k < num_patterns; ++k) {
        std::string P = extract_errory_string(extractor, reposition_rate, pattern_len, rng);
        if (P.empty()) continue;
        if (reposition_rate > 0.0f && T.find(P) == std::string::npos) continue;
        auto a = ms_a(P);
        auto b = ms_b(P);
        if (a.size() != b.size()) return {k, P};
        for (size_t i = 0; i < P.size(); ++i)
            if (a[i] != b[i]) return {k, P};
    }
    return {num_patterns, ""};
}

/** Build an index with given options. */
using IndexBuilder = std::function<std::optional<MSIndexSpillLCP<false>>()>;

static const char* PATTERN_COAL_SPLIT = "AGAAAAGCCCTCTGATTTTTCACCAGGAGGGTAGTGGTAGTGGTATGGCCAAGAGTGGGGAGGTTGGATGCCCAGGCTCTATACATTCTGTGTGTGTGTG";
static const char* PATTERN_MULTISPILL2 = "CACAAAGTACATTCTCAAGAGTGGGGAGAATTACAAAGAATCTTCTTAAGGGTGGGGGAGATTACAAAGTACATTGATCAGTTAGGGTGGGGCAGGAACA";

/** Assert coal_split index matches naive on a specific pattern. */
bool test_coal_split_vs_naive_discrepant(const std::string& data_dir) {
    std::cout << "Testing coal_split vs naive (discrepant pattern)" << std::endl;
    std::string path = data_dir + "/minishred1_20_002_lcp.tsv";
    auto idx_base = ms_io::build_ms_index_spill_from_tsv<false>(path, ms_io::BuildOptions{});
    if (!idx_base) { std::cout << "  DID NOT RUN" << std::endl; return false; }
    std::string T = reconstruct_text(*idx_base);
    ms_io::BuildOptions o; o.coalesce = true; o.split_threshold = 1;
    auto idx = ms_io::build_ms_index_spill_from_tsv<false>(path, o);
    if (!idx) { std::cout << "  DID NOT RUN" << std::endl; return false; }
    std::string P(PATTERN_COAL_SPLIT);
    auto ms_idx = ms_query(*idx, P);
    auto ms_naive = naive_matching_statistics(T, P);
    assert(ms_idx == ms_naive && "coal_split must match naive");
    std::cout << "  PASSED" << std::endl;
    return true;
}

/** Assert multispill index matches naive on a specific pattern. */
bool test_multispill_vs_naive_discrepant(const std::string& data_dir, uchar spill_bits) {
    std::cout << "Testing multispill" << static_cast<unsigned>(spill_bits) << " vs naive (discrepant pattern)" << std::endl;
    std::string path = data_dir + "/minishred1_20_002_lcp.tsv";
    auto idx_base = ms_io::build_ms_index_spill_from_tsv<false>(path, ms_io::BuildOptions{});
    if (!idx_base) { std::cout << "  DID NOT RUN" << std::endl; return false; }
    std::string T = reconstruct_text(*idx_base);
    ms_io::BuildOptions o; o.coalesce = true; o.split_threshold = 1; o.spill_split_bits = spill_bits;
    auto idx = ms_io::build_ms_index_spill_from_tsv<false>(path, o);
    if (!idx) { std::cout << "  DID NOT RUN" << std::endl; return false; }
    std::string P(PATTERN_MULTISPILL2);
    auto ms_idx = ms_query(*idx, P);
    auto ms_naive = naive_matching_statistics(T, P);
    assert(ms_idx == ms_naive && "multispill must match naive");
    std::cout << "  PASSED" << std::endl;
    return true;
}

/**
 * Test that all index types produce identical matching statistics to base.
 * Asserts equality; fails on first mismatch.
 */
bool test_index_types_consistency(const std::string& data_dir) {
    std::cout << "Testing index type consistency (all must match base)" << std::endl;
    std::string path = data_dir + "/minishred1_20_002_lcp.tsv";
    std::vector<uchar> bwt_heads;
    std::vector<ulint> bwt_run_lengths;
    std::vector<std::vector<ulint>> lcps_per_run;
    if (!tsv::load_tsv(path, bwt_heads, bwt_run_lengths, lcps_per_run)) {
        std::cout << "  DID NOT RUN" << std::endl;
        return false;
    }

    const int num_patterns = 500;
    const ulint pattern_len = 100;
    const float reposition_rate = 0.0f;
    std::mt19937 rng(456);

    auto build_base = [&]() { return ms_io::build_ms_index_spill_from_tsv<false>(path, ms_io::BuildOptions{}); };
    auto build_split = [&]() { ms_io::BuildOptions o; o.split_threshold = 1; return ms_io::build_ms_index_spill_from_tsv<false>(path, o); };
    auto build_coal = [&]() { ms_io::BuildOptions o; o.coalesce = true; return ms_io::build_ms_index_spill_from_tsv<false>(path, o); };
    auto build_coal_split = [&]() { ms_io::BuildOptions o; o.coalesce = true; o.split_threshold = 1; return ms_io::build_ms_index_spill_from_tsv<false>(path, o); };
    auto build_multispill = [&](uchar bits) {
        ms_io::BuildOptions o; o.coalesce = true; o.split_threshold = 1; o.spill_split_bits = bits;
        return ms_io::build_ms_index_spill_from_tsv<false>(path, o);
    };

    auto idx_base = build_base();
    if (!idx_base) { std::cout << "  DID NOT RUN" << std::endl; return false; }
    std::string T = reconstruct_text(*idx_base);
    assert(!T.empty());

    auto run_comparison = [&](const char* name, std::optional<MSIndexSpillLCP<false>> other) {
        if (!other) return;
        auto ms_base = [&](const std::string& P) { return ms_query(*idx_base, P); };
        auto ms_other = [&](const std::string& P) { return ms_query(*other, P); };
        compare_index_pair(ms_base, ms_other, *idx_base, T, num_patterns, reposition_rate, pattern_len, rng);
    };

    auto idx_starts = ms_io::build_ms_index_spill_from_tsv<true>(path, ms_io::BuildOptions{});
    if (idx_starts) {
        auto ms_len = [&](const std::string& P) { return ms_query(*idx_base, P); };
        auto ms_starts = [&](const std::string& P) { return ms_query(*idx_starts, P); };
        compare_index_pair(ms_len, ms_starts, *idx_base, T, num_patterns, reposition_rate, pattern_len, rng);
    }
    run_comparison("split", build_split());
    run_comparison("coal", build_coal());
    run_comparison("coal_split", build_coal_split());
    for (uchar bits = 2; bits <= 10; ++bits)
        run_comparison("multispill", build_multispill(bits));

    std::cout << "  PASSED" << std::endl;
    return true;
}

}  // anonymous namespace

namespace ms_test {

bool run_all_tests(const std::string& data_dir) {
    bool all_ran = true;
    test_top_lcp<false>();
    std::cout << std::endl;
    test_top_lcp<true>();
    std::cout << std::endl;
    test_spill_lcp_queries<false>();
    std::cout << std::endl;
    test_spill_lcp_queries<true>();
    std::cout << std::endl;
    test_spillover<false>();
    std::cout << std::endl;
    test_spillover<true>();
    std::cout << std::endl;
    test_ms_query_consistency();
    std::cout << std::endl;
    test_build_spill_data_coalesce();
    std::cout << std::endl;
    test_naive_matching_statistics_known();
    std::cout << std::endl;
    if (!test_ms_query_vs_naive_single_example(data_dir)) all_ran = false;
    std::cout << std::endl;
    if (!test_ms_query_vs_naive_tsv(data_dir)) all_ran = false;
    std::cout << std::endl;
    if (!test_ms_query_vs_naive_tsv_reposition(data_dir)) all_ran = false;
    std::cout << std::endl;
    if (!test_index_types_consistency(data_dir)) all_ran = false;
    std::cout << std::endl;
    if (!test_coal_split_vs_naive_discrepant(data_dir)) all_ran = false;
    std::cout << std::endl;
    if (!test_multispill_vs_naive_discrepant(data_dir, 2)) all_ran = false;
    std::cout << std::endl;
    if (!test_coalesce_spillover_from_tsv(data_dir)) all_ran = false;
    std::cout << std::endl;
    if (all_ran) {
        std::cout << "All ms_test checks PASSED" << std::endl;
    } else {
        std::cout << "SOME TESTS NOT RUN" << std::endl;
    }
    return all_ran;
}

static std::optional<MSIndexSpillLCP<false>> build_index_by_type(const std::string& path, const std::string& index_type) {
    if (index_type == "base") {
        return ms_io::build_ms_index_spill_from_tsv<false>(path, ms_io::BuildOptions{});
    }
    if (index_type == "coal_split") {
        ms_io::BuildOptions o; o.coalesce = true; o.split_threshold = 1;
        return ms_io::build_ms_index_spill_from_tsv<false>(path, o);
    }
    if (index_type.compare(0, 10, "multispill") == 0 && index_type.size() >= 11) {
        uchar bits = static_cast<uchar>(std::stoul(index_type.substr(10)));
        ms_io::BuildOptions o; o.coalesce = true; o.split_threshold = 1; o.spill_split_bits = bits;
        return ms_io::build_ms_index_spill_from_tsv<false>(path, o);
    }
    return std::nullopt;
}

bool run_investigate(const std::string& data_dir, const std::string& pattern,
                    const std::string& index_type, int trace_reposition_at) {
    std::string path = data_dir + "/minishred1_20_002_lcp.tsv";
    auto idx = build_index_by_type(path, index_type);
    if (!idx) {
        std::cerr << "Failed to load index (type=" << index_type << ")" << std::endl;
        return false;
    }
    std::string T = reconstruct_text(*idx);
    if (T.find(pattern) == std::string::npos) {
        std::cerr << "Pattern not found in text" << std::endl;
        return false;
    }
    auto ms_idx = ms_query(*idx, pattern);
    auto ms_naive = naive_matching_statistics(T, pattern);
    bool match = (ms_idx == ms_naive);
    std::cout << "Index type: " << index_type << std::endl;
    std::cout << "ms_query vs naive: " << (match ? "MATCH" : "MISMATCH") << std::endl;
    std::cout << "\n--- Step-by-step trace ---" << std::endl;
    trace_ms_query_steps(*idx, pattern, std::cout, trace_reposition_at);
    return match;
}

/** Compute expected LCP values from raw TSV lcp array (before compression). */
static void expected_from_raw(const std::vector<ulint>& lcp, ulint offset,
                              ulint& boundary, ulint& row_min,
                              ulint& range_min_to_top, ulint& range_min_to_bottom) {
    boundary = lcp.empty() ? LCP_GAP : lcp[0];
    row_min = LCP_GAP;
    range_min_to_top = LCP_GAP;
    range_min_to_bottom = LCP_GAP;
    for (size_t i = 1; i < lcp.size(); ++i) {
        if (lcp[i] != LCP_GAP) row_min = std::min(row_min, lcp[i]);
    }
    for (size_t i = 0; i < lcp.size() && static_cast<ulint>(i) <= offset; ++i) {
        if (lcp[i] != LCP_GAP) range_min_to_top = std::min(range_min_to_top, lcp[i]);
    }
    for (size_t i = offset + 1; i < lcp.size(); ++i) {
        if (lcp[i] != LCP_GAP) range_min_to_bottom = std::min(range_min_to_bottom, lcp[i]);
    }
}

/** Compute expected LCP values from compressed array (after compress_lcps). */
static void expected_from_compressed(const std::vector<ulint>& full, ulint offset,
                                     ulint& boundary, ulint& row_min,
                                     ulint& range_min_to_top, ulint& range_min_to_bottom) {
    boundary = full.empty() ? LCP_GAP : full[0];
    row_min = LCP_GAP;
    range_min_to_top = LCP_GAP;
    range_min_to_bottom = LCP_GAP;
    for (size_t i = 1; i < full.size(); ++i) {
        if (full[i] != LCP_GAP) row_min = std::min(row_min, full[i]);
    }
    for (size_t i = 0; i < full.size() && static_cast<ulint>(i) <= offset; ++i) {
        if (full[i] != LCP_GAP) range_min_to_top = std::min(range_min_to_top, full[i]);
    }
    for (size_t i = offset + 1; i < full.size(); ++i) {
        if (full[i] != LCP_GAP) range_min_to_bottom = std::min(range_min_to_bottom, full[i]);
    }
}

/** Decode spillover for row i into (boundary, row_min, pairs). Pairs are (offset, value). */
static void decode_row_spillover(const MSIndexSpillLCP<false>& idx, ulint i,
                                 ulint& boundary, ulint& row_min,
                                 std::vector<std::pair<ulint, ulint>>& pairs) {
    boundary = LCP_GAP;
    row_min = LCP_GAP;
    pairs.clear();
    ulint top = idx.template get<LCPSpillRunCols::LCP_TOP>(i);
    ulint sub = idx.template get<LCPSpillRunCols::LCP_MIN_SUB>(i);
    ulint so = idx.template get<LCPSpillRunCols::LCP_SPILL>(i);
    const auto& spill = idx.spillover_for_row(i);
    bool jumbo = (top == idx.max_lcp_top() && sub == idx.max_lcp_min_sub());
    size_t p = 0;
    if (jumbo && so != NO_SPILL) {
        p = idx.spill_offset_bytes(so);
        auto [tv, np] = decode_uleb128(spill, p);
        boundary = tv;
        p = np;
        auto [m, nm] = decode_uleb128(spill, p);
        row_min = m;
        p = nm;
    } else {
        boundary = (top != LCP_GAP) ? top : LCP_GAP;
        row_min = (sub != 0 && top != LCP_GAP) ? (top - sub) : boundary;
        if (so == NO_SPILL) return;
        p = idx.spill_offset_bytes(so);
    }
    if (so != NO_SPILL) {
        auto [n, np] = decode_uleb128(spill, p);
        p = np;
        ulint run_len = idx.get_length(i);
        for (ulint k = 0; k < n && p < spill.size(); ++k) {
            if (k > run_len) break;
            auto [o, no] = decode_uleb128(spill, p);
            auto [vv, nv] = decode_uleb128(spill, no);
            p = nv;
            pairs.emplace_back(o, vv);
        }
    }
}

/** From decoded pairs, compute range_min for to_top [0,offset] and to_bottom (offset,run_len). */
static void expected_range_min_from_pairs(const std::vector<std::pair<ulint, ulint>>& pairs,
                                         ulint offset, ulint boundary,
                                         ulint& range_min_up, ulint& range_min_dn) {
    range_min_up = LCP_GAP;
    range_min_dn = LCP_GAP;
    if (boundary != LCP_GAP) range_min_up = boundary;
    for (const auto& [o, v] : pairs) {
        if (v == LCP_GAP) continue;
        if (o <= offset) range_min_up = std::min(range_min_up, v);
        if (offset == 0 ? o > 0 : o > offset) range_min_dn = std::min(range_min_dn, v);
    }
}

bool run_probe_base_internal(const std::string& data_dir, size_t first_row, size_t last_row) {
    std::string path = data_dir + "/minishred1_20_002_lcp.tsv";
    auto idx = ms_io::build_ms_index_spill_from_tsv<false>(path, ms_io::BuildOptions{});
    if (!idx) {
        std::cerr << "Failed to build base index" << std::endl;
        return false;
    }
    const size_t offset = 2;
    std::cout << "Base index internal consistency probe, rows " << first_row << "-" << last_row
              << ", offset=" << offset << std::endl;
    std::cout << "  Decode spillover -> expected boundary, row_min, range_min; compare with actual queries."
              << std::endl;
    std::cout << std::endl;

    std::vector<std::pair<ulint, ulint>> pairs;
    for (size_t r = first_row; r <= last_row && r < idx->move_runs(); ++r) {
        ulint exp_bnd, exp_rm;
        decode_row_spillover(*idx, r, exp_bnd, exp_rm, pairs);
        ulint exp_rt, exp_rb;
        expected_range_min_from_pairs(pairs, offset, exp_bnd, exp_rt, exp_rb);

        ulint act_bnd = boundary_lcp(*idx, r);
        ulint act_rm = row_min_lcp(*idx, r);
        ulint act_rt = range_min(*idx, r, offset, true);
        ulint act_rb = range_min(*idx, r, offset, false);

        bool bad_bnd = (act_bnd != exp_bnd);
        bool bad_rm = (act_rm != exp_rm);
        bool bad_rt = (act_rt != exp_rt);
        bool bad_rb = (act_rb != exp_rb);
        if (!bad_bnd && !bad_rm && !bad_rt && !bad_rb) continue;

        ulint len = idx->get_length(r);
        char c = static_cast<char>(idx->get_character(r));
        std::cout << "row " << r << " len=" << len << " c='" << c << "' INCONSISTENT" << std::endl;
        std::cout << "  boundary:   expected=" << (exp_bnd == LCP_GAP ? ulint(-1) : exp_bnd)
                  << " actual=" << (act_bnd == LCP_GAP ? ulint(-1) : act_bnd);
        if (bad_bnd) std::cout << "  WRONG";
        std::cout << std::endl;
        std::cout << "  row_min:    expected=" << (exp_rm == LCP_GAP ? ulint(-1) : exp_rm)
                  << " actual=" << (act_rm == LCP_GAP ? ulint(-1) : act_rm);
        if (bad_rm) std::cout << "  WRONG";
        std::cout << std::endl;
        std::cout << "  range_min_up:  expected=" << (exp_rt == LCP_GAP ? ulint(-1) : exp_rt)
                  << " actual=" << (act_rt == LCP_GAP ? ulint(-1) : act_rt);
        if (bad_rt) std::cout << "  WRONG";
        std::cout << std::endl;
        std::cout << "  range_min_dn:  expected=" << (exp_rb == LCP_GAP ? ulint(-1) : exp_rb)
                  << " actual=" << (act_rb == LCP_GAP ? ulint(-1) : act_rb);
        if (bad_rb) std::cout << "  WRONG";
        std::cout << std::endl;
        std::cout << std::endl;
    }
    return true;
}

bool run_probe_lcp_queries(const std::string& data_dir, size_t first_row, size_t last_row) {
    std::string path = data_dir + "/minishred1_20_002_lcp.tsv";
    std::vector<uchar> bwt_heads;
    std::vector<ulint> bwt_run_lengths;
    std::vector<std::vector<ulint>> lcps_per_run;
    if (!tsv::load_tsv(path, bwt_heads, bwt_run_lengths, lcps_per_run)) {
        std::cerr << "Failed to load TSV" << std::endl;
        return false;
    }
    auto idx = ms_io::build_ms_index_spill_from_tsv<false>(path, ms_io::BuildOptions{});
    if (!idx) {
        std::cerr << "Failed to build base index" << std::endl;
        return false;
    }
    const size_t offset = 2;  /* the offset from the trace */
    std::cout << "LCP probe: base index, rows " << first_row << "-" << last_row
              << ", offset=" << offset << std::endl;
    std::cout << "  TSV runs=" << lcps_per_run.size() << " index runs=" << idx->move_runs() << std::endl;
    std::cout << "  (raw=from TSV before compress, compressed=after compress_lcps, index=actual query)"
              << std::endl;
    std::cout << std::endl;

    std::vector<ulint> full;
    for (size_t r = first_row; r <= last_row && r < lcps_per_run.size(); ++r) {
        const std::vector<ulint>* next = (r + 1 < lcps_per_run.size()) ? &lcps_per_run[r + 1] : nullptr;
        full.clear();
        compress_lcps(lcps_per_run[r], next, full);

        ulint raw_bnd, raw_rm, raw_rt, raw_rb;
        expected_from_raw(lcps_per_run[r], offset, raw_bnd, raw_rm, raw_rt, raw_rb);

        ulint comp_bnd, comp_rm, comp_rt, comp_rb;
        expected_from_compressed(full, offset, comp_bnd, comp_rm, comp_rt, comp_rb);

        ulint idx_bnd = boundary_lcp(*idx, r);
        ulint idx_rm = row_min_lcp(*idx, r);
        ulint idx_rt = range_min(*idx, r, offset, true);
        ulint idx_rb = range_min(*idx, r, offset, false);

        ulint len = idx->get_length(r);
        char c = static_cast<char>(idx->get_character(r));
        std::cout << "row " << r << " len=" << len << " c='" << c << "'" << std::endl;
        std::cout << "  boundary:      raw=" << (raw_bnd == LCP_GAP ? ulint(-1) : raw_bnd)
                  << " compressed=" << (comp_bnd == LCP_GAP ? ulint(-1) : comp_bnd)
                  << " index=" << (idx_bnd == LCP_GAP ? ulint(-1) : idx_bnd);
        if (idx_bnd != comp_bnd) std::cout << "  MISMATCH";
        std::cout << std::endl;
        std::cout << "  row_min:       raw=" << (raw_rm == LCP_GAP ? ulint(-1) : raw_rm)
                  << " compressed=" << (comp_rm == LCP_GAP ? ulint(-1) : comp_rm)
                  << " index=" << (idx_rm == LCP_GAP ? ulint(-1) : idx_rm);
        if (idx_rm != comp_rm) std::cout << "  MISMATCH";
        std::cout << std::endl;
        std::cout << "  range_min_up:  raw=" << (raw_rt == LCP_GAP ? ulint(-1) : raw_rt)
                  << " compressed=" << (comp_rt == LCP_GAP ? ulint(-1) : comp_rt)
                  << " index=" << (idx_rt == LCP_GAP ? ulint(-1) : idx_rt);
        if (idx_rt != comp_rt) std::cout << "  MISMATCH";
        std::cout << std::endl;
        std::cout << "  range_min_dn:  raw=" << (raw_rb == LCP_GAP ? ulint(-1) : raw_rb)
                  << " compressed=" << (comp_rb == LCP_GAP ? ulint(-1) : comp_rb)
                  << " index=" << (idx_rb == LCP_GAP ? ulint(-1) : idx_rb);
        if (idx_rb != comp_rb) std::cout << "  MISMATCH";
        std::cout << std::endl;
        std::cout << std::endl;
    }
    return true;
}

bool run_discover_failing_patterns(const std::string& data_dir) {
    std::string path = data_dir + "/minishred1_20_002_lcp.tsv";
    std::vector<uchar> bwt_heads;
    std::vector<ulint> bwt_run_lengths;
    std::vector<std::vector<ulint>> lcps_per_run;
    if (!tsv::load_tsv(path, bwt_heads, bwt_run_lengths, lcps_per_run)) {
        std::cerr << "Failed to load TSV" << std::endl;
        return false;
    }
    const int num_patterns = 500;
    const ulint pattern_len = 100;
    const float reposition_rate = 0.0f;
    std::mt19937 rng(456);

    auto build_base = [&]() { return ms_io::build_ms_index_spill_from_tsv<false>(path, ms_io::BuildOptions{}); };
    auto build_coal_split = [&]() {
        ms_io::BuildOptions o; o.coalesce = true; o.split_threshold = 1;
        return ms_io::build_ms_index_spill_from_tsv<false>(path, o);
    };
    auto build_multispill = [&](uchar bits) {
        ms_io::BuildOptions o; o.coalesce = true; o.split_threshold = 1; o.spill_split_bits = bits;
        return ms_io::build_ms_index_spill_from_tsv<false>(path, o);
    };

    auto idx_base = build_base();
    if (!idx_base) return false;
    std::string T = reconstruct_text(*idx_base);

    auto discover = [&](const char* name, std::optional<MSIndexSpillLCP<false>> other,
                       std::function<std::optional<MSIndexSpillLCP<false>>()> build_other) {
        if (!other) return;
        auto ms_base = [&](const std::string& P) { return ms_query(*idx_base, P); };
        auto ms_other = [&](const std::string& P) { return ms_query(*other, P); };
        auto [c, fail_pattern] = compare_find_first_mismatch(ms_base, ms_other, *idx_base, T, num_patterns, reposition_rate, pattern_len, rng);
        if (!fail_pattern.empty()) {
            std::cout << name << " pattern=\"" << fail_pattern << "\"" << std::endl;
        }
    };

    discover("coal_split", build_coal_split(), build_coal_split);
    for (uchar bits = 2; bits <= 10; ++bits) {
        char buf[32];
        snprintf(buf, sizeof(buf), "multispill%u", static_cast<unsigned>(bits));
        discover(buf, build_multispill(bits), [&, bits]() { return build_multispill(bits); });
    }
    return true;
}

}  // namespace ms_test
