/**
 * Interface for MS test suite.
 *
 * Author: Ben Langmead (ben.langmead@gmail.com)
 * Date: Feb 17, 2026
 */

#ifndef _MS_TEST_HPP
#define _MS_TEST_HPP

#include <string>

namespace ms_test {

/**
 * Run all MS tests. Returns true if all pass.
 */
bool run_all_tests(const std::string& data_dir);

/**
 * Trace ms_query step-by-step for a pattern (no assertions). For debugging.
 * Uses minishred1_20_002_lcp.tsv from data_dir.
 * index_type: "base" (default), "coal_split", "multispill2".."multispill10".
 * trace_reposition_at: if >= 0, when we hit that backward step index, print
 *   detailed reposition path (min_up, min_down, LCP values per run).
 */
bool run_investigate(const std::string& data_dir, const std::string& pattern,
                    const std::string& index_type = "base", int trace_reposition_at = -1);

/**
 * Discover and print failing patterns (base vs other index types).
 * For adding to standalone tests. Returns true if run successfully.
 */
bool run_discover_failing_patterns(const std::string& data_dir);

/**
 * Probe LCP queries (boundary, row_min, range_min) for base index rows.
 * Compares index answers vs raw TSV and vs compressed-expected. Use to find
 * which query type yields incorrect answers. first_row and last_row inclusive.
 */
bool run_probe_lcp_queries(const std::string& data_dir, size_t first_row, size_t last_row);

/**
 * Probe base index internal consistency: decode spillover and verify that
 * boundary_lcp, row_min_lcp, and range_min return values consistent with
 * what's stored. Reports which of the 3 query types yield wrong answers.
 * first_row and last_row inclusive.
 */
bool run_probe_base_internal(const std::string& data_dir, size_t first_row, size_t last_row);

}  // namespace ms_test

#endif /* _MS_TEST_HPP */
