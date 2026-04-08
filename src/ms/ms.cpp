/**
 * Command line interface that displays usage and dispatches to the various
 * subcommands.  TODO: slicker command-line parsing.
 *
 * Author: Ben Langmead (ben.langmead@gmail.com)
 * Date: Feb 17, 2026
 */

#include "ms_rlbwt.hpp"
#include "tsv.hpp"
#include "serialize.hpp"
#include "inspect.hpp"
#include "ms_test.hpp"
#include <iostream>
#include <string>
#include <cstring>
#include <cstdlib>

/** Print usage message. */
static void usage(const char* prog) {
    std::cerr << "Usage: " << prog << " <cmd> [args...]\n"
              << "\n"
              << "Main commands:\n"
              << "  build      TSV_PATH INDEX_PATH [--percentile K] [--split-threshold N]\n"
              << "            [--coalesce-spillover] [--spill-align N] [--spill-split-bits X]\n"
              << "              Build index from TSV.\n"
              << "  ms         INDEX_PATH PATTERN\n"
              << "              Compute matching statistics for PATTERN using INDEX_PATH.\n"
              << "  inspect    INDEX_PATH [--spillover-tsv FILE]\n"
              << "              Inspect index and optionally output spillover TSV.\n"
              << "  lcp-list   TSV_PATH\n"
              << "              Print LCP columns per run from TSV.\n"
              << "\n"
              << "Diagnostics and debugging:\n"
              << "  test           [DATA_DIR]\n"
              << "                  Run built-in tests (default DATA_DIR is ./data).\n"
              << "  discover       [DATA_DIR]\n"
              << "                  Find patterns where index types disagree with base.\n"
              << "  investigate    DATA_DIR PATTERN [INDEX_TYPE] [TRACE_STEP]\n"
              << "                  Step-through trace (INDEX_TYPE: base|coal_split|multispill2..10).\n"
              << "  probe          DATA_DIR FIRST_ROW LAST_ROW\n"
              << "                  Compare boundary/row_min/range_min for raw vs compressed vs index.\n"
              << "  probe-internal DATA_DIR FIRST_ROW LAST_ROW\n"
              << "                  Base index: decode spillover, verify 3 query types consistent.\n";
}

/** Main entry point. */
int main(int argc, char** argv) {
    const char* prog = (argc > 0) ? argv[0] : "ms";
    if (argc < 2) {
        usage(prog);
        return 1;
    }
    const std::string cmd = argv[1];
    argc -= 2;
    argv += 2;

    if (cmd == "test") {
        std::string data_dir = (argc > 0 && argv[0][0] != '-') ? argv[0] : "./data";
        return ms_test::run_all_tests(data_dir) ? 0 : 1;
    }

    if (cmd == "discover") {
        std::string data_dir = (argc > 0 && argv[0][0] != '-') ? argv[0] : "./data";
        return ms_test::run_discover_failing_patterns(data_dir) ? 0 : 1;
    }

    if (cmd == "probe") {
        if (argc < 3) {
            std::cerr << "probe requires DATA_DIR FIRST_ROW LAST_ROW\n";
            return 1;
        }
        std::string data_dir = argv[0];
        size_t first_row = static_cast<size_t>(std::stoull(argv[1]));
        size_t last_row = static_cast<size_t>(std::stoull(argv[2]));
        return ms_test::run_probe_lcp_queries(data_dir, first_row, last_row) ? 0 : 1;
    }

    
    if (cmd == "probe-internal") {
        if (argc < 3) {
            std::cerr << "probe-internal requires DATA_DIR FIRST_ROW LAST_ROW\n";
            return 1;
        }
        std::string data_dir = argv[0];
        size_t first_row = static_cast<size_t>(std::stoull(argv[1]));
        size_t last_row = static_cast<size_t>(std::stoull(argv[2]));
        return ms_test::run_probe_base_internal(data_dir, first_row, last_row) ? 0 : 1;
    }

    if (cmd == "investigate") {
        if (argc < 2) {
            std::cerr << "investigate requires DATA_DIR and PATTERN [INDEX_TYPE] [TRACE_STEP]\n";
            return 1;
        }
        std::string data_dir = argv[0], pattern = argv[1];
        std::string index_type = (argc >= 3) ? argv[2] : "base";
        int trace_step = (argc >= 4) ? std::atoi(argv[3]) : -1;
        return ms_test::run_investigate(data_dir, pattern, index_type, trace_step) ? 0 : 1;
    }

    if (cmd == "build") {
        if (argc < 2) {
            std::cerr << "build requires TSV_PATH and INDEX_PATH\n";
            return 1;
        }
        std::string tsv_path = argv[0], idx_path = argv[1];
        double percentile_k = 0.98;
        ulint split_threshold = SPLIT_THRESHOLD_NEVER;
        bool coalesce = false;
        ulint spill_align = 0;
        uchar spill_split_bits = 0;
        for (int i = 2; i < argc; ++i) {
            if (strcmp(argv[i], "--percentile") == 0 && i + 1 < argc) {
                percentile_k = std::stod(argv[++i]);
            } else if (strcmp(argv[i], "--split-threshold") == 0 && i + 1 < argc) {
                split_threshold = std::stoull(argv[++i]);
            } else if (strcmp(argv[i], "--coalesce-spillover") == 0) {
                coalesce = true;
            } else if (strcmp(argv[i], "--spill-align") == 0 && i + 1 < argc) {
                spill_align = std::stoull(argv[++i]);
            } else if (strcmp(argv[i], "--spill-split-bits") == 0 && i + 1 < argc) {
                spill_split_bits = static_cast<uchar>(std::stoul(argv[++i]));
            }
        }
        std::vector<uchar> bwt_heads;
        std::vector<ulint> bwt_run_lengths;
        std::vector<std::vector<ulint>> lcps_per_run;
        if (!tsv::load_tsv(tsv_path, bwt_heads, bwt_run_lengths, lcps_per_run)) {
            std::cerr << "Failed to load TSV: " << tsv_path << "\n";
            return 1;
        }
        apply_lcp_splitting(bwt_heads, bwt_run_lengths, lcps_per_run, split_threshold);
        auto [run_data, spill_vectors, max_top, max_sub, skinny_count, jumbo_count] =
            build_spill_data(lcps_per_run, percentile_k, coalesce, false /* coalesce_lcp_separately */, split_threshold, spill_align, spill_split_bits);
        MSIndexSpillLCP<false> idx(bwt_heads, bwt_run_lengths, run_data, std::move(spill_vectors), max_top, max_sub, spill_align, spill_split_bits);
        if (!ms_serialize::write_index(idx_path, idx)) {
            std::cerr << "Failed to write index: " << idx_path << "\n";
            return 1;
        }
        std::cout << "Built index: " << idx_path << " (skinny=" << skinny_count << ", jumbo=" << jumbo_count << ")\n";
        return 0;
    }

    if (cmd == "ms") {
        if (argc < 2) {
            std::cerr << "ms requires INDEX_PATH and PATTERN\n";
            return 1;
        }
        const char *idx_path = argv[0], *pattern_str = argv[1];
        auto opt = ms_serialize::read_index(idx_path);
        if (!opt) {
            std::cerr << "Failed to load index: " << idx_path << "\n";
            return 1;
        }
        std::string pattern = pattern_str;
        auto ms = ms_query(*opt, pattern);
        std::cout << "ms_query(\"" << pattern << "\") = [";
        for (size_t i = 0; i < ms.size(); ++i) {
            if (i > 0) std::cout << ",";
            std::cout << ms[i];
        }
        std::cout << "]\n";
        return 0;
    }

    if (cmd == "inspect") {
        if (argc < 1) {
            std::cerr << "inspect requires INDEX_PATH\n";
            return 1;
        }
        auto opt = ms_serialize::read_index(argv[0]);
        if (!opt) {
            std::cerr << "Failed to load index: " << argv[0] << "\n";
            return 1;
        }
        ms_inspect::run_inspect(*opt);
        for (int i = 1; i < argc; ++i) {
            if (strcmp(argv[i], "--spillover-tsv") == 0 && i + 1 < argc) {
                std::string path = argv[++i];
                if (ms_inspect::run_spillover_tsv(*opt, path))
                    std::cout << "Spillover TSV written to " << path << "\n";
            }
        }
        return 0;
    }

    if (cmd == "lcp-list") {
        if (argc < 1) {
            std::cerr << "lcp-list requires TSV_PATH\n";
            return 1;
        }
        std::vector<uchar> bwt_heads;
        std::vector<ulint> bwt_run_lengths;
        std::vector<std::vector<ulint>> lcps_per_run;
        if (!tsv::load_tsv(argv[0], bwt_heads, bwt_run_lengths, lcps_per_run)) {
            std::cerr << "Failed to load TSV: " << argv[0] << "\n";
            return 1;
        }
        for (size_t i = 0; i < lcps_per_run.size(); ++i) {
            std::cout << i << ": ";
            for (size_t j = 0; j < lcps_per_run[i].size(); ++j) {
                if (j > 0) std::cout << ",";
                ulint lcp = lcps_per_run[i][j];
                std::cout << (lcp == LCP_GAP ? "-" : std::to_string(lcp));
            }
            std::cout << "\n";
        }
        return 0;
    }

    std::cerr << "Unknown command: " << cmd << "\n";
    usage(prog);
    return 1;
}
