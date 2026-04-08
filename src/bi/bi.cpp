/**
 * Bidirectional BWT search tool — command-line interface.
 *
 * Subcommands:
 *   test [DATA_DIR]   Run built-in tests
 *   build TSV INDEX [--kmer-tables K1,K2,...]  Build from TSV
 *   add-table INDEX K  Add k-mer table to existing index
 *   mem INDEX PATTERN [--kmer-tables K1,K2,...]  Find SMEMs
 *   approx INDEX_PATH PATTERN [--k 1]  Approximate match
 *
 * Author: Ben Langmead (ben.langmead@gmail.com)
 * Date: Feb 23, 2026
 */

#include "bi_fmd.hpp"
#include "bi_io.hpp"
#include "bi_mem.hpp"
#include "bi_scheme.hpp"
#include "bi_test.hpp"
#include "tsv.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <cstring>

/** Print usage message. */
static void usage(const char* prog) {
    std::cerr << "Usage: " << prog << " <cmd> [args...]\n"
              << "\n"
              << "Main commands:\n"
              << "  build      TSV_PATH INDEX_PATH [--kmer-tables K1,K2,...]\n"
              << "              Build index from TSV.\n"
              << "  add-table  INDEX_PATH K\n"
              << "              Add k-mer table of length K to existing index.\n"
              << "  mem        INDEX_PATH PATTERN [--kmer-tables K1,K2,...]\n"
              << "              Find super-maximal exact matches.\n"
              << "  approx     INDEX_PATH PATTERN [--k 1]\n"
              << "              Approximate match with k substitutions (only --k 1 supported).\n"
              << "\n"
              << "Diagnostics and debugging:\n"
              << "  test       [DATA_DIR]\n"
              << "              Run built-in tests (default DATA_DIR is ./data).\n";
}

/** Parse a comma-separated list of k-mer table lengths. */
static std::vector<int> parse_kmer_tables(const char* arg) {
    std::vector<int> out;
    std::istringstream iss(arg);
    std::string tok;
    while (std::getline(iss, tok, ',')) {
        if (!tok.empty())
            out.push_back(std::stoi(tok));
    }
    return out;
}

/** Entry point; parsing and dispatch to subcommand handlers. */
int main(int argc, char** argv) {
    const char* prog = (argc > 0) ? argv[0] : "bi";
    if (argc < 2) {
        usage(prog);
        return 1;
    }
    const std::string cmd = argv[1];
    argc -= 2;
    argv += 2;

    if (cmd == "test") {
        std::string data_dir = (argc > 0 && argv[0][0] != '-') ? argv[0] : "./data";
        return bi_test::run_all_tests(data_dir) ? 0 : 1;
    }

    if (cmd == "build") {
        if (argc < 2) {
            std::cerr << "build requires TSV_PATH and INDEX_PATH\n";
            return 1;
        }
        std::string tsv_path = argv[0], idx_path = argv[1];
        std::vector<int> kmer_ks = {10, 6, 0};
        for (int i = 2; i < argc; ++i) {
            if (strcmp(argv[i], "--kmer-tables") == 0 && i + 1 < argc) {
                kmer_ks = parse_kmer_tables(argv[++i]);
                break;
            }
        }
        std::vector<uchar> bwt_heads;
        std::vector<ulint> bwt_run_lengths;
        std::vector<std::vector<ulint>> lcps_per_run;
        if (!tsv::load_tsv(tsv_path, bwt_heads, bwt_run_lengths, lcps_per_run)) {
            std::cerr << "Failed to load TSV: " << tsv_path << "\n";
            return 1;
        }
        FMDIndex fmd(bwt_heads, bwt_run_lengths);
        if (!bi_io::write_base(idx_path, fmd)) {
            std::cerr << "Failed to write base index: " << idx_path << "\n";
            return 1;
        }
        std::cout << "Built base index: " << idx_path << " (n=" << fmd.domain() << ", r=" << fmd.num_runs() << ")\n";
        for (int K : kmer_ks) {
            KmerTable tab = build_kmer_table(fmd, K);
            std::string path = idx_path + ".k" + std::to_string(K);
            if (!bi_io::write_kmer_table(path, fmd.get_fingerprint(), tab, fmd)) {
                std::cerr << "Failed to write k-mer table: " << path << "\n";
                return 1;
            }
            std::cout << "  Wrote " << path << "\n";
        }
        return 0;
    }

    if (cmd == "add-table") {
        if (argc < 2) {
            std::cerr << "add-table requires INDEX_PATH and K\n";
            return 1;
        }
        std::string idx_path = argv[0];
        int K = std::stoi(argv[1]);
        if (K < 0) {
            std::cerr << "K must be non-negative\n";
            return 1;
        }
        auto opt = bi_io::read_base(idx_path);
        if (!opt) {
            std::cerr << "Failed to load index: " << idx_path << "\n";
            return 1;
        }
        FMDIndex& fmd = *opt;
        KmerTable tab = build_kmer_table(fmd, K);
        std::string path = idx_path + ".k" + std::to_string(K);
        if (!bi_io::write_kmer_table(path, fmd.get_fingerprint(), tab, fmd)) {
            std::cerr << "Failed to write k-mer table: " << path << "\n";
            return 1;
        }
        std::cout << "Wrote " << path << "\n";
        return 0;
    }

    if (cmd == "mem") {
        if (argc < 2) {
            std::cerr << "mem requires INDEX_PATH and PATTERN\n";
            return 1;
        }
        std::string idx_path = argv[0], pattern = argv[1];
        std::vector<int> kmer_ks;
        for (int i = 2; i < argc; ++i) {
            if (strcmp(argv[i], "--kmer-tables") == 0 && i + 1 < argc) {
                kmer_ks = parse_kmer_tables(argv[++i]);
                break;
            }
        }
        auto opt = bi_io::read_base(idx_path);
        if (!opt) {
            std::cerr << "Failed to load index: " << idx_path << "\n";
            return 1;
        }
        FMDIndex& fmd = *opt;
        if (kmer_ks.empty()) {
            auto tables = bi_io::load_available_tables(idx_path, fmd.get_fingerprint(), fmd);
            for (auto& t : tables) fmd.add_kmer_table(std::move(t));
        } else {
            auto tables = bi_io::load_tables(idx_path, fmd.get_fingerprint(), fmd, kmer_ks);
            for (auto& t : tables) fmd.add_kmer_table(std::move(t));
        }
        auto smems = bi_mem::find_smems(fmd, pattern);
        for (const auto& s : smems) {
            std::cout << s.start << "\t" << s.end << "\t" << pattern.substr(s.start, s.end - s.start)
                      << "\t" << s.bi.size << "\n";
        }
        return 0;
    }

    if (cmd == "approx") {
        if (argc < 2) {
            std::cerr << "approx requires INDEX_PATH and PATTERN\n";
            return 1;
        }
        std::string idx_path = argv[0], pattern = argv[1];
        int max_errors = 1;
        for (int i = 2; i < argc; ++i) {
            if (strcmp(argv[i], "--k") == 0 && i + 1 < argc) {
                max_errors = std::stoi(argv[++i]);
                break;
            }
        }
        auto opt = bi_io::read_base(idx_path);
        if (!opt) {
            std::cerr << "Failed to load index: " << idx_path << "\n";
            return 1;
        }
        FMDIndex& fmd = *opt;
        auto tables = bi_io::load_available_tables(idx_path, fmd.get_fingerprint(), fmd);
        for (auto& t : tables) fmd.add_kmer_table(std::move(t));
        std::vector<bi_approx::ApproxMatch> matches;
        if (max_errors == 1) {
            auto scheme = std::vector<bi_approx::Search>{
                {{0, 1}, {0, 0}, {0, 1}},
                {{1, 0}, {0, 0}, {0, 1}}
            };
            matches = bi_approx::approx_match_scheme(fmd, pattern, scheme);
        } else {
            std::cerr << "Only --k 1 supported\n";
            return 1;
        }
        for (const auto& m : matches) {
            std::cout << m.start << "\t" << m.end << "\t" << m.errors << "\t"
                      << pattern.substr(m.start, m.end - m.start) << "\t" << m.bi.size << "\n";
        }
        return 0;
    }

    std::cerr << "Unknown command: " << cmd << "\n";
    usage(prog);
    return 1;
}
