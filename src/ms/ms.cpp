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
#include "rlbwt_io.hpp"
#include "tms_index.hpp"
#include "tms_query.hpp"
#include "tms_smem.hpp"
#include "tms_test.hpp"
#include "perf_counters.hpp"
#include <iostream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cctype>
#include <chrono>
#include <fstream>
#include <memory>
#include <optional>
#include <type_traits>

/** Print usage message. */
static void usage(const char* prog) {
    std::cerr << "Usage: " << prog << " <cmd> [args...]\n"
              << "\n"
              << "Main commands:\n"
              << "  build      TSV_PATH INDEX_PATH [--percentile K] [--split-threshold N]\n"
              << "            [--coalesce-spillover] [--spill-align N] [--spill-split-bits X]\n"
              << "            [--minima-only]\n"
              << "              Build index from TSV.\n"
              << "  build-rlbwt HEADS LENS MINIMA INDEX_PATH [--percentile K]\n"
              << "            [--coalesce-spillover] [--spill-align N] [--spill-split-bits X]\n"
              << "              Build index from an RLBWT (HEADS: one byte per run; LENS:\n"
              << "              fixed-width little-endian lengths) and a TeraLCP -ominima file.\n"
              << "  ms         INDEX_PATH PATTERN\n"
              << "              Compute matching statistics for PATTERN using INDEX_PATH.\n"
              << "  batch      INDEX_PATH READS [-o OUT] [--no-output] [--interleave K]\n"
              << "              Load the index once and compute matching statistics for every\n"
              << "              read in READS (FASTA, FASTQ, or one sequence per line; - for\n"
              << "              stdin).  Writes one line per read: name, tab, space-separated\n"
              << "              values.  --no-output skips writing (for timing).  --interleave K\n"
              << "              keeps K reads in flight, prefetching each one's next row\n"
              << "              (default 32); K = 0 queries one read at a time without\n"
              << "              prefetching (for tms-batch, with tms_query when it reports\n"
              << "              only lengths in psi mode, else the batched engine with one\n"
              << "              read).  Results do not depend on K.  A summary with\n"
              << "              query time per base goes to stderr.\n"
              << "  tms-build  HEADS LENS INDEX_PATH [--minima FILE] [--lf-split B] [--fl-split B]\n"
              << "            [--phi-split B] [--lcp-bin FILE] [--no-phi-inv]\n"
              << "              Build a tms index from an RLBWT: LF and psi, which need no LCP\n"
              << "              input, and with --minima (a TeraLCP -ominima file, of which only\n"
              << "              each run's top LCP is used) or --lcp-bin (one 64-bit LCP per row)\n"
              << "              also phi and phi_inv (phi_inv unless --no-phi-inv; it is needed\n"
              << "              only for tms-batch --report smem-all).  --lf-split, --fl-split\n"
              << "              and --phi-split set Orbit's balancing factor for that structure\n"
              << "              (0 = no splitting; defaults: LF 0, FL and phi Orbit's default\n"
              << "              length capping and balancing).\n"
              << "  tms-build-tsv TSV_PATH INDEX_PATH [--phi] [split options]\n"
              << "              Same, taking the runs and their top LCPs from a TSV.\n"
              << "  tms-batch  INDEX_PATH READS [-o OUT] [--no-output] [--interleave K]\n"
              << "            [--mode psi|phi|phiskip|dual] [--positions]\n"
              << "            [--report ms|smem-one|smem-all] [--min-smem-len T]\n"
              << "              As batch, with a tms index.  --mode sets how repositions\n"
              << "              compute LCEs (default psi; the others need phi).  --positions\n"
              << "              adds a field with an occurrence position for each value (-1\n"
              << "              where it is 0); it needs phi.  --report adds a field with the\n"
              << "              read's SMEMs (super-maximal exact matches), space-separated:\n"
              << "              ms (the default) adds none; smem-one gives each as i:L:p, with\n"
              << "              start i in the read, length L and one text position p where it\n"
              << "              occurs, and needs phi; smem-all gives each as i:L:c:p1,p2,...,\n"
              << "              with all c text positions in BWT row order, and needs phi and\n"
              << "              phi_inv.  --min-smem-len T keeps only SMEMs with L > T (default\n"
              << "              0, all of them) in both SMEM reports.  Results do not depend on\n"
              << "              the mode or K.\n"
              << "  tms-text   INDEX_PATH\n"
              << "              Print the indexed text, read back with LF.\n"
              << "  tms-inspect INDEX_PATH\n"
              << "              Print the structures' sizes and column widths.\n"
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

/**
 * Streaming reader for FASTA, FASTQ, or one sequence per line, chosen by the
 * first non-empty line.  Sequences are upper-cased; FASTA records may span
 * several lines.
 */
class ReadStream {
public:
    explicit ReadStream(std::istream& in) : in_(in) {}

    bool next(std::string& name, std::string& seq) {
        name.clear();
        seq.clear();
        std::string line;
        if (!have_line_) {
            do {
                if (!std::getline(in_, line)) return false;
            } while (line.empty());
            pending_ = line;
            have_line_ = true;
        }
        if (format_ == 0) format_ = (pending_[0] == '>') ? 1 : (pending_[0] == '@') ? 2 : 3;
        ++count_;
        if (format_ == 3) {
            seq = pending_;
            name = "read" + std::to_string(count_);
            have_line_ = false;
        } else if (format_ == 2) {
            name = header_name(pending_);
            if (!std::getline(in_, seq)) return false;
            std::getline(in_, line);  // +
            std::getline(in_, line);  // qualities
            have_line_ = false;
        } else {
            name = header_name(pending_);
            have_line_ = false;
            while (std::getline(in_, line)) {
                if (!line.empty() && line[0] == '>') { pending_ = line; have_line_ = true; break; }
                seq += line;
            }
        }
        if (!seq.empty() && seq.back() == '\r') seq.pop_back();
        for (auto& ch : seq) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
        return true;
    }

private:
    static std::string header_name(const std::string& h) {
        size_t end = h.find_first_of(" \t\r", 1);
        return h.substr(1, end == std::string::npos ? std::string::npos : end - 1);
    }
    std::istream& in_;
    std::string pending_;
    bool have_line_ = false;
    int format_ = 0;  // 1 FASTA, 2 FASTQ, 3 plain
    size_t count_ = 0;
};

static std::vector<ulint> query_one(MSIndexSpillLCP<false>& idx, const std::string& s) { return ms_query(idx, s); }
static void query_many(MSIndexSpillLCP<false>& idx, const std::vector<std::string>& p, size_t k,
                       std::vector<std::vector<ulint>>& out) {
    ms_query_batch(idx, p, k, out);
}
// tms-batch settings, from its command line.
static TmsMode g_tms_mode = TmsMode::PSI;
static bool g_tms_positions = false;
static TmsReport g_tms_report = TmsReport::MS;
static ulint g_tms_min_smem_len = 0;
static std::vector<std::vector<ulint>> g_tms_pos;
static std::vector<TmsSmemHits> g_tms_hits;

// With an SMEM report, the SMEMs are found as part of the query.
static void query_many(TmsIndex& idx, const std::vector<std::string>& p, size_t k,
                       std::vector<std::vector<ulint>>& out) {
    const bool smems = g_tms_report != TmsReport::MS;
    tms_query_batch(idx, p, k, out, g_tms_mode, g_tms_positions || smems ? &g_tms_pos : nullptr);
    if (!smems) return;
    g_tms_hits.resize(p.size());
    for (size_t j = 0; j < p.size(); ++j)
        tms_report_smems(idx, out[j], g_tms_pos[j], g_tms_min_smem_len, g_tms_report, g_tms_hits[j]);
}
// Lengths alone in psi mode use the unbatched tms_query, as ms uses
// ms_query; anything else needs the batched engine, run with one read.
static std::vector<ulint> query_one(TmsIndex& idx, const std::string& s) {
    if (g_tms_mode == TmsMode::PSI && !g_tms_positions && g_tms_report == TmsReport::MS) return tms_query(idx, s);
    std::vector<std::vector<ulint>> len;
    query_many(idx, {s}, 1, len);
    return std::move(len[0]);
}

static std::optional<TmsIndex> read_tms_index(const std::string& path, bool with_phi_inv = true) {
    std::ifstream in(path, std::ios::binary);
    if (!in.good()) return std::nullopt;
    TmsIndex idx;
    try {
        idx.load(in, with_phi_inv);
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
        return std::nullopt;
    }
    return idx;
}

/**
 * ms batch and tms-batch: load the index once, then compute matching
 * statistics for every read.
 */
template <typename Index>
static int run_batch(int argc, char** argv, std::optional<Index> (*read)(const std::string&)) {
    if (argc < 2) {
        std::cerr << "batch requires INDEX_PATH and READS\n";
        return 1;
    }
    std::string idx_path = argv[0], reads_path = argv[1], out_path;
    bool write_output = true;
    size_t interleave = 32;
    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) out_path = argv[++i];
        else if (strcmp(argv[i], "--no-output") == 0) write_output = false;
        else if (strcmp(argv[i], "--interleave") == 0 && i + 1 < argc)
            interleave = static_cast<size_t>(std::stoull(argv[++i]));
        else if (std::is_same_v<Index, TmsIndex> && strcmp(argv[i], "--positions") == 0) g_tms_positions = true;
        else if (std::is_same_v<Index, TmsIndex> && strcmp(argv[i], "--mode") == 0 && i + 1 < argc) {
            const std::string m = argv[++i];
            if (m == "psi") g_tms_mode = TmsMode::PSI;
            else if (m == "phi") g_tms_mode = TmsMode::PHI;
            else if (m == "phiskip") g_tms_mode = TmsMode::PHISKIP;
            else if (m == "dual") g_tms_mode = TmsMode::DUAL;
            else { std::cerr << "Unknown mode: " << m << "\n"; return 1; }
        }
        else if (std::is_same_v<Index, TmsIndex> && strcmp(argv[i], "--report") == 0 && i + 1 < argc) {
            const std::string r = argv[++i];
            if (r == "ms") g_tms_report = TmsReport::MS;
            else if (r == "smem-one") g_tms_report = TmsReport::SMEM_ONE;
            else if (r == "smem-all") g_tms_report = TmsReport::SMEM_ALL;
            else { std::cerr << "Unknown report: " << r << "\n"; return 1; }
        }
        else if (std::is_same_v<Index, TmsIndex> && strcmp(argv[i], "--min-smem-len") == 0 && i + 1 < argc)
            g_tms_min_smem_len = static_cast<ulint>(std::stoull(argv[++i]));
        else { std::cerr << "Unknown batch option: " << argv[i] << "\n"; return 1; }
    }
    using clock = std::chrono::steady_clock;
    auto t0 = clock::now();
    auto opt = read(idx_path);
    if (!opt) {
        std::cerr << "Failed to load index: " << idx_path << "\n";
        return 1;
    }
    const double load_s = std::chrono::duration<double>(clock::now() - t0).count();
    if constexpr (std::is_same_v<Index, TmsIndex>) {
        if ((g_tms_positions || g_tms_report != TmsReport::MS || g_tms_mode != TmsMode::PSI) && !opt->has_phi()) {
            std::cerr << "--positions, --report smem-one/smem-all and modes other than psi need an index with phi\n";
            return 1;
        }
        if (g_tms_report == TmsReport::SMEM_ALL && !opt->has_phi_inv()) {
            std::cerr << "--report smem-all needs an index with phi_inv\n";
            return 1;
        }
    }

    std::ifstream fin;
    std::istream* in = &std::cin;
    if (reads_path != "-") {
        fin.open(reads_path);
        if (!fin.good()) { std::cerr << "Failed to open reads: " << reads_path << "\n"; return 1; }
        in = &fin;
    }
    std::ofstream fout;
    std::ostream* out = &std::cout;
    if (write_output && !out_path.empty()) {
        fout.open(out_path);
        if (!fout.good()) { std::cerr << "Failed to open output: " << out_path << "\n"; return 1; }
        out = &fout;
    }
    std::ios::sync_with_stdio(false);

    ReadStream reads(*in);
    std::string name, seq, line;
    size_t n_reads = 0, n_bases = 0;
    double query_s = 0.0;
    PerfCounters perf;
    auto t_all = clock::now();
    // With --positions, a tab-separated field holds the positions (-1 for
    // none); with an SMEM report, a last one holds the SMEMs.
    auto write_ms = [&](const std::string& nm, const std::vector<ulint>& ms, size_t j) {
        line.clear();
        line += nm;
        line += '\t';
        for (size_t i = 0; i < ms.size(); ++i) {
            if (i > 0) line += ' ';
            line += std::to_string(ms[i]);
        }
        if (std::is_same_v<Index, TmsIndex> && g_tms_positions) {
            line += '\t';
            const auto& pos = g_tms_pos[j];
            for (size_t i = 0; i < pos.size(); ++i) {
                if (i > 0) line += ' ';
                line += pos[i] == TMS_NO_POS ? std::string("-1") : std::to_string(pos[i]);
            }
        }
        if constexpr (std::is_same_v<Index, TmsIndex>) {
            if (g_tms_report != TmsReport::MS) {
                line += '\t';
                tms_format_smems(g_tms_hits[j], g_tms_report, line);
            }
        }
        line += '\n';
        out->write(line.data(), static_cast<std::streamsize>(line.size()));
    };
    if (interleave == 0) {
        // One read at a time with ms_query.
        while (reads.next(name, seq)) {
            auto tq = clock::now();
            perf.start();
            auto ms = query_one(*opt, seq);
            perf.stop();
            query_s += std::chrono::duration<double>(clock::now() - tq).count();
            ++n_reads;
            n_bases += seq.size();
            if (write_output) write_ms(name, ms, 0);
        }
    } else {
        // Blocks of reads with ms_query_batch, interleave reads in flight.
        const size_t block = std::max<size_t>(4096, 64 * interleave);
        std::vector<std::string> names, seqs;
        std::vector<std::vector<ulint>> results;
        bool more = true;
        while (more) {
            names.clear();
            seqs.clear();
            while (seqs.size() < block && (more = reads.next(name, seq))) {
                names.push_back(name);
                seqs.push_back(seq);
                n_bases += seq.size();
            }
            if (seqs.empty()) break;
            auto tq = clock::now();
            perf.start();
            query_many(*opt, seqs, interleave, results);
            perf.stop();
            query_s += std::chrono::duration<double>(clock::now() - tq).count();
            n_reads += seqs.size();
            if (write_output)
                for (size_t j = 0; j < seqs.size(); ++j) write_ms(names[j], results[j], j);
        }
    }
    out->flush();
#ifdef TMS_STATS
    std::cerr << "stats: bases=" << tms_stats.bases << " repositions/base=" << double(tms_stats.repositions) / tms_stats.bases
              << " psi_steps/base=" << double(tms_stats.psi_steps) / tms_stats.bases
              << " scan_rows/rep=" << double(tms_stats.scan_rows) / tms_stats.repositions
              << " len/rep=" << double(tms_stats.len_at_rep) / tms_stats.repositions
              << " phi_steps/base=" << double(tms_stats.phi_steps) / tms_stats.bases
              << " dist/rep=" << double(tms_stats.dist) / tms_stats.repositions
              << " dist1_frac=" << double(tms_stats.dist1) / tms_stats.repositions
              << " lce/rep=" << double(tms_stats.lce) / tms_stats.repositions
              << " capped_frac=" << double(tms_stats.lce_capped) / tms_stats.repositions
              << " scan_visits/rep=" << double(tms_stats.scan_visits) / tms_stats.repositions
              << " lf_ff/step=" << double(tms_stats.lf_ff) / tms_stats.lf_steps
              << " walk_visits/base=" << double(tms_stats.walk_visits) / tms_stats.bases << "\n";
#endif
    perf.report(std::cerr, n_bases);
    const double total_s = std::chrono::duration<double>(clock::now() - t_all).count();
    std::cerr << "batch: reads=" << n_reads << " bases=" << n_bases
              << " index_load_s=" << load_s << " query_s=" << query_s
              << " total_s=" << total_s
              << " query_ns_per_base=" << (n_bases ? query_s * 1e9 / n_bases : 0.0) << "\n";
    return 0;
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

#ifdef NDEBUG
    // The test suites check their results with assert.
    if (cmd == "test" || cmd == "tms-test") {
        std::cerr << cmd << " needs assertions; build with make (not make bench)\n";
        return 1;
    }
#endif

    if (cmd == "test") {
        std::string data_dir = (argc > 0 && argv[0][0] != '-') ? argv[0] : "./data";
        const bool ms_ok = ms_test::run_all_tests(data_dir);
        std::cout << std::endl;
        const bool tms_ok = tms_test::run_all_tests(data_dir);
        return (ms_ok && tms_ok) ? 0 : 1;
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
        bool minima_only = false;
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
            } else if (strcmp(argv[i], "--minima-only") == 0) {
                minima_only = true;
            } else {
                std::cerr << "Unknown build option: " << argv[i] << "\n";
                return 1;
            }
        }
        std::vector<uchar> bwt_heads;
        std::vector<ulint> bwt_run_lengths;
        std::vector<std::vector<ulint>> lcps_per_run;
        if (!tsv::load_tsv(tsv_path, bwt_heads, bwt_run_lengths, lcps_per_run)) {
            std::cerr << "Failed to load TSV: " << tsv_path << "\n";
            return 1;
        }
        apply_lcp_splitting(bwt_heads, bwt_run_lengths, lcps_per_run, split_threshold, minima_only);
        auto [run_data, spill_vectors, max_top, max_sub, skinny_count, jumbo_count] =
            build_spill_data(lcps_per_run, percentile_k, coalesce, false /* coalesce_lcp_separately */, split_threshold, spill_align, spill_split_bits, minima_only);
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

    if (cmd == "build-rlbwt") {
        if (argc < 4) {
            std::cerr << "build-rlbwt requires HEADS, LENS, MINIMA and INDEX_PATH\n";
            return 1;
        }
        const std::string heads_path = argv[0], lens_path = argv[1], minima_path = argv[2], idx_path = argv[3];
        double percentile_k = 0.98;
        bool coalesce = false;
        ulint spill_align = 0;
        uchar spill_split_bits = 0;
        for (int i = 4; i < argc; ++i) {
            if (strcmp(argv[i], "--percentile") == 0 && i + 1 < argc) percentile_k = std::stod(argv[++i]);
            else if (strcmp(argv[i], "--coalesce-spillover") == 0) coalesce = true;
            else if (strcmp(argv[i], "--spill-align") == 0 && i + 1 < argc) spill_align = std::stoull(argv[++i]);
            else if (strcmp(argv[i], "--spill-split-bits") == 0 && i + 1 < argc)
                spill_split_bits = static_cast<uchar>(std::stoul(argv[++i]));
            else { std::cerr << "Unknown build-rlbwt option: " << argv[i] << "\n"; return 1; }
        }
        std::vector<uchar> heads;
        std::vector<ulint> lens;
        std::vector<RunLcpPairs> pairs;
        std::string err;
        if (!rlbwt_io::load_rlbwt(heads_path, lens_path, heads, lens, err)) {
            std::cerr << "Failed to load RLBWT: " << err << "\n";
            return 1;
        }
        ulint n = 0;
        for (ulint l : lens) n += l;
        if (!rlbwt_io::read_minima(minima_path, heads.size(), n, pairs, err)) {
            std::cerr << "Failed to load minima: " << err << "\n";
            return 1;
        }
        auto [run_data, spill_vectors, max_top, max_sub, skinny_count, jumbo_count] =
            build_spill_data_from_pairs(std::move(pairs), percentile_k, coalesce, false, spill_align, spill_split_bits);
        MSIndexSpillLCP<false> idx(heads, lens, run_data, std::move(spill_vectors), max_top, max_sub, spill_align,
                                   spill_split_bits);
        if (!ms_serialize::write_index(idx_path, idx)) {
            std::cerr << "Failed to write index: " << idx_path << "\n";
            return 1;
        }
        std::cout << "Built index: " << idx_path << " (runs=" << heads.size() << ", n=" << n
                  << ", skinny=" << skinny_count << ", jumbo=" << jumbo_count << ")\n";
        return 0;
    }

    if (cmd == "batch") {
        return run_batch<MSIndexSpillLCP<false>>(argc, argv, ms_serialize::read_index);
    }

    if (cmd == "tms-batch") {
        // Only an smem-all report walks phi_inv, so other reports skip loading it.
        return run_batch<TmsIndex>(argc, argv, [](const std::string& path) {
            return read_tms_index(path, g_tms_report == TmsReport::SMEM_ALL);
        });
    }

    if (cmd == "tms-build" || cmd == "tms-build-tsv") {
        const bool from_tsv = cmd == "tms-build-tsv";
        const int npos = from_tsv ? 2 : 3;
        if (argc < npos) {
            std::cerr << cmd << (from_tsv ? " requires TSV_PATH and INDEX_PATH\n" : " requires HEADS, LENS and INDEX_PATH\n");
            return 1;
        }
        TmsBuildOptions opts;
        bool with_phi = false;
        std::string minima_path, lcp_bin_path;
        auto split_arg = [](const char* v) {
            const ulint b = std::stoull(v);
            return b == 0 ? orbit::NO_SPLITTING : orbit::split_params(orbit::DEFAULT_LENGTH_CAPPING, b);
        };
        for (int i = npos; i < argc; ++i) {
            if (strcmp(argv[i], "--lf-split") == 0 && i + 1 < argc) opts.lf_split = split_arg(argv[++i]);
            else if (strcmp(argv[i], "--fl-split") == 0 && i + 1 < argc) opts.fl_split = split_arg(argv[++i]);
            else if (strcmp(argv[i], "--phi-split") == 0 && i + 1 < argc) opts.phi_split = split_arg(argv[++i]);
            else if (strcmp(argv[i], "--no-phi-inv") == 0) opts.phi_inv = false;
            else if (from_tsv && strcmp(argv[i], "--phi") == 0) with_phi = true;
            else if (!from_tsv && strcmp(argv[i], "--minima") == 0 && i + 1 < argc) { minima_path = argv[++i]; with_phi = true; }
            else if (!from_tsv && strcmp(argv[i], "--lcp-bin") == 0 && i + 1 < argc) { lcp_bin_path = argv[++i]; with_phi = true; }
            else { std::cerr << "Unknown " << cmd << " option: " << argv[i] << "\n"; return 1; }
        }
        using clock = std::chrono::steady_clock;
        auto t0 = clock::now();
        std::vector<uchar> heads;
        std::vector<ulint> lens, tops;
        if (from_tsv) {
            std::vector<std::vector<ulint>> lcps;
            if (!tsv::load_tsv(argv[0], heads, lens, lcps)) {
                std::cerr << "Failed to load TSV: " << argv[0] << "\n";
                return 1;
            }
            if (with_phi)
                for (const auto& l : lcps) tops.push_back(l[0]);
        } else {
            std::string err;
            if (!rlbwt_io::load_rlbwt(argv[0], argv[1], heads, lens, err)) {
                std::cerr << "Failed to load RLBWT: " << err << "\n";
                return 1;
            }
            if (with_phi && !lcp_bin_path.empty()) {
                // One little-endian 64-bit LCP per row, optionally after a
                // 64-bit row count (msbench prep's form); keep each run head's.
                std::ifstream lin(lcp_bin_path, std::ios::binary | std::ios::ate);
                if (!lin.good()) { std::cerr << "cannot open " << lcp_bin_path << "\n"; return 1; }
                ulint n = 0;
                for (ulint l : lens) n += l;
                const ulint bytes = static_cast<ulint>(lin.tellg());
                lin.seekg(0);
                if (bytes == 8 * (n + 1)) {
                    uint64_t count = 0;
                    lin.read(reinterpret_cast<char*>(&count), 8);
                    if (count != n) { std::cerr << "lcp file count " << count << " is not n = " << n << "\n"; return 1; }
                } else if (bytes != 8 * n) {
                    std::cerr << "lcp file has " << bytes << " bytes; expected 8 per row for n = " << n << "\n";
                    return 1;
                }
                std::vector<uint64_t> buf(1 << 20);
                size_t have = 0, at = 0;
                tops.reserve(heads.size());
                for (ulint len : lens) {
                    for (ulint j = 0; j < len; ++j) {
                        if (at == have) {
                            lin.read(reinterpret_cast<char*>(buf.data()), static_cast<std::streamsize>(buf.size() * 8));
                            have = static_cast<size_t>(lin.gcount()) / 8;
                            at = 0;
                            if (have == 0) { std::cerr << "lcp file is shorter than the RLBWT\n"; return 1; }
                        }
                        if (j == 0) tops.push_back(buf[at]);
                        ++at;
                    }
                }
            } else if (with_phi) {
                ulint n = 0;
                for (ulint l : lens) n += l;
                std::vector<RunLcpPairs> pairs;
                if (!rlbwt_io::read_minima(minima_path, heads.size(), n, pairs, err)) {
                    std::cerr << "Failed to load minima: " << err << "\n";
                    return 1;
                }
                tops.reserve(pairs.size());
                for (const auto& p : pairs) tops.push_back(p[0].second);
            }
        }
        const double load_s = std::chrono::duration<double>(clock::now() - t0).count();
        auto t1 = clock::now();
        TmsIndex idx(heads, lens, opts, with_phi ? &tops : nullptr);
        const double build_s = std::chrono::duration<double>(clock::now() - t1).count();
        const std::string idx_path = argv[npos - 1];
        std::ofstream out(idx_path, std::ios::binary);
        const size_t bytes = out.good() ? idx.serialize(out) : 0;
        if (!out.good()) {
            std::cerr << "Failed to write index: " << idx_path << "\n";
            return 1;
        }
        idx.describe(std::cout);
        std::cout << "Built tms index: " << idx_path << " (runs=" << heads.size() << ", bytes=" << bytes
                  << ", load_s=" << load_s << ", build_s=" << build_s << ")\n";
        return 0;
    }

    if (cmd == "tms-inspect") {
        if (argc < 1) { std::cerr << "tms-inspect requires INDEX_PATH\n"; return 1; }
        auto opt = read_tms_index(argv[0]);
        if (!opt) { std::cerr << "Failed to load index: " << argv[0] << "\n"; return 1; }
        opt->describe(std::cout);
        return 0;
    }

    if (cmd == "tms-text") {
        if (argc < 1) { std::cerr << "tms-text requires INDEX_PATH\n"; return 1; }
        auto opt = read_tms_index(argv[0]);
        if (!opt) { std::cerr << "Failed to load index: " << argv[0] << "\n"; return 1; }
        // Row 0 holds text position n - 1, so the s-th row LF visits from it
        // has BWT character T[n - 2 - s].
        const ulint n = opt->domain();
        std::string t(n, '\0');
        auto pos = opt->first();
        for (ulint st = 0; st < n; ++st) {
            const uchar c = opt->get_character(pos.interval);
            t[(2 * n - 2 - st) % n] = c == orbit::TERMINATOR ? '$' : c == orbit::SEPARATOR ? '%' : static_cast<char>(c);
            pos = opt->LF_step(pos);
        }
        std::cout << t << "\n";
        return 0;
    }

    if (cmd == "tms-test") {
        std::string data_dir = (argc > 0 && argv[0][0] != '-') ? argv[0] : "./data";
        return tms_test::run_all_tests(data_dir) ? 0 : 1;
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
