/**
 * Command line interface that displays usage and dispatches to the various
 * subcommands.  `ms --help` lists the user commands and options, and
 * `ms --help-advanced` adds the development and advanced ones.
 *
 * Author: Ben Langmead (ben.langmead@gmail.com)
 * Date: Feb 17, 2026
 */

#include "ms_rlbwt.hpp"
#include "ms_build.hpp"
#include "serialize.hpp"
#include "inspect.hpp"
#include "ms_test.hpp"
#include "rlbwt_io.hpp"
#include "mem_report.hpp"
#include "tms_index.hpp"
#include "phi_samples_io.hpp"
#include "tms_query.hpp"
#include "tms_smem.hpp"
#include "tms_test.hpp"
#include "fmd_test.hpp"
#include "perf_counters.hpp"
#include "parallel_blocks.hpp"
#include "read_chunks.hpp"
#include "format_int.hpp"
#include "orient.hpp"
#include "ftab_survey.hpp"
#include "split_report.hpp"
#include "ftab_build.hpp"
#include "ftab_file.hpp"
#include <zlib.h>
#include <iostream>
#include <map>
#include <mutex>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cctype>
#include <atomic>
#include <chrono>
#include <fstream>
#include <functional>
#include <memory>
#include <optional>
#include <type_traits>
#include <filesystem>
#include <iomanip>
#include <thread>
#include <ctime>

// The commit this binary was built from (see the Makefile).
#ifndef MS_COMMIT
#define MS_COMMIT "unknown"
#endif

/**
 * The overview printed by `ms`, `ms --help` and `ms help`: user commands
 * and user options.  With advanced, it leaves off its last line, which
 * points to --help-advanced.
 */
static void usage(std::ostream& os, bool advanced = false) {
    os << "ms: matching statistics over Orbit move structures (commit " << MS_COMMIT << ")\n"
       << "\n"
       << "Usage: ms <command> [args...]\n"
       << "Options take their value as --name VALUE or --name=VALUE.\n"
       << "\n"
       << "Build an index:\n"
       << "  build-rlbwt HEADS LENS INDEX (--minima FILE | --lcp-bin FILE) [--lf-split FACTOR]\n"
       << "      Build an ms index from an RLBWT (HEADS: one byte per run; LENS:\n"
       << "      fixed-width little-endian run lengths) and LCP input: a TeraLCP\n"
       << "      -ominima file, or --lcp-bin, one 64-bit LCP per row, optionally after\n"
       << "      a 64-bit row count.  --lf-split sets Orbit's balancing factor for LF\n"
       << "      (positive; default: Orbit's default).  Results do not depend on it.\n"
       << "  tms-build HEADS LENS INDEX [--phi FILE | --minima FILE | --lcp-bin FILE]\n"
       << "            [--layout full|psi|phi] [--seq-ids FILE] [--no-phi-inv]\n"
       << "            [--lf-split FACTOR] [--phi-split FACTOR] [--lcp-summary]\n"
       << "      Build a tms index.  LF and psi need no LCP input; phi and phi_inv\n"
       << "      need one of --phi (a TeraLCP -ophi file, which allows one terminator\n"
       << "      per string, as in an FMD), --minima or --lcp-bin (which need one\n"
       << "      terminator).  --layout keeps what some queries read: full (default)\n"
       << "      keeps everything; psi keeps LF and psi and serves only --mode psi\n"
       << "      without positions or SMEMs; phi keeps LF, phi and phi_inv and serves\n"
       << "      every query but --mode psi, phiskip and dual.  --no-phi-inv drops\n"
       << "      phi_inv, which only --report smem-all reads.  --seq-ids stores\n"
       << "      ropebwt3's sequence-ID table (the FMD's .len.gz) for smem-all.\n"
       << "      --phi-split 0 leaves phi unsplit.  --lcp-summary adds each LF row's\n"
       << "      LCP at its head and least LCP below it, which --mode phisum reads;\n"
       << "      it needs --minima or --lcp-bin.\n"
       << "\n"
       << "Query:\n"
       << "  query INDEX READS [options]\n"
       << "      Compute matching statistics for every read in READS (FASTA, FASTQ,\n"
       << "      or one sequence per line; - for stdin) with an ms or tms index; its\n"
       << "      header says which.  Writes one line per read: its name, a tab, and\n"
       << "      its values, space-separated (QUERY_OUTPUT.md).  A summary goes to\n"
       << "      stderr.  Results and their order do not depend on the thread count.\n"
       << "    -o FILE                  write to FILE (default stdout); FILE.gz is gzipped\n"
       << "    -p N, --threads N        query on N threads (default 0: every CPU)\n"
       << "    --gzip-level LEVEL       gzip level for .gz output, 0 to 9 (default 1)\n"
       << "    --orient fwd|rc|qual     which strand of each read to match: fwd (default)\n"
       << "                             as given; rc its reverse complement; qual rc only\n"
       << "                             for reads whose last bases have higher quality\n"
       << "                             (FASTQ).  rc and qual equal fwd only on an index\n"
       << "                             holding both strands, as an FMD does.\n"
       << "    --ftab FILE              start reads KLEN steps in from a stored ftab\n"
       << "                             (ftab-build); it must be built from INDEX\n"
       << "    --build-ftab[=KLEN]      build an ftab of KLEN-mers after loading\n"
       << "                             (default KLEN 14)\n"
       << "  tms index only:\n"
       << "    --mode psi|phi|phisum|phiskip|dual\n"
       << "                             how repositions compute LCEs (default psi; psi,\n"
       << "                             phiskip and dual need psi, phi needs phi; phisum\n"
       << "                             walks phi only over the current row, taking whole\n"
       << "                             rows' LCP minima from an index built with\n"
       << "                             --lcp-summary)\n"
       << "    --positions              add a field with an occurrence position for each\n"
       << "                             value (-1 where it is 0); needs phi\n"
       << "    --report ms|smem-one|smem-all\n"
       << "                             ms (default) adds nothing; smem-one adds the\n"
       << "                             read's SMEMs as i:L:p; smem-all writes SMEMs as\n"
       << "                             ropebwt3 mem -p does, one line each\n"
       << "    --min-smem-len L         report SMEMs of length L or more, as ropebwt3's\n"
       << "                             -l L (default 1: all)\n"
       << "    --max-smem-positions MAXPOS\n"
       << "                             smem-all lists at most MAXPOS positions per SMEM\n"
       << "                             (default: all), as ropebwt3's -p\n"
       << "      --orient rc and qual cannot be used with --positions or an SMEM report.\n"
       << "\n"
       << "Ftabs (FTAB_FORMAT.md):\n"
       << "  ftab-build INDEX FTAB [--kmer-len KLEN | --auto] [-p N]\n"
       << "      Build an index's ftab, a table of the query state after each KLEN-mer\n"
       << "      over ACGT, and write it to FTAB.  --auto (the default) takes the\n"
       << "      largest KLEN whose table fits in 5% of the index file, then steps\n"
       << "      down while fewer than 80% of its KLEN-mers occur.  -p as for query.\n"
       << "  ftab-verify INDEX FTAB [--crc]\n"
       << "      Check from the headers that FTAB was built from INDEX; --crc also\n"
       << "      checks its payload CRC-32s.\n"
       << "  ftab-info FTAB\n"
       << "      Print an ftab's header.\n"
       << "\n"
       << "Index files (INDEX_FORMAT.md):\n"
       << "  index-plan INDEX [--mode M] [--positions] [--report R]\n"
       << "      Print an index's header and what each query loads, and its RAM; with\n"
       << "      query options, for that query alone.\n"
       << "  index-verify INDEX\n"
       << "      Check every section's CRC-32 and that the file is whole.\n"
       << "  index-add-seq-ids INDEX SEQ_IDS OUT\n"
       << "      Write the tms index INDEX to OUT with the sequence-ID table SEQ_IDS\n"
       << "      (see tms-build --seq-ids), copying its other sections unchanged.\n";
    if (!advanced) os << "\nRun 'ms --help-advanced' for development and advanced options.\n";
}

/** The overview followed by development and advanced options and commands. */
static void usage_advanced(std::ostream& os) {
    usage(os, true);
    os << "\n"
       << "Advanced query options:\n"
       << "    --interleave DEPTH       keep DEPTH reads in flight per thread, prefetching\n"
       << "                             each one's next row (default 32); 0 queries one\n"
       << "                             read at a time without prefetching\n"
       << "    --block-reads READS      each thread takes about READS reads at a time\n"
       << "                             (default max(512, 16 * DEPTH)), sized in bytes\n"
       << "                             from the reads seen so far\n"
       << "    --block-bytes BYTES      each thread takes BYTES bytes of input at a time\n"
       << "    --no-output              skip writing results (for timing)\n"
       << "    --phi-starts             tms: matching statistics without positions load\n"
       << "                             phi's starts-based form, not its narrower\n"
       << "                             lengths-based form\n"
       << "    --count walk|backward|adaptive\n"
       << "                             how smem-all counts an SMEM's occurrences: walk\n"
       << "                             steps phi and phi_inv over all of them; backward\n"
       << "                             finds its BWT interval by backward search, L steps\n"
       << "                             for length L, and lists positions from the\n"
       << "                             interval's top; adaptive (default) walks at most\n"
       << "                             B + A * L steps, then searches if the walk has not\n"
       << "                             ended\n"
       << "    --count-budget B         adaptive: B (default 8)\n"
       << "    --count-budget-per-len A adaptive: A (default 0.25)\n"
       << "    --count-short-len S      adaptive: search SMEMs of length at most S\n"
       << "                             without walking (default 0)\n"
       << "    --count-row-shift R      backward and adaptive sample every 2^R-th LF row\n"
       << "                             start in a table built after loading and timed\n"
       << "                             with it (default 4)\n"
       << "    --oracle-write FILE      record every reposition's outcome\n"
       << "                             (ORACLE_FORMAT.md); not with an ftab\n"
       << "    --mode oracle|oracle-scan --oracle FILE\n"
       << "                             replay recorded outcomes: oracle moves to the\n"
       << "                             recorded run and length; oracle-scan scans only\n"
       << "                             the recorded side and takes the recorded length.\n"
       << "                             Output equals the recording mode's.  An ms index's\n"
       << "                             other mode is standard, its default.\n"
       << "      Results do not depend on DEPTH, READS, BYTES, --phi-starts or --count.\n"
       << "      On a\n"
       << "      machine with several memory nodes, the index is interleaved over them\n"
       << "      (Linux) unless numactl sets a memory policy.  In the summary,\n"
       << "      query_ns_per_base is the threads' summed query time per base and\n"
       << "      wall_ns_per_base the elapsed time per base; flipped counts the reads\n"
       << "      --orient matched on the other strand.  With smem-all, DEPTH also sets\n"
       << "      how many phi and phi_inv walks list positions at once.\n"
       << "\n"
       << "Advanced build-rlbwt options:\n"
       << "    --minima-only            keep only the minima of --lcp-bin's LCPs, as a\n"
       << "                             --minima file does\n"
       << "    --percentile FRAC        spillover percentile (default 0.98)\n"
       << "    --coalesce-spillover     coalesce spillover\n"
       << "    --spill-align ALIGN      spillover alignment\n"
       << "    --spill-split-bits BITS  spillover split bits\n"
       << "\n"
       << "Development commands:\n"
       << "  pattern INDEX PATTERN\n"
       << "      Print the matching statistics of PATTERN with the ms index INDEX.\n"
       << "  tms-text INDEX\n"
       << "      Print the indexed text, read back with LF; a string collection prints\n"
       << "      each string and its $, one per line, in the BWT order of their\n"
       << "      terminators.  An FMD's N prints as N.\n"
       << "  tms-inspect INDEX [--phi-lengths]\n"
       << "      Print a tms index's structures' sizes and column widths.\n"
       << "      --phi-lengths loads only LF and phi's lengths-based form.\n"
       << "  inspect INDEX [--spillover-tsv FILE]\n"
       << "      Inspect an ms index and optionally write its spillover as TSV.\n"
       << "  ftab-survey INDEX [--kmin KLEN] [--kmax KLEN] [--verify KLEN] [-p N]\n"
       << "      Count the KLEN-mers over ACGT that occur in an ms index's text, for\n"
       << "      KLEN from --kmin (default 6) to --kmax (default 15), and the size of\n"
       << "      a dense table of them; --verify checks counts up to its KLEN by brute\n"
       << "      force.\n"
       << "  split-report INDEX [-p N]\n"
       << "      For each move structure, print its split parameters, intervals,\n"
       << "      longest interval, length cap and balance against the bound 2b - 1.\n"
       << "  index-rewrap LEGACY_INDEX INDEX\n"
       << "      Rewrite an index of an earlier format (MSIX v4, TMSX v8) in the\n"
       << "      current one, copying its structures' bytes.\n"
       << "  test [DATA_DIR]          run the built-in tests (default DATA_DIR ./data)\n"
       << "  tms-test [DATA_DIR]      run only the tms tests\n"
       << "  fmd-test [DATA_DIR]      run only the string-collection (FMD) and large-n tests\n"
       << "  discover [DATA_DIR]      find patterns where index types disagree with base\n"
       << "  investigate DATA_DIR PATTERN [INDEX_TYPE] [TRACE_STEP]\n"
       << "                           step-through trace (INDEX_TYPE:\n"
       << "                           base|coal_split|multispill2..10)\n"
       << "  probe DATA_DIR FIRST_ROW LAST_ROW\n"
       << "                           compare boundary, row_min and range_min for raw,\n"
       << "                           compressed and index\n"
       << "  probe-internal DATA_DIR FIRST_ROW LAST_ROW\n"
       << "                           base index: decode spillover and check that three\n"
       << "                           query types agree\n";
}

/**
 * If argv[i] is the option name, given as "name VALUE" or "name=VALUE", set
 * v to its value, advancing i past it, and return true.
 */
static bool opt_value(int argc, char** argv, int& i, const char* name, std::string& v) {
    const size_t n = strlen(name);
    if (strncmp(argv[i], name, n) != 0) return false;
    if (argv[i][n] == '=') {
        v = argv[i] + n + 1;
        return true;
    }
    if (argv[i][n] != '\0' || i + 1 >= argc) return false;
    v = argv[++i];
    return true;
}

/** As opt_value, for an option with a short alias such as -p. */
static bool opt_value(int argc, char** argv, int& i, const char* name, const char* alias, std::string& v) {
    return opt_value(argc, argv, i, name, v) || opt_value(argc, argv, i, alias, v);
}

/** An unsigned option value, refusing anything that is not a whole number. */
static unsigned long long opt_uint(const char* name, const std::string& v) {
    size_t used = 0;
    unsigned long long x = 0;
    try {
        if (!v.empty() && std::isdigit(static_cast<unsigned char>(v[0]))) x = std::stoull(v, &used);
    } catch (const std::exception&) {
        used = 0;
    }
    if (used == 0 || used != v.size()) throw std::invalid_argument(std::string(name) + " needs a whole number, not '" + v + "'");
    return x;
}

/** A nonnegative decimal option value, refusing anything else. */
static double opt_unonneg(const char* name, const std::string& v) {
    size_t used = 0;
    double x = -1;
    try {
        if (!v.empty() && (std::isdigit(static_cast<unsigned char>(v[0])) || v[0] == '.')) x = std::stod(v, &used);
    } catch (const std::exception&) {
        used = 0;
    }
    if (used == 0 || used != v.size() || !(x >= 0))
        throw std::invalid_argument(std::string(name) + " needs a nonnegative number, not '" + v + "'");
    return x;
}

/**
 * The split parameters for --lf-split or --phi-split FACTOR: Orbit's
 * default length capping and balancing factor FACTOR, or none for 0.
 */
static orbit::split_params split_arg(const char* name, const std::string& v) {
    const ulint b = opt_uint(name, v);
    return b == 0 ? orbit::NO_SPLITTING : orbit::split_params(orbit::DEFAULT_LENGTH_CAPPING, b);
}

/**
 * Parse --lf-split's value into sp.  Returns false, after saying why, for 0:
 * LF is always split.
 */
static bool lf_split_arg(const std::string& v, orbit::split_params& sp) {
    sp = split_arg("--lf-split", v);
    if (sp != orbit::NO_SPLITTING) return true;
    std::cerr << "--lf-split must be positive: LF is always split and balanced\n";
    return false;
}

/** Sets out to in compressed as one gzip member at the given level. */
static void gzip_member(const std::string& in, std::string& out, int level) {
    z_stream z{};
    if (deflateInit2(&z, level, Z_DEFLATED, 15 + 16, 8, Z_DEFAULT_STRATEGY) != Z_OK)
        throw std::runtime_error("deflateInit2 failed");
    out.resize(deflateBound(&z, static_cast<uLong>(in.size())));
    z.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(in.data()));
    z.avail_in = static_cast<uInt>(in.size());
    z.next_out = reinterpret_cast<Bytef*>(&out[0]);
    z.avail_out = static_cast<uInt>(out.size());
    const int r = deflate(&z, Z_FINISH);
    const size_t n = z.total_out;
    deflateEnd(&z);
    if (r != Z_STREAM_END) throw std::runtime_error("deflate failed");
    out.resize(n);
}

// One block's results, owned by one query thread.
struct QueryResults {
    std::vector<std::vector<ulint>> ms;
    std::vector<std::vector<ulint>> pos;  // tms positions, with --positions or an SMEM report
    std::vector<TmsSmemHits> hits;        // tms SMEMs, with an SMEM report
    std::vector<std::vector<RepositionDecision>> decisions;  // reposition outcomes, with --oracle-write
};

// The ftab of a query, from --ftab or --build-ftab, or null.
static const ftab::Ftab* g_ftab = nullptr;
// With --oracle-write, each query records its reposition outcomes; in the
// oracle modes, they are replayed from g_oracle, read from --oracle.
static bool g_record = false;
static RepositionOracle g_oracle;
// The query's --mode on an ms index.
static MsMode g_ms_mode = MsMode::STANDARD;

// first_read is the index of p[0] among all the reads, which says where its
// recorded repositions are in the oracle modes.
static void query_many(MSIndexSpillLCP<false>& idx, const std::vector<std::string>& p, size_t k, QueryResults& r,
                       size_t first_read) {
    RepositionOracleView view{};
    const bool replay = g_ms_mode != MsMode::STANDARD;
    if (replay) view = RepositionOracleView::of(g_oracle, first_read, p.size());
    ms_query_batch(idx, p, k, r.ms, g_ftab, g_ms_mode, replay ? &view : nullptr, g_record ? &r.decisions : nullptr);
}
// Results of one read go to slot 0 of r.  With an ftab or an oracle, the
// batched engine with one read, since ms_query takes neither.
static void query_one(MSIndexSpillLCP<false>& idx, const std::string& s, QueryResults& r, size_t first_read) {
    if (g_ftab || g_record || g_ms_mode != MsMode::STANDARD) return query_many(idx, {s}, 1, r, first_read);
    r.ms.resize(1);
    r.ms[0] = ms_query(idx, s);
}
// A query's settings on a tms index, from its command line.  They are read-only once the
// queries start.
static TmsMode g_tms_mode = TmsMode::PSI;
static bool g_tms_positions = false;
static bool g_tms_phi_starts = false;
static TmsReport g_tms_report = TmsReport::MS;
static ulint g_tms_min_smem_len = 0;
static ulint g_tms_max_smem_positions = TMS_ALL_POSITIONS;
// How smem-all counts occurrences, and what backward search reads, built
// after the index is loaded when the method is not walk.
static TmsCountOptions g_tms_count;
static unsigned g_tms_count_shift = TmsCounter::DEFAULT_SHIFT;
static TmsCounter g_tms_counter;

// With an SMEM report, the SMEMs are found as part of the query.  smem_k
// walks are in flight when smem-all lists positions (0 lists them one SMEM
// at a time).  first_read is the index of p[0] among all the reads, which
// says where its recorded repositions are in the oracle mode.
static void query_many(TmsIndex& idx, const std::vector<std::string>& p, size_t k, QueryResults& r, size_t smem_k,
                       size_t first_read) {
    const bool smems = g_tms_report != TmsReport::MS;
    // Backward search counts and lists smem-all's SMEMs without toeholds.
    const bool toeholds = g_tms_positions || (smems && !(g_tms_report == TmsReport::SMEM_ALL &&
                                                          g_tms_count.method == TmsCount::BACKWARD));
    RepositionOracleView view{};
    const bool replay = g_tms_mode == TmsMode::ORACLE || g_tms_mode == TmsMode::ORACLE_SCAN;
    if (replay) view = RepositionOracleView::of(g_oracle, first_read, p.size());
    tms_query_batch(idx, p, k, r.ms, g_tms_mode, toeholds ? &r.pos : nullptr, true, g_ftab,
                    replay ? &view : nullptr, g_record ? &r.decisions : nullptr);
    if (!toeholds) r.pos.clear();
    if (!smems) return;
    const TmsCounting counting{g_tms_count, &g_tms_counter, &p};
    tms_report_smems_batch(idx, r.ms, r.pos, g_tms_min_smem_len, g_tms_report, smem_k, r.hits,
                           g_tms_max_smem_positions, &counting);
}
static void query_many(TmsIndex& idx, const std::vector<std::string>& p, size_t k, QueryResults& r, size_t first_read) {
    query_many(idx, p, k, r, k, first_read);
}
// One read at a time runs the batched engine with one read: it is faster
// than tms_query, the plain psi reference, whereas ms_query is as fast as
// ms's engine with one read.  smem-all lists positions one SMEM at a time.
static void query_one(TmsIndex& idx, const std::string& s, QueryResults& r, size_t first_read) {
    query_many(idx, {s}, 1, r, 0, first_read);
}

/** Sets m to the tms mode named by name and returns true, or returns false. */
static bool tms_mode_of(const std::string& name, TmsMode& m) {
    if (name == "psi") m = TmsMode::PSI;
    else if (name == "phi") m = TmsMode::PHI;
    else if (name == "phisum") m = TmsMode::PHISUM;
    else if (name == "phiskip") m = TmsMode::PHISKIP;
    else if (name == "dual") m = TmsMode::DUAL;
    else if (name == "oracle") m = TmsMode::ORACLE;
    else if (name == "oracle-scan") m = TmsMode::ORACLE_SCAN;
    else return false;
    return true;
}

/**
 * If argv[i] is a tms query option that chooses what a query loads or how
 * smem-all counts (--positions, --report, --phi-starts, the --count options,
 * and --mode if with_mode), apply it
 * to the tms query settings, advancing i past its value, and return 1;
 * return -1, after printing why, if its value is unknown, and 0 if argv[i]
 * is no such option.
 */
static int tms_query_flag(int argc, char** argv, int& i, bool with_mode = true) {
    std::string v;
    if (strcmp(argv[i], "--positions") == 0) g_tms_positions = true;
    else if (strcmp(argv[i], "--phi-starts") == 0) g_tms_phi_starts = true;
    else if (with_mode && opt_value(argc, argv, i, "--mode", v)) {
        if (!tms_mode_of(v, g_tms_mode)) { std::cerr << "Unknown mode: " << v << "\n"; return -1; }
    }
    else if (opt_value(argc, argv, i, "--report", v)) {
        if (v == "ms") g_tms_report = TmsReport::MS;
        else if (v == "smem-one") g_tms_report = TmsReport::SMEM_ONE;
        else if (v == "smem-all") g_tms_report = TmsReport::SMEM_ALL;
        else { std::cerr << "Unknown report: " << v << " (expected ms, smem-one or smem-all)\n"; return -1; }
    }
    else if (opt_value(argc, argv, i, "--count", v)) {
        if (v == "walk") g_tms_count.method = TmsCount::WALK;
        else if (v == "backward") g_tms_count.method = TmsCount::BACKWARD;
        else if (v == "adaptive") g_tms_count.method = TmsCount::ADAPTIVE;
        else { std::cerr << "Unknown counting method: " << v << " (expected walk, backward or adaptive)\n"; return -1; }
    }
    else if (opt_value(argc, argv, i, "--count-budget", v))
        g_tms_count.budget = static_cast<ulint>(opt_uint("--count-budget", v));
    else if (opt_value(argc, argv, i, "--count-budget-per-len", v))
        g_tms_count.budget_per_len = opt_unonneg("--count-budget-per-len", v);
    else if (opt_value(argc, argv, i, "--count-short-len", v))
        g_tms_count.short_len = static_cast<ulint>(opt_uint("--count-short-len", v));
    else if (opt_value(argc, argv, i, "--count-row-shift", v))
        g_tms_count_shift = static_cast<unsigned>(opt_uint("--count-row-shift", v));
    else return 0;
    return 1;
}

/** The parts a tms query with the current settings loads. */
static TmsParts tms_query_load() {
    TmsParts need = tms_query_parts(g_tms_mode, g_tms_positions, g_tms_report);
    need.phi_starts = need.phi_starts || (need.phi && g_tms_phi_starts);
    return need;
}

// The command line, recorded in the provenance of the indexes it builds.
static std::string g_command_line;

/** Build provenance for an index: this command line and the time, as "key=value" lines. */
static std::string build_info() {
    const std::time_t t = std::time(nullptr);
    char when[32];
    std::strftime(when, sizeof(when), "%Y-%m-%dT%H:%M:%SZ", std::gmtime(&t));
    return "built_by=" + g_command_line + "\nbuilt_at=" + when + "\n";
}

// Read a tms index with the parts in *need, or every part it holds if need is null.
static std::optional<TmsIndex> read_tms_index(const std::string& path, const TmsParts* need = nullptr) {
    std::ifstream in(path, std::ios::binary);
    if (!in.good()) return std::nullopt;
    TmsIndex idx;
    try {
        if (need) idx.load(in, *need);
        else idx.load(in);
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
        return std::nullopt;
    }
    return idx;
}

// A query's command line, parsed before the index's header says which
// kind of index it is.  The tms settings that choose what a query loads go
// to the g_tms_* globals as they are parsed.
struct QueryOptions {
    std::string idx_path, reads_path, out_path;
    bool write_output = true;
    size_t interleave = 32, threads = 0, block_reads = 0, block_bytes = 0;
    int gzip_level = -2;  // -2: from the output name; -1: none
    ReadOrient orient = ReadOrient::FWD;
    std::string mode;  // empty: the index kind's default
    std::string oracle_path, oracle_write_path;
    unsigned ftab_k = 0;    // with --build-ftab, the ftab's k; 0: none built
    std::string ftab_path;  // a stored ftab to load
    std::vector<std::string> tms_only;  // options given that need a tms index
};

/** Parse query's arguments after the command name.  Returns false, after saying why, on a bad one. */
static bool parse_query(int argc, char** argv, QueryOptions& o) {
    if (argc < 2) {
        std::cerr << "query requires INDEX and READS\n";
        return false;
    }
    o.idx_path = argv[0];
    o.reads_path = argv[1];
    std::string v;
    for (int i = 2; i < argc; ++i) {
        const std::string arg = argv[i];
        int q = 0;
        if (opt_value(argc, argv, i, "-o", v)) o.out_path = v;
        else if (opt_value(argc, argv, i, "--ftab", v)) o.ftab_path = v;
        else if (arg == "--build-ftab") o.ftab_k = 14;
        else if (opt_value(argc, argv, i, "--build-ftab", v)) {
            o.ftab_k = static_cast<unsigned>(opt_uint("--build-ftab", v));
            if (o.ftab_k < 1 || o.ftab_k > ftab::max_k) { std::cerr << "--build-ftab's KLEN must be between 1 and 31\n"; return false; }
        }
        else if (arg == "--no-output") o.write_output = false;
        else if (opt_value(argc, argv, i, "--interleave", v)) o.interleave = opt_uint("--interleave", v);
        else if (opt_value(argc, argv, i, "--threads", "-p", v)) o.threads = opt_uint("--threads", v);
        else if (opt_value(argc, argv, i, "--block-reads", v)) o.block_reads = opt_uint("--block-reads", v);
        else if (opt_value(argc, argv, i, "--block-bytes", v)) o.block_bytes = opt_uint("--block-bytes", v);
        else if (opt_value(argc, argv, i, "--gzip-level", v)) o.gzip_level = static_cast<int>(opt_uint("--gzip-level", v));
        else if (opt_value(argc, argv, i, "--orient", v)) {
            if (v == "fwd") o.orient = ReadOrient::FWD;
            else if (v == "rc") o.orient = ReadOrient::RC;
            else if (v == "qual") o.orient = ReadOrient::QUAL;
            else { std::cerr << "Unknown orientation: " << v << " (expected fwd, rc or qual)\n"; return false; }
        }
        else if (opt_value(argc, argv, i, "--mode", v)) o.mode = v;
        else if ((q = tms_query_flag(argc, argv, i, false)) != 0) {
            if (q < 0) return false;
            // --report ms asks for what every query gives.
            if (arg.rfind("--report", 0) != 0 || g_tms_report != TmsReport::MS) o.tms_only.push_back(arg);
        }
        else if (opt_value(argc, argv, i, "--min-smem-len", v)) {
            // ropebwt3's -l L keeps SMEMs of length L or more; the SMEM
            // report keeps those longer than its threshold.
            const ulint l = opt_uint("--min-smem-len", v);
            g_tms_min_smem_len = l > 0 ? l - 1 : 0;
            o.tms_only.push_back("--min-smem-len");
        }
        else if (opt_value(argc, argv, i, "--max-smem-positions", v)) {
            g_tms_max_smem_positions = static_cast<ulint>(opt_uint("--max-smem-positions", v));
            o.tms_only.push_back("--max-smem-positions");
        }
        else if (opt_value(argc, argv, i, "--oracle", v)) o.oracle_path = v;
        else if (opt_value(argc, argv, i, "--oracle-write", v)) o.oracle_write_path = v;
        else { std::cerr << "Unknown query option: " << arg << " (see ms --help)\n"; return false; }
    }
    return true;
}

/**
 * query: load the index once, then compute matching statistics for every
 * read, on one or more threads.  Each thread queries a block of reads at a
 * time with its own DEPTH reads in flight, and output is written in input
 * order, so it does not depend on the thread count.
 */
template <typename Index>
static int run_batch(const QueryOptions& o, std::optional<Index> (*read)(const std::string&)) {
    const std::string &idx_path = o.idx_path, &reads_path = o.reads_path, &out_path = o.out_path;
    const bool write_output = o.write_output;
    size_t interleave = o.interleave, threads = o.threads, block_reads = o.block_reads, block_bytes = o.block_bytes;
    int gzip_level = o.gzip_level;
    const ReadOrient orient = o.orient;
    const std::string &oracle_path = o.oracle_path, &oracle_write_path = o.oracle_write_path;
    unsigned ftab_k = o.ftab_k;
    const std::string& ftab_path = o.ftab_path;
    // Positions found on the reverse complement would have to be mirrored
    // onto the other strand, which needs each string's offset in the text,
    // and the index does not store those.
    if (orient != ReadOrient::FWD && (g_tms_positions || g_tms_report != TmsReport::MS)) {
        std::cerr << "--orient rc and --orient qual give matching statistics only; they cannot be used with "
                     "--positions or an SMEM report\n";
        return 1;
    }
    // smem-all writes one line per SMEM, so there is no line per read to
    // hold the positions of --positions.
    if (g_tms_report == TmsReport::SMEM_ALL && g_tms_positions) {
        std::cerr << "--positions cannot be used with --report smem-all\n";
        return 1;
    }
    if (ftab_k > 0 && !ftab_path.empty()) { std::cerr << "--ftab cannot be used with --build-ftab\n"; return 1; }
    // A stored table is checked against the index's header before either is
    // loaded.  A tms query that walks phi or reports positions needs its
    // toehold, and positions its window positions.
    ftab::FileHeader ftab_header;
    if (!ftab_path.empty()) {
        try {
            ftab_header = ftab::read_header_file(ftab_path);
            const std::string bad = ftab::link_mismatch(ftab_header, orbit_idx::read_header_file(idx_path));
            if (!bad.empty()) {
                std::cerr << ftab_path << " was not built from " << idx_path << ": " << bad << "\n";
                return 1;
            }
        } catch (const std::exception& e) {
            std::cerr << e.what() << "\n";
            return 1;
        }
        const bool positions = g_tms_positions || g_tms_report != TmsReport::MS;
        if (std::is_same_v<Index, TmsIndex> && (g_tms_mode != TmsMode::PSI || positions) && !ftab_header.shape.toehold) {
            std::cerr << ftab_path << " holds no toehold (it was built from a psi layout), which this query needs\n";
            return 1;
        }
        if (std::is_same_v<Index, TmsIndex> && positions && !ftab_header.shape.window) {
            std::cerr << ftab_path << " holds no window positions, which positions and SMEM reports need\n";
            return 1;
        }
    }
    const bool oracle_mode = std::is_same_v<Index, TmsIndex>
                                 ? g_tms_mode == TmsMode::ORACLE || g_tms_mode == TmsMode::ORACLE_SCAN
                                 : g_ms_mode != MsMode::STANDARD;
    if (oracle_mode != !oracle_path.empty()) {
        std::cerr << "--mode oracle or oracle-scan and --oracle FILE go together\n";
        return 1;
    }
    if (oracle_mode && !oracle_write_path.empty()) {
        std::cerr << "--oracle-write records with another mode; the oracle modes replay\n";
        return 1;
    }
    if ((oracle_mode || !oracle_write_path.empty()) && (ftab_k > 0 || !ftab_path.empty())) {
        std::cerr << "--oracle and --oracle-write take no ftab, which skips some repositions\n";
        return 1;
    }
    g_record = !oracle_write_path.empty();
    // The oracle's reads are numbered across the whole input, so blocks are
    // parsed as they are read, in order, to number them.
    const bool number_reads = oracle_mode || g_record;
    threads = resolve_thread_count(threads);
    // Output named *.gz is gzip-compressed, each block as a gzip member by
    // the thread that formats it; members written one after another form one
    // gzip file.
    if (gzip_level == -2)
        gzip_level = out_path.size() > 3 && out_path.compare(out_path.size() - 3, 3, ".gz") == 0 ? 1 : -1;
    if (gzip_level > 9) gzip_level = 9;
    using clock = std::chrono::steady_clock;
    auto t0 = clock::now();
    auto opt = [&] {
        InterleavedAllocation spread;
        return read(idx_path);
    }();
    if (!opt) {
        std::cerr << "Failed to load index: " << idx_path << "\n";
        return 1;
    }
    Index& idx = *opt;
    // smem-all names contigs with the index's sequence-ID table, if it has one.
    if constexpr (std::is_same_v<Index, TmsIndex>) {
        if (g_tms_report == TmsReport::SMEM_ALL) {
            if (const SeqIdTable* t = idx.seq_ids())
                std::cerr << "smem-all: positions as contig:strand:offset, " << t->names.size() << " sequences\n";
            else
                std::cerr << "smem-all: positions as text positions; the index holds no sequence-ID table\n";
        }
    }
    // The counter is part of preparing the index, so it is timed with the load.
    double count_init_s = 0.0;
    if constexpr (std::is_same_v<Index, TmsIndex>) {
        if (g_tms_report == TmsReport::SMEM_ALL && g_tms_count.method != TmsCount::WALK) {
            const auto tc = clock::now();
            try {
                InterleavedAllocation spread;
                g_tms_counter = TmsCounter(idx, g_tms_count_shift);
            } catch (const std::exception& e) {
                std::cerr << "count: " << e.what() << "\n";
                return 1;
            }
            count_init_s = std::chrono::duration<double>(clock::now() - tc).count();
            std::cerr << "count: backward search table " << g_tms_counter.bytes() << " bytes, row shift "
                      << g_tms_counter.rows().shift() << ", built in " << count_init_s << " s\n";
        }
    }
    if (oracle_mode) {
        std::ifstream oin(oracle_path, std::ios::binary);
        if (!oin.good()) { std::cerr << "Failed to open oracle: " << oracle_path << "\n"; return 1; }
        try {
            g_oracle.read(oin);
        } catch (const std::exception& e) {
            std::cerr << oracle_path << ": " << e.what() << "\n";
            return 1;
        }
    }
    std::ofstream oracle_out;
    if (g_record) {
        oracle_out.open(oracle_write_path, std::ios::binary);
        if (!oracle_out.good()) { std::cerr << "Failed to open oracle output: " << oracle_write_path << "\n"; return 1; }
    }
    const double load_s = std::chrono::duration<double>(clock::now() - t0).count();
    // The ftab is loaded, or built on the query's threads, after the index, and is
    // interleaved over the memory nodes as the index is.
    std::optional<ftab::Ftab> tab;
    double ftab_s = 0.0;
    if (ftab_k > 0 || !ftab_path.empty()) {
        const auto tf = clock::now();
        InterleavedAllocation spread;
        try {
            if (!ftab_path.empty()) {
                tab = ftab::read_file(ftab_path, ftab_header);
                ftab_k = tab->k();
            } else {
                if constexpr (std::is_same_v<Index, TmsIndex>) tab = ftab::build_tms(idx, ftab_k, threads).tab;
                else tab = ftab::build_ms(idx, ftab_k, threads).tab;
            }
        } catch (const std::exception& e) {
            std::cerr << "ftab: " << e.what() << "\n";
            return 1;
        }
        ftab_s = std::chrono::duration<double>(clock::now() - tf).count();
        g_ftab = &*tab;
        std::cerr << "ftab: k=" << ftab_k << " entry_bits=" << tab->entry_bits() << " bytes=" << tab->bytes()
                  << (ftab_path.empty() ? " build_s=" : " load_s=") << ftab_s << "\n";
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

    // Query threads take turns reading a block of bytes into their own
    // buffers and cutting it after its last whole record.  Blocks are sized
    // in bytes, from the mean record size seen so far, to hold about
    // block_reads reads (or block_bytes bytes, if given), and at most
    // max_block_bytes unless one record is longer, so that long reads do not
    // make a block's per-base results large.  The default keeps K reads in
    // flight for all but the end of a block, and keeps blocks small enough
    // to spread over many threads.
    if (block_reads == 0) block_reads = std::max<size_t>(512, 16 * interleave);
    constexpr size_t max_block_bytes = size_t(1) << 24;
    ReadCutter cutter(*in);
    std::atomic<size_t> seen_records{0}, seen_bytes{0};
    auto target_bytes = [&]() -> size_t {
        if (block_bytes > 0) return block_bytes;
        const size_t r = seen_records.load(std::memory_order_relaxed), b = seen_bytes.load(std::memory_order_relaxed);
        const size_t per = r > 0 ? std::max<size_t>(b / r, 1) : 160;
        return std::min(std::max<size_t>(block_reads * per, size_t(1) << 16), max_block_bytes);
    };

    // A block's records are found by one thread at a time and parsed by the
    // thread that queries them.
    struct Block {
        ReadBuffer buf;
        ReadChunk chunk;
        std::vector<std::string> names, seqs, quals;
        bool parsed = false;    // names, seqs and quals hold the chunk's reads
        size_t first_read = 0;  // with number_reads, the index of the block's first read
    };
    // Parse a block's reads, timed into *parse_s.
    auto parse_block = [&](Block& b, double* parse_s) {
        auto tp = clock::now();
        if (orient == ReadOrient::QUAL && b.chunk.format != ReadFormat::FASTQ)
            throw std::runtime_error("--orient qual needs FASTQ input, for its base qualities");
        const size_t n =
            parse_chunk(b.chunk, b.names, b.seqs, orient == ReadOrient::QUAL ? &b.quals : nullptr);
        b.names.resize(n);
        b.seqs.resize(n);
        b.parsed = true;
        *parse_s += std::chrono::duration<double>(clock::now() - tp).count();
    };
    // With --oracle-write, each block's recorded repositions, by its first
    // read, until they are written in order at the end.
    std::mutex record_mu;
    std::map<size_t, RepositionOracle> recorded;
    size_t reads_numbered = 0;
    double numbering_parse_s = 0.0;
    // With --positions, a tab-separated field holds the positions (-1 for
    // none); with smem-one, a last one holds the SMEMs.  smem-all writes
    // ropebwt3 mem's lines, one per SMEM, in place of the read's line.
    auto format = [&](std::string& text, const std::string& nm, const QueryResults& r, size_t j) {
        if constexpr (std::is_same_v<Index, TmsIndex>) {
            // ropebwt3 marks the read's unmatched characters only when its
            // length cut keeps every SMEM.
            if (g_tms_report == TmsReport::SMEM_ALL)
                return tms_format_smems_rope(nm, r.hits[j], idx.seq_ids(), text,
                                             g_tms_min_smem_len == 0 ? &r.ms[j] : nullptr);
        }
        const auto& ms = r.ms[j];
        text += nm;
        text += '\t';
        for (size_t i = 0; i < ms.size(); ++i) {
            if (i > 0) text += ' ';
            append_uint(text, ms[i]);
        }
        if constexpr (std::is_same_v<Index, TmsIndex>) {
            if (g_tms_positions) {
                text += '\t';
                const auto& pos = r.pos[j];
                for (size_t i = 0; i < pos.size(); ++i) {
                    if (i > 0) text += ' ';
                    if (pos[i] == TMS_NO_POS) text += "-1";
                    else append_uint(text, pos[i]);
                }
            }
            if (g_tms_report != TmsReport::MS) {
                text += '\t';
                tms_format_smems(r.hits[j], g_tms_report, text);
            }
        }
        text += '\n';
    };
    struct Worker {
        QueryResults r, one;
        PerfCounters perf;
        double query_s = 0.0, parse_s = 0.0, format_s = 0.0, orient_s = 0.0;
        size_t reads = 0, bases = 0, flipped = 0;
        std::string packed;
        std::vector<char> flips;
        std::vector<ulint> scratch;
        std::function<void(Worker&, Block&, std::string&)> body;
        void process(Block& b, std::string& text) { body(*this, b, text); }
    };
    auto body = [&](Worker& w, Block& b, std::string& text) {
        if (!b.parsed) parse_block(b, &w.parse_s);
        b.parsed = false;
        const size_t n = b.seqs.size();
        w.reads += n;
        for (const auto& sq : b.seqs) w.bases += sq.size();
        seen_records.fetch_add(n, std::memory_order_relaxed);
        seen_bytes.fetch_add(b.chunk.size, std::memory_order_relaxed);
        auto timed = [&](auto&& f) {
            auto tq = clock::now();
            w.perf.start();
            f();
            w.perf.stop();
            w.query_s += std::chrono::duration<double>(clock::now() - tq).count();
        };
        // Reads to be matched on the other strand are reverse-complemented
        // in place, and their results converted back after the query.
        const bool oriented = orient != ReadOrient::FWD;
        if (oriented) {
            auto to = clock::now();
            w.flipped += orient_reads(b.seqs, &b.quals, orient, w.flips);
            w.orient_s += std::chrono::duration<double>(clock::now() - to).count();
        }
        if (interleave == 0) {
            // One read at a time, without prefetching.  The block is timed
            // as a whole, so that timing each read does not add to it.
            const size_t n_reads = b.seqs.size();
            timed([&] {
                w.r.ms.resize(n_reads);
                if (g_record) w.r.decisions.resize(n_reads);
                for (size_t j = 0; j < n_reads; ++j) {
                    query_one(idx, b.seqs[j], w.one, b.first_read + j);
                    w.r.ms[j] = std::move(w.one.ms[0]);
                    if (g_record) w.r.decisions[j] = std::move(w.one.decisions[0]);
                    if (!w.one.pos.empty()) {
                        w.r.pos.resize(n_reads);
                        w.r.pos[j] = std::move(w.one.pos[0]);
                    }
                    if (!w.one.hits.empty()) {
                        w.r.hits.resize(n_reads);
                        w.r.hits[j] = std::move(w.one.hits[0]);
                    }
                }
            });
        } else {
            timed([&] { query_many(idx, b.seqs, interleave, w.r, b.first_read); });
        }
        if (g_record) {
            RepositionOracle part;
            part.append(w.r.decisions);
            std::lock_guard<std::mutex> lk(record_mu);
            recorded.emplace(b.first_read, std::move(part));
        }
        if (oriented) {
            auto to = clock::now();
            restore_orientation(w.flips, w.r.ms, w.scratch);
            w.orient_s += std::chrono::duration<double>(clock::now() - to).count();
        }
        if (write_output) {
            auto tf = clock::now();
            for (size_t j = 0; j < b.seqs.size(); ++j) format(text, b.names[j], w.r, j);
            if (gzip_level >= 0) {
                gzip_member(text, w.packed, gzip_level);
                text.swap(w.packed);
            }
            w.format_s += std::chrono::duration<double>(clock::now() - tf).count();
        }
    };
    PerfCounters perf;
    double query_s = 0.0, read_s = 0.0, parse_s = 0.0, format_s = 0.0, orient_s = 0.0;
    size_t n_reads = 0, n_bases = 0, n_flipped = 0;
    TmsCountStats count_total;
#ifdef TMS_STATS
    TmsStats tms_total;
#endif
#ifdef MS_STATS
    MsStats ms_total;
#endif
#ifdef MS_STEP_STATS
    StepAgg step_total;
    FFAgg ff_total;
#endif
    auto t_all = clock::now();
    try {
        run_parallel_blocks<Block>(
            threads, 4 * threads,
            [&](size_t) { Worker w; w.body = body; return w; },
            [&](Block& b) {
                auto tr = clock::now();
                const bool got = cutter.next(b.buf, b.chunk, target_bytes());
                read_s += std::chrono::duration<double>(clock::now() - tr).count();
                // Called under the input lock, in input order.
                if (got && number_reads) {
                    parse_block(b, &numbering_parse_s);
                    b.first_read = reads_numbered;
                    reads_numbered += b.seqs.size();
                }
                return got;
            },
            [&](const std::string& text) {
                if (write_output) out->write(text.data(), static_cast<std::streamsize>(text.size()));
            },
            [&](Worker& w) {
                perf.add(w.perf);
                query_s += w.query_s;
                parse_s += w.parse_s;
                format_s += w.format_s;
                orient_s += w.orient_s;
                n_reads += w.reads;
                n_bases += w.bases;
                n_flipped += w.flipped;
                count_total += tms_count_stats;
#ifdef TMS_STATS
                tms_total += tms_stats;
#endif
#ifdef MS_STATS
                ms_total += ms_stats;
#endif
#ifdef MS_STEP_STATS
                step_total += step_agg;
                ff_total += ff_agg;
#endif
            });
    } catch (const std::exception& e) {
        std::cerr << "query: " << e.what() << "\n";
        return 1;
    }
    out->flush();
    const double total_s = std::chrono::duration<double>(clock::now() - t_all).count();
    parse_s += numbering_parse_s;
    if (oracle_mode && reads_numbered != g_oracle.patterns()) {
        std::cerr << "query: the oracle holds " << g_oracle.patterns() << " reads but the input has "
                  << reads_numbered << "\n";
        return 1;
    }
    if (g_record) {
        RepositionOracle all;
        all.move_runs = idx.move_runs();
        for (const auto& [first, part] : recorded) {
            for (size_t j = 0; j < part.patterns(); ++j)
                all.offsets.push_back(all.offsets.back() + part.offsets[j + 1] - part.offsets[j]);
            all.decisions.insert(all.decisions.end(), part.decisions.begin(), part.decisions.end());
        }
        try {
            all.write(oracle_out);
            oracle_out.close();
        } catch (const std::exception& e) {
            std::cerr << "query: " << e.what() << "\n";
            return 1;
        }
        std::cerr << "oracle: reads=" << all.patterns() << " repositions=" << all.decisions.size()
                  << " bytes=" << 32 + 8 * all.patterns() + 16 * all.decisions.size() << "\n";
    }
#ifdef TMS_STATS
    {
        const TmsStats& tms_stats = tms_total;
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
        if (g_tms_mode == TmsMode::PHISUM)
            std::cerr << "stats: phisum up_settled=" << tms_stats.sum_up_settled << " up_walked=" << tms_stats.sum_up_walked
                      << " down_settled=" << tms_stats.sum_down_settled << " down_walked=" << tms_stats.sum_down_walked
                      << " partial_steps=" << tms_stats.sum_partial_steps << " cont_walks=" << tms_stats.sum_cont_walks
                      << " cont_steps=" << tms_stats.sum_cont_steps
                      << " partial_steps/rep=" << double(tms_stats.sum_partial_steps) / tms_stats.repositions
                      << " cont_steps/rep=" << double(tms_stats.sum_cont_steps) / tms_stats.repositions << "\n";
    }
#endif
#ifdef MS_STATS
    std::cerr << "stats: bases=" << n_bases << " repositions/base=" << double(ms_total.repositions) / n_bases
              << " lf_ff/step=" << double(ms_total.lf_ff) / ms_total.lf_steps
              << " walk_up/rep=" << double(ms_total.walk_up) / ms_total.repositions
              << " walk_down/rep=" << double(ms_total.walk_down) / ms_total.repositions << "\n";
#endif
#ifdef MS_STEP_STATS
    if (interleave == 0) step_total.print(std::cerr);
    else std::cerr << "stepstats: none recorded; steps are recorded with --interleave 0\n";
    ff_total.print(std::cerr);
#endif
    if (std::is_same_v<Index, TmsIndex> && g_tms_report == TmsReport::SMEM_ALL) {
        // Counting's share of query_s (walk_s, search_s and list_s sum the
        // threads' time in each phase) and what it did.
        const TmsCountStats& c = count_total;
        const char* method[] = {"walk", "backward", "adaptive"};
        std::cerr << "count: method=" << method[static_cast<int>(g_tms_count.method)]
                  << " budget=" << g_tms_count.budget << " budget_per_len=" << g_tms_count.budget_per_len
                  << " short_len=" << g_tms_count.short_len
                  << " smems=" << c.smems << " walked=" << c.walked << " switched=" << c.switched
                  << " searched=" << c.searched << " walk_steps=" << c.walk_steps
                  << " search_steps=" << c.search_steps << " top_moves=" << c.top_moves
                  << " bot_moves=" << c.bot_moves << " scan_rows=" << c.scan_rows << " lf_ff=" << c.lf_ff
                  << " far_ends=" << c.far_ends << " sum_rows=" << c.sum_rows << " list_steps=" << c.list_steps
                  << " walk_s=" << c.walk_s << " search_s=" << c.search_s << " list_s=" << c.list_s
                  << " count_init_s=" << count_init_s << "\n";
    }
    perf.report(std::cerr, n_bases);
    // query_s sums the threads' time in queries, parse_s their time
    // parsing reads, and format_s their time formatting (and compressing)
    // output.  read_s is the time the input lock was held for reading; wall_ns_per_base is
    // the elapsed time of the whole read, query and write phase per base.
    std::cerr << "query: reads=" << n_reads << " bases=" << n_bases << " threads=" << threads
              << " index_load_s=" << load_s << " ftab_k=" << ftab_k << " ftab_s=" << ftab_s << " query_s=" << query_s << " read_s=" << read_s
              << " parse_s=" << parse_s << " format_s=" << format_s
              << " orient=" << (orient == ReadOrient::FWD ? "fwd" : orient == ReadOrient::RC ? "rc" : "qual")
              << " flipped=" << n_flipped << " orient_s=" << orient_s
              << " total_s=" << total_s
              << " query_ns_per_base=" << (n_bases ? query_s * 1e9 / n_bases : 0.0)
              << " wall_ns_per_base=" << (n_bases ? total_s * 1e9 / n_bases : 0.0) << "\n";
    return 0;
}

/**
 * query: read the index's header, whose tool field says whether it is an ms
 * or a tms index, then run that kind's engine, refusing options only a tms
 * index serves.
 */
static int run_query(int argc, char** argv) {
    QueryOptions o;
    if (!parse_query(argc, argv, o)) return 1;
    orbit_idx::Header h;
    try {
        h = orbit_idx::read_header_file(o.idx_path);
    } catch (const std::exception& e) {
        std::cerr << "query: " << e.what() << "\n";
        return 1;
    }
    if (h.tool == orbit_idx::Tool::MS) {
        if (!o.tms_only.empty()) {
            std::cerr << "query: " << o.tms_only[0] << " needs a tms index, and " << o.idx_path << " is an ms index\n";
            return 1;
        }
        if (o.mode.empty() || o.mode == "standard") g_ms_mode = MsMode::STANDARD;
        else if (o.mode == "oracle") g_ms_mode = MsMode::ORACLE;
        else if (o.mode == "oracle-scan") g_ms_mode = MsMode::ORACLE_SCAN;
        else if (TmsMode m; tms_mode_of(o.mode, m)) {
            std::cerr << "query: --mode " << o.mode << " needs a tms index, and " << o.idx_path
                      << " is an ms index, whose only mode is standard\n";
            return 1;
        } else {
            std::cerr << "Unknown mode: " << o.mode << "\n";
            return 1;
        }
        return run_batch<MSIndexSpillLCP<false>>(o, ms_serialize::read_index);
    }
    if (o.mode == "standard") {
        std::cerr << "query: --mode standard is for ms indexes; " << o.idx_path
                  << " is a tms index, whose modes are psi (the default), phi, phiskip and dual\n";
        return 1;
    }
    if (!o.mode.empty() && !tms_mode_of(o.mode, g_tms_mode)) {
        std::cerr << "Unknown mode: " << o.mode << "\n";
        return 1;
    }
    // Load only the parts the query reads; the load refuses an index that
    // lacks one.
    return run_batch<TmsIndex>(o, [](const std::string& path) {
        const TmsParts need = tms_query_load();
        return read_tms_index(path, &need);
    });
}

/** Run the command named by argv[1]. */
static int run_command(int argc, char** argv) {
    if (argc < 2) {
        usage(std::cout);
        return 0;
    }
    const std::string cmd = argv[1];
    if (cmd == "--help" || cmd == "-h" || cmd == "help") {
        usage(std::cout);
        return 0;
    }
    if (cmd == "--help-advanced") {
        usage_advanced(std::cout);
        return 0;
    }
    for (int i = 0; i < argc; ++i) g_command_line += (i ? " " : "") + std::string(argv[i]);
    argc -= 2;
    argv += 2;

#ifdef NDEBUG
    // The test suites check their results with assert.
    if (cmd == "test" || cmd == "tms-test" || cmd == "fmd-test") {
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

    if (cmd == "pattern") {
        if (argc < 2) {
            std::cerr << "pattern requires INDEX and PATTERN\n";
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
        if (argc < 3) {
            std::cerr << "build-rlbwt requires HEADS, LENS and INDEX\n";
            return 1;
        }
        const std::string heads_path = argv[0], lens_path = argv[1], idx_path = argv[2];
        std::string minima_path, lcp_bin_path;
        double percentile_k = 0.98;
        bool coalesce = false, minima_only = false;
        ulint spill_align = 0;
        uchar spill_split_bits = 0;
        orbit::split_params lf_split;
        std::string v;
        for (int i = 3; i < argc; ++i) {
            if (opt_value(argc, argv, i, "--minima", v)) minima_path = v;
            else if (opt_value(argc, argv, i, "--lcp-bin", v)) lcp_bin_path = v;
            else if (strcmp(argv[i], "--minima-only") == 0) minima_only = true;
            else if (opt_value(argc, argv, i, "--percentile", v)) percentile_k = std::stod(v);
            else if (strcmp(argv[i], "--coalesce-spillover") == 0) coalesce = true;
            else if (opt_value(argc, argv, i, "--spill-align", v)) spill_align = opt_uint("--spill-align", v);
            else if (opt_value(argc, argv, i, "--spill-split-bits", v))
                spill_split_bits = static_cast<uchar>(opt_uint("--spill-split-bits", v));
            else if (opt_value(argc, argv, i, "--lf-split", v)) {
                if (!lf_split_arg(v, lf_split)) return 1;
            }
            else { std::cerr << "Unknown build-rlbwt option: " << argv[i] << "\n"; return 1; }
        }
        if (minima_path.empty() == lcp_bin_path.empty()) {
            std::cerr << "build-rlbwt needs exactly one of --minima and --lcp-bin\n";
            return 1;
        }
        // A minima file holds only the minima, so it is minima-only by nature.
        if (!minima_path.empty()) minima_only = true;
        using clock = std::chrono::steady_clock;
        auto t0 = clock::now();
        std::vector<uchar> heads;
        std::vector<ulint> lens;
        std::string err;
        if (!rlbwt_io::load_rlbwt(heads_path, lens_path, heads, lens, err)) {
            std::cerr << "Failed to load RLBWT: " << err << "\n";
            return 1;
        }
        ulint n = 0;
        for (ulint l : lens) n += l;
        mem_report::phase("load_rlbwt", heads.size());
        const size_t runs = heads.size();
        const double load_s = std::chrono::duration<double>(clock::now() - t0).count();
        auto t1 = clock::now();
        ms_build::Options bo;
        bo.percentile_k = percentile_k;
        bo.coalesce = coalesce;
        bo.spill_align = spill_align;
        bo.spill_split_bits = spill_split_bits;
        bo.minima_only = minima_only;
        bo.lf_split = lf_split;
        bo.report = true;
        ms_build::Stats st;
        std::optional<MSIndexSpillLCP<false>> idx;
        try {
            if (!minima_path.empty()) {
                idx = ms_build::build(std::move(heads), std::move(lens), ms_build::minima_source(minima_path, runs, n),
                                      bo, &st);
            } else {
                // The LCP file is read by run, so its source keeps the run lengths.
                const std::vector<ulint> run_lens = lens;
                idx = ms_build::build(std::move(heads), std::move(lens),
                                      ms_build::lcp_bin_source(lcp_bin_path, run_lens, minima_only), bo, &st);
            }
        } catch (const std::exception& e) {
            std::cerr << "Failed to build: " << e.what() << "\n";
            return 1;
        }
        const double build_s = std::chrono::duration<double>(clock::now() - t1).count();
        idx->set_info(build_info());
        if (!ms_serialize::write_index(idx_path, *idx)) {
            std::cerr << "Failed to write index: " << idx_path << "\n";
            return 1;
        }
        mem_report::phase("write", runs);
        std::cout << "Built index: " << idx_path << " (runs=" << runs << ", rows=" << idx->move_runs() << ", n=" << n
                  << ", skinny=" << st.skinny << ", jumbo=" << st.jumbo << ", load_s=" << load_s
                  << ", build_s=" << build_s << ")\n";
        return 0;
    }

    if (cmd == "query") return run_query(argc, argv);

    if (cmd == "tms-build") {
        const int npos = 3;
        if (argc < npos) {
            std::cerr << "tms-build requires HEADS, LENS and INDEX\n";
            return 1;
        }
        TmsBuildOptions opts;
        opts.report = true;
        bool with_phi = false, lcp_summary = false;
        std::string minima_path, lcp_bin_path, phi_path, seq_ids_path, v;
        for (int i = npos; i < argc; ++i) {
            if (opt_value(argc, argv, i, "--lf-split", v)) {
                if (!lf_split_arg(v, opts.lf_split)) return 1;
            }
            else if (opt_value(argc, argv, i, "--phi-split", v)) opts.phi_split = split_arg("--phi-split", v);
            else if (strcmp(argv[i], "--no-phi-inv") == 0) opts.phi_inv = false;
            else if (opt_value(argc, argv, i, "--layout", v)) {
                if (v == "full") opts.layout = TmsLayout::FULL;
                else if (v == "psi") opts.layout = TmsLayout::PSI;
                else if (v == "phi") opts.layout = TmsLayout::PHI;
                else { std::cerr << "Unknown layout: " << v << " (expected full, psi or phi)\n"; return 1; }
            }
            else if (opt_value(argc, argv, i, "--minima", v)) { minima_path = v; with_phi = true; }
            else if (opt_value(argc, argv, i, "--lcp-bin", v)) { lcp_bin_path = v; with_phi = true; }
            else if (opt_value(argc, argv, i, "--phi", v)) { phi_path = v; with_phi = true; }
            else if (opt_value(argc, argv, i, "--seq-ids", v)) seq_ids_path = v;
            else if (strcmp(argv[i], "--lcp-summary") == 0) lcp_summary = true;
            else { std::cerr << "Unknown " << cmd << " option: " << argv[i] << "\n"; return 1; }
        }
        if (lcp_summary && minima_path.empty() && lcp_bin_path.empty()) {
            std::cerr << "--lcp-summary needs --minima or --lcp-bin\n";
            return 1;
        }
        if (int(!minima_path.empty()) + int(!lcp_bin_path.empty()) + int(!phi_path.empty()) > 1) {
            std::cerr << "give at most one of --phi, --minima and --lcp-bin\n";
            return 1;
        }
        if (opts.layout == TmsLayout::PSI && with_phi) {
            std::cerr << "--layout psi takes no LCP input (--phi, --minima or --lcp-bin)\n";
            return 1;
        }
        if (opts.layout == TmsLayout::PHI && !with_phi) {
            std::cerr << "--layout phi needs LCP input (--phi, --minima or --lcp-bin)\n";
            return 1;
        }
        // The sequence-ID table is read first, so that a bad one fails
        // before the build.
        std::optional<SeqIdTable> seq_ids;
        if (!seq_ids_path.empty()) {
            try {
                seq_ids = SeqIdTable::read(seq_ids_path);
            } catch (const std::exception& e) {
                std::cerr << e.what() << "\n";
                return 1;
            }
        }
        using clock = std::chrono::steady_clock;
        auto t0 = clock::now();
        std::vector<uchar> heads;
        std::vector<ulint> lens, tops;
        TmsPhiSamples samples;
        std::string err;
        if (!rlbwt_io::load_rlbwt(argv[0], argv[1], heads, lens, err)) {
            std::cerr << "Failed to load RLBWT: " << err << "\n";
            return 1;
        }
        if (!phi_path.empty()) {
            ulint n = 0;
            for (ulint l : lens) n += l;
            if (!phi_samples_io::read_phi_samples(phi_path, n, TmsPhiSamples::records_for(heads, lens), samples, err)) {
                std::cerr << "Failed to load phi samples: " << err << "\n";
                return 1;
            }
        } else if (with_phi && !lcp_bin_path.empty()) {
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
            // Only each run's top is kept, so the file is read a run at a time.
            ulint n = 0;
            for (ulint l : lens) n += l;
            rlbwt_io::MinimaReader rd;
            RunLcpPairs p;
            bool ok = rd.open(minima_path, heads.size(), n, err);
            tops.reserve(heads.size());
            for (size_t i = 0; ok && i < heads.size(); ++i) {
                ok = rd.next(p, err);
                if (ok) tops.push_back(p[0].second);
            }
            if (!ok || !rd.finish(err)) {
                std::cerr << "Failed to load minima: " << err << "\n";
                return 1;
            }
        }
        // The LCP summaries, from two more reads of the LCP input, cut into
        // the LF rows the index will have.
        TmsLcpSummary summary;
        double summary_s = 0;
        if (lcp_summary) {
            const auto ts = clock::now();
            try {
                const orbit::int_vector rows = TmsIndex::lf_row_lengths(heads, lens, opts);
                ulint n = 0;
                for (ulint l : lens) n += l;
                auto sums = minima_path.empty()
                                ? ms_build::row_summaries(lens, rows, ms_build::lcp_bin_source(lcp_bin_path, lens, true))
                                : ms_build::row_summaries(lens, rows, ms_build::minima_source(minima_path, heads.size(), n));
                summary.top = std::move(sums.first);
                summary.int_min = std::move(sums.second);
            } catch (const std::exception& e) {
                std::cerr << "Failed to build LCP summaries: " << e.what() << "\n";
                return 1;
            }
            opts.lcp_summary = &summary;
            summary_s = std::chrono::duration<double>(clock::now() - ts).count();
            mem_report::phase("tms_lcp_summary", heads.size());
        }
        // The run lengths and tops, packed.
        orbit::int_vector run_lens, run_tops;
        {
            ulint max_top = 0;
            for (ulint t : tops) max_top = std::max(max_top, t);
            run_tops = orbit::int_vector(tops.size(), orbit::bit_width(max_top));
            for (size_t i = 0; i < tops.size(); ++i) run_tops.set(i, tops[i]);
            std::vector<ulint>().swap(tops);
        }
        {
            ulint max_len = 0;
            for (ulint l : lens) max_len = std::max(max_len, l);
            run_lens = orbit::int_vector(lens.size(), orbit::bit_width(max_len));
            for (size_t i = 0; i < lens.size(); ++i) run_lens.set(i, lens[i]);
            std::vector<ulint>().swap(lens);
        }
        const double load_s = std::chrono::duration<double>(clock::now() - t0).count();
        mem_report::phase("tms_load", heads.size());
        auto t1 = clock::now();
        TmsIndex idx = phi_path.empty() ? TmsIndex(heads, run_lens, opts, with_phi ? &run_tops : nullptr)
                                        : TmsIndex(heads, run_lens, opts, std::move(samples));
        const double build_s = std::chrono::duration<double>(clock::now() - t1).count();
        const std::string idx_path = argv[npos - 1];
        if (seq_ids) {
            try {
                idx.set_seq_ids(std::move(*seq_ids));
            } catch (const std::exception& e) {
                std::cerr << seq_ids_path << " does not describe the indexed text: " << e.what() << "\n";
                return 1;
            }
        }
        idx.set_info(build_info());
        std::ofstream out(idx_path, std::ios::binary);
        size_t bytes = 0;
        try {
            if (out.good()) bytes = idx.serialize(out);
        } catch (const std::exception& e) {
            std::cerr << e.what() << "\n";
            out.setstate(std::ios::failbit);
        }
        if (!out.good()) {
            std::cerr << "Failed to write index: " << idx_path << "\n";
            return 1;
        }
        mem_report::phase("tms_write", heads.size());
        idx.describe(std::cout);
        std::cout << "Built tms index: " << idx_path << " (runs=" << heads.size() << ", bytes=" << bytes
                  << ", load_s=" << load_s << ", build_s=" << build_s;
        if (lcp_summary) std::cout << ", lcp_summary_s=" << summary_s << " (in load_s)";
        std::cout << ")\n";
        return 0;
    }

    if (cmd == "ftab-build") {
        if (argc < 2) { std::cerr << "ftab-build requires INDEX and FTAB\n"; return 1; }
        const std::string idx_path = argv[0], out_path = argv[1];
        unsigned k = 0;  // 0: --auto
        size_t threads = 0;
        std::string v;
        for (int i = 2; i < argc; ++i) {
            if (opt_value(argc, argv, i, "--kmer-len", v)) {
                k = static_cast<unsigned>(opt_uint("--kmer-len", v));
                if (k < 1 || k > ftab::max_k) { std::cerr << "--kmer-len must be between 1 and 31\n"; return 1; }
            }
            else if (strcmp(argv[i], "--auto") == 0) k = 0;
            else if (opt_value(argc, argv, i, "--threads", "-p", v)) threads = opt_uint("--threads", v);
            else { std::cerr << "Unknown ftab-build option: " << argv[i] << "\n"; return 1; }
        }
        threads = resolve_thread_count(threads);
        using clock = std::chrono::steady_clock;
        try {
            const orbit_idx::Header ih = orbit_idx::read_header_file(idx_path);
            const uint64_t index_bytes = std::filesystem::file_size(idx_path);
            const auto t0 = clock::now();
            std::optional<MSIndexSpillLCP<false>> ms_idx;
            std::optional<TmsIndex> tms_idx;
            if (ih.tool == orbit_idx::Tool::MS) ms_idx = ms_serialize::read_index(idx_path);
            else tms_idx = read_tms_index(idx_path);
            if (!ms_idx && !tms_idx) { std::cerr << "Failed to load index: " << idx_path << "\n"; return 1; }
            const double load_s = std::chrono::duration<double>(clock::now() - t0).count();
            const auto t1 = clock::now();
            const ftab::Shape base = ms_idx ? ftab::ms_shape(*ms_idx, threads) : ftab::tms_shape(*tms_idx, threads);
            auto build = [&](const ftab::Shape& s) {
                return ms_idx ? ftab::build_ms(*ms_idx, s, threads) : ftab::build_tms(*tms_idx, s, threads);
            };
            ftab::AutoChoice choice;
            const ftab::AutoRule rule;
            ftab::Built b;
            if (k == 0) b = ftab::build_auto(base, index_bytes, rule, build, choice);
            else {
                ftab::Shape s = base;
                s.k = k;
                b = build(s);
            }
            const double build_s = std::chrono::duration<double>(clock::now() - t1).count();
            ftab::FileHeader h = ftab::make_header(ih, b.tab);
            h.rule = k == 0 ? ftab::KRule::AUTO : ftab::KRule::GIVEN;
            h.present = b.present[b.tab.k()];
            h.built_at = static_cast<uint64_t>(std::time(nullptr));
            h.index_bytes = index_bytes;
            std::ostringstream info;
            info << build_info() << "commit=" << MS_COMMIT << "\nindex_path=" << idx_path << "\n";
            if (k == 0) {
                info << "k_rule=auto: largest k with table bytes <= " << rule.budget << " of index bytes (k=" << choice.k_fit
                     << "), then down while occupancy < " << rule.min_occupancy << "\noccupancy=";
                const unsigned lo = choice.k_fit > 4 ? choice.k_fit - 4 : 1;
                for (unsigned j = lo; j <= choice.k_fit; ++j)
                    info << (j > lo ? " " : "") << j << ":" << std::setprecision(4) << ftab::occupancy(choice.present, j);
                info << "\n";
            } else {
                info << "k_rule=given\n";
            }
            info << "threads=" << threads << "\nbuild_s=" << build_s << "\n";
            h.info = info.str();
            ftab::write_file(out_path, h, b.tab);
            std::cout << "Built ftab: " << out_path << " (k=" << b.tab.k() << ", entry_bits=" << b.tab.entry_bits()
                      << ", bytes=" << b.tab.bytes() << ", occupancy=" << h.occupancy() << ", load_s=" << load_s
                      << ", build_s=" << build_s << ")\n";
        } catch (const std::exception& e) {
            std::cerr << "ftab-build: " << e.what() << "\n";
            return 1;
        }
        return 0;
    }

    if (cmd == "ftab-verify") {
        if (argc < 2) { std::cerr << "ftab-verify requires INDEX and FTAB\n"; return 1; }
        bool crc = false;
        for (int i = 2; i < argc; ++i) {
            if (strcmp(argv[i], "--crc") == 0) crc = true;
            else { std::cerr << "Unknown ftab-verify option: " << argv[i] << "\n"; return 1; }
        }
        try {
            ftab::FileHeader h = ftab::read_header_file(argv[1]);
            const std::string bad = ftab::link_mismatch(h, orbit_idx::read_header_file(argv[0]));
            if (!bad.empty()) {
                std::cerr << argv[1] << " was not built from " << argv[0] << ": " << bad << "\n";
                return 1;
            }
            if (crc) ftab::read_file(argv[1], h, true);
        } catch (const std::exception& e) {
            std::cerr << e.what() << "\n";
            return 1;
        }
        std::cout << argv[1] << " matches " << argv[0] << (crc ? ", and its payload CRC-32s match" : "") << "\n";
        return 0;
    }

    if (cmd == "ftab-info") {
        if (argc < 1) { std::cerr << "ftab-info requires FTAB\n"; return 1; }
        try {
            ftab::describe(ftab::read_header_file(argv[0]), std::cout);
        } catch (const std::exception& e) {
            std::cerr << e.what() << "\n";
            return 1;
        }
        return 0;
    }

    if (cmd == "ftab-survey") {
        if (argc < 1) { std::cerr << "ftab-survey requires INDEX\n"; return 1; }
        unsigned kmin = 6, kmax = 15, verify = 0;
        size_t threads = 0;
        std::string v;
        for (int i = 1; i < argc; ++i) {
            if (opt_value(argc, argv, i, "--kmin", v)) kmin = static_cast<unsigned>(opt_uint("--kmin", v));
            else if (opt_value(argc, argv, i, "--kmax", v)) kmax = static_cast<unsigned>(opt_uint("--kmax", v));
            else if (opt_value(argc, argv, i, "--threads", "-p", v)) threads = opt_uint("--threads", v);
            else if (opt_value(argc, argv, i, "--verify", v)) verify = static_cast<unsigned>(opt_uint("--verify", v));
            else { std::cerr << "Unknown ftab-survey option: " << argv[i] << "\n"; return 1; }
        }
        if (kmin < 1 || kmax < kmin || kmax > 31) { std::cerr << "ftab-survey needs 1 <= kmin <= kmax <= 31\n"; return 1; }
        threads = resolve_thread_count(threads);
        using clock = std::chrono::steady_clock;
        const auto t0 = clock::now();
        auto opt = ms_serialize::read_index(argv[0]);
        if (!opt) { std::cerr << "Failed to load index: " << argv[0] << "\n"; return 1; }
        auto& idx = *opt;
        const double load_s = std::chrono::duration<double>(clock::now() - t0).count();
        const ulint file_bytes = static_cast<ulint>(std::filesystem::file_size(argv[0]));
        const ulint intervals = idx.move_runs();
        const ulint max_len = ftab_survey::max_interval_length(idx, threads);
        // An entry names an interval and an offset within it, plus one bit
        // marking a k-mer that does not occur.
        const unsigned int_bits = ftab_survey::bits_for(intervals - 1);
        const unsigned off_bits = ftab_survey::bits_for(max_len > 0 ? max_len - 1 : 0);
        const unsigned entry_bits = int_bits + off_bits + 1;
        std::cout << "index: " << argv[0] << " file_bytes=" << file_bytes << " intervals=" << intervals
                  << " n=" << idx.domain() << " max_interval_length=" << max_len << " interval_bits=" << int_bits
                  << " offset_bits=" << off_bits << " entry_bits=" << entry_bits
                  << " len_width=" << int(idx.get_widths()[0]) << " threads=" << threads << " load_s=" << load_s
                  << "\n";
        const auto res = ftab_survey::survey(idx, kmax, threads);
        std::cout << "survey: kmax=" << kmax << " steps=" << res.steps << " seconds=" << res.seconds << "\n";
        std::cout << "k\tpresent\tkmers\toccupancy\tbits_per_entry\ttable_bytes\tpct_of_index\n";
        for (unsigned k = kmin; k <= kmax; ++k) {
            const ulint kmers = ulint{1} << (2 * k);
            const ulint table_bytes = (kmers * entry_bits + 7) / 8;
            std::cout << k << '\t' << res.present[k] << '\t' << kmers << '\t' << std::setprecision(6)
                      << double(res.present[k]) / double(kmers) << '\t' << entry_bits << '\t' << table_bytes << '\t'
                      << 100.0 * double(table_bytes) / double(file_bytes) << "\n";
        }
        for (unsigned k = 1; k <= std::min(verify, kmax); ++k) {
            const ulint bf = ftab_survey::brute_force(idx, k);
            std::cout << "verify: k=" << k << " survey=" << res.present[k] << " brute_force=" << bf
                      << (bf == res.present[k] ? " ok" : " mismatch") << "\n";
            if (bf != res.present[k]) return 1;
        }
        return 0;
    }

    if (cmd == "index-plan") {
        if (argc < 1) { std::cerr << "index-plan requires INDEX\n"; return 1; }
        orbit_idx::Header h;
        try {
            h = orbit_idx::read_header_file(argv[0]);
        } catch (const std::exception& e) {
            std::cerr << e.what() << "\n";
            return 1;
        }
        bool query = false;
        for (int i = 1; i < argc; ++i) {
            const int q = tms_query_flag(argc, argv, i);
            if (q < 0) return 1;
            if (q == 0) { std::cerr << "Unknown index-plan option: " << argv[i] << "\n"; return 1; }
            query = true;
        }
        if (query && h.tool != orbit_idx::Tool::TMS) { std::cerr << "query options apply to tms indexes only\n"; return 1; }
        orbit_idx::describe_header(h, std::cout);
        auto mb = [](size_t bytes) {
            std::ostringstream os;
            os << std::fixed << std::setprecision(bytes < 10'000'000 ? 3 : 1) << bytes / 1e6 << " MB";
            return os.str();
        };
        if (h.tool != orbit_idx::Tool::TMS) {
            size_t total = 0;
            for (const auto& e : h.sections) total += e.bytes;
            std::cout << "load plan of this build: every ms query loads every section = " << mb(total) << "\n";
            return 0;
        }
        // What this build loads for each tms query, from the
        // header's sections: TmsIndex::sections_for, and the start tables
        // load() builds with phi_inv.
        const TmsParts held = TmsIndex::held_parts(h);
        auto serves = [&](const TmsParts& need) {
            return (!need.psi || held.psi) && (!need.phi || held.phi) && (!need.phi_inv || held.phi_inv);
        };
        auto describe_plan = [&](const TmsParts& need) {
            std::ostringstream os;
            const uint32_t mask = TmsIndex::sections_for(need);
            size_t total = 0;
            for (const auto& e : h.sections)
                if (mask & orbit_idx::section_bit(e.id)) {
                    os << " " << orbit_idx::section_name(e.id);
                    total += e.bytes;
                }
            if (need.phi && need.phi_inv) {
                for (auto id : {orbit_idx::Section::PHI_START, orbit_idx::Section::PHI_INV}) {
                    const orbit_idx::SectionEntry* e = h.find(id);
                    if (!e || !(mask & orbit_idx::section_bit(id))) continue;
                    os << " +" << orbit_idx::section_name(id) << "_start_table(shift="
                       << TmsIndex::start_table_shift(e->domain, e->rows) << ")";
                    total += TmsIndex::start_table_bytes(e->domain, e->rows);
                }
            }
            os << "  = " << mb(total);
            return os.str();
        };
        const char* modes[] = {"psi", "phi", "phisum", "phiskip", "dual", "oracle", "oracle-scan"};
        const TmsMode mode_of[] = {TmsMode::PSI,  TmsMode::PHI,    TmsMode::PHISUM,     TmsMode::PHISKIP,
                                   TmsMode::DUAL, TmsMode::ORACLE, TmsMode::ORACLE_SCAN};
        const char* reports[] = {"ms", "smem-one", "smem-all"};
        const TmsReport report_of[] = {TmsReport::MS, TmsReport::SMEM_ONE, TmsReport::SMEM_ALL};
        if (query) {
            const TmsParts need = tms_query_load();
            if (!serves(need)) { std::cout << "this query: refused, the index lacks a part it needs\n"; return 0; }
            std::cout << "this query loads (" << TmsIndex::load_name(need) << "):" << describe_plan(need) << "\n";
            return 0;
        }
        std::vector<std::string> names, plans;
        std::vector<std::string> rows;
        for (int m = 0; m < 7; ++m)
            for (int r = 0; r < 3; ++r)
                for (int pos = 0; pos < 2; ++pos)
                    for (int ps = 0; ps < 2; ++ps) {
                        g_tms_mode = mode_of[m];
                        g_tms_report = report_of[r];
                        g_tms_positions = pos;
                        g_tms_phi_starts = ps;
                        const TmsParts need = tms_query_load();
                        const std::string name = serves(need) ? TmsIndex::load_name(need) : "(refused: the index lacks a part)";
                        if (serves(need) && std::find(names.begin(), names.end(), name) == names.end()) {
                            names.push_back(name);
                            plans.push_back(describe_plan(need));
                        }
                        std::ostringstream row;
                        row << "  " << std::left << std::setw(12) << modes[m] << std::setw(9) << reports[r] << std::setw(12)
                            << (pos ? "--positions" : "") << std::setw(13) << (ps ? "--phi-starts" : "") << std::right
                            << name;
                        rows.push_back(row.str());
                    }
        std::cout << "load plans of this build (sections read, start tables built, and the RAM they take):\n";
        for (size_t i = 0; i < names.size(); ++i) std::cout << "  " << names[i] << ":" << plans[i] << "\n";
        std::cout << "queries on this index (--mode, --report, then --positions and --phi-starts):\n";
        for (const auto& r : rows) std::cout << r << "\n";
        return 0;
    }

    if (cmd == "index-verify") {
        if (argc < 1) { std::cerr << "index-verify requires INDEX\n"; return 1; }
        const std::string bad = orbit_idx::verify_file(argv[0]);
        if (!bad.empty()) { std::cerr << argv[0] << ": " << bad << "\n"; return 1; }
        std::cout << argv[0] << ": every section's CRC-32 matches\n";
        return 0;
    }

    if (cmd == "index-add-seq-ids") {
        if (argc < 3) { std::cerr << "index-add-seq-ids requires INDEX, SEQ_IDS and OUT\n"; return 1; }
        try {
            const orbit_idx::Header h = orbit_idx::read_header_file(argv[0]);
            if (h.tool != orbit_idx::Tool::TMS) throw std::runtime_error(std::string(argv[0]) + " is not a tms index");
            const SeqIdTable t = SeqIdTable::read(argv[1]);
            // The header gives the text length but not the string count,
            // which would take loading LF.
            const std::string bad = TmsIndex::seq_ids_mismatch(t, h.n);
            if (!bad.empty()) throw std::runtime_error(std::string(argv[1]) + " does not describe the indexed text: " + bad);
            std::string info = "seq_ids_from=" + std::string(argv[1]) + "\n" + build_info();
            info.replace(info.find("built_by="), 9, "seq_ids_added_by=");
            info.replace(info.find("built_at="), 9, "seq_ids_added_at=");
            // Earlier provenance of a table this replaces stays, above it.
            const auto t0 = std::chrono::steady_clock::now();
            orbit_idx::SectionSource src = TmsIndex::seq_ids_source(t);
            orbit_idx::add_section(argv[0], argv[2], src, info);
            std::cout << "Wrote " << argv[2] << ": " << argv[0] << " with " << t.names.size() << " sequences from "
                      << argv[1] << " (s=" << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count()
                      << ")\n";
        } catch (const std::exception& e) {
            std::cerr << "index-add-seq-ids: " << e.what() << "\n";
            return 1;
        }
        return 0;
    }
    if (cmd == "index-rewrap") {
        if (argc < 2) { std::cerr << "index-rewrap requires LEGACY_INDEX and INDEX\n"; return 1; }
        std::string info = "rewrapped_from=" + std::string(argv[0]) + "\n" + build_info();
        info.replace(info.find("built_by="), 9, "rewrapped_by=");
        info.replace(info.find("built_at="), 9, "rewrapped_at=");
        const auto t0 = std::chrono::steady_clock::now();
        const std::string err = orbit_idx::rewrap_legacy(argv[0], argv[1], info);
        if (!err.empty()) { std::cerr << "index-rewrap: " << err << "\n"; return 1; }
        std::cout << "Rewrapped " << argv[0] << " as " << argv[1] << " (s="
                  << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << ")\n";
        return 0;
    }

    if (cmd == "split-report") {
        if (argc < 1) { std::cerr << "split-report requires INDEX\n"; return 1; }
        size_t threads = 0;
        std::string v;
        for (int i = 1; i < argc; ++i) {
            if (opt_value(argc, argv, i, "--threads", "-p", v)) threads = opt_uint("--threads", v);
            else { std::cerr << "Unknown split-report option: " << argv[i] << "\n"; return 1; }
        }
        return split_report::run(argv[0], resolve_thread_count(threads), std::cout);
    }

    if (cmd == "tms-inspect") {
        if (argc < 1) { std::cerr << "tms-inspect requires INDEX\n"; return 1; }
        // --phi-lengths loads LF and phi's lengths-based form alone, which
        // takes far less memory than a whole phi or full index.
        const TmsParts lf_phi_len{false, true, false, false};
        const bool phi_len_only = argc > 1 && strcmp(argv[1], "--phi-lengths") == 0;
        auto opt = read_tms_index(argv[0], phi_len_only ? &lf_phi_len : nullptr);
        if (!opt) { std::cerr << "Failed to load index: " << argv[0] << "\n"; return 1; }
        opt->describe(std::cout);
        return 0;
    }

    if (cmd == "tms-text") {
        if (argc < 1) { std::cerr << "tms-text requires INDEX\n"; return 1; }
        // LF alone reads the text back.
        const TmsParts lf_only;
        auto opt = read_tms_index(argv[0], &lf_only);
        if (!opt) { std::cerr << "Failed to load index: " << argv[0] << "\n"; return 1; }
        auto show = [](uchar c) {
            return c == orbit::TERMINATOR ? '$' : c == orbit::SEPARATOR ? '%' : c == ms_alpha::TEXT_N ? 'N'
                                                                                     : static_cast<char>(c);
        };
        const ulint n = opt->domain();
        ulint m = 0;
        for (ulint i = 0; i < opt->move_runs(); ++i)
            if (opt->get_character(i) == orbit::TERMINATOR) m += opt->get_length(i);
        if (m <= 1) {
            // One text, so LF is one cycle.  Row 0 holds text position
            // n - 1, so the s-th row LF visits from it has BWT character
            // T[n - 2 - s].
            std::string t(n, '\0');
            auto pos = opt->first();
            for (ulint st = 0; st < n; ++st) {
                t[(2 * n - 2 - st) % n] = show(opt->get_character(pos.interval));
                pos = opt->LF_step(pos);
            }
            std::cout << t << "\n";
            return 0;
        }
        // A string collection, with one terminator per string, whose LF may
        // have several cycles.  The k-th of the first m rows is the suffix
        // that is only some string's terminator, so LF from it reads that
        // string backward, up to the row whose BWT character is a terminator.
        std::string t, line;
        auto row = opt->first();
        for (ulint k = 0; k < m; ++k) {
            if (k > 0) {
                if (row.offset + 1 < opt->get_length(row.interval)) ++row.offset;
                else row = opt->down(row);  // the next interval's first row
            }
            t.clear();
            for (auto pos = row;;) {
                const uchar c = opt->get_character(pos.interval);
                if (c == orbit::TERMINATOR) break;
                t.push_back(show(c));
                pos = opt->LF_step(pos);
            }
            line.assign(t.rbegin(), t.rend());
            std::cout << line << "$\n";
        }
        return 0;
    }

    if (cmd == "fmd-test") {
        std::string data_dir = (argc > 0 && argv[0][0] != '-') ? argv[0] : "./data";
        return fmd_test::run_all_tests(data_dir) ? 0 : 1;
    }

    if (cmd == "tms-test") {
        std::string data_dir = (argc > 0 && argv[0][0] != '-') ? argv[0] : "./data";
        return tms_test::run_all_tests(data_dir) ? 0 : 1;
    }

    if (cmd == "inspect") {
        if (argc < 1) {
            std::cerr << "inspect requires INDEX\n";
            return 1;
        }
        auto opt = ms_serialize::read_index(argv[0]);
        if (!opt) {
            std::cerr << "Failed to load index: " << argv[0] << "\n";
            return 1;
        }
        ms_inspect::run_inspect(*opt);
        for (int i = 1; i < argc; ++i) {
            if (std::string path; opt_value(argc, argv, i, "--spillover-tsv", path)) {
                if (ms_inspect::run_spillover_tsv(*opt, path))
                    std::cout << "Spillover TSV written to " << path << "\n";
            }
        }
        return 0;
    }

    std::cerr << "Unknown command: " << cmd << " (see ms --help)\n";
    return 1;
}

/** Main entry point. */
int main(int argc, char** argv) {
    try {
        return run_command(argc, argv);
    } catch (const std::invalid_argument& e) {
        // A malformed option value, from opt_uint or a std::sto* call.
        std::cerr << e.what() << "\n";
        return 1;
    } catch (const std::out_of_range& e) {
        std::cerr << "an option value is out of range: " << e.what() << "\n";
        return 1;
    }
}
