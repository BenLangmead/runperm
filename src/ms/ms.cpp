/**
 * Command line interface that displays usage and dispatches to the various
 * subcommands.  TODO: slicker command-line parsing.
 *
 * Author: Ben Langmead (ben.langmead@gmail.com)
 * Date: Feb 17, 2026
 */

#include "ms_rlbwt.hpp"
#include "ms_build.hpp"
#include "tsv.hpp"
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

/** Print usage message. */
static void usage(const char* prog) {
    std::cerr << "Usage: " << prog << " <cmd> [args...]\n"
              << "\n"
              << "Main commands:\n"
              << "  build      TSV_PATH INDEX_PATH [--percentile K] [--split-threshold N]\n"
              << "            [--coalesce-spillover] [--spill-align N] [--spill-split-bits X]\n"
              << "            [--minima-only] [--lf-split B]\n"
              << "              Build index from TSV.  LF is always split: the rows (runs, or\n"
              << "              their pieces after --split-threshold) are cut into LF intervals\n"
              << "              with Orbit's default length capping and balancing, each with\n"
              << "              its own LCP data.  --lf-split B sets the balancing factor\n"
              << "              (default: Orbit's default balancing); it must be positive.\n"
              << "              Results do not depend on B.\n"
              << "  build-rlbwt HEADS LENS INDEX_PATH (--minima FILE | --lcp-bin FILE)\n"
              << "            [--minima-only] [--percentile K] [--coalesce-spillover]\n"
              << "            [--spill-align N] [--spill-split-bits X] [--lf-split B]\n"
              << "              Build index from an RLBWT (HEADS: one byte per run; LENS:\n"
              << "              fixed-width little-endian lengths) and either a TeraLCP\n"
              << "              -ominima file, which keeps only the LCP minima that queries\n"
              << "              need, or --lcp-bin, one 64-bit LCP per row, optionally after a\n"
              << "              64-bit row count.  --minima-only applies to --lcp-bin as to\n"
              << "              build.  --lf-split as for build.\n"
              << "  ms         INDEX_PATH PATTERN\n"
              << "              Compute matching statistics for PATTERN using INDEX_PATH.\n"
              << "  batch      INDEX_PATH READS [-o OUT] [--no-output] [--interleave K]\n"
              << "            [--threads T] [--block-reads B] [--block-bytes N]\n"
              << "            [--orient fwd|rc|qual] [--ftab-file FILE | --ftab | --ftab-k K]\n"
              << "            [--mode standard|oracle|oracle-scan] [--oracle-write FILE | --oracle FILE]\n"
              << "              Load the index once and compute matching statistics for every\n"
              << "              read in READS (FASTA, FASTQ, or one sequence per line, in\n"
              << "              which blank lines are skipped; - for stdin; LF or CRLF line\n"
              << "              endings).  Reads are processed in blocks, so input size is not\n"
              << "              limited by memory.  Writes one line per read: name, tab, space-separated\n"
              << "              values.  --no-output skips writing (for timing).  --interleave K\n"
              << "              keeps K reads in flight, prefetching each one's next row\n"
              << "              (default 32); K = 0 queries one read at a time without\n"
              << "              prefetching (ms_query for batch, the batched engine with one\n"
              << "              read for tms-batch).  --threads T (or -p T) queries on T\n"
              << "              threads, each with its own K reads in flight (default 1; 0 uses\n"
              << "              every CPU the process may run on).  Each thread takes a block\n"
              << "              of about B reads at a time (--block-reads, default\n"
              << "              max(512, 16K)), sized in bytes from the reads seen so far,\n"
              << "              or of N bytes (--block-bytes).  Output named *.gz is written\n"
              << "              gzip-compressed, by all threads (--gzip-level, default 1).\n"
              << "              Results and their order do not depend on K, T or B.  On a\n"
              << "              machine with several memory nodes (sockets), the index's\n"
              << "              memory is interleaved over the nodes (Linux), unless numactl\n"
              << "              sets a memory policy, as when running one process per node on\n"
              << "              part of the reads (numactl -N i -m i).  A summary goes to stderr:\n"
              << "              query_ns_per_base is the threads' summed query time per base,\n"
              << "              and wall_ns_per_base the elapsed time per base.\n"
              << "              --orient sets which strand of each read is matched: fwd (the\n"
              << "              default) matches the read as given; rc matches its reverse\n"
              << "              complement, starting from the read's 5' end, and converts the\n"
              << "              result back to the read's matching statistics; qual does so\n"
              << "              only for reads whose last min(20, len) bases have a higher\n"
              << "              mean base quality than their first, and needs FASTQ input.\n"
              << "              rc and qual give the same output as fwd only for an index that\n"
              << "              holds both strands of every sequence (as an FMD or a ropebwt3\n"
              << "              index does); this cannot be detected from the index, so the\n"
              << "              user must know it.  They refuse tms-batch's --positions and\n"
              << "              SMEM reports.  The summary's flipped field counts the reads\n"
              << "              matched on the other strand, and orient_s is the threads'\n"
              << "              summed time flipping reads and converting their results.\n"
              << "              --ftab-file FILE loads a table (from ftab-build) of the query\n"
              << "              state after each K-mer over ACGT; a read whose last K bases\n"
              << "              occur in the text starts K steps in.  Results are unchanged.\n"
              << "              The table must have been built from this index: its link to\n"
              << "              the index is checked against the index's header, and a table\n"
              << "              of another index is refused.  --ftab instead builds the table\n"
              << "              after loading the index, on every CPU (K = 14, or --ftab-k K,\n"
              << "              which implies --ftab).  tms-batch takes a table in every mode,\n"
              << "              with or without positions; a table built from a psi layout\n"
              << "              serves --mode psi without positions or SMEM reports only.\n"
              << "              --oracle-write FILE (in the standard mode) records every\n"
              << "              reposition's outcome in the format tms-batch's describes, and\n"
              << "              --mode oracle or oracle-scan with --oracle FILE replays them.\n"
              << "              oracle moves to the recorded run with the recorded length,\n"
              << "              without walking for the run or reading LCP bounds; oracle-scan\n"
              << "              walks only the side the recorded run is on, reading each run's\n"
              << "              character but no LCP bounds, and takes the recorded length.\n"
              << "              Their output is the standard mode's.  Neither records nor\n"
              << "              replays with an ftab, which skips some repositions.\n"
              << "  tms-build  HEADS LENS INDEX_PATH [--layout full|psi|phi] [--phi FILE]\n"
              << "            [--minima FILE] [--lcp-bin FILE] [--lf-split B] [--phi-split B]\n"
              << "            [--no-phi-inv] [--seq-ids FILE]\n"
              << "              Build a tms index from an RLBWT: LF and psi, which need no LCP\n"
              << "              input, and with LCP input also phi and phi_inv (phi_inv unless\n"
              << "              --no-phi-inv; it is needed only for tms-batch --report smem-all).\n"
              << "              --phi takes a TeraLCP -ophi file, the suffix array and LCP samples\n"
              << "              of every run, and builds phi from them in O(r) space; the RLBWT\n"
              << "              may then have one terminator per string, as an FMD's runs do.\n"
              << "              --minima (a TeraLCP -ominima file, of which only each run's top\n"
              << "              LCP is used) or --lcp-bin (one 64-bit LCP per row) instead walk LF\n"
              << "              over all n rows to find the suffix array samples, which needs a\n"
              << "              BWT with one terminator.  --layout keeps only what\n"
              << "              some queries read: full (the default) keeps everything; psi keeps\n"
              << "              LF and psi, takes no LCP input, and serves only --mode psi without\n"
              << "              --positions or an SMEM report; phi keeps LF, phi and phi_inv,\n"
              << "              needs LCP input, and serves every query but --mode psi, phiskip\n"
              << "              and dual.  In the psi and full layouts, LF and psi share one\n"
              << "              table of rows, the union of LF's split intervals and their\n"
              << "              images.  --lf-split and --phi-split set Orbit's balancing factor\n"
              << "              for LF and phi, with its default length capping (default: Orbit's\n"
              << "              default balancing).  LF is always split, so --lf-split must be\n"
              << "              positive; --phi-split 0 leaves phi unsplit.  --seq-ids FILE\n"
              << "              stores ropebwt3's sequence-ID table (the .len.gz of the FMD the\n"
              << "              RLBWT came from) in the index, for tms-batch --report smem-all;\n"
              << "              its contig k must be the text's strings 2k and 2k + 1, and both\n"
              << "              its text length and its string count must be the index's.\n"
              << "  tms-build-tsv TSV_PATH INDEX_PATH [--phi] [--layout ...] [split options]\n"
              << "              Same, taking the runs and their top LCPs from a TSV.\n"
              << "  tms-batch  INDEX_PATH READS [-o OUT] [--no-output] [--interleave K]\n"
              << "            [--threads T] [--block-reads B] [--block-bytes N]\n"
              << "            [--mode psi|phi|phiskip|dual|oracle|oracle-scan] [--positions]\n"
              << "            [--report ms|smem-one|smem-all] [--min-smem-len T]\n"
              << "            [--max-smem-positions N] [--phi-starts] [--orient fwd|rc|qual]\n"
              << "            [--ftab-file FILE | --ftab | --ftab-k K]\n"
              << "            [--oracle-write FILE | --oracle FILE]\n"
              << "              As batch, with a tms index.  --mode sets how repositions\n"
              << "              compute LCEs (default psi; psi, phiskip and dual need psi, and\n"
              << "              the others need phi).  Only the index's structures that the\n"
              << "              query reads are loaded, and a query the index's layout cannot\n"
              << "              serve is refused.  Each read gets a line with its name and its\n"
              << "              matching statistics, space-separated.  --positions adds a field\n"
              << "              with an occurrence position for each value (-1 where it is 0);\n"
              << "              it needs phi.  --report ms (the default) adds nothing more;\n"
              << "              smem-one adds a field with the read's SMEMs (super-maximal exact\n"
              << "              matches), space-separated, each as i:L:p, with start i in the\n"
              << "              read, length L and one text position p where it occurs, and\n"
              << "              needs phi.  smem-all writes SMEMs as ropebwt3 mem -p N does, in\n"
              << "              place of the read's line, and needs phi and phi_inv: one line\n"
              << "              per SMEM, tab-separated, with the read's name, start, end\n"
              << "              (exclusive), the count c of its occurrences and, if N > 0,\n"
              << "              min(N, c) and that many positions, the first in BWT row order.\n"
              << "              If the index holds a sequence-ID table (tms-build --seq-ids,\n"
              << "              index-add-seq-ids), each position is contig:strand:offset, where\n"
              << "              a - offset is that of the reverse complement on the forward\n"
              << "              strand; otherwise each is a text position.  With --min-smem-len\n"
              << "              0, each read character that occurs nowhere in the text, such as\n"
              << "              an N, also gets a line with its start, end and a count of 0, as\n"
              << "              ropebwt3 gives with -l 1; reads with neither add no line.\n"
              << "              smem-all cannot be used with --positions.  --min-smem-len T\n"
              << "              keeps only SMEMs with L > T (default 0, all of them) in both\n"
              << "              SMEM reports, so ropebwt3's -l X is --min-smem-len X-1.\n"
              << "              --max-smem-positions N sets N\n"
              << "              (default: list all).  With smem-all, K also sets how many of its\n"
              << "              phi and phi_inv walks list positions at once (K = 0: one SMEM at\n"
              << "              a time).  Matching statistics without positions or SMEMs load\n"
              << "              phi's lengths-based form, whose rows are narrower; --phi-starts\n"
              << "              loads its starts-based form instead, as positions and SMEM\n"
              << "              reports do.  Results do not depend on the mode, K or\n"
              << "              --phi-starts.\n"
              << "              --oracle-write FILE (in a mode other than oracle and\n"
              << "              oracle-scan) records every\n"
              << "              reposition's outcome, and --mode oracle --oracle FILE replays\n"
              << "              them: on a mismatch it moves to the recorded row and length\n"
              << "              without scanning for the candidates or computing an LCE, so it\n"
              << "              reads only LF (and phi with --positions or an SMEM report).  It\n"
              << "              is a lower bound on the time of the other modes, and its output\n"
              << "              is theirs.  --mode oracle-scan replays the same file with less\n"
              << "              given: the recorded interval says only which side (pred above\n"
              << "              or succ below) to scan, and that side is scanned for the\n"
              << "              candidate as the other modes do; the recorded length still\n"
              << "              stands in for the LCE.  Replay needs the same reads, in the\n"
              << "              same order, and an index with the same LF intervals (the same\n"
              << "              build and layout);\n"
              << "              it does not depend on K, T or B.  FILE holds u64s, little-endian:\n"
              << "              the magic \"TMSORC01\" (8 bytes), the index's LF interval count,\n"
              << "              the read count R, the reposition count D, R per-read reposition\n"
              << "              counts in input order, then D (interval, length) pairs: the LF\n"
              << "              interval moved to (its tail row if above the current interval,\n"
              << "              else its head row) and the match length before the LF step.\n"
              << "  tms-text   INDEX_PATH\n"
              << "              Print the indexed text, read back with LF.  For a string\n"
              << "              collection (one terminator per string, as in an FMD), print\n"
              << "              each string and its $, one per line, in the BWT order of\n"
              << "              their terminators.  An FMD's N prints as N.\n"
              << "  ftab-build INDEX_PATH FTAB_PATH [--k K | --auto] [--threads T]\n"
              << "              Build the ftab of an index (see batch --ftab-file) and write it\n"
              << "              to its own file.  --k K sets K; --auto (the default) takes the\n"
              << "              largest K whose table fits in 5% of the index file's bytes, then\n"
              << "              steps K down by one while fewer than 80% of its K-mers occur.\n"
              << "              --threads T builds on T threads (default 0, every CPU).  A tms\n"
              << "              index with phi gets a table with a toehold and window\n"
              << "              positions, which serves every mode, positions and SMEM reports.\n"
              << "  ftab-verify INDEX_PATH FTAB_PATH [--crc]\n"
              << "              Check, from the two headers alone, that the ftab was built from\n"
              << "              this index; on a mismatch, name the field that differs and exit\n"
              << "              nonzero.  --crc also checks the ftab's payload CRC-32s.\n"
              << "  ftab-info  FTAB_PATH\n"
              << "              Print an ftab's header: k, fields, occupancy, link, provenance.\n"
              << "  ftab-survey INDEX_PATH [--kmin K] [--kmax K] [--threads T] [--verify K]\n"
              << "              Development: count the k-mers over ACGT that occur in the\n"
              << "              text, for each k from --kmin (default 6) to --kmax (default\n"
              << "              15), and the size of a dense table holding each k-mer's\n"
              << "              row after its cold-start walk.  --verify K checks the counts\n"
              << "              for k up to K against ms on every k-mer.\n"
              << "  split-report INDEX_PATH [--threads T]\n"
              << "              Development: for each move structure of an ms or tms index,\n"
              << "              print the split parameters stored with it, its\n"
              << "              intervals, longest interval and the length cap Orbit's length\n"
              << "              capping implies, and its balance: how many input interval\n"
              << "              starts lie inside each output interval, against the bound\n"
              << "              2b - 1 that balancing factor b guarantees.  tms parts are\n"
              << "              loaded a group at a time.\n"
              << "  index-plan INDEX_PATH [--mode M] [--positions] [--report R] [--phi-starts]\n"
              << "              Print an index's header, read from the file's first 4 KiB\n"
              << "              alone: its sections with their columns and widths.  Then print\n"
              << "              what this build loads for each query, and the RAM it takes: for\n"
              << "              a tms index, for each tms-batch command line, or with query\n"
              << "              options for that one alone.\n"
              << "  index-verify INDEX_PATH\n"
              << "              Check every section's CRC-32 and that the file is whole.\n"
              << "  index-rewrap LEGACY_PATH INDEX_PATH\n"
              << "              Rewrite an index of the earlier formats (MSIX v4 from build or\n"
              << "              build-rlbwt, TMSX v8 from tms-build) in the current format,\n"
              << "              copying its structures' bytes without loading them.\n"
              << "  index-add-seq-ids INDEX_PATH SEQ_IDS OUT_PATH\n"
              << "              Write the tms index at INDEX_PATH to OUT_PATH with the\n"
              << "              sequence-ID table SEQ_IDS (see tms-build --seq-ids) in it, in\n"
              << "              place of any it holds, copying its structures' bytes without\n"
              << "              loading them.  The table's text length must be the index's.\n"
              << "  tms-inspect INDEX_PATH [--phi-lengths]\n"
              << "              Print the structures' sizes and column widths.  --phi-lengths\n"
              << "              loads only LF and phi's lengths-based form.\n"
              << "  inspect    INDEX_PATH [--spillover-tsv FILE]\n"
              << "              Inspect index and optionally output spillover TSV.\n"
              << "  lcp-list   TSV_PATH\n"
              << "              Print LCP columns per run from TSV.\n"
              << "\n"
              << "Diagnostics and debugging:\n"
              << "  test           [DATA_DIR]\n"
              << "                  Run built-in tests (default DATA_DIR is ./data).\n"
              << "  fmd-test       [DATA_DIR]\n"
              << "                  Run only the string-collection (FMD) and large-n tests.\n"
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
 * The split parameters for --lf-split or --phi-split B: Orbit's default
 * length capping and balancing factor B, or none for B = 0.
 */
static orbit::split_params split_arg(const char* v) {
    const ulint b = std::stoull(v);
    return b == 0 ? orbit::NO_SPLITTING : orbit::split_params(orbit::DEFAULT_LENGTH_CAPPING, b);
}

/**
 * Parse --lf-split's value into sp.  Returns false, after saying why, for 0:
 * LF is always split.
 */
static bool lf_split_arg(const char* v, orbit::split_params& sp) {
    sp = split_arg(v);
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

// The ftab built at index load with --ftab, or null.
static const ftab::Ftab* g_ftab = nullptr;
// With --oracle-write, each query records its reposition outcomes; in the
// oracle modes, they are replayed from g_oracle, read from --oracle.
static bool g_record = false;
static RepositionOracle g_oracle;
// batch's --mode.
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
// tms-batch settings, from its command line.  They are read-only once the
// queries start.
static TmsMode g_tms_mode = TmsMode::PSI;
static bool g_tms_positions = false;
static bool g_tms_phi_starts = false;
static TmsReport g_tms_report = TmsReport::MS;
static ulint g_tms_min_smem_len = 0;
static ulint g_tms_max_smem_positions = TMS_ALL_POSITIONS;

// With an SMEM report, the SMEMs are found as part of the query.  smem_k
// walks are in flight when smem-all lists positions (0 lists them one SMEM
// at a time).  first_read is the index of p[0] among all the reads, which
// says where its recorded repositions are in the oracle mode.
static void query_many(TmsIndex& idx, const std::vector<std::string>& p, size_t k, QueryResults& r, size_t smem_k,
                       size_t first_read) {
    const bool smems = g_tms_report != TmsReport::MS;
    RepositionOracleView view{};
    const bool replay = g_tms_mode == TmsMode::ORACLE || g_tms_mode == TmsMode::ORACLE_SCAN;
    if (replay) view = RepositionOracleView::of(g_oracle, first_read, p.size());
    tms_query_batch(idx, p, k, r.ms, g_tms_mode, g_tms_positions || smems ? &r.pos : nullptr, true, g_ftab,
                    replay ? &view : nullptr, g_record ? &r.decisions : nullptr);
    if (!smems) return;
    tms_report_smems_batch(idx, r.ms, r.pos, g_tms_min_smem_len, g_tms_report, smem_k, r.hits,
                           g_tms_max_smem_positions);
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

/**
 * If argv[i] is a tms query flag that chooses what a query loads (--mode,
 * --positions, --report or --phi-starts), apply it to the tms-batch settings,
 * advancing i past its value, and return 1; return -1, after printing why, if
 * its value is unknown, and 0 if argv[i] is no such flag.
 */
static int tms_query_flag(int argc, char** argv, int& i) {
    if (strcmp(argv[i], "--positions") == 0) g_tms_positions = true;
    else if (strcmp(argv[i], "--phi-starts") == 0) g_tms_phi_starts = true;
    else if (strcmp(argv[i], "--mode") == 0 && i + 1 < argc) {
        const std::string m = argv[++i];
        if (m == "psi") g_tms_mode = TmsMode::PSI;
        else if (m == "phi") g_tms_mode = TmsMode::PHI;
        else if (m == "phiskip") g_tms_mode = TmsMode::PHISKIP;
        else if (m == "dual") g_tms_mode = TmsMode::DUAL;
        else if (m == "oracle") g_tms_mode = TmsMode::ORACLE;
        else if (m == "oracle-scan") g_tms_mode = TmsMode::ORACLE_SCAN;
        else { std::cerr << "Unknown mode: " << m << "\n"; return -1; }
    }
    else if (strcmp(argv[i], "--report") == 0 && i + 1 < argc) {
        const std::string r = argv[++i];
        if (r == "ms") g_tms_report = TmsReport::MS;
        else if (r == "smem-one") g_tms_report = TmsReport::SMEM_ONE;
        else if (r == "smem-all") g_tms_report = TmsReport::SMEM_ALL;
        else { std::cerr << "Unknown report: " << r << "\n"; return -1; }
    }
    else return 0;
    return 1;
}

/** The parts a tms-batch query with the current settings loads. */
static TmsParts tms_batch_parts() {
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

/**
 * ms batch and tms-batch: load the index once, then compute matching
 * statistics for every read, on one or more threads.  Each thread queries a
 * block of reads at a time with its own K reads in flight, and output is
 * written in input order, so it does not depend on the thread count.
 */
template <typename Index>
static int run_batch(int argc, char** argv, std::optional<Index> (*read)(const std::string&)) {
    if (argc < 2) {
        std::cerr << "batch requires INDEX_PATH and READS\n";
        return 1;
    }
    std::string idx_path = argv[0], reads_path = argv[1], out_path;
    bool write_output = true;
    size_t interleave = 32, threads = 1, block_reads = 0, block_bytes = 0;
    int gzip_level = -2;  // -2: from the output name; -1: none
    ReadOrient orient = ReadOrient::FWD;
    std::string oracle_path, oracle_write_path;
    int q = 0;
    unsigned ftab_k = 0;  // 0: no ftab built at load
    std::string ftab_path;  // a stored ftab to load
    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) out_path = argv[++i];
        else if (strcmp(argv[i], "--ftab") == 0) ftab_k = ftab_k ? ftab_k : 14;
        else if (strcmp(argv[i], "--ftab-k") == 0 && i + 1 < argc) ftab_k = static_cast<unsigned>(std::stoul(argv[++i]));
        else if (strcmp(argv[i], "--ftab-file") == 0 && i + 1 < argc) ftab_path = argv[++i];
        else if (strcmp(argv[i], "--no-output") == 0) write_output = false;
        else if (strcmp(argv[i], "--interleave") == 0 && i + 1 < argc)
            interleave = static_cast<size_t>(std::stoull(argv[++i]));
        else if ((strcmp(argv[i], "--threads") == 0 || strcmp(argv[i], "-p") == 0) && i + 1 < argc)
            threads = static_cast<size_t>(std::stoull(argv[++i]));
        else if (strcmp(argv[i], "--block-reads") == 0 && i + 1 < argc)
            block_reads = static_cast<size_t>(std::stoull(argv[++i]));
        else if (strcmp(argv[i], "--gzip-level") == 0 && i + 1 < argc)
            gzip_level = std::stoi(argv[++i]);
        else if (strcmp(argv[i], "--block-bytes") == 0 && i + 1 < argc)
            block_bytes = static_cast<size_t>(std::stoull(argv[++i]));
        else if (strcmp(argv[i], "--orient") == 0 && i + 1 < argc) {
            const std::string o = argv[++i];
            if (o == "fwd") orient = ReadOrient::FWD;
            else if (o == "rc") orient = ReadOrient::RC;
            else if (o == "qual") orient = ReadOrient::QUAL;
            else { std::cerr << "Unknown orientation: " << o << " (expected fwd, rc or qual)\n"; return 1; }
        }
        else if (std::is_same_v<Index, TmsIndex> && (q = tms_query_flag(argc, argv, i)) != 0) {
            if (q < 0) return 1;
        }
        else if (std::is_same_v<Index, TmsIndex> && strcmp(argv[i], "--min-smem-len") == 0 && i + 1 < argc)
            g_tms_min_smem_len = static_cast<ulint>(std::stoull(argv[++i]));
        else if (std::is_same_v<Index, TmsIndex> && strcmp(argv[i], "--max-smem-positions") == 0 && i + 1 < argc)
            g_tms_max_smem_positions = static_cast<ulint>(std::stoull(argv[++i]));
        else if (!std::is_same_v<Index, TmsIndex> && strcmp(argv[i], "--mode") == 0 && i + 1 < argc) {
            const std::string m = argv[++i];
            if (m == "standard") g_ms_mode = MsMode::STANDARD;
            else if (m == "oracle") g_ms_mode = MsMode::ORACLE;
            else if (m == "oracle-scan") g_ms_mode = MsMode::ORACLE_SCAN;
            else { std::cerr << "Unknown mode: " << m << "\n"; return 1; }
        }
        else if (strcmp(argv[i], "--oracle") == 0 && i + 1 < argc)
            oracle_path = argv[++i];
        else if (strcmp(argv[i], "--oracle-write") == 0 && i + 1 < argc)
            oracle_write_path = argv[++i];
        else { std::cerr << "Unknown batch option: " << argv[i] << "\n"; return 1; }
    }
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
    if (ftab_k > 31) { std::cerr << "--ftab-k must be between 1 and 31\n"; return 1; }
    if (ftab_k > 0 && !ftab_path.empty()) { std::cerr << "--ftab-file cannot be used with --ftab or --ftab-k\n"; return 1; }
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
    // The ftab is loaded, or built on every CPU, after the index, and is
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
                const size_t build_threads = resolve_thread_count(0);
                if constexpr (std::is_same_v<Index, TmsIndex>) tab = ftab::build_tms(idx, ftab_k, build_threads).tab;
                else tab = ftab::build_ms(idx, ftab_k, build_threads).tab;
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
        std::cerr << "batch: " << e.what() << "\n";
        return 1;
    }
    out->flush();
    const double total_s = std::chrono::duration<double>(clock::now() - t_all).count();
    parse_s += numbering_parse_s;
    if (oracle_mode && reads_numbered != g_oracle.patterns()) {
        std::cerr << "batch: the oracle holds " << g_oracle.patterns() << " reads but the input has "
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
            std::cerr << "batch: " << e.what() << "\n";
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
    perf.report(std::cerr, n_bases);
    // query_s sums the threads' time in queries, parse_s their time
    // parsing reads, and format_s their time formatting (and compressing)
    // output.  read_s is the time the input lock was held for reading; wall_ns_per_base is
    // the elapsed time of the whole read, query and write phase per base.
    std::cerr << "batch: reads=" << n_reads << " bases=" << n_bases << " threads=" << threads
              << " index_load_s=" << load_s << " ftab_k=" << ftab_k << " ftab_s=" << ftab_s << " query_s=" << query_s << " read_s=" << read_s
              << " parse_s=" << parse_s << " format_s=" << format_s
              << " orient=" << (orient == ReadOrient::FWD ? "fwd" : orient == ReadOrient::RC ? "rc" : "qual")
              << " flipped=" << n_flipped << " orient_s=" << orient_s
              << " total_s=" << total_s
              << " query_ns_per_base=" << (n_bases ? query_s * 1e9 / n_bases : 0.0)
              << " wall_ns_per_base=" << (n_bases ? total_s * 1e9 / n_bases : 0.0) << "\n";
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
        orbit::split_params lf_split;
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
            } else if (strcmp(argv[i], "--lf-split") == 0 && i + 1 < argc) {
                if (!lf_split_arg(argv[++i], lf_split)) return 1;
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
        const size_t runs = bwt_heads.size();
        apply_lcp_splitting(bwt_heads, bwt_run_lengths, lcps_per_run, split_threshold, minima_only);
        ms_build::Options bo;
        bo.percentile_k = percentile_k;
        bo.coalesce = coalesce;
        bo.spill_align = spill_align;
        bo.spill_split_bits = spill_split_bits;
        bo.minima_only = minima_only;
        bo.lf_split = lf_split;
        ms_build::Stats st;
        auto idx = ms_build::build(std::move(bwt_heads), std::move(bwt_run_lengths),
                                   ms_build::lcp_vectors_source(lcps_per_run, minima_only), bo, &st);
        idx.set_info(build_info());
        if (!ms_serialize::write_index(idx_path, idx)) {
            std::cerr << "Failed to write index: " << idx_path << "\n";
            return 1;
        }
        std::cout << "Built index: " << idx_path << " (runs=" << runs << ", rows=" << idx.move_runs()
                  << ", skinny=" << st.skinny << ", jumbo=" << st.jumbo << ")\n";
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
        if (argc < 3) {
            std::cerr << "build-rlbwt requires HEADS, LENS and INDEX_PATH\n";
            return 1;
        }
        const std::string heads_path = argv[0], lens_path = argv[1], idx_path = argv[2];
        std::string minima_path, lcp_bin_path;
        double percentile_k = 0.98;
        bool coalesce = false, minima_only = false;
        ulint spill_align = 0;
        uchar spill_split_bits = 0;
        orbit::split_params lf_split;
        for (int i = 3; i < argc; ++i) {
            if (strcmp(argv[i], "--minima") == 0 && i + 1 < argc) minima_path = argv[++i];
            else if (strcmp(argv[i], "--lcp-bin") == 0 && i + 1 < argc) lcp_bin_path = argv[++i];
            else if (strcmp(argv[i], "--minima-only") == 0) minima_only = true;
            else if (strcmp(argv[i], "--percentile") == 0 && i + 1 < argc) percentile_k = std::stod(argv[++i]);
            else if (strcmp(argv[i], "--coalesce-spillover") == 0) coalesce = true;
            else if (strcmp(argv[i], "--spill-align") == 0 && i + 1 < argc) spill_align = std::stoull(argv[++i]);
            else if (strcmp(argv[i], "--spill-split-bits") == 0 && i + 1 < argc)
                spill_split_bits = static_cast<uchar>(std::stoul(argv[++i]));
            else if (strcmp(argv[i], "--lf-split") == 0 && i + 1 < argc) {
                if (!lf_split_arg(argv[++i], lf_split)) return 1;
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

    if (cmd == "batch") {
        return run_batch<MSIndexSpillLCP<false>>(argc, argv, ms_serialize::read_index);
    }

    if (cmd == "tms-batch") {
        // Load only the parts the query reads; the load refuses an index
        // that lacks one.
        return run_batch<TmsIndex>(argc, argv, [](const std::string& path) {
            const TmsParts need = tms_batch_parts();
            return read_tms_index(path, &need);
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
        opts.report = true;
        bool with_phi = false;
        std::string minima_path, lcp_bin_path, phi_path, seq_ids_path;
        for (int i = npos; i < argc; ++i) {
            if (strcmp(argv[i], "--lf-split") == 0 && i + 1 < argc) {
                if (!lf_split_arg(argv[++i], opts.lf_split)) return 1;
            }
            else if (strcmp(argv[i], "--phi-split") == 0 && i + 1 < argc) opts.phi_split = split_arg(argv[++i]);
            else if (strcmp(argv[i], "--no-phi-inv") == 0) opts.phi_inv = false;
            else if (strcmp(argv[i], "--layout") == 0 && i + 1 < argc) {
                const std::string l = argv[++i];
                if (l == "full") opts.layout = TmsLayout::FULL;
                else if (l == "psi") opts.layout = TmsLayout::PSI;
                else if (l == "phi") opts.layout = TmsLayout::PHI;
                else { std::cerr << "Unknown layout: " << l << "\n"; return 1; }
            }
            else if (from_tsv && strcmp(argv[i], "--phi") == 0) with_phi = true;
            else if (!from_tsv && strcmp(argv[i], "--minima") == 0 && i + 1 < argc) { minima_path = argv[++i]; with_phi = true; }
            else if (!from_tsv && strcmp(argv[i], "--lcp-bin") == 0 && i + 1 < argc) { lcp_bin_path = argv[++i]; with_phi = true; }
            else if (!from_tsv && strcmp(argv[i], "--phi") == 0 && i + 1 < argc) { phi_path = argv[++i]; with_phi = true; }
            else if (!from_tsv && strcmp(argv[i], "--seq-ids") == 0 && i + 1 < argc) seq_ids_path = argv[++i];
            else { std::cerr << "Unknown " << cmd << " option: " << argv[i] << "\n"; return 1; }
        }
        if (int(!minima_path.empty()) + int(!lcp_bin_path.empty()) + int(!phi_path.empty()) > 1) {
            std::cerr << "give at most one of --phi, --minima and --lcp-bin\n";
            return 1;
        }
        if (opts.layout == TmsLayout::PSI && with_phi) {
            std::cerr << "--layout psi takes no LCP input (" << (from_tsv ? "--phi" : "--phi, --minima or --lcp-bin") << ")\n";
            return 1;
        }
        if (opts.layout == TmsLayout::PHI && !with_phi) {
            std::cerr << "--layout phi needs LCP input (" << (from_tsv ? "--phi" : "--phi, --minima or --lcp-bin") << ")\n";
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
                  << ", load_s=" << load_s << ", build_s=" << build_s << ")\n";
        return 0;
    }

    if (cmd == "ftab-build") {
        if (argc < 2) { std::cerr << "ftab-build requires INDEX_PATH and FTAB_PATH\n"; return 1; }
        const std::string idx_path = argv[0], out_path = argv[1];
        unsigned k = 0;  // 0: --auto
        size_t threads = 0;
        for (int i = 2; i < argc; ++i) {
            if (strcmp(argv[i], "--k") == 0 && i + 1 < argc) k = static_cast<unsigned>(std::stoul(argv[++i]));
            else if (strcmp(argv[i], "--auto") == 0) k = 0;
            else if ((strcmp(argv[i], "--threads") == 0 || strcmp(argv[i], "-p") == 0) && i + 1 < argc)
                threads = static_cast<size_t>(std::stoull(argv[++i]));
            else { std::cerr << "Unknown ftab-build option: " << argv[i] << "\n"; return 1; }
        }
        if (k > ftab::max_k) { std::cerr << "--k must be between 1 and 31\n"; return 1; }
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
        if (argc < 2) { std::cerr << "ftab-verify requires INDEX_PATH and FTAB_PATH\n"; return 1; }
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
        if (argc < 1) { std::cerr << "ftab-info requires FTAB_PATH\n"; return 1; }
        try {
            ftab::describe(ftab::read_header_file(argv[0]), std::cout);
        } catch (const std::exception& e) {
            std::cerr << e.what() << "\n";
            return 1;
        }
        return 0;
    }

    if (cmd == "ftab-survey") {
        if (argc < 1) { std::cerr << "ftab-survey requires INDEX_PATH\n"; return 1; }
        unsigned kmin = 6, kmax = 15, verify = 0;
        size_t threads = 1;
        for (int i = 1; i < argc; ++i) {
            if (strcmp(argv[i], "--kmin") == 0 && i + 1 < argc) kmin = static_cast<unsigned>(std::stoul(argv[++i]));
            else if (strcmp(argv[i], "--kmax") == 0 && i + 1 < argc) kmax = static_cast<unsigned>(std::stoul(argv[++i]));
            else if ((strcmp(argv[i], "--threads") == 0 || strcmp(argv[i], "-p") == 0) && i + 1 < argc)
                threads = static_cast<size_t>(std::stoull(argv[++i]));
            else if (strcmp(argv[i], "--verify") == 0 && i + 1 < argc) verify = static_cast<unsigned>(std::stoul(argv[++i]));
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
        if (argc < 1) { std::cerr << "index-plan requires INDEX_PATH\n"; return 1; }
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
        // What this build loads for each tms-batch command line, from the
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
        const char* modes[] = {"psi", "phi", "phiskip", "dual", "oracle", "oracle-scan"};
        const TmsMode mode_of[] = {TmsMode::PSI,  TmsMode::PHI,    TmsMode::PHISKIP,
                                   TmsMode::DUAL, TmsMode::ORACLE, TmsMode::ORACLE_SCAN};
        const char* reports[] = {"ms", "smem-one", "smem-all"};
        const TmsReport report_of[] = {TmsReport::MS, TmsReport::SMEM_ONE, TmsReport::SMEM_ALL};
        if (query) {
            const TmsParts need = tms_batch_parts();
            if (!serves(need)) { std::cout << "this query: refused, the index lacks a part it needs\n"; return 0; }
            std::cout << "this query loads (" << TmsIndex::load_name(need) << "):" << describe_plan(need) << "\n";
            return 0;
        }
        std::vector<std::string> names, plans;
        std::vector<std::string> rows;
        for (int m = 0; m < 6; ++m)
            for (int r = 0; r < 3; ++r)
                for (int pos = 0; pos < 2; ++pos)
                    for (int ps = 0; ps < 2; ++ps) {
                        g_tms_mode = mode_of[m];
                        g_tms_report = report_of[r];
                        g_tms_positions = pos;
                        g_tms_phi_starts = ps;
                        const TmsParts need = tms_batch_parts();
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
        std::cout << "tms-batch queries (--mode, --report, then --positions and --phi-starts):\n";
        for (const auto& r : rows) std::cout << r << "\n";
        return 0;
    }

    if (cmd == "index-verify") {
        if (argc < 1) { std::cerr << "index-verify requires INDEX_PATH\n"; return 1; }
        const std::string bad = orbit_idx::verify_file(argv[0]);
        if (!bad.empty()) { std::cerr << argv[0] << ": " << bad << "\n"; return 1; }
        std::cout << argv[0] << ": every section's CRC-32 matches\n";
        return 0;
    }

    if (cmd == "index-add-seq-ids") {
        if (argc < 3) { std::cerr << "index-add-seq-ids requires INDEX_PATH, SEQ_IDS and OUT_PATH\n"; return 1; }
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
        if (argc < 2) { std::cerr << "index-rewrap requires LEGACY_PATH and INDEX_PATH\n"; return 1; }
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
        if (argc < 1) { std::cerr << "split-report requires INDEX_PATH\n"; return 1; }
        size_t threads = 1;
        for (int i = 1; i < argc; ++i) {
            if ((strcmp(argv[i], "--threads") == 0 || strcmp(argv[i], "-p") == 0) && i + 1 < argc)
                threads = static_cast<size_t>(std::stoull(argv[++i]));
            else { std::cerr << "Unknown split-report option: " << argv[i] << "\n"; return 1; }
        }
        return split_report::run(argv[0], resolve_thread_count(threads), std::cout);
    }

    if (cmd == "tms-inspect") {
        if (argc < 1) { std::cerr << "tms-inspect requires INDEX_PATH\n"; return 1; }
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
        if (argc < 1) { std::cerr << "tms-text requires INDEX_PATH\n"; return 1; }
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
