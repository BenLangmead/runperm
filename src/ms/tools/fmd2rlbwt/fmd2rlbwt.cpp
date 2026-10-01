/**
 * fmd2rlbwt: write the maximal runs of a ropebwt3 FMD as an RLBWT in the
 * .bwt.heads / .bwt.len form that Orbit's ms and tms builds read, and
 * optionally renumber a TeraLCP -f fmd -ominima file to the same runs.
 *
 * The FMD is read with ropebwt3's rld0 code through a read-only memory map,
 * one pass, single-threaded; pages are released once decoded, so resident
 * memory stays small.  The run and head-byte conventions are in FMD_RUNS.md
 * next to this file.
 */

#include "rld0.h"

#include <sys/mman.h>
#include <unistd.h>

#include <array>
#include <cerrno>
#include <chrono>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace {

// Head byte written for each ropebwt3 symbol $ A C G T N (codes 0 to 5).
// Byte order equals ropebwt3's symbol order, so byte-collating readers
// derive the same LF; N is 0x7f, after T, and is not a byte a read holds.
constexpr std::array<uint8_t, 6> HEAD_BYTE = {'$', 'A', 'C', 'G', 'T', 0x7f};
constexpr const char* SYMBOL_NAME[6] = {"$", "A", "C", "G", "T", "N"};
constexpr uint64_t MAX_LEN = (uint64_t{1} << 40) - 1;  // 5-byte lengths

[[noreturn]] void die(const std::string& msg) {
    std::fprintf(stderr, "fmd2rlbwt: %s\n", msg.c_str());
    std::exit(1);
}

/** Buffered binary output. */
class Out {
    FILE* fp_ = nullptr;
    std::vector<uint8_t> buf_;
    std::string path_;

public:
    explicit Out(const std::string& path) : buf_(), path_(path) {
        fp_ = std::fopen(path.c_str(), "wb");
        if (!fp_) die("cannot open " + path + " for writing: " + std::strerror(errno));
        buf_.reserve(1 << 22);
    }
    void put(uint8_t b) {
        buf_.push_back(b);
        if (buf_.size() >= (1 << 22)) flush();
    }
    void put_le(uint64_t v, int bytes) {
        for (int i = 0; i < bytes; ++i) put(static_cast<uint8_t>(v >> (8 * i)));
    }
    void put_uleb(uint64_t v) {
        do {
            uint8_t b = v & 0x7f;
            v >>= 7;
            if (v) b |= 0x80;
            put(b);
        } while (v);
    }
    void flush() {
        if (!buf_.empty() && std::fwrite(buf_.data(), 1, buf_.size(), fp_) != buf_.size())
            die("short write to " + path_);
        buf_.clear();
    }
    /** Overwrite bytes at the start of the file (after flushing). */
    void rewrite_head(const uint8_t* p, size_t len) {
        flush();
        if (std::fseek(fp_, 0, SEEK_SET) != 0 || std::fwrite(p, 1, len, fp_) != len)
            die("cannot rewrite header of " + path_);
    }
    void close() {
        flush();
        if (std::fclose(fp_) != 0) die("error closing " + path_);
        fp_ = nullptr;
    }
};

/** Buffered reader of a TLMSM v1 minima file. */
class MinimaIn {
    FILE* fp_ = nullptr;
    std::vector<uint8_t> buf_;
    size_t have_ = 0, pos_ = 0;
    std::string path_;

public:
    uint64_t r = 0, n = 0, pairs = 0, input_runs = 0;
    unsigned flags = 0;

    explicit MinimaIn(const std::string& path) : buf_(1 << 22), path_(path) {
        fp_ = std::fopen(path.c_str(), "rb");
        if (!fp_) die("cannot open " + path);
        uint8_t h[48];
        if (std::fread(h, 1, 48, fp_) != 48) die("short header in " + path);
        static const uint8_t magic[8] = {0x93, 'T', 'L', 'M', 'S', 'M', 0x00, 0x01};
        if (std::memcmp(h, magic, 8) != 0) die(path + " is not a TLMSM v1 minima file");
        auto u64 = [&](int at) {
            uint64_t v = 0;
            for (int b = 7; b >= 0; --b) v = (v << 8) | h[at + b];
            return v;
        };
        flags = h[8];
        r = u64(16);
        n = u64(24);
        pairs = u64(32);
        input_runs = u64(40);
    }
    bool byte(uint8_t& c) {
        if (pos_ == have_) {
            have_ = std::fread(buf_.data(), 1, buf_.size(), fp_);
            pos_ = 0;
            if (have_ == 0) return false;
        }
        c = buf_[pos_++];
        return true;
    }
    uint64_t uleb() {
        uint64_t v = 0;
        unsigned shift = 0;
        uint8_t c;
        do {
            if (!byte(c) || shift > 63) die("truncated or malformed minima file " + path_);
            v |= static_cast<uint64_t>(c & 0x7f) << shift;
            shift += 7;
        } while (c & 0x80);
        return v;
    }
    bool at_end() {
        uint8_t c;
        return !byte(c);
    }
};

/**
 * Renumbers minima records from TeraLCP's -f fmd numbering (one record per
 * rld run, with every sentinel run split into length-1 runs) to maximal
 * runs.  Pieces arrive in BWT order with their lengths; a maximal run's
 * record is written once the next maximal run's top is known.
 *
 * The merged run's interior values that can be prefix or suffix minima are
 * among its pieces' stored values and the tops of its pieces after the
 * first: a merged prefix (suffix) minimum inside a piece is also below the
 * piece's top (the next piece's top) and below the piece's other values
 * before (after) it, so the piece stored it.  The minima are then
 * recomputed over these candidates in offset order, which gives exactly the
 * merged run's prefix and suffix minima (a value a candidate hides behind is
 * itself no smaller than some candidate on the same side).
 */
class MinimaMerger {
    MinimaIn in_;
    Out out_;
    uint64_t pieces_read_ = 0, runs_written_ = 0, pairs_written_ = 0;
    bool open_ = false;
    uint64_t top_ = 0, run_len_ = 0;
    std::vector<std::pair<uint64_t, uint64_t>> cand_;  // (offset, value), increasing offset
    std::vector<uint8_t> keep_;

    void emit(uint64_t next_top) {
        const size_t k = cand_.size();
        keep_.assign(k, 0);
        uint64_t run_min = top_;
        for (size_t i = 0; i < k; ++i)
            if (cand_[i].second < run_min) { keep_[i] = 1; run_min = cand_[i].second; }
        run_min = next_top;
        for (size_t i = k; i-- > 0;)
            if (cand_[i].second < run_min) { keep_[i] = 1; run_min = cand_[i].second; }
        uint64_t count = 0;
        for (size_t i = 0; i < k; ++i) count += keep_[i];
        out_.put_uleb(top_);
        out_.put_uleb(count);
        uint64_t prev = 0;
        for (size_t i = 0; i < k; ++i) {
            if (!keep_[i]) continue;
            out_.put_uleb(cand_[i].first - prev);
            out_.put_uleb(cand_[i].second);
            prev = cand_[i].first;
        }
        pairs_written_ += count;
        ++runs_written_;
    }

public:
    MinimaMerger(const std::string& in_path, const std::string& out_path, uint64_t n)
        : in_(in_path), out_(out_path) {
        if (in_.n != n)
            die("minima file has n = " + std::to_string(in_.n) + " but the FMD has n = " + std::to_string(n));
        uint8_t zero[48] = {};
        for (uint8_t b : zero) out_.put(b);
    }

    /** One piece in TeraLCP's numbering; starts_run if it begins a maximal run. */
    void piece(uint64_t len, bool starts_run) {
        if (pieces_read_ == in_.r)
            die("minima file has " + std::to_string(in_.r) + " records but the FMD has more runs in TeraLCP's numbering");
        ++pieces_read_;
        const uint64_t top = in_.uleb(), count = in_.uleb();
        if (starts_run) {
            if (open_) emit(top);
            open_ = true;
            top_ = top;
            run_len_ = 0;
            cand_.clear();
        } else {
            cand_.emplace_back(run_len_, top);
        }
        uint64_t off = 0;
        for (uint64_t j = 0; j < count; ++j) {
            const uint64_t gap = in_.uleb(), val = in_.uleb();
            off += gap;
            if (gap == 0 || off >= len) die("bad pair in minima record " + std::to_string(pieces_read_ - 1));
            cand_.emplace_back(run_len_ + off, val);
        }
        run_len_ += len;
    }

    void finish(uint64_t n, uint64_t runs) {
        if (open_) emit(std::numeric_limits<uint64_t>::max());
        if (pieces_read_ != in_.r)
            die("minima file has " + std::to_string(in_.r) + " records but the FMD has " +
                std::to_string(pieces_read_) + " runs in TeraLCP's numbering");
        if (!in_.at_end()) die("trailing bytes after the last minima record");
        if (runs_written_ != runs) die("internal error: merged run count differs from the RLBWT's");
        uint8_t h[48] = {0x93, 'T', 'L', 'M', 'S', 'M', 0x00, 0x01};
        h[8] = 3;  // little-endian, input_runs filled in, numbering equals the RLBWT's
        auto put64 = [&](int at, uint64_t v) { for (int b = 0; b < 8; ++b) h[at + b] = (v >> (8 * b)) & 0xff; };
        put64(16, runs);
        put64(24, n);
        put64(32, pairs_written_);
        put64(40, runs);
        out_.rewrite_head(h, 48);
        out_.close();
        std::fprintf(stderr, "minima: %" PRIu64 " records in (flags %u, input_runs %" PRIu64 "), %" PRIu64
                     " records out, %" PRIu64 " pairs in, %" PRIu64 " pairs out\n",
                     in_.r, in_.flags, in_.input_runs, runs_written_, in_.pairs, pairs_written_);
    }
};

void usage() {
    std::fprintf(stderr,
                 "Usage: fmd2rlbwt [options] IN.fmd OUT_BASE\n"
                 "\n"
                 "Writes the maximal runs of a ropebwt3 FMD (symbols $ACGTN) to\n"
                 "OUT_BASE.bwt.heads (one byte per run: $ A C G T, and 0x7f for N) and\n"
                 "OUT_BASE.bwt.len (5-byte little-endian lengths).  See FMD_RUNS.md.\n"
                 "\n"
                 "Options:\n"
                 "  --minima-in FILE   a TeraLCP -f fmd -ominima file for the same FMD, whose\n"
                 "                     records follow TeraLCP's numbering (rld runs, sentinel\n"
                 "                     runs split); requires --minima-out\n"
                 "  --minima-out FILE  write it renumbered to the maximal runs\n"
                 "  --no-release       keep decoded pages mapped (by default they are\n"
                 "                     released as the scan passes them)\n");
}

}  // namespace

int main(int argc, char** argv) {
    std::string minima_in, minima_out;
    bool release = true;
    std::vector<std::string> pos;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--minima-in" && i + 1 < argc) minima_in = argv[++i];
        else if (a == "--minima-out" && i + 1 < argc) minima_out = argv[++i];
        else if (a == "--no-release") release = false;
        else if (a == "-h" || a == "--help") { usage(); return 0; }
        else if (!a.empty() && a[0] == '-') { usage(); die("unknown option " + a); }
        else pos.push_back(a);
    }
    if (pos.size() != 2) { usage(); return 1; }
    if (minima_in.empty() != minima_out.empty()) die("--minima-in and --minima-out go together");
    const auto t0 = std::chrono::steady_clock::now();

    rld_t* e = rld_restore_mmap(pos[0].c_str());
    if (!e || !e->mem || e->mem == MAP_FAILED) die("cannot load FMD " + pos[0]);
    if (e->asize != 6) die("expected a 6-symbol ($ACGTN) FMD, found asize " + std::to_string(e->asize));
    const uint64_t n = e->mcnt[0];
    madvise(e->mem, ((4 + e->asize) * 8) + e->n_bytes, MADV_SEQUENTIAL);

    Out heads(pos[1] + ".bwt.heads"), lens(pos[1] + ".bwt.len");
    MinimaMerger* merger = minima_in.empty() ? nullptr : new MinimaMerger(minima_in, minima_out, n);

    const long page = sysconf(_SC_PAGESIZE);
    // Release the pages of large block z[b] once the scan has moved past it.
    auto release_block = [&](uint64_t** blk) {
        if (!release) return;
        const uintptr_t lo = reinterpret_cast<uintptr_t>(*blk);
        const uintptr_t hi = lo + static_cast<uintptr_t>(RLD_LSIZE) * 8;
        const uintptr_t a = (lo + page - 1) / page * page, b = hi / page * page;
        if (b > a) madvise(reinterpret_cast<void*>(a), b - a, MADV_DONTNEED);
    };

    std::array<uint64_t, 6> sym_total{}, sym_raw{}, sym_max{};
    uint64_t raw = 0, teralcp = 0, runs = 0, total = 0, max_len = 0;
    int cur_c = -1;
    uint64_t cur_len = 0;
    auto flush_run = [&]() {
        if (cur_c < 0) return;
        if (cur_len > MAX_LEN) die("a run of length " + std::to_string(cur_len) + " does not fit 5 bytes");
        heads.put(HEAD_BYTE[cur_c]);
        lens.put_le(cur_len, 5);
        ++runs;
        ++sym_max[cur_c];
        if (cur_len > max_len) max_len = cur_len;
    };

    rlditr_t itr;
    rld_itr_init(e, &itr, 0);
    uint64_t** blk = itr.i;
    int c = 0;
    int64_t l;
    while ((l = rld_dec(e, &itr, &c, 0)) > 0) {
        if (itr.i != blk) {
            release_block(blk);
            blk = itr.i;
        }
        if (c < 0 || c > 5) die("unexpected symbol code " + std::to_string(c));
        ++raw;
        ++sym_raw[c];
        sym_total[c] += static_cast<uint64_t>(l);
        total += static_cast<uint64_t>(l);
        const bool starts = c != cur_c;
        if (merger) {
            // TeraLCP gives each sentinel its own run.
            if (c == 0) {
                for (int64_t k = 0; k < l; ++k) merger->piece(1, starts && k == 0);
                teralcp += static_cast<uint64_t>(l);
            } else {
                merger->piece(static_cast<uint64_t>(l), starts);
                ++teralcp;
            }
        } else {
            teralcp += c == 0 ? static_cast<uint64_t>(l) : 1;
        }
        if (starts) {
            flush_run();
            cur_c = c;
            cur_len = 0;
        }
        cur_len += static_cast<uint64_t>(l);
    }
    flush_run();
    heads.close();
    lens.close();
    if (total != n) die("decoded " + std::to_string(total) + " symbols but the FMD header says " + std::to_string(n));
    for (int k = 0; k < 6; ++k)
        if (sym_total[k] != e->mcnt[k + 1]) die(std::string("symbol count mismatch for ") + SYMBOL_NAME[k]);
    if (merger) {
        merger->finish(n, runs);
        delete merger;
    }

    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("n\t%" PRIu64 "\n", n);
    std::printf("sequences\t%" PRIu64 "\n", sym_total[0]);
    std::printf("rld_runs\t%" PRIu64 "\n", raw);
    std::printf("teralcp_runs\t%" PRIu64 "\n", teralcp);
    std::printf("maximal_runs\t%" PRIu64 "\n", runs);
    std::printf("max_run_length\t%" PRIu64 "\n", max_len);
    for (int k = 0; k < 6; ++k)
        std::printf("symbol\t%s\tcount\t%" PRIu64 "\trld_runs\t%" PRIu64 "\tmaximal_runs\t%" PRIu64 "\n",
                    SYMBOL_NAME[k], sym_total[k], sym_raw[k], sym_max[k]);
    std::printf("seconds\t%.2f\n", secs);
    return 0;
}
