/**
 * Read a run-length BWT stored as two files: HEADS, one byte per run, and
 * LENS, one little-endian unsigned integer per run of a fixed width (5 bytes
 * in the Movi pipeline's form; the width is inferred from the file sizes).
 *
 * Head bytes are mapped to the codes the index uses: the end-of-text sentinel
 * (0x00 in the Movi form, 0x0a in the TeraLCP form, or # or $) becomes
 * orbit::TERMINATOR,
 * the sequence separator '%' becomes orbit::SEPARATOR, and A, C, G, T are kept.
 * This preserves the sort order sentinel < separator < A < C < G < T that the
 * BWT was built with.
 */

#ifndef _MS_RLBWT_IO_HPP
#define _MS_RLBWT_IO_HPP

#include "orbit/common.hpp"
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

namespace rlbwt_io {

using uchar = orbit::uchar;
using ulint = orbit::ulint;

/** Map one RLBWT head byte to an index code; returns false for unsupported bytes. */
inline bool map_head_byte(uchar b, uchar& code) {
    switch (b) {
        case 0x00: case 0x0a: case '#': case '$': code = orbit::TERMINATOR; return true;
        case '%': code = orbit::SEPARATOR; return true;
        case 'A': case 'C': case 'G': case 'T': code = b; return true;
        default: return false;
    }
}

/**
 * Load HEADS and LENS.  On failure returns false and sets err.  Adjacent runs
 * with the same head are left as they are; the caller's LCP data must use the
 * same run numbering.
 */
inline bool load_rlbwt(const std::string& heads_path, const std::string& lens_path,
                       std::vector<uchar>& heads, std::vector<ulint>& lens, std::string& err) {
    std::ifstream hin(heads_path, std::ios::binary | std::ios::ate);
    std::ifstream lin(lens_path, std::ios::binary | std::ios::ate);
    if (!hin.good()) { err = "cannot open " + heads_path; return false; }
    if (!lin.good()) { err = "cannot open " + lens_path; return false; }
    const size_t r = static_cast<size_t>(hin.tellg());
    const size_t lbytes = static_cast<size_t>(lin.tellg());
    if (r == 0) { err = "empty heads file"; return false; }
    if (lbytes % r != 0 || lbytes / r == 0 || lbytes / r > 8) {
        err = "lens file size " + std::to_string(lbytes) + " is not a whole number (1 to 8) of bytes per run for " +
              std::to_string(r) + " runs";
        return false;
    }
    const size_t width = lbytes / r;
    hin.seekg(0);
    lin.seekg(0);
    heads.resize(r);
    lens.resize(r);
    hin.read(reinterpret_cast<char*>(heads.data()), static_cast<std::streamsize>(r));
    std::vector<uchar> buf(width * 65536);
    size_t done = 0, sentinels = 0;
    while (done < r) {
        const size_t chunk = std::min<size_t>(65536, r - done);
        lin.read(reinterpret_cast<char*>(buf.data()), static_cast<std::streamsize>(chunk * width));
        for (size_t i = 0; i < chunk; ++i) {
            ulint v = 0;
            for (size_t b = 0; b < width; ++b) v |= static_cast<ulint>(buf[i * width + b]) << (8 * b);
            lens[done + i] = v;
        }
        done += chunk;
    }
    if (!hin.good() || !lin.good()) { err = "short read"; return false; }
    for (size_t i = 0; i < r; ++i) {
        uchar code;
        if (!map_head_byte(heads[i], code)) {
            char hex[8];
            std::snprintf(hex, sizeof(hex), "0x%02x", heads[i]);
            err = "unsupported head byte " + std::string(hex) + " at run " + std::to_string(i);
            return false;
        }
        if (code == orbit::TERMINATOR) sentinels += lens[i];
        heads[i] = code;
        if (lens[i] == 0) { err = "zero-length run " + std::to_string(i); return false; }
    }
    if (sentinels != 1) {
        err = "expected exactly one sentinel position, found " + std::to_string(sentinels);
        return false;
    }
    return true;
}

/**
 * Read a TeraLCP -ominima file (format TLMSM v1, documented in TeraTools
 * src/TeraLCP/MS_MINIMA_FORMAT.md) into per-run pairs: element 0 of each run
 * is (0, top) and the rest are the stored interior (offset, value) pairs in
 * increasing offset order.  expected_runs and expected_n are checked against
 * the header.
 */
inline bool read_minima(const std::string& path, size_t expected_runs, ulint expected_n,
                        std::vector<std::vector<std::pair<ulint, ulint>>>& runs, std::string& err) {
    std::ifstream in(path, std::ios::binary);
    if (!in.good()) { err = "cannot open " + path; return false; }
    std::array<unsigned char, 48> h{};
    in.read(reinterpret_cast<char*>(h.data()), 48);
    if (!in.good()) { err = "short header in " + path; return false; }
    static const unsigned char magic[8] = {0x93, 'T', 'L', 'M', 'S', 'M', 0x00, 0x01};
    if (std::memcmp(h.data(), magic, 8) != 0) { err = path + " is not a TLMSM v1 minima file"; return false; }
    auto u64 = [&](size_t at) { ulint v = 0; for (int b = 7; b >= 0; --b) v = (v << 8) | h[at + b]; return v; };
    const unsigned flags = h[8];
    const ulint r = u64(16), n = u64(24), total_pairs = u64(32);
    if (flags & 4) {
        err = "minima file run numbering differs from its input's (sentinel runs were split); "
              "use heads/lens written by TeraLCP for the same input";
        return false;
    }
    if (r != expected_runs) {
        err = "minima file has " + std::to_string(r) + " runs but the RLBWT has " + std::to_string(expected_runs);
        return false;
    }
    if (n != expected_n) {
        err = "minima file has n = " + std::to_string(n) + " but the RLBWT has n = " + std::to_string(expected_n);
        return false;
    }
    std::vector<unsigned char> buf(1 << 20);
    size_t have = 0, pos = 0;
    bool eof = false;
    auto byte = [&](unsigned char& c) {
        if (pos == have) {
            if (eof) return false;
            in.read(reinterpret_cast<char*>(buf.data()), static_cast<std::streamsize>(buf.size()));
            have = static_cast<size_t>(in.gcount());
            pos = 0;
            if (have < buf.size()) eof = true;
            if (have == 0) return false;
        }
        c = buf[pos++];
        return true;
    };
    auto uleb = [&](ulint& v) {
        v = 0;
        unsigned shift = 0;
        unsigned char c;
        do {
            if (!byte(c) || shift > 63) return false;
            v |= static_cast<ulint>(c & 0x7f) << shift;
            shift += 7;
        } while (c & 0x80);
        return true;
    };
    runs.assign(static_cast<size_t>(r), {});
    ulint pairs_seen = 0;
    for (size_t i = 0; i < r; ++i) {
        ulint top, count;
        if (!uleb(top) || !uleb(count)) { err = "truncated record " + std::to_string(i); return false; }
        auto& run = runs[i];
        run.reserve(static_cast<size_t>(count) + 1);
        run.emplace_back(0, top);
        ulint off = 0;
        for (ulint k = 0; k < count; ++k) {
            ulint gap, val;
            if (!uleb(gap) || !uleb(val) || gap == 0) { err = "bad pair in record " + std::to_string(i); return false; }
            off += gap;
            run.emplace_back(off, val);
        }
        pairs_seen += count;
    }
    unsigned char extra;
    if (byte(extra)) { err = "trailing bytes after the last record"; return false; }
    if (pairs_seen != total_pairs) { err = "pair count does not match the header"; return false; }
    return true;
}

/**
 * Read a file of one little-endian 64-bit LCP value per BWT row, optionally
 * after a 64-bit row count (msbench prep's form), and call f(i, values) for
 * each run i in order, with values holding its rows' LCPs.  lens gives the
 * runs' lengths.  On failure returns false and sets err.
 */
template <typename F>
inline bool for_each_run_lcps(const std::string& path, const std::vector<ulint>& lens, F f, std::string& err) {
    std::ifstream in(path, std::ios::binary | std::ios::ate);
    if (!in.good()) { err = "cannot open " + path; return false; }
    ulint n = 0;
    for (ulint l : lens) n += l;
    const ulint bytes = static_cast<ulint>(in.tellg());
    in.seekg(0);
    if (bytes == 8 * (n + 1)) {
        uint64_t count = 0;
        in.read(reinterpret_cast<char*>(&count), 8);
        if (count != n) { err = "lcp file count " + std::to_string(count) + " is not n = " + std::to_string(n); return false; }
    } else if (bytes != 8 * n) {
        err = "lcp file has " + std::to_string(bytes) + " bytes; expected 8 per row for n = " + std::to_string(n);
        return false;
    }
    std::vector<uint64_t> buf(1 << 20);
    size_t have = 0, at = 0;
    std::vector<ulint> values;
    for (size_t i = 0; i < lens.size(); ++i) {
        values.resize(static_cast<size_t>(lens[i]));
        for (auto& v : values) {
            if (at == have) {
                in.read(reinterpret_cast<char*>(buf.data()), static_cast<std::streamsize>(buf.size() * 8));
                have = static_cast<size_t>(in.gcount()) / 8;
                at = 0;
                if (have == 0) { err = "lcp file is shorter than the RLBWT"; return false; }
            }
            v = buf[at++];
        }
        f(i, values);
    }
    return true;
}

/** Write per-run pairs as a TLMSM v1 minima file (used by tests). */
inline bool write_minima(const std::string& path, const std::vector<std::vector<std::pair<ulint, ulint>>>& runs,
                         ulint n) {
    std::ofstream out(path, std::ios::binary);
    if (!out.good()) return false;
    std::array<unsigned char, 48> h{};
    static const unsigned char magic[8] = {0x93, 'T', 'L', 'M', 'S', 'M', 0x00, 0x01};
    std::memcpy(h.data(), magic, 8);
    h[8] = 3;  // little-endian, input_runs filled in
    ulint total = 0;
    for (const auto& run : runs) total += run.size() - 1;
    auto put64 = [&](size_t at, ulint v) { for (int b = 0; b < 8; ++b) h[at + b] = (v >> (8 * b)) & 0xff; };
    put64(16, runs.size()); put64(24, n); put64(32, total); put64(40, runs.size());
    out.write(reinterpret_cast<const char*>(h.data()), 48);
    auto uleb = [&](ulint v) { do { unsigned char b = v & 0x7f; v >>= 7; if (v) b |= 0x80; out.put(static_cast<char>(b)); } while (v); };
    for (const auto& run : runs) {
        uleb(run[0].second);
        uleb(run.size() - 1);
        ulint prev = 0;
        for (size_t k = 1; k < run.size(); ++k) { uleb(run[k].first - prev); uleb(run[k].second); prev = run[k].first; }
    }
    return out.good();
}

/** Write HEADS and LENS in the same form (used by tests). */
inline bool write_rlbwt(const std::string& heads_path, const std::string& lens_path,
                        const std::vector<uchar>& raw_heads, const std::vector<ulint>& lens, size_t width = 5) {
    std::ofstream hout(heads_path, std::ios::binary), lout(lens_path, std::ios::binary);
    if (!hout.good() || !lout.good()) return false;
    hout.write(reinterpret_cast<const char*>(raw_heads.data()), static_cast<std::streamsize>(raw_heads.size()));
    for (ulint v : lens)
        for (size_t b = 0; b < width; ++b) lout.put(static_cast<char>((v >> (8 * b)) & 0xff));
    return hout.good() && lout.good();
}

}  // namespace rlbwt_io

#endif /* _MS_RLBWT_IO_HPP */
