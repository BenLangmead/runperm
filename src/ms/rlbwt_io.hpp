/**
 * Read a run-length BWT stored as two files: HEADS, one byte per run, and
 * LENS, one little-endian unsigned integer per run of a fixed width (5 bytes
 * in the Movi pipeline's form; the width is inferred from the file sizes).
 *
 * Head bytes are mapped to the codes the index uses: the end-of-text sentinel
 * (0x00 in the Movi form, 0x0a in the TeraLCP form) becomes orbit::TERMINATOR,
 * the sequence separator '%' becomes orbit::SEPARATOR, and A, C, G, T are kept.
 * This preserves the sort order sentinel < separator < A < C < G < T that the
 * BWT was built with.
 */

#ifndef _MS_RLBWT_IO_HPP
#define _MS_RLBWT_IO_HPP

#include "orbit/common.hpp"
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

namespace rlbwt_io {

using uchar = orbit::uchar;
using ulint = orbit::ulint;

/** Map one RLBWT head byte to an index code; returns false for unsupported bytes. */
inline bool map_head_byte(uchar b, uchar& code) {
    switch (b) {
        case 0x00: case 0x0a: code = orbit::TERMINATOR; return true;
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
