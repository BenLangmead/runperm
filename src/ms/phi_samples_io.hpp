/**
 * Reading TeraLCP's -ophi file (TLPHI v1, src/TeraLCP/PHI_EXPORT_FORMAT.md in
 * TeraTools) into TmsPhiSamples.
 *
 * The file is a 48-byte header (magic 0x93 "TLPHI" 0x00 0x01, flags, the
 * field width w in bytes, the record count R, the text length n and the
 * number of maximal runs) followed by R records of three little-endian w-byte
 * fields: the text positions of a run's top and bottom rows and the LCP at
 * its top row.  Records follow the BWT's runs in row order, with one record
 * per row of a terminator run.
 */

#ifndef _MS_PHI_SAMPLES_IO_HPP
#define _MS_PHI_SAMPLES_IO_HPP

#include "tms_index.hpp"
#include <cstdint>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

namespace phi_samples_io {

/**
 * Read path into samples, packed.  n is the BWT length the samples must
 * describe and records the number of records the RLBWT needs
 * (TmsPhiSamples::records_for).  The file is read twice: first for the
 * largest top, which sets the width of the packed tops.  On failure returns
 * false and sets err.
 */
inline bool read_phi_samples(const std::string& path, ulint n, ulint records, TmsPhiSamples& samples,
                             std::string& err) {
    std::ifstream in(path, std::ios::binary);
    if (!in.good()) {
        err = "cannot open " + path;
        return false;
    }
    unsigned char h[48];
    in.read(reinterpret_cast<char*>(h), sizeof(h));
    static const unsigned char MAGIC[8] = {0x93, 'T', 'L', 'P', 'H', 'I', 0x00, 0x01};
    if (!in.good() || std::memcmp(h, MAGIC, 8) != 0) {
        err = path + " is not a TeraLCP -ophi file (TLPHI v1)";
        return false;
    }
    auto u64 = [&](size_t at) {
        uint64_t x = 0;
        for (int k = 7; k >= 0; --k) x = (x << 8) | h[at + k];
        return x;
    };
    const unsigned w = h[9];
    const uint64_t R = u64(16), file_n = u64(24);
    if (!(h[8] & 1) || w == 0 || w > 8) {
        err = path + ": unsupported flags or field width";
        return false;
    }
    if (file_n != n) {
        err = path + " describes a text of length " + std::to_string(file_n) + ", but the RLBWT has n = " +
              std::to_string(n);
        return false;
    }
    if (R != records) {
        err = path + " has " + std::to_string(R) + " records, but the RLBWT needs " + std::to_string(records) +
              " (one per run, and one per row of a terminator run)";
        return false;
    }
    const size_t rec = 3 * w, chunk = 1 << 16;
    std::vector<unsigned char> buf(rec * chunk);
    auto field = [&](const unsigned char* p) {
        uint64_t x = 0;
        for (int k = static_cast<int>(w) - 1; k >= 0; --k) x = (x << 8) | p[k];
        return static_cast<ulint>(x);
    };
    // Call f(i, head, tail, top) for every record, in order.
    auto each = [&](auto f) {
        in.clear();
        in.seekg(48);
        for (uint64_t i = 0; i < R;) {
            const size_t take = static_cast<size_t>(std::min<uint64_t>(chunk, R - i));
            in.read(reinterpret_cast<char*>(buf.data()), static_cast<std::streamsize>(take * rec));
            if (static_cast<size_t>(in.gcount()) != take * rec) {
                err = path + " is truncated at record " + std::to_string(i);
                return false;
            }
            for (size_t j = 0; j < take; ++j, ++i) {
                const unsigned char* p = buf.data() + j * rec;
                const ulint head = field(p), tail = field(p + w), top = field(p + 2 * w);
                if (head >= n || tail >= n || top >= n) {
                    err = path + ": record " + std::to_string(i) + " lies outside the text";
                    return false;
                }
                f(i, head, tail, top);
            }
        }
        return true;
    };
    ulint max_top = 0;
    if (!each([&](uint64_t, ulint, ulint, ulint top) { max_top = std::max(max_top, top); })) return false;
    samples = TmsPhiSamples(R, n, max_top);
    if (!each([&](uint64_t i, ulint head, ulint tail, ulint top) {
            samples.head.set(i, head);
            samples.tail.set(i, tail);
            samples.top.set(i, top);
        }))
        return false;
    if (in.peek() != std::char_traits<char>::eof()) {
        err = path + " has bytes after its last record";
        return false;
    }
    return true;
}

}  // namespace phi_samples_io

#endif /* _MS_PHI_SAMPLES_IO_HPP */
