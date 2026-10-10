/**
 * An ftab: a dense table, indexed by every k-mer over ACGT, of the state a
 * matching statistics query reaches after reading that k-mer from a cold
 * start (idx.first() with match length 0).  A query whose read ends in a
 * k-mer that occurs in the text takes its state from the table and starts
 * k steps in.  The statistics in that window are implied: a k-mer that
 * occurs has every suffix occurring, so the statistic d characters from the
 * read's right end is d + 1.
 *
 * An entry holds a bit that is set when the k-mer occurs, the interval and
 * offset of the resolved LF position after the walk, and, in a table with
 * a toehold, the resolved phi point (interval and offset) of that row's
 * text position, which a tms query in a mode that walks phi, or with
 * positions, carries along.  phi's two forms share their intervals, so the
 * toehold serves either.  Entries are bit packed.
 *
 * A tms query with positions also reports, for each statistic in the
 * window, the text position its own walk reached there.  The walk over a
 * k-mer passes through the state of each of its suffixes, so a table with
 * window positions holds, for every j-mer with j < k, the text position its
 * walk ends on; the k-mer's own is its toehold's.
 *
 * Tables are built by ftab_build.hpp and stored in files of their own
 * (ftab_file.hpp, FTAB_FORMAT.md) next to the index they were built from.
 */

#ifndef _FTAB_HPP
#define _FTAB_HPP

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace ftab {

/** Number of bits needed to write any value in [0, v]. */
inline unsigned bits_for(uint64_t v) {
    unsigned b = 0;
    while (v > 0) { ++b; v >>= 1; }
    return b;
}

/** A k-mer's key, or no_key for a window that holds a character other than ACGT. */
constexpr uint64_t no_key = ~uint64_t{0};

/** The largest k a table supports. */
constexpr unsigned max_k = 31;

/** What a table holds beyond the LF position, and the bit width of each field. */
struct Shape {
    unsigned k = 0;
    unsigned interval = 0, offset = 0;                  // LF position
    bool toehold = false;
    unsigned toehold_interval = 0, toehold_offset = 0;  // phi point
    bool window = false;                                // window positions (needs toehold)
    unsigned position = 0;                              // a window position's bits

    unsigned entry_bits() const { return 1 + interval + offset + (toehold ? toehold_interval + toehold_offset : 0); }
    /** Words of the entries, with one spare so that reads need not check the end. */
    uint64_t table_words() const { return ((uint64_t{1} << (2 * k)) * entry_bits() + 63) / 64 + 1; }
    /** Window positions: 4^j for each j from 1 to k - 1. */
    uint64_t window_entries() const { return window && k > 1 ? ((uint64_t{1} << (2 * k)) - 4) / 3 : 0; }
    uint64_t window_words() const { return window && k > 1 ? (window_entries() * position + 63) / 64 + 1 : 0; }
    /** Bytes the table takes in memory and in its file's payloads. */
    uint64_t payload_bytes() const { return 8 * (table_words() + window_words()); }
};

/** A table entry: the state after an occurring k-mer. */
struct Entry {
    uint64_t interval = 0, offset = 0;                  // LF position
    uint64_t toehold_interval = 0, toehold_offset = 0;  // phi point, in a table with a toehold
};

class Ftab {
public:
    Ftab() { init_codes(); }
    /** An empty table (every k-mer absent) of the given shape. */
    explicit Ftab(const Shape& s) : s_(s), eb_(s.entry_bits()) {
        if (s.k < 1 || s.k > max_k) throw std::invalid_argument("ftab k must be between 1 and 31");
        if (s.window && !s.toehold) throw std::invalid_argument("ftab window positions need a toehold");
        for (unsigned w : {s.interval, s.offset, s.toehold_interval, s.toehold_offset, s.position})
            if (w > 64) throw std::invalid_argument("ftab fields are at most 64 bits");
        if (2 * s.k + bits_for(std::max(eb_, s.position)) > 63) throw std::invalid_argument("ftab too large to address");
        words_.assign(s.table_words(), 0);
        win_.assign(s.window_words(), 0);
        init_codes();
    }

    const Shape& shape() const { return s_; }
    unsigned k() const { return s_.k; }
    unsigned entry_bits() const { return eb_; }
    bool has_toehold() const { return s_.toehold; }
    bool has_window() const { return s_.window; }
    size_t bytes() const { return (words_.size() + win_.size()) * sizeof(uint64_t); }

    /**
     * The key of the last k characters of p[0 .. m): the character d places
     * from the right end is in bits 2d and 2d + 1.  no_key if m < k or the
     * window holds a character other than A, C, G and T.
     */
    uint64_t key(const char* p, size_t m) const {
        if (m < s_.k) return no_key;
        uint64_t key = 0;
        for (unsigned d = 0; d < s_.k; ++d) {
            const uint8_t c = code_[static_cast<uint8_t>(p[m - 1 - d])];
            if (c > 3) return no_key;
            key |= uint64_t{c} << (2 * d);
        }
        return key;
    }

    void prefetch(uint64_t key) const {
        if (key != no_key) __builtin_prefetch(&words_[key * eb_ >> 6]);
    }

    /** True, with the entry, if the k-mer with this key occurs. */
    bool lookup(uint64_t key, Entry& e) const {
        if (key == no_key) return false;
        uint64_t bit = key * eb_;
        if (!(get(words_, bit, 1) & 1)) return false;
        ++bit;
        e.interval = get(words_, bit, s_.interval);
        bit += s_.interval;
        e.offset = get(words_, bit, s_.offset);
        bit += s_.offset;
        if (s_.toehold) {
            e.toehold_interval = get(words_, bit, s_.toehold_interval);
            bit += s_.toehold_interval;
            e.toehold_offset = get(words_, bit, s_.toehold_offset);
        }
        return true;
    }

    /** Set the entry of an occurring k-mer.  Safe to call from several threads on different keys. */
    void set(uint64_t key, const Entry& e) {
        uint64_t bit = key * eb_;
        put(words_, bit, 1, 1);
        ++bit;
        put(words_, bit, s_.interval, e.interval);
        bit += s_.interval;
        put(words_, bit, s_.offset, e.offset);
        bit += s_.offset;
        if (s_.toehold) {
            put(words_, bit, s_.toehold_interval, e.toehold_interval);
            bit += s_.toehold_interval;
            put(words_, bit, s_.toehold_offset, e.toehold_offset);
        }
    }

    /**
     * The window position of the j-mer (1 <= j < k) that is the last j
     * characters of the k-mer with this key.
     */
    uint64_t window(unsigned j, uint64_t key) const {
        return get(win_, window_index(j, key) * s_.position, s_.position);
    }
    /** Set the window position of the j-mer with this key.  Safe as set is. */
    void set_window(unsigned j, uint64_t key, uint64_t pos) {
        put(win_, window_index(j, key) * s_.position, s_.position, pos);
    }

    /** Number of occurring k-mers. */
    uint64_t present() const {
        uint64_t n = 0;
        for (uint64_t key = 0; key < (uint64_t{1} << (2 * s_.k)); ++key) n += get(words_, key * eb_, 1);
        return n;
    }

    // The packed words, for reading and writing the table's file.
    std::vector<uint64_t>& table_words() { return words_; }
    const std::vector<uint64_t>& table_words() const { return words_; }
    std::vector<uint64_t>& window_words() { return win_; }
    const std::vector<uint64_t>& window_words() const { return win_; }

private:
    static uint64_t mask(unsigned b) { return b >= 64 ? ~uint64_t{0} : (uint64_t{1} << b) - 1; }
    static uint64_t window_index(unsigned j, uint64_t key) {
        return ((uint64_t{1} << (2 * j)) - 4) / 3 + (key & mask(2 * j));
    }
    static uint64_t get(const std::vector<uint64_t>& w, uint64_t bit, unsigned width) {
        if (width == 0) return 0;
        const unsigned sh = bit & 63;
        uint64_t v = w[bit >> 6] >> sh;
        if (sh + width > 64) v |= w[(bit >> 6) + 1] << (64 - sh);
        return v & mask(width);
    }
    static void put(std::vector<uint64_t>& w, uint64_t bit, unsigned width, uint64_t v) {
        if (width == 0) return;
        v &= mask(width);
        const unsigned sh = bit & 63;
        __atomic_fetch_or(&w[bit >> 6], v << sh, __ATOMIC_RELAXED);
        if (sh + width > 64) __atomic_fetch_or(&w[(bit >> 6) + 1], v >> (64 - sh), __ATOMIC_RELAXED);
    }
    void init_codes() {
        for (int c = 0; c < 256; ++c) code_[c] = 4;
        code_['A'] = 0; code_['C'] = 1; code_['G'] = 2; code_['T'] = 3;
    }

    Shape s_;
    unsigned eb_ = 0;
    std::vector<uint64_t> words_, win_;
    std::array<uint8_t, 256> code_{};
};

}  // namespace ftab

#endif /* _FTAB_HPP */
