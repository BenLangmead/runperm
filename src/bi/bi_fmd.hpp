/**
 * FMD-index bidirectional BWT search built atop the runperm library.
 *
 * Move-structure-native design: all navigation uses Position (interval, offset).
 * No run_starts, run_chars, find_run, or rank — LF(first_a) replaces C[a]+rank.
 *
 * Author: Ben Langmead (ben.langmead@gmail.com)
 * Date: Feb 23, 2026
 */

#ifndef _BI_FMD_HPP
#define _BI_FMD_HPP

#include "orbit/rlbwt.hpp"
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <cassert>

using uchar = orbit::uchar;
using ulint = orbit::ulint;

// Nucleotide mapped values (from alphabet.hpp):
//   TER=0, SEP=1, A=2, C=3, G=4, T=5, N=6
static constexpr uchar NUC_TER = 0;
static constexpr uchar NUC_SEP = 1;
static constexpr uchar NUC_A   = 2;
static constexpr uchar NUC_C   = 3;
static constexpr uchar NUC_G   = 4;
static constexpr uchar NUC_T   = 5;
static constexpr uchar NUC_N   = 6;
static constexpr uchar NUC_SIGMA = 7;

/// Complement in "mapped" space.  Only A/C/G/T are meaningful for extension.
static constexpr uchar complement_map[NUC_SIGMA] = {
    NUC_TER,  // TER → TER (unused for extension)
    NUC_SEP,  // SEP → SEP (unused for extension)
    NUC_T,    // A(2) → T(5)
    NUC_G,    // C(3) → G(4)
    NUC_C,    // G(4) → C(3)
    NUC_A,    // T(5) → A(2)
    NUC_N     // N(6) → N  (unused for extension)
};

/** Complement a character in "mapped" space using complement_map */
inline uchar complement(uchar c) {
    assert(c < NUC_SIGMA);
    return complement_map[c];
}

// LF index type: RunPermLF with empty_data_columns for up()/down() run-walking
using LFIndex = orbit::rlbwt::RunPermLF<orbit::empty_data_columns, false, false, orbit::nucleotide>;
using Position = LFIndex::position;

/**
 * A bi-interval is a pair of suffix array intervals, one for the forward strand
 * and one for the reverse-complement strand.  Uses Position (interval, offset)
 * for both interval starts — no flat BWT offsets.
 */
struct BiInterval {
    Position pos;    // start of forward SA interval
    Position pos_r;  // start of reverse-complement SA interval
    ulint size = 0;  // number of occurrences (same for both strands)

    bool empty() const { return size == 0; }

    bool operator==(const BiInterval& o) const {
        return pos == o.pos && pos_r == o.pos_r && size == o.size;
    }
    bool operator!=(const BiInterval& o) const { return !(*this == o); }
};

/**
 * A table used to jump-start the bidirectional matching process for DNA query
 * strings.  For a given length K, the table contains a pre-computed
 * bi-interval for every k-mer.  That bi-interval can then be extended in
 * either direction to obtain a longer match.
 */
struct KmerTable {
    int K = 0;
    std::vector<BiInterval> table;  // indexed by 2K-bit encoding, size 4^K

    /// Encode a mapped-character k-mer (A=0,C=1,G=2,T=3 in encoding space).
    /// Characters must already be in {NUC_A..NUC_T}.
    static uint64_t encode(const uchar* mapped, int k) {
        uint64_t code = 0;
        for (int i = 0; i < k; ++i)
            code = (code << 2) | (mapped[i] - NUC_A);
        return code;
    }

    BiInterval lookup(uint64_t code) const {
        assert(code < table.size());
        return table[code];
    }
};

/**
 * A bidirectional BWT index for DNA sequences based on the move structure and
 * on the FMD Index.  Position-native: all navigation via Position, no flat
 * run_starts/run_chars.  Uses RunPermLF for up()/down() run-walking.
 */
class FMDIndex {
public:
    FMDIndex() = default;

    /// Construct from RLBWT heads (unmapped ASCII) and run lengths.
    FMDIndex(const std::vector<uchar>& bwt_heads,
             const std::vector<ulint>& bwt_run_lengths);

    // --- Accessors ---
    ulint domain() const { return n; }
    ulint num_runs() const { return r; }
    const std::vector<ulint>& get_C() const { return C; }
    uint64_t get_fingerprint() const { return fp; }

    // --- Core operations ---

    /// Count occurrences of each character in the next `size` BWT positions from `pos`.
    std::array<ulint, NUC_SIGMA> count_all_pos(Position pos, ulint size) const;

    /// Count in BWT[lo..hi). Test shim — uses count_all_pos(advance_position(...)).
    std::array<ulint, NUC_SIGMA> count_all(ulint lo, ulint hi) const;

    /// Find the first Position with character `a` in the next `size` positions from `pos`.
    Position find_first_char_pos(Position pos, ulint size, uchar a) const;

    /// Advance a Position by k BWT positions.
    Position advance_position(Position pos, ulint k) const;

    /// Initialize bi-interval for a single character c (mapped).
    BiInterval init_bi_interval(uchar c) const;

    /// Backward extend: given bi-interval of P, return bi-interval of aP.
    BiInterval backward_extend(BiInterval bi, uchar a) const;

    /// Forward extend: given bi-interval of P, return bi-interval of Pa.
    BiInterval forward_extend(BiInterval bi, uchar a) const;

    // --- K-mer tables ---
    void add_kmer_table(KmerTable table);
    const std::vector<KmerTable>& get_kmer_tables() const { return kmer_tables_; }

    /// Cascading k-mer lookup: try tables in descending K order.
    std::pair<BiInterval, int> lookup_kmer(const uchar* mapped_query,
                                           int pos, int len) const;

    // --- Raw BWT access (for testing / serialization) ---
    const std::vector<uchar>& bwt_heads() const { return heads_; }
    const std::vector<ulint>& bwt_run_lengths() const { return run_lengths_; }

    /// Reconstruct the full BWT string (mapped characters). Expensive — test use only.
    std::vector<uchar> reconstruct_bwt() const;

    /// First Position in BWT (for k=0 table construction).
    Position first_position() const { return lf_.first(); }

private:
    mutable LFIndex lf_;  // mutable: library accessors are non-const

    std::vector<ulint> C;  // C[c] = # chars < c in BWT, size σ+1
    std::array<Position, NUC_SIGMA + 1> C_pos_;  // Position at BWT[C[c]] for each c
    std::vector<KmerTable> kmer_tables_;

    ulint n = 0;  // BWT length
    ulint r = 0;  // number of runs

    std::vector<uchar> heads_;
    std::vector<ulint> run_lengths_;
    uint64_t fp = 0;

    /// Get mapped character for a run (Nucleotide mapped space).
    uchar get_mapped_character(ulint interval) const;

    static uint64_t compute_fingerprint(const std::vector<uchar>& heads,
                                        const std::vector<ulint>& run_lengths);
};

/**
 * Compute a fingerprint from the heads and run lengths of the BWT.  This is
 * used to verify that the k-mer tables were built from the same BWT.
 */
inline uint64_t FMDIndex::compute_fingerprint(
        const std::vector<uchar>& heads,
        const std::vector<ulint>& run_lengths) {
    uint64_t h = 14695981039346656037ULL;
    auto feed = [&](const void* data, size_t len) {
        const auto* p = static_cast<const uint8_t*>(data);
        for (size_t i = 0; i < len; ++i) {
            h ^= p[i];
            h *= 1099511628211ULL;
        }
    };
    for (auto c : heads)       feed(&c, sizeof(c));
    for (auto l : run_lengths) feed(&l, sizeof(l));
    return h;
}

/**
 * Get the mapped character for a run (Nucleotide mapped space).
 */
inline uchar FMDIndex::get_mapped_character(ulint interval) const {
    return orbit::nucleotide::map_char(lf_.get_character(interval));
}

/**
 * Construct the FMD index from the BWT heads and run lengths.
 */
inline FMDIndex::FMDIndex(const std::vector<uchar>& bwt_heads,
                          const std::vector<ulint>& bwt_run_lengths)
    : heads_(bwt_heads), run_lengths_(bwt_run_lengths)
{
    assert(bwt_heads.size() == bwt_run_lengths.size());
    r = bwt_heads.size();

    // No associated run data here
    std::vector<std::array<ulint, 0>> empty_run_data(r);
    lf_ = LFIndex(bwt_heads, bwt_run_lengths, orbit::NO_SPLITTING, empty_run_data);
    n = 0;
    for (ulint L : run_lengths_) n += L;
    r = lf_.intervals();

    // Compute C array from input BWT (runperm first/down order can differ
    // from input order)
    std::array<ulint, NUC_SIGMA + 1> char_count{};
    for (ulint i = 0; i < heads_.size(); ++i) {
        uchar c = orbit::nucleotide::map_char(heads_[i]);
        assert(c < NUC_SIGMA);  // any valid Nucleotide (incl. TER, SEP for double-strand)
        char_count[c] += run_lengths_[i];
    }
    C.assign(NUC_SIGMA + 1, 0);
    ulint sum = 0;
    for (uchar c = 0; c <= NUC_SIGMA; ++c) {
        C[c] = sum;
        sum += char_count[c];
    }
    assert(C[NUC_SIGMA] == n);

    // Compute C_pos_: Position at BWT[C[c]] for each c
    ulint nr = lf_.intervals();
    Position p = lf_.first();
    ulint cumul = 0;
    for (ulint i = 0; i < nr; ++i) {
        ulint len = lf_.get_length(p.interval);
        for (uchar c = 0; c <= NUC_SIGMA; ++c) {
            if (C[c] >= cumul && C[c] < cumul + len) {
                C_pos_[c] = Position{p.interval, C[c] - cumul};
            }
        }
        cumul += len;
        if (i + 1 < nr) p = lf_.down(p);
    }

    // Compute fingerprint, useful for checking whether a k-mer table was truly
    // built over a particular BWT index
    fp = compute_fingerprint(bwt_heads, bwt_run_lengths);
}

/**
 * Count the number of occurrences for each alphabet character in the BWT range
 * [lo, hi).
 */
inline std::array<ulint, NUC_SIGMA>
FMDIndex::count_all(ulint lo, ulint hi) const {
    if (lo >= hi) return {};
    return count_all_pos(advance_position(first_position(), lo), hi - lo);
}

/**
 * Count the number of occurrences for each alphabet character in the next
 * 'size' BWT positions starting from pos.
 */
inline std::array<ulint, NUC_SIGMA>
FMDIndex::count_all_pos(Position pos, ulint size) const {
    std::array<ulint, NUC_SIGMA> cnt{};
    if (size == 0) return cnt;

    ulint remaining = size;
    ulint interval = pos.interval;
    ulint offset = pos.offset;

    while (remaining > 0) {
        ulint run_len = lf_.get_length(interval);
        ulint avail = run_len - offset;
        ulint take = std::min(avail, remaining);
        cnt[get_mapped_character(interval)] += take;
        remaining -= take;
        if (remaining > 0) {
            offset = 0;
            if (interval + 1 < r)
                interval++;
            else
                break;
        }
    }
    return cnt;
}

/**
 * Find the first position with character a that is later than pos.
 */
inline Position FMDIndex::find_first_char_pos(Position pos, ulint size, uchar a) const {
    ulint remaining = size;
    ulint interval = pos.interval;
    ulint offset = pos.offset;

    while (remaining > 0) {
        if (get_mapped_character(interval) == a)
            return Position{interval, offset};
        ulint run_len = lf_.get_length(interval);
        ulint avail = run_len - offset;
        ulint take = std::min(avail, remaining);
        remaining -= take;
        if (remaining > 0) {
            offset = 0;
            interval++;
        }
    }
    assert(false && "find_first_char_pos: char not found");
    return pos;
}

/**
 * Advance a Position by k BWT positions.  Can involve advancing along rows.
 */
inline Position FMDIndex::advance_position(Position pos, ulint k) const {
    if (k == 0) return pos;
    ulint run_avail = lf_.get_length(pos.interval) - pos.offset;
    if (k < run_avail) {
        // No more row advances required
        return Position{pos.interval, pos.offset + k};
    }
    k -= run_avail;
    ulint interval = pos.interval + 1;
    // Advance along rows as needed
    while (k > 0 && interval < r) {
        ulint len = lf_.get_length(interval);
        if (k < len) {
            return Position{interval, k};
        }
        k -= len;
        interval++;
    }
    return Position{interval, 0};
}

/**
 * Initialize a bi-interval for a single character c.
 */
inline BiInterval FMDIndex::init_bi_interval(uchar c) const {
    assert(c < NUC_SIGMA);
    BiInterval bi;
    bi.pos = C_pos_[c];
    bi.pos_r = C_pos_[complement(c)];
    bi.size = C[c + 1] - C[c];
    return bi;
}

/**
 * Backward extend a bi-interval by a single character a.
 */
inline BiInterval FMDIndex::backward_extend(BiInterval bi, uchar a) const {
    if (bi.empty()) return {};

    auto cnt = count_all_pos(bi.pos, bi.size);
    ulint s_new = cnt[a];
    if (s_new == 0) return {};

    // Forward strand: LF(first_a) replaces C[a] + rank(a, lo)
    Position first_a = find_first_char_pos(bi.pos, bi.size, a);
    Position pos_new = lf_.LF(first_a);

    // Reverse complement strand: advance by skip
    uchar a_bar = complement(a);
    ulint skip = 0;
    for (uchar b = 0; b < NUC_SIGMA; ++b) {
        if (complement(b) < a_bar)
            skip += cnt[b];
    }
    Position pos_r_new = advance_position(bi.pos_r, skip);

    return {pos_new, pos_r_new, s_new};
}

/**
 * Forward extend a bi-interval by a single character a.
 */
inline BiInterval FMDIndex::forward_extend(BiInterval bi, uchar a) const {
    if (bi.empty()) return {};
    uchar a_bar = complement(a);
    BiInterval bi_rev = {bi.pos_r, bi.pos, bi.size};
    BiInterval result_rev = backward_extend(bi_rev, a_bar);
    if (result_rev.empty()) return {};
    return {result_rev.pos_r, result_rev.pos, result_rev.size};
}

/**
 * Add a k-mer table to the index.
 * TODO: Where do we ask whether fingerprints match?
 */
inline void FMDIndex::add_kmer_table(KmerTable table) {
    auto it = std::lower_bound(kmer_tables_.begin(), kmer_tables_.end(),
                               table.K,
                               [](const KmerTable& t, int k) { return t.K > k; });
    kmer_tables_.insert(it, std::move(table));
}

/**
 * Lookup a k-mer in the k-mer tables.  Returns the bi-interval for the k-mer
 * as well as the k for the table used.
 */
inline std::pair<BiInterval, int>
FMDIndex::lookup_kmer(const uchar* mapped_query, int pos, int len) const {
    for (const auto& tab : kmer_tables_) {
        if (pos + tab.K > len) continue;
        bool valid = true;
        for (int i = 0; i < tab.K; ++i) {
            uchar c = mapped_query[pos + i];
            if (c < NUC_A || c > NUC_T) { valid = false; break; }
        }
        if (!valid) continue;
        if (tab.K == 0) {
            return {tab.table[0], 0};
        }
        uint64_t code = KmerTable::encode(&mapped_query[pos], tab.K);
        BiInterval bi = tab.lookup(code);
        if (!bi.empty()) return {bi, tab.K};
    }
    return {{}, 0};
}

/**
 * Reconstruct a full, non-RLE BWT string (as "mapped" characters).  Useful for
 * testing.
 */
inline std::vector<uchar> FMDIndex::reconstruct_bwt() const {
    std::vector<uchar> bwt(n);
    Position p = lf_.first();
    ulint idx = 0;
    for (ulint i = 0; i < r && idx < n; ++i) {
        uchar c = get_mapped_character(p.interval);
        ulint len = lf_.get_length(p.interval);
        for (ulint j = 0; j < len && idx < n; ++j)
            bwt[idx++] = c;
        if (i + 1 < r)
            p = lf_.down(p);
    }
    return bwt;
}

/**
 * Decode a k-mer from its 2K-bit encoding.
 */
inline std::vector<uchar> decode_kmer(uint64_t code, int K) {
    std::vector<uchar> kmer(K);
    for (int i = K - 1; i >= 0; --i) {
        kmer[i] = static_cast<uchar>((code & 3) + NUC_A);
        code >>= 2;
    }
    return kmer;
}

/**
 * Build a new k-mer table for a given length K over the BWT index in fmd.
 * TODO: could this be faster by coputing common-prefix bi-intervals once,
 * instead redescending for each k-mer?
 */
inline KmerTable build_kmer_table(const FMDIndex& fmd, int K) {
    KmerTable tab;
    tab.K = K;

    if (K == 0) {
        Position first = fmd.first_position();
        tab.table.push_back({first, first, fmd.domain()});
        return tab;
    }

    uint64_t num_entries = 1ULL << (2 * K);
    tab.table.resize(num_entries);

    Position first = fmd.first_position();
    BiInterval full_range = {first, first, fmd.domain()};

    for (uint64_t code = 0; code < num_entries; ++code) {
        auto kmer = decode_kmer(code, K);
        BiInterval bi = full_range;
        for (int i = K - 1; i >= 0 && !bi.empty(); --i)
            bi = fmd.backward_extend(bi, kmer[i]);
        tab.table[code] = bi;
    }

    return tab;
}

#endif /* _BI_FMD_HPP */
