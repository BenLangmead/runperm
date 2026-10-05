/**
 * Read orientation for batch queries (--orient).  The queries match a read
 * right to left, starting at its 3' end, where Illumina base qualities are
 * lowest.  Querying the reverse complement instead starts at the read's
 * high-quality end.
 *
 * Reversing the per-position output of a query on revcomp(R) does not give
 * the matching statistics of R.  Position j of that output is the length of
 * the longest prefix of revcomp(R)[j..] in the text; when the text holds
 * both strands of every sequence, a string occurs exactly when its reverse
 * complement does, so with e = |R| - 1 - j this is E[e], the length of the
 * longest substring of R ending at e that occurs.  ms_from_match_ends turns
 * E into MS.  Whether an index holds both strands cannot be read from it,
 * so the caller has to know.
 */

#ifndef _ORIENT_HPP
#define _ORIENT_HPP

#include <algorithm>
#include <cstddef>
#include <string>
#include <vector>

/** Which strand of each read the batch queries match. */
enum class ReadOrient { FWD, RC, QUAL };

/**
 * Reverse-complements s in place.  A, C, G and T (either case) map to their
 * complements and every other byte, N included, to itself, so a byte that
 * the queries treat as matching nothing still matches nothing.
 */
inline void revcomp_in_place(std::string& s) {
    std::reverse(s.begin(), s.end());
    for (auto& c : s) {
        switch (c) {
            case 'A': c = 'T'; break;
            case 'C': c = 'G'; break;
            case 'G': c = 'C'; break;
            case 'T': c = 'A'; break;
            case 'a': c = 't'; break;
            case 'c': c = 'g'; break;
            case 'g': c = 'c'; break;
            case 't': c = 'a'; break;
            default: break;
        }
    }
}

/**
 * Matching statistics ms[0 .. n) from match-end lengths end_len[0 .. n),
 * where end_len[e] is the length of the longest substring ending at e that
 * occurs in the text.  Substrings of an occurring string occur, so R[i..e]
 * occurs exactly when start(e) = e - end_len[e] + 1 is at most i, and
 * MS[i] = e - i + 1 for the largest such e >= i (0 if there is none).
 * Since end_len[e + 1] <= end_len[e] + 1, start never decreases as e grows,
 * and the e for each i is found with one pointer that only moves forward.
 * ms and end_len must not overlap.
 */
template <typename T>
inline void ms_from_match_ends(const T* end_len, size_t n, T* ms) {
    size_t p = 0;  // number of e with start(e) <= i
    for (size_t i = 0; i < n; ++i) {
        while (p < n && p + 1 <= i + end_len[p]) ++p;  // start(p) <= i
        ms[i] = p > i ? static_cast<T>(p - i) : 0;
    }
}

/**
 * Replaces v, the output of a query on the reverse complement of a read,
 * with the matching statistics of the read itself, using scratch as space.
 */
template <typename T>
inline void ms_from_revcomp_query(std::vector<T>& v, std::vector<T>& scratch) {
    std::reverse(v.begin(), v.end());
    scratch.resize(v.size());
    ms_from_match_ends(v.data(), v.size(), scratch.data());
    v.swap(scratch);
}

/** Number of bases at each end of a read that --orient qual compares. */
inline constexpr size_t ORIENT_QUAL_BASES = 20;

/**
 * Whether --orient qual queries the reverse complement of a read with base
 * qualities qual: when the mean quality of its last min(20, len) bases is
 * higher than that of its first min(20, len).  Both means are over the same
 * number of bases, so their sums are compared, and the quality encoding's
 * offset does not matter.
 */
inline bool qual_prefers_rc(const std::string& qual) {
    const size_t m = std::min(ORIENT_QUAL_BASES, qual.size());
    long head = 0, tail = 0;
    for (size_t i = 0; i < m; ++i) {
        head += static_cast<unsigned char>(qual[i]);
        tail += static_cast<unsigned char>(qual[qual.size() - 1 - i]);
    }
    return tail > head;
}

/**
 * Reverse-complements in place the reads in seqs that orientation o
 * queries on the other strand (none for FWD, all for RC, and for QUAL those
 * that qual_prefers_rc chooses from quals), sets flipped[j] for each, and
 * returns how many there are.
 */
inline size_t orient_reads(std::vector<std::string>& seqs, const std::vector<std::string>* quals, ReadOrient o,
                           std::vector<char>& flipped) {
    flipped.assign(seqs.size(), 0);
    if (o == ReadOrient::FWD) return 0;
    size_t count = 0;
    for (size_t j = 0; j < seqs.size(); ++j) {
        if (o == ReadOrient::QUAL && !qual_prefers_rc((*quals)[j])) continue;
        revcomp_in_place(seqs[j]);
        flipped[j] = 1;
        ++count;
    }
    return count;
}

/**
 * Turns the query output ms[j] of each flipped read into the matching
 * statistics of the read as given.
 */
template <typename T>
inline void restore_orientation(const std::vector<char>& flipped, std::vector<std::vector<T>>& ms,
                                std::vector<T>& scratch) {
    for (size_t j = 0; j < flipped.size(); ++j)
        if (flipped[j]) ms_from_revcomp_query(ms[j], scratch);
}

#endif /* _ORIENT_HPP */
