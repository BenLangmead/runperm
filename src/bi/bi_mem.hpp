/**
 * MEM and SMEM finding for the bidirectional index.
 *
 * Author: Ben Langmead (ben.langmead@gmail.com)
 * Date: Feb 23, 2026
 */

#ifndef _BI_MEM_HPP
#define _BI_MEM_HPP

#include "bi_fmd.hpp"
#include <vector>
#include <string>

namespace bi_mem {

/**
 * Result of a MEM/SMEM: [start, end) in query, bi-interval.  The number of
 * occurrences is in bi.size.
 */
struct MemResult {
    int start = 0;   // inclusive
    int end = 0;     // exclusive
    BiInterval bi;
};

/**
 * Find super-maximal exact matches (SMEMs) in query Q.
 * Q must be mapped to Nucleotide space (A,C,G,T only for extension; N yields empty).
 */
inline std::vector<MemResult> find_smems(const FMDIndex& fmd,
                                         const uchar* mapped_Q, int m) {
    std::vector<MemResult> results;
    int i = 0;
    while (i < m) {
        auto [bi, chars] = fmd.lookup_kmer(mapped_Q, i, m);
        int j = i + chars;
        if (bi.empty()) {
            i++;
            continue;
        }
        // Forward extend from j
        while (j < m && !bi.empty()) {
            uchar c = mapped_Q[j];
            if (c < NUC_A || c > NUC_T) break;
            BiInterval bi_next = fmd.forward_extend(bi, c);
            if (bi_next.empty()) break;
            bi = bi_next;
            j++;
        }
        int right_end = j;
        // Backward extend from i-1
        int k = i - 1;
        while (k >= 0 && !bi.empty()) {
            uchar c = mapped_Q[k];
            if (c < NUC_A || c > NUC_T) break;
            BiInterval bi_next = fmd.backward_extend(bi, c);
            if (bi_next.empty()) break;
            bi = bi_next;
            k--;
        }
        int left_end = k + 1;
        if (right_end - left_end > 0)
            results.push_back({left_end, right_end, bi});
        i = std::max(i + 1, right_end);
    }
    return results;
}

/** Map ASCII pattern to Nucleotide and find SMEMs. */
inline std::vector<MemResult> find_smems(const FMDIndex& fmd,
                                         const std::string& pattern) {
    std::vector<uchar> mapped(pattern.size());
    for (size_t i = 0; i < pattern.size(); ++i)
        mapped[i] = orbit::nucleotide::map_char(static_cast<uchar>(pattern[i]));
    return find_smems(fmd, mapped.data(), static_cast<int>(mapped.size()));
}

}  // namespace bi_mem

#endif /* _BI_MEM_HPP */
