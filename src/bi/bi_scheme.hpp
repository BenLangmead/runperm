/**
 * Approximate matching with search schemes.
 * Substitution-only (Hamming distance).
 *
 * Author: Ben Langmead (ben.langmead@gmail.com)
 * Date: Feb 23, 2026
 */

#ifndef _BI_SCHEME_HPP
#define _BI_SCHEME_HPP

#include "bi_fmd.hpp"
#include <vector>
#include <string>
#include <set>
#include <algorithm>

namespace bi_approx {

/** A single approximate match: [start,end) in pattern, Hamming errors. */
struct ApproxMatch {
    int start = 0;
    int end = 0;
    int errors = 0;
    BiInterval bi;
    bool operator<(const ApproxMatch& o) const {
        if (start != o.start) return start < o.start;
        if (end != o.end) return end < o.end;
        return errors < o.errors;
    }
};

/** One search: permutation pi, lower bounds L, upper bounds U. */
struct Search {
    std::vector<int> pi, L, U;
};

/** Compute piece boundaries for m-length pattern into p pieces. */
inline std::vector<int> compute_boundaries(int m, int p) {
    std::vector<int> b;
    b.push_back(0);
    int base = m / p, rem = m % p;
    for (int i = 0; i < p; ++i)
        b.push_back(b.back() + base + (i < rem ? 1 : 0));
    return b;
}

/** Extend bi with chars, allowing up to max_errors substitutions. */
inline void extend_with_errors(const FMDIndex& fmd, BiInterval bi,
                               const uchar* chars, int len, int max_errors,
                               bool forward, std::vector<std::pair<BiInterval, int>>& out,
                               int errors_so_far = 0) {
    if (len == 0) {
        if (max_errors >= 0) out.push_back({bi, errors_so_far});
        return;
    }
    uchar c = forward ? chars[0] : chars[len - 1];
    const uchar* rest = forward ? chars + 1 : chars;
    int rest_len = len - 1;
    for (uchar a = NUC_A; a <= NUC_T; ++a) {
        BiInterval bi_ext = forward ? fmd.forward_extend(bi, a) : fmd.backward_extend(bi, a);
        if (bi_ext.empty()) continue;
        int cost = (a == c) ? 0 : 1;
        if (cost <= max_errors)
            extend_with_errors(fmd, bi_ext, rest, rest_len, max_errors - cost, forward, out, errors_so_far + cost);
    }
}

/** Exact match a piece [lo, hi), return bi-interval. */
inline BiInterval exact_match_piece(const FMDIndex& fmd, const uchar* mapped_P, int m,
                                    const std::vector<int>& boundaries, int piece) {
    int lo = boundaries[piece], hi = boundaries[piece + 1];
    if (lo >= hi) return {};
    uchar c = mapped_P[lo];
    if (c < NUC_A || c > NUC_T) return {};
    BiInterval bi = fmd.init_bi_interval(c);
    for (int i = lo + 1; i < hi && !bi.empty(); ++i) {
        c = mapped_P[i];
        if (c < NUC_A || c > NUC_T) return {};
        bi = fmd.forward_extend(bi, c);
    }
    return bi;
}

/** Recursive helper for execute_search. */
inline void execute_search_rec(const FMDIndex& fmd, const uchar* mapped_P, int m,
                               const std::vector<int>& boundaries, const Search& search,
                               BiInterval bi, int step, int errors_so_far,
                               int match_lo, int match_hi, std::vector<ApproxMatch>& results);

/** Execute one search, collect matches. */
inline void execute_search(const FMDIndex& fmd, const uchar* mapped_P, int m,
                           const std::vector<int>& boundaries, const Search& search,
                           std::vector<ApproxMatch>& results) {
    int p = static_cast<int>(search.pi.size());
    if (p == 0) return;
    int first = search.pi[0];
    BiInterval bi = exact_match_piece(fmd, mapped_P, m, boundaries, first);
    if (bi.empty()) return;
    int lo = boundaries[first], hi = boundaries[first + 1];
    execute_search_rec(fmd, mapped_P, m, boundaries, search, bi, 1, 0, lo, hi, results);
}

/** Recursive helper for execute_search. */
inline void execute_search_rec(const FMDIndex& fmd, const uchar* mapped_P, int m,
                               const std::vector<int>& boundaries, const Search& search,
                               BiInterval bi, int step, int errors_so_far,
                               int match_lo, int match_hi, std::vector<ApproxMatch>& results) {
    int p = static_cast<int>(search.pi.size());
    if (step >= p) {
        if (errors_so_far >= search.L[p - 1])
            results.push_back({match_lo, match_hi, errors_so_far, bi});
        return;
    }
    int piece = search.pi[step];
    int piece_lo = boundaries[piece], piece_hi = boundaries[piece + 1];
    int piece_len = piece_hi - piece_lo;
    bool forward = (piece > search.pi[step - 1]);
    int max_new = search.U[step] - errors_so_far;
    int min_new = std::max(0, static_cast<int>(search.L[step]) - errors_so_far);
    int new_lo = forward ? match_lo : match_lo - piece_len;
    int new_hi = forward ? match_hi + piece_len : match_hi;
    std::vector<std::pair<BiInterval, int>> ext;
    extend_with_errors(fmd, bi, mapped_P + piece_lo, piece_len, max_new, forward, ext, 0);
    for (const auto& [bi_ext, e] : ext) {
        if (e >= min_new && e <= max_new)
            execute_search_rec(fmd, mapped_P, m, boundaries, search, bi_ext, step + 1,
                              errors_so_far + e, new_lo, new_hi, results);
    }
}

/** Approximate match with arbitrary scheme. */
inline std::vector<ApproxMatch> approx_match_scheme(const FMDIndex& fmd,
                                                   const uchar* mapped_P, int m,
                                                   const std::vector<Search>& scheme) {
    std::vector<ApproxMatch> results;
    if (m <= 0 || scheme.empty()) return results;
    int p = static_cast<int>(scheme[0].pi.size());
    auto boundaries = compute_boundaries(m, p);
    for (const auto& s : scheme)
        execute_search(fmd, mapped_P, m, boundaries, s, results);
    std::set<ApproxMatch> dedup(results.begin(), results.end());
    return std::vector<ApproxMatch>(dedup.begin(), dedup.end());
}

/** Approximate match from string with custom scheme. */
inline std::vector<ApproxMatch> approx_match_scheme(const FMDIndex& fmd,
                                                   const std::string& pattern,
                                                   const std::vector<Search>& scheme) {
    std::vector<uchar> mapped(pattern.size());
    for (size_t i = 0; i < pattern.size(); ++i)
        mapped[i] = orbit::nucleotide::map_char(static_cast<uchar>(pattern[i]));
    return approx_match_scheme(fmd, mapped.data(), static_cast<int>(mapped.size()), scheme);
}

}  // namespace bi_approx

#endif /* _BI_SCHEME_HPP */
