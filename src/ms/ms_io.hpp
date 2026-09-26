/**
 * Build MSIndexSpillLCP from TSV and related I/O.
 * TSV parsing is in tsv.hpp; this bridges parsed data to the index.
 *
 * Author: Ben Langmead (ben.langmead@gmail.com)
 * Date: Feb 17, 2026
 */

#ifndef _MS_IO_HPP
#define _MS_IO_HPP

#include "ms_rlbwt.hpp"
#include "tsv.hpp"
#include <optional>
#include <string>

namespace ms_io {

/** Options for building an MS index from TSV. */
struct BuildOptions {
    double percentile_k = 0.98;
    ulint split_threshold = SPLIT_THRESHOLD_NEVER;
    bool coalesce = false;
    ulint spill_align = 0;
    uchar spill_split_bits = 0;
    bool minima_only = false;
    // Orbit's LF splitting, applied after LCP splitting (see
    // apply_lf_splitting); it must split.
    orbit::split_params lf_split = orbit::split_params{};
};

/**
 * Build MSIndexSpillLCP from run heads, run lengths and each run's full LCP
 * vector (element 0 its top) with given options.  Throws on bad input.
 */
template <bool StoreAbsolutePositions = false>
inline MSIndexSpillLCP<StoreAbsolutePositions>
build_ms_index_spill(std::vector<uchar> bwt_heads, std::vector<ulint> bwt_run_lengths,
                     std::vector<std::vector<ulint>> lcps_per_run, const BuildOptions& opts) {
    apply_lcp_splitting(bwt_heads, bwt_run_lengths, lcps_per_run, opts.split_threshold, opts.minima_only);
    auto pairs = retained_lcp_pairs(lcps_per_run, opts.minima_only);
    std::vector<std::vector<ulint>>().swap(lcps_per_run);
    apply_lf_splitting(bwt_heads, bwt_run_lengths, pairs, opts.lf_split, opts.minima_only);
    auto [run_data, spill_vectors, max_top, max_sub, skinny_count, jumbo_count] = build_spill_data_from_pairs(
        std::move(pairs), opts.percentile_k, opts.coalesce, false, opts.spill_align, opts.spill_split_bits);
    return MSIndexSpillLCP<StoreAbsolutePositions>(bwt_heads, bwt_run_lengths, run_data, std::move(spill_vectors),
                                                   max_top, max_sub, opts.spill_align, opts.spill_split_bits);
}

/**
 * Build MSIndexSpillLCP from a movify TSV path with given options.
 * Uses spillover storage: no truncation, arbitrary run lengths.
 */
template <bool StoreAbsolutePositions = false>
inline std::optional<MSIndexSpillLCP<StoreAbsolutePositions>>
build_ms_index_spill_from_tsv(const std::string& path, const BuildOptions& opts) {
    std::vector<uchar> bwt_heads;
    std::vector<ulint> bwt_run_lengths;
    std::vector<std::vector<ulint>> lcps_per_run;

    if (!tsv::load_tsv(path, bwt_heads, bwt_run_lengths, lcps_per_run))
        return std::nullopt;

    try {
        return build_ms_index_spill<StoreAbsolutePositions>(std::move(bwt_heads), std::move(bwt_run_lengths),
                                                            std::move(lcps_per_run), opts);
    } catch (...) {
        return std::nullopt;
    }
}

/**
 * Build MSIndexSpillLCP from a movify TSV path (default options).
 */
template <bool StoreAbsolutePositions = false>
inline std::optional<MSIndexSpillLCP<StoreAbsolutePositions>>
build_ms_index_spill_from_tsv(const std::string& path) {
    return build_ms_index_spill_from_tsv<StoreAbsolutePositions>(path, BuildOptions{});
}

}  // namespace ms_io

#endif /* _MS_IO_HPP */
