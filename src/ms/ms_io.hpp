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
};

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
        apply_lcp_splitting(bwt_heads, bwt_run_lengths, lcps_per_run, opts.split_threshold);
        auto [run_data, spill_vectors, max_top, max_sub, skinny_count, jumbo_count] = build_spill_data(
            lcps_per_run, opts.percentile_k, opts.coalesce, false, opts.split_threshold,
            opts.spill_align, opts.spill_split_bits);
        return MSIndexSpillLCP<StoreAbsolutePositions>(bwt_heads, bwt_run_lengths, run_data,
            std::move(spill_vectors), max_top, max_sub, opts.spill_align, opts.spill_split_bits);
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
