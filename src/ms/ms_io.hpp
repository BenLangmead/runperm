/**
 * Build MSIndexSpillLCP from an RLBWT and its per-row LCPs held in memory, or
 * from the files build-rlbwt --lcp-bin reads (see rlbwt_io.hpp).
 *
 * Author: Ben Langmead (ben.langmead@gmail.com)
 * Date: Feb 17, 2026
 */

#ifndef _MS_IO_HPP
#define _MS_IO_HPP

#include "ms_rlbwt.hpp"
#include "ms_build.hpp"
#include "rlbwt_io.hpp"
#include <iostream>
#include <optional>
#include <string>

namespace ms_io {

/** Options for building an MS index. */
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

/** The ms_build options for opts. */
inline ms_build::Options build_options(const BuildOptions& opts) {
    ms_build::Options o;
    o.percentile_k = opts.percentile_k;
    o.coalesce = opts.coalesce;
    o.spill_align = opts.spill_align;
    o.spill_split_bits = opts.spill_split_bits;
    o.minima_only = opts.minima_only;
    o.lf_split = opts.lf_split;
    return o;
}

/**
 * Build MSIndexSpillLCP from run heads, run lengths and each run's full LCP
 * vector (element 0 its top) with given options.  Throws on bad input.
 */
template <bool StoreAbsolutePositions = false>
inline MSIndexSpillLCP<StoreAbsolutePositions>
build_ms_index_spill(std::vector<uchar> bwt_heads, std::vector<ulint> bwt_run_lengths,
                     std::vector<std::vector<ulint>> lcps_per_run, const BuildOptions& opts) {
    apply_lcp_splitting(bwt_heads, bwt_run_lengths, lcps_per_run, opts.split_threshold, opts.minima_only);
    return ms_build::build<StoreAbsolutePositions>(std::move(bwt_heads), std::move(bwt_run_lengths),
                                                   ms_build::lcp_vectors_source(lcps_per_run, opts.minima_only),
                                                   build_options(opts));
}

/**
 * Load an RLBWT and the LCP of each of its rows from PREFIX.bwt.heads,
 * PREFIX.bwt.len and PREFIX.lcp, the files build-rlbwt --lcp-bin reads.  Each
 * run's LCP vector holds its rows' LCPs in order, its top first.  Returns
 * false, after printing why to stderr, if a file is missing or malformed.
 */
inline bool load_rlbwt_lcps(const std::string& prefix, std::vector<uchar>& bwt_heads,
                            std::vector<ulint>& bwt_run_lengths, std::vector<std::vector<ulint>>& lcps_per_run) {
    std::string err;
    if (!rlbwt_io::load_rlbwt(prefix + ".bwt.heads", prefix + ".bwt.len", bwt_heads, bwt_run_lengths, err)) {
        std::cerr << prefix << ": " << err << "\n";
        return false;
    }
    lcps_per_run.assign(bwt_heads.size(), {});
    if (!rlbwt_io::for_each_run_lcps(prefix + ".lcp", bwt_run_lengths,
                                     [&](size_t i, const std::vector<ulint>& v) { lcps_per_run[i] = v; }, err)) {
        std::cerr << prefix << ".lcp: " << err << "\n";
        return false;
    }
    return true;
}

/**
 * Build MSIndexSpillLCP from the files load_rlbwt_lcps reads, with given
 * options.  Uses spillover storage: no truncation, arbitrary run lengths.
 */
template <bool StoreAbsolutePositions = false>
inline std::optional<MSIndexSpillLCP<StoreAbsolutePositions>>
build_ms_index_spill_from_rlbwt(const std::string& prefix, const BuildOptions& opts = BuildOptions{}) {
    std::vector<uchar> bwt_heads;
    std::vector<ulint> bwt_run_lengths;
    std::vector<std::vector<ulint>> lcps_per_run;

    if (!load_rlbwt_lcps(prefix, bwt_heads, bwt_run_lengths, lcps_per_run))
        return std::nullopt;

    try {
        return build_ms_index_spill<StoreAbsolutePositions>(std::move(bwt_heads), std::move(bwt_run_lengths),
                                                            std::move(lcps_per_run), opts);
    } catch (...) {
        return std::nullopt;
    }
}

}  // namespace ms_io

#endif /* _MS_IO_HPP */
