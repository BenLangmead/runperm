/**
 * Move structure with compressed LCPs support (spillover: skinny/jumbo rows).
 * Provides MSIndexSpillLCP and build_spill_data for matching statistics.
 *
 * Data that cannot fit in the main tabular move structure is instead stored in
 * a spillover vector.  The spillover vector is necessary when the row is
 * "jumbo", meaning we lacked enough space to store the top and min LCPs in the
 * main table.  Also, the spillover vector is necessary for any row that has
 * at least one uncompressed interior LCP value.
 *
 * Has facilities for splitting rows so that minimal or near-minimal LCPs can
 * be used as run boundaries to the maximum degree possible.  This in turn
 * allows for more compression of interior LCPs.  The rows are then cut into
 * Orbit's split and balanced LF intervals (apply_lf_splitting), which keeps
 * each LF step's fast-forward short and the length and offset columns narrow.
 *
 * Author: Ben Langmead (ben.langmead@gmail.com)
 * Date: Feb 17, 2026
 */

#ifndef _MS_RLBWT_HPP
#define _MS_RLBWT_HPP

#include "orbit/rlbwt.hpp"
#include "orbit/common.hpp"
#include "ms_constants.hpp"
#include <cassert>
#include <iostream>
#include <utility>
#include <optional>
#include <vector>
#include <string>
#include <algorithm>
#include <limits>
#include <map>
#include <numeric>
#include <array>
#include <tuple>
#include <random>
#include <stdexcept>
#include <type_traits>

using uchar = orbit::uchar;
using ulint = orbit::ulint;

#ifdef MS_STATS
/**
 * Counts from ms_query_batch, for builds with MS_STATS: LF steps and the rows
 * their fast-forward skips, and repositions and the rows their walks up and
 * down read.
 */
struct MsStats {
    ulint lf_steps = 0, lf_ff = 0, repositions = 0, walk_up = 0, walk_down = 0;
};
inline MsStats ms_stats;
#define MS_COUNT(field, v) (ms_stats.field += (v))
#else
#define MS_COUNT(field, v) ((void)0)
#endif

/**
 * Decide which interior LCP values of a run must be stored.  lcp[0] is the
 * run's top (boundary) LCP and lcp_next, if present, is the next run's LCP
 * vector, whose first element is the boundary below this run.
 *
 * Queries only ask for the minimum over rows [0, o] together with the top
 * boundary, or over rows (o, end) followed by the bottom boundary.  With
 * minima_only false, an interior value is kept when it is below the larger of
 * the two boundaries.  With minima_only true, it is kept only when it is a
 * strict new minimum scanning down from the top boundary or scanning up from
 * the bottom boundary; those are the only values that can answer a query.
 * Element 0 of the result is unused.
 */
inline std::vector<bool> lcp_keep_mask(const std::vector<ulint>& lcp,
                                       const std::vector<ulint>* lcp_next,
                                       bool minima_only = false) {
    std::vector<bool> keep(lcp.size(), false);
    if (lcp.size() < 2) return keep;
    const bool has_next = lcp_next && !lcp_next->empty();
    if (!minima_only) {
        ulint m = lcp[0];
        if (has_next) m = std::max(m, (*lcp_next)[0]);
        for (size_t i = 1; i < lcp.size(); ++i) keep[i] = lcp[i] < m;
        return keep;
    }
    ulint m = lcp[0];
    for (size_t i = 1; i < lcp.size(); ++i)
        if (lcp[i] < m) { keep[i] = true; m = lcp[i]; }
    // With no run below, there is no bottom boundary to cap a downward
    // query, so every suffix minimum is needed.
    m = has_next ? (*lcp_next)[0] : std::numeric_limits<ulint>::max();
    for (size_t i = lcp.size() - 1; i >= 1; --i)
        if (lcp[i] < m) { keep[i] = true; m = lcp[i]; }
    return keep;
}

/**
 * compress_lcps: replace interior LCP values that are not needed (see
 *                lcp_keep_mask) with placeholders that will not be stored.
 * Appends to out: [top_lcp, compressed interior...] where gaps are LCP_GAP.
 */
inline void compress_lcps(const std::vector<ulint>& lcp,
                          const std::vector<ulint>* lcp_next,
                          std::vector<ulint>& out,
                          bool minima_only = false) {
    out.clear();
    if (lcp.empty()) return;
    out.push_back(lcp[0]);
    auto keep = lcp_keep_mask(lcp, lcp_next, minima_only);
    for (size_t i = 1; i < lcp.size(); ++i)
        out.push_back(keep[i] ? lcp[i] : LCP_GAP);
}

constexpr ulint NO_SPILL = 0;

enum class LCPSpillRunCols {
    LCP_TOP,
    LCP_MIN_SUB,
    LCP_SPILL,
    COUNT
};

using SpilloverVector = std::vector<uchar>;

/**
 * Append a ULEB128-encoded integer to a vector.  If the remainder of the
 * integer requires more than 7 bits, then a continuation bit is set and we
 * iterate after shifting by 7.
 */
inline void append_uleb128(std::vector<uchar>& out, ulint v) {
    do {
        uchar b = static_cast<uchar>(v & 0x7F);
        v >>= 7;
        if (v) b |= 0x80;
        out.push_back(b);
    } while (v);
}

/**
 * Decode a ULEB128-encoded integer from a vector.  The integer is returned
 * and the position is incremented to the byte just after the parsed number.
 */
inline std::pair<ulint, size_t> decode_uleb128(const std::vector<uchar>& buf, size_t pos) {
    ulint v = 0;
    int shift = 0;
    for (;;) {
        uchar b = buf[pos++];
        v |= static_cast<ulint>(b & 0x7F) << shift;
        if ((b & 0x80) == 0) break;
        shift += 7;
    }
    return {v, pos};
}

/**
 * Skip a ULEB128-encoded integer from a vector.  The position is incremented
 * to the byte just after the parsed number.
 */
inline size_t skip_uleb128(const std::vector<uchar>& buf, size_t pos) {
    return decode_uleb128(buf, pos).second;
}

namespace detail {

/**
 * Convert an LCP vector with some ignored values into a compressed
 * representation with retained LCPs encoded along with their offsets.
 */
inline std::vector<std::pair<ulint, ulint>>
compressed_to_pairs(const std::vector<ulint>& lcp) {
    std::vector<std::pair<ulint, ulint>> pairs;
    for (size_t i = 0; i < lcp.size(); ++i)
        if (lcp[i] != LCP_GAP)
            pairs.emplace_back(static_cast<ulint>(i), lcp[i]);
    return pairs;
}

/**
 * Compute the k-th percentile of a vector of integers.
 */
inline ulint percentile(std::vector<ulint> vals, double k) {
    if (vals.empty()) return 0;
    size_t idx = static_cast<size_t>(k * vals.size());
    if (idx >= vals.size()) idx = vals.size() - 1;
    std::nth_element(vals.begin(), vals.begin() + static_cast<ptrdiff_t>(idx), vals.end());
    return vals[idx];
}

/**
 * Compute the number of bits needed to represent an integer.
 */
inline uchar ceil_log2(ulint v) {
    if (v == 0) return 1;
    ulint b = 0, x = v;
    while (x > 0) { ++b; x >>= 1; }
    return static_cast<uchar>(b);
}

}  // namespace detail

// Constant for when we don't want to split at all
constexpr ulint SPLIT_THRESHOLD_NEVER = std::numeric_limits<ulint>::max();

/**
 * Compute the compression "improvement" of an LCP vector, which is equal to
 * the number of LCPs that are ignored.
 */
inline size_t compression_improvement(const std::vector<ulint>& lcp,
                                      const std::vector<ulint>* lcp_next,
                                      bool minima_only = false)
{
    if (lcp.size() < 2) return 0;
    auto keep = lcp_keep_mask(lcp, lcp_next, minima_only);
    size_t c = 0;
    for (size_t i = 1; i < lcp.size(); ++i)
        if (!keep[i]) ++c;
    return c;
}

/**
 * Attempt to split an LCP vector at its minimum (if not top or bottom) into
 * two parts for better compression. Returns (did_split, top_part, bottom_part);
 * if no split, returns (False, lcp, empty).
 */
inline std::tuple<bool, std::vector<ulint>, std::vector<ulint>>
possibly_split_lcps(const std::vector<ulint>& lcp,
                    const std::vector<ulint>* lcp_next,
                    ulint split_threshold = SPLIT_THRESHOLD_NEVER,
                    bool minima_only = false)
{
    std::vector<ulint> empty;
    if (lcp.size() < 3) return {false, lcp, empty};
    ulint min_val = lcp[1];
    for (size_t i = 2; i + 1 < lcp.size(); ++i)
        min_val = std::min(min_val, lcp[i]);
    std::vector<size_t> candidates;
    for (size_t i = 1; i + 1 < lcp.size(); ++i)
        if (lcp[i] == min_val) candidates.push_back(i);
    size_t mid = lcp.size() / 2;
    size_t min_idx = candidates[0];
    size_t best = (candidates[0] >= mid) ? (candidates[0] - mid) : (mid - candidates[0]);
    for (size_t c : candidates) {
        size_t d = (c >= mid) ? (c - mid) : (mid - c);
        if (d < best) { best = d; min_idx = c; }
    }
    std::vector<ulint> top(lcp.begin(), lcp.begin() + static_cast<ptrdiff_t>(min_idx));
    std::vector<ulint> bot;
    bot.reserve(lcp.size() - min_idx);
    bot.push_back(min_val);
    for (size_t i = min_idx + 1; i < lcp.size(); ++i) bot.push_back(lcp[i]);
    size_t base_imp = compression_improvement(lcp, lcp_next, minima_only);
    size_t top_imp = compression_improvement(top, &bot, minima_only);
    size_t bot_imp = compression_improvement(bot, lcp_next, minima_only);
    size_t gain = (top_imp + bot_imp > base_imp) ? (top_imp + bot_imp - base_imp) : 0;
    if (gain > split_threshold) return {true, std::move(top), std::move(bot)};
    return {false, lcp, empty};
}

namespace detail {

/**
 * Recursively split an LCP vector.  The recursion helps to check if the
 * partitions from one split should themselves be split.
 */
inline void recursive_split_lcps(const std::vector<ulint>& lcp,
                                 const std::vector<ulint>* lcp_next,
                                 size_t orig_i, ulint split_threshold,
                                 std::vector<std::pair<std::vector<ulint>, size_t>>& out,
                                 bool minima_only = false)
{
    if (lcp.size() < 3) { out.emplace_back(lcp, orig_i); return; }
    auto [did, top, bot] = possibly_split_lcps(lcp, lcp_next, split_threshold, minima_only);
    if (!did || bot.empty()) { out.emplace_back(lcp, orig_i); return; }
    recursive_split_lcps(top, &bot, orig_i, split_threshold, out, minima_only);
    recursive_split_lcps(bot, lcp_next, orig_i, split_threshold, out, minima_only);
}

}  // namespace detail

/**
 * Apply LCP splitting to the LCP vectors for each run.
 */
inline void apply_lcp_splitting(std::vector<uchar>& bwt_heads,
                               std::vector<ulint>& bwt_run_lengths,
                               std::vector<std::vector<ulint>>& lcps_per_run,
                               ulint split_threshold,
                               bool minima_only = false) {
    if (split_threshold == SPLIT_THRESHOLD_NEVER || lcps_per_run.empty()) return;
    const size_t n = lcps_per_run.size();
    std::vector<std::pair<std::vector<ulint>, size_t>> new_lcps;
    for (size_t i = 0; i < n; ++i) {
        const std::vector<ulint>* next = (i + 1 < n) ? &lcps_per_run[i + 1] : nullptr;
        detail::recursive_split_lcps(lcps_per_run[i], next, i, split_threshold, new_lcps, minima_only);
    }
    std::vector<uchar> nh;
    std::vector<ulint> nl;
    std::vector<std::vector<ulint>> nlcps;
    for (const auto& [lcp, oi] : new_lcps) {
        nh.push_back(bwt_heads[oi]);
        nl.push_back(static_cast<ulint>(lcp.size()));
        nlcps.push_back(lcp);
    }
    bwt_heads = std::move(nh);
    bwt_run_lengths = std::move(nl);
    lcps_per_run = std::move(nlcps);
}

/**
 * Compute the maximum value for an unsigned integer using the given number of
 * bits.
 */
inline ulint max_for_bits(uchar w) {
    return (w >= 64) ? ~0ULL : ((1ULL << w) - 1);
}

/**
 * Per-run retained LCPs as (offset, value) pairs: element 0 is (0, top LCP)
 * and the rest are the stored interior values in increasing offset order.
 */
using RunLcpPairs = std::vector<std::pair<ulint, ulint>>;

/**
 * Reduce full per-row LCP vectors to the retained pairs of each run (see
 * lcp_keep_mask for which interior values are kept).
 */
inline std::vector<RunLcpPairs>
retained_lcp_pairs(const std::vector<std::vector<ulint>>& lcps_per_run, bool minima_only = false) {
    const size_t r = lcps_per_run.size();
    std::vector<RunLcpPairs> all_pairs(r);
    std::vector<ulint> full;
    for (size_t i = 0; i < r; ++i) {
        const std::vector<ulint>* next = (i + 1 < r) ? &lcps_per_run[i + 1] : nullptr;
        full.clear();
        compress_lcps(lcps_per_run[i], next, full, minima_only);
        all_pairs[i] = detail::compressed_to_pairs(full);
    }
    return all_pairs;
}

/**
 * Keep only the pairs of one row's retained LCPs (element 0 its top) that a
 * query can need, by the rule of lcp_keep_mask, where bottom is the top LCP
 * of the row below, if there is one.  The offsets without a pair count as
 * holding values above every stored one.
 */
inline void prune_lcp_pairs(RunLcpPairs& p, std::optional<ulint> bottom, bool minima_only) {
    if (p.size() < 2) return;
    std::vector<bool> keep(p.size(), false);
    keep[0] = true;
    if (!minima_only) {
        const ulint m = bottom ? std::max(p[0].second, *bottom) : p[0].second;
        for (size_t j = 1; j < p.size(); ++j) keep[j] = p[j].second < m;
    } else {
        ulint m = p[0].second;
        for (size_t j = 1; j < p.size(); ++j)
            if (p[j].second < m) { keep[j] = true; m = p[j].second; }
        m = bottom ? *bottom : std::numeric_limits<ulint>::max();
        for (size_t j = p.size() - 1; j >= 1; --j)
            if (p[j].second < m) { keep[j] = true; m = p[j].second; }
    }
    size_t w = 0;
    for (size_t j = 0; j < p.size(); ++j)
        if (keep[j]) p[w++] = p[j];
    p.resize(w);
}

/**
 * The lengths of the rows into which Orbit's LF splitting sp cuts an RLBWT
 * with the given run heads and lengths, in row order.  Splitting only cuts
 * runs, so each run's rows are consecutive and their lengths sum to its.
 */
inline std::vector<ulint> lf_split_row_lengths(const std::vector<uchar>& heads, const std::vector<ulint>& lens,
                                               const orbit::split_params& sp) {
    using Enc = orbit::rlbwt::rlbwt_interval_encoding<>;
    const Enc enc = Enc::lf_interval_encoding(heads, lens, sp);
    std::vector<ulint> rows(static_cast<size_t>(enc.intervals()));
    for (ulint k = 0; k < enc.intervals(); ++k) rows[k] = enc.get_length(k);
    return rows;
}

/**
 * Throw unless sp splits LF, with length capping or balancing: an ms index's
 * LF is always split.
 */
inline void check_lf_split(const orbit::split_params& sp) {
    const bool caps = sp.length_capping.has_value() && *sp.length_capping > 0;
    const bool balances = sp.balancing.has_value() && *sp.balancing > 0;
    if (!caps && !balances)
        throw std::invalid_argument("LF must be split, with length capping or balancing; unsplit LF is not supported");
}

/**
 * Cut the runs of an RLBWT into the rows of Orbit's LF splitting sp, so that
 * the move structure, built on those rows as they are, gets split and
 * balanced LF intervals, each with LCP data of its own.  sp must split (see
 * check_lf_split).  heads and lens become the rows'
 * and pairs, each run's retained LCP pairs (as from retained_lcp_pairs or a
 * minima file), become each row's.  Rows may be runs already cut by
 * apply_lcp_splitting; they are cut further.
 *
 * A row cut from inside a run at offset a needs a top LCP, which the
 * retained pairs need not hold.  It gets max(u, d), where u is the minimum of
 * the run's values at offsets up to a, top included, and d the minimum of
 * those from a on together with the next run's top.  Both are exact from
 * the retained pairs, and u, d <= LCP[a].  This leaves every result the
 * same.  A reposition walks through rows whose character is not the
 * target's, so a walk up that crosses offset a of a run goes on to the run's
 * top, and its minimum includes that of [0, a], at most u.  Likewise a walk
 * down that crosses a crosses the rest of the run and the top of the next,
 * with a minimum of at most d.  Taking the minimum with a value no smaller
 * changes nothing.
 * Each row's pairs are then pruned with prune_lcp_pairs.
 */
inline void apply_lf_splitting(std::vector<uchar>& heads, std::vector<ulint>& lens, std::vector<RunLcpPairs>& pairs,
                               const orbit::split_params& sp, bool minima_only = false) {
    check_lf_split(sp);
    if (pairs.size() != heads.size() || lens.size() != heads.size())
        throw std::invalid_argument("heads, lens and LCP pairs must have one entry per run");
    const std::vector<ulint> rows = lf_split_row_lengths(heads, lens, sp);
    std::vector<uchar> nh;
    std::vector<ulint> nl;
    std::vector<RunLcpPairs> np;
    nh.reserve(rows.size());
    nl.reserve(rows.size());
    np.reserve(rows.size());
    std::vector<ulint> starts, tops, suffix_min;
    size_t k = 0;
    for (size_t i = 0; i < heads.size(); ++i) {
        const std::optional<ulint> next_top =
            i + 1 < pairs.size() ? std::optional<ulint>(pairs[i + 1][0].second) : std::nullopt;
        // The offsets in run i at which its rows start.
        starts.clear();
        for (ulint covered = 0; covered < lens[i]; covered += rows[k++]) {
            if (k >= rows.size()) throw std::logic_error("LF split rows do not cover the runs");
            starts.push_back(covered);
        }
        if (std::accumulate(rows.begin() + static_cast<ptrdiff_t>(k - starts.size()),
                            rows.begin() + static_cast<ptrdiff_t>(k), ulint{0}) != lens[i])
            throw std::logic_error("an LF split row crosses a run boundary");
        RunLcpPairs& p = pairs[i];
        if (starts.size() == 1) {
            nh.push_back(heads[i]);
            nl.push_back(lens[i]);
            np.push_back(std::move(p));
            continue;
        }
        // suffix_min[j]: minimum of the values of pairs j onward.
        suffix_min.assign(p.size() + 1, LCP_GAP);
        for (size_t j = p.size(); j-- > 0;) suffix_min[j] = std::min(suffix_min[j + 1], p[j].second);
        tops.assign(starts.size(), 0);
        tops[0] = p[0].second;
        ulint prefix_min = p[0].second;
        size_t j = 1;  // the first pair not yet in prefix_min
        for (size_t s = 1; s < starts.size(); ++s) {
            size_t first_at = j;  // the first pair at or after starts[s]
            while (first_at < p.size() && p[first_at].first < starts[s]) ++first_at;
            for (; j < p.size() && p[j].first <= starts[s]; ++j) prefix_min = std::min(prefix_min, p[j].second);
            ulint d = suffix_min[first_at];
            if (next_top) d = std::min(d, *next_top);
            tops[s] = d == LCP_GAP ? prefix_min : std::max(prefix_min, d);
        }
        size_t at = 1;  // the next interior pair to place
        for (size_t s = 0; s < starts.size(); ++s) {
            const ulint lo = starts[s], hi = s + 1 < starts.size() ? starts[s + 1] : lens[i];
            RunLcpPairs row{{0, tops[s]}};
            for (; at < p.size() && p[at].first < hi; ++at)
                if (p[at].first > lo) row.emplace_back(p[at].first - lo, p[at].second);
            prune_lcp_pairs(row, s + 1 < starts.size() ? std::optional<ulint>(tops[s + 1]) : next_top, minima_only);
            nh.push_back(heads[i]);
            nl.push_back(hi - lo);
            np.push_back(std::move(row));
        }
    }
    if (k != rows.size()) throw std::logic_error("LF split rows do not match the runs");
    heads = std::move(nh);
    lens = std::move(nl);
    pairs = std::move(np);
}

/**
 * Build the full spillover data for the index from each run's retained LCP
 * pairs (as produced by retained_lcp_pairs, or read from a TeraLCP minima
 * file).
 * spill_align: 0=none, else align chunks to multiples of spill_align bytes.
 * spill_split_bits: X LSBs of row index select which spillover array (0=single array).
 * spill_off is stored as byte_offset / spill_align to save bits.
 */
inline std::tuple<
    std::vector<std::array<ulint, static_cast<size_t>(LCPSpillRunCols::COUNT)>>,
    std::vector<SpilloverVector>, ulint, ulint, size_t, size_t>
build_spill_data_from_pairs(std::vector<RunLcpPairs> all_pairs,
                            double percentile_k = 0.98,
                            bool coalesce_spillover = false,
                            bool coalesce_lcp_separately = false,
                            ulint spill_align = 0,
                            uchar spill_split_bits = 0)
{
    if (coalesce_lcp_separately && spill_align > 0)
        throw std::invalid_argument("coalescing LCPs separately is not compatible with spill-align");
    const size_t r = all_pairs.size();
    size_t skinny_count = 0, jumbo_count = 0;
    std::vector<ulint> all_top, all_sub;

    for (size_t i = 0; i < r; ++i) {
        const auto& pairs = all_pairs[i];
        assert(!pairs.empty());
        if (pairs.size() == 1) {
            all_top.push_back(pairs[0].second);
            all_sub.push_back(0);
        } else {
            ulint min_s = LCP_GAP;
            for (size_t j = 1; j < pairs.size(); ++j)
                min_s = std::min(min_s, pairs[j].second);
            ulint sub = (min_s != LCP_GAP && min_s < pairs[0].second)
                ? (pairs[0].second - min_s) : 0;
            all_top.push_back(pairs[0].second);
            all_sub.push_back(sub);
        }
    }

    ulint p_top = detail::percentile(all_top, percentile_k);
    ulint p_sub = detail::percentile(all_sub, percentile_k);
    uchar w_top = detail::ceil_log2(p_top);
    uchar w_sub = detail::ceil_log2(p_sub);
    ulint max_top = max_for_bits(w_top);
    ulint max_sub = max_for_bits(w_sub);

    auto is_skinny = [max_top, max_sub](ulint t, ulint s) {
        if (t > max_top || s > max_sub) return false;
        return (t < max_top || s < max_sub);
    };

    const size_t align_x = (spill_align > 0) ? static_cast<size_t>(spill_align) : 1;
    const size_t num_arrays = (spill_split_bits > 0) ? (size_t{1} << spill_split_bits) : 1;

    std::vector<std::array<ulint, static_cast<size_t>(LCPSpillRunCols::COUNT)>> run_data(r);
    std::vector<SpilloverVector> spill_vectors(num_arrays);
    for (size_t b = 0; b < num_arrays; ++b) {
        spill_vectors[b].reserve(4096);
        for (size_t i = 0; i < align_x; ++i) spill_vectors[b].push_back(0);
    }

    // Per-bucket coalesce maps (when coalesce_spillover).
    std::vector<std::map<SpilloverVector, size_t>> jumbo_maps(num_arrays);
    std::vector<std::map<SpilloverVector, size_t>> lcp_maps(num_arrays);

    auto append_payload = [&](size_t bucket, SpilloverVector& p) -> size_t {
        auto& spill = spill_vectors[bucket];
        while (spill.size() % align_x != 0) spill.push_back(0);
        size_t off = spill.size();
        spill.insert(spill.end(), p.begin(), p.end());
        return off;
    };

    for (size_t i = 0; i < r; ++i) {
        const size_t bucket = (spill_split_bits > 0) ? (i & ((1ULL << spill_split_bits) - 1)) : 0;
        const auto& pairs = all_pairs[i];
        ulint top_val = all_top[i], sub_val = all_sub[i];

        if (pairs.empty()) {
            run_data[i][0] = run_data[i][1] = 0;
            run_data[i][2] = NO_SPILL;
            ++skinny_count;
        } else if (is_skinny(top_val, sub_val)) {
            run_data[i][0] = top_val;
            run_data[i][1] = sub_val;
            if (pairs.size() == 1) {
                run_data[i][2] = NO_SPILL;
            } else {
                SpilloverVector p; // temporary holder of spillover data
                append_uleb128(p, static_cast<ulint>(pairs.size() - 1));
                for (size_t j = 1; j < pairs.size(); ++j) {
                    assert(pairs[j].first > 0 && pairs[j].second != LCP_GAP);
                    append_uleb128(p, pairs[j].first);
                    append_uleb128(p, pairs[j].second);
                }
                if (coalesce_spillover) {
                    auto* lcp_p = coalesce_lcp_separately ? &lcp_maps[bucket] : &jumbo_maps[bucket];
                    auto it = lcp_p->find(p);
                    if (it != lcp_p->end())
                        run_data[i][2] = static_cast<ulint>(it->second);
                    else {
                        size_t off = append_payload(bucket, p);
                        ulint stored = static_cast<ulint>(off / align_x);
                        run_data[i][2] = stored;
                        (*lcp_p)[std::move(p)] = stored;
                    }
                } else {
                    size_t off = append_payload(bucket, p);
                    run_data[i][2] = static_cast<ulint>(off / align_x);
                }
            }
            ++skinny_count;
        } else {
            run_data[i][0] = max_top;
            run_data[i][1] = max_sub;
            ulint row_min = LCP_GAP;
            for (const auto& [o, v] : pairs) row_min = std::min(row_min, v);
            SpilloverVector p_head, p_lcp, p_full; // temporary holders
            append_uleb128(p_full, top_val);
            append_uleb128(p_full, row_min);
            append_uleb128(p_full, static_cast<ulint>(pairs.size() - 1));
            if(coalesce_spillover) {
                append_uleb128(p_head, top_val);
                append_uleb128(p_head, row_min);
                append_uleb128(p_lcp, static_cast<ulint>(pairs.size() - 1));
            }
            for (size_t j = 1; j < pairs.size(); ++j) {
                assert(pairs[j].first > 0 && pairs[j].second != LCP_GAP);
                append_uleb128(p_full, pairs[j].first);
                append_uleb128(p_full, pairs[j].second);
                if(coalesce_spillover) {
                    append_uleb128(p_lcp, pairs[j].first);
                    append_uleb128(p_lcp, pairs[j].second);
                }
            }
            if (coalesce_spillover) {
                auto* jumbo_p = &jumbo_maps[bucket];
                auto* lcp_p = coalesce_lcp_separately ? &lcp_maps[bucket] : jumbo_p;
                auto it = jumbo_p->find(p_full);
                if (it != jumbo_p->end())
                    run_data[i][2] = static_cast<ulint>(it->second);
                else {
                    // Readers decode the head and then the LCP record at the
                    // next byte, so the two must be contiguous with no
                    // alignment padding between them.
                    size_t off_head = append_payload(bucket, p_full);
                    size_t off_lcp = off_head + p_head.size();
                    ulint stored = static_cast<ulint>(off_head / align_x);
                    run_data[i][2] = stored;
                    (*jumbo_p)[std::move(p_full)] = stored;
                    if (coalesce_lcp_separately)
                        (*lcp_p)[std::move(p_lcp)] = static_cast<size_t>(off_lcp / align_x);
                }
            } else {
                size_t off = append_payload(bucket, p_full);
                run_data[i][2] = static_cast<ulint>(off / align_x);
            }
            ++jumbo_count;
        }
    }
    return {run_data, spill_vectors, max_top, max_sub, skinny_count, jumbo_count};
}

/**
 * Build the full spillover data for the index from full per-row LCP vectors.
 * split_threshold is accepted for interface compatibility; splitting is
 * applied beforehand by apply_lcp_splitting.
 */
inline std::tuple<
    std::vector<std::array<ulint, static_cast<size_t>(LCPSpillRunCols::COUNT)>>,
    std::vector<SpilloverVector>, ulint, ulint, size_t, size_t>
build_spill_data(const std::vector<std::vector<ulint>>& lcps_per_run,
                 double percentile_k = 0.98,
                 bool coalesce_spillover = false,
                 bool coalesce_lcp_separately = false,
                 ulint split_threshold = SPLIT_THRESHOLD_NEVER,
                 ulint spill_align = 0,
                 uchar spill_split_bits = 0,
                 bool minima_only = false)
{
    (void)split_threshold;
    return build_spill_data_from_pairs(retained_lcp_pairs(lcps_per_run, minima_only), percentile_k,
                                       coalesce_spillover, coalesce_lcp_separately, spill_align,
                                       spill_split_bits);
}

enum class LCPRunCols { TOP_LCP, COUNT };
template <bool SP> using MSIndexTopLCP = orbit::rlbwt::lf_permutation<LCPRunCols, true, SP>;

/**
 * The main index class that combines the run permutation and spillover data.
 */
template <bool StoreAbsolutePositions = false>
class MSIndexSpillLCP {
    using IndexImpl = orbit::rlbwt::lf_permutation<LCPSpillRunCols, true, StoreAbsolutePositions>;
    IndexImpl idx_;
    std::vector<SpilloverVector> spill_vectors_;
    ulint max_lcp_top_ = 0, max_lcp_min_sub_ = 0;
    ulint spill_align_ = 1;
    uchar spill_split_bits_ = 0;
    std::array<bool, 256> occurs_{};  // occurs_[c]: byte c appears somewhere in the BWT
    std::array<uchar, 256> code_{};   // code_[c]: alphabet code of byte c if it occurs, else no_code

    void compute_codes() {
        for (size_t c = 0; c < 256; ++c) {
            code_[c] = no_code;
            if (!occurs_[c]) continue;
            const auto k = idx_.character_code(static_cast<uchar>(c));
            if (k) code_[c] = *k;
        }
    }
    void compute_occurs() {
        occurs_.fill(false);
        for (ulint i = 0; i < idx_.intervals(); ++i) occurs_[idx_.get_character(i)] = true;
        compute_codes();
    }

public:
    using Position = typename IndexImpl::position;
    using position = Position;
    MSIndexSpillLCP(const std::vector<uchar>& chars,
                    const std::vector<ulint>& lens,
                    const std::vector<std::array<ulint, static_cast<size_t>(LCPSpillRunCols::COUNT)>>& run_data,
                    std::vector<SpilloverVector> spill_vectors, ulint max_top, ulint max_sub,
                    ulint spill_align = 0, uchar spill_split_bits = 0)
        : idx_(chars, lens, orbit::NO_SPLITTING, run_data)
        , spill_vectors_(std::move(spill_vectors))
        , max_lcp_top_(max_top)
        , max_lcp_min_sub_(max_sub)
        , spill_align_(spill_align > 0 ? spill_align : 1)
        , spill_split_bits_(spill_split_bits)
    {
        for (size_t i = 0; i < chars.size(); ++i)
            if (lens[i] > 0) occurs_[chars[i]] = true;
        compute_codes();
    }

    /** An empty index, to be filled by load(). */
    MSIndexSpillLCP() = default;

    /**
     * Write the index: the Orbit move structure (including the integrated
     * LCP columns) in Orbit's own packed serialization, then the spillover
     * arrays and the scalar settings.  Integers are written in host byte
     * order, as Orbit does.
     */
    size_t serialize(std::ostream& out) {
        size_t bytes = idx_.serialize(out);
        auto put = [&](ulint v) { out.write(reinterpret_cast<const char*>(&v), sizeof(v)); bytes += sizeof(v); };
        put(max_lcp_top_);
        put(max_lcp_min_sub_);
        put(spill_align_);
        put(spill_split_bits_);
        put(spill_vectors_.size());
        for (const auto& v : spill_vectors_) {
            put(v.size());
            out.write(reinterpret_cast<const char*>(v.data()), static_cast<std::streamsize>(v.size()));
            bytes += v.size();
        }
        return bytes;
    }

    /** Read an index written by serialize().  Throws on malformed input. */
    void load(std::istream& in) {
        idx_.load(in);
        auto get = [&]() {
            ulint v = 0;
            in.read(reinterpret_cast<char*>(&v), sizeof(v));
            if (!in.good()) throw std::runtime_error("truncated ms index");
            return v;
        };
        max_lcp_top_ = get();
        max_lcp_min_sub_ = get();
        spill_align_ = get();
        spill_split_bits_ = static_cast<uchar>(get());
        const ulint arrays = get();
        const ulint expected = (spill_split_bits_ > 0) ? (ulint{1} << spill_split_bits_) : 1;
        if (arrays != expected || spill_align_ == 0) throw std::runtime_error("inconsistent ms index settings");
        spill_vectors_.assign(static_cast<size_t>(arrays), {});
        for (auto& v : spill_vectors_) {
            v.resize(static_cast<size_t>(get()));
            in.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(v.size()));
            if (!in.good() && !v.empty()) throw std::runtime_error("truncated ms spillover");
        }
        compute_occurs();
    }

    /** True if byte c occurs in the indexed text. */
    bool occurs(uchar c) const { return occurs_[c]; }

    ulint get_length(ulint i) const { return idx_.get_length(i); }
    ulint get_length(Position p) const { return idx_.get_length(p); }
    uchar get_character(ulint i) { return idx_.get_character(i); }
    uchar get_character(Position p) { return idx_.get_character(p); }
    Position LF(Position p) { return idx_.LF(p); }
    /**
     * LF in two halves (see Orbit's start_next and finish_next):
     * finish_LF(start_LF(p)) equals LF(p), and start_LF reads only p's row,
     * so the landing row can be prefetched in between.
     */
    Position start_LF(Position p) const { return idx_.start_next(p); }
    Position finish_LF(Position p) const { return idx_.finish_next(p); }
    /** Hint that the row of interval i will be read soon. */
    void prefetch(ulint i) const { idx_.prefetch(i); }
    /** Hint that the rows of intervals lo to hi, lo <= hi, will be read soon. */
    void prefetch_rows(ulint lo, ulint hi) const { idx_.prefetch_rows(lo, hi); }
    /**
     * Hint that the spillover record of row i will be read soon, given its
     * LCP_SPILL column so.  A row without a record prefetches the start of its
     * spillover array, which is cheaper than testing for it.
     */
    void prefetch_spill(ulint i, ulint so) const {
        const uchar* p = spillover_for_row(i).data() + spill_offset_bytes(so);
        ORBIT_PREFETCH(p);
        ORBIT_PREFETCH(p + 64);
    }
    Position first() { return idx_.first(); }
    Position last() { return idx_.last(); }
    ulint move_runs() const { return idx_.intervals(); }
    ulint domain() const { return idx_.domain(); }
    Position down(Position p) { return idx_.down(p); }
    Position up(Position p) { return idx_.up(p); }

    template <LCPSpillRunCols Col>
    ulint get(ulint i) const { return idx_.template get<Col>(i); }
    template <LCPSpillRunCols Col>
    ulint get(Position p) const { return idx_.template get<Col>(p); }

    /** The alphabet code of byte c, or no_code if c does not occur in the text. */
    uchar code(uchar c) const { return code_[c]; }
    static constexpr uchar no_code = 0xFF;

    /**
     * Row access through a local copy of the packed layout (see Orbit's
     * packed_matrix::reader), for relative positions.  A row is named by its
     * first bit, row(i).  Its character, LCP_TOP and LCP_MIN_SUB columns are
     * adjacent and read together with one load, tail(r), as are its pointer
     * and offset; fits() says whether the index's column widths allow that.
     */
    struct PackedAccess {
        decltype(std::declval<const IndexImpl&>().get_reader()) rows;
        static constexpr size_t len_col = IndexImpl::length_column();
        static constexpr size_t ptr_col = IndexImpl::pointer_column();
        static constexpr size_t off_col = IndexImpl::offset_column();
        static constexpr size_t chr_col = IndexImpl::character_column();
        static constexpr size_t top_col = IndexImpl::template data_column<LCPSpillRunCols::LCP_TOP>();
        static constexpr size_t sub_col = IndexImpl::template data_column<LCPSpillRunCols::LCP_MIN_SUB>();
        static constexpr size_t spill_col = IndexImpl::template data_column<LCPSpillRunCols::LCP_SPILL>();
        static_assert(ptr_col < off_col && chr_col < top_col && top_col < sub_col, "unexpected column order");
        using Handle = ulint;
        using Tail = ulint;
        bool fits() const {
            return rows.template span_fits<ptr_col, off_col>() && rows.template span_fits<chr_col, sub_col>();
        }
        Handle row(ulint i) const { return rows.row_start(i); }
        void prefetch(ulint i) const { rows.prefetch(i); }
        void prefetch_rows(ulint lo, ulint hi) const { rows.prefetch_rows(lo, hi); }
        ulint length(Handle r) const { return rows.template get_at<len_col>(r); }
        /** The pointer and offset columns. */
        std::pair<ulint, ulint> pointer_offset(Handle r) const {
            const ulint span = rows.template get_span<ptr_col>(r);
            return {rows.template extract_span<ptr_col, ptr_col>(span), rows.template extract_span<ptr_col, off_col>(span)};
        }
        Tail tail(Handle r) const { return rows.template get_span<chr_col>(r); }
        uchar code(Tail t) const { return static_cast<uchar>(rows.template extract_span<chr_col, chr_col>(t)); }
        ulint top(Tail t) const { return rows.template extract_span<chr_col, top_col>(t); }
        ulint sub(Tail t) const { return rows.template extract_span<chr_col, sub_col>(t); }
        ulint spill(Handle r) const { return rows.template get_at<spill_col>(r); }
    };
    /** Row access through the index's own column reads, for any index. */
    struct ColumnAccess {
        MSIndexSpillLCP* idx;  // get_character is not const in Orbit
        using Handle = ulint;
        using Tail = ulint;
        Handle row(ulint i) const { return i; }
        void prefetch(ulint i) const { idx->prefetch(i); }
        void prefetch_rows(ulint lo, ulint hi) const { idx->prefetch_rows(lo, hi); }
        ulint length(Handle i) const { return idx->get_length(i); }
        Tail tail(Handle i) const { return i; }
        uchar code(Tail i) const { return idx->code(idx->get_character(i)); }
        ulint top(Tail i) const { return idx->template get<LCPSpillRunCols::LCP_TOP>(i); }
        ulint sub(Tail i) const { return idx->template get<LCPSpillRunCols::LCP_MIN_SUB>(i); }
        ulint spill(Handle i) const { return idx->template get<LCPSpillRunCols::LCP_SPILL>(i); }
    };
    static constexpr bool packed_access_supported = !StoreAbsolutePositions;
    PackedAccess packed_access() const {
        static_assert(packed_access_supported, "packed access needs relative positions");
        return PackedAccess{idx_.get_reader()};
    }
    ColumnAccess column_access() const { return ColumnAccess{const_cast<MSIndexSpillLCP*>(this)}; }

    /** Spillover array for row i (uses X LSBs of i when spill_split_bits>0). */
    const SpilloverVector& spillover_for_row(ulint i) const {
        size_t b = (spill_split_bits_ > 0) ? (i & ((1ULL << spill_split_bits_) - 1)) : 0;
        return spill_vectors_[b];
    }
    /** Total spillover bytes (sum over all arrays). */
    size_t spillover_total_bytes() const {
        size_t tot = 0;
        for (const auto& v : spill_vectors_) tot += v.size();
        return tot;
    }
    const std::vector<SpilloverVector>& spillover_vectors() const { return spill_vectors_; }
    const SpilloverVector& spillover() const { return spill_vectors_[0]; }  // legacy; use spillover_for_row
    ulint max_lcp_top() const { return max_lcp_top_; }
    ulint max_lcp_min_sub() const { return max_lcp_min_sub_; }
    ulint spill_align() const { return spill_align_; }
    uchar spill_split_bits() const { return spill_split_bits_; }
    /** Convert stored spill_off to byte offset (multiply by spill_align); NO_SPILL gives 0. */
    size_t spill_offset_bytes(ulint so) const {
        static_assert(NO_SPILL == 0, "NO_SPILL must scale to offset 0");
        return static_cast<size_t>(so) * spill_align_;
    }

    /** Return per-column bit widths for the underlying PackedVector (length, pointer, offset, character, lcp_top, lcp_min_sub, lcp_spill). */
    const std::array<uchar, 7>& get_widths() const { return idx_.get_widths(); }

    /** Return the alphabet characters (unmapped) in index order. */
    std::vector<uchar> get_alphabet() const { return idx_.get_alphabet(); }
};

/**
 * The top LCP value of run i, whose row is r, with tail t, in the row access
 * acc (see MSIndexSpillLCP::PackedAccess).
 */
template <bool SP, class Access>
inline ulint boundary_lcp(const MSIndexSpillLCP<SP>& idx, const Access& acc, typename Access::Handle r,
                          typename Access::Tail t, ulint i) {
    const ulint top = acc.top(t);
    const ulint sub = acc.sub(t);
    if (top == idx.max_lcp_top() && sub == idx.max_lcp_min_sub()) {
        // Jumbo row: the value is the first of its spillover record.
        const ulint so = acc.spill(r);
        assert(so != NO_SPILL);
        return decode_uleb128(idx.spillover_for_row(i), idx.spill_offset_bytes(so)).first;
    }
    return top;
}

/** Compute the top LCP value for a run. */
template <bool SP>
inline ulint boundary_lcp(const MSIndexSpillLCP<SP>& idx, ulint i) {
    return boundary_lcp(idx, idx.column_access(), i, i, i);
}

/**
 * The minimum LCP value of run i, whose row is r, with tail t, in the row
 * access acc.
 */
template <bool SP, class Access>
inline ulint row_min_lcp(const MSIndexSpillLCP<SP>& idx, const Access& acc, typename Access::Handle r,
                         typename Access::Tail t, ulint i) {
    const ulint top = acc.top(t);
    const ulint sub = acc.sub(t);
    if (top == idx.max_lcp_top() && sub == idx.max_lcp_min_sub()) {
        // Jumbo row: the value is the second of its spillover record.
        const ulint so = acc.spill(r);
        assert(so != NO_SPILL);
        const auto& spill = idx.spillover_for_row(i);
        return decode_uleb128(spill, skip_uleb128(spill, idx.spill_offset_bytes(so))).first;
    }
    return top - sub;
}

/** Compute the minimum LCP value in the row. */
template <bool SP>
inline ulint row_min_lcp(const MSIndexSpillLCP<SP>& idx, ulint i) {
    return row_min_lcp(idx, idx.column_access(), i, i, i);
}

/**
 * Compute the minimum LCP value in a range within the run.
 * to_top: true  -> range [0, offset] (from top down to offset); include top boundary.
 *         false -> range [offset, run_len) (from offset down to bottom); do NOT include
 *                  top boundary (it is the LCP with the previous run, irrelevant when moving down).
 */
template <bool SP>
inline ulint range_min(const MSIndexSpillLCP<SP>& idx, ulint interval, ulint offset, bool to_top) {
    ulint m = LCP_GAP;
    ulint top = idx.template get<LCPSpillRunCols::LCP_TOP>(interval);
    ulint sub = idx.template get<LCPSpillRunCols::LCP_MIN_SUB>(interval);
    ulint so = idx.template get<LCPSpillRunCols::LCP_SPILL>(interval);
    const auto& spill = idx.spillover_for_row(interval);
    bool jumbo = (top == idx.max_lcp_top() && sub == idx.max_lcp_min_sub());
    bool include_boundary = to_top;  /* only for upward range; top LCP irrelevant when moving down */
    size_t p = 0;
    bool do_spill = false;
    if (jumbo && so != NO_SPILL) {
        p = idx.spill_offset_bytes(so);
        auto [tv, np] = decode_uleb128(spill, p);
        p = skip_uleb128(spill, np);
        if (include_boundary) m = std::min(m, tv);
        do_spill = true;
    } else {
        if (include_boundary && top != LCP_GAP) m = std::min(m, top);
        if (so == NO_SPILL) return m;
        p = idx.spill_offset_bytes(so);
        do_spill = true;
    }
    if (do_spill) {
        auto [n, np] = decode_uleb128(spill, p);
        p = np;
        ulint run_len = idx.get_length(interval);
        for (ulint k = 0; k < n && k <= run_len; ++k) {
            if (p >= spill.size()) break;
            auto [o, no] = decode_uleb128(spill, p);
            auto [vv, nv] = decode_uleb128(spill, no);
            p = nv;
            /* Downward: exclude o=0 (top LCP, boundary with previous run).
             * Exclude o=offset: LCP[o] is between position o and o-1; when at offset we traverse
             * o+1, o+2, ... so we need LCP at o>offset. */
            bool in_range = to_top ? (o <= offset) : (offset == 0 ? o > 0 : o > offset);
            if (in_range && vv != LCP_GAP) m = std::min(m, vv);
        }
    }
    return m;
}

template <bool SP>
inline ulint range_min_to_top(const MSIndexSpillLCP<SP>& idx, ulint interval, ulint offset) {
    return range_min(idx, interval, offset, true);
}
template <bool SP>
inline ulint range_min_to_bottom(const MSIndexSpillLCP<SP>& idx, ulint interval, ulint offset) {
    return range_min(idx, interval, offset, false);
}

/**
 * Minimum LCP values a reposition crosses inside its starting run, as a pair:
 * first, going up from offset, the values at offsets [0, offset] including
 * the top boundary; second, going down, the values at offsets past offset.
 * Equal to (range_min(..., true), range_min(..., false)) but decodes the
 * run's spillover record once.  r is the run's row in the row access acc.
 */
template <bool SP, class Access>
inline std::pair<ulint, ulint> range_min_both(const MSIndexSpillLCP<SP>& idx, const Access& acc, typename Access::Handle r,
                                              ulint interval, ulint offset) {
    ulint up = LCP_GAP, down = LCP_GAP;
    const auto t = acc.tail(r);
    const ulint top = acc.top(t);
    const ulint sub = acc.sub(t);
    const ulint so = acc.spill(r);
    const bool jumbo = (top == idx.max_lcp_top() && sub == idx.max_lcp_min_sub());
    size_t p;
    const auto& spill = idx.spillover_for_row(interval);
    if (jumbo && so != NO_SPILL) {
        p = idx.spill_offset_bytes(so);
        auto [tv, np] = decode_uleb128(spill, p);
        p = skip_uleb128(spill, np);
        up = tv;
    } else {
        if (top != LCP_GAP) up = top;
        if (so == NO_SPILL) return {up, down};
        p = idx.spill_offset_bytes(so);
    }
    auto [n, np] = decode_uleb128(spill, p);
    p = np;
    const ulint run_len = acc.length(r);
    for (ulint k = 0; k < n && k <= run_len; ++k) {
        if (p >= spill.size()) break;
        auto [o, no] = decode_uleb128(spill, p);
        auto [v, nv] = decode_uleb128(spill, no);
        p = nv;
        if (v == LCP_GAP) continue;
        if (o <= offset) up = std::min(up, v);
        else down = std::min(down, v);
    }
    return {up, down};
}

template <bool SP>
inline std::pair<ulint, ulint> range_min_both(const MSIndexSpillLCP<SP>& idx, ulint interval, ulint offset) {
    return range_min_both(idx, idx.column_access(), interval, interval, offset);
}

/**
 * Find the row to reposition to for character c: the nearest run above or
 * below (interval, offset) whose character is c, choosing the side with the
 * larger range-minimum LCP.  Returns that position (before the LF step) and
 * the LCP bound, or nullopt if c does not occur in the text.  The walks stop
 * at the first and last runs rather than wrapping.
 *
 * Going up to the last row of run u crosses the LCP values of the starting
 * run from its top down to offset and every value of the runs strictly
 * between, including their tops.  Going down to the first row of run d
 * crosses the starting run's values below offset, every value of the runs
 * in between, and d's top.  A run's row minimum includes its top, so the
 * runs in between need only row_min_lcp.  Each walk reads a run's row once,
 * for its character and, unless it is the target, its row minimum, through
 * the row access acc (see MSIndexSpillLCP::PackedAccess).
 */
template <bool SP, class Access>
inline std::optional<std::pair<typename MSIndexSpillLCP<SP>::position, ulint>>
reposition_target_impl(MSIndexSpillLCP<SP>& idx, const Access acc, ulint interval, ulint offset, uchar c) {
    using Position = typename MSIndexSpillLCP<SP>::position;
    assert(idx.get_character(interval) != c);
    // A character absent from the text can never be found by walking; say so
    // at once instead of scanning the whole index.
    if (!idx.occurs(c)) return std::nullopt;
    const ulint last_run = idx.move_runs() - 1;
    // The starting run's spillover record is read first, so that a cache
    // miss on it overlaps the walks.  LCP_GAP is above every LCP value, so
    // std::min skips it.
    auto [min_up, min_down] = range_min_both(idx, acc, acc.row(interval), interval, offset);
    const uchar code = idx.code(c);
    ulint u = interval, d = interval;
    bool found_up = false, found_down = false;
    MS_COUNT(repositions, 1);
    while (u > 0) {
        MS_COUNT(walk_up, 1);
        const auto r = acc.row(--u);
        const auto t = acc.tail(r);
        if (acc.code(t) == code) { found_up = true; break; }
        min_up = std::min(min_up, row_min_lcp(idx, acc, r, t, u));
    }
    while (d < last_run) {
        MS_COUNT(walk_down, 1);
        const auto r = acc.row(++d);
        const auto t = acc.tail(r);
        if (acc.code(t) == code) {
            found_down = true;
            min_down = std::min(min_down, boundary_lcp(idx, acc, r, t, d));
            break;
        }
        min_down = std::min(min_down, row_min_lcp(idx, acc, r, t, d));
    }
    if (!found_up && !found_down) return std::nullopt;

    /* For comparison and match_len: treat LCP_GAP as high (no cap); use domain as sentinel */
    const ulint eff_high = idx.domain();
    const ulint eff_min_up = (min_up != LCP_GAP) ? min_up : eff_high;
    const ulint eff_min_down = (min_down != LCP_GAP) ? min_down : eff_high;
    const bool go_up = found_up && (!found_down || eff_min_up >= eff_min_down);
    // up and down fill in the offset (and absolute position, if stored) of
    // the row next to the one given.
    Position from{};
    if (go_up) {
        from.interval = u + 1;
        return std::make_pair(idx.up(from), eff_min_up);
    }
    from.interval = d - 1;
    return std::make_pair(idx.down(from), eff_min_down);
}

/**
 * reposition_target for one read at a time, as in ms_query.  It reads
 * columns through the index rather than through packed access: with nothing
 * prefetched, more of the walks' cache misses overlap that way.
 */
template <bool SP>
inline std::optional<std::pair<typename MSIndexSpillLCP<SP>::position, ulint>>
reposition_target(MSIndexSpillLCP<SP>& idx, ulint interval, ulint offset, uchar c) {
    return reposition_target_impl(idx, idx.column_access(), interval, offset, c);
}

/**
 * Reposition with LCP: find the next position and LCP value for the given
 * character.  The position is the LF image of reposition_target's row.
 */
template <typename Index>
inline std::optional<std::pair<typename Index::position, ulint>>
reposition_with_lcp(Index& idx, ulint interval, ulint offset, uchar c) {
    auto t = reposition_target(idx, interval, offset, c);
    if (t) t->first = idx.LF(t->first);
    return t;
}

/**
 * Matching statistics: MS[i] = length of longest prefix of pattern[i..]
 * that is present in the index. Uses a single (row, offset) pair; at each
 * step matches or repositions (up/down by min LCP) then LF.
 */
template <bool SP>
inline std::vector<ulint> ms_query(MSIndexSpillLCP<SP>& idx, const std::string& pattern) {
    using Position = typename MSIndexSpillLCP<SP>::position;
    std::vector<ulint> out;
    if (pattern.empty()) return out;
    out.resize(pattern.size(), 0);
    Position pos = idx.first();
    ulint match_len = 0;
    for (size_t i = pattern.size(); i > 0; --i) {
        uchar c = static_cast<uchar>(pattern[i - 1]);
        if (idx.get_character(pos.interval) == c) {
            // Simple "case 1" LF step
            pos = idx.LF(pos);
            out[i - 1] = ++match_len;
            continue;
        }
        // Case 2; note that reposition_with_lcp does the LF
        auto opt = reposition_with_lcp(idx, pos.interval, static_cast<ulint>(pos.offset), c);
        if (!opt) {
            // c does not occur in the text (for example N): no match can
            // include this position, so the statistic is 0 and matching
            // restarts from the current row.
            out[i - 1] = match_len = 0;
            continue;
        }
        assert(opt->second != LCP_GAP);
        out[i - 1] = match_len = std::min(match_len, opt->second) + 1;
        pos = opt->first;
    }
    return out;
}

/**
 * Matching statistics for many patterns at once, giving the same results as
 * ms_query on each.  Up to k patterns are in flight.  Each LF step is split
 * so that after computing a pattern's landing interval the loop prefetches
 * that row and moves on to the next pattern, which lets the row reads of
 * different patterns overlap instead of waiting on one another.  A reposition
 * takes two visits: the first prefetches the rows around the current one
 * and the current row's spillover record, and the second runs
 * reposition_target and starts its LF step.  Rows are read through the row
 * access acc (see MSIndexSpillLCP::PackedAccess); with packed access, a
 * visit resolves the LF step and reads the landing row's columns itself.
 * out[j] receives the statistics for patterns[j].
 */
template <bool SP, class Access>
inline void ms_query_batch_impl(MSIndexSpillLCP<SP>& idx, const Access acc, const std::vector<std::string>& patterns,
                                size_t k, std::vector<std::vector<ulint>>& out) {
    constexpr bool packed = std::is_same_v<Access, typename MSIndexSpillLCP<SP>::PackedAccess>;
    using Position = typename MSIndexSpillLCP<SP>::position;
    struct Slot {
        const char* pat;       // the pattern being matched
        ulint* ms;             // its output array
        size_t i;              // pattern[i - 1] is the next character to match
        Position pos;          // unresolved position from start_LF, or a resolved one
        ulint match_len;
        bool repositioning;    // pos is resolved and the next visit repositions for pattern[i - 1]
    };
    // Rows prefetched on each side of the current row before a reposition.
    // Most repositions stop within this many rows; a wider window costs more
    // in prefetch instructions than it saves.
    constexpr ulint rep_window = 8;
    const ulint last_run = idx.move_runs() - 1;
    out.resize(patterns.size());
    if (k == 0) k = 1;
    // The alphabet codes of all bytes, kept local so that stores to the
    // output cannot force them to be reloaded.
    std::array<uchar, 256> codes;
    for (size_t c = 0; c < 256; ++c) codes[c] = idx.code(static_cast<uchar>(c));
    std::vector<Slot> slots(k);
    size_t next = 0;
    // Start the next nonempty pattern; returns false when none are left.
    auto start = [&](Slot& s) {
        while (next < patterns.size()) {
            const size_t j = next++;
            out[j].assign(patterns[j].size(), 0);
            if (patterns[j].empty()) continue;
            s = Slot{patterns[j].data(), out[j].data(), patterns[j].size(), idx.first(), 0, false};
            return true;
        }
        return false;
    };
    // Slots [slots.data(), end) are in flight.
    Slot* end = slots.data();
    while (end < slots.data() + k && start(*end)) ++end;
    while (end != slots.data()) {
        for (Slot* sp = slots.data(); sp < end;) {
            Slot& s = *sp;
            const uchar c = static_cast<uchar>(s.pat[s.i - 1]);
            if (s.repositioning) {
                auto opt = reposition_target_impl(idx, acc, s.pos.interval, static_cast<ulint>(s.pos.offset), c);
                assert(opt && opt->second != LCP_GAP);
                s.ms[s.i - 1] = s.match_len = std::min(s.match_len, opt->second) + 1;
                s.pos = idx.start_LF(opt->first);
                s.repositioning = false;
            } else {
                // Resolve the pending LF step and read the landing row.
                Position pos;
                typename Access::Handle r;
                if constexpr (packed) {
                    ulint i = s.pos.interval, o = s.pos.offset;
                    r = acc.row(i);
                    MS_COUNT(lf_steps, 1);
                    for (ulint len = acc.length(r); o >= len; len = acc.length(r)) {
                        MS_COUNT(lf_ff, 1);
                        o -= len;
                        r = acc.row(++i);
                    }
                    pos.interval = i;
                    pos.offset = o;
                } else {
                    pos = idx.finish_LF(s.pos);
                    r = acc.row(pos.interval);
                }
                if (acc.code(acc.tail(r)) == codes[c]) {
                    s.ms[s.i - 1] = ++s.match_len;
                    if constexpr (packed) {
                        const auto [ptr, off] = acc.pointer_offset(r);
                        s.pos.interval = ptr;
                        s.pos.offset = pos.offset + off;
                    } else {
                        s.pos = idx.start_LF(pos);
                    }
                } else if (codes[c] == MSIndexSpillLCP<SP>::no_code) {
                    // c does not occur in the text: the statistic is 0 and
                    // matching restarts from the current row.
                    s.ms[s.i - 1] = s.match_len = 0;
                    s.pos = pos;
                } else {
                    // Prefetch what the reposition reads first and do it on
                    // the next visit; this character is not consumed yet.
                    const ulint cur = pos.interval;
                    const ulint lo = cur > rep_window ? cur - rep_window : 0;
                    const ulint hi = std::min(cur + rep_window, last_run);
                    acc.prefetch_rows(lo, hi);
                    idx.prefetch_spill(cur, acc.spill(r));
                    s.pos = pos;
                    s.repositioning = true;
                    ++sp;
                    continue;
                }
            }
            acc.prefetch(s.pos.interval);
            if (--s.i > 0 || start(s)) ++sp;
            else s = *--end;
        }
    }
}

template <bool SP>
inline void ms_query_batch(MSIndexSpillLCP<SP>& idx, const std::vector<std::string>& patterns, size_t k,
                           std::vector<std::vector<ulint>>& out) {
    if constexpr (MSIndexSpillLCP<SP>::packed_access_supported) {
        const auto acc = idx.packed_access();
        if (acc.fits()) return ms_query_batch_impl(idx, acc, patterns, k, out);
    }
    ms_query_batch_impl(idx, idx.column_access(), patterns, k, out);
}

/**
 * Extract a string of the given length and given repositioning (case 2) rate
 * by taking a random walk in the index.  All pseudo-randomness comes from the
 * provided RNG.  We start at a uniformly random (run, offset) pair.  At each
 * step, with probability reposition_rate, repositions to an adjacent run (via
 * one of the other characters in the index alphabet) using LCP-guided up/down;
 * then takes an LF step.  The letters we LF on are recorded in a string and
 * returned.
 */
template <bool SP, typename RNG>
inline std::string extract_errory_string(
    MSIndexSpillLCP<SP>& idx,
    float reposition_rate,
    ulint length,
    RNG& rng)
{
    using Position = typename MSIndexSpillLCP<SP>::position;
    const ulint r = idx.move_runs();
    if (r == 0) return {};
    std::uniform_real_distribution<float> unif(0.0f, 1.0f);
    std::vector<uchar> sigma = idx.get_alphabet();
    ulint interval = static_cast<ulint>(rng()) % r;
    ulint len = idx.get_length(interval);
    ulint offset = (len > 0) ? (static_cast<ulint>(rng()) % len) : 0;
    Position pos{interval, offset};

    std::string ret;
    ret.reserve(length);
    while (ret.size() < length) {
        if (sigma.size() > 1 && unif(rng) < reposition_rate) {
            uchar cur_c = idx.get_character(pos.interval);
            uchar c;
            do {
                c = sigma[static_cast<size_t>(rng()) % sigma.size()];
            } while (c == cur_c);
            auto opt = reposition_with_lcp(idx, pos.interval,
                                           static_cast<ulint>(pos.offset), c);
            if (opt) pos = opt->first;
        }
        ret += static_cast<char>(idx.get_character(pos.interval));
        pos = idx.LF(pos);
    }
    std::reverse(ret.begin(), ret.end());
    return ret;
}

#endif /* _MS_RLBWT_HPP */
