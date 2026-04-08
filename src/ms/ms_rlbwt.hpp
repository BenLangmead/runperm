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
 * allows for more compression of interior LCPs.
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
#include <tuple>
#include <random>
#include <stdexcept>

using uchar = orbit::uchar;
using ulint = orbit::ulint;

/**
 * compress_lcps: replace interior LCP values that are >= max(boundary)
 *                with placeholders that will not be stored.
 * Appends to out: [top_lcp, compressed interior...] where gaps are LCP_GAP.
 */
inline void compress_lcps(const std::vector<ulint>& lcp,
                          const std::vector<ulint>* lcp_next,
                          std::vector<ulint>& out) {
    out.clear();
    if (lcp.empty()) return;
    ulint top_lcp = lcp[0];
    out.push_back(top_lcp);
    ulint lcp_top_bot_max = top_lcp;
    if (lcp_next && !lcp_next->empty())
        lcp_top_bot_max = std::max(lcp_top_bot_max, (*lcp_next)[0]);
    for (size_t i = 1; i < lcp.size(); ++i)
        out.push_back(lcp[i] < lcp_top_bot_max ? lcp[i] : LCP_GAP);
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
                                      const std::vector<ulint>* lcp_next)
{
    if (lcp.size() < 2) return 0;
    ulint m = lcp[0];
    if (lcp_next && !lcp_next->empty())
        m = std::max(m, (*lcp_next)[0]);
    size_t c = 0;
    for (size_t i = 1; i < lcp.size(); ++i)
        if (lcp[i] >= m) ++c;
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
                    ulint split_threshold = SPLIT_THRESHOLD_NEVER)
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
    size_t base_imp = compression_improvement(lcp, lcp_next);
    size_t top_imp = compression_improvement(top, &bot);
    size_t bot_imp = compression_improvement(bot, lcp_next);
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
                                 std::vector<std::pair<std::vector<ulint>, size_t>>& out)
{
    if (lcp.size() < 3) { out.emplace_back(lcp, orig_i); return; }
    auto [did, top, bot] = possibly_split_lcps(lcp, lcp_next, split_threshold);
    if (!did || bot.empty()) { out.emplace_back(lcp, orig_i); return; }
    recursive_split_lcps(top, &bot, orig_i, split_threshold, out);
    recursive_split_lcps(bot, lcp_next, orig_i, split_threshold, out);
}

}  // namespace detail

/**
 * Apply LCP splitting to the LCP vectors for each run.
 */
inline void apply_lcp_splitting(std::vector<uchar>& bwt_heads,
                               std::vector<ulint>& bwt_run_lengths,
                               std::vector<std::vector<ulint>>& lcps_per_run,
                               ulint split_threshold) {
    if (split_threshold == SPLIT_THRESHOLD_NEVER || lcps_per_run.empty()) return;
    const size_t n = lcps_per_run.size();
    std::vector<std::pair<std::vector<ulint>, size_t>> new_lcps;
    for (size_t i = 0; i < n; ++i) {
        const std::vector<ulint>* next = (i + 1 < n) ? &lcps_per_run[i + 1] : nullptr;
        detail::recursive_split_lcps(lcps_per_run[i], next, i, split_threshold, new_lcps);
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
 * Build the full spillover data for the index.
 * spill_align: 0=none, else align chunks to multiples of spill_align bytes.
 * spill_split_bits: X LSBs of row index select which spillover array (0=single array).
 * spill_off is stored as byte_offset / spill_align to save bits.
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
                 uchar spill_split_bits = 0)
{
    if (coalesce_lcp_separately && spill_align > 0)
        throw std::invalid_argument("coalescing LCPs seprately is not compatible with coalesce-spillover");
    const size_t r = lcps_per_run.size();
    size_t skinny_count = 0, jumbo_count = 0;
    std::vector<ulint> all_top, all_sub;
    std::vector<std::vector<std::pair<ulint, ulint>>> all_pairs(r);
    std::vector<ulint> full;

    for (size_t i = 0; i < r; ++i) {
        const std::vector<ulint>* next = (i + 1 < r) ? &lcps_per_run[i + 1] : nullptr;
        full.clear();
        compress_lcps(lcps_per_run[i], next, full);
        auto pairs = detail::compressed_to_pairs(full);
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
        all_pairs[i] = std::move(pairs);
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
                    size_t off_head = append_payload(bucket, p_head);
                    size_t off_lcp = append_payload(bucket, p_lcp);
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

enum class LCPRunCols { TOP_LCP, COUNT };
template <bool SP> using MSIndexTopLCP = orbit::rlbwt::RunPermLF<LCPRunCols, true, SP>;

/**
 * The main index class that combines the run permutation and spillover data.
 */
template <bool StoreAbsolutePositions = false>
class MSIndexSpillLCP {
    using IndexImpl = orbit::rlbwt::RunPermLF<LCPSpillRunCols, true, StoreAbsolutePositions>;
    IndexImpl idx_;
    std::vector<SpilloverVector> spill_vectors_;
    ulint max_lcp_top_, max_lcp_min_sub_;
    ulint spill_align_;
    uchar spill_split_bits_;

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
    {}

    ulint get_length(ulint i) const { return idx_.get_length(i); }
    ulint get_length(Position p) const { return idx_.get_length(p); }
    uchar get_character(ulint i) { return idx_.get_character(i); }
    uchar get_character(Position p) { return idx_.get_character(p); }
    Position LF(Position p) { return idx_.LF(p); }
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
    /** Convert stored spill_off to byte offset (multiply by spill_align). */
    size_t spill_offset_bytes(ulint so) const {
        if (so == NO_SPILL) return 0;
        return static_cast<size_t>(so) * spill_align_;
    }

    /** Return per-column bit widths for the underlying PackedVector (length, pointer, offset, character, lcp_top, lcp_min_sub, lcp_spill). */
    const std::array<uchar, 7>& get_widths() const { return idx_.get_widths(); }

    /** Return the alphabet characters (unmapped) in index order. */
    std::vector<uchar> get_alphabet() const { return idx_.get_alphabet(); }
};

/**
 * Compute the top LCP value for a run.
 */
template <bool SP>
inline ulint boundary_lcp(const MSIndexSpillLCP<SP>& idx, ulint i) {
    ulint top = idx.template get<LCPSpillRunCols::LCP_TOP>(i);
    ulint sub = idx.template get<LCPSpillRunCols::LCP_MIN_SUB>(i);
    bool jumbo = (top == idx.max_lcp_top() && sub == idx.max_lcp_min_sub());
    if (jumbo) {
        ulint so = idx.template get<LCPSpillRunCols::LCP_SPILL>(i);
        assert(so != NO_SPILL);
        const auto& spill = idx.spillover_for_row(i);
        size_t p = idx.spill_offset_bytes(so);
        auto [v, _] = decode_uleb128(spill, p);
        return v;
    }
    return top;
}

/**
 * Compute the minimum LCP value in the row.
 */
template <bool SP>
inline ulint row_min_lcp(const MSIndexSpillLCP<SP>& idx, ulint i) {
    ulint top = idx.template get<LCPSpillRunCols::LCP_TOP>(i);
    ulint sub = idx.template get<LCPSpillRunCols::LCP_MIN_SUB>(i);
    bool jumbo = (top == idx.max_lcp_top() && sub == idx.max_lcp_min_sub());
    if (jumbo) {
        ulint so = idx.template get<LCPSpillRunCols::LCP_SPILL>(i);
        assert(so != NO_SPILL);
        const auto& spill = idx.spillover_for_row(i);
        size_t p = idx.spill_offset_bytes(so);
        p = skip_uleb128(spill, p);
        // Grab the min LCP
        auto [m, _] = decode_uleb128(spill, p);
        return m;
    }
    return top - sub;
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
 * Reposition with LCP: find the next position and LCP value for the given
 * character.
 */
template <typename Index>
inline std::optional<std::pair<typename Index::position, ulint>>
reposition_with_lcp(Index& idx, ulint interval, ulint offset, uchar c) {
    using Position = typename Index::position;
    assert(idx.get_character(interval) != c);
    Position cur{interval, offset};
    const Position first_pos = idx.first();
    const Position last_pos = idx.last();
    const ulint max_walk = idx.domain();

    // Reposition up
    ulint min_up = LCP_GAP;
    Position pos_up = cur;
    bool found_up = false;
    for (ulint k = 0; k < max_walk; ++k) {
        if (idx.get_character(pos_up.interval) == c) {
            // Arrived at destination row
            ulint lcp_arrive = range_min(idx, pos_up.interval, pos_up.offset, false);
            if (lcp_arrive != LCP_GAP) min_up = std::min(min_up, lcp_arrive);
            found_up = true;
            break;
        }
        ulint lcp_val = k ? row_min_lcp(idx, pos_up.interval) :
                            range_min(idx, pos_up.interval, pos_up.offset, true);
        if (lcp_val != LCP_GAP) min_up = std::min(min_up, lcp_val);
        pos_up = idx.up(pos_up);
        if (pos_up.interval == last_pos.interval && pos_up.offset == last_pos.offset) break; /* wrapped */
    }

    // Reposition down
    ulint min_down = LCP_GAP;
    Position pos_down = cur;
    bool found_down = false;
    for (ulint k = 0; k < max_walk; ++k) {
        if (idx.get_character(pos_down.interval) == c) {
            // Arrived at destination row
            ulint lcp_arrive = range_min(idx, pos_down.interval, pos_down.offset, true);
            if (lcp_arrive != LCP_GAP) min_down = std::min(min_down, lcp_arrive);
            found_down = true;
            break;
        }
        ulint lcp_val = k ? row_min_lcp(idx, pos_down.interval) :
                            range_min(idx, pos_down.interval, pos_down.offset, false);
        if (lcp_val != LCP_GAP) min_down = std::min(min_down, lcp_val);
        pos_down = idx.down(pos_down);
        if (pos_down.interval == first_pos.interval && pos_down.offset == first_pos.offset) break; /* wrapped */
        ulint lcp_boundary = boundary_lcp(idx, pos_down.interval);
        if (lcp_boundary != LCP_GAP) min_down = std::min(min_down, lcp_boundary);
    }

    if (!found_up && !found_down) return std::nullopt;

    /* For comparison and match_len: treat LCP_GAP as high (no cap); use domain as sentinel */
    const ulint eff_high = idx.domain();
    ulint eff_min_up = (min_up != LCP_GAP) ? min_up : eff_high;
    ulint eff_min_down = (min_down != LCP_GAP) ? min_down : eff_high;

    ulint match_len = 0;
    Position out_pos;
    if (!found_up) {
        out_pos = pos_down;
        match_len = eff_min_down;
    } else if (!found_down) {
        out_pos = pos_up;
        match_len = eff_min_up;
    } else if (eff_min_up >= eff_min_down) {
        out_pos = pos_up;
        match_len = eff_min_up;
    } else {
        out_pos = pos_down;
        match_len = eff_min_down;
    }
    return std::make_pair(idx.LF(out_pos), match_len);
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
        if (!opt) return out;
        assert(opt->second != LCP_GAP);
        out[i - 1] = match_len = std::min(match_len, opt->second) + 1;
        pos = opt->first;
    }
    return out;
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
