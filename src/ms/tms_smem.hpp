/**
 * Super-maximal exact matches (SMEMs) from matching statistics on a
 * TmsIndex, and their suffix array entries.
 *
 * P[i .. i + L) is an SMEM when it occurs in the text and no other substring
 * of P that occurs in the text contains it.  With MS the matching statistics
 * of P, position i starts an SMEM, of length L = MS[i], exactly when
 * MS[i] > 0 and (i = 0 or MS[i - 1] <= MS[i]); otherwise P[i - 1 ..] matches
 * past i + MS[i] and contains it.
 *
 * The suffix array (SA) entries of an SMEM are the text positions where it
 * occurs, one per row of its BWT interval.  tms_query_batch with positions
 * gives one of them for each i (the toehold).  The rest come from walking
 * from that row up with phi while PLCP >= L and down with phi_inv while
 * PLCPB >= L, since the interval is the maximal block of rows around it whose
 * adjacent LCPs are all at least L.  Each walk starts from the toehold's
 * text position, located in phi or phi_inv by binary search.
 *
 * tms_report_smems_batch lists the SA entries of a whole block of reads at
 * once.  It treats every up walk and every down walk as an independent chain
 * of dependent phi or phi_inv steps and keeps up to K of them in flight,
 * visiting them round-robin: each visit finishes one step, reads PLCP or
 * PLCPB, starts the next step and prefetches the row it will read, so that
 * the K chains' cache misses overlap.
 */

#ifndef _TMS_SMEM_HPP
#define _TMS_SMEM_HPP

#include "tms_index.hpp"
#include "tms_query.hpp"
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

/** What tms-batch reports for each read. */
enum class TmsReport {
    MS,        // matching statistics only
    SMEM_ONE,  // plus one SA entry per SMEM
    SMEM_ALL,  // plus every SA entry of each SMEM
};

/** A max_listed that lists every position. */
constexpr ulint TMS_ALL_POSITIONS = std::numeric_limits<ulint>::max();

/** Starts i of the SMEMs of the pattern with statistics ms, keeping those with MS[i] > min_len. */
inline std::vector<size_t> tms_smems(const std::vector<ulint>& ms, ulint min_len = 0) {
    std::vector<size_t> out;
    for (size_t i = 0; i < ms.size(); ++i)
        if (ms[i] > min_len && (i == 0 || ms[i - 1] <= ms[i])) out.push_back(i);
    return out;
}

/**
 * Append to out the text positions of every row in the BWT interval of a
 * length-len match whose row holds text position x, in row order (top
 * first).  Needs an index with phi and phi_inv.  Returns how many it
 * appended.
 */
inline ulint tms_interval_positions(TmsIndex& idx, ulint x, ulint len, std::vector<ulint>& out) {
    const size_t first = out.size();
    // Append the positions of a walk from x while the LCP stays at least len.
    auto walk = [&](const auto& wk) {
        const auto p = wk.locate(x);
        if (wk.lcp(p) < len) return;
        auto q = wk.start_step(p);
        ulint lcp;
        do out.push_back(wk.step(q, lcp));
        while (lcp >= len);
    };
    // Up, collected bottom to top, then reversed.
    walk(idx.phi_walker());
    std::reverse(out.begin() + static_cast<std::ptrdiff_t>(first), out.end());
    out.push_back(x);
    walk(idx.phi_inv_walker());
    return out.size() - first;
}

/**
 * One read's SMEMs: for each, its start, length, count of positions, and
 * how many of those are listed (all of them unless a report caps it).
 */
struct TmsSmemHits {
    struct Smem { ulint start, len, count, listed; };
    std::vector<Smem> smems;
    std::vector<ulint> pos;  // the listed positions of each SMEM in turn
};

/**
 * The SMEMs longer than min_len of a read with statistics ms and toehold
 * positions pos (from tms_query_batch), with one position each for
 * SMEM_ONE and all of them for SMEM_ALL.  SMEM_ALL lists at most max_listed
 * positions per SMEM, the first ones in row order, but counts all of them.
 */
inline void tms_report_smems(TmsIndex& idx, const std::vector<ulint>& ms, const std::vector<ulint>& pos,
                             ulint min_len, TmsReport report, TmsSmemHits& hits,
                             ulint max_listed = TMS_ALL_POSITIONS) {
    hits.smems.clear();
    hits.pos.clear();
    for (size_t i : tms_smems(ms, min_len)) {
        ulint count = 1, listed = 1;
        if (report == TmsReport::SMEM_ALL) {
            count = tms_interval_positions(idx, pos[i], ms[i], hits.pos);
            listed = std::min(count, max_listed);
            hits.pos.resize(hits.pos.size() - (count - listed));
        } else {
            hits.pos.push_back(pos[i]);
        }
        hits.smems.push_back({i, ms[i], count, listed});
    }
}

namespace tms_smem_detail {

/** One SMEM of a block, with what its up and down walks found. */
struct Walks {
    ulint x, len;                       // toehold text position and SMEM length
    ulint up_count = 0, down_count = 0; // rows above and below x in the interval
    size_t up_at = 0, down_at = 0;      // where the kept positions start in found
    ulint up_kept = 0, down_kept = 0;   // how many positions each walk kept
    // The current walk's first step, from start_phi or start_phi_inv, if
    // it has one.
    TmsIndex::PhiPos first{};
    bool has_step = false;
};

/**
 * Run the up walks (phi, Down false) or the down walks (phi_inv, Down true)
 * of every SMEM in w, appending the positions each keeps to found.
 *
 * First, k walks at a time, locate each x among the interval starts: in
 * three passes over the k, prefetch the start table entries that bound x's
 * interval, read them and prefetch the rows between them, and binary search
 * those rows.  Each search ends by reading the LCP at x and starting the
 * walk's first step, if the LCP is at least the SMEM's length.
 *
 * Then walk, k at a time, round-robin: each visit takes a step with the
 * walker, which also reads the LCP there and starts the next step, and
 * prefetches the next step's row while the LCP is at least the length.  The
 * prefetch covers the row's cache line and the next one, where most
 * fast-forwards end.  A slot whose walk ends takes the next walk with a
 * step.
 *
 * An up walk keeps the max_listed positions nearest the top of the
 * interval, which it finds last.  A down walk runs after the up walks and
 * keeps the positions nearest x that fit in max_listed after the up walk's
 * and x.
 */
template <bool Down>
void run_walks(TmsIndex& idx, std::vector<Walks>& w, size_t k, ulint max_listed, std::vector<ulint>& found) {
    using Pos = TmsIndex::PhiPos;
    static_assert(std::is_same_v<TmsIndex::PhiPos, TmsIndex::PhiInvPos>, "phi and phi_inv share a position type");
    const auto wk = [&] {
        if constexpr (Down) return idx.phi_inv_walker();
        else return idx.phi_walker();
    }();

    // Search, k at a time: prefetch each x's start table entries, then
    // read them and prefetch the rows between them, then search those rows.
    constexpr size_t MAX_LINES = 4;  // rows prefetched per search, in cache lines
    std::vector<std::pair<ulint, ulint>> range(std::min(k, w.size()));
    for (size_t g = 0; g < w.size(); g += range.size()) {
        const size_t m = std::min(range.size(), w.size() - g);
        for (size_t t = 0; t < m; ++t) wk.prefetch_table(w[g + t].x);
        for (size_t t = 0; t < m; ++t) {
            auto& [lo, hi] = range[t];
            wk.range(w[g + t].x, lo, hi);
            wk.prefetch_range(lo, hi, MAX_LINES);
        }
        for (size_t t = 0; t < m; ++t) {
            Walks& e = w[g + t];
            const Pos p = wk.locate(e.x, range[t].first, range[t].second);
            e.has_step = wk.lcp(p) >= e.len;
            if (e.has_step) e.first = wk.start_step(p);
        }
    }

    // Walk.
    struct Slot {
        Pos p;        // the pending step
        size_t smem;  // which entry of w
        ulint len, count, budget;
        std::vector<ulint> buf;
    };
    // Record a finished walk in w and found.
    auto flush = [&](Slot& s) {
        Walks& e = w[s.smem];
        if constexpr (Down) {
            e.down_count = s.count;
            e.down_at = found.size();
            e.down_kept = s.buf.size();
            found.insert(found.end(), s.buf.begin(), s.buf.end());
        } else {
            e.up_count = s.count;
            e.up_kept = std::min<ulint>(s.count, max_listed);
            e.up_at = found.size();
            found.insert(found.end(), s.buf.end() - static_cast<std::ptrdiff_t>(e.up_kept), s.buf.end());
        }
    };
    size_t next = 0;
    // Load the next walk with a step to take into s.  Walks without one
    // keep their counts of 0.  Returns false once no walks are left.
    auto begin = [&](Slot& s) {
        while (next < w.size() && !w[next].has_step) ++next;
        if (next == w.size()) return false;
        s.smem = next++;
        const Walks& e = w[s.smem];
        s.p = e.first;
        wk.prefetch(s.p.interval);
        s.len = e.len;
        s.count = 0;
        s.buf.clear();
        if constexpr (Down) {
            const ulint up_listed = std::min(e.up_count, max_listed);
            s.budget = max_listed > up_listed ? max_listed - up_listed - 1 : 0;
        } else {
            s.budget = TMS_ALL_POSITIONS;
        }
        return true;
    };
    // Take s's pending step and start its next one.  Returns false once its
    // walk and every walk after it are done.
    auto visit = [&](Slot& s) {
        ulint lcp;
        const ulint x = wk.step(s.p, lcp);
        ++s.count;
        if (s.buf.size() < s.budget) s.buf.push_back(x);
        if (lcp >= s.len) {
            wk.prefetch(s.p.interval);
            return true;
        }
        flush(s);
        return begin(s);
    };
    std::vector<Slot> slots(std::min(k, w.size()));
    size_t live = 0;
    for (auto& s : slots)
        if (begin(s)) ++live;
        else break;
    // Round-robin over the live slots, moving each finished one past the end.
    while (live > 0)
        for (size_t s = 0; s < live;) {
            if (visit(slots[s])) ++s;
            else std::swap(slots[s], slots[--live]);
        }
}

}  // namespace tms_smem_detail

/**
 * tms_report_smems for reads 0 to ms.size() - 1 of a block, into hits[j].
 * For SMEM_ALL with k > 0, the up and down walks of all the block's SMEMs
 * run up to k at a time; the results are the same as with k = 0, which runs
 * tms_report_smems on each read in turn.
 */
inline void tms_report_smems_batch(TmsIndex& idx, const std::vector<std::vector<ulint>>& ms,
                                   const std::vector<std::vector<ulint>>& pos, ulint min_len, TmsReport report,
                                   size_t k, std::vector<TmsSmemHits>& hits,
                                   ulint max_listed = TMS_ALL_POSITIONS) {
    hits.resize(ms.size());
    if (report != TmsReport::SMEM_ALL || k == 0) {
        for (size_t j = 0; j < ms.size(); ++j) tms_report_smems(idx, ms[j], pos[j], min_len, report, hits[j], max_listed);
        return;
    }
    using tms_smem_detail::Walks;
    std::vector<Walks> w;
    for (size_t j = 0; j < ms.size(); ++j) {
        hits[j].smems.clear();
        hits[j].pos.clear();
        for (size_t i : tms_smems(ms[j], min_len)) {
            hits[j].smems.push_back({i, ms[j][i], 0, 0});
            w.push_back({pos[j][i], ms[j][i]});
        }
    }
    std::vector<ulint> found;
    tms_smem_detail::run_walks<false>(idx, w, k, max_listed, found);
    tms_smem_detail::run_walks<true>(idx, w, k, max_listed, found);
    // Row order: the up walk's positions reversed, then x, then the down
    // walk's.
    size_t at = 0;
    for (auto& h : hits)
        for (auto& s : h.smems) {
            const Walks& e = w[at++];
            s.count = e.up_count + 1 + e.down_count;
            s.listed = std::min(s.count, max_listed);
            const auto up = found.begin() + static_cast<std::ptrdiff_t>(e.up_at);
            h.pos.insert(h.pos.end(), std::make_reverse_iterator(up + static_cast<std::ptrdiff_t>(e.up_kept)),
                         std::make_reverse_iterator(up));
            if (s.listed > e.up_kept) h.pos.push_back(e.x);
            const auto down = found.begin() + static_cast<std::ptrdiff_t>(e.down_at);
            h.pos.insert(h.pos.end(), down, down + static_cast<std::ptrdiff_t>(e.down_kept));
        }
}

/**
 * Append hits to line in tms-batch's format: space-separated i:L:p, or
 * i:L:c:p1,p2,... for SMEM_ALL, with count c and the listed positions.
 */
inline void tms_format_smems(const TmsSmemHits& hits, TmsReport report, std::string& line) {
    size_t at = 0;
    for (size_t s = 0; s < hits.smems.size(); ++s) {
        const auto& h = hits.smems[s];
        if (s > 0) line += ' ';
        line += std::to_string(h.start);
        line += ':';
        line += std::to_string(h.len);
        line += ':';
        if (report == TmsReport::SMEM_ALL) {
            line += std::to_string(h.count);
            line += ':';
        }
        for (ulint c = 0; c < h.listed; ++c) {
            if (c > 0) line += ',';
            line += std::to_string(hits.pos[at++]);
        }
    }
}

#endif /* _TMS_SMEM_HPP */
