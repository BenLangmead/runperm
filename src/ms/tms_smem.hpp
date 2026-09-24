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
 */

#ifndef _TMS_SMEM_HPP
#define _TMS_SMEM_HPP

#include "tms_index.hpp"
#include "tms_query.hpp"
#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

/** What tms-batch reports for each read. */
enum class TmsReport {
    MS,        // matching statistics only
    SMEM_ONE,  // plus one SA entry per SMEM
    SMEM_ALL,  // plus every SA entry of each SMEM
};

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
    // Up, collected bottom to top, then reversed.
    auto p = idx.phi_at(x);
    while (idx.plcp(p) >= len) {
        p = idx.phi(p);
        out.push_back(p.idx);
    }
    std::reverse(out.begin() + static_cast<std::ptrdiff_t>(first), out.end());
    out.push_back(x);
    auto q = idx.phi_inv_at(x);
    while (idx.plcpb(q) >= len) {
        q = idx.phi_inv(q);
        out.push_back(q.idx);
    }
    return out.size() - first;
}

/** One read's SMEMs: for each, its start, length and count of positions. */
struct TmsSmemHits {
    struct Smem { ulint start, len, count; };
    std::vector<Smem> smems;
    std::vector<ulint> pos;  // the positions of each SMEM in turn
};

/**
 * The SMEMs longer than min_len of a read with statistics ms and toehold
 * positions pos (from tms_query_batch), with one position each for
 * SMEM_ONE and all of them for SMEM_ALL.
 */
inline void tms_report_smems(TmsIndex& idx, const std::vector<ulint>& ms, const std::vector<ulint>& pos,
                             ulint min_len, TmsReport report, TmsSmemHits& hits) {
    hits.smems.clear();
    hits.pos.clear();
    for (size_t i : tms_smems(ms, min_len)) {
        ulint count = 1;
        if (report == TmsReport::SMEM_ALL) count = tms_interval_positions(idx, pos[i], ms[i], hits.pos);
        else hits.pos.push_back(pos[i]);
        hits.smems.push_back({i, ms[i], count});
    }
}

/** Append hits to line in tms-batch's format: space-separated i:L:p, or i:L:c:p1,p2,... for SMEM_ALL. */
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
        for (ulint c = 0; c < h.count; ++c) {
            if (c > 0) line += ',';
            line += std::to_string(hits.pos[at++]);
        }
    }
}

#endif /* _TMS_SMEM_HPP */
