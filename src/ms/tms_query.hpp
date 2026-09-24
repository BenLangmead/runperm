/**
 * Matching statistics on a TmsIndex (see tms_index.hpp), computing each
 * repositioning LCE with psi.
 *
 * The state is the current LF position (row j) and the match length len,
 * with the invariant that the pattern characters just consumed,
 * P[i .. i + len), are the first len characters of row j's suffix.  On a
 * mismatch at P[i - 1] = c, the candidates are the tail of the nearest run of
 * c above j and the head of the nearest run of c below j.  The new length is
 * min(LCE(j, j'), len) + 1 for the candidate j' with the larger LCE, ties
 * going to the one above.  Because only min(LCE, len) matters, the LCE is
 * computed by comparing P[i ..] with row j''s suffix read forward with psi,
 * stopping at a mismatch or after len characters.
 */

#ifndef _TMS_QUERY_HPP
#define _TMS_QUERY_HPP

#include "tms_index.hpp"
#include <algorithm>
#include <optional>
#include <string>
#include <vector>

/**
 * Number of leading characters of pat[0 .. m) that the suffix at FL point q
 * shares, capped at cap.  At most cap - 1 psi steps.
 */
inline ulint tms_psi_lce(TmsIndex& idx, TmsIndex::FLPos q, const char* pat, size_t m, ulint cap) {
    ulint lce = 0;
    while (lce < cap && lce < m && idx.psi_character(q) == static_cast<uchar>(pat[lce])) {
        if (++lce == cap) break;
        q = idx.psi(q);
    }
    return lce;
}

/**
 * Reposition from LF position (interval, offset) for character c, with the
 * matched pattern characters rest[0 .. len) following it.  Returns the row
 * to take the LF step from and the capped LCE, or nullopt if c does not occur
 * in the text.  The scans stop at the first and last intervals rather than
 * wrapping.
 */
inline std::optional<std::pair<TmsIndex::LFPos, ulint>>
tms_reposition(TmsIndex& idx, ulint interval, uchar c, const char* rest, ulint len) {
    if (!idx.occurs(c)) return std::nullopt;
    const ulint last_run = idx.move_runs() - 1;
    ulint u = interval, d = interval;
    bool found_up = false, found_down = false;
    while (u > 0)
        if (idx.get_character(--u) == c) { found_up = true; break; }
    while (d < last_run)
        if (idx.get_character(++d) == c) { found_down = true; break; }
    if (!found_up && !found_down) return std::nullopt;

    ulint up_lce = 0, down_lce = 0;
    if (found_up && len > 0) up_lce = tms_psi_lce(idx, idx.psi_at_tail(u), rest, len, len);
    // succ can only win with a strictly larger LCE, which is impossible once
    // pred has reached the cap.
    if (found_down && len > 0 && (!found_up || up_lce < len))
        down_lce = tms_psi_lce(idx, idx.psi_at_head(d), rest, len, len);
    TmsIndex::LFPos from{};
    if (found_up && (!found_down || up_lce >= down_lce)) {
        from.interval = u + 1;
        return std::make_pair(idx.up(from), up_lce);
    }
    from.interval = d - 1;
    return std::make_pair(idx.down(from), down_lce);
}

/** Matching statistics of pattern against the index; same result as ms_query. */
inline std::vector<ulint> tms_query(TmsIndex& idx, const std::string& pattern) {
    std::vector<ulint> out(pattern.size(), 0);
    if (pattern.empty()) return out;
    const char* pat = pattern.data();
    const size_t m = pattern.size();
    TmsIndex::LFPos pos = idx.first();
    ulint len = 0;
    for (size_t i = m; i > 0; --i) {
        const uchar c = static_cast<uchar>(pat[i - 1]);
        if (idx.get_character(pos.interval) != c) {
            auto t = tms_reposition(idx, pos.interval, c, pat + i, std::min<ulint>(len, m - i));
            if (!t) {
                // c does not occur in the text: the statistic is 0 and
                // matching restarts from the current row.
                out[i - 1] = len = 0;
                continue;
            }
            pos = t->first;
            len = t->second;
        }
        pos = idx.LF_step(pos);
        out[i - 1] = ++len;
    }
    return out;
}

/** Matching statistics for many patterns; out[j] receives those of patterns[j]. */
inline void tms_query_batch(TmsIndex& idx, const std::vector<std::string>& patterns, size_t k,
                            std::vector<std::vector<ulint>>& out) {
    (void)k;
    out.resize(patterns.size());
    for (size_t j = 0; j < patterns.size(); ++j) out[j] = tms_query(idx, patterns[j]);
}

#endif /* _TMS_QUERY_HPP */
