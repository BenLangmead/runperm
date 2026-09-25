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
#include <array>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <optional>
#include <string>
#include <vector>

#ifndef TMS_REP_WINDOW
#define TMS_REP_WINDOW 8
#endif
#ifndef TMS_SCAN_EXTEND
#define TMS_SCAN_EXTEND 16
#endif

#ifdef TMS_STATS
struct TmsStats { ulint bases = 0, repositions = 0, psi_steps = 0, phi_steps = 0, scan_rows = 0, dist = 0, len_at_rep = 0, lce = 0, lce_capped = 0, dist1 = 0,
                  lf_ff = 0, lf_steps = 0, scan_visits = 0, walk_visits = 0; };
inline TmsStats tms_stats;
#define TMS_COUNT(field, v) (tms_stats.field += (v))
#else
#define TMS_COUNT(field, v) ((void)0)
#endif

/**
 * Number of leading characters of pat[0 .. m) that the suffix at FL point q
 * shares, capped at cap.  At most cap - 1 psi steps.
 */
inline ulint tms_psi_lce(TmsIndex& idx, TmsIndex::FLPos q, const char* pat, size_t m, ulint cap) {
    ulint lce = 0;
    while (lce < cap && lce < m) {
        uchar c;
        const TmsIndex::FLPos next = idx.psi_step(q, c);
        if (c != static_cast<uchar>(pat[lce]) || ++lce == cap) break;
        q = next;
        TMS_COUNT(psi_steps, 1);
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
    TMS_COUNT(repositions, 1);
    TMS_COUNT(scan_rows, d - u);
    TMS_COUNT(len_at_rep, len);

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

/** Matching statistics of pattern against an index with psi; same result as ms_query. */
inline std::vector<ulint> tms_query(TmsIndex& idx, const std::string& pattern) {
    if (!idx.has_psi()) throw std::invalid_argument("tms_query needs an index with psi");
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
        TMS_COUNT(bases, 1);
    }
    return out;
}

/** How a reposition computes the LCE of the current row with a candidate. */
enum class TmsMode {
    PSI,      // compare the pattern with the candidate's suffix, read with psi
    PHI,      // walk phi across the rows between, reading PLCP
    PHISKIP,  // psi if the match length is below the row distance, else phi
    DUAL,     // psi if the match length is below the row distance, else both, first to finish wins
};

/** The position reported with a statistic of 0 (an absent character). */
constexpr ulint TMS_NO_POS = std::numeric_limits<ulint>::max();

/**
 * Matching statistics for many patterns at once.  Up to k patterns are in
 * flight, and every dependent memory access is split so that the row it
 * lands on is prefetched and the loop moves on to the next pattern before
 * reading it.  Each visit to a pattern does one of:
 *
 *  - finish its LF step and compare the next character, then start the next
 *    LF step (a match), restart (an absent character), or prefetch the rows
 *    around the current one for a reposition;
 *  - scan for the pred and succ candidates and set up their LCE walks,
 *    prefetching the rows where they start;
 *  - advance every running LCE walk by one step, prefetching the rows they
 *    land on;
 *  - for a pred reposition with positions, finish the phi step that gives
 *    the text position of the pred tail row.
 *
 * The pred and succ walks run side by side, and in DUAL each side can run a
 * psi and a phi walk at once.  A psi walk's count so far is a lower bound on
 * its side's capped LCE and a phi walk's minimum so far is an upper bound, so
 * a side stops as soon as it can no longer win (ties go to pred).  Every
 * mode therefore chooses the same candidate, and the lengths and positions
 * do not depend on the mode or on k.
 *
 * Rows are read through the row access acc (TmsIndex::PackedAccess or
 * ColumnAccess).  out_len[j] receives the statistics for patterns[j]; if positions is true,
 * out_pos[j] receives for each one a text position where it occurs
 * (TMS_NO_POS with a statistic of 0).  PHI, PHISKIP, DUAL and positions need
 * an index with phi, and PSI, PHISKIP and DUAL one with psi.
 */
template <TmsMode Mode, bool Positions, class Access>
inline void tms_query_batch_impl(TmsIndex& idx, const Access acc, const std::vector<std::string>& patterns, size_t k,
                                 std::vector<std::vector<ulint>>& out_len,
                                 std::vector<std::vector<ulint>>* out_pos) {
    using LFPos = TmsIndex::LFPos;
    using FLPos = TmsIndex::FLPos;
    using PhiPos = TmsIndex::PhiPos;
    constexpr bool Toehold = Positions || Mode != TmsMode::PSI;
    constexpr ulint INF = std::numeric_limits<ulint>::max();
    enum : uint8_t { STEP, SCAN, WALK, PRED_POS };
    struct Side {
        FLPos q;          // psi walk: unresolved FL point of the next character
        PhiPos p;         // phi walk: unresolved phi point of the next LCP value
        ulint psi_lce;    // psi walk: characters matched so far
        ulint phi_min;    // phi walk: minimum LCP so far
        ulint phi_left;   // phi walk: LCP values still to read
        ulint value;      // the capped LCE, once done
        bool psi_on, phi_on, done, lost;
        bool active() const { return !done && !lost; }
    };
    struct Slot {
        const char* pat;
        ulint* ms;
        ulint* ms_pos;
        size_t i;         // pat[i - 1] is the next character to match
        LFPos pos;        // unresolved after start_LF; resolved in the other states
        PhiPos ph;        // the toehold: the text position of pos's row
        ulint len;
        uint8_t state;
        ulint u, d;       // pred and succ candidate intervals, or how far their scans have got
        ulint lo, hi;     // the LF rows prefetched for the scans
        ulint dist_up, dist_down;  // rows from the current row to u's tail and d's head
        bool found_up, found_down, scanned_up, scanned_down;
        Side up, down;
    };
    // LF rows prefetched on each side of the current row before a
    // reposition's scan, and how many more a scan prefetches on a side
    // before yielding when it has not found the character.
    constexpr ulint rep_window = TMS_REP_WINDOW;
    constexpr ulint scan_extend = TMS_SCAN_EXTEND;
    const ulint last_run = idx.move_runs() - 1;
    // What the access's LF and FL characters are for each byte, kept local
    // so that stores to the output cannot force them to be reloaded.
    constexpr bool packed = std::is_same_v<Access, TmsIndex::PackedAccess>;
    std::array<uchar, 256> lf_sym, fl_sym;
    for (size_t c = 0; c < 256; ++c) {
        lf_sym[c] = packed ? idx.lf_code(static_cast<uchar>(c)) : static_cast<uchar>(c);
        fl_sym[c] = packed ? idx.fl_code(static_cast<uchar>(c)) : static_cast<uchar>(c);
    }
    out_len.resize(patterns.size());
    if constexpr (Positions) out_pos->resize(patterns.size());
    if (k == 0) k = 1;
    PhiPos init_ph{};
    if constexpr (Toehold) init_ph = idx.resolve_phi(idx.phi_at_head(0));
    std::vector<Slot> slots(k);
    size_t next = 0;
    auto start = [&](Slot& s) {
        while (next < patterns.size()) {
            const size_t j = next++;
            out_len[j].assign(patterns[j].size(), 0);
            if constexpr (Positions) (*out_pos)[j].assign(patterns[j].size(), TMS_NO_POS);
            if (patterns[j].empty()) continue;
            s = Slot{};
            s.pat = patterns[j].data();
            s.ms = out_len[j].data();
            if constexpr (Positions) s.ms_pos = (*out_pos)[j].data();
            s.i = patterns[j].size();
            s.pos = idx.first();
            s.ph = init_ph;
            s.state = STEP;
            return true;
        }
        return false;
    };
    // Slots [slots.data(), end) are in flight.
    Slot* end = slots.data();
    while (end < slots.data() + k && start(*end)) ++end;
    auto lower = [](const Side& x) { return x.done ? x.value : x.psi_on ? x.psi_lce : 0; };
    auto upper = [](const Side& x, ulint cap) { return x.done ? x.value : x.phi_on ? std::min(x.phi_min, cap) : cap; };
    auto setup = [&](Side& x, ulint dist, ulint cap) {
        if constexpr (Mode == TmsMode::PSI) x.psi_on = true;
        else if constexpr (Mode == TmsMode::PHI) x.phi_on = true;
        else if constexpr (Mode == TmsMode::PHISKIP) { x.psi_on = cap < dist; x.phi_on = !x.psi_on; }
        else { x.psi_on = true; x.phi_on = cap >= dist; }
        x.psi_lce = 0;
        x.phi_min = INF;
        x.phi_left = dist;
    };
    constexpr bool PsiWalks = Mode != TmsMode::PHI, PhiWalks = Mode != TmsMode::PSI;
    auto step_side = [&](Side& x, const char* rest, ulint cap) {
        if (PsiWalks && x.psi_on) {
            TMS_COUNT(psi_steps, 1);
            uchar ch;
            const FLPos next = acc.psi_step(x.q, ch);
            if (ch == fl_sym[static_cast<uchar>(rest[x.psi_lce])] && ++x.psi_lce < cap) {
                x.q = next;
                acc.prefetch_psi(x.q.interval);
            } else {
                x.done = true;
                x.value = x.psi_lce;
                return;
            }
        }
        if (PhiWalks && x.phi_on) {
            TMS_COUNT(phi_steps, 1);
            const PhiPos p = idx.finish_phi(x.p);
            x.phi_min = std::min(x.phi_min, idx.plcp(p));
            if (x.phi_min == 0 || --x.phi_left == 0) {
                x.done = true;
                x.value = std::min(x.phi_min, cap);
            } else {
                x.p = idx.start_phi(p);
                idx.prefetch_phi(x.p.interval);
            }
        }
    };
    // The LF step from resolved position p, unresolved.
    auto start_LF = [&](LFPos p) { return acc.start_LF(acc.row(p.interval), p); };
    // Record pat[i - 1]'s statistic, take the LF step started at next, and
    // move the toehold one text position left with it.
    auto consume = [&](Slot& s, LFPos next) {
        TMS_COUNT(bases, 1);
        s.ms[s.i - 1] = ++s.len;
        s.pos = next;
        acc.prefetch(s.pos.interval);
        if constexpr (Toehold) {
            s.ph = idx.phi_left(s.ph);
            // The next move left from offset 0 reads the row above.
            if (s.ph.offset == 0) idx.prefetch_phi(s.ph.interval > 0 ? s.ph.interval - 1 : idx.phi_intervals() - 1);
        }
        if constexpr (Positions) s.ms_pos[s.i - 1] = s.ph.idx;
        s.state = STEP;
    };
    while (end != slots.data()) {
        for (Slot* sp = slots.data(); sp < end;) {
            Slot& s = *sp;
            const uchar c = static_cast<uchar>(s.pat[s.i - 1]);
            bool advanced = false;  // pat[i - 1] is done
            if (s.state == STEP) {
                LFPos pos = s.pos;
                const auto r = acc.resolve_LF(pos);
                TMS_COUNT(lf_steps, 1);
                TMS_COUNT(lf_ff, pos.interval - s.pos.interval);
                if (acc.code(r) == lf_sym[c]) {
                    consume(s, acc.start_LF(r, pos));
                    advanced = true;
                } else if (!idx.occurs(c)) {
                    // Absent character: the statistic is 0 and matching
                    // restarts from the current row.
                    s.ms[s.i - 1] = s.len = 0;
                    s.pos = pos;
                    advanced = true;
                } else {
                    const ulint cur = pos.interval;
                    s.u = s.d = cur;
                    s.lo = cur > rep_window ? cur - rep_window : 0;
                    s.hi = std::min(cur + rep_window, last_run);
                    acc.prefetch_rows(s.lo, s.hi);
                    s.dist_up = pos.offset + 1;
                    s.dist_down = acc.length(r) - pos.offset;
                    s.found_up = s.found_down = s.scanned_up = s.scanned_down = false;
                    s.pos = pos;
                    s.state = SCAN;
                }
            } else if (s.state == PRED_POS) {
                s.ph = idx.finish_phi(s.ph);
                consume(s, start_LF(s.pos));
                advanced = true;
            } else {
                const ulint cap = s.len;
                if (s.state == SCAN) {
                    // Scan only rows prefetched on an earlier visit; a side
                    // that runs out of them prefetches more and yields.
                    TMS_COUNT(scan_visits, 1);
                    while (!s.scanned_up && s.u > s.lo) {
                        const auto r = acc.row(--s.u);
                        if (acc.code(r) == lf_sym[c]) s.found_up = s.scanned_up = true;
                        else s.dist_up += acc.length(r);
                    }
                    if (s.u == 0) s.scanned_up = true;
                    while (!s.scanned_down && s.d < s.hi) {
                        const auto r = acc.row(++s.d);
                        if (acc.code(r) == lf_sym[c]) s.found_down = s.scanned_down = true;
                        else s.dist_down += acc.length(r);
                    }
                    if (s.d == last_run) s.scanned_down = true;
                    if (!s.scanned_up || !s.scanned_down) {
                        if (!s.scanned_up) {
                            const ulint lo = s.lo > scan_extend ? s.lo - scan_extend : 0;
                            acc.prefetch_rows(lo, s.lo - 1);
                            s.lo = lo;
                        }
                        if (!s.scanned_down) {
                            const ulint hi = std::min(s.hi + scan_extend, last_run);
                            acc.prefetch_rows(s.hi + 1, hi);
                            s.hi = hi;
                        }
                        ++sp;
                        continue;
                    }
                    const bool found_up = s.found_up, found_down = s.found_down;
                    const ulint dist_up = s.dist_up, dist_down = s.dist_down;
                    TMS_COUNT(repositions, 1);
                    TMS_COUNT(scan_rows, s.d - s.u);
                    TMS_COUNT(len_at_rep, cap);
                    TMS_COUNT(dist, std::min(found_up ? dist_up : INF, found_down ? dist_down : INF));
                    TMS_COUNT(dist1, (found_up && dist_up == 1) || (found_down && dist_down == 1));
                    s.up = Side{};
                    s.down = Side{};
                    s.up.lost = !found_up;
                    s.down.lost = !found_down;
                    if (cap == 0) {
                        // Both LCEs cap at 0; take pred if there is one.
                        s.up.done = found_up;
                        s.down.done = found_down;
                        s.up.value = s.down.value = 0;
                    } else {
                        if (found_up) {
                            setup(s.up, dist_up, cap);
                            if (s.up.psi_on) { s.up.q = acc.psi_at_tail(acc.row(s.u)); acc.prefetch_psi(s.up.q.interval); }
                            if (s.up.phi_on) s.up.p = s.ph;  // the current row is the lower one
                        }
                        if (found_down) {
                            setup(s.down, dist_down, cap);
                            if (s.down.psi_on) {
                                // The head's FL point is one past the tail
                                // of the interval above (d > 0 here),
                                // unresolved: the first psi step may move on
                                // to the next FL row.
                                s.down.q = acc.psi_at_tail(acc.row(s.d - 1));
                                ++s.down.q.offset;
                                acc.prefetch_psi(s.down.q.interval);
                                if (s.down.q.interval + 1 < idx.psi_intervals()) acc.prefetch_psi(s.down.q.interval + 1);
                            }
                            if (s.down.phi_on) { s.down.p = idx.phi_at_head(s.d); idx.prefetch_phi(s.down.p.interval); }
                        }
                    }
                    if constexpr (Toehold) {
                        // Rows the new toehold is read from, whichever side wins.
                        if (found_up) idx.prefetch_phi(idx.phi_at_head(s.u + 1).interval);
                        if (found_down) idx.prefetch_phi(idx.phi_at_head(s.d).interval);
                    }
                    s.state = WALK;
                } else {
                    const char* rest = s.pat + s.i;
                    TMS_COUNT(walk_visits, 1);
                    if (s.up.active()) step_side(s.up, rest, cap);
                    if constexpr (Mode == TmsMode::PSI) {
                        // With psi counts alone, the only early decision is
                        // pred reaching the cap, since succ needs a strictly
                        // larger LCE; every other case waits for both counts.
                        if (s.up.done && s.up.value == cap) s.down.lost = true;
                        else if (s.down.active()) step_side(s.down, rest, cap);
                    } else {
                        if (s.down.active()) step_side(s.down, rest, cap);
                    }
                }
                if constexpr (Mode != TmsMode::PSI) {
                    if (!s.up.lost && !s.down.lost) {
                        if (upper(s.down, cap) <= lower(s.up)) s.down.lost = true;
                        else if (upper(s.up, cap) < lower(s.down)) s.up.lost = true;
                    }
                }
                if (!s.up.active() && !s.down.active()) {
                    TMS_COUNT(lce, std::max(s.up.done ? s.up.value : 0, s.down.done ? s.down.value : 0));
                    TMS_COUNT(lce_capped, (s.up.done && s.up.value == cap) || (s.down.done && s.down.value == cap));
                    LFPos from{};
                    if (s.up.done && (s.down.lost || s.up.value >= s.down.value)) {
                        from.interval = s.u + 1;
                        from = idx.up(from);
                        s.len = s.up.value;
                        if constexpr (Toehold) {
                            // The tail's text position is one phi step from
                            // the next interval head's.
                            s.ph = idx.start_phi(idx.phi_at_head(s.u + 1));
                            idx.prefetch_phi(s.ph.interval);
                            s.pos = from;
                            s.state = PRED_POS;
                        } else {
                            consume(s, start_LF(from));
                            advanced = true;
                        }
                    } else {
                        from.interval = s.d - 1;
                        from = idx.down(from);
                        s.len = s.down.value;
                        if constexpr (Toehold) s.ph = idx.resolve_phi(idx.phi_at_head(s.d));
                        consume(s, start_LF(from));
                        advanced = true;
                    }
                }
            }
            if (!advanced || --s.i > 0 || start(s)) ++sp;
            else s = *--end;
        }
    }
}

/**
 * tms_query_batch_impl with the mode and positions chosen at run time.  Rows
 * are read through packed access when the index's column widths allow it
 * and packed is true, else through column access.
 */
inline void tms_query_batch(TmsIndex& idx, const std::vector<std::string>& patterns, size_t k,
                            std::vector<std::vector<ulint>>& out_len,
                            TmsMode mode = TmsMode::PSI, std::vector<std::vector<ulint>>* out_pos = nullptr,
                            bool packed = true) {
    if ((mode != TmsMode::PSI || out_pos) && !idx.has_phi())
        throw std::invalid_argument("this mode or positions need an index built with phi");
    if (mode != TmsMode::PHI && !idx.has_psi()) throw std::invalid_argument("this mode needs an index with psi");
    const TmsIndex::PackedAccess pa = idx.packed_access();
    const bool use_packed = packed && pa.fits();
    auto run = [&](auto m) {
        constexpr TmsMode M = decltype(m)::value;
        if (use_packed) {
            if (out_pos) tms_query_batch_impl<M, true>(idx, pa, patterns, k, out_len, out_pos);
            else tms_query_batch_impl<M, false>(idx, pa, patterns, k, out_len, nullptr);
        } else {
            const TmsIndex::ColumnAccess ca = idx.column_access();
            if (out_pos) tms_query_batch_impl<M, true>(idx, ca, patterns, k, out_len, out_pos);
            else tms_query_batch_impl<M, false>(idx, ca, patterns, k, out_len, nullptr);
        }
    };
    switch (mode) {
        case TmsMode::PSI: run(std::integral_constant<TmsMode, TmsMode::PSI>{}); break;
        case TmsMode::PHI: run(std::integral_constant<TmsMode, TmsMode::PHI>{}); break;
        case TmsMode::PHISKIP: run(std::integral_constant<TmsMode, TmsMode::PHISKIP>{}); break;
        case TmsMode::DUAL: run(std::integral_constant<TmsMode, TmsMode::DUAL>{}); break;
    }
}

#endif /* _TMS_QUERY_HPP */
