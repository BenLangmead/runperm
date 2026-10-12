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
#include "step_stats.hpp"
#include "ftab.hpp"
#include <algorithm>
#include <array>
#include <cstdint>
#include <istream>
#include <ostream>
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
                  lf_ff = 0, lf_steps = 0, scan_visits = 0, walk_visits = 0;
                  TmsStats& operator+=(const TmsStats& o) {
                      bases += o.bases; repositions += o.repositions; psi_steps += o.psi_steps; phi_steps += o.phi_steps;
                      scan_rows += o.scan_rows; dist += o.dist; len_at_rep += o.len_at_rep; lce += o.lce;
                      lce_capped += o.lce_capped; dist1 += o.dist1; lf_ff += o.lf_ff; lf_steps += o.lf_steps;
                      scan_visits += o.scan_visits; walk_visits += o.walk_visits;
                      return *this;
                  } };
// Each query thread counts into its own copy; the batch driver adds them up.
inline thread_local TmsStats tms_stats;
#define TMS_COUNT(field, v) (tms_stats.field += (v))
#else
#define TMS_COUNT(field, v) ((void)0)
#endif
#if defined(MS_STEP_STATS) && !defined(TMS_STATS)
#error "MS_STEP_STATS reads the TMS_STATS counters"
#endif

/**
 * Number of leading characters of pat[0 .. m) that the suffix at psi point q
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
    ORACLE,   // replay recorded outcomes: no candidate scan and no LCE (see TmsOracle)
    ORACLE_SCAN,  // replay recorded outcomes' side and LCE, but scan that side for the candidate
};

/**
 * The outcome of one reposition: the LF interval the match moves to (the
 * pred candidate's interval, whose tail row is taken, if it is above the
 * current interval, else the succ candidate's, whose head row is taken) and
 * the capped LCE, the match length before the LF step.
 */
struct TmsDecision {
    uint64_t interval;
    uint64_t len;
};

/**
 * Recorded reposition outcomes for a sequence of patterns, for TmsMode::ORACLE:
 * pattern j's, in the order its query makes them, are
 * decisions[offsets[j] .. offsets[j + 1]).  move_runs is the LF interval
 * count of the index they were recorded on, since intervals are numbered
 * per index (and per layout).
 *
 * File format (all integers u64, little-endian): the magic "TMSORC01"
 * (8 bytes), move_runs, the pattern count p, the decision count d, then p
 * per-pattern decision counts, then d (interval, len) pairs.
 */
struct TmsOracle {
    uint64_t move_runs = 0;
    std::vector<uint64_t> offsets{0};
    std::vector<TmsDecision> decisions;

    size_t patterns() const { return offsets.size() - 1; }

    /** Append the decisions of the patterns in rec, in order. */
    void append(const std::vector<std::vector<TmsDecision>>& rec) {
        for (const auto& v : rec) {
            decisions.insert(decisions.end(), v.begin(), v.end());
            offsets.push_back(decisions.size());
        }
    }

    void write(std::ostream& out) const {
        static_assert(sizeof(TmsDecision) == 16, "decisions are written as two u64s");
        auto put = [&](uint64_t x) { out.write(reinterpret_cast<const char*>(&x), sizeof(x)); };
        out.write("TMSORC01", 8);
        put(move_runs);
        put(patterns());
        put(decisions.size());
        for (size_t j = 0; j < patterns(); ++j) put(offsets[j + 1] - offsets[j]);
        out.write(reinterpret_cast<const char*>(decisions.data()),
                  static_cast<std::streamsize>(decisions.size() * sizeof(TmsDecision)));
        if (!out) throw std::runtime_error("failed writing the oracle");
    }

    void read(std::istream& in) {
        char magic[8];
        auto get = [&] {
            uint64_t x = 0;
            in.read(reinterpret_cast<char*>(&x), sizeof(x));
            return x;
        };
        in.read(magic, 8);
        if (!in || std::string(magic, 8) != "TMSORC01") throw std::runtime_error("not a tms oracle file");
        move_runs = get();
        const uint64_t p = get(), d = get();
        if (!in) throw std::runtime_error("truncated oracle file");
        offsets.assign(1, 0);
        offsets.reserve(p + 1);
        for (uint64_t j = 0; j < p; ++j) offsets.push_back(offsets.back() + get());
        if (!in || offsets.back() != d) throw std::runtime_error("oracle file's per-pattern counts do not add up");
        decisions.resize(d);
        in.read(reinterpret_cast<char*>(decisions.data()), static_cast<std::streamsize>(d * sizeof(TmsDecision)));
        if (!in) throw std::runtime_error("truncated oracle file");
    }
};

/**
 * Part of a TmsOracle that a batch replays: its pattern j's decisions are
 * decisions[offsets[j] .. offsets[j + 1]).
 */
struct TmsOracleView {
    const TmsDecision* decisions;
    const uint64_t* offsets;
    uint64_t move_runs;
    /** Patterns [first, first + count) of oracle, which must hold them. */
    static TmsOracleView of(const TmsOracle& o, size_t first, size_t count) {
        if (first + count > o.patterns())
            throw std::runtime_error("the oracle holds fewer patterns than the query has");
        return TmsOracleView{o.decisions.data(), o.offsets.data() + first, o.move_runs};
    }
};

namespace tms_detail {
// Per-pattern oracle state of the batched engine, as a base of its slot so
// that it takes no room when neither recording nor replaying.
struct NoOracleSlot {};
struct ReplaySlot { const TmsDecision *orc, *orc_end; };
// ORACLE_SCAN also keeps the outcome being replayed while its side is scanned.
struct ReplayScanSlot : ReplaySlot { TmsDecision orc_now; };
struct RecordSlot { std::vector<TmsDecision>* rec; };
}  // namespace tms_detail

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
 * an index with phi, and PSI, PHISKIP and DUAL one with psi.  PhiStarts
 * says which form of phi the walks and the toehold use; positions need the
 * starts-based form.  With an ftab (ftab.hpp) built for idx, a pattern
 * whose last k-mer occurs starts from the table's state, k steps in; a
 * mode that walks phi, or positions, needs a table with a toehold, and
 * positions need its window positions.
 *
 * With Record, (*record)[j] receives every reposition outcome of patterns[j]
 * in order.  ORACLE replays such outcomes from *replay instead of scanning
 * for the candidates and computing their LCEs: on a mismatch it takes the
 * next recorded interval and length, prefetches that interval's LF rows,
 * and on the next visit moves there as the other modes do once they have
 * chosen, so it reads only LF (and phi for positions).  ORACLE_SCAN replays
 * the same outcomes with less given: the side is whichever of pred and succ
 * the recorded interval lies on, and that side alone is scanned for the
 * candidate as the other modes scan it.  The recorded length still stands in
 * for the LCE.  The interval the scan finds must be the recorded one.
 * Recording and replaying take no ftab, since a table skips the
 * repositions of a pattern's last k characters.
 */
template <TmsMode Mode, bool Positions, bool PhiStarts, bool Record, class Access>
inline void tms_query_batch_impl(TmsIndex& idx, const Access acc, const std::vector<std::string>& patterns, size_t k,
                                 std::vector<std::vector<ulint>>& out_len,
                                 std::vector<std::vector<ulint>>* out_pos, const ftab::Ftab* ft = nullptr,
                                 const TmsOracleView* replay = nullptr,
                                 std::vector<std::vector<TmsDecision>>* record = nullptr) {
    using LFPos = TmsIndex::LFPos;
    using FLPos = TmsIndex::FLPos;
    using PhiPos = TmsIndex::PhiPosOf<PhiStarts>;
    static_assert(PhiStarts || !Positions, "positions need phi's starts-based form");
    constexpr bool ReplayScan = Mode == TmsMode::ORACLE_SCAN;
    constexpr bool Replay = Mode == TmsMode::ORACLE || ReplayScan;
    constexpr bool Scans = !Replay || ReplayScan;
    static_assert(!(Replay && Record), "the oracle mode replays outcomes rather than recording them");
    constexpr bool Toehold = Positions || (Mode != TmsMode::PSI && !Replay);
    constexpr ulint INF = std::numeric_limits<ulint>::max();
    enum : uint8_t { STEP, SCAN, WALK, PRED_POS };
    struct Side {
        FLPos q;          // psi walk: unresolved psi point of the next character
        PhiPos p;         // phi walk: unresolved phi point of the next LCP value
        ulint psi_lce;    // psi walk: characters matched so far
        ulint phi_min;    // phi walk: minimum LCP so far
        ulint phi_left;   // phi walk: LCP values still to read
        ulint value;      // the capped LCE, once done
        bool psi_on, phi_on, done, lost;
        bool active() const { return !done && !lost; }
    };
    using OracleSlot = std::conditional_t<
        ReplayScan, tms_detail::ReplayScanSlot,
        std::conditional_t<Replay, tms_detail::ReplaySlot,
                           std::conditional_t<Record, tms_detail::RecordSlot, tms_detail::NoOracleSlot>>>;
    struct Slot : OracleSlot {
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
#ifdef MS_STEP_STATS
        size_t m;         // the pattern's length
#endif
    };
    // LF rows prefetched on each side of the current row before a
    // reposition's scan, and how many more a scan prefetches on a side
    // before yielding when it has not found the character.
    constexpr ulint rep_window = TMS_REP_WINDOW;
    constexpr ulint scan_extend = TMS_SCAN_EXTEND;
    const ulint last_run = idx.move_runs() - 1;
    // What the access's LF and FL characters are for each byte, kept local
    // so that stores to the output cannot force them to be reloaded.
    constexpr bool packed = !std::is_same_v<Access, TmsIndex::ColumnAccess>;
    std::array<uchar, 256> lf_sym, fl_sym;
    for (size_t c = 0; c < 256; ++c) {
        lf_sym[c] = packed ? idx.lf_code(static_cast<uchar>(c)) : static_cast<uchar>(c);
        fl_sym[c] = packed ? idx.fl_code(static_cast<uchar>(c)) : static_cast<uchar>(c);
    }
    out_len.resize(patterns.size());
    if constexpr (Positions) out_pos->resize(patterns.size());
    if constexpr (Record) record->resize(patterns.size());
    if (k == 0) k = 1;
    PhiPos init_ph{};
    if constexpr (Toehold) init_ph = idx.resolve_phi(idx.template phi_at_head<PhiStarts>(0));
    std::vector<Slot> slots(k);
    size_t next = 0;
    // With an ftab, each pattern's key, and how many patterns ahead of the
    // one starting its table entry is prefetched.
    std::vector<uint64_t> keys;
    const size_t ahead = k;
    if (ft) {
        keys.resize(patterns.size());
        for (size_t j = 0; j < patterns.size(); ++j) keys[j] = ft->key(patterns[j].data(), patterns[j].size());
        for (size_t j = 0; j < std::min(ahead, patterns.size()); ++j) ft->prefetch(keys[j]);
    }
#ifdef MS_STEP_STATS
    // Steps are recorded only with one pattern in flight, so that each
    // visit's time belongs to one step.
    const bool record_steps = k == 1;
    const ulint no_rows = 0;
#endif
    auto start = [&](Slot& s) {
        while (next < patterns.size()) {
            const size_t j = next++;
            out_len[j].assign(patterns[j].size(), 0);
            if constexpr (Positions) (*out_pos)[j].assign(patterns[j].size(), TMS_NO_POS);
            if constexpr (Record) (*record)[j].clear();
            if (patterns[j].empty()) continue;
            s = Slot{};
            if constexpr (Replay) {
                s.orc = replay->decisions + replay->offsets[j];
                s.orc_end = replay->decisions + replay->offsets[j + 1];
            }
            if constexpr (Record) s.rec = &(*record)[j];
            s.pat = patterns[j].data();
            s.ms = out_len[j].data();
            if constexpr (Positions) s.ms_pos = (*out_pos)[j].data();
            s.i = patterns[j].size();
            s.pos = idx.first();
            s.ph = init_ph;
            s.state = STEP;
            if (ft) {
                if (j + ahead < patterns.size()) ft->prefetch(keys[j + ahead]);
                ftab::Entry e;
                if (ft->lookup(keys[j], e)) {
                    const unsigned kk = ft->k();
                    const size_t m = patterns[j].size();
                    for (unsigned d = 0; d < kk; ++d) out_len[j][m - 1 - d] = d + 1;
                    if constexpr (Toehold) {
                        s.ph.interval = e.toehold_interval;
                        s.ph.offset = e.toehold_offset;
                        if constexpr (PhiStarts) s.ph.idx = idx.phi_start(e.toehold_interval) + e.toehold_offset;
                    }
                    if constexpr (Positions) {
                        // The window's positions are those the walk over
                        // each suffix of the k-mer ends on.
                        for (unsigned d = 0; d + 1 < kk; ++d) s.ms_pos[m - 1 - d] = ft->window(d + 1, keys[j]);
                        s.ms_pos[m - kk] = s.ph.idx;
                    }
                    if (m == kk) continue;
                    s.pos.interval = e.interval;
                    s.pos.offset = e.offset;
                    s.len = kk;
                    s.i = m - kk;
                    acc.prefetch(e.interval);
                    if constexpr (Toehold) idx.template prefetch_phi<PhiStarts>(s.ph.interval);
                }
            }
#ifdef MS_STEP_STATS
            s.m = patterns[j].size();
            if (record_steps) step_rec.begin(s.pat, s.m, out_len[j].data());
#endif
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
    constexpr bool PsiWalks = Mode != TmsMode::PHI && !Replay, PhiWalks = Mode != TmsMode::PSI && !Replay;
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
            MS_FF((PhiStarts ? FF_PHI_STARTS : FF_PHI_LEN), p.interval - x.p.interval);
            x.phi_min = std::min(x.phi_min, idx.plcp(p));
            if (x.phi_min == 0 || --x.phi_left == 0) {
                x.done = true;
                x.value = std::min(x.phi_min, cap);
            } else {
                x.p = idx.start_phi(p);
                idx.template prefetch_phi<PhiStarts>(x.p.interval);
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
            if (s.ph.offset == 0)
                idx.template prefetch_phi<PhiStarts>(s.ph.interval > 0 ? s.ph.interval - 1 : idx.phi_intervals() - 1);
        }
        if constexpr (Positions) s.ms_pos[s.i - 1] = s.ph.idx;
        s.state = STEP;
    };
    while (end != slots.data()) {
        for (Slot* sp = slots.data(); sp < end;) {
            Slot& s = *sp;
#ifdef MS_STEP_STATS
            // Rows read are the reposition scans' rows.
            StepProbe probe(s.m - s.i, &tms_stats.repositions, &tms_stats.scan_rows, &no_rows);
#endif
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
                } else if constexpr (Replay) {
                    // The recorded outcome stands in for the scan and the
                    // LCE walks; the move to it is finished on the next visit.
                    if (s.orc == s.orc_end)
                        throw std::runtime_error("the oracle has fewer repositions for a pattern than its query makes");
                    const TmsDecision o = *s.orc++;
                    const ulint cur = pos.interval, v = o.interval;
                    if (v == cur || v > last_run || o.len > s.len)
                        throw std::runtime_error("an oracle decision does not fit the query; was it recorded "
                                                 "on this index and these patterns?");
                    if constexpr (ReplayScan) {
                        // Scan only the side the recorded interval is on,
                        // as the SCAN state does for both sides otherwise.
                        s.orc_now = o;
                        s.u = s.d = cur;
                        s.found_up = s.found_down = false;
                        s.dist_up = pos.offset + 1;
                        s.dist_down = acc.length(r) - pos.offset;
                        if (v < cur) {
                            s.lo = cur > rep_window ? cur - rep_window : 0;
                            s.hi = cur;
                            acc.prefetch_rows(s.lo, cur - 1);
                            s.scanned_up = false;
                            s.scanned_down = true;
                        } else {
                            s.lo = cur;
                            s.hi = std::min(cur + rep_window, last_run);
                            acc.prefetch_rows(cur + 1, s.hi);
                            s.scanned_up = true;
                            s.scanned_down = false;
                        }
                        s.pos = pos;
                        s.state = SCAN;
                    } else {
                        TMS_COUNT(repositions, 1);
                        TMS_COUNT(len_at_rep, s.len);
                        s.up = Side{};
                        s.down = Side{};
                        if (v < cur) {
                            s.u = v;
                            s.up.done = true;
                            s.up.value = o.len;
                            s.down.lost = true;
                            // Row v's length gives its tail, and with positions
                            // row v + 1's phi point gives the tail's.
                            acc.prefetch_rows(v, Toehold ? v + 1 : v);
                        } else {
                            s.d = v;
                            s.down.done = true;
                            s.down.value = o.len;
                            s.up.lost = true;
                            acc.prefetch(v);
                    }
                    s.pos = pos;
                    s.state = WALK;
                    }
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
#ifdef MS_STEP_STATS
                const ulint from = s.ph.interval;
#endif
                s.ph = idx.finish_phi(s.ph);
                MS_FF(FF_PHI_TOEHOLD, s.ph.interval - from);
                consume(s, start_LF(s.pos));
                advanced = true;
            } else {
                const ulint cap = s.len;
                bool chosen_now = false;  // ORACLE_SCAN: the scan just ended; the move waits a visit
                if (Scans && s.state == SCAN) {
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
                    if constexpr (ReplayScan) {
                        // The scanned side's candidate takes the recorded
                        // length, and the move to it is finished on the next
                        // visit, as in ORACLE.
                        const TmsDecision o = s.orc_now;
                        const bool up = o.interval < s.pos.interval;
                        if (!(up ? found_up && s.u == o.interval : found_down && s.d == o.interval))
                            throw std::runtime_error("the scan found another candidate than the oracle recorded; "
                                                     "was it recorded on this index and these patterns?");
                        s.up = Side{};
                        s.down = Side{};
                        if (up) {
                            s.up.done = true;
                            s.up.value = o.len;
                            s.down.lost = true;
                            acc.prefetch_rows(s.u, Toehold ? s.u + 1 : s.u);
                        } else {
                            s.down.done = true;
                            s.down.value = o.len;
                            s.up.lost = true;
                        }
                        s.state = WALK;
                        chosen_now = true;
                    } else {
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
                                if (s.up.psi_on) { s.up.q = acc.psi_at_tail(s.u); acc.prefetch_psi(s.up.q.interval); }
                                if (s.up.phi_on) s.up.p = s.ph;  // the current row is the lower one
                            }
                            if (found_down) {
                                setup(s.down, dist_down, cap);
                                if (s.down.psi_on) {
                                    // Row d and the row below it, which the
                                    // reads that follow from row d often reach.
                                    s.down.q = acc.psi_at_head(s.d);
                                    acc.prefetch_psi(s.down.q.interval);
                                    if (s.d < last_run) acc.prefetch_psi(s.d + 1);
                                }
                                if (s.down.phi_on) {
                                    s.down.p = idx.template phi_at_head<PhiStarts>(s.d);
                                    idx.template prefetch_phi<PhiStarts>(s.down.p.interval);
                                }
                            }
                        }
                        if constexpr (Toehold) {
                            // Rows the new toehold is read from, whichever side wins.
                            if (found_up)
                                idx.template prefetch_phi<PhiStarts>(idx.template phi_at_head<PhiStarts>(s.u + 1).interval);
                            if (found_down)
                                idx.template prefetch_phi<PhiStarts>(idx.template phi_at_head<PhiStarts>(s.d).interval);
                        }
                        s.state = WALK;
                    }
                } else if constexpr (!Replay) {
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
                if constexpr (Mode != TmsMode::PSI && !Replay) {
                    if (!s.up.lost && !s.down.lost) {
                        if (upper(s.down, cap) <= lower(s.up)) s.down.lost = true;
                        else if (upper(s.up, cap) < lower(s.down)) s.up.lost = true;
                    }
                }
                if (!chosen_now && !s.up.active() && !s.down.active()) {
                    TMS_COUNT(lce, std::max(s.up.done ? s.up.value : 0, s.down.done ? s.down.value : 0));
                    TMS_COUNT(lce_capped, (s.up.done && s.up.value == cap) || (s.down.done && s.down.value == cap));
                    LFPos from{};
                    if (s.up.done && (s.down.lost || s.up.value >= s.down.value)) {
                        from.interval = s.u + 1;
                        from = idx.up(from);
                        s.len = s.up.value;
                        if constexpr (Record) s.rec->push_back(TmsDecision{s.u, s.len});
                        if constexpr (Toehold) {
                            // The tail's text position is one phi step from
                            // the next interval head's.
                            s.ph = idx.start_phi(idx.template phi_at_head<PhiStarts>(s.u + 1));
                            idx.template prefetch_phi<PhiStarts>(s.ph.interval);
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
                        if constexpr (Record) s.rec->push_back(TmsDecision{s.d, s.len});
                        if constexpr (Toehold) s.ph = idx.resolve_phi(idx.template phi_at_head<PhiStarts>(s.d));
                        consume(s, start_LF(from));
                        advanced = true;
                    }
                }
            }
#ifdef MS_STEP_STATS
            // Before start() files this pattern and begins the next.
            probe.record();
#endif
            if constexpr (Replay) {
                if (advanced && s.i == 1 && s.orc != s.orc_end)
                    throw std::runtime_error("the oracle has more repositions for a pattern than its query makes");
            }
            if (!advanced || --s.i > 0 || start(s)) ++sp;
            else s = *--end;
        }
    }
#ifdef MS_STEP_STATS
    if (record_steps) step_rec.flush();
#endif
}

/**
 * tms_query_batch_impl with the mode and positions chosen at run time.  Rows
 * are read through packed access when the index's column widths allow it
 * and packed is true, else through column access.  With record, the
 * reposition outcomes are recorded there (any mode but ORACLE); ORACLE
 * replays them from replay, which must be recorded on an index with the
 * same LF intervals, for the same patterns.
 */
inline void tms_query_batch(TmsIndex& idx, const std::vector<std::string>& patterns, size_t k,
                            std::vector<std::vector<ulint>>& out_len,
                            TmsMode mode = TmsMode::PSI, std::vector<std::vector<ulint>>* out_pos = nullptr,
                            bool packed = true, const ftab::Ftab* ft = nullptr, const TmsOracleView* replay = nullptr,
                            std::vector<std::vector<TmsDecision>>* record = nullptr) {
    const bool oracle = mode == TmsMode::ORACLE || mode == TmsMode::ORACLE_SCAN;
    if (((mode != TmsMode::PSI && !oracle) || out_pos) && !idx.has_phi())
        throw std::invalid_argument("this mode or positions need an index built with phi");
    if (out_pos && !idx.has_phi_starts())
        throw std::invalid_argument("positions need phi's starts-based form, which this index was loaded without");
    if (mode != TmsMode::PHI && !oracle && !idx.has_psi()) throw std::invalid_argument("this mode needs an index with psi");
    if (ft && (oracle || record))
        throw std::invalid_argument("recording or replaying repositions takes no ftab, which skips some of them");
    if (ft && (mode != TmsMode::PSI || out_pos) && !ft->has_toehold())
        throw std::invalid_argument("this mode or positions need an ftab with a toehold");
    if (ft && out_pos && !ft->has_window()) throw std::invalid_argument("positions need an ftab with window positions");
    if (oracle && !replay) throw std::invalid_argument("the oracle mode needs recorded repositions to replay");
    if (!oracle && replay) throw std::invalid_argument("only the oracle mode replays recorded repositions");
    if (oracle && record) throw std::invalid_argument("the oracle mode cannot record repositions");
    if (oracle && replay->move_runs != idx.move_runs())
        throw std::invalid_argument("the oracle was recorded on an index with other LF intervals");
    auto run = [&](auto m, auto rec) {
        constexpr TmsMode M = decltype(m)::value;
        constexpr bool R = decltype(rec)::value;
        auto with = [&](const auto acc) {
            // Without positions, phi's lengths-based form when present: its
            // rows are narrower.  The oracle mode reads phi only for positions.
            if (out_pos) tms_query_batch_impl<M, true, true, R>(idx, acc, patterns, k, out_len, out_pos, ft, replay, record);
            else if constexpr (M == TmsMode::ORACLE || M == TmsMode::ORACLE_SCAN)
                tms_query_batch_impl<M, false, true, R>(idx, acc, patterns, k, out_len, nullptr, ft, replay, record);
            else if (idx.has_phi_lengths())
                tms_query_batch_impl<M, false, false, R>(idx, acc, patterns, k, out_len, nullptr, ft, replay, record);
            else tms_query_batch_impl<M, false, true, R>(idx, acc, patterns, k, out_len, nullptr, ft, replay, record);
        };
        const TmsIndex::PackedAccess pa = idx.packed_access();
        if (packed && pa.fits()) with(pa);
        else with(idx.column_access());
    };
    auto run_rec = [&](auto m) {
        if (record) run(m, std::true_type{});
        else run(m, std::false_type{});
    };
    switch (mode) {
        case TmsMode::PSI: run_rec(std::integral_constant<TmsMode, TmsMode::PSI>{}); break;
        case TmsMode::PHI: run_rec(std::integral_constant<TmsMode, TmsMode::PHI>{}); break;
        case TmsMode::PHISKIP: run_rec(std::integral_constant<TmsMode, TmsMode::PHISKIP>{}); break;
        case TmsMode::DUAL: run_rec(std::integral_constant<TmsMode, TmsMode::DUAL>{}); break;
        case TmsMode::ORACLE: run(std::integral_constant<TmsMode, TmsMode::ORACLE>{}, std::false_type{}); break;
        case TmsMode::ORACLE_SCAN:
            run(std::integral_constant<TmsMode, TmsMode::ORACLE_SCAN>{}, std::false_type{});
            break;
    }
}

#endif /* _TMS_QUERY_HPP */
