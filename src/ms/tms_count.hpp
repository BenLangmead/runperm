/**
 * Counting the occurrences of an SMEM by backward search on a TmsIndex, and
 * listing its first positions in row order from the top of its BWT interval.
 *
 * tms_smem.hpp counts an SMEM of length L by walking phi up and phi_inv down
 * from the toehold while the LCP stays at least L, which takes one step per
 * occurrence.  Backward search instead finds the SMEM's BWT interval
 * [top, bot] in L extension steps, and the count is bot - top + 1.
 *
 * An extension by character c maps both ends with LF.  An end whose LF row
 * holds c is stepped as it is.  Otherwise top moves to the head of the
 * nearest run of c below it and bot to the tail of the nearest run of c above
 * it, found by scanning LF's rows as a matching statistics reposition does,
 * and those are stepped.  The interval starts as the F block of the SMEM's
 * last character, whose ends are the LF images of the head of the first run
 * of c and of the tail of the last one (TmsCounter keeps them per character).
 *
 * LF positions are (interval, offset) pairs, since LF's rows store lengths.
 * The count needs the rows' absolute numbers, or at least their difference:
 * when top and bot lie at most a block of intervals apart, the lengths of
 * the intervals between them are summed; otherwise each row's number comes
 * from a sampled table of interval starts (TmsRowStarts) plus the lengths of
 * at most a block of intervals before it.
 *
 * The text position of top comes almost free (the r-index toehold).  Where
 * top steps without moving, SA[LF(top)] = SA[top] - 1.  Where it moves to
 * the head of a run, SA there is the text position of that row's phi point,
 * which LF's PHI columns give at every run head.  So only the last move
 * matters: with h the run head top last moved to and s the extension steps
 * taken since, inclusive, top's text position is SA[h] - s, one phi row read
 * per SMEM.  The rest of a capped listing of N positions comes from N - 1
 * phi_inv steps from there, with no LCP tests, since the count is known.
 * Those are the first N positions in row order, the same as the walk lists.
 *
 * The batched forms keep up to k searches (or listings) in flight, visited
 * round-robin; each visit takes one extension (or phi_inv step) and
 * prefetches the rows the next one reads, so that cache misses overlap.
 */

#ifndef _TMS_COUNT_HPP
#define _TMS_COUNT_HPP

#include "tms_index.hpp"
#include <array>
#include <chrono>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

/** How tms_report_smems_batch counts the occurrences of an SMEM for SMEM_ALL. */
enum class TmsCount {
    WALK,      // walk phi and phi_inv from the toehold while the LCP is at least L
    BACKWARD,  // backward search for the SMEM's BWT interval
    ADAPTIVE,  // walk up to a budget of steps, then backward search if the walk has not ended
};

/** The counting method and the adaptive method's thresholds. */
struct TmsCountOptions {
    TmsCount method = TmsCount::WALK;
    // ADAPTIVE: walk steps (up and down together) before switching to
    // backward search.
    ulint budget = 64;
    // ADAPTIVE: SMEMs at most this long are backward searched without
    // walking first.
    ulint short_len = 0;
};

/** What counting did, summed over a thread's blocks. */
struct TmsCountStats {
    ulint smems = 0;            // SMEMs counted
    ulint walked = 0;           // counted by a walk
    ulint switched = 0;         // ADAPTIVE: walked to the budget, then searched
    ulint searched = 0;         // counted by backward search
    ulint walk_steps = 0;       // phi and phi_inv steps of the counting walks
    ulint search_steps = 0;     // backward extensions
    ulint top_moves = 0;        // extensions whose top moved to another run
    ulint bot_moves = 0;        // extensions whose bot moved to another run
    ulint scan_rows = 0;        // LF rows read by those moves' scans
    ulint lf_ff = 0;            // fast-forward rows of the ends' LF steps
    ulint far_ends = 0;         // counts that read the row-start table
    ulint sum_rows = 0;         // interval lengths summed for counts
    ulint list_steps = 0;       // phi_inv steps listing positions after a search
    double walk_s = 0, search_s = 0, list_s = 0;
    TmsCountStats& operator+=(const TmsCountStats& o) {
        smems += o.smems; walked += o.walked; switched += o.switched; searched += o.searched;
        walk_steps += o.walk_steps; search_steps += o.search_steps; top_moves += o.top_moves;
        bot_moves += o.bot_moves; scan_rows += o.scan_rows; lf_ff += o.lf_ff; far_ends += o.far_ends;
        sum_rows += o.sum_rows; list_steps += o.list_steps;
        walk_s += o.walk_s; search_s += o.search_s; list_s += o.list_s;
        return *this;
    }
};
// Each query thread counts into its own copy; the batch driver adds them up.
inline thread_local TmsCountStats tms_count_stats;

/**
 * The first row of every 2^shift-th LF interval: absolute at every
 * 2^SUPER-th interval, and relative to that in 32 bits at the others.
 * row(p) adds the lengths of the at most 2^shift - 1 intervals between the
 * sample and p.
 */
class TmsRowStarts {
public:
    static constexpr unsigned SUPER = 16;

    TmsRowStarts() = default;
    TmsRowStarts(const TmsIndex& idx, unsigned shift) : shift_(shift) {
        if (shift > SUPER) throw std::invalid_argument("the row-start sampling shift must be at most 16");
        const auto acc = idx.packed_access();
        const ulint m = idx.move_runs();
        rel_.reserve((m >> shift) + 1);
        abs_.reserve((m >> SUPER) + 1);
        ulint row = 0, base = 0;
        for (ulint i = 0; i < m; ++i) {
            if ((i & ((ulint(1) << SUPER) - 1)) == 0) {
                abs_.push_back(row);
                base = row;
            }
            if ((i & ((ulint(1) << shift) - 1)) == 0) {
                if (row - base > UINT32_MAX) throw std::runtime_error("LF rows too long for the row-start table");
                rel_.push_back(static_cast<uint32_t>(row - base));
            }
            row += acc.length(acc.row(i));
        }
    }
    unsigned shift() const { return shift_; }
    size_t bytes() const { return rel_.size() * sizeof(uint32_t) + abs_.size() * sizeof(ulint); }
    /** The first row of LF interval i. */
    ulint start(const TmsIndex::PackedAccess& acc, ulint i) const {
        const ulint b = i >> shift_;
        ulint r = abs_[i >> SUPER] + rel_[b];
        for (ulint j = b << shift_; j < i; ++j) r += acc.length(acc.row(j));
        tms_count_stats.sum_rows += i - (b << shift_);
        return r;
    }

private:
    unsigned shift_ = 0;
    std::vector<uint32_t> rel_;
    std::vector<ulint> abs_;
};

/**
 * What backward search needs beyond the index: the row-start table and,
 * for each character, the unresolved LF images of the head of its first run
 * and the tail of its last run, and the text position of the first row of
 * its F block.  Needs an index with phi's starts-based form.  Valid while
 * the index is alive and unchanged.
 */
class TmsCounter {
public:
    using LFPos = TmsIndex::LFPos;
    static constexpr unsigned DEFAULT_SHIFT = 4;

    TmsCounter() = default;
    explicit TmsCounter(const TmsIndex& idx, unsigned shift = DEFAULT_SHIFT) : rows_(idx, shift) {
        if (!idx.has_phi_starts()) throw std::invalid_argument("backward search counting needs phi's starts-based form");
        const auto acc = idx.packed_access();
        const ulint m = idx.move_runs();
        auto& ix = const_cast<TmsIndex&>(idx);  // Orbit's character reads are not const
        for (size_t c = 1; c < 256; ++c) {
            if (!idx.occurs(static_cast<uchar>(c))) continue;
            const uchar code = ix.lf_code(static_cast<uchar>(c));
            if (code == TmsIndex::no_code) continue;
            ulint d = 0, u = m - 1;
            while (d < m && acc.code(acc.row(d)) != code) ++d;
            while (u > 0 && acc.code(acc.row(u)) != code) --u;
            if (d == m) continue;
            Start& s = start_[code];
            s.present = true;
            LFPos h{};
            h.interval = d;
            s.top = acc.start_LF(acc.row(d), h);
            s.top_x = head_text_pos(idx, d) - 1;
            LFPos t{};
            t.interval = u;
            t.offset = acc.length(acc.row(u)) - 1;
            s.bot = acc.start_LF(acc.row(u), t);
        }
    }
    const TmsRowStarts& rows() const { return rows_; }
    size_t bytes() const { return rows_.bytes(); }

    /** The text position of the head row of LF interval d, which must start a BWT run. */
    static ulint head_text_pos(const TmsIndex& idx, ulint d) {
        return idx.resolve_phi(idx.template phi_at_head<true>(d)).idx;
    }

    struct Start {
        LFPos top{}, bot{};
        ulint top_x = 0;
        bool present = false;
    };
    const Start& start(uchar code) const { return start_[code]; }

private:
    TmsRowStarts rows_;
    std::array<Start, 256> start_{};
};

/** One SMEM to count by backward search: P[0 .. len) with P = pat. */
struct TmsSearchItem {
    const char* pat;
    ulint len;
    // Results: the count and the text position of the interval's top row.
    ulint count = 0, top_x = 0;
};

/**
 * Count each item by backward search, up to k at a time, and set its top_x.
 * Every item must occur in the text.
 */
inline void tms_backward_search_batch(const TmsIndex& idx, const TmsCounter& ctr, std::vector<TmsSearchItem>& items,
                                      size_t k) {
    using LFPos = TmsIndex::LFPos;
    if (items.empty()) return;
    auto& ix = const_cast<TmsIndex&>(idx);
    const auto acc = idx.packed_access();
    const ulint last = idx.move_runs() - 1;
    std::array<uchar, 256> code;
    for (size_t c = 0; c < 256; ++c) code[c] = ix.lf_code(static_cast<uchar>(c));
    TmsCountStats& st = tms_count_stats;
    struct Slot {
        LFPos top, bot;    // unresolved
        ulint left;        // characters still to extend by: pat[0 .. left)
        ulint head;        // the run head top last moved to, or last + 1 for none
        ulint since;       // LF steps of top since that run head (or the first run of the last character)
        size_t item;
    };
    size_t next = 0;
    auto begin = [&](Slot& s) {
        if (next == items.size()) return false;
        s.item = next++;
        const TmsSearchItem& it = items[s.item];
        const uchar c = code[static_cast<uchar>(it.pat[it.len - 1])];
        if (c == TmsIndex::no_code || !ctr.start(c).present)
            throw std::invalid_argument("backward search of a pattern that does not occur");
        const auto& b = ctr.start(c);
        s.top = b.top;
        s.bot = b.bot;
        s.left = it.len - 1;
        s.head = last + 1;
        s.since = 1;
        acc.prefetch(s.top.interval);
        acc.prefetch(s.bot.interval);
        return true;
    };
    // Finish s's search: the count and top's text position.
    auto finish = [&](Slot& s) {
        TmsSearchItem& it = items[s.item];
        LFPos top = s.top, bot = s.bot;
        acc.resolve_LF(top);
        acc.resolve_LF(bot);
        const TmsRowStarts& rs = ctr.rows();
        ulint diff;
        if (bot.interval - top.interval <= (ulint(1) << rs.shift())) {
            diff = 0;
            for (ulint i = top.interval; i < bot.interval; ++i) diff += acc.length(acc.row(i));
            st.sum_rows += bot.interval - top.interval;
        } else {
            diff = rs.start(acc, bot.interval) - rs.start(acc, top.interval);
            ++st.far_ends;
        }
        it.count = diff + bot.offset - top.offset + 1;
        // The F block's first row is one LF step from its run head, as is
        // the row after a move from the run head the move found.
        const ulint x = s.head > last ? ctr.start(code[static_cast<uchar>(it.pat[it.len - 1])]).top_x + 1
                                      : TmsCounter::head_text_pos(idx, s.head);
        it.top_x = x - s.since;
    };
    // One extension of s, or its finish.  Returns false once its search and
    // every search after it are done.
    auto visit = [&](Slot& s) {
        if (s.left == 0) {
            finish(s);
            return begin(s);
        }
        const TmsSearchItem& it = items[s.item];
        const uchar c = code[static_cast<uchar>(it.pat[s.left - 1])];
        ++st.search_steps;
        LFPos top = s.top, bot = s.bot;
        const auto rt = acc.resolve_LF(top);
        const auto rb = acc.resolve_LF(bot);
        st.lf_ff += (top.interval - s.top.interval) + (bot.interval - s.bot.interval);
        if (acc.code(rt) == c) {
            s.top = acc.start_LF(rt, top);
        } else {
            ulint d = top.interval;
            do {
                if (d == last) throw std::logic_error("backward search left the BWT");
                ++d;
            } while (acc.code(acc.row(d)) != c);
            ++st.top_moves;
            st.scan_rows += d - top.interval;
            LFPos h{};
            h.interval = d;
            s.top = acc.start_LF(acc.row(d), h);
            s.head = d;
            s.since = 0;
        }
        if (acc.code(rb) == c) {
            s.bot = acc.start_LF(rb, bot);
        } else {
            ulint u = bot.interval;
            do {
                if (u == 0) throw std::logic_error("backward search left the BWT");
                --u;
            } while (acc.code(acc.row(u)) != c);
            ++st.bot_moves;
            st.scan_rows += bot.interval - u;
            LFPos t{};
            t.interval = u;
            const auto ru = acc.row(u);
            t.offset = acc.length(ru) - 1;
            s.bot = acc.start_LF(ru, t);
        }
        ++s.since;
        --s.left;
        acc.prefetch(s.top.interval);
        acc.prefetch(s.bot.interval);
        return true;
    };
    std::vector<Slot> slots(std::min(k == 0 ? size_t(1) : k, items.size()));
    size_t live = 0;
    for (auto& s : slots)
        if (begin(s)) ++live;
        else break;
    while (live > 0)
        for (size_t s = 0; s < live;) {
            if (visit(slots[s])) ++s;
            else std::swap(slots[s], slots[--live]);
        }
}

/** A listing of the n positions in row order from text position x, the first of them. */
struct TmsListItem {
    ulint x, n;
    ulint* out;  // receives the n positions
};

/**
 * Fill each item's positions by n - 1 phi_inv steps from x, up to k
 * listings at a time.  Needs an index with phi_inv.
 */
inline void tms_list_batch(const TmsIndex& idx, std::vector<TmsListItem>& items, size_t k) {
    using Pos = TmsIndex::PhiInvPos;
    const auto wk = idx.phi_inv_walker();
    TmsCountStats& st = tms_count_stats;
    struct Slot {
        Pos p;
        ulint at;
        size_t item;
    };
    size_t next = 0;
    auto begin = [&](Slot& s) {
        while (next < items.size()) {
            TmsListItem& it = items[next++];
            it.out[0] = it.x;
            if (it.n <= 1) continue;
            s.item = static_cast<size_t>(&it - items.data());
            s.p = wk.start_step(wk.locate(it.x));
            wk.prefetch(s.p.interval);
            s.at = 1;
            return true;
        }
        return false;
    };
    auto visit = [&](Slot& s) {
        TmsListItem& it = items[s.item];
        ulint lcp;
        it.out[s.at++] = wk.step(s.p, lcp);
        ++st.list_steps;
        if (s.at < it.n) {
            wk.prefetch(s.p.interval);
            return true;
        }
        return begin(s);
    };
    std::vector<Slot> slots(std::min(k == 0 ? size_t(1) : k, std::max<size_t>(items.size(), 1)));
    size_t live = 0;
    for (auto& s : slots)
        if (begin(s)) ++live;
        else break;
    while (live > 0)
        for (size_t s = 0; s < live;) {
            if (visit(slots[s])) ++s;
            else std::swap(slots[s], slots[--live]);
        }
}

#endif /* _TMS_COUNT_HPP */
