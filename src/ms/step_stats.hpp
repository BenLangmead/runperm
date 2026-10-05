/**
 * Per-step counts for builds with MS_STEP_STATS (make ms-steps): where in a
 * read the repositions, the rows they read and the query time fall.  They
 * measure how much of a query an ftab, a table that jumps over a read's
 * first k steps, would skip.  Normal builds leave this header empty.
 *
 * A step is one character, counted from where the walk starts, the read's
 * right end.  The query records each step's repositions, rows read and
 * cycles (from the time stamp counter) as it goes, and finish() files the
 * read's steps by their index and by their distance from the last absent
 * character (such as N), and works out which of the read's windows an ftab
 * could jump: the last j bases, and the j bases after each run of absent
 * characters, can be jumped when they occur in the text, which the read's
 * matching statistics say.  Reads are recorded only when one is in flight
 * at a time (ms batch and tms-batch with --interleave 0).
 */

#ifndef _STEP_STATS_HPP
#define _STEP_STATS_HPP

#ifdef MS_STEP_STATS

#include "orbit/common.hpp"
#include <algorithm>
#include <array>
#include <cstdint>
#include <ostream>
#include <vector>
#if defined(__x86_64__) || defined(__i386__)
#include <x86intrin.h>
#endif

using step_ulint = orbit::ulint;

inline uint64_t step_ticks() {
#if defined(__x86_64__) || defined(__i386__)
    return __rdtsc();
#else
    return 0;
#endif
}

/** Counts summed over reads. */
struct StepAgg {
    // Window lengths and step indexes 0 to J - 1 are kept apart; later steps
    // share bucket J.
    static constexpr size_t J = 32;
    using Arr = std::array<step_ulint, J + 1>;
    // By step index from the read's right end.
    Arr steps{}, reps{}, rows{}, cyc{};
    // By steps since the end of the last run of absent characters, for steps
    // after one.
    Arr n_steps{}, n_reps{}, n_rows{}, n_cyc{};
    // For j = 1 to J: reads with at least j bases, those whose last j bases
    // occur, and the repositions, rows and cycles of the first j steps of
    // those reads.
    Arr start_fit{}, start_jump{}, start_jump_reps{}, start_jump_rows{}, start_jump_cyc{};
    // The same for the j bases after each run of absent characters.
    Arr postn_fit{}, postn_jump{}, postn_jump_reps{}, postn_jump_rows{}, postn_jump_cyc{};
    step_ulint reads = 0, bases = 0, reads_absent = 0, absent_chars = 0, absent_runs = 0;
    step_ulint tot_reps = 0, tot_rows = 0, tot_cyc = 0;

    StepAgg& operator+=(const StepAgg& o) {
        auto add = [](Arr& a, const Arr& b) { for (size_t i = 0; i <= J; ++i) a[i] += b[i]; };
        add(steps, o.steps); add(reps, o.reps); add(rows, o.rows); add(cyc, o.cyc);
        add(n_steps, o.n_steps); add(n_reps, o.n_reps); add(n_rows, o.n_rows); add(n_cyc, o.n_cyc);
        add(start_fit, o.start_fit); add(start_jump, o.start_jump); add(start_jump_reps, o.start_jump_reps);
        add(start_jump_rows, o.start_jump_rows); add(start_jump_cyc, o.start_jump_cyc);
        add(postn_fit, o.postn_fit); add(postn_jump, o.postn_jump); add(postn_jump_reps, o.postn_jump_reps);
        add(postn_jump_rows, o.postn_jump_rows); add(postn_jump_cyc, o.postn_jump_cyc);
        reads += o.reads; bases += o.bases; reads_absent += o.reads_absent;
        absent_chars += o.absent_chars; absent_runs += o.absent_runs;
        tot_reps += o.tot_reps; tot_rows += o.tot_rows; tot_cyc += o.tot_cyc;
        return *this;
    }

    /** One line of totals, then one line per j of the counts above. */
    void print(std::ostream& os) const {
        os << "stepstats: reads=" << reads << " bases=" << bases << " reads_absent=" << reads_absent
           << " absent_chars=" << absent_chars << " absent_runs=" << absent_runs << " reps=" << tot_reps
           << " rows=" << tot_rows << " cycles=" << tot_cyc << "\n";
        os << "stepstats-cols: j steps reps rows cyc n_steps n_reps n_rows n_cyc start_fit start_jump start_jump_reps "
              "start_jump_rows start_jump_cyc postn_fit postn_jump postn_jump_reps postn_jump_rows postn_jump_cyc\n";
        for (size_t j = 0; j <= J; ++j) {
            os << "stepstats-j " << j << ' ' << steps[j] << ' ' << reps[j] << ' ' << rows[j] << ' ' << cyc[j] << ' '
               << n_steps[j] << ' ' << n_reps[j] << ' ' << n_rows[j] << ' ' << n_cyc[j] << ' ' << start_fit[j] << ' '
               << start_jump[j] << ' ' << start_jump_reps[j] << ' ' << start_jump_rows[j] << ' ' << start_jump_cyc[j]
               << ' ' << postn_fit[j] << ' ' << postn_jump[j] << ' ' << postn_jump_reps[j] << ' '
               << postn_jump_rows[j] << ' ' << postn_jump_cyc[j] << "\n";
        }
    }
};

inline thread_local StepAgg step_agg;

/** The steps of the read in flight. */
struct StepRecorder {
    std::vector<step_ulint> reps, rows, cyc;
    const step_ulint* ms = nullptr;  // the read's statistics, complete when it is finished
    size_t m = 0;
    bool pending = false;

    /** Start recording a read of length m whose statistics go to ms, filing the previous one. */
    void begin(size_t len, const step_ulint* out) {
        flush();
        m = len;
        ms = out;
        reps.assign(m, 0);
        rows.assign(m, 0);
        cyc.assign(m, 0);
        pending = true;
    }
    void add(size_t t, step_ulint r, step_ulint w, step_ulint c) {
        if (!pending || t >= m) return;
        reps[t] += r;
        rows[t] += w;
        cyc[t] += c;
    }
    /** File the read in flight, if any, into step_agg. */
    void flush() {
        if (!pending) return;
        pending = false;
        StepAgg& a = step_agg;
        constexpr size_t J = StepAgg::J;
        ++a.reads;
        a.bases += m;
        // Step t reads pattern position m - 1 - t, whose statistic is 0
        // exactly when that character does not occur in the text.
        auto absent = [&](size_t t) { return ms[m - 1 - t] == 0; };
        // Prefix sums over steps, for the windows.
        std::vector<step_ulint> pr(m + 1, 0), pw(m + 1, 0), pc(m + 1, 0);
        for (size_t t = 0; t < m; ++t) {
            pr[t + 1] = pr[t] + reps[t];
            pw[t + 1] = pw[t] + rows[t];
            pc[t + 1] = pc[t] + cyc[t];
        }
        a.tot_reps += pr[m];
        a.tot_rows += pw[m];
        a.tot_cyc += pc[m];
        auto window = [&](size_t s, StepAgg::Arr& fit, StepAgg::Arr& jump,
                          StepAgg::Arr& jr, StepAgg::Arr& jw, StepAgg::Arr& jc) {
            for (size_t j = 1; j <= J && s + j <= m; ++j) {
                ++fit[j];
                // The window is steps s to s + j - 1, pattern positions
                // m - s - j to m - s - 1; it occurs when the statistic at its
                // left end covers it.
                if (ms[m - s - j] >= j) {
                    ++jump[j];
                    jr[j] += pr[s + j] - pr[s];
                    jw[j] += pw[s + j] - pw[s];
                    jc[j] += pc[s + j] - pc[s];
                }
            }
        };
        window(0, a.start_fit, a.start_jump, a.start_jump_reps, a.start_jump_rows, a.start_jump_cyc);
        bool any_absent = false;
        size_t seg = SIZE_MAX;  // first step after the last run of absent characters
        for (size_t t = 0; t < m; ++t) {
            const size_t b = std::min(t, J);
            a.steps[b] += 1;
            a.reps[b] += reps[t];
            a.rows[b] += rows[t];
            a.cyc[b] += cyc[t];
            if (absent(t)) {
                ++a.absent_chars;
                any_absent = true;
                if (t + 1 == m || !absent(t + 1)) {
                    ++a.absent_runs;
                    seg = t + 1;
                    if (seg < m)
                        window(seg, a.postn_fit, a.postn_jump, a.postn_jump_reps, a.postn_jump_rows,
                               a.postn_jump_cyc);
                }
                continue;
            }
            if (seg != SIZE_MAX) {
                const size_t u = std::min(t - seg, J);
                a.n_steps[u] += 1;
                a.n_reps[u] += reps[t];
                a.n_rows[u] += rows[t];
                a.n_cyc[u] += cyc[t];
            }
        }
        if (any_absent) ++a.reads_absent;
    }
};

inline thread_local StepRecorder step_rec;

/**
 * Measures one visit to step t: the repositions and rows counted in
 * *reps and *rows_a + *rows_b, and the cycles, between its construction and
 * record(), which the destructor calls if nothing else has.
 */
struct StepProbe {
    size_t t;
    const step_ulint *reps, *rows_a, *rows_b;
    step_ulint r0, w0;
    uint64_t c0;
    bool done = false;
    StepProbe(size_t step, const step_ulint* r, const step_ulint* wa, const step_ulint* wb)
        : t(step), reps(r), rows_a(wa), rows_b(wb), r0(*r), w0(*wa + *wb), c0(step_ticks()) {}
    void record() {
        if (done) return;
        done = true;
        const uint64_t c1 = step_ticks();
        step_rec.add(t, *reps - r0, *rows_a + *rows_b - w0, c1 - c0);
    }
    ~StepProbe() { record(); }
};

#endif /* MS_STEP_STATS */

#endif /* _STEP_STATS_HPP */
