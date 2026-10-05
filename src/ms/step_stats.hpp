/**
 * Per-step counts for builds with MS_STEP_STATS (make ms-steps): where in a
 * read the repositions, the rows they read and the query time fall.  They
 * measure how much of a query an ftab, a table that jumps over a read's
 * first k steps, would skip.  Normal builds leave this header empty.
 *
 * A step is one character, counted from where the walk starts, the read's
 * right end.  The query records each step's repositions, rows read and
 * cycles (from the time stamp counter) as it goes, and flush() files the
 * read's steps by their index and by their distance from the last character
 * other than A, C, G and T (such as N), and works out which of the read's
 * windows an ftab over ACGT could jump: the last j bases, and the j bases
 * after each run of other characters, can be jumped when they occur in the
 * text, which the read's matching statistics say.  An index whose text holds
 * N (an FMD's, for one) matches a read's N rather than restarting there, but
 * an ftab over ACGT has no entry for a window holding one either way.  Reads are recorded only when one is in flight
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
    // By steps since the end of the last run of characters other than ACGT,
    // for steps after one.
    Arr n_steps{}, n_reps{}, n_rows{}, n_cyc{};
    // For j = 1 to J: reads with at least j bases, those whose last j bases
    // occur, and the repositions, rows and cycles of the first j steps of
    // those reads.
    Arr start_fit{}, start_jump{}, start_jump_reps{}, start_jump_rows{}, start_jump_cyc{};
    // The same for the j bases after each run of characters other than ACGT.
    Arr postn_fit{}, postn_jump{}, postn_jump_reps{}, postn_jump_rows{}, postn_jump_cyc{};
    step_ulint reads = 0, bases = 0, reads_other = 0, other_chars = 0, other_runs = 0;
    step_ulint tot_reps = 0, tot_rows = 0, tot_cyc = 0;

    StepAgg& operator+=(const StepAgg& o) {
        auto add = [](Arr& a, const Arr& b) { for (size_t i = 0; i <= J; ++i) a[i] += b[i]; };
        add(steps, o.steps); add(reps, o.reps); add(rows, o.rows); add(cyc, o.cyc);
        add(n_steps, o.n_steps); add(n_reps, o.n_reps); add(n_rows, o.n_rows); add(n_cyc, o.n_cyc);
        add(start_fit, o.start_fit); add(start_jump, o.start_jump); add(start_jump_reps, o.start_jump_reps);
        add(start_jump_rows, o.start_jump_rows); add(start_jump_cyc, o.start_jump_cyc);
        add(postn_fit, o.postn_fit); add(postn_jump, o.postn_jump); add(postn_jump_reps, o.postn_jump_reps);
        add(postn_jump_rows, o.postn_jump_rows); add(postn_jump_cyc, o.postn_jump_cyc);
        reads += o.reads; bases += o.bases; reads_other += o.reads_other;
        other_chars += o.other_chars; other_runs += o.other_runs;
        tot_reps += o.tot_reps; tot_rows += o.tot_rows; tot_cyc += o.tot_cyc;
        return *this;
    }

    /** One line of totals, then one line per j of the counts above. */
    void print(std::ostream& os) const {
        os << "stepstats: reads=" << reads << " bases=" << bases << " reads_other=" << reads_other
           << " other_chars=" << other_chars << " other_runs=" << other_runs << " reps=" << tot_reps
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
    const char* pat = nullptr;       // the read
    size_t m = 0;
    bool pending = false;

    /** Start recording read p of length len whose statistics go to out, filing the previous one. */
    void begin(const char* p, size_t len, const step_ulint* out) {
        flush();
        pat = p;
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
        // Step t reads pattern position m - 1 - t.
        auto other = [&](size_t t) {
            const char c = pat[m - 1 - t];
            return c != 'A' && c != 'C' && c != 'G' && c != 'T';
        };
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
                // left end covers it and it holds only ACGT.
                bool acgt = true;
                for (size_t t = s; t < s + j && acgt; ++t) acgt = !other(t);
                if (acgt && ms[m - s - j] >= j) {
                    ++jump[j];
                    jr[j] += pr[s + j] - pr[s];
                    jw[j] += pw[s + j] - pw[s];
                    jc[j] += pc[s + j] - pc[s];
                }
            }
        };
        window(0, a.start_fit, a.start_jump, a.start_jump_reps, a.start_jump_rows, a.start_jump_cyc);
        bool any_other = false;
        size_t seg = SIZE_MAX;  // first step after the last run of characters other than ACGT
        for (size_t t = 0; t < m; ++t) {
            const size_t b = std::min(t, J);
            a.steps[b] += 1;
            a.reps[b] += reps[t];
            a.rows[b] += rows[t];
            a.cyc[b] += cyc[t];
            if (other(t)) {
                ++a.other_chars;
                any_other = true;
                if (t + 1 == m || !other(t + 1)) {
                    ++a.other_runs;
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
        if (any_other) ++a.reads_other;
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

/**
 * Fast-forward lengths: for each move structure a query steps through, a
 * histogram of how many intervals each fast-forward crosses after a move,
 * that is, how far the landing interval lies past the one the pointer names.
 * A move resolved with exponential search (starts-based phi) reads about
 * 2 log2(k) + 2 interval starts to cross k intervals and a linear one k + 1,
 * so each histogram also sums the starts an exponential search would read
 * (exp_probes), which is the real count for starts-based phi and what it
 * would be for the others.  Counts are kept for every move, whatever
 * --interleave is.
 */
enum FFKind : int {
    FF_LF,          // LF steps
    FF_PSI,         // psi steps, over the same rows as LF in the psi and full layouts
    FF_PHI_LEN,     // phi steps of matching statistics, lengths-based (linear)
    FF_PHI_STARTS,  // phi steps of matching statistics, starts-based (exponential)
    FF_PHI_TOEHOLD, // the phi step from a pred candidate's next row head to its tail, either form
    FF_PHI_WALK,    // phi steps of SMEM position walks, starts-based (linear)
    FF_PHI_INV,     // phi_inv steps of SMEM position walks, starts-based (linear)
    FF_KINDS
};

inline const char* ff_kind_name(int k) {
    static const char* names[] = {"lf", "psi", "phi_len", "phi_starts", "phi_toehold", "phi_walk", "phi_inv"};
    return names[k];
}

/**
 * The interval starts fast_forward_exponential (move_structure_impl.hpp)
 * reads to cross k intervals, away from the table's end: one to test the
 * next start, then the doubling probes, then the binary search.
 */
inline step_ulint exp_probes(step_ulint k) {
    if (k == 0) return 1;
    step_ulint probes = 1, hi = 1;
    while (hi <= k) { ++probes; hi *= 2; }
    ++probes;  // the doubling probe that fails
    step_ulint lo = hi / 2;
    while (lo + 1 < hi) {
        ++probes;
        const step_ulint mid = lo + (hi - lo) / 2;
        if (k >= mid) lo = mid;
        else hi = mid;
    }
    return probes;
}

struct FFHist {
    // Crossings 0 to B - 1 are counted apart, longer ones in bucket B.
    static constexpr size_t B = 1024;
    std::array<step_ulint, B + 1> by_k{};
    std::array<step_ulint, 64> by_probes{};
    step_ulint moves = 0, crossed = 0, probes = 0, max = 0;

    void add(step_ulint k) {
        ++moves;
        crossed += k;
        max = std::max(max, k);
        ++by_k[std::min<step_ulint>(k, B)];
        const step_ulint p = exp_probes(k);
        probes += p;
        ++by_probes[std::min<step_ulint>(p, 63)];
    }
    FFHist& operator+=(const FFHist& o) {
        for (size_t i = 0; i <= B; ++i) by_k[i] += o.by_k[i];
        for (size_t i = 0; i < 64; ++i) by_probes[i] += o.by_probes[i];
        moves += o.moves; crossed += o.crossed; probes += o.probes; max = std::max(max, o.max);
        return *this;
    }
    /** The smallest k with at least fraction q of the moves crossing k or fewer. */
    step_ulint quantile(double q) const {
        const double want = q * double(moves);
        step_ulint acc = 0;
        for (size_t i = 0; i <= B; ++i) {
            acc += by_k[i];
            if (double(acc) >= want) return i;
        }
        return B;
    }
    step_ulint range(size_t lo, size_t hi) const {
        step_ulint s = 0;
        for (size_t i = lo; i <= std::min(hi, B); ++i) s += by_k[i];
        return s;
    }
};

struct FFAgg {
    std::array<FFHist, FF_KINDS> h;
    FFAgg& operator+=(const FFAgg& o) {
        for (int k = 0; k < FF_KINDS; ++k) h[k] += o.h[k];
        return *this;
    }
    /**
     * One line per structure with moves: the mean, p50, p99 and largest
     * crossing count, the shares crossing 0, 1, 2-3, 4-7 and 8 or more, the
     * mean starts read by linear (k + 1) and exponential search, and then
     * the histogram of crossings as k:count pairs.
     */
    void print(std::ostream& os) const {
        for (int k = 0; k < FF_KINDS; ++k) {
            const FFHist& x = h[k];
            if (x.moves == 0) continue;
            const double m = double(x.moves);
            os << "ffstats: kind=" << ff_kind_name(k) << " moves=" << x.moves << " mean=" << double(x.crossed) / m
               << " p50=" << x.quantile(0.5) << " p99=" << x.quantile(0.99) << " p999=" << x.quantile(0.999)
               << " max=" << x.max << " share0=" << double(x.range(0, 0)) / m
               << " share1=" << double(x.range(1, 1)) / m << " share2_3=" << double(x.range(2, 3)) / m
               << " share4_7=" << double(x.range(4, 7)) / m << " share8plus=" << double(x.range(8, FFHist::B)) / m
               << " linear_reads=" << double(x.crossed + x.moves) / m << " exp_reads=" << double(x.probes) / m
               << "\n";
            os << "ffhist: kind=" << ff_kind_name(k);
            for (size_t i = 0; i <= FFHist::B; ++i)
                if (x.by_k[i]) os << ' ' << i << (i == FFHist::B ? "+" : "") << ':' << x.by_k[i];
            os << "\n";
        }
    }
};

inline thread_local FFAgg ff_agg;
#define MS_FF(kind, k) (ff_agg.h[(kind)].add(static_cast<step_ulint>(k)))

#endif /* MS_STEP_STATS */

#ifndef MS_FF
#define MS_FF(kind, k) ((void)0)
#endif

#endif /* _STEP_STATS_HPP */
