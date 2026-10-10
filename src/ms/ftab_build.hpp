/**
 * Building an ftab (ftab.hpp) for an ms index or a tms index.
 *
 * As in ftab_survey.hpp, the walk over a k-mer reads it right to left from
 * idx.first() with match length 0, exactly as the index's query does, and
 * k-mers that share a suffix share the walk over it.  The builder walks the
 * trie of suffixes depth first, records the state at every occurring k-mer,
 * and counts the occurring j-mers at every depth j up to k.
 *
 * A tms reposition compares the pattern to the right of the current
 * character, so each trie node carries its suffix.  On an index with phi,
 * each node also carries the toehold the batched engine carries: the phi
 * point of the current row's text position, moved left with each LF step
 * and set from the chosen candidate at each reposition.  The LCE of a
 * reposition is read with psi where the index has it, and otherwise (the
 * phi layout) from the PLCP values phi walks over between the rows; both
 * choose the candidate the query chooses.  A suffix occurs only if its
 * capped LCE reaches the cap, so a walk stops as soon as it falls short.
 */

#ifndef _FTAB_BUILD_HPP
#define _FTAB_BUILD_HPP

#include "ftab.hpp"
#include "ftab_survey.hpp"
#include "ms_rlbwt.hpp"
#include "tms_index.hpp"
#include "tms_query.hpp"
#include <algorithm>
#include <atomic>
#include <stdexcept>
#include <thread>
#include <vector>

namespace ftab {

/** A table, and how many j-mers occur for each j from 0 to its k, as its walk counted them. */
struct Built {
    Ftab tab;
    std::vector<uint64_t> present;
};

/**
 * Walk the suffix trie to depth k on the given number of threads.
 * step(state, depth, c, rest) reads c to the left of the occurring suffix
 * rest (of length depth) whose state is state, and returns true, with the
 * new state, if c + rest occurs.  record(state, depth, key) is called, from
 * any thread, for every occurring suffix of depth 1 to k, whose key holds
 * the character d places from its right end in bits 2d and 2d + 1.  Returns
 * the number of occurring suffixes at each depth from 0 to k.
 */
template <class State, class Step, class Record>
inline std::vector<uint64_t> walk(unsigned k, const State& root, Step step, Record record, size_t threads) {
    if (k < 1 || k > max_k) throw std::invalid_argument("ftab k must be between 1 and 31");
    static constexpr char acgt[4] = {'A', 'C', 'G', 'T'};
    struct Node {
        State st;
        uint64_t key;
        unsigned depth;
        char buf[max_k];  // the suffix, in buf[max_k - depth .. max_k)
    };
    // Children of an occurring node go to out if their depth is below k.
    auto expand = [&](const Node& nd, std::vector<Node>& out, std::vector<uint64_t>& present) {
        for (unsigned c = 0; c < 4; ++c) {
            Node ch = nd;
            if (!step(ch.st, nd.depth, static_cast<uchar>(acgt[c]), nd.buf + max_k - nd.depth)) continue;
            ch.depth = nd.depth + 1;
            ch.key = nd.key | (uint64_t{c} << (2 * nd.depth));
            ch.buf[max_k - ch.depth] = acgt[c];
            ++present[ch.depth];
            record(ch.st, ch.depth, ch.key);
            if (ch.depth < k) out.push_back(ch);
        }
    };
    if (threads == 0) threads = 1;
    std::vector<std::vector<uint64_t>> counts(threads + 1, std::vector<uint64_t>(k + 1, 0));
    Node root_node{};
    root_node.st = root;
    std::vector<Node> frontier{root_node};
    const unsigned seed_depth = std::min<unsigned>(k - 1, 6);
    for (unsigned d = 0; d < seed_depth; ++d) {
        std::vector<Node> next;
        for (const auto& nd : frontier) expand(nd, next, counts[threads]);
        frontier.swap(next);
    }
    std::atomic<size_t> next_task{0};
    auto work = [&](size_t t) {
        std::vector<Node> stack;
        for (size_t j; (j = next_task.fetch_add(1)) < frontier.size();) {
            stack.push_back(frontier[j]);
            while (!stack.empty()) {
                const Node nd = stack.back();
                stack.pop_back();
                expand(nd, stack, counts[t]);
            }
        }
    };
    std::vector<std::thread> pool;
    for (size_t t = 1; t < threads; ++t) pool.emplace_back(work, t);
    work(0);
    for (auto& th : pool) th.join();
    std::vector<uint64_t> present(k + 1, 0);
    present[0] = 1;
    for (const auto& c : counts)
        for (unsigned j = 1; j <= k; ++j) present[j] += c[j];
    return present;
}

/** The largest of f(i) for i in [0, n), on the given number of threads. */
template <class F>
inline uint64_t max_over(uint64_t n, F f, size_t threads) {
    if (threads == 0) threads = 1;
    std::vector<uint64_t> best(threads, 0);
    auto work = [&](size_t t) {
        const uint64_t lo = n / threads * t, hi = t + 1 == threads ? n : n / threads * (t + 1);
        uint64_t m = 0;
        for (uint64_t i = lo; i < hi; ++i) m = std::max<uint64_t>(m, f(i));
        best[t] = m;
    };
    std::vector<std::thread> pool;
    for (size_t t = 1; t < threads; ++t) pool.emplace_back(work, t);
    work(0);
    for (auto& th : pool) th.join();
    return *std::max_element(best.begin(), best.end());
}

/** Bits for an interval of a structure with n intervals, and for an offset below max_len. */
inline unsigned interval_bits(uint64_t n) { return bits_for(n > 0 ? n - 1 : 0); }
inline unsigned offset_bits(uint64_t max_len) { return bits_for(max_len > 0 ? max_len - 1 : 0); }

// --- ms ------------------------------------------------------------------

/** The shape of an ms index's tables, with k = 0. */
template <bool SP>
inline Shape ms_shape(MSIndexSpillLCP<SP>& idx, size_t threads) {
    Shape s;
    s.interval = interval_bits(idx.move_runs());
    s.offset = offset_bits(max_over(idx.move_runs(), [&](uint64_t i) { return idx.get_length(i); }, threads));
    return s;
}

/** The ftab of an ms index: states as ms_query and ms_query_batch reach them. */
template <bool SP>
inline Built build_ms(MSIndexSpillLCP<SP>& idx, const Shape& shape, size_t threads) {
    using Position = typename MSIndexSpillLCP<SP>::position;
    Built b{Ftab(shape), {}};
    auto step = [&idx](Position& p, unsigned depth, uchar c, const char*) {
        ulint len = depth;
        return ftab_survey::step(idx, p, len, c) && len == depth + 1;
    };
    auto record = [&b, k = shape.k](const Position& p, unsigned depth, uint64_t key) {
        if (depth == k) b.tab.set(key, Entry{p.interval, p.offset, 0, 0});
    };
    b.present = walk(shape.k, idx.first(), step, record, threads);
    return b;
}

template <bool SP>
inline Built build_ms(MSIndexSpillLCP<SP>& idx, unsigned k, size_t threads) {
    Shape s = ms_shape(idx, threads);
    s.k = k;
    return build_ms(idx, s, threads);
}

// --- tms -----------------------------------------------------------------

/**
 * The shape of a tms index's tables, with k = 0: with a toehold if the
 * index has phi, and with window positions too if it has phi's
 * starts-based form, which gives a phi point's text position.
 */
inline Shape tms_shape(TmsIndex& idx, size_t threads) {
    Shape s;
    s.interval = interval_bits(idx.move_runs());
    s.offset = offset_bits(max_over(idx.move_runs(), [&](uint64_t i) { return idx.get_length(i); }, threads));
    if (idx.has_phi()) {
        s.toehold = true;
        s.toehold_interval = interval_bits(idx.phi_intervals());
        s.toehold_offset = offset_bits(max_over(idx.phi_intervals(), [&](uint64_t i) { return idx.phi_length(i); }, threads));
        s.window = idx.has_phi_starts();
        s.position = s.window ? bits_for(idx.domain() - 1) : 0;
    }
    return s;
}

namespace detail {
template <bool Toehold, bool Starts>
inline Built build_tms(TmsIndex& idx, const Shape& shape, size_t threads) {
    using LFPos = TmsIndex::LFPos;
    using PhiPos = TmsIndex::PhiPosOf<Starts>;
    struct State {
        LFPos pos;
        PhiPos ph;  // the toehold, with Toehold
    };
    const bool psi = idx.has_psi();
    if (!psi && !Toehold) throw std::invalid_argument("a tms ftab needs an index with psi or phi");
    Built b{Ftab(shape), {}};
    const ulint last_run = idx.move_runs() - 1;
    // Whether the dist PLCP values from phi point p on are all at least cap.
    auto phi_reaches = [&idx](PhiPos p, ulint dist, ulint cap) {
        for (ulint t = 0; t < dist; ++t) {
            if (idx.plcp(p) < cap) return false;
            if (t + 1 < dist) p = idx.finish_phi(idx.start_phi(p));
        }
        return true;
    };
    auto step = [&](State& s, unsigned depth, uchar c, const char* rest) {
        if (idx.get_character(s.pos.interval) != c) {
            if (!idx.occurs(c)) return false;
            const ulint cur = s.pos.interval;
            ulint u = cur, d = cur;
            ulint dist_up = s.pos.offset + 1, dist_down = idx.get_length(cur) - s.pos.offset;
            bool found_up = false, found_down = false;
            while (u > 0) {
                if (idx.get_character(--u) == c) { found_up = true; break; }
                dist_up += idx.get_length(u);
            }
            while (d < last_run) {
                if (idx.get_character(++d) == c) { found_down = true; break; }
                dist_down += idx.get_length(d);
            }
            // The suffix occurs if a candidate's capped LCE reaches the
            // cap, depth; ties go to pred, as in the query.
            auto reaches_up = [&] {
                if (depth == 0) return true;
                if (psi) return tms_psi_lce(idx, idx.psi_at_tail(u), rest, depth, depth) == depth;
                if constexpr (Toehold) return phi_reaches(s.ph, dist_up, depth);
                return false;
            };
            auto reaches_down = [&] {
                if (depth == 0) return true;
                if (psi) return tms_psi_lce(idx, idx.psi_at_head(d), rest, depth, depth) == depth;
                if constexpr (Toehold) return phi_reaches(idx.resolve_phi(idx.template phi_at_head<Starts>(d)), dist_down, depth);
                return false;
            };
            LFPos from{};
            if (found_up && reaches_up()) {
                from.interval = u + 1;
                from = idx.up(from);
                if constexpr (Toehold) s.ph = idx.finish_phi(idx.start_phi(idx.template phi_at_head<Starts>(u + 1)));
            } else if (found_down && reaches_down()) {
                from.interval = d - 1;
                from = idx.down(from);
                if constexpr (Toehold) s.ph = idx.resolve_phi(idx.template phi_at_head<Starts>(d));
            } else {
                return false;
            }
            s.pos = from;
        }
        s.pos = idx.LF_step(s.pos);
        if constexpr (Toehold) s.ph = idx.phi_left(s.ph);
        return true;
    };
    const unsigned k = shape.k;
    const bool window = shape.window;
    auto record = [&b, k, window](const State& s, unsigned depth, uint64_t key) {
        if (depth == k) {
            Entry e{s.pos.interval, s.pos.offset, 0, 0};
            if constexpr (Toehold) {
                e.toehold_interval = s.ph.interval;
                e.toehold_offset = s.ph.offset;
            }
            b.tab.set(key, e);
        } else if constexpr (Starts) {
            if (window) b.tab.set_window(depth, key, s.ph.idx);
        }
    };
    State root{};
    root.pos = idx.first();
    if constexpr (Toehold) root.ph = idx.resolve_phi(idx.template phi_at_head<Starts>(0));
    b.present = walk(k, root, step, record, threads);
    return b;
}
}  // namespace detail

/**
 * The ftab of a tms index: states as tms_query_batch reaches them in every
 * mode.  The shape (from tms_shape) says whether it holds a toehold and
 * window positions; those need the index loaded with phi, and window
 * positions phi's starts-based form.
 */
inline Built build_tms(TmsIndex& idx, const Shape& shape, size_t threads) {
    if (shape.toehold && !idx.has_phi()) throw std::invalid_argument("a tms ftab with a toehold needs phi");
    if (shape.window && !idx.has_phi_starts())
        throw std::invalid_argument("a tms ftab with window positions needs phi's starts-based form");
    if (!shape.toehold) return detail::build_tms<false, false>(idx, shape, threads);
    if (idx.has_phi_starts()) return detail::build_tms<true, true>(idx, shape, threads);
    return detail::build_tms<true, false>(idx, shape, threads);
}

inline Built build_tms(TmsIndex& idx, unsigned k, size_t threads) {
    Shape s = tms_shape(idx, threads);
    s.k = k;
    return build_tms(idx, s, threads);
}

// --- Choosing k ----------------------------------------------------------

/**
 * The rule `ftab-build --auto` follows: the largest k whose table fits in
 * budget of the index file's bytes; then, while fewer than min_occupancy of
 * that k's k-mers occur, one smaller (but at least 1).
 */
struct AutoRule {
    double budget = 0.05;
    double min_occupancy = 0.8;
};

/** The largest k (at least 1) whose table, of shape s otherwise, takes at most budget_bytes. */
inline unsigned largest_fitting_k(Shape s, double budget_bytes) {
    unsigned best = 1;
    for (unsigned k = 1; k <= 24; ++k) {
        s.k = k;
        if (double(s.payload_bytes()) <= budget_bytes) best = k;
        else break;
    }
    return best;
}

/** Fraction of the 4^j j-mers that occur, given present from a walk. */
inline double occupancy(const std::vector<uint64_t>& present, unsigned j) {
    return double(present[j]) / double(uint64_t{1} << (2 * j));
}

/** From k0, step down while fewer than min_occupancy of the k-mers occur, stopping at 1. */
inline unsigned step_down(unsigned k0, const std::vector<uint64_t>& present, double min_occupancy) {
    unsigned k = k0;
    while (k > 1 && occupancy(present, k) < min_occupancy) --k;
    return k;
}

/** What the auto rule chose, and the counts of the walk it chose from. */
struct AutoChoice {
    unsigned k_fit = 0;             // the largest k that fits
    unsigned k = 0;                 // the k chosen
    std::vector<uint64_t> present;  // occurring j-mers for j up to k_fit
};

/**
 * Build with the auto rule for an index of index_bytes bytes, whose tables
 * have shape base (apart from k), with build(shape) returning a Built.  The
 * table at the largest fitting k is built first, counting every depth; if
 * the rule then steps down, the table at the chosen k is built in a second
 * walk.
 */
template <class Build>
inline Built build_auto(const Shape& base, uint64_t index_bytes, const AutoRule& rule, Build build, AutoChoice& choice) {
    Shape s = base;
    s.k = largest_fitting_k(base, rule.budget * double(index_bytes));
    Built b = build(s);
    choice.k_fit = s.k;
    choice.present = b.present;
    choice.k = step_down(s.k, b.present, rule.min_occupancy);
    if (choice.k == s.k) return b;
    s.k = choice.k;
    return build(s);
}

}  // namespace ftab

#endif /* _FTAB_BUILD_HPP */
