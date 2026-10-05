/**
 * Survey of k-mer presence for sizing an ftab, a table indexed by every
 * k-mer over ACGT that holds the state an MS query reaches after its
 * cold-start walk over the k-mer.  This is a development tool: it measures
 * how many k-mers occur in the text and how large a dense table would be,
 * and does not build the table.
 *
 * The walk over a k-mer reads it right to left from idx.first() with match
 * length 0, exactly as ms_query does, and the k-mer occurs in the text when
 * its statistic at the left end equals k.  The state after reading a suffix
 * depends only on that suffix, so k-mers that share a suffix share the walk
 * over it.  The survey therefore walks a trie of suffixes depth first,
 * extending each one by the four characters to its left, and every depth up
 * to kmax is counted in one pass.  A suffix that does not occur has no
 * occurring extension, so its subtree is skipped.
 */

#ifndef _FTAB_SURVEY_HPP
#define _FTAB_SURVEY_HPP

#include "ms_rlbwt.hpp"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <string>
#include <thread>
#include <vector>

namespace ftab_survey {

/** Number of bits needed to write any value in [0, v]. */
inline unsigned bits_for(ulint v) {
    unsigned b = 0;
    while (v > 0) { ++b; v >>= 1; }
    return b;
}

/**
 * One step of ms_query's walk: read character c to the left of the
 * suffix whose state is (pos, match_len).  Returns false if c does not
 * occur in the text.
 */
template <bool SP>
inline bool step(MSIndexSpillLCP<SP>& idx, typename MSIndexSpillLCP<SP>::position& pos, ulint& match_len, uchar c) {
    if (idx.get_character(pos.interval) == c) {
        pos = idx.LF(pos);
        ++match_len;
        return true;
    }
    auto opt = reposition_with_lcp(idx, pos.interval, static_cast<ulint>(pos.offset), c);
    if (!opt) {
        match_len = 0;
        return false;
    }
    match_len = std::min(match_len, opt->second) + 1;
    pos = opt->first;
    return true;
}

/** present[L]: how many L-mers over ACGT occur in the text, for L = 0 to kmax. */
struct Result {
    std::vector<ulint> present;
    ulint steps = 0;  // walk steps taken
    double seconds = 0.0;
};

/**
 * Count the occurring L-mers for every L up to kmax, on the given number of
 * threads.  The trie is expanded breadth first to a shallow depth on one
 * thread, and the subtrees below are shared among the threads.
 */
template <bool SP>
inline Result survey(MSIndexSpillLCP<SP>& idx, unsigned kmax, size_t threads) {
    using Position = typename MSIndexSpillLCP<SP>::position;
    static constexpr uchar acgt[4] = {'A', 'C', 'G', 'T'};
    struct Node {
        Position pos;
        ulint depth;  // also the match length, since the suffix occurs
    };
    const auto t0 = std::chrono::steady_clock::now();
    Result res;
    res.present.assign(kmax + 1, 0);
    res.present[0] = 1;
    // Extend node by each character; occurring children go to out and are
    // counted in present.
    auto expand = [&](const Node& nd, std::vector<Node>& out, std::vector<ulint>& present, ulint& steps) {
        for (uchar c : acgt) {
            Position p = nd.pos;
            ulint ml = nd.depth;
            ++steps;
            if (!step(idx, p, ml, c) || ml != nd.depth + 1) continue;
            ++present[nd.depth + 1];
            out.push_back(Node{p, nd.depth + 1});
        }
    };
    const unsigned seed_depth = std::min<unsigned>(kmax, 6);
    std::vector<Node> frontier{Node{idx.first(), 0}};
    for (unsigned d = 0; d < seed_depth; ++d) {
        std::vector<Node> next;
        for (const auto& nd : frontier) expand(nd, next, res.present, res.steps);
        frontier.swap(next);
    }
    if (threads == 0) threads = 1;
    std::atomic<size_t> next_task{0};
    std::vector<std::vector<ulint>> counts(threads, std::vector<ulint>(kmax + 1, 0));
    std::vector<ulint> steps(threads, 0);
    auto work = [&](size_t t) {
        std::vector<Node> stack;
        for (size_t j; (j = next_task.fetch_add(1)) < frontier.size();) {
            stack.push_back(frontier[j]);
            while (!stack.empty()) {
                const Node nd = stack.back();
                stack.pop_back();
                if (nd.depth < kmax) expand(nd, stack, counts[t], steps[t]);
            }
        }
    };
    std::vector<std::thread> pool;
    for (size_t t = 1; t < threads; ++t) pool.emplace_back(work, t);
    work(0);
    for (auto& th : pool) th.join();
    for (size_t t = 0; t < threads; ++t) {
        for (unsigned L = 0; L <= kmax; ++L) res.present[L] += counts[t][L];
        res.steps += steps[t];
    }
    res.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return res;
}

/**
 * The same count for one k by running ms_query from a cold start on every
 * k-mer, as a check of survey.  Practical only for small k.
 */
template <bool SP>
inline ulint brute_force(MSIndexSpillLCP<SP>& idx, unsigned k) {
    static constexpr char acgt[4] = {'A', 'C', 'G', 'T'};
    ulint present = 0;
    std::string kmer(k, 'A');
    for (ulint code = 0; code < (ulint{1} << (2 * k)); ++code) {
        for (unsigned i = 0; i < k; ++i) kmer[i] = acgt[(code >> (2 * i)) & 3];
        if (ms_query(idx, kmer)[0] == k) ++present;
    }
    return present;
}

/** The longest LF interval of idx, found by scanning every interval's length. */
template <bool SP>
inline ulint max_interval_length(MSIndexSpillLCP<SP>& idx, size_t threads) {
    const ulint n = idx.move_runs();
    if (threads == 0) threads = 1;
    std::vector<ulint> best(threads, 0);
    auto work = [&](size_t t) {
        const ulint lo = n / threads * t, hi = t + 1 == threads ? n : n / threads * (t + 1);
        ulint m = 0;
        for (ulint i = lo; i < hi; ++i) m = std::max(m, idx.get_length(i));
        best[t] = m;
    };
    std::vector<std::thread> pool;
    for (size_t t = 1; t < threads; ++t) pool.emplace_back(work, t);
    work(0);
    for (auto& th : pool) th.join();
    return *std::max_element(best.begin(), best.end());
}

}  // namespace ftab_survey

#endif /* _FTAB_SURVEY_HPP */
