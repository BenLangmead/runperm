/**
 * Recorded reposition outcomes, which the oracle modes of the ms engine
 * (ms_rlbwt.hpp) and the tms engine (tms_query.hpp) replay instead of
 * computing them.  Each engine replays its own recordings: when both sides'
 * LCEs reach the match length, the ms engine takes the side with the larger
 * uncapped LCP bound and the tms engine the side above, so their outcomes
 * can differ.
 */

#ifndef _REPOSITION_ORACLE_HPP
#define _REPOSITION_ORACLE_HPP

#include <cstddef>
#include <cstdint>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

/**
 * The outcome of one reposition: the LF interval the match moves to (the
 * nearest run of the character above the current interval, whose tail row
 * is taken, if it is above, else the nearest run below, whose head row is
 * taken) and the capped LCE, the match length before the LF step.
 */
struct RepositionDecision {
    uint64_t interval;
    uint64_t len;
};

/**
 * Recorded reposition outcomes for a sequence of patterns, for the oracle
 * modes of the ms and tms engines: pattern j's, in the order its query makes them, are
 * decisions[offsets[j] .. offsets[j + 1]).  move_runs is the LF interval
 * count of the index they were recorded on, since intervals are numbered
 * per index (and per layout).
 *
 * File format (all integers u64, little-endian): the magic "TMSORC01"
 * (8 bytes), move_runs, the pattern count p, the decision count d, then p
 * per-pattern decision counts, then d (interval, len) pairs.
 */
struct RepositionOracle {
    uint64_t move_runs = 0;
    std::vector<uint64_t> offsets{0};
    std::vector<RepositionDecision> decisions;

    size_t patterns() const { return offsets.size() - 1; }

    /** Append the decisions of the patterns in rec, in order. */
    void append(const std::vector<std::vector<RepositionDecision>>& rec) {
        for (const auto& v : rec) {
            decisions.insert(decisions.end(), v.begin(), v.end());
            offsets.push_back(decisions.size());
        }
    }

    void write(std::ostream& out) const {
        static_assert(sizeof(RepositionDecision) == 16, "decisions are written as two u64s");
        auto put = [&](uint64_t x) { out.write(reinterpret_cast<const char*>(&x), sizeof(x)); };
        out.write("TMSORC01", 8);
        put(move_runs);
        put(patterns());
        put(decisions.size());
        for (size_t j = 0; j < patterns(); ++j) put(offsets[j + 1] - offsets[j]);
        out.write(reinterpret_cast<const char*>(decisions.data()),
                  static_cast<std::streamsize>(decisions.size() * sizeof(RepositionDecision)));
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
        if (!in || std::string(magic, 8) != "TMSORC01") throw std::runtime_error("not an oracle file");
        move_runs = get();
        const uint64_t p = get(), d = get();
        if (!in) throw std::runtime_error("truncated oracle file");
        offsets.assign(1, 0);
        offsets.reserve(p + 1);
        for (uint64_t j = 0; j < p; ++j) offsets.push_back(offsets.back() + get());
        if (!in || offsets.back() != d) throw std::runtime_error("oracle file's per-pattern counts do not add up");
        decisions.resize(d);
        in.read(reinterpret_cast<char*>(decisions.data()), static_cast<std::streamsize>(d * sizeof(RepositionDecision)));
        if (!in) throw std::runtime_error("truncated oracle file");
    }
};

/**
 * Part of a RepositionOracle that a batch replays: its pattern j's decisions are
 * decisions[offsets[j] .. offsets[j + 1]).
 */
struct RepositionOracleView {
    const RepositionDecision* decisions;
    const uint64_t* offsets;
    uint64_t move_runs;
    /** Patterns [first, first + count) of oracle, which must hold them. */
    static RepositionOracleView of(const RepositionOracle& o, size_t first, size_t count) {
        if (first + count > o.patterns())
            throw std::runtime_error("the oracle holds fewer patterns than the query has");
        return RepositionOracleView{o.decisions.data(), o.offsets.data() + first, o.move_runs};
    }
};


#endif /* _REPOSITION_ORACLE_HPP */
