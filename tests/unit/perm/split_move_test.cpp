// Unit tests for the two-half move step (start_next and finish_next), the
// whole-row and reader reads of a permutation's packed rows, and its prefetch
// hints, checked against next(), get_length, get_start and get.
// Simple assert-based tests, no external framework.

#include "orbit/permutation.hpp"
#include "orbit/rlbwt.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <random>
#include <vector>

using std::size_t;
using std::vector;

using namespace orbit;

enum class TestRunCols {
    VAL1,
    VAL2,
    COUNT
};

using TestRunData = columns_tuple<TestRunCols>;

// A permutation of [0, n) made of runs: [0, n) is cut into blocks of random
// lengths, which are then shuffled, so that next() often fast-forwards.
static vector<ulint> make_runny_permutation(size_t n, size_t max_run, std::mt19937 &rng) {
    vector<std::pair<ulint, ulint>> blocks;
    for (ulint start = 0; start < n;) {
        ulint len = std::min<ulint>(1 + rng() % max_run, n - start);
        blocks.push_back({start, len});
        start += len;
    }
    std::shuffle(blocks.begin(), blocks.end(), rng);
    vector<ulint> perm;
    perm.reserve(n);
    for (auto [start, len] : blocks) {
        for (ulint i = 0; i < len; ++i) perm.push_back(start + i);
    }
    return perm;
}

// The position of global index idx, from the intervals' lengths.
template <bool Absolute, typename RP>
static typename RP::position make_pos(const RP &rp, ulint idx) {
    typename RP::position pos{};
    ulint prefix = 0;
    for (ulint interval = 0; interval < rp.intervals(); ++interval) {
        ulint len = rp.get_length(interval);
        if (idx < prefix + len) {
            pos.interval = interval;
            pos.offset = idx - prefix;
            if constexpr (Absolute) pos.idx = idx;
            return pos;
        }
        prefix += len;
    }
    assert(false && "index out of range");
    return pos;
}

template <bool Absolute, typename Pos>
static void assert_same_position(const Pos &a, const Pos &b) {
    assert(a.interval == b.interval);
    assert(a.offset == b.offset);
    if constexpr (Absolute) assert(a.idx == b.idx);
}

// The global index of a position, from the intervals' lengths.
template <typename RP>
static ulint global_index(const RP &rp, ulint interval, ulint offset) {
    ulint prefix = 0;
    for (ulint i = 0; i < interval; ++i) prefix += rp.get_length(i);
    return prefix + offset;
}

/*
 * start_next followed by finish_next must give next(), whether each
 * position is finished right away or a batch of positions is started,
 * their landing rows prefetched, and then finished.
 */
template <bool Absolute, typename RP>
static void check_two_half_move(RP &rp, const vector<ulint> &perm) {
    using position = typename RP::position;
    const ulint domain = rp.domain();
    for (ulint idx = 0; idx < domain; ++idx) {
        const position pos = make_pos<Absolute>(rp, idx);
        const position expected = rp.next(pos);
        assert(global_index(rp, expected.interval, expected.offset) == perm[idx]);

        const position started = rp.start_next(pos);
        assert(started.interval < rp.intervals());
        rp.prefetch(started.interval);
        assert_same_position<Absolute>(rp.finish_next(started), expected);
    }

    // Several positions in flight, as a batched caller keeps them.
    constexpr size_t batch = 5;
    for (ulint base = 0; base < domain; base += batch) {
        const size_t k = std::min<size_t>(batch, domain - base);
        std::array<position, batch> started{};
        for (size_t j = 0; j < k; ++j) {
            started[j] = rp.start_next(make_pos<Absolute>(rp, base + j));
            rp.prefetch(started[j].interval);
        }
        for (size_t j = 0; j < k; ++j) {
            const position expected = rp.next(make_pos<Absolute>(rp, base + j));
            assert_same_position<Absolute>(rp.finish_next(started[j]), expected);
        }
    }

    // Following a chain of steps, each split in two, matches next() steps.
    position a = rp.first();
    position b = rp.first();
    for (ulint step = 0; step < 2 * domain; ++step) {
        a = rp.next(a);
        position u = rp.start_next(b);
        rp.prefetch(u.interval);
        b = rp.finish_next(u);
        assert_same_position<Absolute>(a, b);
    }
}

/*
 * The reader's and the whole-row reads of each interval's row must agree
 * with get_length or get_start, with get for integrated data columns, and
 * with the pointer and offset that start_next follows from the interval's
 * first position.  Prefetches on every row, the last included, change none
 * of it.
 */
template <bool Absolute, bool Integrated, typename RP>
static void check_row_reads(RP &rp) {
    using position = typename RP::position;
    const auto r = rp.get_reader();
    const ulint last = rp.intervals() - 1;
    for (ulint i = 0; i < rp.intervals(); ++i) {
        rp.prefetch(i);
        rp.prefetch_rows(i, last);
        rp.prefetch_rows(i, i);
        r.prefetch(i);
        r.prefetch_rows(i, last);

        position head{};
        head.interval = i;
        head.offset = 0;
        if constexpr (Absolute) {
            head.idx = rp.get_start(i);
            assert(r.template get<RP::start_column()>(i) == rp.get_start(i));
        } else {
            assert(r.template get<RP::length_column()>(i) == rp.get_length(i));
        }
        const position started = rp.start_next(head);
        const ulint pointer = r.template get<RP::pointer_column()>(i);
        const ulint offset = r.template get<RP::offset_column()>(i);
        assert(pointer == started.interval);
        assert(offset == started.offset);

        if (rp.row_fits_word()) {
            const ulint bits = rp.row_bits(i);
            assert(rp.pointer_of(bits) == pointer);
            assert(rp.offset_of(bits) == offset);
            if constexpr (!Absolute) assert(rp.length_of(bits) == rp.get_length(i));
        }

        if constexpr (Integrated) {
            assert(r.template get<RP::template data_column<TestRunCols::VAL1>()>(i) ==
                   rp.template get<TestRunCols::VAL1>(i));
            assert(r.template get<RP::template data_column<TestRunCols::VAL2>()>(i) ==
                   rp.template get<TestRunCols::VAL2>(i));
        }
    }
    rp.prefetch_rows(0, last);
}

// A permutation with run data, built with or without splitting.
template <typename RP>
static RP build_with_data(const vector<ulint> &perm, const split_params &sp) {
    using Enc = typename RP::interval_encoding_t;
    Enc enc = Enc::from_permutation(perm, sp);
    vector<TestRunData> run_data(enc.intervals());
    for (size_t i = 0; i < run_data.size(); ++i) {
        run_data[i] = {static_cast<ulint>(i), static_cast<ulint>(3 * i + 1000)};
    }
    return RP(enc, run_data);
}

// A permutation with empty data columns, built with or without splitting.
template <typename RP>
static RP build_without_data(const vector<ulint> &perm, const split_params &sp) {
    using Enc = typename RP::interval_encoding_t;
    Enc enc = Enc::from_permutation(perm, sp);
    return RP(enc);
}

static void test_permutations() {
    std::mt19937 rng(23);
    const vector<split_params> splits = {NO_SPLITTING, split_params(2.0, 2), DEFAULT_SPLITTING};
    for (size_t n : {1, 2, 17, 64, 300}) {
        for (size_t max_run : {1, 4, 40}) {
            const vector<ulint> perm = make_runny_permutation(n, max_run, rng);
            for (const auto &sp : splits) {
                {
                    auto rp = build_with_data<permutation_separated<TestRunCols>>(perm, sp);
                    check_two_half_move<false>(rp, perm);
                    check_row_reads<false, false>(rp);
                }
                {
                    auto rp = build_with_data<permutation_integrated<TestRunCols>>(perm, sp);
                    check_two_half_move<false>(rp, perm);
                    check_row_reads<false, true>(rp);
                }
                {
                    auto rp = build_with_data<permutation_separated_absolute<TestRunCols>>(perm, sp);
                    check_two_half_move<true>(rp, perm);
                    check_row_reads<true, false>(rp);
                }
                {
                    auto rp = build_with_data<permutation_integrated_absolute<TestRunCols>>(perm, sp);
                    check_two_half_move<true>(rp, perm);
                    check_row_reads<true, true>(rp);
                }
                // With empty data columns the separate data table has no
                // rows, and prefetch covers only the move structure's rows.
                {
                    auto rp = build_without_data<permutation_separated<>>(perm, sp);
                    check_two_half_move<false>(rp, perm);
                    check_row_reads<false, false>(rp);
                }
                {
                    auto rp = build_without_data<permutation_separated_absolute<>>(perm, sp);
                    check_two_half_move<true>(rp, perm);
                    check_row_reads<true, false>(rp);
                }
                {
                    auto rp = build_without_data<move_permutation_relative>(perm, sp);
                    check_two_half_move<false>(rp, perm);
                    check_row_reads<false, false>(rp);
                }
                {
                    auto rp = build_without_data<move_permutation_absolute>(perm, sp);
                    check_two_half_move<true>(rp, perm);
                    check_row_reads<true, false>(rp);
                }
            }
        }
    }
}

// LF over a small RLBWT: the two-half move matches LF, and the character
// read from a whole row or the reader matches get_character.
static void test_lf_rows() {
    using namespace orbit::rlbwt;
    const vector<uchar> bwt_heads =       {'T','C','G','A','T', 1 ,'A','T','A'};
    const vector<ulint> bwt_run_lengths = { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };

    lf_move<> lf(bwt_heads, bwt_run_lengths);
    using position = typename lf_move<>::position;

    const auto r = lf.get_reader();
    for (ulint i = 0; i < lf.intervals(); ++i) {
        lf.prefetch(i);
        const uchar c = lf.get_character(i);
        if (lf.row_fits_word()) {
            assert(lf.character_of(lf.row_bits(i)) == c);
        }
        const auto code = lf.character_code(c);
        assert(code.has_value());
        assert(r.template get<lf_move<>::character_column()>(i) == *code);
    }
    assert(!lf.character_code('Z').has_value());

    position a = lf.first();
    position b = lf.first();
    for (ulint step = 0; step < 3 * lf.domain(); ++step) {
        a = lf.LF(a);
        position u = lf.start_next(b);
        lf.prefetch(u.interval);
        b = lf.finish_next(u);
        assert(a.interval == b.interval);
        assert(a.offset == b.offset);
    }
}

int main() {
    test_permutations();
    test_lf_rows();

    std::cout << "split move unit tests passed" << std::endl;
    return 0;
}
