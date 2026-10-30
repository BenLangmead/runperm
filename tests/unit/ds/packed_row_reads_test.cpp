// Unit tests for packed_matrix's whole-row reads, its reader and its prefetch
// hints, checked against get().
// Simple assert-based tests, no external framework.

#include "orbit/internal/ds/packed_vector.hpp"

#include <array>
#include <cassert>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>

using std::array;
using std::size_t;
using std::vector;

using namespace orbit;

constexpr size_t word_bits = num_bits_type(packed_matrix<1>::word_t);

// Fill every column of every row with a random value of its width.
template <size_t NumCols>
static void fill_random(packed_matrix<NumCols> &m, std::mt19937_64 &rng) {
    const auto &widths = m.get_widths();
    for (size_t row = 0; row < m.rows(); ++row) {
        array<ulint, NumCols> values{};
        for (size_t c = 0; c < NumCols; ++c) {
            values[c] = rng() & mask(widths[c]);
        }
        m.set_row(row, values);
    }
}

// Every read path for one row and column must agree with get().
template <size_t NumCols, size_t Col>
static void assert_column_reads_match(const packed_matrix<NumCols> &m,
                                      const typename packed_matrix<NumCols>::reader &r,
                                      size_t row) {
    const ulint expected = m.template get<Col>(row);
    assert(r.template get<Col>(row) == expected);
    assert(r.template get_at<Col>(r.row_start(row)) == expected);
    if (m.row_fits_word()) {
        assert(m.template extract<Col>(m.get_row_bits(row)) == expected);
    }
    // A span from column 0 to Col read with one load, when it fits.
    if (r.template span_fits<0, Col>()) {
        const ulint span = r.template get_span<0>(r.row_start(row));
        assert((r.template extract_span<0, Col>(span)) == expected);
    }
    // A span starting at Col itself always fits, since one column does.
    assert((r.template span_fits<Col, Col>()));
    const ulint own = r.template get_span<Col>(r.row_start(row));
    assert((r.template extract_span<Col, Col>(own)) == expected);
}

template <size_t NumCols, size_t... Cols>
static void assert_row_reads_match(const packed_matrix<NumCols> &m, size_t row,
                                   std::index_sequence<Cols...>) {
    const auto r = m.get_reader();
    (assert_column_reads_match<NumCols, Cols>(m, r, row), ...);
}

/*
 * get(), get_row_bits() and the reader's reads copy a whole ulint from the
 * byte holding a column's or row's first bit, so a read of the last row's
 * columns can touch up to sizeof(ulint) - 1 bytes past the last byte holding
 * row data.  These reads stay inside the buffer only because packed_matrix
 * allocates sizeof(ulint) bytes of padding past the last row (data_size()).
 * The prefetch hints likewise name the byte sizeof(ulint) - 1 past the start
 * of the row after the one prefetched, which for the last row lies in that
 * padding.  This checks, from the reader's copy of the layout, that every
 * such byte is below data_size(), so a change to the padding fails here
 * rather than as an out-of-bounds read.
 */
template <size_t NumCols>
static void assert_last_row_within_padding(const packed_matrix<NumCols> &m) {
    if (m.rows() == 0) return;
    const auto r = m.get_reader();
    const size_t last = m.rows() - 1;
    const size_t limit = m.data_size();
    for (size_t c = 0; c < NumCols; ++c) {
        const size_t first_byte = (r.row_start(last) + r.offsets[c]) / word_bits;
        assert(first_byte + sizeof(ulint) <= limit);
    }
    // get_row_bits() and get_span<0>() load from the row's first byte.
    assert(r.row_start(last) / word_bits + sizeof(ulint) <= limit);
    // The last byte prefetch(last) and prefetch_rows(.., last) name.
    const size_t end = r.row_start(last) + r.row_width;
    assert(end / word_bits + sizeof(ulint) - 1 < limit);
}

template <size_t NumCols>
static void check_matrix(const packed_matrix<NumCols> &m) {
    assert_last_row_within_padding(m);
    for (size_t row = 0; row < m.rows(); ++row) {
        assert_row_reads_match(m, row, std::make_index_sequence<NumCols>{});
    }
    if (m.rows() == 0) return;
    // Prefetches never change results; issue them on every row, the last
    // row included, and on ranges ending at the last row, then read again.
    const auto r = m.get_reader();
    const size_t last = m.rows() - 1;
    for (size_t row = 0; row < m.rows(); ++row) {
        m.prefetch(row);
        r.prefetch(row);
        m.prefetch_rows(row, last);
        r.prefetch_rows(row, last);
        m.prefetch_rows(row, row);
    }
    m.prefetch_rows(0, last);
    assert_row_reads_match(m, last, std::make_index_sequence<NumCols>{});
}

// Rows of many widths and row counts, so that the last row starts at every
// bit offset within a byte and rows straddle cache lines.
static void test_row_reads_match_get() {
    std::mt19937_64 rng(7);
    const vector<array<uchar, 3>> widths_list = {
        {1, 1, 1}, {3, 5, 7}, {7, 7, 7}, {13, 1, 9}, {20, 20, 17},
        {31, 13, 13}, {40, 9, 8}, {57, 57, 57}, {1, 56, 1}, {33, 33, 33},
    };
    for (const auto &widths : widths_list) {
        for (size_t rows : {1, 2, 3, 7, 8, 9, 63, 64, 65, 200}) {
            packed_matrix<3> m(rows, widths);
            fill_random(m, rng);
            check_matrix(m);
        }
    }
}

// Whole-row reads apply only to rows of at most max_width bits.
static void test_row_fits_word() {
    const uchar max = packed_matrix<2>::max_width;
    packed_matrix<2> fits(4, array<uchar, 2>{static_cast<uchar>(max - 10), 10});
    assert(fits.row_fits_word());
    packed_matrix<2> wide(4, array<uchar, 2>{static_cast<uchar>(max - 10), 11});
    assert(!wide.row_fits_word());

    std::mt19937_64 rng(11);
    fill_random(fits, rng);
    fill_random(wide, rng);
    check_matrix(fits);
    check_matrix(wide);
}

// Columns of width 0 take no bits and read as 0 through every path, and a
// matrix whose columns are all of width 0 has rows of width 0, whose
// prefetches name the padding.
static void test_zero_width_columns() {
    std::mt19937_64 rng(13);
    for (const auto &widths : vector<array<uchar, 3>>{
             {0, 9, 4}, {9, 0, 4}, {9, 4, 0}, {0, 0, 5}, {0, 0, 0}}) {
        for (size_t rows : {1, 5, 64}) {
            packed_matrix<3> m(rows, widths);
            fill_random(m, rng);
            for (size_t row = 0; row < rows; ++row) {
                for (size_t c = 0; c < 3; ++c) {
                    if (widths[c] != 0) continue;
                    const ulint v = c == 0 ? m.get<0>(row) : c == 1 ? m.get<1>(row) : m.get<2>(row);
                    assert(v == 0);
                }
            }
            check_matrix(m);
        }
    }
}

// A matrix with no rows keeps sizeof(ulint) bytes of padding.
static void test_empty_matrix() {
    packed_matrix<2> m(0, array<uchar, 2>{5, 7});
    assert(m.rows() == 0);
    assert(m.data_size() == sizeof(ulint));
    check_matrix(m);
}

// The padding survives serialization: a loaded matrix sizes its buffer the
// same way, so its last row reads the same.
static void test_padding_after_load() {
    std::mt19937_64 rng(17);
    packed_matrix<3> m(37, array<uchar, 3>{11, 23, 5});
    fill_random(m, rng);

    std::stringstream ss;
    m.serialize(ss);
    packed_matrix<3> loaded;
    loaded.load(ss);
    assert(loaded.rows() == m.rows());
    assert(loaded.data_size() == m.data_size());
    check_matrix(loaded);
    assert(loaded.get_row(m.rows() - 1) == m.get_row(m.rows() - 1));
}

int main() {
    test_row_reads_match_get();
    test_row_fits_word();
    test_zero_width_columns();
    test_empty_matrix();
    test_padding_after_load();

    std::cout << "packed row read tests passed" << std::endl;
    return 0;
}
