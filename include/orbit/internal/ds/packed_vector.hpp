#ifndef _PACKED_VECTOR_HPP
#define _PACKED_VECTOR_HPP

#include "orbit/common.hpp"
#include <cassert>
#include <cstring>
#include <cstdint>
#include <cstddef>
#include <vector>
#include <array>
#include <iostream>
#include <algorithm>
#include <iterator>

namespace orbit {

template <size_t num_cols>
class packed_matrix {
public:
    // uses ulint for unpacked data
    using word_t = uchar;

    // We read ulint at a time, this ensures we never need to read more than one ulint
    // should be 57 bits for 64 bit ulint and 8 bit word_t
    constexpr static uchar max_width = num_bits_type(ulint) - (num_bits_type(word_t) - 1);
    // Bytes per cache line, the stride of prefetch_rows.
    constexpr static size_t cache_line_bytes = 64;

    packed_matrix() = default;
    packed_matrix(const ulint rows, const std::array<uchar, num_cols>& widths) {
        packed_matrix::num_rows = rows;
        packed_matrix::widths = widths;

        init();
    }

    packed_matrix(packed_matrix&& other) noexcept = default;
    packed_matrix& operator=(packed_matrix&& other) noexcept = default;
    packed_matrix(const packed_matrix& other) = default;
    packed_matrix& operator=(const packed_matrix& other) = default;
    ~packed_matrix() = default;

    template<size_t col>
    ulint get(size_t row) const {
        static_assert(col < num_cols, "Column out of bounds");
        assert(row < num_rows);

        bit_pos pos(get_row_start(row) + offsets[col]);
        ulint bits = 0;
        std::memcpy(&bits, &data[pos.chunk], sizeof(ulint));
        return extract_bits(bits, pos.offset, masks_extract[col]);
    } 

    /** Whether get_row_bits can read a whole row with one load. */
    bool row_fits_word() const { return row_width <= max_width; }

    /**
     * The bits of a row, its first column lowest, read with one unaligned
     * load; only when row_fits_word().  Bits above the row are unspecified,
     * and extract() takes a column out of the result.
     */
    ulint get_row_bits(size_t row) const {
        assert(row < num_rows && row_fits_word());
        bit_pos pos(get_row_start(row));
        ulint bits = 0;
        std::memcpy(&bits, &data[pos.chunk], sizeof(ulint));
        return bits >> pos.offset;
    }

    /**
     * A copy of what reading columns needs: the buffer's address, the row
     * width, and each column's bit offset and mask.  A reader kept in a local
     * variable lets the compiler hold these in registers, where reads through
     * the matrix reload them after any store that might alias them.  It is
     * valid while the matrix is alive and not resized.  Rows are named by
     * their first bit, row_start(row), so that several columns of one row
     * share the multiplication.
     */
    struct reader {
        const word_t* data;
        size_t row_width;
        std::array<uint16_t, num_cols> offsets;
        std::array<ulint, num_cols> masks;

        size_t row_start(size_t row) const { return row * row_width; }
        template<size_t col>
        ulint get_at(size_t start) const {
            static_assert(col < num_cols, "Column out of bounds");
            const size_t bit = start + offsets[col];
            ulint bits = 0;
            std::memcpy(&bits, &data[bit / num_bits_type(word_t)], sizeof(ulint));
            return (bits >> (bit % num_bits_type(word_t))) & masks[col];
        }
        template<size_t col>
        ulint get(size_t row) const { return get_at<col>(row_start(row)); }
    };
    reader get_reader() const { return reader{data.data(), row_width, offsets, masks_extract}; }

    /** Column col of a row read with get_row_bits. */
    template<size_t col>
    ulint extract(ulint row_bits) const {
        static_assert(col < num_cols, "Column out of bounds");
        return extract_bits(row_bits, static_cast<uchar>(offsets[col]), masks_extract[col]);
    }

    /**
     * Hint that row will be read soon.  get() and get_row_bits() read a
     * whole ulint from the byte where a column or row starts, so the reads
     * of a row, and a get_row_bits() of the row after it (as a move's
     * fast-forward may do), span from the row's first byte to sizeof(ulint)
     * - 1 bytes past the start of the next row.  This prefetches the
     * lines holding the first and last bytes.  They are often the same line,
     * but testing for that costs more in mispredicted branches than the
     * second prefetch does.
     */
    void prefetch(size_t row) const {
        const size_t start = get_row_start(row);
        const word_t* first = &data[start / num_bits_type(word_t)];
        const word_t* last = &data[(start + row_width) / num_bits_type(word_t) + sizeof(ulint) - 1];
        ORBIT_PREFETCH(first);
        ORBIT_PREFETCH(last);
    }

    /**
     * Hint that rows lo to hi, lo <= hi, will be read soon: prefetches each
     * cache line from the first byte of row lo to the last byte a read of
     * row hi touches, once per line.
     */
    void prefetch_rows(size_t lo, size_t hi) const {
        const word_t* p = &data[get_row_start(lo) / num_bits_type(word_t)];
        const word_t* last = &data[(get_row_start(hi) + row_width) / num_bits_type(word_t) + sizeof(ulint) - 1];
        for (; p < last; p += cache_line_bytes) ORBIT_PREFETCH(p);
        ORBIT_PREFETCH(last);
    }

    template<size_t col>
    void set(size_t row, ulint val) {
        static_assert(col < num_cols, "Column out of bounds");
        assert(row < num_rows);
        assert(val < pow2(widths[col]));

        bit_pos pos(get_row_start(row) + offsets[col]);

        ulint bits = 0;
        std::memcpy(&bits, &data[pos.chunk], sizeof(ulint));
        write_bits(bits, pos.offset, masks_write[pos.offset][col], val);
        std::memcpy(&data[pos.chunk], &bits, sizeof(ulint));
    }

    template<size_t... col>
    void set_row(size_t row, const std::array<ulint, num_cols>& values, std::index_sequence<col...>) {
        (set<col>(row, values[col]), ...);
    }
    void set_row(size_t row, const std::array<ulint, num_cols>& values) {
        set_row(row, values, std::make_index_sequence<num_cols>{});
    }

    template<size_t... col>
    std::array<ulint, num_cols> get_row(size_t row, std::index_sequence<col...>) const { 
        return {get<col>(row)...};
    }
    std::array<ulint, num_cols> get_row(size_t row) const { 
        return get_row(row, std::make_index_sequence<num_cols>{});
    }


    [[nodiscard]] size_t size() const noexcept { return num_rows; }
    [[nodiscard]] size_t rows() const noexcept { return num_rows; }
    [[nodiscard]] static constexpr size_t cols() noexcept { return num_cols; }
    [[nodiscard]] size_t data_size() const noexcept {
        return ceil_div(vector_width, num_bits_type(word_t)) + sizeof(ulint)/sizeof(word_t);
    }
    [[nodiscard]] const std::array<uchar, num_cols>& get_widths() const noexcept { return widths; }

    size_t serialize(std::ostream &out) {
        size_t written_bytes = 0;

        written_bytes += write_magic(out, MAGIC);
        written_bytes += serialize_version(out);

        out.write((char *)&num_rows, sizeof(num_rows));
        written_bytes += sizeof(num_rows);

        // Serialize column widths (may be zero columns)
        size_t widths_bytes = widths.size() * sizeof(uchar);
        if (widths_bytes > 0) {
            out.write((char *)widths.data(), widths_bytes);
            written_bytes += widths_bytes;
        }

        // Serialize packed data buffer
        size_t data_bytes = data.size() * sizeof(word_t);
        if (data_bytes > 0) {
            out.write((char *)data.data(), data_bytes);
            written_bytes += data_bytes;
        }

        return written_bytes;
    }

    void load(std::istream &in) {
        check_magic(in, MAGIC);
        auto [serialized_major, serialized_minor, serialized_patch] = load_version(in);
        if (serialized_major != VERSION_MAJOR || serialized_minor != VERSION_MINOR || serialized_patch != VERSION_PATCH) {
            // TODO handle version mismatches
        }

        in.read((char *)&num_rows, sizeof(num_rows));
        // Load column widths (may be zero columns)
        size_t widths_bytes = widths.size() * sizeof(uchar);
        if (widths_bytes > 0) {
            in.read((char *)widths.data(), widths_bytes);
        }
        init();
        // Load packed data buffer
        size_t data_bytes = data_size() * sizeof(word_t);
        if (data_bytes > 0) {
            in.read((char *)data.data(), data_bytes);
        }
    }

private:
    // OrBit Packed Mector
    static constexpr std::array<char, MAGIC_BYTES> MAGIC = {'O', 'B', 'P', 'M'};

    size_t num_rows;
    size_t vector_width; // bit width of stored data (actual data size might be larger due to padding)
    size_t row_width; // width of each row in bits

    std::array<uchar, num_cols> widths; // Bit width of each column
    std::array<uint16_t, num_cols> offsets; // Offset of the first bit of each column in the vector
    std::array<ulint, num_cols> masks_extract; // Mask of the bits of each column for get
    std::array<std::array<ulint, num_cols>, num_bits_type(word_t)> masks_write; // The mask for each offset to reset the bits

    std::vector<word_t> data;

    void init() {
        size_t bit_pos = 0;
        for (size_t i = 0; i < num_cols; i++) {
            assert(widths[i] <= max_width);
            offsets[i] = bit_pos;
            masks_extract[i] = mask(widths[i]);
            for (size_t j = 0; j < num_bits_type(word_t); j++) {
                masks_write[j][i] = ~(masks_extract[i] << j);
            }
            bit_pos += widths[i];
        }
        row_width = bit_pos;
        vector_width = num_rows * row_width;
        data.resize(data_size());
    }

    struct bit_pos {
        ulint chunk;
        ulint offset;

        bit_pos(size_t bit) {
            chunk = bit / num_bits_type(word_t);
            offset = bit % num_bits_type(word_t);
        }
    };

    inline size_t get_row_start(size_t row) const {
        return row*row_width;
    }

    inline ulint extract_bits(ulint bits, uchar start, ulint mask) const {
        return (bits >> start) & mask;
    }

    inline void write_bits(ulint& bits, uchar start, ulint mask, ulint val) {
        // clear old bits
        bits &= mask;
        // set new bits
        bits |= (val << start);
    }

};

template<class columns>
class packed_vector : public packed_matrix<static_cast<size_t>(columns::COUNT)> {
    using base = packed_matrix<static_cast<size_t>(columns::COUNT)>;

public:
    packed_vector() = default;
    packed_vector(size_t rows, const std::array<uchar, static_cast<size_t>(columns::COUNT)>& widths)
        : base(rows, widths) {}

    template<columns col>
    ulint get(size_t row) const {
        return base::template get<static_cast<size_t>(col)>(row);
    }

    template<columns col>
    void set(size_t row, ulint val) {
        base::template set<static_cast<size_t>(col)>(row, val);
    }

    template<columns col>
    ulint extract(ulint row_bits) const {
        return base::template extract<static_cast<size_t>(col)>(row_bits);
    }
};

class int_vector : public packed_matrix<1> {
    using base = packed_matrix<1>;
public:
    using value_type = ulint;
    using size_type = size_t;
    using difference_type = std::ptrdiff_t;

    class reference {
    public:
        reference(int_vector* v, size_type i) : vec(v), idx(i) {}

        reference(const reference&) = default;
        reference& operator=(const reference& other) {
            return *this = static_cast<ulint>(other);
        }

        reference& operator=(ulint value) {
            vec->set(idx, value);
            return *this;
        }

        operator ulint() const {
            return vec->get(idx);
        }

    private:
        int_vector* vec;
        size_type idx;
    };

    class const_reference {
    public:
        const_reference(const int_vector* v, size_type i) : vec(v), idx(i) {}

        operator ulint() const {
            return vec->get(idx);
        }

    private:
        const int_vector* vec;
        size_type idx;
    };

    class iterator {
    public:
        using iterator_category = std::random_access_iterator_tag;
        using value_type        = int_vector::value_type;
        using difference_type   = int_vector::difference_type;
        using reference         = int_vector::reference;
        using pointer           = void;

        iterator() : vec(nullptr), idx(0) {}
        iterator(int_vector* v, size_type i) : vec(v), idx(i) {}

        reference operator*() const { return reference(vec, idx); }

        iterator& operator++() { ++idx; return *this; }
        iterator operator++(int) { iterator tmp(*this); ++(*this); return tmp; }

        iterator& operator--() { --idx; return *this; }
        iterator operator--(int) { iterator tmp(*this); --(*this); return tmp; }

        iterator& operator+=(difference_type n) { idx += n; return *this; }
        iterator& operator-=(difference_type n) { idx -= n; return *this; }

        iterator operator+(difference_type n) const { return iterator(vec, idx + n); }
        iterator operator-(difference_type n) const { return iterator(vec, idx - n); }

        difference_type operator-(const iterator& other) const {
            return static_cast<difference_type>(idx) - static_cast<difference_type>(other.idx);
        }

        bool operator==(const iterator& other) const { return vec == other.vec && idx == other.idx; }
        bool operator!=(const iterator& other) const { return !(*this == other); }
        bool operator<(const iterator& other) const { return idx < other.idx; }
        bool operator>(const iterator& other) const { return other < *this; }
        bool operator<=(const iterator& other) const { return !(other < *this); }
        bool operator>=(const iterator& other) const { return !(*this < other); }

        int_vector* vec;
        size_type idx;
    };

    class const_iterator {
    public:
        using iterator_category = std::random_access_iterator_tag;
        using value_type        = int_vector::value_type;
        using difference_type   = int_vector::difference_type;
        using reference         = int_vector::const_reference;
        using pointer           = void;

        const_iterator() : vec(nullptr), idx(0) {}
        const_iterator(const int_vector* v, size_type i) : vec(v), idx(i) {}
        const_iterator(const iterator& it) : vec(it.vec), idx(it.idx) {}

        reference operator*() const { return reference(vec, idx); }

        const_iterator& operator++() { ++idx; return *this; }
        const_iterator operator++(int) { const_iterator tmp(*this); ++(*this); return tmp; }

        const_iterator& operator--() { --idx; return *this; }
        const_iterator operator--(int) { const_iterator tmp(*this); --(*this); return tmp; }

        const_iterator& operator+=(difference_type n) { idx += n; return *this; }
        const_iterator& operator-=(difference_type n) { idx -= n; return *this; }

        const_iterator operator+(difference_type n) const { return const_iterator(vec, idx + n); }
        const_iterator operator-(difference_type n) const { return const_iterator(vec, idx - n); }

        difference_type operator-(const const_iterator& other) const {
            return static_cast<difference_type>(idx) - static_cast<difference_type>(other.idx);
        }

        bool operator==(const const_iterator& other) const { return vec == other.vec && idx == other.idx; }
        bool operator!=(const const_iterator& other) const { return !(*this == other); }
        bool operator<(const const_iterator& other) const { return idx < other.idx; }
        bool operator>(const const_iterator& other) const { return other < *this; }
        bool operator<=(const const_iterator& other) const { return !(other < *this); }
        bool operator>=(const const_iterator& other) const { return !(*this < other); }

        const int_vector* vec;
        size_type idx;
    };

    int_vector() = default;
    int_vector(size_t rows, uchar width) : base(rows, {width}) {}
    int_vector(std::vector<ulint> data) 
    : base(data.size(), {(data.empty() ? static_cast<uchar>(0) : bit_width(*std::max_element(data.begin(), data.end())))}) {
        for (size_t i = 0; i < data.size(); i++) {
            set(i, data[i]);
        }
    }
    int_vector(std::vector<ulint> data, uchar width)
    : base(data.size(), {width}) {
        for (size_t i = 0; i < data.size(); i++) {
            set(i, data[i]);
        }
    }
    
    ulint get(size_t row) const {
        return base::template get<0>(row);
    }

    void set(size_t row, ulint val) {
        base::template set<0>(row, val);
    }

    /// Bits per stored element (same as column 0 width in the underlying packed matrix).
    [[nodiscard]] uchar get_width() const noexcept {
        return base::get_widths()[0];
    }

    reference operator[](size_t row) {
        return reference(this, row);
    }

    ulint operator[](size_t row) const {
        return get(row);
    }

    iterator begin() {
        return iterator(this, 0);
    }

    iterator end() {
        return iterator(this, this->size());
    }

    const_iterator begin() const {
        return cbegin();
    }

    const_iterator end() const {
        return cend();
    }

    const_iterator cbegin() const {
        return const_iterator(this, 0);
    }

    const_iterator cend() const {
        return const_iterator(this, this->size());
    }
};

/** Prvalues from operator* are proxies; libc++ std::iter_swap calls swap(*a,*b) with rvalues. */
inline void swap(int_vector::reference&& x, int_vector::reference&& y) noexcept {
    ulint t = static_cast<ulint>(x);
    x = static_cast<ulint>(y);
    y = t;
}

} // namespace orbit

#endif // end of include guard: _PACKED_VECTOR_HPP