/**
 * The alphabet of the ms and tms indexes: Orbit's nucleotide codes
 * (terminator, separator, A, C, G, T, N) and one more, TEXT_N, for N in a
 * text read from a ropebwt3 FMD.
 *
 * Orbit orders the BWT's characters by byte value when it builds LF, so a
 * run head's byte must sort where its symbol sorts in the BWT.  ropebwt3
 * orders its symbols $ < A < C < G < T < N, with N last, while the byte 'N'
 * sorts before 'T'; TEXT_N (0x7f) sorts after 'T' as ropebwt3's N does.  It
 * is also not a byte a read holds (reads are upper-cased letters), so N in
 * such a text never matches a read character, as in ropebwt3, where N
 * matches nothing.  The byte 'N' keeps its own code for texts built with
 * byte order, where an 'N' in a read matches an 'N' in the text.
 *
 * Like orbit::nucleotide, the mapping is fixed, so nothing is serialized and
 * indexes written with orbit::nucleotide load unchanged.
 */

#ifndef _MS_ALPHABET_HPP
#define _MS_ALPHABET_HPP

#include "orbit/common.hpp"
#include "orbit/internal/ds/alphabet.hpp"
#include <array>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ms_alpha {

using uchar = orbit::uchar;
using ulint = orbit::ulint;

/** The index byte of N in a text read from an FMD (see the file comment). */
inline constexpr uchar TEXT_N = 0x7f;

namespace detail {
/** Byte to code: TER, SEP, A, C, G, T, N (either case), TEXT_N; others unmapped. */
constexpr std::array<uchar, 256> make_map() {
    std::array<uchar, 256> m{};
    for (auto& x : m) x = orbit::UNMAPPED;
    m[orbit::TERMINATOR] = 0;
    m[orbit::SEPARATOR] = 1;
    m['A'] = 2; m['C'] = 3; m['G'] = 4; m['T'] = 5; m['N'] = 6;
    m['a'] = 2; m['c'] = 3; m['g'] = 4; m['t'] = 5; m['n'] = 6;
    m[TEXT_N] = 7;
    return m;
}
}  // namespace detail

class alphabet {
public:
    /** Number of codes: TER, SEP, A, C, G, T, N and TEXT_N. */
    static constexpr uchar SIGMA = 8;
    static constexpr uchar TER = 0, SEP = 1, A = 2, C = 3, G = 4, T = 5, N = 6, XN = 7;

    alphabet() = default;

    /** Check that every byte with a nonzero count has a code. */
    explicit alphabet(const std::vector<ulint>& char_counts) {
        for (size_t i = 0; i < char_counts.size(); ++i)
            if (char_counts[i] > 0 && map_[i] == NUL)
                throw std::runtime_error("ms alphabet has no code for byte " + std::to_string(i));
    }

    static constexpr uchar map_char(uchar c) { return map_[c]; }
    static constexpr uchar unmap_char(uchar c) { return unmap_[c]; }
    static constexpr uchar size() { return SIGMA; }

    size_t serialize(std::ostream&) { return 0; }
    void load(std::istream&) {}

private:
    static constexpr uchar NUL = orbit::UNMAPPED;
    static constexpr std::array<uchar, 256> map_ = detail::make_map();
    static constexpr std::array<uchar, SIGMA> unmap_ = {orbit::TERMINATOR, orbit::SEPARATOR, 'A', 'C', 'G', 'T', 'N',
                                                        TEXT_N};
};

}  // namespace ms_alpha

#endif /* _MS_ALPHABET_HPP */
