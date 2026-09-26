/**
 * Appending integers to output text without a temporary string per value.
 */

#ifndef _FORMAT_INT_HPP
#define _FORMAT_INT_HPP

#include <charconv>
#include <cstdint>
#include <string>

/** Appends the decimal digits of v to s. */
inline void append_uint(std::string& s, uint64_t v) {
    char buf[20];
    const auto r = std::to_chars(buf, buf + sizeof(buf), v);
    s.append(buf, r.ptr);
}

#endif /* _FORMAT_INT_HPP */
