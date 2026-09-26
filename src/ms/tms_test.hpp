/**
 * Tests for TmsIndex and tms_query, and small-text helpers they share.
 */

#ifndef _TMS_TEST_HPP
#define _TMS_TEST_HPP

#include "orbit/common.hpp"
#include <random>
#include <string>
#include <vector>

namespace tms_test {

using uchar = orbit::uchar;
using ulint = orbit::ulint;

/**
 * A small text with its suffix array, LCP array and run-length BWT.  str is
 * the text as given plus a final '$'; text holds index codes, with '$' as
 * orbit::TERMINATOR and '%' as orbit::SEPARATOR.
 */
struct TextBwt {
    std::string str;
    std::vector<uchar> text;
    std::vector<ulint> sa, isa, lcp;
    std::vector<uchar> heads;
    std::vector<ulint> lens;
    std::vector<std::vector<ulint>> lcps_per_run;
};

/** Build a TextBwt from s, which must not contain '$'.  Quadratic worst case. */
TextBwt make_text_bwt(const std::string& s);

/** Texts from the families that found bugs in TeraMS, plus haplotype sets. */
std::vector<std::string> fuzz_texts(std::mt19937& rng);

/** Matching statistics of P against T by substring search. */
std::vector<ulint> naive_ms(const std::string& T, const std::string& P);

bool run_all_tests(const std::string& data_dir);

}  // namespace tms_test

#endif /* _TMS_TEST_HPP */
