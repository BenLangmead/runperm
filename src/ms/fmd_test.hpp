/**
 * Tests for indexes over string collections with one terminator per string,
 * as read from a ropebwt3 FMD (whose LF may have several cycles, and whose N
 * sorts after T and matches nothing), and for BWTs longer than 2^41.
 */

#ifndef _FMD_TEST_HPP
#define _FMD_TEST_HPP

#include <string>

namespace fmd_test {

/** Run the tests; returns true if all ran (failures assert). */
bool run_all_tests(const std::string& data_dir);

}  // namespace fmd_test

#endif /* _FMD_TEST_HPP */
