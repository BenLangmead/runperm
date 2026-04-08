/**
 * Bidirectional index test suite.
 *
 * Author: Ben Langmead (ben.langmead@gmail.com)
 * Date: Feb 23, 2026
 */

#ifndef _BI_TEST_HPP
#define _BI_TEST_HPP

#include <string>

namespace bi_test {

/// Run all built-in tests.  Returns true if all pass.
bool run_all_tests(const std::string& data_dir = "");

}  // namespace bi_test

#endif /* _BI_TEST_HPP */
