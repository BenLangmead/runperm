/**
 * Shared constants for MS index and TSV parsing.
 *
 * Author: Ben Langmead (ben.langmead@gmail.com)
 * Date: Feb 17, 2026
 */

#ifndef _MS_CONSTANTS_HPP
#define _MS_CONSTANTS_HPP

#include "orbit/common.hpp"

using ulint = orbit::ulint;

/** Sentinel for LCP values that are gaps (e.g. "-" in TSV, compressed-out values). */
constexpr ulint LCP_GAP = (1ULL << 57) - 1;

#endif /* _MS_CONSTANTS_HPP */
