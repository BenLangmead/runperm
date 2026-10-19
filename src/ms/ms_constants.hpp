/**
 * Shared constants for the MS index.
 *
 * Author: Ben Langmead (ben.langmead@gmail.com)
 * Date: Feb 17, 2026
 */

#ifndef _MS_CONSTANTS_HPP
#define _MS_CONSTANTS_HPP

#include "orbit/common.hpp"

using ulint = orbit::ulint;

/** Sentinel for LCP values that are gaps (compressed-out values). */
constexpr ulint LCP_GAP = (1ULL << 57) - 1;

#endif /* _MS_CONSTANTS_HPP */
