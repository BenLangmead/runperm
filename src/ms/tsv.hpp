/**
 * Parse movify.py TSV format.
 * TSV columns: id, off, len, c, sa, lcp (lcp is comma-separated list).
 * Maps movify alphabet (#ACGT) to runperm Nucleotide (SEP, A, C, G, T).
 *
 * Author: Ben Langmead (ben.langmead@gmail.com)
 * Date: Feb 17, 2026
 */

#ifndef _MS_TSV_HPP
#define _MS_TSV_HPP

#include "orbit/common.hpp"
#include <string>
#include <vector>

using uchar = orbit::uchar;
using ulint = orbit::ulint;

namespace tsv {

/** Map movify '#' to SEP (1), keep 'A','C','G','T' as-is */
uchar map_movify_char(char c);

/** Parse LCP row "0,0,2465,..." or "0,-,2465" (gap as -) -> vector<ulint> */
std::vector<ulint> parse_lcp_row(const std::string& s);

/** Load TSV and populate bwt_heads, bwt_run_lengths, lcps_per_run. Returns true on success. */
bool load_tsv(
    const std::string& path,
    std::vector<uchar>& bwt_heads,
    std::vector<ulint>& bwt_run_lengths,
    std::vector<std::vector<ulint>>& lcps_per_run);

}  // namespace tsv

#endif /* _MS_TSV_HPP */
