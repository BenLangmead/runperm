/**
 * Bit-usage inspection for serialized MS index.
 * 
 * Author: Ben Langmead (ben.langmead@gmail.com)
 * Date: Feb 17, 2026
 */

#ifndef _MS_INSPECT_HPP
#define _MS_INSPECT_HPP

#include "ms_rlbwt.hpp"
#include <string>

namespace ms_inspect {

/**
 * Compute bit usage for different components of the index and print a brief
 * report to stdout.
 */
void run_inspect(MSIndexSpillLCP<false>& index);

/**
 * Dump spillover table to TSV: jumbo (T/F), num_interior, rest.
 * rest is "-" when num_interior=0, else "offset0,value0,offset1,value1,...".
 * Returns true on success.
 */
bool run_spillover_tsv(MSIndexSpillLCP<false>& index, const std::string& path);

}  // namespace ms_inspect

#endif /* _MS_INSPECT_HPP */
