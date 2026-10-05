/**
 * Serialization for MS spill index, as an ORBITIDX file (see
 * index_format.hpp).
 *
 * Author: Ben Langmead (ben.langmead@gmail.com)
 * Date: Feb 17, 2026
 */

#ifndef _MS_SERIALIZE_HPP
#define _MS_SERIALIZE_HPP

#include "ms_rlbwt.hpp"
#include <optional>
#include <string>

namespace ms_serialize {

/**
 * Write MSIndexSpillLCP<false> (lengths-based) to a binary file.
 * Returns true on success.
 */
bool write_index(const std::string& path, MSIndexSpillLCP<false>& index);

/**
 * Read MS spill index from binary file.  Returns nullopt on error.
 */
std::optional<MSIndexSpillLCP<false>> read_index(const std::string& path);

}  // namespace ms_serialize

#endif /* _MS_SERIALIZE_HPP */
