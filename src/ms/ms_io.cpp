/**
 * Serialization for the compressed-LCP MS index: an ORBITIDX file (see
 * index_format.hpp) written by MSIndexSpillLCP::serialize.
 *
 * Author: Ben Langmead (ben.langmead@gmail.com)
 * Date: Feb 17, 2026
 */

#include "serialize.hpp"
#include "ms_rlbwt.hpp"
#include <fstream>
#include <iostream>

namespace ms_serialize {

/**
 * Write the index to a single binary file.  Loading needs no rebuild, and the
 * file is about the size of the index in memory.
 */
bool write_index(const std::string& path, MSIndexSpillLCP<false>& index) {
    std::ofstream out(path, std::ios::binary);
    if (!out.good()) return false;
    try {
        index.serialize(out);
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
        return false;
    }
    return out.good();
}

/**
 * Read an index written by write_index.  Returns nullopt, after printing why
 * to stderr, if the file cannot be opened, is not an ms index, or is
 * malformed.
 */
std::optional<MSIndexSpillLCP<false>> read_index(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in.good()) return std::nullopt;
    MSIndexSpillLCP<false> idx;
    try {
        idx.load(in);
    } catch (const std::exception& e) {
        std::cerr << path << ": " << e.what() << "\n";
        return std::nullopt;
    }
    return idx;
}

}  // namespace ms_serialize
