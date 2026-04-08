/**
 * Bidirectional index serialization.
 * Base index (.bi) + k-mer tables (.bi.kN), fingerprint-linked.
 *
 * Author: Ben Langmead (ben.langmead@gmail.com)
 * Date: Feb 23, 2026
 */

#ifndef _BI_IO_HPP
#define _BI_IO_HPP

#include "bi_fmd.hpp"
#include <optional>
#include <string>
#include <vector>

namespace bi_io {

/**
 * Write a base index to a file.
 * Returns true on success.
 */
bool write_base(const std::string& path, const FMDIndex& index);

/**
 * Read a base index from a file. Checks magic and version to ensure it matches expected.
 * Returns std::nullopt on error.
 */
std::optional<FMDIndex> read_base(const std::string& path);

/**
 * Write a k-mer table to a file. Converts Position-based BiIntervals to flat
 * (lo, lo_r, size) for storage. Requires fmd for Position->flat conversion.
 * Returns true on success.
 */
bool write_kmer_table(const std::string& path, uint64_t fingerprint,
                      const KmerTable& table, const FMDIndex& fmd);

/**
 * Read a k-mer table from a file. Checks fingerprint to ensure it matches
 * expected_fingerprint from the fmd. Converts flat BiIntervals to Position-based.
 * Returns std::nullopt if fingerprint mismatch or I/O error.
 */
std::optional<KmerTable> read_kmer_table(const std::string& path,
                                         uint64_t expected_fingerprint,
                                         const FMDIndex& fmd);

/**
 * Load all available k-mer tables from a base path. Checks fingerprints on each
 * to ensure they match fmd. Sorts tables descending by K. If no K=0 table was loaded, synthesizes a trivial one.
 */
std::vector<KmerTable> load_available_tables(const std::string& base_path,
                                             uint64_t fingerprint,
                                             const FMDIndex& fmd);

/**
 * Load specific k-mer tables from a base path. Checks fingerprints on each to
 * ensure they match fmd. If no K=0 table was loaded, synthesizes a trivial one.
 */
std::vector<KmerTable> load_tables(const std::string& base_path,
                                   uint64_t fingerprint,
                                   const FMDIndex& fmd,
                                   const std::vector<int>& desired_Ks);

}  // namespace bi_io

#endif /* _BI_IO_HPP */
