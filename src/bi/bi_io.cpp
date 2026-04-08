/**
 * Bidirectional index serialization implementation.
 * Base: BIIX magic, fingerprint, BWT heads/lengths, C array.
 * K-mer tables: BIKT magic, fingerprint, K, 4^K BiIntervals (flat lo, lo_r, size).
 *
 * Author: Ben Langmead (ben.langmead@gmail.com)
 * Date: Feb 23, 2026
 */

#include "bi_io.hpp"
#include <fstream>
#include <cstring>
#include <algorithm>

namespace bi_io {

static const char BASE_MAGIC[] = "BIIX";
static const char KMER_MAGIC[] = "BIKT";
static constexpr uint32_t BASE_VERSION = 1;
static constexpr uint32_t KMER_VERSION = 1;
static constexpr uchar SIGMA = 7;

// Position <-> flat BWT offset conversion using run lengths
static ulint position_to_flat(const FMDIndex& fmd, Position pos) {
    const auto& lens = fmd.bwt_run_lengths();
    ulint flat = 0;
    for (ulint i = 0; i < pos.interval && i < lens.size(); ++i)
        flat += lens[i];
    return flat + pos.offset;
}

/**
 * Convert a flat BWT offset to a Position.
 */
static Position flat_to_position(const FMDIndex& fmd, ulint flat) {
    const auto& lens = fmd.bwt_run_lengths();
    ulint cumul = 0;
    for (ulint i = 0; i < lens.size(); ++i) {
        if (flat < cumul + lens[i])
            return Position{i, flat - cumul};
        cumul += lens[i];
    }
    return Position{static_cast<ulint>(lens.size()), 0};
}

/**
 * Write a base index to a file.
 */
bool write_base(const std::string& path, const FMDIndex& index) {
    std::ofstream out(path, std::ios::binary);
    if (!out.good()) return false;

    out.write(BASE_MAGIC, 4);
    uint32_t v = BASE_VERSION;
    out.write(reinterpret_cast<const char*>(&v), sizeof(v));

    uint64_t fp = index.get_fingerprint();
    out.write(reinterpret_cast<const char*>(&fp), sizeof(fp));

    ulint r = index.num_runs();
    ulint n = index.domain();
    out.write(reinterpret_cast<const char*>(&r), sizeof(r));
    out.write(reinterpret_cast<const char*>(&n), sizeof(n));
    out.write(reinterpret_cast<const char*>(&SIGMA), sizeof(SIGMA));

    const auto& heads = index.bwt_heads();
    const auto& lens = index.bwt_run_lengths();
    out.write(reinterpret_cast<const char*>(heads.data()), static_cast<std::streamsize>(r * sizeof(uchar)));
    out.write(reinterpret_cast<const char*>(lens.data()), static_cast<std::streamsize>(r * sizeof(ulint)));

    const auto& C = index.get_C();
    out.write(reinterpret_cast<const char*>(C.data()), static_cast<std::streamsize>((SIGMA + 1) * sizeof(ulint)));

    return out.good();
}

/**
 * Read a base index from a file.  Checks magic and version to ensure it matches
 * expected.
 */
std::optional<FMDIndex> read_base(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in.good()) return std::nullopt;

    char magic[5] = {};
    in.read(magic, 4);
    if (std::strncmp(magic, BASE_MAGIC, 4) != 0) return std::nullopt;

    uint32_t v = 0;
    in.read(reinterpret_cast<char*>(&v), sizeof(v));
    if (v != BASE_VERSION) return std::nullopt;

    uint64_t fp = 0;
    ulint r = 0, n = 0;
    uchar sigma = 0;
    in.read(reinterpret_cast<char*>(&fp), sizeof(fp));
    in.read(reinterpret_cast<char*>(&r), sizeof(r));
    in.read(reinterpret_cast<char*>(&n), sizeof(n));
    in.read(reinterpret_cast<char*>(&sigma), sizeof(sigma));
    if (sigma != SIGMA || r == 0) return std::nullopt;

    std::vector<uchar> heads(r);
    std::vector<ulint> lens(r);
    in.read(reinterpret_cast<char*>(heads.data()), static_cast<std::streamsize>(r * sizeof(uchar)));
    in.read(reinterpret_cast<char*>(lens.data()), static_cast<std::streamsize>(r * sizeof(ulint)));

    std::vector<ulint> C(SIGMA + 1);
    in.read(reinterpret_cast<char*>(C.data()), static_cast<std::streamsize>((SIGMA + 1) * sizeof(ulint)));

    if (!in.good()) return std::nullopt;

    FMDIndex index(heads, lens);
    return index;
}

/**
 * Write a k-mer table to a file.  Converts Position-based BiIntervals to flat
 * (lo, lo_r, size) for storage. Requires fmd for Position->flat conversion.
 */
bool write_kmer_table(const std::string& path, uint64_t fingerprint,
                      const KmerTable& table, const FMDIndex& fmd) {
    std::ofstream out(path, std::ios::binary);
    if (!out.good()) return false;

    out.write(KMER_MAGIC, 4);
    uint32_t v = KMER_VERSION;
    out.write(reinterpret_cast<const char*>(&v), sizeof(v));
    out.write(reinterpret_cast<const char*>(&fingerprint), sizeof(fingerprint));

    uint32_t K = static_cast<uint32_t>(table.K);
    out.write(reinterpret_cast<const char*>(&K), sizeof(K));

    for (const auto& bi : table.table) {
        ulint lo = position_to_flat(fmd, bi.pos);
        ulint lo_r = position_to_flat(fmd, bi.pos_r);
        out.write(reinterpret_cast<const char*>(&lo), sizeof(lo));
        out.write(reinterpret_cast<const char*>(&lo_r), sizeof(lo_r));
        out.write(reinterpret_cast<const char*>(&bi.size), sizeof(bi.size));
    }

    return out.good();
}

/**
 * Read a k-mer table from a file.  Checks fingerprint to ensure it matches
 * expected_fingerprint from the fmd.
 */
std::optional<KmerTable> read_kmer_table(const std::string& path,
                                          uint64_t expected_fingerprint,
                                          const FMDIndex& fmd) {
    std::ifstream in(path, std::ios::binary);
    if (!in.good()) return std::nullopt;

    char magic[5] = {};
    in.read(magic, 4);
    if (std::strncmp(magic, KMER_MAGIC, 4) != 0) return std::nullopt;

    uint32_t v = 0;
    in.read(reinterpret_cast<char*>(&v), sizeof(v));
    if (v != KMER_VERSION) return std::nullopt;

    uint64_t fp = 0;
    in.read(reinterpret_cast<char*>(&fp), sizeof(fp));
    if (fp != expected_fingerprint) return std::nullopt;

    uint32_t K = 0;
    in.read(reinterpret_cast<char*>(&K), sizeof(K));

    uint64_t num_entries = (K == 0) ? 1 : (1ULL << (2 * K));
    KmerTable tab;
    tab.K = static_cast<int>(K);
    tab.table.resize(num_entries);

    for (uint64_t i = 0; i < num_entries; ++i) {
        ulint lo = 0, lo_r = 0, size = 0;
        in.read(reinterpret_cast<char*>(&lo), sizeof(lo));
        in.read(reinterpret_cast<char*>(&lo_r), sizeof(lo_r));
        in.read(reinterpret_cast<char*>(&size), sizeof(size));
        if (!in.good()) return std::nullopt;
        tab.table[i] = {flat_to_position(fmd, lo), flat_to_position(fmd, lo_r), size};
    }

    return tab;
}

/**
 * Check if a file exists.
 */
static bool file_exists(const std::string& path) {
    std::ifstream f(path);
    return f.good();
}

/**
 * Load all available k-mer tables from a base path.  Checks fingerprints on each
 * to ensure they match fmd.
 */
std::vector<KmerTable> load_available_tables(const std::string& base_path,
                                              uint64_t fingerprint,
                                              const FMDIndex& fmd) {
    std::vector<KmerTable> tables;
    for (int K = 0; K <= 20; ++K) {
        std::string path = base_path + ".k" + std::to_string(K);
        if (!file_exists(path)) continue;
        auto opt = read_kmer_table(path, fingerprint, fmd);
        if (opt) tables.push_back(std::move(*opt));
    }
    std::sort(tables.begin(), tables.end(),
              [](const KmerTable& a, const KmerTable& b) { return a.K > b.K; });
    if (tables.empty() || tables.back().K != 0) {
        // If no K=0 table was loaded, synthesize a trivial one
        KmerTable k0;
        k0.K = 0;
        k0.table.push_back({fmd.first_position(), fmd.first_position(), fmd.domain()});
        tables.push_back(std::move(k0));
    }
    return tables;
}

/**
 * Load specific k-mer tables from a base path.  Checks fingerprints on each to
 * ensure they match fmd.
 */
std::vector<KmerTable> load_tables(const std::string& base_path,
                                    uint64_t fingerprint,
                                    const FMDIndex& fmd,
                                    const std::vector<int>& desired_Ks) {
    std::vector<KmerTable> tables;
    for (int K : desired_Ks) {
        std::string path = base_path + ".k" + std::to_string(K);
        auto opt = read_kmer_table(path, fingerprint, fmd);
        if (opt) tables.push_back(std::move(*opt));
    }
    std::sort(tables.begin(), tables.end(),
              [](const KmerTable& a, const KmerTable& b) { return a.K > b.K; });
    bool has_k0 = false;
    for (const auto& t : tables) if (t.K == 0) { has_k0 = true; break; }
    if (!has_k0) {
        // If no K=0 table was loaded, synthesize a trivial one
        KmerTable k0;
        k0.K = 0;
        k0.table.push_back({fmd.first_position(), fmd.first_position(), fmd.domain()});
        tables.push_back(std::move(k0));
    }
    return tables;
}

}  // namespace bi_io
