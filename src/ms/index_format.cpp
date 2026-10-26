/**
 * ORBITIDX container: header description, CRC verification, and the rewrap
 * of legacy MSIX v4 and TMSX v8 indexes (see index_format.hpp).
 */

#include "index_format.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <memory>

namespace orbit_idx {

Header read_header_file(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in.good()) throw std::runtime_error("cannot open " + path);
    std::vector<char> b(HEADER_BYTES);
    in.read(b.data(), HEADER_BYTES);
    if (in.gcount() < 8 || !has_magic(b.data(), 8)) throw std::runtime_error(path + " is not an Orbit index (no ORBITIDX magic)");
    if (!in.good()) throw std::runtime_error("truncated Orbit index header in " + path);
    return decode_header(b.data());
}

namespace {

std::string mb(uint64_t bytes) {
    std::ostringstream os;
    os << std::fixed << std::setprecision(bytes < 10'000'000 ? 3 : 1) << bytes / 1e6 << " MB";
    return os.str();
}

const char* kind_name(Kind k) {
    switch (k) {
        case Kind::MOVE_TABLE: return "move_table";
        case Kind::U64_ARRAY: return "u64_array";
        case Kind::SPILL_ARRAYS: return "spill_arrays";
        case Kind::SEQ_TABLE: return "seq_table";
    }
    return "?";
}

}  // namespace

void describe_header(const Header& h, std::ostream& os) {
    os << "format: ORBITIDX v" << FORMAT_VERSION << "  tool: " << (h.tool == Tool::MS ? "ms" : "tms")
       << "  layout: " << layout_name(h.layout) << "  n=" << h.n << "  bwt_runs=" << h.bwt_runs << "\n";
    if (h.tool == Tool::MS)
        os << "scalars: max_lcp_top=" << h.scalars[0] << " max_lcp_min_sub=" << h.scalars[1]
           << " spill_align=" << h.scalars[2] << " spill_split_bits=" << h.scalars[3] << "\n";
    os << "sections:\n";
    for (const auto& e : h.sections) {
        os << "  " << std::left << std::setw(10) << section_name(e.id) << std::right << " " << kind_name(e.kind)
           << " offset=" << e.offset << " bytes=" << e.bytes << " crc=" << std::hex << std::setw(8)
           << std::setfill('0') << e.crc << std::dec << std::setfill(' ');
        if (e.kind == Kind::MOVE_TABLE) {
            os << " rows=" << e.rows << " domain=" << e.domain << " row_bits=" << e.row_bits << " cols=";
            for (size_t c = 0; c < e.cols.size(); ++c)
                os << (c ? "," : "") << col_name(e.cols[c].first) << ":" << int(e.cols[c].second);
        } else {
            os << " " << (e.kind == Kind::U64_ARRAY ? "values=" : e.kind == Kind::SEQ_TABLE ? "sequences=" : "arrays=")
               << e.rows;
        }
        os << "\n";
    }
    if (!h.info.empty()) {
        os << "info:\n";
        std::istringstream is(h.info);
        for (std::string line; std::getline(is, line);) os << "  " << line << "\n";
    }
}

std::string verify_file(const std::string& path) {
    try {
        const Header h = read_header_file(path);
        std::ifstream in(path, std::ios::binary);
        std::vector<char> buf(1 << 22);
        std::string bad;
        for (const auto& e : h.sections) {
            in.seekg(static_cast<std::streamoff>(e.offset));
            uint32_t crc = 0;
            for (uint64_t left = e.bytes; left > 0;) {
                const auto k = static_cast<std::streamsize>(std::min<uint64_t>(left, buf.size()));
                in.read(buf.data(), k);
                if (in.gcount() != k) return std::string("section ") + section_name(e.id) + " is truncated";
                crc = crc32(crc, buf.data(), static_cast<size_t>(k));
                left -= static_cast<uint64_t>(k);
            }
            if (crc != e.crc) bad += std::string(bad.empty() ? "" : "; ") + "section " + section_name(e.id) + " CRC mismatch";
        }
        const uint64_t end = h.sections.empty() ? HEADER_BYTES : h.sections.back().offset + h.sections.back().bytes;
        const uint64_t size = std::filesystem::file_size(path);
        if (size != end)
            bad += std::string(bad.empty() ? "" : "; ") + "file is " + std::to_string(size) + " bytes, sections end at " +
                   std::to_string(end);
        return bad;
    } catch (const std::exception& e) {
        return e.what();
    }
}

// --- Rewrap of legacy indexes ---------------------------------------------

namespace {

/** Reads byte ranges of a legacy file. */
struct LegacyFile {
    std::ifstream in;
    uint64_t size = 0;
    explicit LegacyFile(const std::string& path) : in(path, std::ios::binary) {
        if (!in.good()) throw std::runtime_error("cannot open " + path);
        size = std::filesystem::file_size(path);
    }
    void read(uint64_t at, void* p, size_t n) {
        if (at + n > size) throw std::runtime_error("legacy index is truncated");
        in.clear();
        in.seekg(static_cast<std::streamoff>(at));
        in.read(static_cast<char*>(p), static_cast<std::streamsize>(n));
        if (!in.good()) throw std::runtime_error("failed reading legacy index");
    }
    template <typename T>
    T get(uint64_t at) {
        T v;
        read(at, &v, sizeof(v));
        return v;
    }
};

/** An Orbit permutation's serialization in a legacy file, as its preamble describes it. */
struct OrbitPerm {
    uint64_t begin = 0, data = 0, end = 0;  // its first byte, its packed rows' first byte, and one past it
    uint64_t n = 0, rows = 0;
    std::vector<uint8_t> widths;
};

/**
 * Parse the preamble of the Orbit permutation with ncols columns at at:
 * "OBPT", a version, split parameters, "OBMS", a version, n and r, "OBPM", a
 * version, the row count and one width per column, then the packed rows,
 * ceil(rows * row_bits / 8) + 8 bytes.  The alphabet, of ms_alpha's fixed
 * kind, adds no bytes.
 */
OrbitPerm parse_perm(LegacyFile& f, uint64_t at, size_t ncols) {
    auto magic = [&](uint64_t p, const char* m) {
        char b[4];
        f.read(p, b, 4);
        if (std::memcmp(b, m, 4) != 0) throw std::runtime_error(std::string("legacy index lacks the ") + m + " magic where expected");
    };
    OrbitPerm o;
    o.begin = at;
    magic(at, "OBPT");
    uint64_t p = at + 4 + 24;
    for (int k = 0; k < 2; ++k) p += 1 + (f.get<uint8_t>(p) ? 8 : 0);  // length capping, then balancing
    magic(p, "OBMS");
    p += 4 + 24;
    o.n = f.get<uint64_t>(p);
    p += 16;
    magic(p, "OBPM");
    p += 4 + 24;
    o.rows = f.get<uint64_t>(p);
    p += 8;
    o.widths.resize(ncols);
    f.read(p, o.widths.data(), ncols);
    p += ncols;
    uint64_t bits = 0;
    for (auto w : o.widths) bits += w;
    o.data = p;
    o.end = p + (o.rows * bits + 7) / 8 + 8;
    if (o.end > f.size) throw std::runtime_error("legacy index is truncated inside a move table");
    return o;
}

/**
 * Copy bytes [from, to) of f to out.  With runs, also count the maximal runs
 * of one value in column col of the move table perm (whose rows lie inside
 * the range) into *runs.
 */
void copy_range(LegacyFile& f, uint64_t from, uint64_t to, std::ostream& out, const OrbitPerm* perm = nullptr,
                size_t col = 0, uint64_t* runs = nullptr) {
    std::vector<char> buf(1 << 23);
    // Decoding: rows are counted in groups of eight, which span row_bits
    // whole bytes, held in carry until a group is complete.
    unsigned row_bits = 0, col_off = 0, col_w = 0;
    if (perm && runs) {
        for (size_t c = 0; c < perm->widths.size(); ++c) {
            if (c < col) col_off += perm->widths[c];
            row_bits += perm->widths[c];
        }
        col_w = perm->widths[col];
        *runs = 0;
    }
    std::vector<unsigned char> carry;
    uint64_t row = 0, prev = ~uint64_t(0);
    auto decode = [&](const char* p, size_t n, uint64_t at) {
        // The part of [at, at + n) inside the rows' bytes.
        const uint64_t lo = std::max(at, perm->data), hi = std::min(at + n, perm->end - 8);
        if (lo >= hi) return;
        carry.insert(carry.end(), p + (lo - at), p + (hi - at));
        const size_t group = row_bits;  // bytes per eight rows
        size_t used = 0;
        auto take = [&](const unsigned char* g, size_t rows_in) {
            for (size_t j = 0; j < rows_in && row < perm->rows; ++j, ++row) {
                const uint64_t bit = j * row_bits + col_off;
                uint64_t w = 0;
                const size_t avail = carry.size() - used;
                if (bit / 8 >= avail) throw std::runtime_error("legacy move table has fewer row bytes than rows");
                std::memcpy(&w, g + bit / 8, std::min<size_t>(8, avail - bit / 8));
                const uint64_t v = (w >> (bit % 8)) & ((uint64_t(1) << col_w) - 1);
                *runs += v != prev;
                prev = v;
            }
        };
        while (carry.size() - used >= group + 8) {
            take(carry.data() + used, 8);
            used += group;
        }
        // The last rows, once every row byte is in; a group may be short.
        if (hi == perm->end - 8)
            while (row < perm->rows) {
                if (used >= carry.size()) throw std::runtime_error("legacy move table has fewer row bytes than rows");
                take(carry.data() + used, 8);
                used = std::min(used + group, carry.size());
            }
        carry.erase(carry.begin(), carry.begin() + static_cast<std::ptrdiff_t>(used));
    };
    for (uint64_t at = from; at < to;) {
        const size_t k = static_cast<size_t>(std::min<uint64_t>(to - at, buf.size()));
        f.read(at, buf.data(), k);
        out.write(buf.data(), static_cast<std::streamsize>(k));
        if (perm && runs && row_bits > 0 && col_w > 0) decode(buf.data(), k, at);
        at += k;
    }
}

template <size_t N>
SectionEntry legacy_move_entry(Section id, const OrbitPerm& o, const std::array<Col, N>& codes) {
    SectionEntry e;
    e.id = id;
    e.kind = Kind::MOVE_TABLE;
    e.rows = o.rows;
    e.domain = o.n;
    for (size_t c = 0; c < N; ++c) {
        e.cols.emplace_back(static_cast<uint8_t>(codes[c]), o.widths[c]);
        e.row_bits = static_cast<uint16_t>(e.row_bits + o.widths[c]);
    }
    return e;
}

// The column codes of each legacy structure, in row order.
constexpr std::array<Col, 8> TMS_LF_COLS{Col::LENGTH, Col::POINTER, Col::OFFSET, Col::CHARACTER,
                                         Col::PSI_INT, Col::PSI_OFF, Col::PHI_INT, Col::PHI_OFF};
constexpr std::array<Col, 4> TMS_PHI_LEN_COLS{Col::LENGTH, Col::POINTER, Col::OFFSET, Col::PLCP};
constexpr std::array<Col, 4> TMS_PHI_COLS{Col::START, Col::POINTER, Col::OFFSET, Col::PLCP};
constexpr std::array<Col, 4> TMS_PHI_INV_COLS{Col::START, Col::POINTER, Col::OFFSET, Col::PLCPB};
constexpr std::array<Col, 7> MS_LF_COLS{Col::LENGTH, Col::POINTER, Col::OFFSET, Col::CHARACTER,
                                        Col::LCP_TOP, Col::LCP_MIN_SUB, Col::LCP_SPILL};

}  // namespace

std::string rewrap_legacy(const std::string& from, const std::string& to, const std::string& info) {
    try {
        LegacyFile f(from);
        char magic[4];
        f.read(0, magic, 4);
        const uint32_t version = f.get<uint32_t>(4);
        Header h;
        std::vector<SectionSource> srcs;
        uint64_t runs = 0;
        OrbitPerm lf;
        if (std::memcmp(magic, "MSIX", 4) == 0) {
            if (version != 4) return "MSIX version " + std::to_string(version) + " is not rewrapped (only 4)";
            h.tool = Tool::MS;
            lf = parse_perm(f, 8, MS_LF_COLS.size());
            uint64_t p = lf.end;
            for (size_t i = 0; i < SCALARS; ++i, p += 8) h.scalars[i] = f.get<uint64_t>(p);
            const uint64_t arrays = f.get<uint64_t>(p);
            p += 8;
            const uint64_t spill_begin = p;
            for (uint64_t i = 0; i < arrays; ++i) p += 8 + f.get<uint64_t>(p);
            if (p != f.size) return "MSIX file has bytes after its spillover arrays";
            srcs.push_back({legacy_move_entry(Section::LF, lf, MS_LF_COLS), {}});
            SectionEntry sp;
            sp.id = Section::MS_SPILL;
            sp.kind = Kind::SPILL_ARRAYS;
            sp.rows = arrays;
            srcs.push_back({sp, [&f, spill_begin](std::ostream& o) { copy_range(f, spill_begin, f.size, o); }});
        } else if (std::memcmp(magic, "TMSX", 4) == 0) {
            if (version != 8) return "TMSX version " + std::to_string(version) + " is not rewrapped (only 8)";
            h.tool = Tool::TMS;
            const uint32_t flags = f.get<uint32_t>(8);
            const bool phi = flags & 1, phi_inv = flags & 2, psi = flags & 4;
            h.layout = psi && phi ? Layout::FULL : psi ? Layout::PSI : phi ? Layout::PHI : Layout::NONE;
            lf = parse_perm(f, 12, TMS_LF_COLS.size());
            srcs.push_back({legacy_move_entry(Section::LF, lf, TMS_LF_COLS), {}});
            uint64_t p = lf.end;
            if (psi) {
                SectionEntry e;
                e.id = Section::F_STARTS;
                e.kind = Kind::U64_ARRAY;
                e.rows = 7;
                const uint64_t at = p;
                srcs.push_back({e, [&f, at](std::ostream& o) { copy_range(f, at, at + 56, o); }});
                p += 56;
            }
            // A structure written after its size.
            auto sized = [&](Section id, const auto& codes) {
                const uint64_t size = f.get<uint64_t>(p);
                const OrbitPerm o = parse_perm(f, p + 8, codes.size());
                if (o.end != p + 8 + size) throw std::runtime_error("legacy structure's size disagrees with its preamble");
                srcs.push_back({legacy_move_entry(id, o, codes),
                                [&f, o](std::ostream& out) { copy_range(f, o.begin, o.end, out); }});
                p = o.end;
            };
            if (phi) {
                sized(Section::PHI_LEN, TMS_PHI_LEN_COLS);
                sized(Section::PHI_START, TMS_PHI_COLS);
            }
            if (phi_inv) sized(Section::PHI_INV, TMS_PHI_INV_COLS);
            if (p != f.size) return "TMSX file has bytes after its last structure";
        } else {
            return from + " is neither an MSIX nor a TMSX index";
        }
        h.n = lf.n;
        // LF's payload; its first writing also counts the BWT runs in its character column.
        bool counted = false;
        srcs[0].write = [&f, &lf, &runs, &counted](std::ostream& o) {
            copy_range(f, lf.begin, lf.end, o, &lf, 3, counted ? nullptr : &runs);
            counted = true;
        };
        h.info = info;
        std::ofstream out(to, std::ios::binary);
        if (!out.good()) return "cannot write " + to;
        write_file(out, h, srcs, [&](Header& hh) { hh.bwt_runs = runs; });
        out.close();
        if (!out) return "failed writing " + to;
        return "";
    } catch (const std::exception& e) {
        return e.what();
    }
}

namespace {
/**
 * Lay out sections after the header: set each entry's offset, and return the
 * offsets they had before (where they are read from).
 */
std::vector<uint64_t> lay_out(std::vector<SectionEntry>& sections) {
    std::vector<uint64_t> from_offsets;
    uint64_t at = HEADER_BYTES;
    for (auto& e : sections) {
        from_offsets.push_back(e.offset);
        e.offset = align_up(at);
        at = e.offset + e.bytes;
    }
    return from_offsets;
}

/**
 * Copy section e, at from_offset in in (the file from), to out, checking it
 * against its CRC-32.  buf is the copy buffer.
 */
void copy_section(std::ifstream& in, const std::string& from, std::ostream& out, const SectionEntry& e,
                  uint64_t from_offset, std::vector<char>& buf) {
    in.seekg(static_cast<std::streamoff>(from_offset));
    uint32_t crc = 0;
    for (uint64_t left = e.bytes; left > 0;) {
        const size_t k = static_cast<size_t>(std::min<uint64_t>(left, buf.size()));
        in.read(buf.data(), static_cast<std::streamsize>(k));
        if (!in.good()) throw std::runtime_error(from + ": section " + section_name(e.id) + " is truncated");
        crc = crc32(crc, buf.data(), k);
        out.write(buf.data(), static_cast<std::streamsize>(k));
        left -= k;
    }
    if (crc != e.crc) throw std::runtime_error(from + ": section " + section_name(e.id) + " fails its CRC-32");
}

/** Open from for reading and check that to is another file. */
std::ifstream open_copy(const std::string& from, const std::string& to) {
    if (std::filesystem::exists(to) && std::filesystem::equivalent(from, to))
        throw std::runtime_error("the output must be another file than the input");
    std::ifstream in(from, std::ios::binary);
    if (!in.good()) throw std::runtime_error("cannot open " + from);
    return in;
}

/** Open to and write the header of h with sections and info appended. */
std::ofstream write_copy_header(const std::string& to, Header h, const std::vector<SectionEntry>& sections,
                                const std::string& info) {
    h.sections = sections;
    h.info += info;
    const std::vector<char> hb = encode_header(h);
    std::ofstream out(to, std::ios::binary | std::ios::trunc);
    if (!out.good()) throw std::runtime_error("cannot open " + to + " for writing");
    out.write(hb.data(), static_cast<std::streamsize>(hb.size()));
    return out;
}
}  // namespace

void add_section(const std::string& from, const std::string& to, SectionSource& s, const std::string& info) {
    std::ifstream in = open_copy(from, to);
    Header h = read_header_file(from);
    // Size and checksum the new payload, then lay out every section, the
    // new one where a section it replaces was, else last.
    detail::CrcCountBuf sink;
    {
        std::ostream os(&sink);
        s.write(os);
    }
    SectionEntry added = s.entry;
    added.bytes = sink.count;
    added.crc = sink.crc;
    std::vector<SectionEntry> sections = h.sections;
    size_t slot = sections.size();
    for (size_t i = 0; i < sections.size(); ++i)
        if (sections[i].id == added.id) slot = i;
    if (slot == sections.size()) sections.push_back(added);
    else sections[slot] = added;
    const std::vector<uint64_t> from_offsets = lay_out(sections);
    std::ofstream out = write_copy_header(to, h, sections, info);
    uint64_t pos = HEADER_BYTES;
    std::vector<char> buf(size_t(1) << 22);
    for (size_t i = 0; i < sections.size(); ++i) {
        const SectionEntry& e = sections[i];
        detail::write_zeros(out, e.offset - pos);
        pos = e.offset + e.bytes;
        if (i != slot) {
            copy_section(in, from, out, e, from_offsets[i], buf);
            continue;
        }
        detail::CrcCountBuf tee(out.rdbuf());
        {
            std::ostream os(&tee);
            s.write(os);
            os.flush();
        }
        if (tee.count != added.bytes || tee.crc != added.crc)
            throw std::logic_error(std::string("section ") + section_name(added.id) + " serialized differently twice");
    }
    out.flush();
    if (!out.good()) throw std::runtime_error("failed writing " + to);
}

void drop_sections(const std::string& from, const std::string& to, uint32_t drop, Layout layout,
                   const std::string& info) {
    std::ifstream in = open_copy(from, to);
    Header h = read_header_file(from);
    std::vector<SectionEntry> kept;
    for (const auto& e : h.sections)
        if (!(drop & section_bit(e.id))) kept.push_back(e);
    const std::vector<uint64_t> from_offsets = lay_out(kept);
    h.layout = layout;
    std::ofstream out = write_copy_header(to, h, kept, info);
    uint64_t pos = HEADER_BYTES;
    std::vector<char> buf(size_t(1) << 22);
    for (size_t i = 0; i < kept.size(); ++i) {
        detail::write_zeros(out, kept[i].offset - pos);
        copy_section(in, from, out, kept[i], from_offsets[i], buf);
        pos = kept[i].offset + kept[i].bytes;
    }
    out.flush();
    if (!out.good()) throw std::runtime_error("failed writing " + to);
}

}  // namespace orbit_idx
