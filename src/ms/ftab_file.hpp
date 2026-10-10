/**
 * The ftab file: an ftab (ftab.hpp) stored next to the ORBITIDX index it was
 * built from.  See FTAB_FORMAT.md for the byte layout.
 *
 * The first HEADER_BYTES of the file describe the table (k, its fields' bit
 * widths, how k was chosen, its occupancy, its provenance) and link it to
 * its index: copies of the index header's tool, layout, n, BWT runs and
 * format version, and the full section-table entry of each section the
 * table's walk read.  A table is used with an index only if the link
 * matches the index's header, which needs the header alone.  The rest of the
 * index header (its info text, its other sections) does not enter the link,
 * so that an index rewrapped from the same content still matches.
 */

#ifndef _FTAB_FILE_HPP
#define _FTAB_FILE_HPP

#include "ftab.hpp"
#include "index_format.hpp"
#include <cstdint>
#include <cstring>
#include <fstream>
#include <istream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

namespace ftab {

inline constexpr char FILE_MAGIC[8] = {'O', 'R', 'B', 'I', 'F', 'T', 'A', 'B'};
inline constexpr uint32_t FILE_VERSION = 1;
inline constexpr uint32_t HEADER_BYTES = 4096;
inline constexpr size_t FIXED_BYTES = 128;

/** How a table's k was chosen. */
enum class KRule : uint8_t { GIVEN = 1, AUTO = 2 };

/** One payload of the file: where it starts, its bytes and their CRC-32. */
struct Payload {
    uint64_t offset = 0, bytes = 0;
    uint32_t crc = 0;
};

struct FileHeader {
    Shape shape;
    KRule rule = KRule::GIVEN;
    uint64_t present = 0;      // occurring k-mers
    uint64_t built_at = 0;     // seconds since the Unix epoch, UTC
    uint64_t index_bytes = 0;  // the index file's size when the table was built
    // The link: copies of the index header's fixed fields and section entries.
    orbit_idx::Tool tool{};
    orbit_idx::Layout layout = orbit_idx::Layout::NONE;
    uint64_t n = 0, bwt_runs = 0;
    uint32_t index_version = 0;
    std::vector<orbit_idx::SectionEntry> link;
    Payload table, window;
    std::string info;  // provenance, "key=value" lines

    double occupancy() const { return double(present) / double(uint64_t{1} << (2 * shape.k)); }
};

/**
 * The sections a table's walk reads, among those index header h holds: LF
 * and MS_SPILL for ms; LF and, with psi, F_STARTS for tms, plus both forms
 * of phi for a table with a toehold, whose phi points either form reads.
 */
inline std::vector<orbit_idx::Section> link_sections(const orbit_idx::Header& h, bool toehold) {
    using orbit_idx::Section;
    std::vector<Section> want;
    if (h.tool == orbit_idx::Tool::MS) want = {Section::LF, Section::MS_SPILL};
    else {
        want = {Section::LF, Section::F_STARTS};
        if (toehold) want.insert(want.end(), {Section::PHI_LEN, Section::PHI_START});
    }
    std::vector<Section> out;
    for (Section s : want)
        if (h.find(s)) out.push_back(s);
    return out;
}

/** A header for table tab built from the index whose header is ih, with its link filled in. */
inline FileHeader make_header(const orbit_idx::Header& ih, const Ftab& tab) {
    FileHeader h;
    h.shape = tab.shape();
    h.tool = ih.tool;
    h.layout = ih.layout;
    h.n = ih.n;
    h.bwt_runs = ih.bwt_runs;
    h.index_version = orbit_idx::FORMAT_VERSION;
    for (auto s : link_sections(ih, tab.has_toehold())) {
        orbit_idx::SectionEntry e = *ih.find(s);
        e.offset = 0;  // where a section lies is not part of its content
        h.link.push_back(e);
    }
    return h;
}

/** The first HEADER_BYTES bytes of an index file, decoded. */
inline orbit_idx::Header index_header_of(std::istream& in) {
    std::vector<char> b(orbit_idx::HEADER_BYTES);
    in.read(b.data(), orbit_idx::HEADER_BYTES);
    if (in.gcount() != static_cast<std::streamsize>(orbit_idx::HEADER_BYTES))
        throw std::runtime_error("truncated Orbit index header");
    return orbit_idx::decode_header(b.data());
}

namespace detail {
using orbit_idx::detail::get;
using orbit_idx::detail::put;

inline std::vector<char> encode(const FileHeader& h) {
    if (FIXED_BYTES + h.link.size() * orbit_idx::SECTION_ENTRY_BYTES + h.info.size() > HEADER_BYTES)
        throw std::length_error("ftab header does not fit in its block");
    std::vector<char> b(HEADER_BYTES, 0);
    const Shape& s = h.shape;
    std::memcpy(b.data(), FILE_MAGIC, 8);
    put<uint32_t>(b, 8, FILE_VERSION);
    put<uint32_t>(b, 12, orbit_idx::ENDIAN_CHECK);
    put<uint32_t>(b, 16, HEADER_BYTES);
    put<uint8_t>(b, 20, static_cast<uint8_t>(s.k));
    put<uint8_t>(b, 21, static_cast<uint8_t>(h.rule));
    put<uint8_t>(b, 22, static_cast<uint8_t>((s.toehold ? 1 : 0) | (s.window ? 2 : 0)));
    put<uint8_t>(b, 23, static_cast<uint8_t>(s.entry_bits()));
    put<uint8_t>(b, 24, static_cast<uint8_t>(s.interval));
    put<uint8_t>(b, 25, static_cast<uint8_t>(s.offset));
    put<uint8_t>(b, 26, static_cast<uint8_t>(s.toehold_interval));
    put<uint8_t>(b, 27, static_cast<uint8_t>(s.toehold_offset));
    put<uint8_t>(b, 28, static_cast<uint8_t>(s.position));
    put<uint16_t>(b, 30, static_cast<uint16_t>(h.link.size()));
    put<uint16_t>(b, 32, static_cast<uint16_t>(h.info.size()));
    put<uint64_t>(b, 40, h.present);
    put<uint64_t>(b, 48, h.built_at);
    put<uint64_t>(b, 56, h.index_bytes);
    put<uint16_t>(b, 64, static_cast<uint16_t>(h.tool));
    put<uint16_t>(b, 66, static_cast<uint16_t>(h.layout));
    put<uint32_t>(b, 68, h.index_version);
    put<uint64_t>(b, 72, h.n);
    put<uint64_t>(b, 80, h.bwt_runs);
    put<uint64_t>(b, 88, h.table.offset);
    put<uint64_t>(b, 96, h.table.bytes);
    put<uint64_t>(b, 104, h.window.offset);
    put<uint64_t>(b, 112, h.window.bytes);
    put<uint32_t>(b, 120, h.table.crc);
    put<uint32_t>(b, 124, h.window.crc);
    size_t at = FIXED_BYTES;
    for (const auto& e : h.link) {
        orbit_idx::encode_entry(b, at, e);
        at += orbit_idx::SECTION_ENTRY_BYTES;
    }
    std::memcpy(b.data() + at, h.info.data(), h.info.size());
    return b;
}

inline FileHeader decode(const char* b) {
    if (std::memcmp(b, FILE_MAGIC, 8) != 0) throw std::runtime_error("not an ftab file (no ORBIFTAB magic)");
    if (get<uint32_t>(b, 12) != orbit_idx::ENDIAN_CHECK) throw std::runtime_error("ftab written with another byte order");
    const uint32_t v = get<uint32_t>(b, 8);
    if (v != FILE_VERSION) throw std::runtime_error("unsupported ftab format version " + std::to_string(v));
    if (get<uint32_t>(b, 16) != HEADER_BYTES) throw std::runtime_error("unsupported ftab header size");
    FileHeader h;
    Shape& s = h.shape;
    s.k = get<uint8_t>(b, 20);
    h.rule = static_cast<KRule>(get<uint8_t>(b, 21));
    const uint8_t flags = get<uint8_t>(b, 22);
    s.toehold = flags & 1;
    s.window = flags & 2;
    s.interval = get<uint8_t>(b, 24);
    s.offset = get<uint8_t>(b, 25);
    s.toehold_interval = get<uint8_t>(b, 26);
    s.toehold_offset = get<uint8_t>(b, 27);
    s.position = get<uint8_t>(b, 28);
    if (s.k < 1 || s.k > max_k || get<uint8_t>(b, 23) != s.entry_bits() || (flags & ~3) != 0 ||
        (h.rule != KRule::GIVEN && h.rule != KRule::AUTO))
        throw std::runtime_error("malformed ftab header: bad k, flags, rule or entry bits");
    bool reserved_zero = b[29] == 0;
    for (size_t i = 34; i < 40; ++i) reserved_zero = reserved_zero && b[i] == 0;
    if (!reserved_zero) throw std::runtime_error("malformed ftab header: reserved bytes are not 0");
    const size_t nl = get<uint16_t>(b, 30), ni = get<uint16_t>(b, 32);
    h.present = get<uint64_t>(b, 40);
    h.built_at = get<uint64_t>(b, 48);
    h.index_bytes = get<uint64_t>(b, 56);
    h.tool = static_cast<orbit_idx::Tool>(get<uint16_t>(b, 64));
    h.layout = static_cast<orbit_idx::Layout>(get<uint16_t>(b, 66));
    h.index_version = get<uint32_t>(b, 68);
    h.n = get<uint64_t>(b, 72);
    h.bwt_runs = get<uint64_t>(b, 80);
    h.table = {get<uint64_t>(b, 88), get<uint64_t>(b, 96), get<uint32_t>(b, 120)};
    h.window = {get<uint64_t>(b, 104), get<uint64_t>(b, 112), get<uint32_t>(b, 124)};
    if (FIXED_BYTES + nl * orbit_idx::SECTION_ENTRY_BYTES + ni > HEADER_BYTES)
        throw std::runtime_error("malformed ftab header: its tables overrun it");
    if (h.table.bytes != 8 * s.table_words() || h.window.bytes != 8 * s.window_words())
        throw std::runtime_error("malformed ftab header: payload sizes disagree with the shape");
    size_t at = FIXED_BYTES;
    for (size_t i = 0; i < nl; ++i, at += orbit_idx::SECTION_ENTRY_BYTES)
        h.link.push_back(orbit_idx::decode_entry(b, at, "ftab"));
    h.info.assign(b + at, ni);
    return h;
}

inline void write_zeros(std::ostream& out, uint64_t n) {
    static const char zeros[4096] = {};
    while (n > 0) {
        const uint64_t k = std::min<uint64_t>(n, sizeof(zeros));
        out.write(zeros, static_cast<std::streamsize>(k));
        n -= k;
    }
}

inline uint32_t crc_of(const std::vector<uint64_t>& w) { return orbit_idx::crc32(0, w.data(), w.size() * 8); }
}  // namespace detail

/**
 * Write table tab with header h (whose shape, payloads and CRCs are filled
 * in here): the header, then the entries from byte HEADER_BYTES, then the
 * window positions, if any, from the next multiple of 4096.
 */
inline void write(std::ostream& out, FileHeader h, const Ftab& tab) {
    h.shape = tab.shape();
    h.table = {HEADER_BYTES, tab.table_words().size() * 8, detail::crc_of(tab.table_words())};
    h.window = {0, tab.window_words().size() * 8, 0};
    if (h.window.bytes > 0) {
        h.window.offset = orbit_idx::align_up(h.table.offset + h.table.bytes);
        h.window.crc = detail::crc_of(tab.window_words());
    }
    const std::vector<char> hb = detail::encode(h);
    out.write(hb.data(), static_cast<std::streamsize>(hb.size()));
    out.write(reinterpret_cast<const char*>(tab.table_words().data()), static_cast<std::streamsize>(h.table.bytes));
    if (h.window.bytes > 0) {
        detail::write_zeros(out, h.window.offset - (h.table.offset + h.table.bytes));
        out.write(reinterpret_cast<const char*>(tab.window_words().data()), static_cast<std::streamsize>(h.window.bytes));
    }
    if (!out.good()) throw std::runtime_error("failed writing ftab");
}

inline void write_file(const std::string& path, const FileHeader& h, const Ftab& tab) {
    std::ofstream out(path, std::ios::binary);
    if (!out.good()) throw std::runtime_error("cannot open " + path + " for writing");
    write(out, h, tab);
    out.close();
    if (!out) throw std::runtime_error("failed writing " + path);
}

/** The header at the start of in. */
inline FileHeader read_header(std::istream& in) {
    std::vector<char> b(HEADER_BYTES);
    in.read(b.data(), HEADER_BYTES);
    if (in.gcount() < 8 || std::memcmp(b.data(), FILE_MAGIC, 8) != 0)
        throw std::runtime_error("not an ftab file (no ORBIFTAB magic)");
    if (in.gcount() != HEADER_BYTES) throw std::runtime_error("truncated ftab header");
    return detail::decode(b.data());
}

inline FileHeader read_header_file(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in.good()) throw std::runtime_error("cannot open " + path);
    return read_header(in);
}

/**
 * Read the table whose header h has just been read from in.  With check_crc,
 * also check each payload's CRC-32.  Throws if the file is short or a CRC
 * differs.
 */
inline Ftab read_payload(std::istream& in, const FileHeader& h, bool check_crc) {
    Ftab tab(h.shape);
    auto read_at = [&](const Payload& p, std::vector<uint64_t>& w, uint64_t at, const char* what) {
        if (p.offset < at) throw std::runtime_error("malformed ftab: payloads out of order");
        in.ignore(static_cast<std::streamsize>(p.offset - at));
        in.read(reinterpret_cast<char*>(w.data()), static_cast<std::streamsize>(p.bytes));
        if (!in.good()) throw std::runtime_error(std::string("ftab ") + what + " is truncated");
        if (check_crc && detail::crc_of(w) != p.crc) throw std::runtime_error(std::string("ftab ") + what + " fails its CRC-32");
        return p.offset + p.bytes;
    };
    uint64_t at = read_at(h.table, tab.table_words(), HEADER_BYTES, "table");
    if (h.window.bytes > 0) read_at(h.window, tab.window_words(), at, "window positions");
    return tab;
}

/** Read the table file at path, with its header in h. */
inline Ftab read_file(const std::string& path, FileHeader& h, bool check_crc = false) {
    std::ifstream in(path, std::ios::binary);
    if (!in.good()) throw std::runtime_error("cannot open " + path);
    h = read_header(in);
    return read_payload(in, h, check_crc);
}

/**
 * How the link of table header f disagrees with index header ih, naming
 * the first field that differs, or an empty string if it matches.
 */
inline std::string link_mismatch(const FileHeader& f, const orbit_idx::Header& ih) {
    std::ostringstream os;
    auto differs = [&](const char* field, auto want, auto got) {
        if (want == got) return false;
        os << field << ": the ftab was built from " << want << ", the index has " << got;
        return true;
    };
    auto tool = [](orbit_idx::Tool t) { return t == orbit_idx::Tool::MS ? "ms" : t == orbit_idx::Tool::TMS ? "tms" : "?"; };
    if (differs("tool", std::string(tool(f.tool)), std::string(tool(ih.tool))) || differs("layout", std::string(orbit_idx::layout_name(f.layout)), std::string(orbit_idx::layout_name(ih.layout))) ||
        differs("n", f.n, ih.n) || differs("BWT runs", f.bwt_runs, ih.bwt_runs) ||
        differs("index format version", f.index_version, orbit_idx::FORMAT_VERSION))
        return os.str();
    for (const auto& e : f.link) {
        const orbit_idx::SectionEntry* g = ih.find(e.id);
        const std::string name = std::string("section ") + orbit_idx::section_name(e.id) + " ";
        if (!g) return name + "is missing from the index";
        if (differs((name + "kind").c_str(), int(e.kind), int(g->kind)) ||
            differs((name + "rows").c_str(), e.rows, g->rows) || differs((name + "domain").c_str(), e.domain, g->domain) ||
            differs((name + "payload bytes").c_str(), e.bytes, g->bytes) ||
            differs((name + "row bits").c_str(), e.row_bits, g->row_bits) ||
            differs((name + "columns").c_str(), e.cols.size(), g->cols.size()))
            return os.str();
        for (size_t c = 0; c < e.cols.size(); ++c)
            if (differs((name + "column " + std::to_string(c) + " code").c_str(), int(e.cols[c].first), int(g->cols[c].first)) ||
                differs((name + "column " + std::to_string(c) + " width").c_str(), int(e.cols[c].second), int(g->cols[c].second)))
                return os.str();
        if (e.crc != g->crc) {
            os << name << "CRC-32: the ftab was built from " << std::hex << e.crc << ", the index has " << g->crc;
            return os.str();
        }
    }
    // A section the walk would read that the table's index lacked.
    for (auto id : link_sections(ih, f.shape.toehold)) {
        bool linked = false;
        for (const auto& e : f.link) linked = linked || e.id == id;
        if (!linked) return std::string("section ") + orbit_idx::section_name(id) + " is in the index but not in the ftab's link";
    }
    return "";
}

/** A description of h: its shape, sizes, link and info. */
inline void describe(const FileHeader& h, std::ostream& os) {
    const Shape& s = h.shape;
    os << "ftab format version " << FILE_VERSION << "\n";
    os << "k=" << s.k << " (" << (h.rule == KRule::AUTO ? "auto" : "given") << ")"
       << " present=" << h.present << " of " << (uint64_t{1} << (2 * s.k)) << " occupancy=" << h.occupancy() << "\n";
    os << "entry: " << s.entry_bits() << " bits = present 1 + interval " << s.interval << " + offset " << s.offset;
    if (s.toehold) os << " + toehold interval " << s.toehold_interval << " + toehold offset " << s.toehold_offset;
    os << (s.toehold ? "" : " (no toehold)") << "\n";
    if (s.window) os << "window positions: " << s.window_entries() << " of " << s.position << " bits\n";
    os << "payloads: table " << h.table.bytes << " bytes at " << h.table.offset << " (crc " << std::hex << h.table.crc
       << std::dec << ")";
    if (h.window.bytes > 0)
        os << ", window " << h.window.bytes << " bytes at " << h.window.offset << " (crc " << std::hex << h.window.crc << std::dec << ")";
    os << "; total " << s.payload_bytes() << " bytes";
    if (h.index_bytes > 0) os << " = " << 100.0 * double(s.payload_bytes()) / double(h.index_bytes) << "% of the index's " << h.index_bytes;
    os << "\n";
    os << "index: tool=" << (h.tool == orbit_idx::Tool::MS ? "ms" : "tms") << " layout=" << orbit_idx::layout_name(h.layout)
       << " n=" << h.n << " bwt_runs=" << h.bwt_runs << " format_version=" << h.index_version << "\n";
    for (const auto& e : h.link) {
        os << "  link " << orbit_idx::section_name(e.id) << ": rows=" << e.rows << " domain=" << e.domain
           << " bytes=" << e.bytes << " crc=" << std::hex << e.crc << std::dec;
        if (!e.cols.empty()) {
            os << " cols=";
            for (size_t c = 0; c < e.cols.size(); ++c)
                os << (c ? "," : "") << orbit_idx::col_name(e.cols[c].first) << ":" << int(e.cols[c].second);
        }
        os << "\n";
    }
    if (!h.info.empty()) os << h.info << (h.info.back() == '\n' ? "" : "\n");
}

}  // namespace ftab

#endif /* _FTAB_FILE_HPP */
