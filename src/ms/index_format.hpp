/**
 * The ORBITIDX container that Orbit-ms and Orbit-tms indexes are written in.
 * See INDEX_FORMAT.md for the byte layout.
 *
 * The first HEADER_BYTES of a file describe what it holds: the tool and
 * layout, the text length, and a table of sections, each with its offset,
 * size and checksum and, for a move table, its rows, domain and columns (a
 * column code and a bit width each, in the order the rows hold them).  Which
 * sections a query loads is the software's choice, not the file's; with the
 * section table, the memory that choice takes follows from the header alone
 * (see `ms index-plan`).  Section payloads are the structures' own Orbit
 * serializations, each starting on a SECTION_ALIGN boundary.
 *
 * A writer serializes each section twice, first into a sink that counts its
 * bytes and checksums them, so that the header can go first and the output
 * need not seek.  A reader reads the header, then seeks (or, on a stream
 * that cannot seek, reads) past the sections it does not load, and checks
 * each section it loads against its entry: the bytes it took, and for a move
 * table its column codes, widths, rows and domain.
 *
 * Integers are little-endian; a header written on a big-endian host is
 * refused through its endian check.
 */

#ifndef _INDEX_FORMAT_HPP
#define _INDEX_FORMAT_HPP

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <functional>
#include <istream>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <vector>

namespace orbit_idx {

inline constexpr char MAGIC[8] = {'O', 'R', 'B', 'I', 'T', 'I', 'D', 'X'};
inline constexpr uint32_t FORMAT_VERSION = 1;
inline constexpr uint32_t ENDIAN_CHECK = 0x01020304;
inline constexpr uint32_t HEADER_BYTES = 4096;
inline constexpr uint64_t SECTION_ALIGN = 4096;
// Byte sizes of the header's parts (see INDEX_FORMAT.md).
inline constexpr size_t FIXED_BYTES = 128, SECTION_ENTRY_BYTES = 80;
inline constexpr size_t MAX_COLS = 12, SCALARS = 4;

enum class Tool : uint16_t { MS = 1, TMS = 2 };

/** A tms index's layout, or NONE for an ms index. */
enum class Layout : uint16_t { NONE = 0, FULL = 1, PSI = 2, PHI = 3 };

/** Section ids.  A set of sections is a mask with bit id set for each (see section_bit). */
enum class Section : uint16_t { LF = 1, F_STARTS = 2, PHI_LEN = 3, PHI_START = 4, PHI_INV = 5, MS_SPILL = 6, SEQ_IDS = 7 };

enum class Kind : uint16_t {
    MOVE_TABLE = 1,    // an Orbit permutation's serialization
    U64_ARRAY = 2,     // rows little-endian uint64_t values
    SPILL_ARRAYS = 3,  // rows byte arrays, each a uint64_t size and its bytes
    SEQ_TABLE = 4,     // rows sequences, each a uint64_t length, a uint32_t name size and the name
};

/** Column codes, one registry for every move table. */
enum class Col : uint8_t {
    LENGTH = 1, START = 2, POINTER = 3, OFFSET = 4, CHARACTER = 5,
    PSI_INT = 16, PSI_OFF = 17, PHI_INT = 18, PHI_OFF = 19, PLCP = 20, PLCPB = 21,
    TOP_LCP = 22, INTERIOR_MIN = 23,
    LCP_TOP = 32, LCP_MIN_SUB = 33, LCP_SPILL = 34,
};

inline const char* section_name(Section s) {
    switch (s) {
        case Section::LF: return "LF";
        case Section::F_STARTS: return "F_STARTS";
        case Section::PHI_LEN: return "PHI_LEN";
        case Section::PHI_START: return "PHI_START";
        case Section::PHI_INV: return "PHI_INV";
        case Section::MS_SPILL: return "MS_SPILL";
        case Section::SEQ_IDS: return "SEQ_IDS";
    }
    return "?";
}

inline const char* col_name(uint8_t c) {
    switch (static_cast<Col>(c)) {
        case Col::LENGTH: return "len";
        case Col::START: return "start";
        case Col::POINTER: return "ptr";
        case Col::OFFSET: return "off";
        case Col::CHARACTER: return "chr";
        case Col::PSI_INT: return "psi_int";
        case Col::PSI_OFF: return "psi_off";
        case Col::PHI_INT: return "phi_int";
        case Col::PHI_OFF: return "phi_off";
        case Col::PLCP: return "plcp";
        case Col::PLCPB: return "plcpb";
        case Col::TOP_LCP: return "top_lcp";
        case Col::INTERIOR_MIN: return "int_min";
        case Col::LCP_TOP: return "lcp_top";
        case Col::LCP_MIN_SUB: return "lcp_min_sub";
        case Col::LCP_SPILL: return "lcp_spill";
    }
    return "?";
}

inline const char* layout_name(Layout l) {
    switch (l) {
        case Layout::NONE: return "none";
        case Layout::FULL: return "full";
        case Layout::PSI: return "psi";
        case Layout::PHI: return "phi";
    }
    return "?";
}

constexpr uint32_t section_bit(Section s) { return uint32_t(1) << static_cast<unsigned>(s); }

struct SectionEntry {
    Section id{};
    Kind kind{};
    uint16_t row_bits = 0;  // a move table's bits per row
    uint64_t offset = 0, bytes = 0;
    uint64_t rows = 0;      // a move table's rows, or the values or arrays of the other kinds
    uint64_t domain = 0;    // a move table's domain
    uint32_t crc = 0;       // CRC-32 of the payload
    std::vector<std::pair<uint8_t, uint8_t>> cols;  // (code, width) in row order
};

struct Header {
    Tool tool{};
    Layout layout = Layout::NONE;
    uint64_t n = 0;         // text length, terminators included
    uint64_t bwt_runs = 0;  // maximal runs of one character in LF's rows
    // Tool scalars.  ms: max LCP top, max LCP min-sub, spill alignment, spill split bits.
    std::array<uint64_t, SCALARS> scalars{};
    std::vector<SectionEntry> sections;
    std::string info;  // build provenance, "key=value" lines

    const SectionEntry* find(Section s) const {
        for (const auto& e : sections)
            if (e.id == s) return &e;
        return nullptr;
    }
    uint32_t held_mask() const {
        uint32_t m = 0;
        for (const auto& e : sections) m |= section_bit(e.id);
        return m;
    }
};

// --- CRC-32 (IEEE, reflected), slicing by eight -------------------------

namespace detail {
struct CrcTables {
    uint32_t t[8][256];
    CrcTables() {
        for (uint32_t i = 0; i < 256; ++i) {
            uint32_t c = i;
            for (int k = 0; k < 8; ++k) c = (c & 1) ? 0xEDB88320u ^ (c >> 1) : c >> 1;
            t[0][i] = c;
        }
        for (uint32_t i = 0; i < 256; ++i)
            for (int s = 1; s < 8; ++s) t[s][i] = (t[s - 1][i] >> 8) ^ t[0][t[s - 1][i] & 0xFF];
    }
};
inline const CrcTables& crc_tables() {
    static const CrcTables tables;
    return tables;
}
}  // namespace detail

/** Extend crc, a CRC-32 of earlier bytes (0 for none), with n bytes at p. */
inline uint32_t crc32(uint32_t crc, const void* p, size_t n) {
    const auto& t = detail::crc_tables().t;
    const auto* b = static_cast<const unsigned char*>(p);
    crc = ~crc;
    while (n >= 8) {
        uint32_t lo, hi;
        std::memcpy(&lo, b, 4);
        std::memcpy(&hi, b + 4, 4);
        lo ^= crc;
        crc = t[7][lo & 0xFF] ^ t[6][(lo >> 8) & 0xFF] ^ t[5][(lo >> 16) & 0xFF] ^ t[4][lo >> 24] ^
              t[3][hi & 0xFF] ^ t[2][(hi >> 8) & 0xFF] ^ t[1][(hi >> 16) & 0xFF] ^ t[0][hi >> 24];
        b += 8;
        n -= 8;
    }
    while (n--) crc = t[0][(crc ^ *b++) & 0xFF] ^ (crc >> 8);
    return ~crc;
}

// --- Header bytes --------------------------------------------------------

namespace detail {
template <typename T>
void put(std::vector<char>& b, size_t at, T v) {
    std::memcpy(b.data() + at, &v, sizeof(v));
}
template <typename T>
T get(const char* b, size_t at) {
    T v;
    std::memcpy(&v, b + at, sizeof(v));
    return v;
}
}  // namespace detail

/** Write section entry e as the SECTION_ENTRY_BYTES bytes at b[at]. */
inline void encode_entry(std::vector<char>& b, size_t at, const SectionEntry& e) {
    using detail::put;
    if (e.cols.size() > MAX_COLS) throw std::length_error("a move table has more columns than the header holds");
    put<uint16_t>(b, at, static_cast<uint16_t>(e.id));
    put<uint16_t>(b, at + 2, static_cast<uint16_t>(e.kind));
    put<uint8_t>(b, at + 4, static_cast<uint8_t>(e.cols.size()));
    put<uint16_t>(b, at + 6, e.row_bits);
    put<uint64_t>(b, at + 8, e.offset);
    put<uint64_t>(b, at + 16, e.bytes);
    put<uint64_t>(b, at + 24, e.rows);
    put<uint64_t>(b, at + 32, e.domain);
    put<uint32_t>(b, at + 40, e.crc);
    for (size_t c = 0; c < e.cols.size(); ++c) {
        put<uint8_t>(b, at + 48 + 2 * c, e.cols[c].first);
        put<uint8_t>(b, at + 49 + 2 * c, e.cols[c].second);
    }
}

/**
 * The section entry in the SECTION_ENTRY_BYTES bytes at b + at.  Throws,
 * naming the file kind what, if its columns are malformed.
 */
inline SectionEntry decode_entry(const char* b, size_t at, const char* what) {
    using detail::get;
    SectionEntry e;
    e.id = static_cast<Section>(get<uint16_t>(b, at));
    e.kind = static_cast<Kind>(get<uint16_t>(b, at + 2));
    const size_t nc = get<uint8_t>(b, at + 4);
    e.row_bits = get<uint16_t>(b, at + 6);
    e.offset = get<uint64_t>(b, at + 8);
    e.bytes = get<uint64_t>(b, at + 16);
    e.rows = get<uint64_t>(b, at + 24);
    e.domain = get<uint64_t>(b, at + 32);
    e.crc = get<uint32_t>(b, at + 40);
    if (nc > MAX_COLS) throw std::runtime_error(std::string("malformed ") + what + " header: too many columns");
    unsigned bits = 0;
    for (size_t c = 0; c < nc; ++c) {
        e.cols.emplace_back(get<uint8_t>(b, at + 48 + 2 * c), get<uint8_t>(b, at + 49 + 2 * c));
        bits += e.cols.back().second;
    }
    if (bits != e.row_bits)
        throw std::runtime_error(std::string("malformed ") + what + " header: row bits disagree with widths");
    return e;
}

/** Bytes the header's tables and info take; at most HEADER_BYTES. */
inline size_t header_used(const Header& h) {
    return FIXED_BYTES + h.sections.size() * SECTION_ENTRY_BYTES + h.info.size();
}

/** The HEADER_BYTES bytes of h.  Throws if its tables and info do not fit. */
inline std::vector<char> encode_header(const Header& h) {
    using detail::put;
    if (header_used(h) > HEADER_BYTES) throw std::length_error("index header does not fit in its block");
    std::vector<char> b(HEADER_BYTES, 0);
    std::memcpy(b.data(), MAGIC, 8);
    put<uint32_t>(b, 8, FORMAT_VERSION);
    put<uint32_t>(b, 12, ENDIAN_CHECK);
    put<uint32_t>(b, 16, HEADER_BYTES);
    put<uint16_t>(b, 20, static_cast<uint16_t>(h.tool));
    put<uint16_t>(b, 22, static_cast<uint16_t>(h.layout));
    put<uint64_t>(b, 24, h.n);
    put<uint64_t>(b, 32, h.bwt_runs);
    put<uint16_t>(b, 40, static_cast<uint16_t>(h.sections.size()));
    put<uint16_t>(b, 44, static_cast<uint16_t>(h.info.size()));
    for (size_t i = 0; i < SCALARS; ++i) put<uint64_t>(b, 48 + 8 * i, h.scalars[i]);
    size_t at = FIXED_BYTES;
    for (const auto& e : h.sections) {
        encode_entry(b, at, e);
        at += SECTION_ENTRY_BYTES;
    }
    std::memcpy(b.data() + at, h.info.data(), h.info.size());
    return b;
}

/** True if the first bytes of a file, at least 8 of them, are an ORBITIDX magic. */
inline bool has_magic(const char* b, size_t n) { return n >= 8 && std::memcmp(b, MAGIC, 8) == 0; }

/** The header in the HEADER_BYTES bytes at b.  Throws if they are not a valid one. */
inline Header decode_header(const char* b) {
    using detail::get;
    if (!has_magic(b, 8)) throw std::runtime_error("not an Orbit index (no ORBITIDX magic)");
    if (get<uint32_t>(b, 12) != ENDIAN_CHECK) throw std::runtime_error("Orbit index written with another byte order");
    const uint32_t v = get<uint32_t>(b, 8);
    if (v != FORMAT_VERSION) throw std::runtime_error("unsupported Orbit index format version " + std::to_string(v));
    if (get<uint32_t>(b, 16) != HEADER_BYTES) throw std::runtime_error("unsupported Orbit index header size");
    Header h;
    h.tool = static_cast<Tool>(get<uint16_t>(b, 20));
    if (h.tool != Tool::MS && h.tool != Tool::TMS) throw std::runtime_error("Orbit index of an unknown tool");
    h.layout = static_cast<Layout>(get<uint16_t>(b, 22));
    h.n = get<uint64_t>(b, 24);
    h.bwt_runs = get<uint64_t>(b, 32);
    const size_t ns = get<uint16_t>(b, 40), ni = get<uint16_t>(b, 44);
    for (size_t i = 0; i < SCALARS; ++i) h.scalars[i] = get<uint64_t>(b, 48 + 8 * i);
    // Reserved bytes are 0, so that a later version can give them meaning.
    bool reserved_zero = get<uint16_t>(b, 42) == 0 && get<uint16_t>(b, 46) == 0;
    for (size_t i = 48 + 8 * SCALARS; i < FIXED_BYTES; ++i) reserved_zero = reserved_zero && b[i] == 0;
    if (!reserved_zero) throw std::runtime_error("malformed Orbit index header: reserved bytes are not 0");
    if (FIXED_BYTES + ns * SECTION_ENTRY_BYTES + ni > HEADER_BYTES)
        throw std::runtime_error("malformed Orbit index header: its tables overrun it");
    size_t at = FIXED_BYTES;
    uint64_t end = HEADER_BYTES;
    for (size_t i = 0; i < ns; ++i, at += SECTION_ENTRY_BYTES) {
        SectionEntry e = decode_entry(b, at, "Orbit index");
        if (e.offset % SECTION_ALIGN != 0 || e.offset < end)
            throw std::runtime_error("malformed Orbit index header: sections out of order or unaligned");
        if (h.find(e.id)) throw std::runtime_error("malformed Orbit index header: a section appears twice");
        end = e.offset + e.bytes;
        h.sections.push_back(std::move(e));
    }
    h.info.assign(b + at, ni);
    return h;
}

// --- Writing -------------------------------------------------------------

namespace detail {
/** A stream buffer that counts and checksums what is written to it, and forwards it to next if any. */
struct CrcCountBuf : std::streambuf {
    std::streambuf* next = nullptr;
    uint64_t count = 0;
    uint32_t crc = 0;
    explicit CrcCountBuf(std::streambuf* n = nullptr) : next(n) {}
    int_type overflow(int_type c) override {
        if (traits_type::eq_int_type(c, traits_type::eof())) return traits_type::not_eof(c);
        const char ch = traits_type::to_char_type(c);
        return xsputn(&ch, 1) == 1 ? c : traits_type::eof();
    }
    std::streamsize xsputn(const char* s, std::streamsize n) override {
        crc = crc32(crc, s, static_cast<size_t>(n));
        count += static_cast<uint64_t>(n);
        if (next && next->sputn(s, n) != n) return 0;
        return n;
    }
};

inline void write_zeros(std::ostream& out, uint64_t n) {
    static const char zeros[4096] = {};
    while (n > 0) {
        const uint64_t k = std::min<uint64_t>(n, sizeof(zeros));
        out.write(zeros, static_cast<std::streamsize>(k));
        n -= k;
    }
}

/**
 * The offset of the column widths in a serialized Orbit permutation whose
 * first have bytes are b: past the OBPT preamble, the split parameters (two
 * flags, each followed by 8 bytes when set), and the OBMS and OBPM
 * preambles.  0 if more bytes are needed to tell.
 */
inline size_t perm_widths_at(const char* b, size_t have) {
    size_t p = 4 + 24;
    for (int k = 0; k < 2; ++k) {
        if (have <= p) return 0;
        p += 1 + (b[p] ? 8 : 0);
    }
    return p + 4 + 24 + 16 + 4 + 24 + 8;
}

/**
 * Passes a serialized Orbit permutation to dst with the drop column widths
 * after the first keep left out.  The packed rows are unchanged when the
 * columns left out have width 0 and come last.
 */
struct DropWidthsBuf : std::streambuf {
    std::streambuf* dst;
    size_t keep, drop;
    std::string head;
    bool passing = false;
    DropWidthsBuf(std::streambuf* d, size_t k, size_t n) : dst(d), keep(k), drop(n) {}
    int_type overflow(int_type c) override {
        if (traits_type::eq_int_type(c, traits_type::eof())) return traits_type::not_eof(c);
        const char ch = traits_type::to_char_type(c);
        return xsputn(&ch, 1) == 1 ? c : traits_type::eof();
    }
    std::streamsize xsputn(const char* s, std::streamsize n) override {
        if (passing) return dst->sputn(s, n);
        head.append(s, static_cast<size_t>(n));
        const size_t p = perm_widths_at(head.data(), head.size());
        if (p && head.size() >= p + keep + drop) {
            passing = true;
            const auto first = static_cast<std::streamsize>(p + keep);
            const auto rest = static_cast<std::streamsize>(head.size() - (p + keep + drop));
            if (dst->sputn(head.data(), first) != first ||
                dst->sputn(head.data() + p + keep + drop, rest) != rest)
                return 0;
            head.clear();
        }
        return n;
    }
};

/** Serves the bytes of pre, then those of src. */
struct PrefixBuf : std::streambuf {
    std::string pre;
    std::streambuf* src;
    char ch = 0;
    PrefixBuf(std::string p, std::streambuf* s) : pre(std::move(p)), src(s) {
        setg(pre.data(), pre.data(), pre.data() + pre.size());
    }
    int_type underflow() override {
        const int_type c = src->sbumpc();
        if (traits_type::eq_int_type(c, traits_type::eof())) return c;
        ch = traits_type::to_char_type(c);
        setg(&ch, &ch, &ch + 1);
        return c;
    }
    std::streamsize xsgetn(char* s, std::streamsize n) override {
        const std::streamsize take = std::min<std::streamsize>(n, egptr() - gptr());
        std::memcpy(s, gptr(), static_cast<size_t>(take));
        gbump(static_cast<int>(take));
        return take + (n > take ? src->sgetn(s + take, n - take) : 0);
    }
};
}  // namespace detail

inline uint64_t align_up(uint64_t x) { return (x + SECTION_ALIGN - 1) / SECTION_ALIGN * SECTION_ALIGN; }

/** A section to write: its entry, apart from offset, bytes and crc, and how to write its payload. */
struct SectionSource {
    SectionEntry entry;
    std::function<void(std::ostream&)> write;
};

/**
 * Write a file with header h, whose sections are those of srcs, filled in.
 * Each payload is written twice: once to count and checksum it, and once to
 * out.  finish, if given, completes h between the two passes.  Returns the
 * bytes written.  Throws if a payload's second writing differs from its first
 * or out fails.
 */
inline uint64_t write_file(std::ostream& out, Header h, std::vector<SectionSource>& srcs,
                           const std::function<void(Header&)>& finish = {}) {
    h.sections.clear();
    uint64_t at = HEADER_BYTES;
    for (auto& s : srcs) {
        detail::CrcCountBuf sink;
        std::ostream os(&sink);
        s.write(os);
        s.entry.offset = align_up(at);
        s.entry.bytes = sink.count;
        s.entry.crc = sink.crc;
        at = s.entry.offset + s.entry.bytes;
        h.sections.push_back(s.entry);
    }
    if (finish) finish(h);
    const std::vector<char> hb = encode_header(h);
    out.write(hb.data(), static_cast<std::streamsize>(hb.size()));
    uint64_t pos = HEADER_BYTES;
    for (const auto& s : srcs) {
        detail::write_zeros(out, s.entry.offset - pos);
        detail::CrcCountBuf tee(out.rdbuf());
        std::ostream os(&tee);
        s.write(os);
        os.flush();
        if (tee.count != s.entry.bytes || tee.crc != s.entry.crc)
            throw std::logic_error(std::string("section ") + section_name(s.entry.id) + " serialized differently twice");
        pos = s.entry.offset + s.entry.bytes;
    }
    if (!out.good()) throw std::runtime_error("failed writing Orbit index");
    return pos;
}

/** A move-table section entry for a structure with the given column codes. */
template <typename Structure, size_t N>
SectionEntry move_table_entry(Section id, const Structure& x, const std::array<Col, N>& codes) {
    static_assert(N <= MAX_COLS);
    SectionEntry e;
    e.id = id;
    e.kind = Kind::MOVE_TABLE;
    e.rows = x.intervals();
    e.domain = x.domain();
    const auto& w = x.get_widths();
    static_assert(std::tuple_size_v<std::remove_cv_t<std::remove_reference_t<decltype(w)>>> == N,
                  "column codes and the structure's columns disagree in number");
    for (size_t c = 0; c < N; ++c) {
        e.cols.emplace_back(static_cast<uint8_t>(codes[c]), w[c]);
        e.row_bits = static_cast<uint16_t>(e.row_bits + w[c]);
    }
    return e;
}

/** A move-table section source that serializes x, which must outlive the write. */
template <typename Structure, size_t N>
SectionSource move_table_source(Section id, Structure& x, const std::array<Col, N>& codes) {
    return SectionSource{move_table_entry(id, x, codes), [&x](std::ostream& o) { x.serialize(o); }};
}

/**
 * As move_table_source, but the last optional columns are left out of the
 * section, entry and payload alike, where they and every column after them
 * have width 0, so that the section reads the same as one written by a
 * build without them.  Reader::load_move_table with the same optional reads
 * it back.
 */
template <typename Structure, size_t N>
SectionSource move_table_source(Section id, Structure& x, const std::array<Col, N>& codes, size_t optional) {
    SectionSource s = move_table_source(id, x, codes);
    const auto& w = x.get_widths();
    size_t keep = N;
    while (keep > N - optional && w[keep - 1] == 0) --keep;
    if (keep == N) return s;
    s.entry.cols.resize(keep);
    s.write = [&x, keep](std::ostream& o) {
        detail::DropWidthsBuf d(o.rdbuf(), keep, N - keep);
        std::ostream os(&d);
        x.serialize(os);
        if (!os.good() || !d.passing) o.setstate(std::ios::badbit);
    };
    return s;
}

// --- Reading -------------------------------------------------------------

/**
 * Reads a file's sections in order from a stream, counting the bytes taken
 * so that each section's size can be checked and the sections not loaded
 * skipped: by seeking where the stream can, else by reading past them.
 */
class Reader {
public:
    explicit Reader(std::istream& in) : buf_(in.rdbuf()), in_(&buf_) {
        std::vector<char> b(HEADER_BYTES);
        in_.read(b.data(), HEADER_BYTES);
        if (in_.gcount() < 8 || !has_magic(b.data(), 8)) throw std::runtime_error("not an Orbit index (no ORBITIDX magic)");
        if (!in_.good()) throw std::runtime_error("truncated Orbit index header");
        header_ = decode_header(b.data());
    }
    const Header& header() const { return header_; }

    /**
     * Load the move table of section id into x, which has the given column
     * codes, and check it against its entry.  The section may lack up to
     * the last optional columns, which are then loaded with width 0.
     * Throws if the section is missing or disagrees.
     */
    template <typename Structure, size_t N>
    void load_move_table(Section id, Structure& x, const std::array<Col, N>& codes, size_t optional = 0) {
        const SectionEntry& e = entry(id, Kind::MOVE_TABLE);
        const size_t held = e.cols.size();
        if (held > N || held + optional < N) fail(e, "has another number of columns than this build reads");
        for (size_t c = 0; c < held; ++c)
            if (e.cols[c].first != static_cast<uint8_t>(codes[c]))
                fail(e, std::string("has column ") + col_name(e.cols[c].first) + " where this build reads " +
                            col_name(static_cast<uint8_t>(codes[c])));
        seek(e.offset);
        try {
            if (held == N) {
                x.load(in_);
            } else {
                // The preamble and the widths held, then a zero width for
                // each column the section lacks, then the rest as stored.
                std::string head;
                for (size_t p = 0; !p || head.size() < p + held; p = detail::perm_widths_at(head.data(), head.size())) {
                    const int c = in_.get();
                    if (c == std::char_traits<char>::eof()) fail(e, "is truncated");
                    head.push_back(static_cast<char>(c));
                }
                head.append(N - held, '\0');
                detail::PrefixBuf pb(std::move(head), &buf_);
                std::istream is(&pb);
                x.load(is);
                if (!is.good()) in_.setstate(std::ios::failbit);
            }
        } catch (const std::exception&) {
            // A structure cut short fails a check of the bytes it did not get.
            if (!in_.good()) fail(e, "is truncated");
            throw;
        }
        check_end(e);
        const auto& w = x.get_widths();
        for (size_t c = 0; c < N; ++c)
            if (w[c] != (c < held ? e.cols[c].second : 0)) fail(e, "has widths other than its header entry");
        if (x.intervals() != e.rows || x.domain() != e.domain) fail(e, "has rows or domain other than its header entry");
    }

    /** The values of a U64_ARRAY section. */
    std::vector<uint64_t> load_u64s(Section id) {
        const SectionEntry& e = entry(id, Kind::U64_ARRAY);
        if (e.bytes != e.rows * 8) fail(e, "has a size other than 8 bytes per value");
        seek(e.offset);
        std::vector<uint64_t> v(e.rows);
        in_.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(e.bytes));
        check_end(e);
        return v;
    }

    /** The byte arrays of a SPILL_ARRAYS section. */
    std::vector<std::vector<unsigned char>> load_arrays(Section id) {
        const SectionEntry& e = entry(id, Kind::SPILL_ARRAYS);
        seek(e.offset);
        std::vector<std::vector<unsigned char>> v(e.rows);
        for (auto& a : v) {
            uint64_t size = 0;
            in_.read(reinterpret_cast<char*>(&size), sizeof(size));
            if (!in_.good() || size > e.bytes) fail(e, "is truncated or malformed");
            a.resize(size);
            in_.read(reinterpret_cast<char*>(a.data()), static_cast<std::streamsize>(size));
        }
        check_end(e);
        return v;
    }

    /** The payload of section id, which must be present and of kind k, as bytes. */
    std::string load_bytes(Section id, Kind k) {
        const SectionEntry& e = entry(id, k);
        seek(e.offset);
        std::string v(e.bytes, '\0');
        in_.read(v.data(), static_cast<std::streamsize>(e.bytes));
        check_end(e);
        if (crc32(0, v.data(), v.size()) != e.crc) fail(e, "fails its CRC-32");
        return v;
    }
    /** The entry of section id, which must be present and of kind k. */
    const SectionEntry& entry(Section id, Kind k) const {
        const SectionEntry* e = header_.find(id);
        if (!e) throw std::runtime_error(std::string("Orbit index lacks section ") + section_name(id));
        if (e->kind != k) fail(*e, "is of another kind than expected");
        return *e;
    }

private:
    /** Counts the bytes read through it from src, buffering one byte for underflow. */
    struct CountingBuf : std::streambuf {
        std::streambuf* src;
        uint64_t count = 0;
        char ch = 0;
        explicit CountingBuf(std::streambuf* s) : src(s) {}
        int_type underflow() override {
            const int_type c = src->sbumpc();
            if (traits_type::eq_int_type(c, traits_type::eof())) return c;
            ++count;
            ch = traits_type::to_char_type(c);
            setg(&ch, &ch, &ch + 1);
            return c;
        }
        std::streamsize xsgetn(char* s, std::streamsize n) override {
            std::streamsize got = 0;
            if (n > 0 && gptr() < egptr()) {
                *s++ = *gptr();
                gbump(1);
                got = 1;
                --n;
            }
            const std::streamsize k = n > 0 ? src->sgetn(s, n) : 0;
            count += static_cast<uint64_t>(k);
            return got + k;
        }
        /** Position: bytes taken by the reader, not counting one buffered and unread. */
        uint64_t pos() const { return count - static_cast<uint64_t>(egptr() - gptr()); }
        /** Move forward to offset to, seeking src if it can, else reading. */
        bool skip_to(uint64_t to) {
            const uint64_t at = pos();
            if (to < at) return false;
            if (gptr() < egptr()) {
                // The buffered byte is at position at; past it, src is at count.
                if (to == at) return true;
                gbump(1);
            }
            uint64_t gap = to - count;
            if (gap == 0) return true;
            if (src->pubseekoff(static_cast<std::streamoff>(gap), std::ios::cur, std::ios::in) != std::streampos(-1)) {
                count = to;
                return true;
            }
            char tmp[1 << 16];
            while (gap > 0) {
                const std::streamsize k = src->sgetn(tmp, static_cast<std::streamsize>(std::min<uint64_t>(gap, sizeof(tmp))));
                if (k <= 0) return false;
                gap -= static_cast<uint64_t>(k);
                count += static_cast<uint64_t>(k);
            }
            return true;
        }
    };

    [[noreturn]] static void fail(const SectionEntry& e, const std::string& what) {
        throw std::runtime_error(std::string("Orbit index section ") + section_name(e.id) + " " + what);
    }
    void seek(uint64_t offset) {
        if (!in_.good() || !buf_.skip_to(offset)) throw std::runtime_error("Orbit index sections out of reach (truncated?)");
    }
    void check_end(const SectionEntry& e) {
        if (!in_.good()) fail(e, "is truncated");
        if (buf_.pos() != e.offset + e.bytes) fail(e, "took a number of bytes other than its header entry");
    }

    CountingBuf buf_;
    std::istream in_;
    Header header_;
};

/** The header of the file at path, read from its first HEADER_BYTES alone. */
Header read_header_file(const std::string& path);

/** A description of h: its sections, with bytes and columns, and its info. */
void describe_header(const Header& h, std::ostream& os);

/**
 * Check every section's CRC-32 against its entry and that the file ends with
 * the last section.  Returns an empty string if all agree, else what differs.
 */
std::string verify_file(const std::string& path);

/**
 * Rewrite a legacy index (MSIX v4 or TMSX v8) at from as an ORBITIDX file at
 * to, copying each structure's bytes without loading it.  info is appended to
 * the provenance.  Returns an empty string on success, else the error.
 */
std::string rewrap_legacy(const std::string& from, const std::string& to, const std::string& info);
/**
 * Copy the index at from to to with section s added last, or put in place
 * of a section with the same id, and info appended to the provenance.  The
 * other sections' bytes are copied once each, checked against their CRC-32s,
 * and those before s keep their offsets.  Throws on
 * failure.
 */
void add_section(const std::string& from, const std::string& to, SectionSource& s, const std::string& info);

}  // namespace orbit_idx

#endif
