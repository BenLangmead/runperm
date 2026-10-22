/**
 * Contig names and lengths for an index built from a ropebwt3 FMD, and the
 * map from a text position to a contig, strand and forward-strand offset.
 *
 * The table is read from ropebwt3's sequence-ID file (the .len.gz next to an
 * .fmd, plain or gzipped): one line per contig, its name and length
 * separated by whitespace, with any further fields ignored.  As in
 * ropebwt3's rb3_sid_read, a line without a positive length is skipped.
 *
 * Contig k, of length len_k, is strings 2k (its forward strand) and 2k + 1
 * (its reverse complement) of the FMD, so the indexed text is
 *
 *     C_0 $ rc(C_0) $ C_1 $ rc(C_1) $ ...
 *
 * and strand s of contig k starts at text position 2 * sum_{j<k} (len_j + 1)
 * + s * (len_k + 1).  The text length is 2 * sum_k (len_k + 1), which is
 * checked against the index.
 *
 * In a tms index the table is the SEQ_IDS section (see INDEX_FORMAT.md): per
 * contig, in order, its length as a uint64_t, its name's byte count as a
 * uint32_t and the name, all little-endian.
 *
 * Locating a position, done once per listed SMEM position, reads a bucket
 * table over the high bits of text positions (build_lookup), which narrows
 * the search for its strand to the few strands that start in its bucket,
 * and appends the strand's "name:+:" or "name:-:" label in one copy.
 */

#ifndef _MS_SEQ_IDS_HPP
#define _MS_SEQ_IDS_HPP

#include "format_int.hpp"
#include <algorithm>
#include <charconv>
#include <cstdint>
#include <cstring>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <zlib.h>

struct SeqIdTable {
    std::vector<std::string> names;
    std::vector<uint64_t> lens;
    // starts[2k + s]: the text position of strand s of contig k.
    std::vector<uint64_t> starts;
    // label_text[label_at[j] .. label_at[j + 1]): "name:+:" or "name:-:" for
    // strand j = 2k + s, strand s of contig k, all in one buffer.
    std::string label_text;
    std::vector<uint32_t> label_at{0};
    // The longest label.
    size_t max_label = 0;
    // bucket[b]: the strand holding text position b << bucket_shift; empty
    // until build_lookup, and then the last entry is the last strand.
    std::vector<uint32_t> bucket;
    unsigned bucket_shift = 0;

    /** The text length the table describes: both strands of every contig, each with its terminator. */
    uint64_t text_length() const { return starts.empty() ? 0 : starts.back() + lens.back() + 1; }

    /** Append a contig of length len named name. */
    void add(const std::string& name, uint64_t len) {
        const uint64_t at = text_length();
        names.push_back(name);
        lens.push_back(len);
        starts.push_back(at);
        starts.push_back(at + len + 1);
        for (const char* tail : {":+:", ":-:"}) {
            label_text += name;
            label_text += tail;
            max_label = std::max(max_label, name.size() + 3);
            if (label_text.size() > UINT32_MAX) throw std::runtime_error("sequence names too long for the table");
            label_at.push_back(static_cast<uint32_t>(label_text.size()));
        }
        bucket.clear();
    }

    /**
     * Build the bucket table that narrows locate's search, with about four
     * buckets per strand.  read and from_payload build it; a table extended
     * with add afterwards searches all strands until it is built again.
     */
    void build_lookup() {
        bucket.clear();
        const uint64_t n = text_length();
        if (n == 0 || starts.size() > UINT32_MAX) return;
        bucket_shift = 0;
        while ((n - 1) >> bucket_shift >= 4 * starts.size()) ++bucket_shift;
        const uint64_t buckets = ((n - 1) >> bucket_shift) + 1;
        bucket.resize(buckets + 1);
        size_t j = 0;
        for (uint64_t b = 0; b < buckets; ++b) {
            const uint64_t x = b << bucket_shift;
            while (j + 1 < starts.size() && starts[j + 1] <= x) ++j;
            bucket[b] = static_cast<uint32_t>(j);
        }
        bucket[buckets] = static_cast<uint32_t>(starts.size() - 1);
    }

    /** The strand holding text position p, or starts.size() if p precedes every strand. */
    size_t strand_of(uint64_t p) const {
        if (bucket.empty() || p >= text_length())
            return static_cast<size_t>(std::upper_bound(starts.begin(), starts.end(), p) - starts.begin()) - 1;
        const uint64_t b = p >> bucket_shift;
        // The strand lies between the ones holding the bucket's first
        // position and the next bucket's.
        size_t lo = bucket[b], hi = bucket[b + 1];
        while (lo < hi) {
            const size_t mid = (lo + hi + 1) / 2;
            if (starts[mid] <= p) lo = mid;
            else hi = mid - 1;
        }
        return lo;
    }

    /** Write the table as a SEQ_IDS section payload. */
    void serialize(std::ostream& out) const {
        for (size_t k = 0; k < names.size(); ++k) {
            const uint64_t len = lens[k];
            const uint32_t nb = static_cast<uint32_t>(names[k].size());
            out.write(reinterpret_cast<const char*>(&len), sizeof(len));
            out.write(reinterpret_cast<const char*>(&nb), sizeof(nb));
            out.write(names[k].data(), nb);
        }
    }

    /** The table of n contigs in the SEQ_IDS payload b; throws if it is malformed. */
    static SeqIdTable from_payload(const std::string& b, uint64_t n) {
        SeqIdTable t;
        size_t at = 0;
        for (uint64_t k = 0; k < n; ++k) {
            uint64_t len;
            uint32_t nb;
            if (b.size() - at < sizeof(len) + sizeof(nb)) throw std::runtime_error("malformed SEQ_IDS section");
            std::memcpy(&len, b.data() + at, sizeof(len));
            std::memcpy(&nb, b.data() + at + sizeof(len), sizeof(nb));
            at += sizeof(len) + sizeof(nb);
            if (b.size() - at < nb || len == 0) throw std::runtime_error("malformed SEQ_IDS section");
            t.add(b.substr(at, nb), len);
            at += nb;
        }
        if (at != b.size()) throw std::runtime_error("malformed SEQ_IDS section");
        t.build_lookup();
        return t;
    }

    /** Read path, plain or gzipped; throws std::runtime_error on failure or an empty table. */
    static SeqIdTable read(const std::string& path) {
        gzFile f = gzopen(path.c_str(), "rb");
        if (!f) throw std::runtime_error("cannot open " + path);
        SeqIdTable t;
        std::string line;
        char buf[1 << 16];
        auto take = [&] {
            std::istringstream in(line);
            std::string name, len;
            if (!(in >> name >> len)) return;
            uint64_t v = 0;
            try {
                if (len.empty() || len[0] == '-') return;
                size_t used = 0;
                v = std::stoull(len, &used);
                if (used != len.size()) return;
            } catch (const std::exception&) {
                return;
            }
            if (v == 0) return;
            t.add(name, v);
        };
        while (gzgets(f, buf, sizeof(buf))) {
            line += buf;
            if (!line.empty() && line.back() == '\n') {
                take();
                line.clear();
            }
        }
        const bool err = !gzeof(f);
        gzclose(f);
        if (err) throw std::runtime_error("error reading " + path);
        if (!line.empty()) take();
        if (t.names.empty()) throw std::runtime_error(path + " lists no sequences");
        t.build_lookup();
        return t;
    }

    /**
     * The contig k, strand (0 for +, 1 for -) and forward-strand offset of
     * the match of length len at text position p.  On the - strand the
     * offset is that of the match's reverse complement on the forward
     * strand, len_k - (local + len), where local is p's offset in the
     * strand.  Throws std::runtime_error if the match leaves its strand,
     * which a match of a read never does in an index of this table's text.
     */
    void locate(uint64_t p, uint64_t len, size_t& k, int& strand, uint64_t& off) const {
        if (starts.empty() || p < starts[0])
            throw std::runtime_error("text position " + std::to_string(p) + " precedes every sequence");
        const size_t j = strand_of(p);
        k = j / 2;
        strand = static_cast<int>(j % 2);
        const uint64_t local = p - starts[j];
        if (local + len > lens[k])
            throw std::runtime_error("match at text position " + std::to_string(p) + " of length " +
                                     std::to_string(len) + " crosses the end of " + names[k] + " (" +
                                     (strand ? "-" : "+") + " strand, length " + std::to_string(lens[k]) + ")");
        off = strand ? lens[k] - (local + len) : local;
    }

    /**
     * The strand j (2k + s) and forward-strand offset off of the match of
     * length len at text position p, as locate gives them.
     */
    void locate_strand(uint64_t p, uint64_t len, size_t& j, uint64_t& off) const {
        size_t k;
        int strand;
        locate(p, len, k, strand, off);
        j = 2 * k + static_cast<size_t>(strand);
    }

    /** Write strand j's label and offset off at out, as put_location does. */
    char* put_strand_offset(char* out, size_t j, uint64_t off) const {
        const size_t n = label_at[j + 1] - label_at[j];
        std::memcpy(out, label_text.data() + label_at[j], n);
        return std::to_chars(out + n, out + n + 20, off).ptr;
    }

    /** The most bytes put_location writes. */
    size_t max_location() const { return max_label + 20; }

    /**
     * Write "name:strand:offset" for the match of length len at text
     * position p at out, which has room for max_location() bytes, and return
     * the end of what it wrote.
     */
    char* put_location(char* out, uint64_t p, uint64_t len) const {
        size_t j;
        uint64_t off;
        locate_strand(p, len, j, off);
        return put_strand_offset(out, j, off);
    }

    /** Append "name:strand:offset" for the match of length len at text position p. */
    void append_location(std::string& line, uint64_t p, uint64_t len) const {
        const size_t at = line.size();
        line.resize(at + max_location());
        line.resize(static_cast<size_t>(put_location(&line[at], p, len) - line.data()));
    }
};

#endif /* _MS_SEQ_IDS_HPP */
