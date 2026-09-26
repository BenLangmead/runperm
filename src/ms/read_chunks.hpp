/**
 * Reads for batch queries, split into chunks of whole records so that
 * several threads can parse them.  ReadChunker::next, called by one thread
 * at a time, only finds record boundaries (a scan for line ends) and copies
 * the chunk's bytes; parse_chunk, which any thread may call, turns a chunk
 * into names and upper-cased sequences.
 *
 * The input is FASTA, FASTQ, or one sequence per line, chosen by the first
 * non-empty line ('>' FASTA, '@' FASTQ, else one per line).  Lines may end
 * in LF or CRLF.
 * - One per line: every non-empty line is a read, named read<i> with i
 *   counting reads from 1.
 * - FASTQ: blank lines before a record are skipped; a record is a header
 *   line, a sequence line, and a '+' line and a quality line, which may be
 *   missing at the end of the input.  A header with no sequence line after
 *   it is dropped.
 * - FASTA: a record is a header line and the lines up to the next line
 *   starting with '>', concatenated.
 * Names are the header after its first character, up to the first space,
 * tab or carriage return.
 */

#ifndef _READ_CHUNKS_HPP
#define _READ_CHUNKS_HPP

#include <algorithm>
#include <cstddef>
#include <cstring>
#include <istream>
#include <memory>
#include <string>
#include <vector>

enum class ReadFormat { UNKNOWN, FASTA, FASTQ, PLAIN };

/** The bytes of whole records, and the number of records before them. */
struct ReadChunk {
    std::string bytes;
    size_t first = 0;
    ReadFormat format = ReadFormat::UNKNOWN;
};

namespace read_chunks_detail {

// Length of the line [b, e) without a trailing carriage return.
inline size_t content_len(const char* b, const char* e) {
    return (e > b && e[-1] == '\r') ? size_t(e - b - 1) : size_t(e - b);
}

inline void header_name(const char* b, size_t len, std::string& name) {
    size_t end = 1;
    while (end < len && b[end] != ' ' && b[end] != '\t' && b[end] != '\r') ++end;
    name.assign(len > 0 ? b + 1 : b, len > 0 ? end - 1 : 0);
}

inline void append_upper(std::string& s, const char* b, size_t len) {
    const size_t at = s.size();
    s.append(b, len);
    for (size_t i = at; i < s.size(); ++i) {
        const char c = s[i];
        if (c >= 'a' && c <= 'z') s[i] = static_cast<char>(c - ('a' - 'A'));
    }
}

}  // namespace read_chunks_detail

/**
 * Parses the records of c into names[0 .. n) and seqs[0 .. n), reusing the
 * vectors' strings, and returns n.
 */
inline size_t parse_chunk(const ReadChunk& c, std::vector<std::string>& names, std::vector<std::string>& seqs) {
    using namespace read_chunks_detail;
    const char* p = c.bytes.data();
    const char* const end = p + c.bytes.size();
    size_t n = 0;
    // Next line [lb, le), with p moved past its line ending.
    const char *lb = nullptr, *le = nullptr;
    auto next_line = [&]() {
        if (p >= end) return false;
        lb = p;
        const void* nl = std::memchr(p, '\n', size_t(end - p));
        le = nl ? static_cast<const char*>(nl) : end;
        p = nl ? le + 1 : end;
        return true;
    };
    auto slot = [&]() {
        if (names.size() <= n) {
            names.resize(n + 1);
            seqs.resize(n + 1);
        }
        names[n].clear();
        seqs[n].clear();
    };
    if (c.format == ReadFormat::PLAIN) {
        while (next_line()) {
            const size_t len = content_len(lb, le);
            if (len == 0) continue;
            slot();
            names[n] = "read" + std::to_string(c.first + n + 1);
            append_upper(seqs[n], lb, len);
            ++n;
        }
    } else if (c.format == ReadFormat::FASTQ) {
        while (next_line()) {
            if (content_len(lb, le) == 0) continue;
            const char* hb = lb;
            const size_t hlen = content_len(lb, le);
            if (!next_line()) break;
            slot();
            header_name(hb, hlen, names[n]);
            append_upper(seqs[n], lb, content_len(lb, le));
            ++n;
            next_line();  // +
            next_line();  // qualities
        }
    } else if (c.format == ReadFormat::FASTA) {
        bool have = false;
        while (next_line()) {
            const size_t len = content_len(lb, le);
            if (!have && len == 0) continue;
            if (!have || (len > 0 && lb[0] == '>')) {
                slot();
                header_name(lb, len, names[n]);
                ++n;
                have = true;
            } else {
                append_upper(seqs[n - 1], lb, len);
            }
        }
    }
    return n;
}

/** Splits an input stream into chunks of whole records. */
class ReadChunker {
public:
    explicit ReadChunker(std::istream& in, size_t read_size = size_t(1) << 22) : in_(in), read_size_(read_size) {}

    /**
     * Moves the next records into c: at most max_records of them, and no
     * more once the chunk holds max_bytes bytes.  Returns false when the
     * input holds no more records.
     */
    bool next(ReadChunk& c, size_t max_records, size_t max_bytes) {
        using namespace read_chunks_detail;
        // Drop the bytes of earlier chunks once they are most of the
        // buffer, so each byte is moved at most about once.
        if (pos_ > size_ / 2) {
            std::memmove(buf_.get(), buf_.get() + pos_, size_ - pos_);
            size_ -= pos_;
            pos_ = 0;
        }
        if (format_ == ReadFormat::UNKNOWN && !detect()) return false;
        const size_t start = pos_;
        size_t at = pos_, records = 0;
        size_t lb, le;
        auto blank = [&](size_t b, size_t e) { return content_len(buf_.get() + b, buf_.get() + e) == 0; };
        while (records < max_records && at - start < max_bytes) {
            if (format_ == ReadFormat::PLAIN) {
                if (!line(at, lb, le)) break;
                at = after(le);
                if (!blank(lb, le)) ++records;
            } else if (format_ == ReadFormat::FASTQ) {
                size_t at2 = at;
                bool found = false;
                while (line(at2, lb, le)) {
                    at2 = after(le);
                    if (!blank(lb, le)) { found = true; break; }
                }
                if (!found) { at = at2; break; }
                if (!line(at2, lb, le)) { at = at2; break; }  // a header without a sequence is dropped
                at2 = after(le);
                for (int i = 0; i < 2 && line(at2, lb, le); ++i) at2 = after(le);
                at = at2;
                ++records;
            } else {
                size_t at2 = at;
                bool found = false;
                while (line(at2, lb, le)) {
                    at2 = after(le);
                    if (!blank(lb, le)) { found = true; break; }
                }
                if (!found) { at = at2; break; }
                // Lines up to the next header.
                while (line(at2, lb, le) && !(le > lb && buf_[lb] == '>')) at2 = after(le);
                at = at2;
                ++records;
            }
        }
        if (records == 0) {
            pos_ = at;
            return false;
        }
        c.bytes.assign(buf_.get() + start, at - start);
        c.first = count_;
        c.format = format_;
        count_ += records;
        pos_ = at;
        return true;
    }

private:
    // Finds the line starting at b: [lb, le) with le at its '\n' or at the
    // end of the input.  Reads more input as needed.  False at the end.
    bool line(size_t b, size_t& lb, size_t& le) {
        for (;;) {
            if (b < size_) {
                const void* nl = std::memchr(buf_.get() + b, '\n', size_ - b);
                if (nl) {
                    lb = b;
                    le = size_t(static_cast<const char*>(nl) - buf_.get());
                    return true;
                }
            }
            if (eof_) {
                if (b >= size_) return false;
                lb = b;
                le = size_;
                return true;
            }
            fill();
        }
    }
    size_t after(size_t le) const { return le < size_ ? le + 1 : le; }

    // Appends up to read_size bytes of input to the buffer.
    void fill() {
        if (size_ + read_size_ > cap_) {
            const size_t cap = std::max(2 * cap_, size_ + read_size_);
            std::unique_ptr<char[]> b(new char[cap]);
            if (size_ > 0) std::memcpy(b.get(), buf_.get(), size_);
            buf_ = std::move(b);
            cap_ = cap;
        }
        in_.read(buf_.get() + size_, static_cast<std::streamsize>(read_size_));
        const size_t got = static_cast<size_t>(in_.gcount());
        size_ += got;
        if (got < read_size_) eof_ = true;
    }

    // Sets the format from the first non-empty line.
    bool detect() {
        size_t at = pos_, lb, le;
        while (line(at, lb, le)) {
            if (read_chunks_detail::content_len(buf_.get() + lb, buf_.get() + le) > 0) {
                const char c = buf_[lb];
                format_ = c == '>' ? ReadFormat::FASTA : c == '@' ? ReadFormat::FASTQ : ReadFormat::PLAIN;
                return true;
            }
            at = after(le);
        }
        return false;
    }

    std::istream& in_;
    size_t read_size_;
    std::unique_ptr<char[]> buf_;  // input bytes [0, size_), of which [pos_, size_) are unread
    size_t size_ = 0, cap_ = 0;
    size_t pos_ = 0;
    size_t count_ = 0;
    bool eof_ = false;
    ReadFormat format_ = ReadFormat::UNKNOWN;
};

#endif /* _READ_CHUNKS_HPP */
