/**
 * Reads for batch queries, split into chunks of whole records so that
 * several threads can parse them.  ReadChunker::next only finds record
 * boundaries (a scan for line ends); ReadPrefetcher runs it on a thread of
 * its own; and parse_chunk, which any thread may call, turns a chunk into
 * names and upper-cased sequences.
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
#include <chrono>
#include <condition_variable>
#include <deque>
#include <exception>
#include <mutex>
#include <thread>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <istream>
#include <memory>
#include <string>
#include <vector>

enum class ReadFormat { UNKNOWN, FASTA, FASTQ, PLAIN };

/** A block of input bytes, shared by the chunks that point into it. */
struct ReadSegment {
    std::unique_ptr<char[]> data;
    size_t size = 0, capacity = 0;
};

/**
 * Segments no longer referenced by any chunk, for reuse, so that reading
 * does not allocate and fault in fresh memory for every segment.  Chunks
 * are released on the query threads, so the pool is locked.
 */
class ReadSegmentPool : public std::enable_shared_from_this<ReadSegmentPool> {
public:
    /** A segment of at least `capacity` bytes, returned here once unused. */
    std::shared_ptr<ReadSegment> take(size_t capacity) {
        std::unique_ptr<ReadSegment> seg;
        {
            std::lock_guard<std::mutex> lk(mu_);
            for (size_t i = 0; i < free_.size(); ++i)
                if (free_[i]->capacity >= capacity) {
                    seg = std::move(free_[i]);
                    free_[i] = std::move(free_.back());
                    free_.pop_back();
                    break;
                }
        }
        if (!seg) {
            seg.reset(new ReadSegment);
            seg->data.reset(new char[capacity]);
            seg->capacity = capacity;
        }
        seg->size = 0;
        std::weak_ptr<ReadSegmentPool> pool = shared_from_this();
        return std::shared_ptr<ReadSegment>(seg.release(), [pool](ReadSegment* s) {
            if (auto p = pool.lock()) p->give(s);
            else delete s;
        });
    }

private:
    // Keeps a bounded number of segments; the rest are freed.
    void give(ReadSegment* s) {
        std::unique_ptr<ReadSegment> seg(s);
        std::lock_guard<std::mutex> lk(mu_);
        if (free_.size() < 64) free_.push_back(std::move(seg));
    }
    std::mutex mu_;
    std::vector<std::unique_ptr<ReadSegment>> free_;
};

/** The bytes of whole records, and the number of records before them. */
struct ReadChunk {
    std::shared_ptr<const ReadSegment> seg;  // keeps data alive
    const char* data = nullptr;
    size_t size = 0;
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
    const char* p = c.data;
    const char* const end = p + c.size;
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

/**
 * Splits an input stream into chunks of whole records.  The input is read
 * into segments, and a chunk is a range of one segment, shared with the
 * chunk rather than copied.  A record left incomplete at the end of a
 * segment is copied to the start of the next one.
 */
class ReadChunker {
public:
    explicit ReadChunker(std::istream& in, size_t read_size = size_t(1) << 22) : in_(in), read_size_(read_size) {}

    /**
     * Sets c to the next records: at most max_records of them, and no more
     * once the chunk holds max_bytes bytes.  Returns false when the input
     * holds no more records.
     */
    bool next(ReadChunk& c, size_t max_records, size_t max_bytes) {
        using namespace read_chunks_detail;
        for (;;) {
            if (format_ == ReadFormat::UNKNOWN) {
                const int d = detect();
                if (d == END) return false;
                if (d == MORE) { refill(); continue; }
            }
            const char* const buf = seg_ ? seg_->data.get() : nullptr;
            auto blank = [&](size_t b, size_t e) { return content_len(buf + b, buf + e) == 0; };
            const size_t start = pos_;
            size_t at = pos_, records = 0, lb, le;
            int r = FOUND;
            // The byte cap ends a chunk only once it holds a record, since
            // blank lines use bytes without adding records.
            while (records < max_records && (records == 0 || at - start < max_bytes)) {
                size_t at2 = at;
                if (format_ == ReadFormat::PLAIN) {
                    if ((r = line(at2, lb, le)) != FOUND) break;
                    at = after(le);
                    if (!blank(lb, le)) ++records;
                    continue;
                }
                // FASTQ and FASTA: skip blank lines, then a header.
                while ((r = line(at2, lb, le)) == FOUND && blank(lb, le)) at2 = after(le);
                if (r != FOUND) {
                    if (r == END) at = at2;
                    break;
                }
                at2 = after(le);
                if (format_ == ReadFormat::FASTQ) {
                    // Sequence, then '+' and qualities if present.
                    if ((r = line(at2, lb, le)) != FOUND) {
                        if (r == END) at = at2;  // a header without a sequence is dropped
                        break;
                    }
                    at2 = after(le);
                    for (int i = 0; i < 2 && (r = line(at2, lb, le)) == FOUND; ++i) at2 = after(le);
                    if (r == MORE) break;
                } else {
                    // Lines up to the next header.
                    while ((r = line(at2, lb, le)) == FOUND && !(le > lb && buf[lb] == '>')) at2 = after(le);
                    if (r == MORE) break;
                }
                at = at2;
                ++records;
                r = FOUND;
            }
            if (records > 0) {
                c.seg = seg_;
                c.data = buf + start;
                c.size = at - start;
                c.first = count_;
                c.format = format_;
                count_ += records;
                pos_ = at;
                return true;
            }
            pos_ = at;
            if (r != MORE) return false;
            refill();
        }
    }

private:
    enum { FOUND, MORE, END };

    // Finds the line starting at b: [lb, le) with le at its '\n', or at the
    // end of the input.  MORE when the segment ends before a line end and
    // input remains, END at the end of the input.
    int line(size_t b, size_t& lb, size_t& le) const {
        const size_t size = seg_ ? seg_->size : 0;
        if (b < size) {
            const char* buf = seg_->data.get();
            const void* nl = std::memchr(buf + b, '\n', size - b);
            if (nl) {
                lb = b;
                le = size_t(static_cast<const char*>(nl) - buf);
                return FOUND;
            }
        }
        if (!eof_) return MORE;
        if (b >= size) return END;
        lb = b;
        le = size;
        return FOUND;
    }
    size_t after(size_t le) const { return le < seg_->size ? le + 1 : le; }

    // Starts a new segment with the unread bytes of the current one and
    // more input.  Chunks of the old segment keep it alive.
    void refill() {
        const size_t carry = seg_ ? seg_->size - pos_ : 0;
        // Reading at least as much as is carried keeps the copying of a
        // record longer than read_size linear in its length.
        const size_t want = std::max(read_size_, carry);
        auto seg = pool_->take(carry + want);
        if (carry > 0) std::memcpy(seg->data.get(), seg_->data.get() + pos_, carry);
        in_.read(seg->data.get() + carry, static_cast<std::streamsize>(want));
        const size_t got = static_cast<size_t>(in_.gcount());
        seg->size = carry + got;
        if (got < want) eof_ = true;
        seg_ = std::move(seg);
        pos_ = 0;
    }

    // Sets the format from the first non-empty line.
    int detect() {
        size_t at = pos_, lb, le;
        int r;
        while ((r = line(at, lb, le)) == FOUND) {
            const char* buf = seg_->data.get();
            if (read_chunks_detail::content_len(buf + lb, buf + le) > 0) {
                const char c = buf[lb];
                format_ = c == '>' ? ReadFormat::FASTA : c == '@' ? ReadFormat::FASTQ : ReadFormat::PLAIN;
                return FOUND;
            }
            at = after(le);
        }
        pos_ = at;
        return r;
    }

    std::istream& in_;
    size_t read_size_;
    std::shared_ptr<ReadSegmentPool> pool_ = std::make_shared<ReadSegmentPool>();
    std::shared_ptr<ReadSegment> seg_;  // bytes [pos_, seg_->size) are unread
    size_t pos_ = 0;
    size_t count_ = 0;
    bool eof_ = false;
    ReadFormat format_ = ReadFormat::UNKNOWN;
};

/**
 * Runs a ReadChunker on its own thread, which keeps up to `capacity`
 * chunks ready.  Keeping the input buffer on one thread avoids moving it
 * between the cores of the threads that take chunks.
 */
class ReadPrefetcher {
public:
    ReadPrefetcher(std::istream& in, size_t max_records, size_t max_bytes, size_t capacity)
        : chunker_(in), max_records_(max_records), max_bytes_(max_bytes), capacity_(std::max<size_t>(capacity, 1)),
          thread_([this] { run(); }) {}
    ~ReadPrefetcher() {
        {
            std::lock_guard<std::mutex> lk(mu_);
            stop_ = true;
        }
        cv_.notify_all();
        thread_.join();
    }

    /** Takes the next chunk, waiting for one; false at the end of the input. */
    bool next(ReadChunk& c) {
        std::unique_lock<std::mutex> lk(mu_);
        cv_.wait(lk, [&] { return !ready_.empty() || done_; });
        if (!ready_.empty()) {
            c = std::move(ready_.front());
            ready_.pop_front();
            cv_.notify_all();
            return true;
        }
        if (error_) std::rethrow_exception(error_);
        return false;
    }

    /** Seconds the reading thread spent reading and splitting input. */
    double busy_seconds() const {
        std::lock_guard<std::mutex> lk(mu_);
        return busy_s_;
    }

private:
    void run() {
        using clock = std::chrono::steady_clock;
        double busy = 0.0;
        try {
            for (;;) {
                {
                    std::unique_lock<std::mutex> lk(mu_);
                    cv_.wait(lk, [&] { return stop_ || ready_.size() < capacity_; });
                    if (stop_) break;
                }
                auto t0 = clock::now();
                ReadChunk c;
                const bool got = chunker_.next(c, max_records_, max_bytes_);
                busy += std::chrono::duration<double>(clock::now() - t0).count();
                if (!got) break;
                std::lock_guard<std::mutex> lk(mu_);
                ready_.push_back(std::move(c));
                cv_.notify_all();
            }
        } catch (...) {
            std::lock_guard<std::mutex> lk(mu_);
            error_ = std::current_exception();
        }
        std::lock_guard<std::mutex> lk(mu_);
        done_ = true;
        busy_s_ = busy;
        cv_.notify_all();
    }

    ReadChunker chunker_;
    size_t max_records_, max_bytes_, capacity_;
    mutable std::mutex mu_;
    std::condition_variable cv_;
    std::deque<ReadChunk> ready_;
    bool done_ = false, stop_ = false;
    double busy_s_ = 0.0;
    std::exception_ptr error_;
    std::thread thread_;  // last, so it starts after the members it uses
};

/** A buffer owned by one query thread, which ReadCutter reads blocks into. */
struct ReadBuffer {
    std::unique_ptr<char[]> data;
    size_t capacity = 0;
    // Grows to at least n bytes, keeping the first `keep` bytes.
    void reserve(size_t n, size_t keep) {
        if (n <= capacity) return;
        const size_t cap = std::max(n, 2 * capacity);
        std::unique_ptr<char[]> d(new char[cap]);
        if (keep > 0) std::memcpy(d.get(), data.get(), keep);
        data = std::move(d);
        capacity = cap;
    }
};

namespace read_chunks_detail {

// Bytes of an 8-byte word equal to c, as the high bit of each such byte.
inline uint64_t byte_eq_mask(uint64_t w, uint64_t c) {
    const uint64_t lo7 = 0x7F7F7F7F7F7F7F7FULL;
    const uint64_t x = w ^ (c * 0x0101010101010101ULL);
    return ~(((x & lo7) + lo7) | x | lo7);
}
// Number of bytes flagged in a byte_eq_mask result.
inline unsigned mask_count(uint64_t m) { return unsigned(((m >> 7) * 0x0101010101010101ULL) >> 56); }
inline uint64_t load64(const char* p) {
    uint64_t w;
    std::memcpy(&w, p, 8);
    return w;
}

/**
 * Number of non-empty lines in [b, b + n), which starts at a line start,
 * counting a last line without a line ending.  Counts line ends and empty
 * lines (a '\n' right after a '\n' or at the start) a word at a time; a
 * block with a carriage return is counted line by line.
 */
inline size_t count_nonempty_lines(const char* b, size_t n) {
    size_t lines = 0, empty = 0, cr = 0, i = 0;
    for (; i + 9 <= n; i += 8) {
        const uint64_t m0 = byte_eq_mask(load64(b + i), '\n');
        const uint64_t m1 = byte_eq_mask(load64(b + i + 1), '\n');
        lines += mask_count(m0);
        empty += mask_count(m0 & m1);  // '\n' at i + k + 1 right after one at i + k
        cr |= byte_eq_mask(load64(b + i), '\r');
    }
    for (; i < n; ++i) {
        lines += b[i] == '\n';
        empty += b[i] == '\n' && i + 1 < n && b[i + 1] == '\n';
        cr |= b[i] == '\r';
    }
    if (cr) {
        size_t count = 0;
        for (const char *p = b, *end = b + n; p < end;) {
            const void* nl = std::memchr(p, '\n', size_t(end - p));
            const char* le = nl ? static_cast<const char*>(nl) : end;
            count += content_len(p, le) > 0;
            p = nl ? le + 1 : end;
        }
        return count;
    }
    // Each '\n' ends a line; pairs counted above end empty lines, as does a
    // '\n' at the start.  A last line without a line ending is non-empty.
    if (n > 0 && b[0] == '\n') ++empty;
    return lines - empty + (n > 0 && b[n - 1] != '\n');
}

}  // namespace read_chunks_detail

/**
 * Splits an input stream into blocks under the caller's lock, doing as
 * little as possible there: it reads bytes into the calling thread's own
 * buffer, after the partial record left by the previous block, and cuts the
 * block after its last whole record.  For one sequence per line the cut is
 * the last line end, and for FASTA the last line starting with '>'; both are
 * found from the end.  FASTQ records are found line by line.  The bytes after
 * the cut are kept for the next block.  For one sequence per line, the lines
 * are counted a word at a time, to number the reads.
 */
class ReadCutter {
public:
    explicit ReadCutter(std::istream& in) : in_(in) {}

    /**
     * Reads a block of about target bytes (more if one record is longer)
     * into buf and sets c to its whole records.  Returns false at the end of
     * the input.  A block may hold only blank lines.
     */
    bool next(ReadBuffer& buf, ReadChunk& c, size_t target) {
        using namespace read_chunks_detail;
        target = std::max<size_t>(target, 1);
        size_t size = carry_.size();
        buf.reserve(size + target, 0);
        if (size > 0) std::memcpy(buf.data.get(), carry_.data(), size);
        carry_.clear();
        size_t cut = 0, want = target, searched = 0;
        for (;;) {
            if (!eof_) {
                // Reading at least as much again as the block holds keeps a
                // record longer than the target linear to read and search.
                want = std::max(want, size);
                buf.reserve(size + want, size);
                in_.read(buf.data.get() + size, static_cast<std::streamsize>(want));
                const size_t got = static_cast<size_t>(in_.gcount());
                size += got;
                if (got < want) eof_ = true;
            }
            if (format_ == ReadFormat::UNKNOWN) detect(buf.data.get(), size);
            if (format_ == ReadFormat::UNKNOWN && !eof_) continue;  // only blank lines so far
            if (eof_) {
                cut = size;
                break;
            }
            cut = find_cut(buf.data.get(), size, searched);
            if (cut > 0) break;
            searched = size;
        }
        carry_.assign(buf.data.get() + cut, size - cut);
        if (cut == 0) return false;
        c.seg.reset();
        c.data = buf.data.get();
        c.size = cut;
        c.first = count_;
        c.format = format_;
        if (format_ == ReadFormat::PLAIN) count_ += count_nonempty_lines(c.data, cut);
        return true;
    }

private:
    // Sets the format from the first non-empty line of [b, b + n), if it has
    // a whole one or the input has ended.
    void detect(const char* b, size_t n) {
        for (const char *p = b, *end = b + n; p < end;) {
            const void* nl = std::memchr(p, '\n', size_t(end - p));
            if (!nl && !eof_) return;
            const char* le = nl ? static_cast<const char*>(nl) : end;
            if (read_chunks_detail::content_len(p, le) > 0) {
                format_ = *p == '>' ? ReadFormat::FASTA : *p == '@' ? ReadFormat::FASTQ : ReadFormat::PLAIN;
                return;
            }
            p = nl ? le + 1 : end;
        }
    }

    // End of the last whole record in [b, b + n), which starts at a record
    // start (or blank lines before one); 0 if there is none.  [b, b + done)
    // is known to hold no cut, so the search from the end stops there.
    size_t find_cut(const char* b, size_t n, size_t done) const {
        using namespace read_chunks_detail;
        const size_t stop = done > 0 ? done - 1 : 0;
        if (format_ == ReadFormat::PLAIN) {
            for (size_t i = n; i > stop; --i)
                if (b[i - 1] == '\n') return i;
            return 0;
        }
        if (format_ == ReadFormat::FASTA) {
            // The last header after the first one starts the partial record.
            size_t first = 0;
            while (first < n && b[first] != '>') ++first;  // blank lines before the first header
            for (size_t i = n; i > std::max(first + 1, stop); --i)
                if (b[i - 1] == '>' && b[i - 2] == '\n') return i - 1;
            return 0;
        }
        // FASTQ: blank lines, then a header, a sequence, '+' and qualities.
        size_t at = 0, last = 0;
        auto line_end = [&](size_t from, size_t& le) {
            const void* nl = from < n ? std::memchr(b + from, '\n', n - from) : nullptr;
            if (!nl) return false;
            le = size_t(static_cast<const char*>(nl) - b);
            return true;
        };
        for (;;) {
            size_t le;
            if (!line_end(at, le)) return last;
            if (content_len(b + at, b + le) == 0) {
                at = le + 1;
                last = at;  // blank lines before a header can end a block
                continue;
            }
            size_t p = le + 1;
            int lines = 1;
            for (; lines < 4 && line_end(p, le); ++lines) p = le + 1;
            if (lines < 4) return last;
            at = last = p;
        }
    }

    std::istream& in_;
    std::string carry_;  // bytes after the last cut: the start of a partial record
    size_t count_ = 0;
    bool eof_ = false;
    ReadFormat format_ = ReadFormat::UNKNOWN;
};

#endif /* _READ_CHUNKS_HPP */
