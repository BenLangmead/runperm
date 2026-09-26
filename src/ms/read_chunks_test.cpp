/**
 * Differential test of ReadCutter and parse_chunk.  Each input is split into
 * blocks of every size from 1 to 300 bytes, so that block boundaries fall at
 * every byte offset of its records, and every split must give the records
 * that a plain line-by-line reader (Reference, below) gives, including the
 * numbering of reads given one per line.  The inputs cover FASTA, FASTQ and
 * one sequence per line, with CRLF line endings, blank lines, lower case,
 * quality strings of '@', a record longer than the blocks, and odd records
 * at the end.  Build with the address and undefined behavior sanitizers
 * (make test-reader).
 */

#include "read_chunks.hpp"
#include <cctype>
#include <iostream>
#include <random>
#include <sstream>

using Records = std::vector<std::pair<std::string, std::string>>;

/** A reader of whole records one line at a time, for comparison. */
class Reference {
public:
    explicit Reference(std::istream& in) : in_(in) {}

    bool next(std::string& name, std::string& seq) {
        name.clear();
        seq.clear();
        std::string line;
        if (!have_line_) {
            do {
                if (!get_line(line)) return false;
            } while (line.empty());
            pending_ = line;
            have_line_ = true;
        }
        if (format_ == 0) format_ = (pending_[0] == '>') ? 1 : (pending_[0] == '@') ? 2 : 3;
        ++count_;
        if (format_ == 3) {
            seq = pending_;
            name = "read" + std::to_string(count_);
            have_line_ = false;
        } else if (format_ == 2) {
            name = header_name(pending_);
            if (!get_line(seq)) return false;
            get_line(line);  // +
            get_line(line);  // qualities
            have_line_ = false;
        } else {
            name = header_name(pending_);
            have_line_ = false;
            while (get_line(line)) {
                if (!line.empty() && line[0] == '>') {
                    pending_ = line;
                    have_line_ = true;
                    break;
                }
                seq += line;
            }
        }
        for (auto& ch : seq) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
        return true;
    }

private:
    bool get_line(std::string& line) {
        if (!std::getline(in_, line)) return false;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        return true;
    }
    static std::string header_name(const std::string& h) {
        size_t end = h.find_first_of(" \t\r", 1);
        return h.substr(1, end == std::string::npos ? std::string::npos : end - 1);
    }
    std::istream& in_;
    std::string pending_;
    bool have_line_ = false;
    int format_ = 0;  // 1 FASTA, 2 FASTQ, 3 plain
    size_t count_ = 0;
};

static Records reference(const std::string& data) {
    std::istringstream in(data);
    Reference r(in);
    Records out;
    std::string name, seq;
    while (r.next(name, seq)) out.emplace_back(name, seq);
    return out;
}

// Splits data with ReadCutter into blocks of about `target` bytes.
static bool cut(const std::string& data, size_t target, Records& out) {
    std::istringstream in(data);
    ReadCutter cutter(in);
    ReadBuffer buf;
    ReadChunk c;
    std::vector<std::string> names, seqs;
    out.clear();
    size_t numbered = 0;
    while (cutter.next(buf, c, target)) {
        if (c.format == ReadFormat::PLAIN && c.first != numbered) return false;  // reads numbered in order
        const size_t n = parse_chunk(c, names, seqs);
        numbered += n;
        for (size_t i = 0; i < n; ++i) out.emplace_back(names[i], seqs[i]);
    }
    return true;
}

static std::vector<std::pair<std::string, std::string>> inputs() {
    std::mt19937 rng(3);
    auto pick = [&](size_t n) { return size_t(rng() % n); };
    auto nl = [&]() { return pick(2) ? std::string("\n") : std::string("\r\n"); };
    std::vector<std::string> seqs;
    const size_t lens[] = {0, 1, 2, 5, 40, 90};
    for (int i = 0; i < 60; ++i) {
        std::string s;
        for (size_t j = lens[pick(6)]; j > 0; --j) s += "acgtnACGTN"[pick(10)];
        seqs.push_back(s);
    }
    std::string fa = nl() + nl(), fq, txt = nl();
    const char* suffix[] = {"", " d", "\tt", "\r"};
    const size_t widths[] = {1, 3, 60};
    for (size_t i = 0; i < seqs.size(); ++i) {
        const std::string& s = seqs[i];
        fa += ">r" + std::to_string(i) + suffix[pick(4)] + nl();
        const size_t w = widths[pick(3)];
        for (size_t j = 0; j < s.size(); j += w) fa += s.substr(j, w) + nl();
        if (pick(5) == 0) fa += nl();
        if (pick(5) == 0) fq += nl();
        fq += "@q" + std::to_string(i) + " x" + nl() + s + nl() + "+" + nl() + std::string(s.size(), '@') + nl();
        txt += s + nl();
        if (pick(5) == 0) txt += pick(2) ? nl() : nl() + nl();
    }
    fa += ">end";
    fq += "@t\nAC";
    txt += "acgt";
    // One sequence per line with no carriage returns or blank lines, which
    // ReadCutter counts a word at a time.
    std::string lf;
    for (const auto& q : seqs)
        if (!q.empty()) lf += q + "\n";
    std::string big = ">big\n" + std::string(5000, 'A') + "\n>small\nC\n";
    return {{"fasta", fa}, {"fastq", fq}, {"plain", txt}, {"plain lf", lf}, {"plain lf blank", "\n" + lf + "\n\nAC"},
            {"long record", big},
            {"blank lines only", "\n\r\n\n"}, {"fastq header only", "@only\n"}, {"empty", ""}};
}

int main() {
    size_t cases = 0, failures = 0;
    for (const auto& [name, data] : inputs()) {
        const Records ref = reference(data);
        Records got;
        for (size_t target : {size_t(1) << 26}) {
            ++cases;
            if (!cut(data, target, got) || got != ref) {
                ++failures;
                std::cerr << name << ": one block differs\n";
            }
        }
        for (size_t target = 1; target <= 300; ++target) {
            ++cases;
            if (!cut(data, target, got) || got != ref) {
                if (++failures <= 10) std::cerr << name << ": blocks of " << target << " bytes differ\n";
            }
        }
        std::cout << name << ": " << ref.size() << " records\n";
    }
    std::cout << cases << " splits, " << failures << " failures\n";
    if (failures == 0) std::cout << "All read chunking checks PASSED\n";
    return failures != 0;
}
