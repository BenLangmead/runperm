/**
 * Differential test of ReadChunker and parse_chunk.  Each input is split
 * with segment sizes from 1 byte up, so that segment boundaries fall at
 * every byte offset of its records, and with several caps on records and
 * bytes per chunk; every split must give the same records as one segment
 * holding the whole input.  ReadCutter is checked the same way, with
 * blocks of every size from 1 byte up, including its numbering of reads
 * given one per line.  The inputs cover FASTA, FASTQ and one sequence
 * per line, with CRLF line endings, blank lines, lower case, quality
 * strings of '@', and odd records at the end.  Build with the address and
 * undefined behavior sanitizers (make test-reader).
 */

#include "read_chunks.hpp"
#include <iostream>
#include <random>
#include <sstream>

using Records = std::vector<std::pair<std::string, std::string>>;

static bool split(const std::string& data, size_t read_size, size_t max_records, size_t max_bytes, Records& out) {
    std::istringstream in(data);
    ReadChunker chunker(in, read_size);
    ReadChunk c;
    std::vector<std::string> names, seqs;
    out.clear();
    while (chunker.next(c, max_records, max_bytes)) {
        const size_t n = parse_chunk(c, names, seqs);
        if (n == 0) return false;  // a chunk must hold a record
        for (size_t i = 0; i < n; ++i) out.emplace_back(names[i], seqs[i]);
    }
    return true;
}

// The same with ReadCutter, blocks of about `target` bytes.
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
        Records ref, got;
        if (!split(data, size_t(1) << 26, size_t(-1), size_t(-1), ref)) {
            std::cerr << name << ": reference split failed\n";
            return 1;
        }
        for (size_t rs = 1; rs <= 300; ++rs)
            for (size_t mr : {size_t(1), size_t(2), size_t(7), size_t(1000)})
                for (size_t mb : {size_t(1), size_t(100), size_t(1) << 24}) {
                    ++cases;
                    if (!split(data, rs, mr, mb, got) || got != ref) {
                        if (++failures <= 10)
                            std::cerr << name << ": read_size=" << rs << " max_records=" << mr << " max_bytes=" << mb
                                      << " differs\n";
                    }
                }
        for (size_t target = 1; target <= 300; ++target) {
            ++cases;
            if (!cut(data, target, got) || got != ref) {
                if (++failures <= 10) std::cerr << name << ": cutter target=" << target << " differs\n";
            }
        }
        std::cout << name << ": " << ref.size() << " records\n";
    }
    std::cout << cases << " splits, " << failures << " failures\n";
    if (failures == 0) std::cout << "All read chunking checks PASSED\n";
    return failures != 0;
}
