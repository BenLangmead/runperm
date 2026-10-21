// Benchmark of the ways ms query --report smem-all counts SMEM occurrences
// (see tms_count.hpp), on one load of a tms index.
//
// Usage: tms_count_bench INDEX OUT_PREFIX [options] NAME=READS.fa[:MAX_READS] ...
//   --k K            searches, walks or listings in flight (default 32)
//   --mode M         psi (default), dual or phi: how matching statistics reposition
//   --row-shift R    backward search's row-start sampling (default 4)
//   --caps a,b,...   positions listed per SMEM (default 1,100)
//   --methods LIST   comma-separated walk, backward, adaptive:B or adaptive:B:S
//                    (default walk,backward,adaptive:16,adaptive:64,adaptive:256)
//   --rounds R       rounds over the configurations, in rotated order (default 3)
//   --block-bases N  bases per query block (default 4,000,000)
//   --max-walk-steps N  skip walk at a cap once a run has shown that it would
//                    take more than N phi and phi_inv steps (default: never)
//
// For each read set, cap, method and round, the reads are queried in blocks
// as ms query does with one thread: matching statistics (with toeholds
// unless the method is backward), then the SMEM report, then formatting as
// ropebwt3 mem -p N.  OUT_PREFIX.tsv gets one row per run with the time of
// each phase (load time is reported apart from all of them), what counting
// did, and a hash of the formatted output, which must be the same for every
// method at a cap; the program fails if it is not.  backward runs first in
// the first round, so that its counts say how long the walk would take.  OUT_PREFIX.dist.tsv gets
// each set's SMEM length and count distribution, and OUT_PREFIX.hist.tsv the
// joint histogram of length and count in powers of two, and
// OUT_PREFIX.skipped.tsv the walks skipped and the steps they would take.
#include "tms_smem.hpp"
#include "tms_query.hpp"
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {
using clock_type = std::chrono::steady_clock;
double since(clock_type::time_point t) { return std::chrono::duration<double>(clock_type::now() - t).count(); }

struct ReadSet {
    std::string name;
    std::vector<std::string> names, seqs;
    size_t bases = 0;
};

ReadSet read_fasta(const std::string& name, const std::string& path, size_t max_reads) {
    ReadSet r;
    r.name = name;
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (r.names.size() == max_reads) break;
            r.names.push_back(line.substr(1, line.find_first_of(" \t") - 1));
            r.seqs.emplace_back();
        } else if (!r.seqs.empty()) {
            r.seqs.back() += line;
        }
    }
    for (const auto& s : r.seqs) r.bases += s.size();
    return r;
}

struct Method {
    std::string label;
    TmsCountOptions opts;
};

Method parse_method(const std::string& m) {
    Method x;
    x.label = m;
    if (m == "walk") x.opts.method = TmsCount::WALK;
    else if (m == "backward") x.opts.method = TmsCount::BACKWARD;
    else if (m.rfind("adaptive", 0) == 0) {
        x.opts.method = TmsCount::ADAPTIVE;
        std::stringstream ss(m);
        std::string part;
        std::getline(ss, part, ':');
        if (std::getline(ss, part, ':')) x.opts.budget = std::stoull(part);
        if (std::getline(ss, part, ':')) x.opts.short_len = std::stoull(part);
    } else {
        throw std::invalid_argument("unknown method " + m);
    }
    return x;
}

std::vector<std::string> split(const std::string& s, char sep) {
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string part;
    while (std::getline(ss, part, sep)) out.push_back(part);
    return out;
}

uint64_t fnv(uint64_t h, const std::string& s) {
    for (unsigned char c : s) h = (h ^ c) * 1099511628211ULL;
    return h;
}

ulint quantile(std::vector<ulint>& v, double q) {
    if (v.empty()) return 0;
    const size_t at = std::min(v.size() - 1, static_cast<size_t>(q * double(v.size())));
    std::nth_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(at), v.end());
    return v[at];
}

unsigned log2_bucket(ulint x) {
    unsigned b = 0;
    while ((ulint(2) << b) <= x) ++b;
    return b;
}
}  // namespace

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "usage: tms_count_bench INDEX OUT_PREFIX [options] NAME=READS.fa[:MAX_READS] ...\n";
        return 1;
    }
    const std::string idx_path = argv[1], out_prefix = argv[2];
    size_t k = 32, rounds = 3, block_bases = 4000000;
    ulint max_walk_steps = ~ulint(0);
    unsigned row_shift = TmsCounter::DEFAULT_SHIFT;
    TmsMode mode = TmsMode::PSI;
    std::vector<ulint> caps{1, 100};
    std::vector<std::string> method_names{"walk", "backward", "adaptive:16", "adaptive:64", "adaptive:256"};
    std::vector<std::pair<std::string, std::string>> sets;
    for (int i = 3; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--k" && i + 1 < argc) k = std::stoul(argv[++i]);
        else if (a == "--rounds" && i + 1 < argc) rounds = std::stoul(argv[++i]);
        else if (a == "--row-shift" && i + 1 < argc) row_shift = static_cast<unsigned>(std::stoul(argv[++i]));
        else if (a == "--block-bases" && i + 1 < argc) block_bases = std::stoul(argv[++i]);
        else if (a == "--max-walk-steps" && i + 1 < argc) max_walk_steps = std::stoull(argv[++i]);
        else if (a == "--mode" && i + 1 < argc) {
            const std::string m = argv[++i];
            mode = m == "psi" ? TmsMode::PSI : m == "dual" ? TmsMode::DUAL : m == "phi" ? TmsMode::PHI
                                                                                       : throw std::invalid_argument(m);
        } else if (a == "--caps" && i + 1 < argc) {
            caps.clear();
            for (const auto& c : split(argv[++i], ',')) caps.push_back(std::stoull(c));
        } else if (a == "--methods" && i + 1 < argc) {
            method_names = split(argv[++i], ',');
        } else {
            const size_t eq = a.find('=');
            if (eq == std::string::npos) {
                std::cerr << "unknown argument " << a << "\n";
                return 1;
            }
            sets.emplace_back(a.substr(0, eq), a.substr(eq + 1));
        }
    }
    std::vector<Method> methods;
    for (const auto& m : method_names) methods.push_back(parse_method(m));
    std::stable_partition(methods.begin(), methods.end(),
                          [](const Method& m) { return m.opts.method == TmsCount::BACKWARD; });

    // Load, timed apart from everything else.
    const auto t_load = clock_type::now();
    TmsIndex idx;
    {
        std::ifstream in(idx_path, std::ios::binary);
        if (!in) {
            std::cerr << "cannot open " << idx_path << "\n";
            return 1;
        }
        const TmsParts need = tms_query_parts(mode, false, TmsReport::SMEM_ALL);
        idx.load(in, need);
    }
    const double load_s = since(t_load);
    const auto t_ctr = clock_type::now();
    const TmsCounter counter(idx, row_shift);
    const double counter_s = since(t_ctr);
    std::cerr << "load_s=" << load_s << " counter_s=" << counter_s << " counter_bytes=" << counter.bytes()
              << " lf_intervals=" << idx.move_runs() << " n=" << idx.domain() << "\n";

    std::ofstream tsv(out_prefix + ".tsv"), dist(out_prefix + ".dist.tsv"), hist(out_prefix + ".hist.tsv"),
        skipped(out_prefix + ".skipped.tsv");
    tsv << "set\treads\tbases\tmethod\tbudget\tshort_len\tcap\tround\tms_s\treport_s\twalk_s\tsearch_s\tlist_s"
           "\tformat_s\tsmems\twalked\tswitched\tsearched\twalk_steps\tsearch_steps\ttop_moves\tbot_moves"
           "\tscan_rows\tlf_ff\tfar_ends\tsum_rows\tlist_steps\tlisted\trepositions\thash\n";
    dist << "set\treads\tbases\tsmems\tsmems_per_kb\trepositions_per_base\tlen_mean\tlen_q10\tlen_q25\tlen_q50"
            "\tlen_q75\tlen_q90\tlen_max\tcount_mean\tcount_q10\tcount_q25\tcount_q50\tcount_q75\tcount_q90"
            "\tcount_q99\tcount_max\tsum_count\tsum_len\n";
    hist << "set\tlen_log2\tcount_log2\tsmems\n";
    skipped << "set\tcap\twalk_steps\n";
    const SeqIdTable* ids = idx.seq_ids();

    for (const auto& [name, spec] : sets) {
        std::string path = spec;
        size_t max_reads = SIZE_MAX;
        const size_t colon = spec.rfind(':');
        if (colon != std::string::npos && spec.find('/', colon) == std::string::npos) {
            path = spec.substr(0, colon);
            max_reads = std::stoull(spec.substr(colon + 1));
        }
        const ReadSet rs = read_fasta(name, path, max_reads);
        std::cerr << "set " << name << ": " << rs.seqs.size() << " reads, " << rs.bases << " bases\n";
        // Blocks of reads, by bases.
        std::vector<std::pair<size_t, size_t>> blocks;
        for (size_t a = 0; a < rs.seqs.size();) {
            size_t b = a, bases = 0;
            while (b < rs.seqs.size() && (b == a || bases + rs.seqs[b].size() <= block_bases)) bases += rs.seqs[b++].size();
            blocks.emplace_back(a, b);
            a = b;
        }
        bool have_dist = false;
        for (ulint cap : caps) {
            std::map<std::string, uint64_t> hashes;
            // The walk's steps, from any run's counts: each SMEM's count less one.
            ulint walk_steps = 0;
            bool walk_known = false;
            for (size_t round = 0; round < rounds; ++round)
                for (size_t mi = 0; mi < methods.size(); ++mi) {
                    const Method& m = methods[(mi + round) % methods.size()];
                    if (m.opts.method == TmsCount::WALK && walk_known && walk_steps > max_walk_steps) {
                        if (round == 0)
                            std::cerr << name << " cap " << cap << ": walk skipped, it would take " << walk_steps
                                      << " steps\n";
                        if (round == 0) skipped << name << '\t' << cap << '\t' << walk_steps << '\n';
                        continue;
                    }
                    tms_count_stats = TmsCountStats{};
#ifdef TMS_STATS
                    tms_stats = TmsStats{};
#endif
                    const bool toeholds = m.opts.method != TmsCount::BACKWARD;
                    double ms_s = 0, report_s = 0, format_s = 0;
                    ulint listed = 0;
                    uint64_t h = 1469598103934665603ULL;
                    std::vector<ulint> lens, counts;
                    std::map<std::pair<unsigned, unsigned>, ulint> joint;
                    ulint sum_count = 0, sum_len = 0, run_steps = 0;
                    std::vector<std::string> block;
                    std::vector<std::vector<ulint>> ms, pos;
                    std::vector<TmsSmemHits> hits;
                    std::string text;
                    for (const auto& [a, b] : blocks) {
                        block.assign(rs.seqs.begin() + static_cast<std::ptrdiff_t>(a),
                                     rs.seqs.begin() + static_cast<std::ptrdiff_t>(b));
                        auto t = clock_type::now();
                        tms_query_batch(idx, block, k, ms, mode, toeholds ? &pos : nullptr, true);
                        ms_s += since(t);
                        if (!toeholds) pos.clear();
                        t = clock_type::now();
                        const TmsCounting counting{m.opts, &counter, &block};
                        tms_report_smems_batch(idx, ms, pos, 0, TmsReport::SMEM_ALL, k, hits, cap, &counting);
                        report_s += since(t);
                        t = clock_type::now();
                        text.clear();
                        for (size_t j = 0; j < block.size(); ++j)
                            tms_format_smems_rope(rs.names[a + j], hits[j], ids, text, &ms[j]);
                        format_s += since(t);
                        h = fnv(h, text);
                        for (const auto& hh : hits)
                            for (const auto& s : hh.smems) {
                                listed += s.listed;
                                run_steps += s.count - 1;
                                if (!have_dist) {
                                    lens.push_back(s.len);
                                    counts.push_back(s.count);
                                    sum_count += s.count;
                                    sum_len += s.len;
                                    ++joint[{log2_bucket(s.len), log2_bucket(s.count)}];
                                }
                            }
                    }
                    walk_steps = run_steps;
                    walk_known = true;
                    ulint reps = 0;
#ifdef TMS_STATS
                    reps = tms_stats.repositions;
#endif
                    const TmsCountStats& c = tms_count_stats;
                    tsv << name << '\t' << rs.seqs.size() << '\t' << rs.bases << '\t' << m.label << '\t'
                        << m.opts.budget << '\t' << m.opts.short_len << '\t' << cap << '\t' << round << '\t' << ms_s
                        << '\t' << report_s << '\t' << c.walk_s << '\t' << c.search_s << '\t' << c.list_s << '\t'
                        << format_s << '\t' << c.smems << '\t' << c.walked << '\t' << c.switched << '\t'
                        << c.searched << '\t' << c.walk_steps << '\t' << c.search_steps << '\t' << c.top_moves
                        << '\t' << c.bot_moves << '\t' << c.scan_rows << '\t' << c.lf_ff << '\t' << c.far_ends
                        << '\t' << c.sum_rows << '\t' << c.list_steps << '\t' << listed << '\t' << reps << '\t'
                        << std::hex << h << std::dec << '\n';
                    tsv.flush();
                    std::cerr << name << " cap " << cap << " " << m.label << " round " << round << ": ms_s=" << ms_s
                              << " report_s=" << report_s << " format_s=" << format_s << "\n";
                    if (!have_dist) {
                        have_dist = true;
                        const size_t n = lens.size();
                        dist << name << '\t' << rs.seqs.size() << '\t' << rs.bases << '\t' << n << '\t'
                             << 1000.0 * double(n) / double(rs.bases) << '\t' << double(reps) / double(rs.bases)
                             << '\t' << (n ? double(sum_len) / double(n) : 0.0);
                        for (double q : {0.10, 0.25, 0.50, 0.75, 0.90}) dist << '\t' << quantile(lens, q);
                        dist << '\t' << (n ? *std::max_element(lens.begin(), lens.end()) : 0) << '\t'
                             << (n ? double(sum_count) / double(n) : 0.0);
                        for (double q : {0.10, 0.25, 0.50, 0.75, 0.90, 0.99}) dist << '\t' << quantile(counts, q);
                        dist << '\t' << (n ? *std::max_element(counts.begin(), counts.end()) : 0) << '\t'
                             << sum_count << '\t' << sum_len << '\n';
                        dist.flush();
                        for (const auto& [key, v] : joint)
                            hist << name << '\t' << key.first << '\t' << key.second << '\t' << v << '\n';
                        hist.flush();
                    }
                    auto [it, fresh] = hashes.emplace(m.label, h);
                    if (!fresh && it->second != h) {
                        std::cerr << "FAILED: " << m.label << " output differs between rounds\n";
                        return 2;
                    }
                    for (const auto& [other, oh] : hashes)
                        if (oh != h) {
                            std::cerr << "FAILED: " << name << " cap " << cap << ": " << m.label << " and " << other
                                      << " give different output\n";
                            return 2;
                        }
                }
        }
    }
    std::cerr << "all methods agree\n";
    return 0;
}
