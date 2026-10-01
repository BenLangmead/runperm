/**
 * Streaming build of an MSIndexSpillLCP from run heads, run lengths and a
 * source of each run's retained LCP pairs, such as a TeraLCP minima file.
 *
 * The result is the index that apply_lf_splitting, build_spill_data_from_pairs
 * and the MSIndexSpillLCP constructor make (see ms_rlbwt.hpp), byte for byte,
 * built with memory near the size of the index itself:
 *
 *  - The runs' pairs are never all held.  The source is read twice, a run
 *    at a time, and each run is cut into its LF split rows on the fly.  The
 *    first pass counts the bit widths of the rows' top and min-sub LCPs,
 *    which is all the percentiles that set the column widths need, and
 *    bounds the spillover offsets.  The second writes each row's columns
 *    into a packed table and its spillover record into its spillover
 *    array.
 *  - Coalescing finds equal records with a hash table of offsets into the
 *    spillover arrays, so a record's bytes are stored once.
 *  - The move structure is built with its data columns written in place
 *    (Orbit's in_place_data), with no table of run data beside it.
 *
 * Author: Ben Langmead (ben.langmead@gmail.com)
 */

#ifndef _MS_BUILD_HPP
#define _MS_BUILD_HPP

#include "ms_rlbwt.hpp"
#include "mem_report.hpp"
#include "rlbwt_io.hpp"
#include <array>
#include <cstring>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace ms_build {

/** Build settings, as for build_spill_data_from_pairs and apply_lf_splitting. */
struct Options {
    double percentile_k = 0.98;
    bool coalesce = false;
    ulint spill_align = 0;
    uchar spill_split_bits = 0;
    bool minima_only = false;
    orbit::split_params lf_split = orbit::split_params{};
    // Print memory at the phase boundaries (see mem_report).
    bool report = false;
};

/** Counts from a build. */
struct Stats {
    size_t runs = 0, rows = 0, skinny = 0, jumbo = 0;
};

/**
 * The lengths of the rows into which Orbit's LF splitting sp cuts an RLBWT,
 * in row order, packed (see lf_split_row_lengths).
 */
template <typename Heads, typename Lens>
inline orbit::int_vector split_row_lengths(const Heads& heads, const Lens& lens, const orbit::split_params& sp) {
    using Enc = orbit::rlbwt::rlbwt_interval_encoding<>;
    const Enc enc = Enc::lf_interval_encoding(heads, lens, sp);
    ulint max_len = 0;
    for (ulint k = 0; k < enc.intervals(); ++k) max_len = std::max(max_len, enc.get_length(k));
    orbit::int_vector rows(enc.intervals(), orbit::bit_width(max_len));
    for (ulint k = 0; k < enc.intervals(); ++k) rows.set(k, enc.get_length(k));
    return rows;
}

/** Reused buffers of RowCutter. */
struct CutScratch {
    std::vector<ulint> starts, tops, suffix_min;
    RunLcpPairs row;
};

/**
 * Cut one run, of length len with retained pairs p, into its LF split rows,
 * whose lengths are rows[k] onward, exactly as apply_lf_splitting does, and
 * call emit(row_len, row_pairs) for each in order.  next_top is the next
 * run's top, if there is one.  Returns the index of the run's last row plus
 * one.  p may be left changed.
 */
template <typename Emit>
inline size_t cut_run(RunLcpPairs& p, ulint len, const orbit::int_vector& rows, size_t k,
                      std::optional<ulint> next_top, bool minima_only, CutScratch& s, Emit&& emit) {
    auto& starts = s.starts;
    starts.clear();
    ulint covered = 0;
    while (covered < len) {
        if (k >= rows.size()) throw std::logic_error("LF split rows do not cover the runs");
        starts.push_back(covered);
        covered += rows[k++];
    }
    if (covered != len) throw std::logic_error("an LF split row crosses a run boundary");
    if (starts.size() == 1) {
        emit(len, p);
        return k;
    }
    // As in apply_lf_splitting: suffix_min[j] is the minimum of the values of
    // pairs j onward, and each row after the first gets the top max(u, d).
    auto& suffix_min = s.suffix_min;
    suffix_min.assign(p.size() + 1, LCP_GAP);
    for (size_t j = p.size(); j-- > 0;) suffix_min[j] = std::min(suffix_min[j + 1], p[j].second);
    auto& tops = s.tops;
    tops.assign(starts.size(), 0);
    tops[0] = p[0].second;
    ulint prefix_min = p[0].second;
    size_t j = 1;
    for (size_t t = 1; t < starts.size(); ++t) {
        size_t first_at = j;
        while (first_at < p.size() && p[first_at].first < starts[t]) ++first_at;
        for (; j < p.size() && p[j].first <= starts[t]; ++j) prefix_min = std::min(prefix_min, p[j].second);
        ulint d = suffix_min[first_at];
        if (next_top) d = std::min(d, *next_top);
        tops[t] = d == LCP_GAP ? prefix_min : std::max(prefix_min, d);
    }
    size_t at = 1;
    auto& row = s.row;
    for (size_t t = 0; t < starts.size(); ++t) {
        const ulint lo = starts[t], hi = t + 1 < starts.size() ? starts[t + 1] : len;
        row.clear();
        row.emplace_back(0, tops[t]);
        for (; at < p.size() && p[at].first < hi; ++at)
            if (p[at].first > lo) row.emplace_back(p[at].first - lo, p[at].second);
        prune_lcp_pairs(row, t + 1 < starts.size() ? std::optional<ulint>(tops[t + 1]) : next_top, minima_only);
        emit(hi - lo, row);
    }
    return k;
}

/**
 * Call f(row, run, row_pairs) for every LF split row in order, where row is
 * its index, run the index of its run and row_pairs its retained pairs (see
 * apply_lf_splitting).  source(g) must call g(pairs) for each run's
 * retained pairs in run order; g may change the pairs it is given.
 */
template <typename Lens, typename Source, typename F>
inline void for_each_row(const Lens& lens, const orbit::int_vector& rows, bool minima_only, Source&& source, F&& f) {
    RunLcpPairs cur;
    bool have = false;
    size_t run = 0, k = 0, row = 0;
    CutScratch scratch;
    auto process = [&](std::optional<ulint> next_top) {
        if (run >= lens.size()) throw std::invalid_argument("the LCP source has more runs than the RLBWT");
        k = cut_run(cur, lens[run], rows, k, next_top, minima_only, scratch,
                    [&](ulint, RunLcpPairs& rp) { f(row++, run, rp); });
        ++run;
    };
    source([&](RunLcpPairs& pairs) {
        if (pairs.empty()) throw std::invalid_argument("a run's LCP pairs must hold its top");
        if (have) process(pairs[0].second);
        std::swap(cur, pairs);
        have = true;
    });
    if (have) process(std::nullopt);
    if (run != lens.size()) throw std::invalid_argument("heads, lens and LCP pairs must have one entry per run");
    if (k != rows.size()) throw std::logic_error("LF split rows do not match the runs");
}

/** A row's top and min-sub LCP values, as build_spill_data_from_pairs computes them. */
inline std::pair<ulint, ulint> top_and_sub(const RunLcpPairs& pairs) {
    if (pairs.size() == 1) return {pairs[0].second, 0};
    ulint min_s = LCP_GAP;
    for (size_t j = 1; j < pairs.size(); ++j) min_s = std::min(min_s, pairs[j].second);
    const ulint sub = (min_s != LCP_GAP && min_s < pairs[0].second) ? (pairs[0].second - min_s) : 0;
    return {pairs[0].second, sub};
}

/** Bytes of v in ULEB128. */
inline size_t uleb128_size(ulint v) {
    size_t b = 1;
    while (v >>= 7) ++b;
    return b;
}

/**
 * The bit width of the k-th percentile of the values whose bit widths
 * (detail::ceil_log2) have the histogram hist, the same as
 * detail::ceil_log2(detail::percentile(values, k)): ceil_log2 does not
 * decrease, so the k-th smallest value's width is the k-th smallest width.
 */
inline uchar percentile_width(const std::array<size_t, 65>& hist, size_t count, double k) {
    if (count == 0) return detail::ceil_log2(0);
    size_t idx = static_cast<size_t>(k * count);
    if (idx >= count) idx = count - 1;
    size_t seen = 0;
    for (size_t b = 0; b < hist.size(); ++b) {
        seen += hist[b];
        if (seen > idx) return static_cast<uchar>(b);
    }
    throw std::logic_error("percentile histogram is short");
}

/**
 * A set of spillover records, each stored once in a spillover array, for
 * coalescing: an open-addressing hash table of the records' offsets.  A
 * slot holds 0 when empty, else 1 + (offset << 1 | jumbo), and the record's
 * length is decoded from its bytes, which are self-delimiting given whether
 * it is a jumbo record (top, row minimum, count, pairs) or not (count,
 * pairs).  Records are equal when their bytes are, whatever their kinds.
 */
class RecordSet {
public:
    /**
     * The offset of a stored record equal to bytes rec, or of rec appended
     * to spill by append(), which returns its offset, if there is none.
     */
    template <typename Append>
    size_t find_or_add(const SpilloverVector& spill, const SpilloverVector& rec, bool jumbo, Append&& append) {
        if (slots_.empty() || 10 * (count_ + 1) > 7 * slots_.size()) grow(spill);
        const ulint h = hash(rec.data(), rec.size());
        const size_t mask = slots_.size() - 1;
        for (size_t at = h & mask;; at = (at + 1) & mask) {
            const ulint slot = slots_[at];
            if (slot == 0) {
                const size_t off = append();
                slots_[at] = 1 + ((static_cast<ulint>(off) << 1) | (jumbo ? 1 : 0));
                ++count_;
                return off;
            }
            const size_t off = static_cast<size_t>((slot - 1) >> 1);
            const size_t len = record_size(spill, off, ((slot - 1) & 1) != 0);
            if (len == rec.size() && std::memcmp(spill.data() + off, rec.data(), len) == 0) return off;
        }
    }

    /** Free the table. */
    void clear() { std::vector<ulint>().swap(slots_); count_ = 0; }

private:
    std::vector<ulint> slots_;
    size_t count_ = 0;

    static ulint hash(const uchar* p, size_t len) {
        ulint h = 0x9E3779B97F4A7C15ULL ^ len;
        auto mix = [&](ulint w) {
            h ^= w;
            h *= 0xBF58476D1CE4E5B9ULL;
            h ^= h >> 31;
        };
        size_t i = 0;
        for (; i + 8 <= len; i += 8) {
            ulint w;
            std::memcpy(&w, p + i, 8);
            mix(w);
        }
        ulint w = 0;
        std::memcpy(&w, p + i, len - i);
        mix(w);
        h *= 0x94D049BB133111EBULL;
        return h ^ (h >> 29);
    }

    static size_t record_size(const SpilloverVector& spill, size_t off, bool jumbo) {
        size_t pos = off;
        if (jumbo) pos = skip_uleb128(spill, skip_uleb128(spill, pos));
        auto [count, after] = decode_uleb128(spill, pos);
        pos = after;
        for (ulint j = 0; j < 2 * count; ++j) pos = skip_uleb128(spill, pos);
        return pos - off;
    }

    void grow(const SpilloverVector& spill) {
        std::vector<ulint> old;
        old.swap(slots_);
        slots_.assign(old.empty() ? 1024 : 2 * old.size(), 0);
        const size_t mask = slots_.size() - 1;
        for (ulint slot : old) {
            if (slot == 0) continue;
            const size_t off = static_cast<size_t>((slot - 1) >> 1);
            const size_t len = record_size(spill, off, ((slot - 1) & 1) != 0);
            size_t at = hash(spill.data() + off, len) & mask;
            while (slots_[at] != 0) at = (at + 1) & mask;
            slots_[at] = slot;
        }
    }
};

/**
 * Build the index.  heads and lens are the runs' (as from load_rlbwt);
 * source(g) calls g(pairs) with each run's retained pairs in run order (see
 * for_each_row) and is called twice.  heads and lens are freed on the way.
 */
template <bool SP = false, typename Source>
inline MSIndexSpillLCP<SP> build(std::vector<uchar>&& heads_in, std::vector<ulint>&& lens_in, Source&& source,
                                    const Options& o, Stats* stats = nullptr) {
    check_lf_split(o.lf_split);
    std::vector<uchar> heads = std::move(heads_in);
    const size_t runs = heads.size();
    if (lens_in.size() != runs) throw std::invalid_argument("heads and lens must have one entry per run");
    auto phase = [&](const char* name) { if (o.report) mem_report::phase(name, runs); };
    orbit::int_vector rows = split_row_lengths(heads, lens_in, o.lf_split);
    // The run lengths, packed.
    orbit::int_vector lens;
    {
        ulint max_len = 0;
        for (ulint l : lens_in) max_len = std::max(max_len, l);
        lens = orbit::int_vector(runs, orbit::bit_width(max_len));
        for (size_t i = 0; i < runs; ++i) lens.set(i, lens_in[i]);
        std::vector<ulint>().swap(lens_in);
    }
    const size_t R = rows.size();
    phase("ms_split_rows");

    const size_t align_x = o.spill_align > 0 ? static_cast<size_t>(o.spill_align) : 1;
    const size_t num_arrays = o.spill_split_bits > 0 ? (size_t{1} << o.spill_split_bits) : 1;
    const ulint bucket_mask = o.spill_split_bits > 0 ? (ulint{1} << o.spill_split_bits) - 1 : 0;

    // Pass 1: the widths of the rows' tops and min-subs, and a bound on each
    // spillover array's size (every row's record as a jumbo one, padded).
    std::array<size_t, 65> hist_top{}, hist_sub{};
    std::vector<size_t> spill_bound(num_arrays, align_x);
    for_each_row(lens, rows, o.minima_only, source, [&](size_t i, size_t, const RunLcpPairs& pairs) {
        const auto [top, sub] = top_and_sub(pairs);
        ++hist_top[detail::ceil_log2(top)];
        ++hist_sub[detail::ceil_log2(sub)];
        ulint row_min = LCP_GAP;
        size_t bytes = 0;
        for (size_t j = 0; j < pairs.size(); ++j) {
            row_min = std::min(row_min, pairs[j].second);
            if (j > 0) bytes += uleb128_size(pairs[j].first) + uleb128_size(pairs[j].second);
        }
        bytes += uleb128_size(top) + uleb128_size(row_min) + uleb128_size(pairs.size() - 1);
        spill_bound[i & bucket_mask] += bytes + align_x - 1;
    });
    const uchar w_top = percentile_width(hist_top, R, o.percentile_k);
    const uchar w_sub = percentile_width(hist_sub, R, o.percentile_k);
    const ulint max_top = max_for_bits(w_top), max_sub = max_for_bits(w_sub);
    ulint max_stored_spill = 0;
    for (size_t b : spill_bound) max_stored_spill = std::max<ulint>(max_stored_spill, b / align_x);
    phase("ms_pass1");

    auto is_skinny = [max_top, max_sub](ulint t, ulint s) {
        if (t > max_top || s > max_sub) return false;
        return (t < max_top || s < max_sub);
    };

    // Pass 2: each row's columns, in a packed table, and its spillover record.
    enum { TOP, SUB, SPILL };
    orbit::packed_matrix<3> cols(R, {w_top, w_sub, orbit::bit_width(max_stored_spill)});
    std::array<ulint, 3> col_max{};
    std::vector<SpilloverVector> spill(num_arrays);
    for (auto& v : spill) {
        v.reserve(4096);
        v.assign(align_x, 0);
    }
    std::vector<RecordSet> sets(o.coalesce ? num_arrays : 0);
    SpilloverVector rec;
    size_t skinny = 0, jumbo = 0;
    auto store = [&](size_t i, ulint top, ulint sub, ulint so) {
        cols.set<TOP>(i, top);
        cols.set<SUB>(i, sub);
        cols.set<SPILL>(i, so);
        col_max[TOP] = std::max(col_max[TOP], top);
        col_max[SUB] = std::max(col_max[SUB], sub);
        col_max[SPILL] = std::max(col_max[SPILL], so);
    };
    // The stored offset of record rec in array b, appended or coalesced.
    auto place = [&](size_t b, bool is_jumbo) -> ulint {
        auto& sv = spill[b];
        auto append = [&]() {
            while (sv.size() % align_x != 0) sv.push_back(0);
            const size_t off = sv.size();
            sv.insert(sv.end(), rec.begin(), rec.end());
            return off;
        };
        const size_t off = o.coalesce ? sets[b].find_or_add(sv, rec, is_jumbo, append) : append();
        return static_cast<ulint>(off / align_x);
    };
    for_each_row(lens, rows, o.minima_only, source, [&](size_t i, size_t, const RunLcpPairs& pairs) {
        const size_t b = static_cast<size_t>(i & bucket_mask);
        const auto [top, sub] = top_and_sub(pairs);
        if (is_skinny(top, sub)) {
            ulint so = NO_SPILL;
            if (pairs.size() > 1) {
                rec.clear();
                append_uleb128(rec, static_cast<ulint>(pairs.size() - 1));
                for (size_t j = 1; j < pairs.size(); ++j) {
                    append_uleb128(rec, pairs[j].first);
                    append_uleb128(rec, pairs[j].second);
                }
                so = place(b, false);
            }
            store(i, top, sub, so);
            ++skinny;
        } else {
            ulint row_min = LCP_GAP;
            for (const auto& pr : pairs) row_min = std::min(row_min, pr.second);
            rec.clear();
            append_uleb128(rec, top);
            append_uleb128(rec, row_min);
            append_uleb128(rec, static_cast<ulint>(pairs.size() - 1));
            for (size_t j = 1; j < pairs.size(); ++j) {
                append_uleb128(rec, pairs[j].first);
                append_uleb128(rec, pairs[j].second);
            }
            store(i, max_top, max_sub, place(b, true));
            ++jumbo;
        }
    });
    std::vector<RecordSet>().swap(sets);
    for (auto& v : spill) v.shrink_to_fit();
    SpilloverVector().swap(rec);
    phase("ms_pass2");

    // The rows' heads, and their encoding as runs of their own.
    using Impl = typename MSIndexSpillLCP<SP>::Impl;
    using Enc = orbit::rlbwt::rlbwt_interval_encoding<>;
    Enc enc;
    {
        std::vector<uchar> row_heads(R);
        for (size_t r = 0, k = 0; r < runs; ++r)
            for (ulint covered = 0, len = lens[r]; covered < len; covered += rows[k++]) row_heads[k] = heads[r];
        std::vector<uchar>().swap(heads);
        lens = orbit::int_vector();
        enc = Enc::lf_interval_encoding(row_heads, rows, orbit::NO_SPLITTING);
        rows = orbit::int_vector();
    }
    phase("ms_encoding");
    std::array<uchar, 3> widths{};
    for (size_t c = 0; c < 3; ++c) widths[c] = orbit::bit_width(col_max[c]);
    Impl impl(enc, widths, Impl::in_place_data);
    enc = Enc();
    for (size_t i = 0; i < R; ++i) {
        impl.template set_data<LCPSpillRunCols::LCP_TOP>(i, cols.get<TOP>(i));
        impl.template set_data<LCPSpillRunCols::LCP_MIN_SUB>(i, cols.get<SUB>(i));
        impl.template set_data<LCPSpillRunCols::LCP_SPILL>(i, cols.get<SPILL>(i));
    }
    cols = orbit::packed_matrix<3>();
    phase("ms_assemble");
    if (stats) *stats = Stats{runs, R, skinny, jumbo};
    return MSIndexSpillLCP<SP>(std::move(impl), std::move(spill), max_top, max_sub, o.spill_align, o.spill_split_bits);
}

/** A source (see build) over runs' pairs held in memory. */
inline auto pairs_source(const std::vector<RunLcpPairs>& all) {
    return [&all](auto&& g) {
        RunLcpPairs p;
        for (const auto& run : all) {
            p = run;
            g(p);
        }
    };
}

/**
 * A source (see build) over each run's full LCP vector (element 0 its top),
 * reduced to its retained pairs as retained_lcp_pairs does.
 */
inline auto lcp_vectors_source(const std::vector<std::vector<ulint>>& lcps, bool minima_only) {
    return [&lcps, minima_only](auto&& g) {
        std::vector<ulint> full;
        RunLcpPairs p;
        for (size_t i = 0; i < lcps.size(); ++i) {
            compress_lcps(lcps[i], i + 1 < lcps.size() ? &lcps[i + 1] : nullptr, full, minima_only);
            p = detail::compressed_to_pairs(full);
            g(p);
        }
    };
}

/** A source (see build) over a TeraLCP minima file; throws std::runtime_error on a bad file. */
inline auto minima_source(const std::string& path, size_t runs, ulint n) {
    return [path, runs, n](auto&& g) {
        rlbwt_io::MinimaReader rd;
        std::string err;
        if (!rd.open(path, runs, n, err)) throw std::runtime_error(err);
        RunLcpPairs p;
        for (size_t i = 0; i < runs; ++i) {
            if (!rd.next(p, err)) throw std::runtime_error(err);
            g(p);
        }
        if (!rd.finish(err)) throw std::runtime_error(err);
    };
}

/**
 * A source (see build) over a file of one 64-bit LCP per row (see
 * for_each_run_lcps), reduced to each run's retained pairs; lens are the
 * runs' lengths.  Throws std::runtime_error on a bad file.
 */
template <typename Lens>
inline auto lcp_bin_source(const std::string& path, const Lens& lens, bool minima_only) {
    return [path, &lens, minima_only](auto&& g) {
        std::vector<ulint> prev, kept;
        RunLcpPairs p;
        auto finish = [&](const std::vector<ulint>* next) {
            compress_lcps(prev, next, kept, minima_only);
            p = detail::compressed_to_pairs(kept);
            g(p);
        };
        std::string err;
        const bool ok = rlbwt_io::for_each_run_lcps(path, lens, [&](size_t i, const std::vector<ulint>& v) {
            if (i > 0) finish(&v);
            prev = v;
        }, err);
        if (!ok) throw std::runtime_error(err);
        if (!lens.empty()) finish(nullptr);
    };
}

}  // namespace ms_build

#endif /* _MS_BUILD_HPP */
