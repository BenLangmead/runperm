/**
 * Index for TeraMS-style matching statistics on Orbit move structures.
 *
 * Matching statistics are computed as in ms_query: the pattern is scanned
 * right to left with LF, and on a mismatch the scan repositions to the
 * nearest run of the wanted character above (pred) or below (succ).  The
 * difference is where the longest common extension (LCE) of the current row
 * with a candidate row comes from.  ms stores per-run LCP minima; this index
 * computes each LCE on demand, in one of two ways:
 *
 *  - psi: compare the pattern with the candidate row's suffix, read forward
 *    with psi (Orbit's FL), stopping at a mismatch or at the match length.
 *    Needs no LCP information.
 *  - phi: walk phi from the lower of the two rows to the upper one, taking
 *    the minimum of the LCP values crossed.  phi works on text positions, so
 *    the scan carries the text position of its current row (the toehold);
 *    LCP at a text position comes from one PLCP sample per phi interval,
 *    since inside a phi interval PLCP drops by one per position.  The
 *    toehold also gives each matching statistic an occurrence position.
 *
 * The index holds:
 *
 *  - LF over the BWT runs, runs-based (rows store lengths), with integrated
 *    columns per LF interval: PSI_INT and PSI_OFF, the FL point of the
 *    interval's head row, and PHI_INT and PHI_OFF, the phi point of that
 *    row's text position.  The FL point of a tail row is the one just before
 *    the next interval's head, and its phi point is one phi step from the
 *    next interval's head's phi point.
 *  - FL over the F runs, runs-based, with no data columns.  Its character at
 *    a position is the first character of that row's suffix.
 *  - Optionally, phi over text positions, starts-based (rows store absolute
 *    starts), with an integrated PLCP column holding PLCP at each interval's
 *    start.
 *
 * LF, FL and the PSI columns come from the run heads and lengths alone in
 * O(r) time and space.  phi also needs the LCP value at each run head and a
 * walk over all n rows.  Each structure takes its own splitting parameters.
 */

#ifndef _TMS_INDEX_HPP
#define _TMS_INDEX_HPP

#include "orbit/rlbwt.hpp"
#include "orbit/common.hpp"
#include <algorithm>
#include <array>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using uchar = orbit::uchar;
using ulint = orbit::ulint;

enum class TmsLFCols { PSI_INT, PSI_OFF, PHI_INT, PHI_OFF, COUNT };
enum class TmsPhiCols { PLCP, COUNT };

/** Splitting parameters for each move structure of a TmsIndex. */
struct TmsBuildOptions {
    // LF is unsplit by default, as in ms, so that the two tools take the
    // same LF steps and differ only in how they compute LCEs.
    orbit::split_params lf_split = orbit::NO_SPLITTING;
    // FL and phi are walked one dependent step at a time, so balancing them
    // keeps each step's fast-forward short.
    orbit::split_params fl_split = orbit::split_params{};
    orbit::split_params phi_split = orbit::split_params{};
};

class TmsIndex {
public:
    using LF = orbit::rlbwt::lf_permutation<TmsLFCols, true, false>;
    using FL = orbit::rlbwt::fl_permutation<orbit::empty_data_columns, false, false>;
    using Phi = orbit::rlbwt::phi_permutation_impl<TmsPhiCols, true, true, orbit::move_vector>;
    using LFPos = typename LF::position;
    using FLPos = typename FL::position;
    using PhiPos = typename Phi::position;
    using position = LFPos;

    TmsIndex() = default;

    /**
     * Build from run heads (index codes: orbit::TERMINATOR, orbit::SEPARATOR,
     * or a nucleotide) and run lengths.  If run_tops is given, it holds the
     * LCP value at each run's head row (0 for row 0), and the index also gets
     * phi; this walks LF over all n rows.
     */
    TmsIndex(const std::vector<uchar>& heads, const std::vector<ulint>& lens,
             const TmsBuildOptions& opts = TmsBuildOptions{}, const std::vector<ulint>* run_tops = nullptr) {
        using Enc = orbit::rlbwt::rlbwt_interval_encoding<>;
        // FL first, so that its encoding is freed before LF's is built.
        fl_ = FL(Enc::fl_interval_encoding(heads, lens, opts.fl_split));
        Enc lf_enc = Enc::lf_interval_encoding(heads, lens, opts.lf_split);
        const ulint lf_count = lf_enc.intervals();
        std::vector<typename LF::data_tuple> cols(lf_count);
        {
            // Merge the two partitions of the rows: for each LF interval,
            // the FL interval holding its head row and the offset within it.
            ulint l_pos = 0, f_pos = 0, f_int = 0;
            const ulint f_count = fl_.intervals();
            for (ulint k = 0; k < lf_count; ++k) {
                while (f_int < f_count && f_pos + fl_.get_length(f_int) <= l_pos)
                    f_pos += fl_.get_length(f_int++);
                cols[k][col(TmsLFCols::PSI_INT)] = f_int;
                cols[k][col(TmsLFCols::PSI_OFF)] = l_pos - f_pos;
                l_pos += lf_enc.get_length(k);
            }
        }
        if (run_tops) {
            if (run_tops->size() != heads.size()) throw std::invalid_argument("run_tops must have one value per run");
            build_phi(heads, lens, lf_enc, *run_tops, opts.phi_split, cols);
            has_phi_ = true;
        }
        lf_ = LF(lf_enc, cols);
        compute_occurs();
    }

    /** Write the index in Orbit's packed serialization, host byte order. */
    size_t serialize(std::ostream& out) {
        size_t bytes = 0;
        out.write(MAGIC, 4);
        const uint32_t v = VERSION, flags = has_phi_ ? 1 : 0;
        out.write(reinterpret_cast<const char*>(&v), sizeof(v));
        out.write(reinterpret_cast<const char*>(&flags), sizeof(flags));
        bytes += 12;
        bytes += lf_.serialize(out);
        bytes += fl_.serialize(out);
        if (has_phi_) bytes += phi_.serialize(out);
        return bytes;
    }

    /** Read an index written by serialize().  Throws on malformed input. */
    void load(std::istream& in) {
        char magic[4] = {};
        uint32_t v = 0, flags = 0;
        in.read(magic, 4);
        in.read(reinterpret_cast<char*>(&v), sizeof(v));
        in.read(reinterpret_cast<char*>(&flags), sizeof(flags));
        if (!in.good() || std::memcmp(magic, MAGIC, 4) != 0) throw std::runtime_error("not a tms index");
        if (v != VERSION) throw std::runtime_error("unsupported tms index version " + std::to_string(v));
        lf_.load(in);
        fl_.load(in);
        has_phi_ = (flags & 1) != 0;
        if (has_phi_) phi_.load(in);
        if (!in.good()) throw std::runtime_error("truncated tms index");
        compute_occurs();
    }

    bool has_phi() const { return has_phi_; }

    /** True if byte c occurs in the indexed text. */
    bool occurs(uchar c) const { return occurs_[c]; }

    // LF side, with the same names MSIndexSpillLCP uses.
    uchar get_character(ulint i) { return lf_.get_character(i); }
    ulint get_length(ulint i) const { return lf_.get_length(i); }
    LFPos LF_step(LFPos p) { return lf_.LF(p); }
    LFPos start_LF(LFPos p) const { return lf_.start_next(p); }
    LFPos finish_LF(LFPos p) const { return lf_.finish_next(p); }
    void prefetch(ulint i) const { lf_.prefetch(i); }
    LFPos first() { return lf_.first(); }
    LFPos last() { return lf_.last(); }
    LFPos up(LFPos p) { return lf_.up(p); }
    LFPos down(LFPos p) { return lf_.down(p); }
    ulint move_runs() const { return lf_.intervals(); }
    ulint domain() const { return lf_.domain(); }

    /** FL point of the head row of LF interval k. */
    FLPos psi_at_head(ulint k) const {
        FLPos q;
        q.interval = lf_.template get<TmsLFCols::PSI_INT>(k);
        q.offset = lf_.template get<TmsLFCols::PSI_OFF>(k);
        return q;
    }
    /** FL point of the tail row of LF interval k, which must not be the last. */
    FLPos psi_at_tail(ulint k) const {
        FLPos q = psi_at_head(k + 1);
        if (q.offset > 0) {
            --q.offset;
        } else {
            --q.interval;
            q.offset = fl_.get_length(q.interval) - 1;
        }
        return q;
    }

    // FL side.
    uchar psi_character(FLPos q) { return fl_.get_character(q.interval); }
    FLPos psi(FLPos q) { return fl_.FL(q); }
    FLPos start_psi(FLPos q) const { return fl_.start_next(q); }
    FLPos finish_psi(FLPos q) const { return fl_.finish_next(q); }
    void prefetch_psi(ulint i) const { fl_.prefetch(i); }
    ulint psi_intervals() const { return fl_.intervals(); }

    // phi side.  The text position of a phi point is its idx.
    /** phi point of the text position of LF interval k's head row. */
    PhiPos phi_at_head(ulint k) const {
        PhiPos p;
        p.interval = lf_.template get<TmsLFCols::PHI_INT>(k);
        p.offset = lf_.template get<TmsLFCols::PHI_OFF>(k);
        p.idx = 0;  // set by resolve_phi
        return p;
    }
    /** Fill in the text position of a phi point from its interval and offset. */
    PhiPos resolve_phi(PhiPos p) const { return phi_.finish_next(p); }
    PhiPos phi(PhiPos p) { return phi_.phi(p); }
    PhiPos start_phi(PhiPos p) const { return phi_.start_next(p); }
    PhiPos finish_phi(PhiPos p) const { return phi_.finish_next(p); }
    void prefetch_phi(ulint i) const { phi_.prefetch(i); }
    /** The phi point one text position to the left, wrapping at position 0. */
    PhiPos phi_left(PhiPos p) {
        if (p.offset > 0) {
            --p.offset;
            --p.idx;
            return p;
        }
        return phi_.up(p);
    }
    /** PLCP at a phi point: the LCP of its row with the row above. */
    ulint plcp(PhiPos p) const { return phi_.template get<TmsPhiCols::PLCP>(p.interval) - p.offset; }
    ulint phi_intervals() const { return phi_.intervals(); }

    const LF& lf() const { return lf_; }
    const FL& fl() const { return fl_; }

    /** One line per structure: intervals and column widths in bits. */
    void describe(std::ostream& os) const {
        auto row_bits = [](const auto& w) { size_t b = 0; for (auto x : w) b += x; return b; };
        auto widths = [&](const auto& w) { for (size_t i = 0; i < w.size(); ++i) os << (i ? "," : "") << int(w[i]); };
        const auto& w = lf_.get_widths();
        os << "LF: intervals=" << lf_.intervals() << " runs=" << lf_.runs() << " n=" << lf_.domain()
           << " widths(len,ptr,off,chr,psi_int,psi_off,phi_int,phi_off)=";
        widths(w);
        os << " row_bits=" << row_bits(w) << "\n";
        const auto& fw = fl_.get_widths();
        os << "FL: intervals=" << fl_.intervals() << " widths(len,ptr,off,chr)=";
        widths(fw);
        os << " row_bits=" << row_bits(fw) << "\n";
        double total = double(row_bits(w)) * lf_.intervals() + double(row_bits(fw)) * fl_.intervals();
        if (has_phi_) {
            const auto& pw = phi_.get_widths();
            os << "phi: intervals=" << phi_.intervals() << " widths(start,ptr,off,plcp)=";
            widths(pw);
            os << " row_bits=" << row_bits(pw) << "\n";
            total += double(row_bits(pw)) * phi_.intervals();
        }
        os << "total: " << total / 8 / 1e6 << " MB, " << total / lf_.runs() << " bits/run\n";
    }

private:
    static constexpr char MAGIC[4] = {'T', 'M', 'S', 'X'};
    static constexpr uint32_t VERSION = 2;

    LF lf_;
    FL fl_;
    Phi phi_;
    bool has_phi_ = false;
    std::array<bool, 256> occurs_{};

    static constexpr size_t col(TmsLFCols c) { return static_cast<size_t>(c); }

    void compute_occurs() {
        occurs_.fill(false);
        for (ulint i = 0; i < lf_.intervals(); ++i) occurs_[lf_.get_character(i)] = true;
    }

    /**
     * Build phi and fill the PHI columns.  Walking LF from row 0, whose text
     * position is n - 1, gives the text position of every row; the walk
     * records it at each LF interval's head and tail.  phi intervals start at
     * the text positions of true run heads (rows whose BWT character differs
     * from the row above's, and row 0), and the interval starting at SA[h]
     * maps to SA[h - 1] (SA[0] to SA[n - 1]).  The PLCP sample at SA[h] is
     * LCP[h], the top of the run whose head is h.
     */
    template <typename Enc>
    void build_phi(const std::vector<uchar>& heads, const std::vector<ulint>& lens, const Enc& lf_enc,
                   const std::vector<ulint>& run_tops, const orbit::split_params& sp,
                   std::vector<typename LF::data_tuple>& cols) {
        const ulint lf_count = lf_enc.intervals();
        std::vector<ulint> sa_head(lf_count), sa_tail(lf_count), top(lf_count, NOT_A_RUN_HEAD);
        {
            // The LCP at each LF interval that starts an original run.
            ulint k = 0;
            for (size_t r = 0; r < lens.size(); ++r) {
                top[k] = run_tops[r];
                ulint covered = 0;
                while (covered < lens[r]) covered += lf_enc.get_length(k++);
            }
            orbit::rlbwt::lf_move<false> lf(heads, lens, lf_enc.get_split_params());
            const ulint n = lf.domain();
            auto pos = lf.first();
            ulint sa = n - 1;
            std::vector<bool> seen(lf_count, false);
            for (ulint i = 0; i < n; ++i, --sa) {
                if (pos.offset == 0) { sa_head[pos.interval] = sa; seen[pos.interval] = true; }
                if (pos.offset + 1 == lf.get_length(pos.interval)) sa_tail[pos.interval] = sa;
                pos = lf.LF(pos);
            }
            if (pos != lf.first() || std::find(seen.begin(), seen.end(), false) != seen.end())
                throw std::runtime_error("LF is not a single cycle; phi needs one terminator");
        }
        // True run heads, sorted by text position.
        std::vector<ulint> order;
        for (ulint k = 0; k < lf_count; ++k)
            if (k == 0 || lf_enc.get_heads()[k] != lf_enc.get_heads()[k - 1]) order.push_back(k);
        std::sort(order.begin(), order.end(), [&](ulint a, ulint b) { return sa_head[a] < sa_head[b]; });
        const ulint n = lf_enc.domain();
        std::vector<ulint> starts(order.size()), images(order.size()), samples(order.size()), lengths(order.size());
        ulint max_length = 0;
        for (size_t x = 0; x < order.size(); ++x) {
            const ulint k = order[x];
            if (top[k] == NOT_A_RUN_HEAD) throw std::logic_error("a true run head must start an original run");
            starts[x] = sa_head[k];
            images[x] = sa_tail[k == 0 ? lf_count - 1 : k - 1];
            samples[x] = top[k];
        }
        if (starts.empty() || starts[0] != 0) throw std::logic_error("phi must start at text position 0");
        for (size_t x = 0; x < order.size(); ++x) {
            lengths[x] = (x + 1 < order.size() ? starts[x + 1] : n) - starts[x];
            max_length = std::max(max_length, lengths[x]);
        }
        auto enc = orbit::interval_encoding_impl<>::from_lengths_and_images(lengths, images, n, max_length, sp);
        // PLCP at the start of each split interval, and those starts.
        std::vector<orbit::columns_tuple<TmsPhiCols>> phi_cols(enc.intervals());
        std::vector<ulint> split_starts(enc.intervals());
        ulint s = 0;
        size_t x = 0;
        for (ulint i = 0; i < enc.intervals(); ++i) {
            while (x + 1 < starts.size() && starts[x + 1] <= s) ++x;
            const ulint into = s - starts[x];
            if (samples[x] < into) throw std::runtime_error("PLCP samples are inconsistent with the RLBWT");
            phi_cols[i][0] = samples[x] - into;
            split_starts[i] = s;
            s += enc.get_length(i);
        }
        phi_ = Phi(enc, phi_cols);
        for (ulint k = 0; k < lf_count; ++k) {
            const ulint p = sa_head[k];
            const ulint i = static_cast<ulint>(std::upper_bound(split_starts.begin(), split_starts.end(), p) -
                                               split_starts.begin()) - 1;
            cols[k][col(TmsLFCols::PHI_INT)] = i;
            cols[k][col(TmsLFCols::PHI_OFF)] = p - split_starts[i];
        }
    }

    static constexpr ulint NOT_A_RUN_HEAD = std::numeric_limits<ulint>::max();
};

#endif /* _TMS_INDEX_HPP */
