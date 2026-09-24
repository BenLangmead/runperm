/**
 * Index for TeraMS-style matching statistics on Orbit move structures.
 *
 * Matching statistics are computed as in ms_query: the pattern is scanned
 * right to left with LF, and on a mismatch the scan repositions to the
 * nearest run of the wanted character above (pred) or below (succ).  The
 * difference is where the longest common extension (LCE) of the current row
 * with a candidate row comes from.  ms stores per-run LCP minima; this index
 * stores no LCP information and instead compares the pattern with the
 * candidate row's suffix, read forward with psi (Orbit's FL).  So the index
 * holds:
 *
 *  - LF over the BWT runs, runs-based (rows store lengths), with two
 *    integrated columns per LF interval: PSI_INT and PSI_OFF, the FL point
 *    (interval, offset) of the interval's head row.  The FL point of a tail
 *    row is the one just before the next interval's head.
 *  - FL over the F runs, runs-based, with no data columns.  Its character at
 *    a position is the first character of that row's suffix.
 *
 * Both are built from the run heads and lengths alone, in O(r) time and
 * space beyond Orbit's own encodings.  Each structure takes its own
 * splitting parameters.
 */

#ifndef _TMS_INDEX_HPP
#define _TMS_INDEX_HPP

#include "orbit/rlbwt.hpp"
#include "orbit/common.hpp"
#include <array>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using uchar = orbit::uchar;
using ulint = orbit::ulint;

enum class TmsLFCols { PSI_INT, PSI_OFF, COUNT };

/** Splitting parameters for each move structure of a TmsIndex. */
struct TmsBuildOptions {
    // LF is unsplit by default, as in ms, so that the two tools take the
    // same LF steps and differ only in how they compute LCEs.
    orbit::split_params lf_split = orbit::NO_SPLITTING;
    // FL is walked one dependent step per matched character, so balancing
    // it keeps each step's fast-forward short.
    orbit::split_params fl_split = orbit::split_params{};
};

class TmsIndex {
public:
    using LF = orbit::rlbwt::lf_permutation<TmsLFCols, true, false>;
    using FL = orbit::rlbwt::fl_permutation<orbit::empty_data_columns, false, false>;
    using LFPos = typename LF::position;
    using FLPos = typename FL::position;
    using position = LFPos;

    TmsIndex() = default;

    /**
     * Build from run heads (index codes: orbit::TERMINATOR, orbit::SEPARATOR,
     * or a nucleotide) and run lengths.
     */
    TmsIndex(const std::vector<uchar>& heads, const std::vector<ulint>& lens,
             const TmsBuildOptions& opts = TmsBuildOptions{}) {
        using Enc = orbit::rlbwt::rlbwt_interval_encoding<>;
        std::vector<typename LF::data_tuple> cols;
        {
            Enc fl_enc = Enc::fl_interval_encoding(heads, lens, opts.fl_split);
            Enc lf_enc = Enc::lf_interval_encoding(heads, lens, opts.lf_split);
            // Merge the two partitions of the rows: for each LF interval,
            // the FL interval holding its head row and the offset within it.
            cols.resize(lf_enc.intervals());
            ulint l_pos = 0, f_pos = 0, f_int = 0;
            const ulint f_count = fl_enc.intervals();
            for (ulint k = 0; k < lf_enc.intervals(); ++k) {
                while (f_int < f_count && f_pos + fl_enc.get_length(f_int) <= l_pos)
                    f_pos += fl_enc.get_length(f_int++);
                cols[k][static_cast<size_t>(TmsLFCols::PSI_INT)] = f_int;
                cols[k][static_cast<size_t>(TmsLFCols::PSI_OFF)] = l_pos - f_pos;
                l_pos += lf_enc.get_length(k);
            }
            fl_ = FL(fl_enc);
            lf_ = LF(lf_enc, cols);
        }
        compute_occurs();
    }

    /** Write the index in Orbit's packed serialization, host byte order. */
    size_t serialize(std::ostream& out) {
        size_t bytes = 0;
        out.write(MAGIC, 4);
        const uint32_t v = VERSION;
        out.write(reinterpret_cast<const char*>(&v), sizeof(v));
        bytes += 8;
        bytes += lf_.serialize(out);
        bytes += fl_.serialize(out);
        return bytes;
    }

    /** Read an index written by serialize().  Throws on malformed input. */
    void load(std::istream& in) {
        char magic[4] = {};
        uint32_t v = 0;
        in.read(magic, 4);
        in.read(reinterpret_cast<char*>(&v), sizeof(v));
        if (!in.good() || std::memcmp(magic, MAGIC, 4) != 0) throw std::runtime_error("not a tms index");
        if (v != VERSION) throw std::runtime_error("unsupported tms index version " + std::to_string(v));
        lf_.load(in);
        fl_.load(in);
        if (!in.good()) throw std::runtime_error("truncated tms index");
        compute_occurs();
    }

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

    const LF& lf() const { return lf_; }
    const FL& fl() const { return fl_; }

    /** One line per structure: intervals and column widths in bits. */
    void describe(std::ostream& os) const {
        auto w = lf_.get_widths();
        size_t lf_bits = 0;
        for (auto b : w) lf_bits += b;
        os << "LF: intervals=" << lf_.intervals() << " runs=" << lf_.runs() << " n=" << lf_.domain()
           << " widths(len,ptr,off,chr,psi_int,psi_off)=";
        for (size_t i = 0; i < w.size(); ++i) os << (i ? "," : "") << int(w[i]);
        os << " row_bits=" << lf_bits << "\n";
        auto fw = fl_.get_widths();
        size_t fl_bits = 0;
        for (auto b : fw) fl_bits += b;
        os << "FL: intervals=" << fl_.intervals() << " widths(len,ptr,off,chr)=";
        for (size_t i = 0; i < fw.size(); ++i) os << (i ? "," : "") << int(fw[i]);
        os << " row_bits=" << fl_bits << "\n";
        const double total_bits = double(lf_bits) * lf_.intervals() + double(fl_bits) * fl_.intervals();
        os << "total: " << total_bits / 8 / 1e6 << " MB, " << total_bits / lf_.runs() << " bits/run\n";
    }

private:
    static constexpr char MAGIC[4] = {'T', 'M', 'S', 'X'};
    static constexpr uint32_t VERSION = 1;

    LF lf_;
    FL fl_;
    std::array<bool, 256> occurs_{};

    void compute_occurs() {
        occurs_.fill(false);
        for (ulint i = 0; i < lf_.intervals(); ++i) occurs_[lf_.get_character(i)] = true;
    }
};

#endif /* _TMS_INDEX_HPP */
