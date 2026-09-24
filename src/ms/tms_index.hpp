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
 *    interval's tail row, and PHI_INT and PHI_OFF, the phi point of its head
 *    row's text position.  The FL point of a head row is the one just after
 *    the previous interval's tail, and the phi point of a tail row is one phi
 *    step from the next interval's head's phi point.  So a reposition reads
 *    both candidates' start points from rows its scan has already read.
 *  - FL over the F runs, runs-based, with no data columns.  Its character at
 *    a position is the first character of that row's suffix.
 *  - Optionally, phi over text positions, starts-based (rows store absolute
 *    starts), with an integrated PLCP column holding PLCP at each interval's
 *    start.
 *  - Optionally, with phi, phi_inv over text positions (SA[j] to SA[j + 1]),
 *    starts-based, with an integrated PLCPB column: the LCP of the row with
 *    the row below, at each interval's start.  Inside a phi_inv interval
 *    PLCPB drops by one per position, as PLCP does inside a phi interval.
 *    phi and phi_inv together enumerate every row of a BWT interval from one
 *    of its text positions (see tms_smem.hpp).
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
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using uchar = orbit::uchar;
using ulint = orbit::ulint;

enum class TmsLFCols { PSI_INT, PSI_OFF, PHI_INT, PHI_OFF, COUNT };
enum class TmsPhiCols { PLCP, COUNT };
enum class TmsPhiInvCols { PLCPB, COUNT };

/** Splitting parameters for each move structure of a TmsIndex. */
struct TmsBuildOptions {
    // LF is unsplit by default, as in ms, so that the two tools take the
    // same LF steps and differ only in how they compute LCEs.
    orbit::split_params lf_split = orbit::NO_SPLITTING;
    // FL and phi are walked one dependent step at a time, so balancing them
    // keeps each step's fast-forward short.
    orbit::split_params fl_split = orbit::split_params{};
    orbit::split_params phi_split = orbit::split_params{};
    // Build phi_inv along with phi (same splitting as phi).
    bool phi_inv = true;
};

class TmsIndex {
public:
    using LF = orbit::rlbwt::lf_permutation<TmsLFCols, true, false>;
    using FL = orbit::rlbwt::fl_permutation<orbit::empty_data_columns, false, false>;
    using Phi = orbit::rlbwt::phi_permutation_impl<TmsPhiCols, true, true, orbit::move_vector>;
    using PhiInv = orbit::rlbwt::phi_inv_permutation_impl<TmsPhiInvCols, true, true, orbit::move_vector>;
    using LFPos = typename LF::position;
    using FLPos = typename FL::position;
    using PhiPos = typename Phi::position;
    using PhiInvPos = typename PhiInv::position;
    using position = LFPos;

    TmsIndex() = default;

    /**
     * Build from run heads (index codes: orbit::TERMINATOR, orbit::SEPARATOR,
     * or a nucleotide) and run lengths.  If run_tops is given, it holds the
     * LCP value at each run's head row (0 for row 0), and the index also gets
     * phi (and phi_inv if opts.phi_inv); this walks LF over all n rows.
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
            // the FL interval holding its tail row and the offset within it.
            ulint l_pos = 0, f_pos = 0, f_int = 0;
            const ulint f_count = fl_.intervals();
            for (ulint k = 0; k < lf_count; ++k) {
                l_pos += lf_enc.get_length(k);
                const ulint tail = l_pos - 1;
                while (f_int < f_count && f_pos + fl_.get_length(f_int) <= tail)
                    f_pos += fl_.get_length(f_int++);
                cols[k][col(TmsLFCols::PSI_INT)] = f_int;
                cols[k][col(TmsLFCols::PSI_OFF)] = tail - f_pos;
            }
        }
        if (run_tops) {
            if (run_tops->size() != heads.size()) throw std::invalid_argument("run_tops must have one value per run");
            build_phi(heads, lens, lf_enc, *run_tops, opts.phi_split, opts.phi_inv, cols);
            has_phi_ = true;
            has_phi_inv_ = opts.phi_inv;
        }
        lf_ = LF(lf_enc, cols);
        compute_occurs();
    }

    /** Write the index in Orbit's packed serialization, host byte order. */
    size_t serialize(std::ostream& out) {
        size_t bytes = 0;
        out.write(MAGIC, 4);
        const uint32_t v = VERSION, flags = (has_phi_ ? 1 : 0) | (has_phi_inv_ ? 2 : 0);
        out.write(reinterpret_cast<const char*>(&v), sizeof(v));
        out.write(reinterpret_cast<const char*>(&flags), sizeof(flags));
        bytes += 12;
        bytes += lf_.serialize(out);
        bytes += fl_.serialize(out);
        if (has_phi_) bytes += phi_.serialize(out);
        if (has_phi_inv_) bytes += phi_inv_.serialize(out);
        return bytes;
    }

    /**
     * Read an index written by serialize().  Throws on malformed input.
     * phi_inv, which only SMEM enumeration uses, is read only if
     * with_phi_inv; it is the last structure, so skipping it reads nothing
     * more.
     */
    void load(std::istream& in, bool with_phi_inv = true) {
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
        has_phi_inv_ = with_phi_inv && (flags & 2) != 0;
        if (has_phi_inv_) phi_inv_.load(in);
        if (!in.good()) throw std::runtime_error("truncated tms index");
        compute_occurs();
    }

    bool has_phi() const { return has_phi_; }
    bool has_phi_inv() const { return has_phi_inv_; }

    /** True if byte c occurs in the indexed text. */
    bool occurs(uchar c) const { return occurs_[c]; }

    // LF side, with the same names MSIndexSpillLCP uses.
    uchar get_character(ulint i) { return lf_.get_character(i); }
    ulint get_length(ulint i) const { return lf_.get_length(i); }
    LFPos LF_step(LFPos p) { return lf_.LF(p); }
    LFPos start_LF(LFPos p) const { return lf_.start_next(p); }
    LFPos finish_LF(LFPos p) const { return lf_.finish_next(p); }
    void prefetch(ulint i) const { lf_.prefetch(i); }
    /** Prefetch LF rows lo to hi inclusive, lo <= hi, once per cache line. */
    void prefetch_rows(ulint lo, ulint hi) const { lf_.prefetch_rows(lo, hi); }
    LFPos first() { return lf_.first(); }
    LFPos last() { return lf_.last(); }
    LFPos up(LFPos p) { return lf_.up(p); }
    LFPos down(LFPos p) { return lf_.down(p); }
    ulint move_runs() const { return lf_.intervals(); }
    ulint domain() const { return lf_.domain(); }

    /** FL point of the tail row of LF interval k. */
    FLPos psi_at_tail(ulint k) const {
        FLPos q;
        q.interval = lf_.template get<TmsLFCols::PSI_INT>(k);
        q.offset = lf_.template get<TmsLFCols::PSI_OFF>(k);
        return q;
    }
    /**
     * FL point of the head row of LF interval k, unresolved: its offset may
     * run past its interval's end, and finish_psi or resolve_psi resolves it.
     * Reads LF row k - 1 only.
     */
    FLPos psi_at_head_unresolved(ulint k) const {
        if (k == 0) return FLPos{};
        FLPos q = psi_at_tail(k - 1);
        ++q.offset;
        return q;
    }
    /** FL point of the head row of LF interval k. */
    FLPos psi_at_head(ulint k) const { return resolve_psi(psi_at_head_unresolved(k)); }

    // FL side.
    uchar psi_character(FLPos q) { return fl_.get_character(q.interval); }
    FLPos psi(FLPos q) { return fl_.FL(q); }
    FLPos start_psi(FLPos q) const { return fl_.start_next(q); }
    FLPos finish_psi(FLPos q) const { return fl_.finish_next(q); }
    /**
     * One psi step on an FL point q, resolved or not: resolve it, set c to
     * the first character of its suffix, and return the unresolved FL point
     * of the next character.  When FL's rows fit in a word, each row is read
     * with one load.
     */
    FLPos psi_step(FLPos q, uchar& c) {
        if (!fl_rows_fit_word_) {
            q = finish_psi(q);
            c = psi_character(q);
            return start_psi(q);
        }
        ulint w = fl_.row_bits(q.interval);
        ulint len = fl_.length_of(w);
        while (q.offset >= len) {
            q.offset -= len;
            w = fl_.row_bits(++q.interval);
            len = fl_.length_of(w);
        }
        c = fl_.character_of(w);
        FLPos next{};
        next.interval = fl_.pointer_of(w);
        next.offset = q.offset + fl_.offset_of(w);
        return next;
    }
    /** Resolve an FL point whose offset may run past its interval's end. */
    FLPos resolve_psi(FLPos q) const { return fl_.finish_next(q); }
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
    /** The first text position of phi interval i. */
    ulint phi_start(ulint i) const { return phi_.get_start(i); }
    /** The phi point of text position x, by binary search over interval starts. */
    PhiPos phi_at(ulint x) const { return locate<PhiPos>(phi_, x); }

    // phi_inv side, like the phi side.
    /** The phi_inv point of text position x, by binary search over interval starts. */
    PhiInvPos phi_inv_at(ulint x) const { return locate<PhiInvPos>(phi_inv_, x); }
    PhiInvPos phi_inv(PhiInvPos p) { return phi_inv_.phi_inv(p); }
    PhiInvPos start_phi_inv(PhiInvPos p) const { return phi_inv_.start_next(p); }
    PhiInvPos finish_phi_inv(PhiInvPos p) const { return phi_inv_.finish_next(p); }
    void prefetch_phi_inv(ulint i) const { phi_inv_.prefetch(i); }
    /** PLCPB at a phi_inv point: the LCP of its row with the row below. */
    ulint plcpb(PhiInvPos p) const { return phi_inv_.template get<TmsPhiInvCols::PLCPB>(p.interval) - p.offset; }
    ulint phi_inv_intervals() const { return phi_inv_.intervals(); }
    /** The first text position of phi_inv interval i. */
    ulint phi_inv_start(ulint i) const { return phi_inv_.get_start(i); }

    const LF& lf() const { return lf_; }
    const FL& fl() const { return fl_; }

    /** The alphabet code of byte c in LF's and FL's character columns, or no_code if c does not occur. */
    uchar lf_code(uchar c) { return occurs_[c] ? lf_.character_code(c).value_or(no_code) : no_code; }
    uchar fl_code(uchar c) { return occurs_[c] ? fl_.character_code(c).value_or(no_code) : no_code; }
    static constexpr uchar no_code = 0xFF;

    /**
     * Row access for the batched engine through local copies of LF's and
     * FL's packed layouts (see Orbit's packed_matrix::reader).  Kept in a
     * local variable, it lets the compiler hold the layouts in registers,
     * where reads through the index reload them after any store that might
     * alias them.  An LF row is named by its first bit, row(i).  Its pointer
     * and offset columns are read with one load, as are its character and
     * PSI columns, and an FL row is read whole with one load; fits() says
     * whether the index's column widths allow that.  Characters are alphabet
     * codes (lf_code, fl_code).
     */
    struct PackedAccess {
        using LFRows = decltype(std::declval<const LF&>().get_reader());
        using FLRows = decltype(std::declval<const FL&>().get_reader());
        static constexpr size_t len_col = LF::length_column();
        static constexpr size_t ptr_col = LF::pointer_column();
        static constexpr size_t off_col = LF::offset_column();
        static constexpr size_t chr_col = LF::character_column();
        static constexpr size_t psi_int_col = LF::template data_column<TmsLFCols::PSI_INT>();
        static constexpr size_t psi_off_col = LF::template data_column<TmsLFCols::PSI_OFF>();
        static constexpr size_t fl_len_col = FL::length_column();
        static constexpr size_t fl_ptr_col = FL::pointer_column();
        static constexpr size_t fl_off_col = FL::offset_column();
        static constexpr size_t fl_chr_col = FL::character_column();
        static_assert(ptr_col < off_col && chr_col < psi_int_col && psi_int_col < psi_off_col, "unexpected LF column order");
        static_assert(fl_len_col == 0 && fl_ptr_col <= 3 && fl_off_col <= 3 && fl_chr_col <= 3, "unexpected FL columns");
        LFRows lf;
        FLRows fl;
        using Row = ulint;

        bool fits() const {
            return lf.template span_fits<ptr_col, off_col>() && lf.template span_fits<chr_col, psi_off_col>() &&
                   fl.template span_fits<0, 3>();
        }
        Row row(ulint i) const { return lf.row_start(i); }
        ulint length(Row r) const { return lf.template get_at<len_col>(r); }
        uchar code(Row r) const { return static_cast<uchar>(lf.template get_at<chr_col>(r)); }
        /** Resolve an LF position whose offset may run past its interval's end; returns its row. */
        Row resolve_LF(LFPos& p) const {
            Row r = row(p.interval);
            for (ulint len = length(r); p.offset >= len; len = length(r)) {
                p.offset -= len;
                r = row(++p.interval);
            }
            return r;
        }
        /** start_LF of position p, resolved, whose row is r. */
        LFPos start_LF(Row r, LFPos p) const {
            const ulint span = lf.template get_span<ptr_col>(r);
            p.interval = lf.template extract_span<ptr_col, ptr_col>(span);
            p.offset += lf.template extract_span<ptr_col, off_col>(span);
            return p;
        }
        /** psi_at_tail of the LF interval whose row is r. */
        FLPos psi_at_tail(Row r) const {
            const ulint span = lf.template get_span<chr_col>(r);
            FLPos q;
            q.interval = lf.template extract_span<chr_col, psi_int_col>(span);
            q.offset = lf.template extract_span<chr_col, psi_off_col>(span);
            return q;
        }
        /** As TmsIndex::psi_step, with c an alphabet code. */
        FLPos psi_step(FLPos q, uchar& c) const {
            ulint start = fl.row_start(q.interval);
            ulint w = fl.template get_span<0>(start);
            ulint len = fl.template extract_span<0, fl_len_col>(w);
            while (q.offset >= len) {
                q.offset -= len;
                ++q.interval;
                start += fl.row_width;
                w = fl.template get_span<0>(start);
                len = fl.template extract_span<0, fl_len_col>(w);
            }
            c = static_cast<uchar>(fl.template extract_span<0, fl_chr_col>(w));
            FLPos next;
            next.interval = fl.template extract_span<0, fl_ptr_col>(w);
            next.offset = q.offset + fl.template extract_span<0, fl_off_col>(w);
            return next;
        }
        void prefetch(ulint i) const { lf.prefetch(i); }
        void prefetch_rows(ulint lo, ulint hi) const { lf.prefetch_rows(lo, hi); }
        void prefetch_psi(ulint i) const { fl.prefetch(i); }
    };
    PackedAccess packed_access() const { return PackedAccess{lf_.get_reader(), fl_.get_reader()}; }

    /** The same row access through the index's own column reads, with bytes for characters. */
    struct ColumnAccess {
        TmsIndex* idx;  // Orbit's character reads are not const
        using Row = ulint;
        Row row(ulint i) const { return i; }
        ulint length(Row i) const { return idx->get_length(i); }
        uchar code(Row i) const { return idx->get_character(i); }
        Row resolve_LF(LFPos& p) const {
            p = idx->finish_LF(p);
            return p.interval;
        }
        LFPos start_LF(Row, LFPos p) const { return idx->start_LF(p); }
        FLPos psi_at_tail(Row i) const { return idx->psi_at_tail(i); }
        FLPos psi_step(FLPos q, uchar& c) const { return idx->psi_step(q, c); }
        void prefetch(ulint i) const { idx->prefetch(i); }
        void prefetch_rows(ulint lo, ulint hi) const { idx->prefetch_rows(lo, hi); }
        void prefetch_psi(ulint i) const { idx->prefetch_psi(i); }
    };
    ColumnAccess column_access() { return ColumnAccess{this}; }

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
        if (has_phi_inv_) {
            const auto& pw = phi_inv_.get_widths();
            os << "phi_inv: intervals=" << phi_inv_.intervals() << " widths(start,ptr,off,plcpb)=";
            widths(pw);
            os << " row_bits=" << row_bits(pw) << "\n";
            total += double(row_bits(pw)) * phi_inv_.intervals();
        }
        os << "total: " << total / 8 / 1e6 << " MB, " << total / lf_.runs() << " bits/run\n";
    }

private:
    static constexpr char MAGIC[4] = {'T', 'M', 'S', 'X'};
    static constexpr uint32_t VERSION = 4;

    LF lf_;
    FL fl_;
    Phi phi_;
    PhiInv phi_inv_;
    bool has_phi_ = false, has_phi_inv_ = false;
    std::array<bool, 256> occurs_{};
    bool fl_rows_fit_word_ = false;

    static constexpr size_t col(TmsLFCols c) { return static_cast<size_t>(c); }

    // Derived fields: whole-row FL reads, and which characters occur.
    void compute_occurs() {
        fl_rows_fit_word_ = fl_.row_fits_word();
        occurs_.fill(false);
        for (ulint i = 0; i < lf_.intervals(); ++i) occurs_[lf_.get_character(i)] = true;
    }

    /** The point of text position x in a starts-based permutation. */
    template <typename Pos, typename Perm>
    static Pos locate(const Perm& perm, ulint x) {
        ulint lo = 0, hi = perm.intervals();  // perm.get_start(lo) <= x < perm.get_start(hi)
        while (hi - lo > 1) {
            const ulint mid = lo + (hi - lo) / 2;
            if (perm.get_start(mid) <= x) lo = mid;
            else hi = mid;
        }
        Pos p;
        p.interval = lo;
        p.offset = x - perm.get_start(lo);
        p.idx = x;
        return p;
    }

    /**
     * A starts-based permutation over text positions whose intervals start at
     * starts (sorted, starting at 0) with the given images, and one integrated
     * column holding an LCP value that drops by one per position inside an
     * interval, sampled at each start.  Split intervals get the sample minus
     * their distance from the original start.  Also returns the starts of the
     * split intervals.
     */
    template <typename Perm, typename Cols>
    static std::pair<Perm, std::vector<ulint>> sampled_permutation(const std::vector<ulint>& starts,
                                                                   const std::vector<ulint>& images,
                                                                   const std::vector<ulint>& samples, ulint n,
                                                                   const orbit::split_params& sp) {
        if (starts.empty() || starts[0] != 0) throw std::logic_error("a permutation over text positions must start at 0");
        std::vector<ulint> lengths(starts.size());
        ulint max_length = 0;
        for (size_t x = 0; x < starts.size(); ++x) {
            lengths[x] = (x + 1 < starts.size() ? starts[x + 1] : n) - starts[x];
            max_length = std::max(max_length, lengths[x]);
        }
        auto enc = orbit::interval_encoding_impl<>::from_lengths_and_images(lengths, images, n, max_length, sp);
        std::vector<orbit::columns_tuple<Cols>> cols(enc.intervals());
        std::vector<ulint> split_starts(enc.intervals());
        ulint s = 0;
        size_t x = 0;
        for (ulint i = 0; i < enc.intervals(); ++i) {
            while (x + 1 < starts.size() && starts[x + 1] <= s) ++x;
            const ulint into = s - starts[x];
            if (samples[x] < into) throw std::runtime_error("LCP samples are inconsistent with the RLBWT");
            cols[i][0] = samples[x] - into;
            split_starts[i] = s;
            s += enc.get_length(i);
        }
        return {Perm(enc, cols), std::move(split_starts)};
    }

    /**
     * Build phi (and phi_inv if with_inv) and fill the PHI columns.  Walking
     * LF from row 0, whose text position is n - 1, gives the text position of
     * every row; the walk records it at each LF interval's head and tail.
     *
     * phi intervals start at the text positions of true run heads (rows whose
     * BWT character differs from the row above's, and row 0), and the
     * interval starting at SA[h] maps to SA[h - 1] (SA[0] to SA[n - 1]).  The
     * PLCP sample at SA[h] is LCP[h], the top of the run whose head is h.
     *
     * phi_inv intervals start at the text positions of true run tails (rows
     * whose BWT character differs from the row below's, and row n - 1), and
     * the interval starting at SA[t] maps to SA[t + 1] (SA[n - 1] to SA[0] =
     * n - 1).  The PLCPB sample at SA[t] is LCP[t + 1], the top of the next
     * run, or 0 for row n - 1.
     */
    template <typename Enc>
    void build_phi(const std::vector<uchar>& heads, const std::vector<ulint>& lens, const Enc& lf_enc,
                   const std::vector<ulint>& run_tops, const orbit::split_params& sp, bool with_inv,
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
        const ulint n = lf_enc.domain();
        const auto& lf_heads = lf_enc.get_heads();
        std::vector<ulint> order, starts, images, samples;
        auto fill = [&](auto start_of, auto image_of, auto sample_of) {
            std::sort(order.begin(), order.end(), [&](ulint a, ulint b) { return start_of(a) < start_of(b); });
            starts.resize(order.size());
            images.resize(order.size());
            samples.resize(order.size());
            for (size_t x = 0; x < order.size(); ++x) {
                starts[x] = start_of(order[x]);
                images[x] = image_of(order[x]);
                samples[x] = sample_of(order[x]);
            }
        };
        // phi from true run heads.
        order.clear();
        for (ulint k = 0; k < lf_count; ++k)
            if (k == 0 || lf_heads[k] != lf_heads[k - 1]) order.push_back(k);
        fill([&](ulint k) { return sa_head[k]; },
             [&](ulint k) { return sa_tail[k == 0 ? lf_count - 1 : k - 1]; },
             [&](ulint k) {
                 if (top[k] == NOT_A_RUN_HEAD) throw std::logic_error("a true run head must start an original run");
                 return top[k];
             });
        std::vector<ulint> split_starts;
        std::tie(phi_, split_starts) = sampled_permutation<Phi, TmsPhiCols>(starts, images, samples, n, sp);
        for (ulint k = 0; k < lf_count; ++k) {
            const ulint p = sa_head[k];
            const ulint i = static_cast<ulint>(std::upper_bound(split_starts.begin(), split_starts.end(), p) -
                                               split_starts.begin()) - 1;
            cols[k][col(TmsLFCols::PHI_INT)] = i;
            cols[k][col(TmsLFCols::PHI_OFF)] = p - split_starts[i];
        }
        if (!with_inv) return;
        // phi_inv from true run tails.
        order.clear();
        for (ulint k = 0; k < lf_count; ++k)
            if (k + 1 == lf_count || lf_heads[k] != lf_heads[k + 1]) order.push_back(k);
        fill([&](ulint k) { return sa_tail[k]; },
             [&](ulint k) { return sa_head[k + 1 == lf_count ? 0 : k + 1]; },
             [&](ulint k) {
                 if (k + 1 == lf_count) return ulint(0);
                 if (top[k + 1] == NOT_A_RUN_HEAD) throw std::logic_error("a true run head must start an original run");
                 return top[k + 1];
             });
        std::tie(phi_inv_, split_starts) = sampled_permutation<PhiInv, TmsPhiInvCols>(starts, images, samples, n, sp);
    }

    static constexpr ulint NOT_A_RUN_HEAD = std::numeric_limits<ulint>::max();
};

#endif /* _TMS_INDEX_HPP */
