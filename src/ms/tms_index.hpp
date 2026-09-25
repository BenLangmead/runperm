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
 *  - LF over the BWT, runs-based (rows store lengths), split and balanced
 *    with Orbit's splitting.  In the psi and full layouts its rows are the
 *    union of the split LF intervals and their images, so that each row lies
 *    inside one LF interval and inside one image, and psi (Orbit's FL) is a
 *    move structure over the same rows.  LF has integrated columns per row:
 *    PSI_INT and PSI_OFF, psi's pointer and offset for the row's head, and
 *    PHI_INT and PHI_OFF, the phi point of its head row's text position.
 *    The phi point of a tail row is one phi step from the next row's head's
 *    phi point, so a reposition reads both candidates' start points from
 *    rows its scan has already read.
 *  - With psi, a table of F block starts.  An F block starts at a row start
 *    (the image of the first run of its character), so a row's F character
 *    is the block its index falls in.
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
 * LF positions and psi points are the same (interval, offset) pairs, so a
 * candidate row's psi walk starts at the row itself: its tail is
 * (u, length - 1) and its head (d, 0).
 *
 * LF, its rows and the PSI columns come from the run heads and lengths alone
 * in O(r) time and space.  phi also needs the LCP value at each run head and
 * a walk over all n rows.  LF and phi each take their own splitting
 * parameters.
 *
 * An index is built with one of three layouts (TmsLayout), so that it holds
 * only what its queries read:
 *
 *  - full: everything above.  Built without LCP values, it holds what psi
 *    holds.
 *  - psi: LF with the PSI columns and the F block starts.  Enough for
 *    matching statistics with psi, without positions or SMEMs.
 *  - phi: LF with the PHI columns, phi and optionally phi_inv.  LF's rows are
 *    then the split LF intervals alone, since nothing walks psi.  Enough for
 *    every query that does not walk psi.
 *
 * The LF data columns a layout omits have width 0 (see Orbit's
 * permutation_impl), so they take no bits and read as 0, and the columns
 * kept stay where the full layout has them.  The header records which parts
 * the index holds (TmsParts), and load() reads only the parts a query needs,
 * seeking past the rest.
 */

#ifndef _TMS_INDEX_HPP
#define _TMS_INDEX_HPP

#include "orbit/rlbwt.hpp"
#include "orbit/common.hpp"
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using uchar = orbit::uchar;
using ulint = orbit::ulint;

// Inline a function the batched engine calls in its hot loop, where GCC
// otherwise may keep it out of line.
#if defined(__GNUC__)
#define TMS_ALWAYS_INLINE inline __attribute__((always_inline))
#else
#define TMS_ALWAYS_INLINE inline
#endif

enum class TmsLFCols { PSI_INT, PSI_OFF, PHI_INT, PHI_OFF, COUNT };
enum class TmsPhiCols { PLCP, COUNT };
enum class TmsPhiInvCols { PLCPB, COUNT };

/** Which structures a TmsIndex is built with (see the file comment). */
enum class TmsLayout { FULL, PSI, PHI };

/**
 * Parts of a TmsIndex beyond LF, which every query reads: psi is LF's PSI
 * columns and the F block starts, phi is phi with LF's PHI columns, and
 * phi_inv is phi_inv.
 * An index holds some of them, and a query needs some of them.
 */
struct TmsParts {
    bool psi = false, phi = false, phi_inv = false;
};

/** Splitting parameters for each move structure of a TmsIndex. */
struct TmsBuildOptions {
    // LF and phi are walked one dependent step at a time, so balancing them
    // keeps each step's fast-forward short.  LF must be split: its balancing
    // also bounds row lengths and offsets, and in the psi and full layouts
    // the union with the images makes the rows psi's too.
    orbit::split_params lf_split = orbit::split_params{};
    orbit::split_params phi_split = orbit::split_params{};
    // Build phi_inv along with phi (same splitting as phi).
    bool phi_inv = true;
    TmsLayout layout = TmsLayout::FULL;
};

class TmsIndex {
public:
    using LF = orbit::rlbwt::lf_permutation<TmsLFCols, true, false>;
    using Phi = orbit::rlbwt::phi_permutation_impl<TmsPhiCols, true, true, orbit::move_vector>;
    using PhiInv = orbit::rlbwt::phi_inv_permutation_impl<TmsPhiInvCols, true, true, orbit::move_vector>;
    using LFPos = typename LF::position;
    // A psi point: an LF position whose F character and psi step psi reads.
    using FLPos = LFPos;
    using PhiPos = typename Phi::position;
    using PhiInvPos = typename PhiInv::position;
    using position = LFPos;
    // The most F blocks after the first: one per character of the
    // alphabet but the smallest.
    static constexpr size_t F_BLOCKS = orbit::nucleotide::SIGMA - 1;

    TmsIndex() = default;

    /**
     * Build from run heads (index codes: orbit::TERMINATOR, orbit::SEPARATOR,
     * or a nucleotide) and run lengths.  If run_tops is given, it holds the
     * LCP value at each run's head row (0 for row 0), and the index also gets
     * phi (and phi_inv if opts.phi_inv); this walks LF over all n rows.
     * opts.layout says which parts to build: the psi layout takes no
     * run_tops, and the phi layout needs them.  opts.lf_split must split.
     */
    TmsIndex(const std::vector<uchar>& heads, const std::vector<ulint>& lens,
             const TmsBuildOptions& opts = TmsBuildOptions{}, const std::vector<ulint>* run_tops = nullptr) {
        using Enc = orbit::rlbwt::rlbwt_interval_encoding<>;
        if (opts.layout == TmsLayout::PSI && run_tops) throw std::invalid_argument("a psi layout takes no LCP values");
        if (opts.layout == TmsLayout::PHI && !run_tops) throw std::invalid_argument("a phi layout needs LCP values");
        check_lf_split(opts.lf_split);
        has_psi_ = opts.layout != TmsLayout::PHI;
        Enc lf_enc;
        std::vector<typename LF::data_tuple> cols;
        if (has_psi_) {
            std::vector<uchar> row_heads;
            std::vector<ulint> row_lens;
            union_rows(heads, lens, opts.lf_split, row_heads, row_lens);
            cols.resize(row_lens.size());
            fill_psi(row_heads, row_lens, cols);
            lf_enc = Enc::lf_interval_encoding(row_heads, row_lens, orbit::NO_SPLITTING);
        } else {
            lf_enc = Enc::lf_interval_encoding(heads, lens, opts.lf_split);
            cols.resize(lf_enc.intervals());
        }
        if (run_tops) {
            if (run_tops->size() != heads.size()) throw std::invalid_argument("run_tops must have one value per run");
            build_phi(lens, lf_enc, *run_tops, opts.phi_split, opts.phi_inv, cols);
            has_phi_ = true;
            has_phi_inv_ = opts.phi_inv;
        }
        // The columns of parts the index lacks hold zeros and take no bits.
        lf_ = LF(lf_enc, cols, {!has_psi_, !has_psi_, !has_phi_, !has_phi_});
        compute_occurs();
        if (has_phi_inv_) build_start_tables();
    }

    /**
     * Throw unless sp splits LF: TmsIndex's LF is always split, with length
     * capping or balancing (see TmsBuildOptions).
     */
    static void check_lf_split(const orbit::split_params& sp) {
        const bool caps = sp.length_capping.has_value() && *sp.length_capping > 0;
        const bool balances = sp.balancing.has_value() && *sp.balancing > 0;
        if (!caps && !balances)
            throw std::invalid_argument("LF must be split, with length capping or balancing; unsplit LF is not supported");
    }

    /**
     * Write the index in Orbit's packed serialization, host byte order: a
     * header with the parts the index holds, LF, with psi the table of F
     * block starts, and then each part's structure present, phi and phi_inv,
     * after its size in bytes so that load() can skip it.
     */
    size_t serialize(std::ostream& out) {
        size_t bytes = 0;
        out.write(MAGIC, 4);
        const uint32_t v = VERSION, flags = parts_flags(parts());
        out.write(reinterpret_cast<const char*>(&v), sizeof(v));
        out.write(reinterpret_cast<const char*>(&flags), sizeof(flags));
        bytes += 12;
        bytes += lf_.serialize(out);
        if (has_psi_) {
            for (ulint s : f_rows_) {
                const uint64_t x = s;
                out.write(reinterpret_cast<const char*>(&x), sizeof(x));
            }
            bytes += F_BLOCKS * sizeof(uint64_t);
        }
        if (has_phi_) bytes += serialize_sized(phi_, out);
        if (has_phi_inv_) bytes += serialize_sized(phi_inv_, out);
        return bytes;
    }

    /** Read an index written by serialize(), with every part it holds.  Throws on malformed input. */
    void load(std::istream& in) { load(in, nullptr); }

    /**
     * Read an index written by serialize(), with only LF and the parts in
     * need, seeking past the rest.  Throws, before reading LF, if the index
     * lacks a part in need.
     */
    void load(std::istream& in, const TmsParts& need) { load(in, &need); }

    bool has_psi() const { return has_psi_; }
    bool has_phi() const { return has_phi_; }
    bool has_phi_inv() const { return has_phi_inv_; }
    /** The parts the index holds, or those loaded if load() skipped some. */
    TmsParts parts() const { return TmsParts{has_psi_, has_phi_, has_phi_inv_}; }
    /** The name of the layout with the given parts: full, psi, phi, or none (LF alone). */
    static std::string layout_name(const TmsParts& p) {
        return p.psi && p.phi ? "full" : p.psi ? "psi" : p.phi ? "phi" : "none";
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
    /** Prefetch LF rows lo to hi inclusive, lo <= hi, once per cache line. */
    void prefetch_rows(ulint lo, ulint hi) const { lf_.prefetch_rows(lo, hi); }
    LFPos first() { return lf_.first(); }
    LFPos last() { return lf_.last(); }
    LFPos up(LFPos p) { return lf_.up(p); }
    LFPos down(LFPos p) { return lf_.down(p); }
    ulint move_runs() const { return lf_.intervals(); }
    ulint domain() const { return lf_.domain(); }

    // psi side.  psi points are LF positions, and psi reads LF's rows.
    /** psi point of the tail row of LF interval k: (k, length - 1). */
    FLPos psi_at_tail(ulint k) const {
        FLPos q;
        q.interval = k;
        q.offset = lf_.get_length(k) - 1;
        return q;
    }
    /** psi point of the head row of LF interval k: (k, 0). */
    FLPos psi_at_head(ulint k) const {
        FLPos q;
        q.interval = k;
        q.offset = 0;
        return q;
    }
    /** The F character of a resolved psi point: the first character of its row's suffix. */
    uchar psi_character(FLPos q) const { return f_byte_[f_block(q.interval)]; }
    FLPos psi(FLPos q) const { return finish_psi(start_psi(q)); }
    FLPos start_psi(FLPos q) const {
        FLPos next;
        next.interval = lf_.template get<TmsLFCols::PSI_INT>(q.interval);
        next.offset = lf_.template get<TmsLFCols::PSI_OFF>(q.interval) + q.offset;
        return next;
    }
    /** Resolve a psi point whose offset may run past its row's end. */
    FLPos finish_psi(FLPos q) const { return lf_.finish_next(q); }
    FLPos resolve_psi(FLPos q) const { return finish_psi(q); }
    /**
     * One psi step on a psi point q, resolved or not: resolve it, set c to
     * the first character of its suffix, and return the unresolved psi point
     * of the next character.
     */
    FLPos psi_step(FLPos q, uchar& c) const {
        q = finish_psi(q);
        c = psi_character(q);
        return start_psi(q);
    }
    void prefetch_psi(ulint i) const { lf_.prefetch(i); }
    ulint psi_intervals() const { return lf_.intervals(); }

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
    /** The phi point of text position x, by binary search over interval starts narrowed by the start table. */
    PhiPos phi_at(ulint x) const { return phi_walker().locate(x); }

    // phi_inv side, like the phi side.
    /** The phi_inv point of text position x, as phi_at. */
    PhiInvPos phi_inv_at(ulint x) const { return phi_inv_walker().locate(x); }
    PhiInvPos phi_inv(PhiInvPos p) { return phi_inv_.phi_inv(p); }
    PhiInvPos start_phi_inv(PhiInvPos p) const { return phi_inv_.start_next(p); }
    PhiInvPos finish_phi_inv(PhiInvPos p) const { return phi_inv_.finish_next(p); }
    void prefetch_phi_inv(ulint i) const { phi_inv_.prefetch(i); }
    /** PLCPB at a phi_inv point: the LCP of its row with the row below. */
    ulint plcpb(PhiInvPos p) const { return phi_inv_.template get<TmsPhiInvCols::PLCPB>(p.interval) - p.offset; }
    ulint phi_inv_intervals() const { return phi_inv_.intervals(); }
    /** The first text position of phi_inv interval i. */
    ulint phi_inv_start(ulint i) const { return phi_inv_.get_start(i); }

    /**
     * A table over the high bits of text positions that narrows the search
     * for a text position's phi or phi_inv interval: entry b holds the
     * interval of text position b << shift, so the interval of x lies
     * between entries x >> shift and (x >> shift) + 1.  The last entry, past
     * every bucket, holds the last interval.  Built from the interval starts,
     * not stored in the index.
     */
    class StartTable {
    public:
        using Rows = orbit::packed_matrix<1>;

        StartTable() = default;
        /** The table of a starts-based permutation with 2^shift positions per entry. */
        template <typename Perm>
        StartTable(const Perm& perm, unsigned shift) : shift_(shift), built_(true) {
            const ulint n = perm.domain(), last = perm.intervals() - 1;
            const ulint buckets = ((n - 1) >> shift) + 1;
            rows_ = Rows(buckets + 1, {std::max<uchar>(1, orbit::bit_width(last))});
            // Interval i gets the buckets whose first position lies in it.
            const auto rd = perm.get_reader();
            ulint b = 0;
            for (ulint i = 0; i <= last; ++i) {
                const ulint end = i < last ? rd.template get<Perm::start_column()>(i + 1) : n;
                for (const ulint to = ((end - 1) >> shift) + 1; b < to; ++b) rows_.template set<0>(b, i);
            }
            rows_.template set<0>(buckets, last);
        }
        bool empty() const { return !built_; }
        unsigned shift() const { return shift_; }
        /** Bytes the table takes. */
        size_t bytes() const { return built_ ? rows_.data_size() : 0; }
        const Rows& rows() const { return rows_; }

    private:
        unsigned shift_ = 0;
        bool built_ = false;
        Rows rows_{};  // zero-initialized, so that an unbuilt table's reader is well defined
    };

    /**
     * The table shift for a permutation over n text positions with the given
     * intervals: the smallest whose entries cover, on average, at least
     * START_TABLE_DENSITY intervals each.
     */
    static unsigned start_table_shift(ulint n, ulint intervals) {
        unsigned s = 0;
        while ((ulint(1) << s) * intervals < START_TABLE_DENSITY * n && s < 62) ++s;
        return s;
    }
    static constexpr ulint START_TABLE_DENSITY = 4;

    /**
     * Build the start tables of phi and phi_inv, which phi_at, phi_inv_at and
     * the walkers' locate use, with 2^shift text positions per entry;
     * shift = -1 picks start_table_shift.  load() builds them when it reads
     * phi_inv, and so does the constructor that builds phi_inv.
     */
    void build_start_tables(int shift = -1) {
        auto pick = [&](const auto& perm) {
            return shift < 0 ? start_table_shift(perm.domain(), perm.intervals()) : unsigned(shift);
        };
        if (has_phi_) phi_table_ = StartTable(phi_, pick(phi_));
        if (has_phi_inv_) phi_inv_table_ = StartTable(phi_inv_, pick(phi_inv_));
    }
    /** Bytes the start tables take. */
    size_t start_table_bytes() const { return phi_table_.bytes() + phi_inv_table_.bytes(); }

    /**
     * Steps of phi or phi_inv walks, and locating text positions in them,
     * from local copies of what reading the rows and the start table needs
     * (see packed_matrix::reader).  A step reads its row's start, the next
     * row's start, and the pointer, offset and LCP columns, the last three
     * with one load when they fit in one.  It fast-forwards one interval at
     * a time.  Valid while the index is alive and unchanged.
     */
    template <typename Perm, auto LcpCol>
    class TextWalker {
        using Reader = decltype(std::declval<const Perm&>().get_reader());
        using TableReader = decltype(std::declval<const StartTable::Rows&>().get_reader());
        static constexpr size_t START = Perm::start_column(), PTR = Perm::pointer_column(),
                                OFF = Perm::offset_column(), LCP = Perm::template data_column<LcpCol>();
        static_assert(PTR + 1 == OFF && OFF + 1 == LCP, "a step reads the pointer, offset and LCP columns together");

    public:
        using Pos = typename Perm::position;

        TextWalker(const Perm& perm, const StartTable& table)
            : rd_(perm.get_reader()), tr_(table.rows().get_reader()), last_(perm.intervals() - 1),
              n_(perm.domain()), shift_(table.shift()), has_table_(!table.empty()),
              span_(rd_.template span_fits<PTR, LCP>()) {}

        /** Hint that interval i's row, and the rows in the next cache line, will be read soon. */
        void prefetch(ulint i) const {
            const auto* a = rd_.data + rd_.row_start(i) / 8;
            ORBIT_PREFETCH(a);
            ORBIT_PREFETCH(a + 64);
        }
        /** Hint that the start table entries for text position x will be read soon. */
        void prefetch_table(ulint x) const {
            if (has_table_) tr_.prefetch(x >> shift_);
        }
        /** Intervals lo and hi such that x's interval lies between them. */
        void range(ulint x, ulint& lo, ulint& hi) const {
            if (!has_table_) {
                lo = 0;
                hi = last_;
                return;
            }
            const ulint b = x >> shift_;
            lo = tr_.template get<0>(b);
            hi = tr_.template get<0>(b + 1);
        }
        /** Hint that the rows of intervals lo to hi will be read soon, at most max_lines cache lines. */
        void prefetch_range(ulint lo, ulint hi, size_t max_lines) const {
            const size_t a = rd_.row_start(lo) / 8, b = rd_.row_start(hi) / 8;
            const auto* p = rd_.data + a;
            const auto* end = rd_.data + std::min<size_t>(b + 16, a + 64 * max_lines);
            for (; p < end; p += 64) ORBIT_PREFETCH(p);
        }
        /** The resolved point of text position x, whose interval lies between lo and hi. */
        Pos locate(ulint x, ulint lo, ulint hi) const {
            ulint base = lo, len = hi - lo + 1;
            while (len > 1) {
                const ulint half = len / 2;
                base = start(base + half) <= x ? base + half : base;
                len -= half;
            }
            Pos p;
            p.interval = base;
            p.offset = x - start(base);
            p.idx = x;
            return p;
        }
        /** The resolved point of text position x. */
        Pos locate(ulint x) const {
            ulint lo, hi;
            range(x, lo, hi);
            return locate(x, lo, hi);
        }
        /** The LCP column at a resolved point, minus its offset: PLCP or PLCPB there. */
        ulint lcp(Pos p) const { return rd_.template get<LCP>(p.interval) - p.offset; }
        /** The unresolved point one step on from a resolved point. */
        Pos start_step(Pos p) const {
            const size_t row = rd_.row_start(p.interval);
            Pos q;
            q.interval = rd_.template get_at<PTR>(row);
            q.offset = rd_.template get_at<OFF>(row) + p.offset;
            q.idx = 0;
            return q;
        }
        /**
         * One step of a walk.  p is unresolved, from start_step or a
         * previous step.  Resolves it, sets lcp to the LCP there (PLCP or
         * PLCPB) and returns its text position; p becomes the unresolved
         * point one step on.
         */
        ulint step(Pos& p, ulint& lcp) const {
            ulint i = p.interval;
            size_t row = rd_.row_start(i);
            ulint first = rd_.template get_at<START>(row);
            const ulint x = first + p.offset;
            // Row last + 1 lies in the matrix's padding, so reading its
            // start is safe; the value is replaced by n.
            ulint next = rd_.template get_at<START>(row + rd_.row_width);
            if (i == last_) next = n_;
            while (x >= next) {
                ++i;
                row += rd_.row_width;
                first = next;
                next = rd_.template get_at<START>(row + rd_.row_width);
                if (i == last_) next = n_;
            }
            const ulint off = x - first;
            ulint ptr, poff, l;
            if (span_) {
                const ulint s = rd_.template get_span<PTR>(row);
                ptr = rd_.template extract_span<PTR, PTR>(s);
                poff = rd_.template extract_span<PTR, OFF>(s);
                l = rd_.template extract_span<PTR, LCP>(s);
            } else {
                ptr = rd_.template get_at<PTR>(row);
                poff = rd_.template get_at<OFF>(row);
                l = rd_.template get_at<LCP>(row);
            }
            lcp = l - off;
            p.interval = ptr;
            p.offset = poff + off;
            return x;
        }

    private:
        Reader rd_;
        TableReader tr_;
        ulint last_, n_;
        unsigned shift_;
        bool has_table_, span_;

        ulint start(ulint i) const { return rd_.template get<START>(i); }
    };
    using PhiWalker = TextWalker<Phi, TmsPhiCols::PLCP>;
    using PhiInvWalker = TextWalker<PhiInv, TmsPhiInvCols::PLCPB>;
    PhiWalker phi_walker() const { return PhiWalker(phi_, phi_table_); }
    PhiInvWalker phi_inv_walker() const { return PhiInvWalker(phi_inv_, phi_inv_table_); }

    const LF& lf() const { return lf_; }
    /**
     * The alphabet code of byte c in LF's character column, or no_code if c
     * does not occur; and the code of c as an F character, its F block, or
     * no_code if c does not occur or the index lacks psi.
     */
    uchar lf_code(uchar c) { return occurs_[c] ? lf_.character_code(c).value_or(no_code) : no_code; }
    uchar fl_code(uchar c) const { return occurs_[c] && has_psi_ ? f_block_of_[c] : no_code; }
    static constexpr uchar no_code = 0xFF;

    /**
     * Row access for the batched engine through a local copy of LF's packed
     * layout (see Orbit's packed_matrix::reader) and of the F block starts.
     * Kept in a local variable, it lets the compiler hold them in registers,
     * where reads through the index reload them after any store that might
     * alias them.  A row is named by its first bit, row(i).  Its pointer and
     * offset columns are read with one load, as are its character and PSI
     * columns; fits() says whether the index's column widths allow that.  A
     * psi step reads a row's length and then its PSI columns, and takes the
     * F character from the block starts.  Characters are alphabet codes:
     * lf_code for LF's, fl_code for F's.  Without psi, only the LF reads are
     * valid.
     */
    struct PackedAccess {
        using LFRows = decltype(std::declval<const LF&>().get_reader());
        static constexpr size_t len_col = LF::length_column();
        static constexpr size_t ptr_col = LF::pointer_column();
        static constexpr size_t off_col = LF::offset_column();
        static constexpr size_t chr_col = LF::character_column();
        static constexpr size_t psi_int_col = LF::template data_column<TmsLFCols::PSI_INT>();
        static constexpr size_t psi_off_col = LF::template data_column<TmsLFCols::PSI_OFF>();
        static_assert(ptr_col < off_col && chr_col < psi_int_col && psi_int_col < psi_off_col, "unexpected LF column order");
        LFRows lf;
        // The first row of F blocks 1 to F_BLOCKS, or a value past every row.
        std::array<ulint, F_BLOCKS> f_rows;
        using Row = ulint;

        bool fits() const {
            return lf.template span_fits<ptr_col, off_col>() && lf.template span_fits<chr_col, psi_off_col>();
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
        /** The F block, and so the F character's code, of row i. */
        uchar f_code(ulint i) const {
            uchar b = 0;
            for (size_t k = 0; k < F_BLOCKS; ++k) b += i >= f_rows[k];
            return b;
        }
        /** psi_at_tail of LF interval k: (k, length - 1). */
        FLPos psi_at_tail(ulint k) const {
            FLPos q;
            q.interval = k;
            q.offset = length(row(k)) - 1;
            return q;
        }
        /** psi_at_head of LF interval k: (k, 0). */
        FLPos psi_at_head(ulint k) const {
            FLPos q;
            q.interval = k;
            q.offset = 0;
            return q;
        }
        /** As TmsIndex::psi_step, with c the F character's code. */
        TMS_ALWAYS_INLINE FLPos psi_step(FLPos q, uchar& c) const {
            Row r = row(q.interval);
            for (ulint len = length(r); q.offset >= len; len = length(r)) {
                q.offset -= len;
                r = row(++q.interval);
            }
            c = f_code(q.interval);
            const ulint span = lf.template get_span<chr_col>(r);
            FLPos next;
            next.interval = lf.template extract_span<chr_col, psi_int_col>(span);
            next.offset = q.offset + lf.template extract_span<chr_col, psi_off_col>(span);
            return next;
        }
        void prefetch(ulint i) const { lf.prefetch(i); }
        void prefetch_rows(ulint lo, ulint hi) const { lf.prefetch_rows(lo, hi); }
        void prefetch_psi(ulint i) const { lf.prefetch(i); }
    };
    PackedAccess packed_access() const { return PackedAccess{lf_.get_reader(), f_rows_}; }

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
        FLPos psi_at_tail(ulint k) const { return idx->psi_at_tail(k); }
        FLPos psi_at_head(ulint k) const { return idx->psi_at_head(k); }
        FLPos psi_step(FLPos q, uchar& c) const { return idx->psi_step(q, c); }
        void prefetch(ulint i) const { idx->prefetch(i); }
        void prefetch_rows(ulint lo, ulint hi) const { idx->prefetch_rows(lo, hi); }
        void prefetch_psi(ulint i) const { idx->prefetch_psi(i); }
    };
    ColumnAccess column_access() { return ColumnAccess{this}; }

    /** The layout, then one line per structure: intervals and column widths in bits. */
    void describe(std::ostream& os) const {
        auto row_bits = [](const auto& w) { size_t b = 0; for (auto x : w) b += x; return b; };
        auto widths = [&](const auto& w) { for (size_t i = 0; i < w.size(); ++i) os << (i ? "," : "") << int(w[i]); };
        os << "layout: " << layout_name(parts()) << (has_phi_ && !has_phi_inv_ ? " (without phi_inv)" : "")
           << "\n";
        const auto& w = lf_.get_widths();
        os << "LF: intervals=" << lf_.intervals() << " runs=" << bwt_runs_ << " n=" << lf_.domain()
           << " widths(len,ptr,off,chr,psi_int,psi_off,phi_int,phi_off)=";
        widths(w);
        os << " row_bits=" << row_bits(w) << "\n";
        double total = double(row_bits(w)) * lf_.intervals();
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
        os << "total: " << total / 8 / 1e6 << " MB, " << total / bwt_runs_ << " bits/run\n";
    }

private:
    static constexpr char MAGIC[4] = {'T', 'M', 'S', 'X'};
    static constexpr uint32_t VERSION = 6;
    // Header flags for the parts an index holds.
    static constexpr uint32_t FLAG_PHI = 1, FLAG_PHI_INV = 2, FLAG_PSI = 4;

    LF lf_;
    Phi phi_;
    PhiInv phi_inv_;
    bool has_psi_ = false, has_phi_ = false, has_phi_inv_ = false;
    // Start tables of phi and phi_inv, built when phi_inv is present.
    StartTable phi_table_, phi_inv_table_;
    std::array<bool, 256> occurs_{};
    // With psi, the first row of each F block after the first (as
    // PackedAccess::f_rows), each block's byte, and each byte's block.
    std::array<ulint, F_BLOCKS> f_rows_{};
    std::array<uchar, F_BLOCKS + 1> f_byte_{};
    std::array<uchar, 256> f_block_of_{};
    // Maximal runs of equal characters in LF's rows: the BWT runs.
    ulint bwt_runs_ = 0;

    static constexpr size_t col(TmsLFCols c) { return static_cast<size_t>(c); }

    static uint32_t parts_flags(const TmsParts& p) {
        return (p.phi ? FLAG_PHI : 0) | (p.phi_inv ? FLAG_PHI_INV : 0) | (p.psi ? FLAG_PSI : 0);
    }

    /** A stream buffer that counts the bytes written to it and keeps none. */
    struct CountingBuf : std::streambuf {
        size_t count = 0;
        int_type overflow(int_type c) override {
            if (!traits_type::eq_int_type(c, traits_type::eof())) ++count;
            return traits_type::not_eof(c);
        }
        std::streamsize xsputn(const char*, std::streamsize n) override {
            count += static_cast<size_t>(n);
            return n;
        }
    };

    /**
     * Write x after its serialized size in bytes, as a uint64_t.  Returns the
     * bytes written.  The size is counted from a first serialization rather
     * than taken from serialize()'s return value, which for some Orbit
     * structures leaves out a few bytes.
     */
    template <typename Structure>
    static size_t serialize_sized(Structure& x, std::ostream& out) {
        CountingBuf counter;
        std::ostream sizing(&counter);
        x.serialize(sizing);
        const uint64_t size = counter.count;
        out.write(reinterpret_cast<const char*>(&size), sizeof(size));
        x.serialize(out);
        return sizeof(size) + size;
    }

    /**
     * Read a structure written by serialize_sized into x if want, else seek
     * past it.  Throws if a read structure's size disagrees with the one
     * recorded.
     */
    template <typename Structure>
    static void load_sized(Structure& x, std::istream& in, bool want) {
        uint64_t size = 0;
        in.read(reinterpret_cast<char*>(&size), sizeof(size));
        if (!in.good()) throw std::runtime_error("truncated tms index");
        if (!want) {
            in.seekg(static_cast<std::streamoff>(size), std::ios::cur);
            if (!in.good()) {
                // A stream that cannot seek reads past the structure instead.
                in.clear();
                in.ignore(static_cast<std::streamsize>(size));
            }
            return;
        }
        const std::streampos at = in.tellg();
        x.load(in);
        if (at != std::streampos(-1) && in.good() && static_cast<uint64_t>(in.tellg() - at) != size)
            throw std::runtime_error("malformed tms index: a structure's size disagrees with its header");
    }

    /**
     * load() with the parts in *need, or every part the index holds if need
     * is null.  phi_inv, whose start tables take time to build, is read only
     * with phi.
     */
    void load(std::istream& in, const TmsParts* need) {
        char magic[4] = {};
        uint32_t v = 0, flags = 0;
        in.read(magic, 4);
        in.read(reinterpret_cast<char*>(&v), sizeof(v));
        in.read(reinterpret_cast<char*>(&flags), sizeof(flags));
        if (!in.good() || std::memcmp(magic, MAGIC, 4) != 0) throw std::runtime_error("not a tms index");
        if (v != VERSION) throw std::runtime_error("unsupported tms index version " + std::to_string(v));
        const TmsParts held{(flags & FLAG_PSI) != 0, (flags & FLAG_PHI) != 0, (flags & FLAG_PHI_INV) != 0};
        const TmsParts want = need ? *need : held;
        check_parts(held, want);
        lf_.load(in);
        if (held.psi) {
            for (ulint& s : f_rows_) {
                uint64_t x = 0;
                in.read(reinterpret_cast<char*>(&x), sizeof(x));
                s = x;
            }
        }
        has_psi_ = want.psi;
        has_phi_ = want.phi;
        has_phi_inv_ = want.phi && want.phi_inv;
        if (held.phi) load_sized(phi_, in, has_phi_);
        if (held.phi_inv && has_phi_inv_) load_sized(phi_inv_, in, true);
        if (!in.good()) throw std::runtime_error("truncated tms index");
        compute_occurs();
        if (has_phi_inv_) build_start_tables();
    }

    /** Throw, naming what is missing and how to build it, if an index with parts held lacks one in want. */
    static void check_parts(const TmsParts& held, const TmsParts& want) {
        std::string missing;
        if (want.psi && !held.psi) missing = "psi (the PSI columns)";
        else if (want.phi && !held.phi) missing = "phi (phi and the PHI columns)";
        else if (want.phi_inv && !held.phi_inv) missing = "phi_inv";
        else return;
        // The smallest layout with every part the query needs.
        std::string build = want.psi && want.phi ? "full" : want.psi ? "psi or full" : "phi or full";
        if (want.phi) build += ", with --minima or --lcp-bin";
        if (want.phi_inv) build += " and without --no-phi-inv";
        throw std::runtime_error("this query needs " + missing + ", which this index (layout " + layout_name(held) +
                                 (held.phi && !held.phi_inv ? ", without phi_inv" : "") +
                                 ") lacks; build one with tms-build --layout " + build);
    }

    // Derived fields: which characters occur, the BWT runs, and each F
    // block's byte and each byte's block.
    void compute_occurs() {
        occurs_.fill(false);
        bwt_runs_ = 0;
        for (ulint i = 0, prev = 256; i < lf_.intervals(); ++i) {
            const uchar c = lf_.get_character(i);
            bwt_runs_ += c != prev;
            occurs_[c] = true;
            prev = c;
        }
        // F blocks are in byte order.
        f_block_of_.fill(no_code);
        f_byte_.fill(0);
        uchar b = 0;
        for (size_t c = 0; c < 256; ++c)
            if (occurs_[c]) {
                if (b > F_BLOCKS) throw std::runtime_error("more characters than the F block table holds");
                f_byte_[b] = static_cast<uchar>(c);
                f_block_of_[c] = b++;
            }
    }

    /**
     * The rows of the psi and full layouts (see the file comment) of the
     * RLBWT with the given run heads and lengths: the union of the starts of
     * LF's intervals, split with sp, and of their images.  Each row gets the
     * byte of the run it lies in.  O(r) time and space.
     */
    static void union_rows(const std::vector<uchar>& heads, const std::vector<ulint>& lens,
                           const orbit::split_params& sp, std::vector<uchar>& row_heads,
                           std::vector<ulint>& row_lens) {
        using Enc = orbit::rlbwt::rlbwt_interval_encoding<>;
        const Enc split = Enc::lf_interval_encoding(heads, lens, sp);
        const ulint m = split.intervals();
        // Call f(k, c) for each split interval k in row order, with c its byte.
        auto each = [&](auto f) {
            ulint k = 0;
            for (size_t r = 0; r < lens.size(); ++r)
                for (ulint covered = 0; covered < lens[r]; ++k) {
                    f(k, heads[r]);
                    covered += split.get_length(k);
                }
        };
        // Byte c's F block holds the images of c's intervals in their order,
        // so counting rows and intervals per byte places every image.
        std::array<ulint, 256> rows_of{}, intervals_of{}, next_image{}, next_slot{};
        each([&](ulint k, uchar c) {
            rows_of[c] += split.get_length(k);
            ++intervals_of[c];
        });
        for (size_t c = 0, a = 0, b = 0; c < 256; ++c) {
            next_image[c] = a;
            a += rows_of[c];
            next_slot[c] = b;
            b += intervals_of[c];
        }
        // The images' starts, in row order.
        std::vector<ulint> images(m);
        each([&](ulint k, uchar c) {
            images[next_slot[c]++] = next_image[c];
            next_image[c] += split.get_length(k);
        });
        // Merge the two sorted lists of starts, calling emit(c, length) for
        // each row.
        auto merge = [&](auto emit) {
            ulint start = 0;
            size_t j = 0;
            each([&](ulint k, uchar c) {
                const ulint end = start + split.get_length(k);
                ulint at = start;
                while (j < m && images[j] <= start) ++j;
                for (; j < m && images[j] < end; ++j) {
                    emit(c, images[j] - at);
                    at = images[j];
                }
                emit(c, end - at);
                start = end;
            });
        };
        size_t rows = 0;
        merge([&](uchar, ulint) { ++rows; });
        row_heads.clear();
        row_lens.clear();
        row_heads.reserve(rows);
        row_lens.reserve(rows);
        merge([&](uchar c, ulint len) {
            row_heads.push_back(c);
            row_lens.push_back(len);
        });
    }

    /**
     * For LF rows with the given heads and lengths, fill the PSI columns
     * with psi's pointer and offset for each row's head, and the table of F
     * block starts.  Row i with byte c holds the occurrences of c ranked R
     * to R + length - 1 in the BWT, where R counts those in earlier rows, so
     * psi maps F positions C[c] + R onward into it; the rows of c's F block
     * that start there are found by walking the block alongside.
     */
    void fill_psi(const std::vector<uchar>& row_heads, const std::vector<ulint>& row_lens,
                  std::vector<typename LF::data_tuple>& cols) {
        const size_t rows = row_lens.size();
        std::array<ulint, 256> count{}, f_pos{}, f_row{};
        for (size_t i = 0; i < rows; ++i) count[row_heads[i]] += row_lens[i];
        for (size_t c = 0, a = 0; c < 256; ++c) {
            f_pos[c] = a;
            a += count[c];
        }
        // Each F block's first row, and the table of blocks after the first.
        f_rows_.fill(std::numeric_limits<ulint>::max());
        {
            ulint start = 0;
            size_t i = 0, b = 0;
            for (size_t c = 0; c < 256; ++c) {
                if (count[c] == 0) continue;
                while (start < f_pos[c]) start += row_lens[i++];
                if (start != f_pos[c]) throw std::logic_error("an F block must start at a row start");
                f_row[c] = i;
                if (b > F_BLOCKS) throw std::runtime_error("more characters than the F block table holds");
                if (b > 0) f_rows_[b - 1] = i;
                ++b;
            }
        }
        std::array<ulint, 256> rank{};
        std::array<ulint, 256> next_row = f_row, next_start = f_pos;
        for (size_t i = 0; i < rows; ++i) {
            const uchar c = row_heads[i];
            const ulint lo = f_pos[c] + rank[c], hi = lo + row_lens[i];
            for (; next_start[c] < hi; next_start[c] += row_lens[next_row[c]++]) {
                cols[next_row[c]][col(TmsLFCols::PSI_INT)] = i;
                cols[next_row[c]][col(TmsLFCols::PSI_OFF)] = next_start[c] - lo;
            }
            rank[c] += row_lens[i];
        }
    }

    /** The F block of row i: how many F blocks after the first start at or before it. */
    uchar f_block(ulint i) const {
        uchar b = 0;
        for (ulint s : f_rows_) b += i >= s;
        return b;
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
    void build_phi(const std::vector<ulint>& lens, const Enc& lf_enc,
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
            orbit::rlbwt::lf_move<false> lf(lf_enc);
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
