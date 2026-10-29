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
 *    Queries read phi points only at rows that start a BWT run: row 0, and
 *    the candidates of a reposition, the head of the nearest run of the
 *    wanted character below and the row after the tail of the nearest one
 *    above.  phi's intervals start at the text positions of run heads, so
 *    such a phi point has offset 0: PHI_INT holds the phi interval that
 *    starts there (0 at rows that start no run), and PHI_OFF has width 0.
 *    The phi point of a tail row is one phi step from the next row's head's
 *    phi point, so a reposition reads both candidates' start points from
 *    rows its scan has already read.
 *  - Optionally, with phi, two LCP summary columns per LF row (see
 *    TmsLcpSummary): TOP_LCP, the LCP of the row's head with the row above,
 *    and INT_MIN, the least LCP of the row's other rows.  They let a phi
 *    reposition (--mode phisum) take whole rows' minima from the rows its
 *    scan reads, walking phi only over the current row.
 *  - With psi, a table of F block starts.  An F block starts at a row start
 *    (the image of the first run of its character), so a row's F character
 *    is the block its index falls in.
 *  - Optionally, phi over text positions, with an integrated PLCP column
 *    holding PLCP at each interval's start, in two forms over the same
 *    intervals: lengths-based (rows store lengths) and starts-based (rows
 *    store absolute starts).  A phi point of the lengths-based form is an
 *    interval and offset, enough for matching statistics; the starts-based
 *    form also gives its text position, which positions and SMEM reports
 *    need, at the cost of a wider row.  The file holds both, and a load
 *    reads only the one its query needs.
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
 * in O(r) time and space.  phi and phi_inv come either from per-run suffix
 * array and LCP samples (TmsPhiSamples, as TeraLCP -ophi writes them), in
 * O(r log r) time and O(r) space, or from the LCP value at each run head and
 * a walk of LF over all n rows, which needs a BWT with one terminator.  Either
 * way the starts-based form of phi is built first, and the lengths-based form
 * is copied from its rows at the end of the build.  LF and phi each take their
 * own splitting parameters.
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
 * kept stay where the full layout has them.  An index is written as an
 * ORBITIDX file (see index_format.hpp), with a section per structure whose
 * header entry gives its columns and their widths.  load() reads only the
 * sections holding the parts a query needs (TmsParts, see sections_for),
 * seeking past the rest.  LF is loaded whole, with every column the layout
 * keeps, whichever parts a query needs.
 *
 * An index built from a ropebwt3 FMD may also hold the FMD's sequence-ID
 * table (SeqIdTable, the SEQ_IDS section), which names the contigs its
 * text's strings come from.  It is not a part a query needs: a query that
 * asks for it gets it if the index holds it.
 */

#ifndef _TMS_INDEX_HPP
#define _TMS_INDEX_HPP

#include "orbit/rlbwt.hpp"
#include "orbit/common.hpp"
#include "ms_alphabet.hpp"
#include "mem_report.hpp"
#include "step_stats.hpp"
#include "index_format.hpp"
#include "seq_ids.hpp"
#include <algorithm>
#include <array>
#include <cstdint>
#include <numeric>
#include <cstring>
#include <fstream>
#include <iostream>
#include <type_traits>
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

enum class TmsLFCols { PSI_INT, PSI_OFF, PHI_INT, PHI_OFF, TOP_LCP, INT_MIN, COUNT };
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
    // With phi, phi's starts-based form rather than its lengths-based form.
    // An index that holds phi holds both forms.
    bool phi_starts = false;
    // The sequence-ID table, if the index holds one.
    bool seq_ids = false;
};

/**
 * Suffix array and LCP samples from which a TmsIndex builds phi and phi_inv
 * without walking LF, in the record order of TeraLCP's -ophi output
 * (src/TeraLCP/PHI_EXPORT_FORMAT.md in TeraTools): one record per BWT run,
 * except that a run of the terminator (orbit::TERMINATOR, the end of each
 * string) has one record per row.  Record i holds, for its rows s to e, the
 * text positions head[i] = SA[s] and tail[i] = SA[e], and top[i] = LCP[s], the
 * LCP with the row above (0 for row 0), where a common prefix never extends
 * over a terminator.  The three are packed.
 */
struct TmsPhiSamples {
    orbit::int_vector head, tail, top;
    TmsPhiSamples() = default;
    /** Room for records records of a text of length n whose tops are at most max_top. */
    TmsPhiSamples(size_t records, ulint n, ulint max_top)
        : head(records, orbit::bit_width(n - 1)), tail(records, orbit::bit_width(n - 1)),
          top(records, orbit::bit_width(max_top)) {}
    size_t size() const { return head.size(); }
    /** The number of records an RLBWT with these run heads and lengths needs. */
    template <typename Lens>
    static ulint records_for(const std::vector<uchar>& heads, const Lens& lens) {
        ulint records = 0;
        for (size_t r = 0; r < heads.size(); ++r) records += heads[r] == orbit::TERMINATOR ? ulint(lens[r]) : 1;
        return records;
    }
};

/**
 * LCP summaries of each LF row R of a TmsIndex, in row order: top[R] is the
 * LCP at R's head row and int_min[R] the least LCP of its other rows, or a
 * value at least every LCP in the index if R has one row.  They come from
 * each run's retained LCP pairs (as a TeraLCP minima file holds them), cut
 * into the LF rows as ms's build cuts them (ms_build::row_summaries), so a
 * row inside a run gets the top max(u, d) of apply_lf_splitting and an
 * interior minimum over the pairs it keeps.  Every value is an LCP inside a
 * reposition's range or at least that range's minimum, and the range's
 * minimum is among them or the LCPs the walk reads, so the minima a
 * reposition takes are exact.
 */
struct TmsLcpSummary {
    orbit::int_vector top, int_min;
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
    // With phi, the LCP summary columns' values, one per LF row (see
    // lf_row_lengths); null builds the index without them.
    const TmsLcpSummary* lcp_summary = nullptr;
    // Print memory at the phase boundaries of the build (see mem_report).
    bool report = false;
};

class TmsIndex {
public:
    using LF = orbit::rlbwt::lf_permutation<TmsLFCols, true, false, ms_alpha::alphabet>;
    using Phi = orbit::rlbwt::phi_permutation_impl<TmsPhiCols, true, true, orbit::move_vector>;
    using PhiInv = orbit::rlbwt::phi_inv_permutation_impl<TmsPhiInvCols, true, true, orbit::move_vector>;
    // phi's lengths-based form, over the same intervals as Phi.
    using PhiLen = orbit::permutation_impl<TmsPhiCols, true, false, false, orbit::move_columns, orbit::move_structure,
                                           orbit::move_vector>;
    using LFPos = typename LF::position;
    // A psi point: an LF position whose F character and psi step psi reads.
    using FLPos = LFPos;
    using PhiPos = typename Phi::position;
    using PhiInvPos = typename PhiInv::position;
    // A phi point of the lengths-based form: an interval and offset, without a text position.
    using PhiLenPos = typename PhiLen::position;
    template <bool Starts>
    using PhiPosOf = std::conditional_t<Starts, PhiPos, PhiLenPos>;
    using position = LFPos;
    // The most F blocks after the first: one per character of the
    // alphabet but the smallest.
    static constexpr size_t F_BLOCKS = ms_alpha::alphabet::SIGMA - 1;

    TmsIndex() = default;

    /**
     * Build from run heads (index codes: orbit::TERMINATOR, orbit::SEPARATOR,
     * a nucleotide or ms_alpha::TEXT_N) and run lengths.  There may be many
     * terminators (one per string of a collection, as in an FMD), so that
     * LF has several cycles; LF and psi do not depend on that.  If run_tops
     * is given, it holds the LCP value at each run's head row (0 for row 0),
     * and the index also gets phi (and phi_inv if opts.phi_inv); this walks
     * LF over all n rows from row 0 and so needs a single terminator.
     * opts.layout says which parts to build: the psi layout takes no
     * run_tops, and the phi layout needs them.  opts.lf_split must split.
     * lens and run_tops may be std::vectors or packed orbit::int_vectors.
     */
    template <typename Lens = std::vector<ulint>, typename Tops = std::vector<ulint>>
    TmsIndex(const std::vector<uchar>& heads, const Lens& lens,
             const TmsBuildOptions& opts = TmsBuildOptions{}, const Tops* run_tops = nullptr) {
        build(heads, lens, opts, run_tops, nullptr);
    }

    /**
     * Build as above, with phi (and phi_inv if opts.phi_inv) from samples,
     * which must have TmsPhiSamples::records_for(heads, lens) records, instead
     * of from run tops and a walk of LF.  The BWT may have many terminators.
     * opts.layout must not be psi.  samples is freed during the build.
     */
    template <typename Lens = std::vector<ulint>>
    TmsIndex(const std::vector<uchar>& heads, const Lens& lens, const TmsBuildOptions& opts,
             TmsPhiSamples&& samples) {
        build(heads, lens, opts, static_cast<const std::vector<ulint>*>(nullptr), &samples);
    }

private:
    template <typename Lens, typename Tops>
    void build(const std::vector<uchar>& heads, const Lens& lens, const TmsBuildOptions& opts,
               const Tops* run_tops, TmsPhiSamples* samples) {
        using Enc = orbit::rlbwt::rlbwt_interval_encoding<ms_alpha::alphabet>;
        const bool lcp = run_tops || samples;
        if (opts.layout == TmsLayout::PSI && lcp) throw std::invalid_argument("a psi layout takes no LCP values");
        if (opts.layout == TmsLayout::PHI && !lcp) throw std::invalid_argument("a phi layout needs LCP values");
        check_lf_split(opts.lf_split);
        auto phase = [&](const char* name) { if (opts.report) mem_report::phase(name, heads.size()); };
        has_psi_ = opts.layout != TmsLayout::PHI;
        Enc lf_enc;
        // With psi, the rows' bytes, which fill_psi reads.
        std::vector<uchar> row_heads;
        orbit::int_vector row_lens;
        if (has_psi_) {
            union_rows(heads, lens, opts.lf_split, row_heads, row_lens);
            phase("tms_union_rows");
            lf_enc = Enc::lf_interval_encoding(row_heads, row_lens, orbit::NO_SPLITTING);
        } else {
            lf_enc = Enc::lf_interval_encoding(heads, lens, opts.lf_split);
        }
        phase("tms_lf_enc");
        // With phi, the PHI_INT column of each LF row (see the file comment).
        orbit::packed_matrix<1> phi_cols;
        if (run_tops) {
            if (run_tops->size() != heads.size()) throw std::invalid_argument("run_tops must have one value per run");
            // build_phi frees LF's encoding once its walk is done, and it is
            // made again afterwards, from the rows with psi (which is cheap)
            // or by splitting the runs again, when less is held.
            auto release = [&]() { lf_enc = Enc(); };
            build_phi(lens, lf_enc, *run_tops, opts.phi_split, opts.phi_inv, phi_cols, phase, release);
            lf_enc = has_psi_ ? Enc::lf_interval_encoding(row_heads, row_lens, orbit::NO_SPLITTING)
                              : Enc::lf_interval_encoding(heads, lens, opts.lf_split);
            phase("tms_build_phi");
            has_phi_ = has_phi_starts_ = has_phi_len_ = true;
            has_phi_inv_ = opts.phi_inv;
        } else if (samples) {
            // As above, without a walk: LF's encoding is needed only to find
            // the LF intervals that start runs.
            auto release = [&]() { lf_enc = Enc(); };
            build_phi_from_samples(heads, lens, lf_enc, *samples, opts.phi_split, opts.phi_inv, phi_cols, phase,
                                   release);
            lf_enc = has_psi_ ? Enc::lf_interval_encoding(row_heads, row_lens, orbit::NO_SPLITTING)
                              : Enc::lf_interval_encoding(heads, lens, opts.lf_split);
            phase("tms_build_phi");
            has_phi_ = has_phi_starts_ = has_phi_len_ = true;
            has_phi_inv_ = opts.phi_inv;
        }
        // The data columns take the widths of their largest values, and the
        // columns of parts the index lacks take no bits and read as 0.  LF is
        // built with them and its data written in place.
        std::array<uchar, static_cast<size_t>(TmsLFCols::COUNT)> widths{};
        if (has_psi_) {
            ulint max_int = 0, max_off = 0;
            fill_psi(row_heads, [&](ulint i) { return lf_enc.get_length(i); }, [&](ulint, ulint to, ulint off) {
                max_int = std::max(max_int, to);
                max_off = std::max(max_off, off);
            });
            widths[col(TmsLFCols::PSI_INT)] = orbit::bit_width(max_int);
            widths[col(TmsLFCols::PSI_OFF)] = orbit::bit_width(max_off);
        }
        if (has_phi_) {
            ulint max_int = 0;
            for (size_t k = 0; k < phi_cols.size(); ++k) max_int = std::max(max_int, phi_cols.get<0>(k));
            widths[col(TmsLFCols::PHI_INT)] = orbit::bit_width(max_int);
            // Every phi point read has offset 0, so PHI_OFF takes no bits.
            widths[col(TmsLFCols::PHI_OFF)] = 0;
        }
        if (const TmsLcpSummary* sum = opts.lcp_summary) {
            if (!has_phi_) throw std::invalid_argument("LCP summaries need an index with phi");
            if (sum->top.size() != lf_enc.intervals() || sum->int_min.size() != lf_enc.intervals())
                throw std::invalid_argument("LCP summaries must have one value per LF row");
            // At least one bit each, so that a loaded index can tell it has them.
            ulint max_top = 0, max_int = 0;
            for (size_t k = 0; k < sum->top.size(); ++k) {
                max_top = std::max<ulint>(max_top, sum->top[k]);
                max_int = std::max<ulint>(max_int, sum->int_min[k]);
            }
            widths[col(TmsLFCols::TOP_LCP)] = std::max<uchar>(1, orbit::bit_width(max_top));
            widths[col(TmsLFCols::INT_MIN)] = std::max<uchar>(1, orbit::bit_width(max_int));
        }
        lf_ = LF(lf_enc, widths, LF::in_place_data);
        lf_enc = Enc();
        row_lens = orbit::int_vector();
        if (has_psi_) {
            fill_psi(row_heads, [&](ulint i) { return lf_.get_length(i); }, [&](ulint k, ulint to, ulint off) {
                lf_.template set_data<TmsLFCols::PSI_INT>(k, to);
                lf_.template set_data<TmsLFCols::PSI_OFF>(k, off);
            });
            std::vector<uchar>().swap(row_heads);
        }
        if (has_phi_) {
            for (size_t k = 0; k < phi_cols.size(); ++k) lf_.template set_data<TmsLFCols::PHI_INT>(k, phi_cols.get<0>(k));
            phi_cols = orbit::packed_matrix<1>();
        }
        if (const TmsLcpSummary* sum = opts.lcp_summary) {
            for (size_t k = 0; k < sum->top.size(); ++k) {
                lf_.template set_data<TmsLFCols::TOP_LCP>(k, sum->top[k]);
                lf_.template set_data<TmsLFCols::INT_MIN>(k, sum->int_min[k]);
            }
        }
        phase("tms_assemble");
        if (has_phi_len_) {
            // phi's lengths-based form comes from the finished starts-based
            // form, once LF is built and the build's other data are freed.
            phi_len_ = lengths_form(phi_);
            phase("tms_phi_len");
        }
        compute_occurs();
        if (has_phi_inv_) build_start_tables();
    }

public:

    /**
     * The lengths of the LF rows that an index built with heads, lens and
     * opts has, in row order: the union rows with psi, else LF's split
     * intervals.  An LCP summary (TmsBuildOptions::lcp_summary) has one value
     * per row.
     */
    template <typename Lens>
    static orbit::int_vector lf_row_lengths(const std::vector<uchar>& heads, const Lens& lens,
                                            const TmsBuildOptions& opts) {
        check_lf_split(opts.lf_split);
        if (opts.layout != TmsLayout::PHI) {
            std::vector<uchar> row_heads;
            orbit::int_vector row_lens;
            union_rows(heads, lens, opts.lf_split, row_heads, row_lens);
            return row_lens;
        }
        using Enc = orbit::rlbwt::rlbwt_interval_encoding<ms_alpha::alphabet>;
        const Enc enc = Enc::lf_interval_encoding(heads, lens, opts.lf_split);
        ulint max_len = 0;
        for (ulint k = 0; k < enc.intervals(); ++k) max_len = std::max(max_len, enc.get_length(k));
        orbit::int_vector rows(enc.intervals(), orbit::bit_width(max_len));
        for (ulint k = 0; k < enc.intervals(); ++k) rows.set(k, enc.get_length(k));
        return rows;
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
     * Write the index as an ORBITIDX file (see index_format.hpp): a header
     * with a section per structure, then LF, with psi the F block starts, and
     * with phi both of its forms and phi_inv if present.  Throws if the index was loaded with only
     * some of its parts.
     */
    size_t serialize(std::ostream& out) {
        if (!complete_)
            throw std::logic_error("an index loaded with only some of its parts cannot be written");
        using namespace orbit_idx;
        Header h;
        h.tool = Tool::TMS;
        h.layout = has_psi_ && has_phi_ ? Layout::FULL : has_psi_ ? Layout::PSI : has_phi_ ? Layout::PHI : Layout::NONE;
        h.n = lf_.domain();
        h.bwt_runs = bwt_runs_;
        h.info = info_;
        std::vector<SectionSource> srcs;
        srcs.push_back(move_table_source(Section::LF, lf_, LF_COLS));
        if (has_psi_) {
            SectionEntry e;
            e.id = Section::F_STARTS;
            e.kind = Kind::U64_ARRAY;
            e.rows = F_BLOCKS;
            srcs.push_back({e, [this](std::ostream& o) {
                                for (ulint s : f_rows_) {
                                    const uint64_t x = s;
                                    o.write(reinterpret_cast<const char*>(&x), sizeof(x));
                                }
                            }});
        }
        if (has_phi_) {
            srcs.push_back(move_table_source(Section::PHI_LEN, phi_len_, PHI_LEN_COLS));
            srcs.push_back(move_table_source(Section::PHI_START, phi_, PHI_COLS));
        }
        if (has_phi_inv_) srcs.push_back(move_table_source(Section::PHI_INV, phi_inv_, PHI_INV_COLS));
        if (seq_ids_) srcs.push_back(seq_ids_source(*seq_ids_));
        return write_file(out, h, srcs);
    }

    /** A SEQ_IDS section source for t, which must outlive the write. */
    static orbit_idx::SectionSource seq_ids_source(const SeqIdTable& t) {
        orbit_idx::SectionEntry e;
        e.id = orbit_idx::Section::SEQ_IDS;
        e.kind = orbit_idx::Kind::SEQ_TABLE;
        e.rows = t.names.size();
        return {e, [&t](std::ostream& o) { t.serialize(o); }};
    }

    /**
     * An empty string if t can describe a text of length n with the given
     * number of terminators, both strands of each contig, each ending in a
     * terminator; else what differs.  Without the count, terminators is 0
     * and only the length is checked.
     */
    static std::string seq_ids_mismatch(const SeqIdTable& t, ulint n, ulint terminators = 0) {
        if (t.text_length() != n)
            return "the table describes a text of length " + std::to_string(t.text_length()) + " (" +
                   std::to_string(t.names.size()) + " sequences, both strands) but the index's has length " +
                   std::to_string(n);
        if (terminators != 0 && terminators != 2 * t.names.size())
            return "the table lists " + std::to_string(t.names.size()) + " sequences, but the index's text has " +
                   std::to_string(terminators) + " strings rather than two per sequence";
        return "";
    }

    /** The terminators in the BWT, one per string of the text, counted from LF. */
    ulint terminators() {
        ulint t = 0;
        for (ulint i = 0; i < lf_.intervals(); ++i)
            if (lf_.get_character(i) == orbit::TERMINATOR) t += lf_.get_length(i);
        return t;
    }

    /** Hold t as the index's sequence-ID table; throws if it does not describe the text. */
    void set_seq_ids(SeqIdTable t) {
        const std::string bad = seq_ids_mismatch(t, lf_.domain(), terminators());
        if (!bad.empty()) throw std::invalid_argument(bad);
        seq_ids_ = std::move(t);
    }

    /** The sequence-ID table, or null if the index holds none or it was not loaded. */
    const SeqIdTable* seq_ids() const { return seq_ids_ ? &*seq_ids_ : nullptr; }

    /** Read an index written by serialize(), with every part it holds.  Throws on malformed input. */
    void load(std::istream& in) { load(in, nullptr); }

    /**
     * Read an index written by serialize(), with only LF and the parts in
     * need: the sections sections_for(need), skipping the rest.  Throws,
     * before reading LF, if the index lacks a part in need.
     */
    void load(std::istream& in, const TmsParts& need) { load(in, &need); }

    /**
     * The sections load() reads for need, as a mask of orbit_idx::section_bit:
     * LF; with psi the F block starts; with phi one form of it, starts-based
     * if need.phi_starts, else lengths-based; and with phi and phi_inv,
     * phi_inv, from which load() also builds start tables (see
     * start_table_bytes).
     */
    static uint32_t sections_for(const TmsParts& need) {
        using orbit_idx::Section;
        using orbit_idx::section_bit;
        uint32_t m = section_bit(Section::LF);
        if (need.psi) m |= section_bit(Section::F_STARTS);
        if (need.phi) m |= section_bit(need.phi_starts ? Section::PHI_START : Section::PHI_LEN);
        if (need.phi && need.phi_inv) m |= section_bit(Section::PHI_INV);
        if (need.seq_ids) m |= section_bit(Section::SEQ_IDS);
        return m;
    }

    /** A short name for what sections_for(need) loads: lf, or its parts joined with +. */
    static std::string load_name(const TmsParts& need) {
        std::string s;
        auto add = [&](const char* p) { s += s.empty() ? p : std::string("+") + p; };
        if (need.psi) add("psi");
        if (need.phi) add(need.phi_starts ? "phi-starts" : "phi-len");
        if (need.phi && need.phi_inv) add("phi-inv");
        if (need.seq_ids) add("seq-ids");
        return s.empty() ? "lf" : s;
    }

    /** Build provenance written with the index, as "key=value" lines. */
    const std::string& info() const { return info_; }
    void set_info(std::string info) { info_ = std::move(info); }

    bool has_psi() const { return has_psi_; }
    bool has_phi() const { return has_phi_; }
    bool has_phi_inv() const { return has_phi_inv_; }
    /** Whether LF holds the LCP summary columns (see TmsLcpSummary). */
    bool has_lcp_summary() const { return lf_.get_widths()[LF::template data_column<TmsLFCols::TOP_LCP>()] > 0; }
    /** The LCP summaries of LF row k: the LCP at its head row, and the least of its other rows'. */
    ulint top_lcp(ulint k) const { return lf_.template get<TmsLFCols::TOP_LCP>(k); }
    ulint int_min(ulint k) const { return lf_.template get<TmsLFCols::INT_MIN>(k); }
    /** Whether phi's starts-based form, and its lengths-based form, are present. */
    bool has_phi_starts() const { return has_phi_starts_; }
    bool has_phi_lengths() const { return has_phi_len_; }
    /**
     * The parts the index holds, or those loaded if load() skipped some;
     * phi_starts says whether the starts-based form of phi is present.
     */
    TmsParts parts() const { return TmsParts{has_psi_, has_phi_, has_phi_inv_, has_phi_starts_}; }
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
    LFPos finish_LF(LFPos p) const {
        const LFPos q = lf_.finish_next(p);
        MS_FF(FF_LF, q.interval - p.interval);
        return q;
    }
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
    FLPos finish_psi(FLPos q) const {
        const FLPos r = lf_.finish_next(q);
        MS_FF(FF_PSI, r.interval - q.interval);
        return r;
    }
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

    // phi side.  The text position of a phi point of the starts-based form
    // is its idx.  The functions that take a phi point serve either form,
    // and the others take the form as Starts.
    /** phi point of the text position of LF interval k's head row, unresolved. */
    template <bool Starts = true>
    PhiPosOf<Starts> phi_at_head(ulint k) const {
        PhiPosOf<Starts> p;
        p.interval = lf_.template get<TmsLFCols::PHI_INT>(k);
        p.offset = lf_.template get<TmsLFCols::PHI_OFF>(k);
        if constexpr (Starts) p.idx = 0;  // set by resolve_phi
        return p;
    }
    /** Resolve a phi point whose offset may run past its interval, filling in its text position. */
    PhiPos resolve_phi(PhiPos p) const { return phi_.finish_next(p); }
    PhiLenPos resolve_phi(PhiLenPos p) const { return phi_len_.finish_next(p); }
    PhiPos phi(PhiPos p) { return phi_.phi(p); }
    PhiPos start_phi(PhiPos p) const { return phi_.start_next(p); }
    PhiLenPos start_phi(PhiLenPos p) const { return phi_len_.start_next(p); }
    PhiPos finish_phi(PhiPos p) const { return phi_.finish_next(p); }
    PhiLenPos finish_phi(PhiLenPos p) const { return phi_len_.finish_next(p); }
    template <bool Starts = true>
    void prefetch_phi(ulint i) const {
        if constexpr (Starts) phi_.prefetch(i);
        else phi_len_.prefetch(i);
    }
    /** The phi point one text position to the left, wrapping at position 0. */
    PhiPos phi_left(PhiPos p) {
        if (p.offset > 0) {
            --p.offset;
            --p.idx;
            return p;
        }
        return phi_.up(p);
    }
    PhiLenPos phi_left(PhiLenPos p) {
        if (p.offset > 0) {
            --p.offset;
            return p;
        }
        return phi_len_.up(p);
    }
    /** PLCP at a phi point: the LCP of its row with the row above. */
    ulint plcp(PhiPos p) const { return phi_.template get<TmsPhiCols::PLCP>(p.interval) - p.offset; }
    ulint plcp(PhiLenPos p) const { return phi_len_.template get<TmsPhiCols::PLCP>(p.interval) - p.offset; }
    /** phi's intervals, the same in both forms. */
    ulint phi_intervals() const { return has_phi_starts_ ? phi_.intervals() : phi_len_.intervals(); }
    /** The length of phi interval i, in whichever form is present. */
    ulint phi_length(ulint i) const { return has_phi_starts_ ? phi_.get_length(i) : phi_len_.get_length(i); }
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
     * Bytes of the start table that build_start_tables makes, with the
     * default shift, for a starts-based permutation over n text positions in
     * the given intervals: one entry per 2^shift positions, plus one, of the
     * bits an interval takes.
     */
    static size_t start_table_bytes(ulint n, ulint intervals) {
        const unsigned shift = start_table_shift(n, intervals);
        const size_t entries = ((n - 1) >> shift) + 2;
        const size_t width = std::max<uchar>(1, orbit::bit_width(intervals - 1));
        // As packed_matrix: whole bytes of packed entries, and a word of slack.
        return (entries * width + 7) / 8 + sizeof(ulint);
    }

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
        if (has_phi_starts_) phi_table_ = StartTable(phi_, pick(phi_));
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
#ifdef MS_STEP_STATS
            const ulint from = i;
#endif
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
            MS_FF((std::is_same_v<Perm, PhiInv> ? FF_PHI_INV : FF_PHI_WALK), i - from);
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
    // The other structures, for development tools that read their rows
    // directly; each is empty unless it was loaded.
    const Phi& phi_perm() const { return phi_; }
    const PhiLen& phi_len_perm() const { return phi_len_; }
    const PhiInv& phi_inv_perm() const { return phi_inv_; }
    /** Maximal runs of equal characters in LF's rows: the BWT runs. */
    ulint bwt_runs() const { return bwt_runs_; }

    /** The parts the index file at path holds, from its header.  Throws if it is not a tms index. */
    static TmsParts held_parts(const std::string& path) { return held_parts(orbit_idx::read_header_file(path)); }
    static TmsParts held_parts(const orbit_idx::Header& h) {
        using orbit_idx::Section;
        if (h.tool != orbit_idx::Tool::TMS) throw std::runtime_error("not a tms index");
        const bool phi = h.find(Section::PHI_LEN) && h.find(Section::PHI_START);
        return TmsParts{h.find(Section::F_STARTS) != nullptr, phi, phi && h.find(Section::PHI_INV) != nullptr, phi};
    }
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
        static constexpr size_t top_lcp_col = LF::template data_column<TmsLFCols::TOP_LCP>();
        static constexpr size_t int_min_col = LF::template data_column<TmsLFCols::INT_MIN>();
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
        ulint top_lcp(Row r) const { return lf.template get_at<top_lcp_col>(r); }
        ulint int_min(Row r) const { return lf.template get_at<int_min_col>(r); }
        /** Resolve an LF position whose offset may run past its interval's end; returns its row. */
        Row resolve_LF(LFPos& p) const {
#ifdef MS_STEP_STATS
            const ulint from = p.interval;
#endif
            Row r = row(p.interval);
            for (ulint len = length(r); p.offset >= len; len = length(r)) {
                p.offset -= len;
                r = row(++p.interval);
            }
            MS_FF(FF_LF, p.interval - from);
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
#ifdef MS_STEP_STATS
            const ulint from = q.interval;
#endif
            Row r = row(q.interval);
            for (ulint len = length(r); q.offset >= len; len = length(r)) {
                q.offset -= len;
                r = row(++q.interval);
            }
            MS_FF(FF_PSI, q.interval - from);
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
        ulint top_lcp(Row i) const { return idx->top_lcp(i); }
        ulint int_min(Row i) const { return idx->int_min(i); }
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
           << " widths(len,ptr,off,chr,psi_int,psi_off,phi_int,phi_off,top_lcp,int_min)=";
        widths(w);
        os << " row_bits=" << row_bits(w) << "\n";
        double total = double(row_bits(w)) * lf_.intervals();
        if (has_phi_len_) {
            const auto& pw = phi_len_.get_widths();
            os << "phi (lengths): intervals=" << phi_len_.intervals() << " widths(len,ptr,off,plcp)=";
            widths(pw);
            os << " row_bits=" << row_bits(pw) << "\n";
            total += double(row_bits(pw)) * phi_len_.intervals();
        }
        if (has_phi_starts_) {
            const auto& pw = phi_.get_widths();
            os << "phi (starts): intervals=" << phi_.intervals() << " widths(start,ptr,off,plcp)=";
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
    // Each structure's column codes, in the order its rows hold them.
    using Col = orbit_idx::Col;
    static constexpr std::array<Col, 10> LF_COLS{Col::LENGTH, Col::POINTER, Col::OFFSET, Col::CHARACTER,
                                                 Col::PSI_INT, Col::PSI_OFF, Col::PHI_INT, Col::PHI_OFF,
                                                 Col::TOP_LCP, Col::INT_MIN};
    static constexpr std::array<Col, 4> PHI_LEN_COLS{Col::LENGTH, Col::POINTER, Col::OFFSET, Col::PLCP};
    static constexpr std::array<Col, 4> PHI_COLS{Col::START, Col::POINTER, Col::OFFSET, Col::PLCP};
    static constexpr std::array<Col, 4> PHI_INV_COLS{Col::START, Col::POINTER, Col::OFFSET, Col::PLCPB};

    LF lf_;
    Phi phi_;
    PhiLen phi_len_;
    PhiInv phi_inv_;
    // has_phi_: the PHI columns and at least one form of phi.
    bool has_psi_ = false, has_phi_ = false, has_phi_inv_ = false, has_phi_starts_ = false, has_phi_len_ = false;
    // Every part the index holds is present: built, or loaded whole.
    bool complete_ = true;
    std::optional<SeqIdTable> seq_ids_;
    std::string info_;
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

    /**
     * load() with the parts in *need, or every part the index holds if need
     * is null, when both forms of phi are read.  phi_inv, whose start tables
     * take time to build, is read only with phi.
     */
    void load(std::istream& in, const TmsParts* need) {
        using namespace orbit_idx;
        Reader rd(in);
        const Header& h = rd.header();
        const TmsParts held = held_parts(h);
        const TmsParts want = need ? *need : held;
        check_parts(held, want);
        const uint32_t mask = need ? sections_for(want) : h.held_mask();
        auto wants = [&](Section s) { return (mask & section_bit(s)) != 0; };
        rd.load_move_table(Section::LF, lf_, LF_COLS);
        if (wants(Section::F_STARTS)) {
            const std::vector<uint64_t> v = rd.load_u64s(Section::F_STARTS);
            if (v.size() != F_BLOCKS) throw std::runtime_error("tms index has another number of F block starts");
            for (size_t i = 0; i < F_BLOCKS; ++i) f_rows_[i] = v[i];
        }
        has_psi_ = want.psi;
        has_phi_ = want.phi;
        has_phi_inv_ = want.phi && want.phi_inv;
        has_phi_len_ = want.phi && wants(Section::PHI_LEN);
        has_phi_starts_ = want.phi && wants(Section::PHI_START);
        if (has_phi_len_) rd.load_move_table(Section::PHI_LEN, phi_len_, PHI_LEN_COLS);
        if (has_phi_starts_) rd.load_move_table(Section::PHI_START, phi_, PHI_COLS);
        if (has_phi_inv_) rd.load_move_table(Section::PHI_INV, phi_inv_, PHI_INV_COLS);
        seq_ids_.reset();
        if (const SectionEntry* e = h.find(Section::SEQ_IDS); e && wants(Section::SEQ_IDS)) {
            seq_ids_ = SeqIdTable::from_payload(rd.load_bytes(Section::SEQ_IDS, Kind::SEQ_TABLE), e->rows);
            const std::string bad = seq_ids_mismatch(*seq_ids_, h.n);
            if (!bad.empty()) throw std::runtime_error("tms index's SEQ_IDS section does not fit it: " + bad);
        }
        complete_ = !need || (want.psi == held.psi && want.phi == held.phi && has_phi_inv_ == held.phi_inv &&
                              (!held.phi || (has_phi_len_ && has_phi_starts_)) &&
                              (seq_ids_.has_value() || !h.find(Section::SEQ_IDS)));
        info_ = h.info;
        compute_occurs();
        if (h.bwt_runs != bwt_runs_) throw std::runtime_error("tms index's BWT runs disagree with its header");
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
    template <typename Lens>
    static void union_rows(const std::vector<uchar>& heads, const Lens& lens,
                           const orbit::split_params& sp, std::vector<uchar>& row_heads,
                           orbit::int_vector& row_lens) {
        using Enc = orbit::rlbwt::rlbwt_interval_encoding<ms_alpha::alphabet>;
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
        // The images' starts, in row order, packed.
        orbit::int_vector images(m, orbit::bit_width(split.domain() - 1));
        each([&](ulint k, uchar c) {
            images.set(next_slot[c]++, next_image[c]);
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
        ulint max_len = 0;
        merge([&](uchar, ulint len) {
            ++rows;
            max_len = std::max(max_len, len);
        });
        row_heads.assign(rows, 0);
        row_lens = orbit::int_vector(rows, orbit::bit_width(max_len));
        size_t i = 0;
        merge([&](uchar c, ulint len) {
            row_heads[i] = c;
            row_lens.set(i++, len);
        });
    }

    /**
     * For LF rows with the given heads and lengths (len(i) the length of row
     * i), find psi's pointer and offset for each row's head, calling
     * put(k, pointer, offset) for each row k that psi maps a row head into,
     * and fill the table of F block starts.  Row i with byte c holds the
     * occurrences of c ranked R to R + length - 1 in the BWT, where R counts
     * those in earlier rows, so psi maps F positions C[c] + R onward into it;
     * the rows of c's F block that start there are found by walking the
     * block alongside.
     */
    template <typename Len, typename Put>
    void fill_psi(const std::vector<uchar>& row_heads, Len&& len, Put&& put) {
        const size_t rows = row_heads.size();
        std::array<ulint, 256> count{}, f_pos{}, f_row{};
        for (size_t i = 0; i < rows; ++i) count[row_heads[i]] += len(i);
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
                while (start < f_pos[c]) start += len(i++);
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
            const ulint li = len(i);
            const ulint lo = f_pos[c] + rank[c], hi = lo + li;
            for (; next_start[c] < hi; next_start[c] += len(next_row[c]++)) put(next_row[c], i, next_start[c] - lo);
            rank[c] += li;
        }
    }

    /** The F block of row i: how many F blocks after the first start at or before it. */
    uchar f_block(ulint i) const {
        uchar b = 0;
        for (ulint s : f_rows_) b += i >= s;
        return b;
    }

    /**
     * A starts-based permutation over text positions with m intervals, the
     * x-th starting at start(x) (increasing, from 0) with image image(x),
     * and one integrated column holding an LCP value that drops by one per
     * position inside an interval, sampled at each start as sample(x).
     * Split intervals get the sample minus their distance from the original
     * start.  The column is written in place, at the width of its largest
     * value.
     */
    template <typename Perm, typename Cols, typename Start, typename Image, typename Sample>
    static Perm sampled_permutation(size_t m, Start&& start, Image&& image, Sample&& sample, ulint n,
                                    const orbit::split_params& sp) {
        if (m == 0 || start(0) != 0) throw std::logic_error("a permutation over text positions must start at 0");
        // Views of the intervals' lengths and images, as Orbit's encoding reads them.
        struct Lengths {
            size_t m;
            ulint n;
            const Start* start;
            size_t size() const { return m; }
            ulint operator[](size_t x) const { return (x + 1 < m ? (*start)(x + 1) : n) - (*start)(x); }
        } lengths{m, n, &start};
        struct Images {
            size_t m;
            const Image* image;
            size_t size() const { return m; }
            ulint operator[](size_t x) const { return (*image)(x); }
        } images{m, &image};
        ulint max_length = 0;
        for (size_t x = 0; x < m; ++x) max_length = std::max(max_length, lengths[x]);
        auto enc = orbit::interval_encoding_impl<>::from_lengths_and_images(lengths, images, n, max_length, sp);
        // Call put(i, value) with the column's value for each split interval i.
        auto each = [&](auto put) {
            ulint s = 0;
            size_t x = 0;
            ulint sx = start(0);
            for (ulint i = 0; i < enc.intervals(); ++i) {
                while (x + 1 < m && start(x + 1) <= s) sx = start(++x);
                const ulint into = s - sx, v = sample(x);
                if (v < into) throw std::runtime_error("LCP samples are inconsistent with the RLBWT");
                put(i, v - into);
                s += enc.get_length(i);
            }
        };
        ulint max_value = 0;
        each([&](ulint, ulint v) { max_value = std::max(max_value, v); });
        Perm perm(enc, {orbit::bit_width(max_value)}, Perm::in_place_data);
        each([&](ulint i, ulint v) { perm.template set_data<static_cast<Cols>(0)>(i, v); });
        return perm;
    }

    /**
     * phi's lengths-based form over the intervals of phi, its starts-based
     * form: each row's start becomes its interval's length, and its pointer,
     * offset and PLCP columns are copied.  It reads phi's rows alone, so no
     * list of intervals is held a second time.  Phi points of the two forms
     * then have the same interval and offset.
     */
    static PhiLen lengths_form(const Phi& phi) {
        // PhiLen's row type and fields, which only derived classes reach.
        struct Access : PhiLen {
            using Rows = orbit::packed_vector<typename PhiLen::columns>;
            void assign(Rows&& rows, ulint n, ulint runs, const orbit::split_params& sp) {
                this->move_structure = typename PhiLen::move_structure_perm(std::move(rows), n, runs);
                this->split_params_ = sp;
            }
        };
        using Rows = typename Access::Rows;
        constexpr size_t START = Phi::start_column(), PTR = Phi::pointer_column(), OFF = Phi::offset_column(),
                         PLCP = Phi::template data_column<TmsPhiCols::PLCP>();
        constexpr size_t LEN = PhiLen::length_column(), LPTR = PhiLen::pointer_column(),
                         LOFF = PhiLen::offset_column(), LPLCP = PhiLen::template data_column<TmsPhiCols::PLCP>();
        const ulint m = phi.intervals(), n = phi.domain();
        const auto rd = phi.get_reader();
        auto length = [&](ulint i) { return (i + 1 < m ? rd.template get<START>(i + 1) : n) - rd.template get<START>(i); };
        ulint max_length = 0;
        for (ulint i = 0; i < m; ++i) max_length = std::max(max_length, length(i));
        const auto& pw = phi.get_widths();
        std::array<uchar, 4> w{};
        w[LEN] = orbit::bit_width(max_length);
        w[LPTR] = pw[PTR];
        w[LOFF] = pw[OFF];
        w[LPLCP] = pw[PLCP];
        Rows rows(m, w);
        // Both forms' rows have four columns.
        orbit::packed_matrix<4>& mat = rows;
        for (ulint i = 0; i < m; ++i) {
            std::array<ulint, 4> v{};
            v[LEN] = length(i);
            v[LPTR] = rd.template get<PTR>(i);
            v[LOFF] = rd.template get<OFF>(i);
            v[LPLCP] = rd.template get<PLCP>(i);
            mat.set_row(i, v);
        }
        Access out;
        out.assign(std::move(rows), n, phi.runs(), phi.get_split_params());
        return PhiLen(std::move(out));
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
     *
     * The text positions and tops are kept packed, and the intervals of phi
     * and phi_inv are read through views of them in text order.
     */
    template <typename Lens, typename Enc, typename Tops, typename Phase, typename Release>
    void build_phi(const Lens& lens, const Enc& lf_enc,
                   const Tops& run_tops, const orbit::split_params& sp, bool with_inv,
                   orbit::packed_matrix<1>& phi_cols, Phase&& phase, Release&& release_enc) {
        const ulint lf_count = lf_enc.intervals();
        const ulint n = lf_enc.domain();
        // top holds 1 + the LCP at each LF interval that starts an original
        // run, and 0 at the others.
        orbit::int_vector sa_head(lf_count, orbit::bit_width(n - 1)), sa_tail(lf_count, orbit::bit_width(n - 1)), top;
        {
            ulint max_top = 0;
            for (ulint t : run_tops) max_top = std::max(max_top, t);
            top = orbit::int_vector(lf_count, orbit::bit_width(max_top + 1));
        }
        {
            // The LCP at each LF interval that starts an original run.
            ulint k = 0;
            for (size_t r = 0; r < lens.size(); ++r) {
                top.set(k, run_tops[r] + 1);
                ulint covered = 0;
                while (covered < lens[r]) covered += lf_enc.get_length(k++);
            }
            orbit::rlbwt::lf_move<false, ms_alpha::alphabet> lf(lf_enc);
            auto pos = lf.first();
            ulint sa = n - 1;
            std::vector<bool> seen(lf_count, false);
            for (ulint i = 0; i < n; ++i, --sa) {
                if (pos.offset == 0) { sa_head[pos.interval] = sa; seen[pos.interval] = true; }
                if (pos.offset + 1 == lf.get_length(pos.interval)) sa_tail[pos.interval] = sa;
                pos = lf.LF(pos);
            }
            if (pos != lf.first() || std::find(seen.begin(), seen.end(), false) != seen.end())
                throw std::runtime_error("LF is not a single cycle (a string collection, such as an FMD, with one "
                                         "terminator per string); building phi by an LF walk needs one terminator, "
                                         "so build it from TeraLCP -ophi samples (tms-build --phi) instead");
        }
        // Which LF intervals' bytes differ from the one before's (k = 0
        // counts as differing), after which release_enc may free lf_enc.
        std::vector<bool> differs(lf_count);
        {
            const auto& lf_heads = lf_enc.get_heads();
            for (ulint k = 0; k < lf_count; ++k) differs[k] = k == 0 || lf_heads[k] != lf_heads[k - 1];
        }
        release_enc();
        phase("tms_phi_walk");
        auto top_of = [&](ulint k) {
            const ulint t = top[k];
            if (t == 0) throw std::logic_error("a true run head must start an original run");
            return t - 1;
        };
        // The LF intervals that give the intervals, in text order.
        std::vector<ulint> order;
        // phi from true run heads.
        for (ulint k = 0; k < lf_count; ++k)
            if (differs[k]) order.push_back(k);
        std::sort(order.begin(), order.end(), [&](ulint a, ulint b) { return sa_head[a] < sa_head[b]; });
        phi_ = sampled_permutation<Phi, TmsPhiCols>(
            order.size(), [&](size_t x) -> ulint { return sa_head[order[x]]; },
            [&](size_t x) -> ulint { const ulint k = order[x]; return sa_tail[k == 0 ? lf_count - 1 : k - 1]; },
            [&](size_t x) { return top_of(order[x]); }, n, sp);
        std::vector<ulint>().swap(order);
        phase("tms_phi");
        // The phi interval that starts at each run head's text position,
        // located with a start table (see the file comment); rows that start
        // no run get 0.
        phi_table_ = StartTable(phi_, start_table_shift(phi_.domain(), phi_.intervals()));
        phi_cols = orbit::packed_matrix<1>(lf_count, {orbit::bit_width(phi_.intervals() - 1)});
        {
            const PhiWalker w = phi_walker();
            for (ulint k = 0; k < lf_count; ++k) {
                if (!differs[k]) {
                    phi_cols.set<0>(k, 0);
                    continue;
                }
                phi_cols.set<0>(k, run_head_phi_interval(w, sa_head[k]));
            }
        }
        phi_table_ = StartTable();
        phase("tms_phi_cols");
        if (!with_inv) return;
        // phi_inv from true run tails.
        for (ulint k = 0; k < lf_count; ++k)
            if (k + 1 == lf_count || differs[k + 1]) order.push_back(k);
        std::sort(order.begin(), order.end(), [&](ulint a, ulint b) { return sa_tail[a] < sa_tail[b]; });
        // phi_inv's intervals, packed, so that the walk's data is freed before
        // phi_inv's encoding is made.
        const size_t m = order.size();
        orbit::int_vector starts(m, orbit::bit_width(n - 1)), images(m, orbit::bit_width(n - 1)),
            samples(m, top.get_width());
        for (size_t x = 0; x < m; ++x) {
            const ulint k = order[x];
            starts.set(x, sa_tail[k]);
            images.set(x, sa_head[k + 1 == lf_count ? 0 : k + 1]);
            samples.set(x, k + 1 == lf_count ? ulint(0) : top_of(k + 1));
        }
        std::vector<ulint>().swap(order);
        sa_head = orbit::int_vector();
        sa_tail = orbit::int_vector();
        top = orbit::int_vector();
        phi_inv_ = sampled_permutation<PhiInv, TmsPhiInvCols>(
            m, [&](size_t x) -> ulint { return starts[x]; }, [&](size_t x) -> ulint { return images[x]; },
            [&](size_t x) -> ulint { return samples[x]; }, n, sp);
    }

    /**
     * Build phi (and phi_inv if with_inv) from samples (see TmsPhiSamples)
     * and fill the PHI columns, without walking LF; samples is freed.
     *
     * The phi interval of record i starts at head[i] and maps to tail[i - 1]
     * (tail[R - 1] for i = 0, so SA[0] maps to SA[n - 1]), with PLCP sample
     * top[i].  The phi_inv interval of record i starts at tail[i] and maps to
     * head[i + 1] (head[0] for i = R - 1), with PLCPB sample top[i + 1] (0 for
     * i = R - 1).  Since each terminator row is a record, no interval holds a
     * text position whose predecessor is a terminator, which is where the LCP
     * of the next position need not be one less.
     *
     * The PHI columns of the first LF interval of each BWT run hold the phi
     * point of the run's head row, which is where queries read them; the PHI
     * columns of other LF intervals hold 0.  They are filled last, when only
     * the heads of the samples are still held.  LF's encoding is needed only
     * to find the first LF interval of each run, after which release_enc may
     * free it.
     */
    template <typename Lens, typename Enc, typename Phase, typename Release>
    void build_phi_from_samples(const std::vector<uchar>& heads, const Lens& lens, const Enc& lf_enc,
                                TmsPhiSamples& samples, const orbit::split_params& sp, bool with_inv,
                                orbit::packed_matrix<1>& phi_cols, Phase&& phase, Release&& release_enc) {
        const ulint R = samples.size(), n = lf_enc.domain(), lf_count = lf_enc.intervals();
        if (samples.tail.size() != R || samples.top.size() != R)
            throw std::invalid_argument("phi samples must have as many tails and tops as heads");
        const ulint want = TmsPhiSamples::records_for(heads, lens);
        if (R != want)
            throw std::invalid_argument("phi samples have " + std::to_string(R) + " records, but the RLBWT needs " +
                                        std::to_string(want) + " (one per run, and one per row of a terminator run)");
        const auto& head = samples.head;
        const auto& tail = samples.tail;
        const auto& top = samples.top;
        for (ulint i = 0; i < R; ++i)
            if (head[i] >= n || tail[i] >= n || top[i] >= n)
                throw std::runtime_error("phi sample " + std::to_string(i) + " lies outside the text");
        // The first LF interval of each run.
        std::vector<bool> run_start(lf_count, false);
        {
            ulint k = 0;
            for (size_t r = 0; r < heads.size(); ++r) {
                run_start[k] = true;
                for (ulint covered = 0; covered < ulint(lens[r]);) covered += lf_enc.get_length(k++);
            }
        }
        release_enc();
        // The records in the text order of key, packed.
        auto text_order = [&](const orbit::int_vector& key) {
            std::vector<ulint> order(R);
            std::iota(order.begin(), order.end(), ulint(0));
            std::sort(order.begin(), order.end(), [&](ulint a, ulint b) { return key[a] < key[b]; });
            if (key[order[0]] != 0) throw std::runtime_error("phi samples lack text position 0");
            for (ulint x = 1; x < R; ++x)
                if (key[order[x]] == key[order[x - 1]]) throw std::runtime_error("phi samples repeat a text position");
            orbit::int_vector packed(R, orbit::bit_width(R - 1));
            for (ulint x = 0; x < R; ++x) packed.set(x, order[x]);
            return packed;
        };
        // A structure's interval starts, images and LCP samples in text
        // order, packed, so that building it reads them in order.
        orbit::int_vector starts, images, lcps;
        auto lay_out = [&](const orbit::int_vector& order, auto start, auto image, auto sample) {
            starts = orbit::int_vector(R, orbit::bit_width(n - 1));
            images = orbit::int_vector(R, orbit::bit_width(n - 1));
            lcps = orbit::int_vector(R, top.get_width());
            for (ulint x = 0; x < R; ++x) {
                const ulint i = order[x];
                starts.set(x, start(i));
                images.set(x, image(i));
                lcps.set(x, sample(i));
            }
        };
        auto build_from_lists = [&](auto& perm, auto cols_tag) {
            using P = std::remove_reference_t<decltype(perm)>;
            perm = sampled_permutation<P, decltype(cols_tag)>(
                R, [&](size_t x) -> ulint { return starts[x]; }, [&](size_t x) -> ulint { return images[x]; },
                [&](size_t x) -> ulint { return lcps[x]; }, n, sp);
            starts = images = lcps = orbit::int_vector();
        };
        {
            orbit::int_vector order = text_order(head);
            lay_out(order, [&](ulint i) -> ulint { return head[i]; },
                    [&](ulint i) -> ulint { return tail[i == 0 ? R - 1 : i - 1]; },
                    [&](ulint i) -> ulint { return top[i]; });
        }
        phase("tms_phi_order");
        build_from_lists(phi_, TmsPhiCols{});
        phase("tms_phi");
        if (with_inv) {
            {
                orbit::int_vector order = text_order(tail);
                lay_out(order, [&](ulint i) -> ulint { return tail[i]; },
                        [&](ulint i) -> ulint { return head[i + 1 == R ? 0 : i + 1]; },
                        [&](ulint i) -> ulint { return i + 1 == R ? ulint(0) : ulint(top[i + 1]); });
            }
            // Only the heads are needed after this.
            samples.tail = orbit::int_vector();
            samples.top = orbit::int_vector();
            phase("tms_phi_inv_order");
            build_from_lists(phi_inv_, TmsPhiInvCols{});
            phase("tms_phi_inv");
        }
        samples.tail = orbit::int_vector();
        samples.top = orbit::int_vector();
        // The phi interval that starts at each run head's text position,
        // located with a start table (see the file comment); rows that start
        // no run get 0.
        phi_table_ = StartTable(phi_, start_table_shift(phi_.domain(), phi_.intervals()));
        phi_cols = orbit::packed_matrix<1>(lf_count, {orbit::bit_width(phi_.intervals() - 1)});
        {
            const PhiWalker w = phi_walker();
            ulint r = 0, x = 0;
            for (ulint k = 0; k < lf_count; ++k) {
                if (!run_start[k]) {
                    phi_cols.set<0>(k, 0);
                    continue;
                }
                phi_cols.set<0>(k, run_head_phi_interval(w, head[x]));
                x += heads[r] == orbit::TERMINATOR ? ulint(lens[r]) : 1;
                ++r;
            }
        }
        phi_table_ = StartTable();
        std::vector<bool>().swap(run_start);
        phase("tms_phi_cols");
        samples = TmsPhiSamples();
    }

    static constexpr ulint NOT_A_RUN_HEAD = std::numeric_limits<ulint>::max();

    /**
     * The phi interval that starts at text position x, a run head's.  Throws
     * if none starts there: LF's phi columns store no offset (see the file
     * comment).
     */
    static ulint run_head_phi_interval(const PhiWalker& w, ulint x) {
        const PhiPos q = w.locate(x);
        if (q.offset != 0) throw std::logic_error("no phi interval starts at a run head's text position");
        return q.interval;
    }
};

#endif /* _TMS_INDEX_HPP */
