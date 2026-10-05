/**
 * ms split-report: the splitting an index's move structures were built
 * with and the balance they reach.  This is a development tool; it reads an
 * index and changes nothing.
 *
 * For each move structure an MSIX (Orbit-ms) or TMSX (Orbit-tms) file
 * holds, it prints:
 *
 *  - the split_params stored with it.  A structure built with NO_SPLITTING
 *    on rows that were already split (Orbit-ms's LF, and Orbit-tms's LF and
 *    psi in the psi and full layouts) stores NO_SPLITTING, so its stored
 *    parameters do not say how its rows were cut; the report says so.
 *  - its intervals, its longest interval, n, the runs it stores, the BWT
 *    runs, and the length cap that Orbit's length capping (factor 8 unless
 *    stored otherwise) implies for n over the runs that were split.  Orbit
 *    computes the cap as the largest value of bit_width(ceil(n / r * L))
 *    bits (move_splitting::split_by_length_capping), so it is one less than
 *    a power of two.
 *  - its static balance: for each output interval (the image of an input
 *    interval), how many input interval starts lie strictly inside it.  A
 *    move lands in an output interval and fast-forwards across exactly
 *    those starts at most, so the maximum bounds every fast-forward.
 *    Orbit's balancing with factor b (move_splitting::balance) splits an
 *    interval of either side once 2b starts of the other side lie strictly
 *    inside it, so a balanced structure has at most 2b - 1 in every
 *    interval.  The report gives the maximum, a histogram, and whether the
 *    maximum is within 2b - 1 for the default b of 16 (or the stored b).
 *
 * The counts take one pass over the rows per structure, shared among
 * threads.  TMSX parts are loaded one group at a time (LF with psi, then
 * phi's lengths-based form, then phi's starts-based form with phi_inv), as
 * tms-batch loads them, so that less of a large index is held at once.
 */

#ifndef _SPLIT_REPORT_HPP
#define _SPLIT_REPORT_HPP

#include "ms_rlbwt.hpp"
#include "serialize.hpp"
#include "tms_index.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

namespace split_report {

/** Counts of input starts strictly inside each output interval, and interval lengths. */
struct Balance {
    static constexpr size_t B = 256;  // counts 0 to B - 1 apart, larger ones in bucket B
    std::array<ulint, B + 1> by_count{};
    ulint intervals = 0, total = 0, max = 0, max_at = 0, max_len = 0;
    bool bad = false;  // a pointer or offset that does not describe a valid image

    void add(ulint i, ulint c, ulint len) {
        ++intervals;
        total += c;
        if (c > max || intervals == 1) { max = c; max_at = i; }
        ++by_count[std::min<ulint>(c, B)];
        max_len = std::max(max_len, len);
    }
    Balance& operator+=(const Balance& o) {
        for (size_t k = 0; k <= B; ++k) by_count[k] += o.by_count[k];
        if (o.intervals && (o.max > max || intervals == 0)) { max = o.max; max_at = o.max_at; }
        intervals += o.intervals;
        total += o.total;
        max_len = std::max(max_len, o.max_len);
        bad = bad || o.bad;
        return *this;
    }
    ulint quantile(double q) const {
        const double want = q * double(intervals);
        ulint acc = 0;
        for (size_t k = 0; k <= B; ++k) {
            acc += by_count[k];
            if (double(acc) >= want) return k;
        }
        return B;
    }
};

/** Run f(lo, hi, part) over [0, m) cut into one range per thread, and add up the parts. */
template <typename F>
inline Balance in_parallel(ulint m, size_t threads, F&& f) {
    threads = std::max<size_t>(1, std::min<size_t>(threads, std::max<ulint>(1, m / 4096)));
    std::vector<Balance> parts(threads);
    std::vector<std::thread> pool;
    for (size_t t = 0; t < threads; ++t) {
        const ulint lo = m * t / threads, hi = m * (t + 1) / threads;
        pool.emplace_back([&, lo, hi, t] { f(lo, hi, parts[t]); });
    }
    for (auto& th : pool) th.join();
    Balance all;
    for (const auto& p : parts) all += p;
    return all;
}

/**
 * Balance of a lengths-based structure with m intervals: interval i has
 * length len(i) and image starting at offset off(i) into interval ptr(i).
 * The input starts strictly inside the image are those of intervals
 * ptr(i) + 1 onward that start before the image ends.
 */
template <typename Len, typename Ptr, typename Off>
inline Balance balance_lengths(ulint m, Len&& len, Ptr&& ptr, Off&& off, size_t threads) {
    return in_parallel(m, threads, [&](ulint lo, ulint hi, Balance& b) {
        for (ulint i = lo; i < hi; ++i) {
            const ulint L = len(i), p = ptr(i), o = off(i);
            if (p >= m || o >= len(p)) { b.bad = true; continue; }
            // Distance from the start of interval p to the image's end, and
            // to the start of the next interval after p + c.
            const ulint end = o + L;
            ulint at = len(p), c = 0;
            while (at < end && p + c + 1 < m) {
                ++c;
                at += len(p + c);
            }
            if (at < end) b.bad = true;
            b.add(i, c, L);
        }
    });
}

/**
 * Balance of a starts-based structure with m intervals over n positions:
 * interval i starts at start(i) and its image starts at start(ptr(i)) +
 * off(i).
 */
template <typename Start, typename Ptr, typename Off>
inline Balance balance_starts(ulint m, ulint n, Start&& start, Ptr&& ptr, Off&& off, size_t threads) {
    return in_parallel(m, threads, [&](ulint lo, ulint hi, Balance& b) {
        for (ulint i = lo; i < hi; ++i) {
            const ulint L = (i + 1 < m ? start(i + 1) : n) - start(i), p = ptr(i);
            if (p >= m) { b.bad = true; continue; }
            const ulint img = start(p) + off(i), end = img + L;
            if (end > n || img >= (p + 1 < m ? start(p + 1) : n)) { b.bad = true; continue; }
            ulint c = 0;
            while (p + c + 1 < m && start(p + c + 1) < end) ++c;
            b.add(i, c, L);
        }
    });
}

/** split_params as text. */
inline std::string params_text(const orbit::split_params& sp) {
    std::ostringstream os;
    os << "length_capping=";
    if (sp.length_capping) os << *sp.length_capping;
    else os << "none";
    os << " balancing=";
    if (sp.balancing) os << *sp.balancing;
    else os << "none";
    return os.str();
}

/** The length cap Orbit's length capping with factor L gives n positions in r runs. */
inline ulint implied_cap(ulint n, ulint r, double L) {
    const double avg = double(n) / double(r);
    const ulint want = static_cast<ulint>(std::ceil(avg * L));
    return orbit::max_val(orbit::bit_width(want));
}

/** What the report knows about one structure besides its balance. */
struct Info {
    std::string name;            // lf, psi, phi_len, phi_starts, phi_inv
    orbit::split_params stored;  // as stored with the structure
    std::string note;            // how its rows were cut, when the stored parameters do not say
    ulint n = 0, intervals = 0, runs_stored = 0, runs_split = 0;
    std::string runs_split_what;  // what runs_split counts
};

/**
 * Print one structure's report: a line of key=value fields, then its
 * histogram as count:intervals pairs.  The cap and the balance guarantee
 * use the stored parameters when they split, else Orbit's defaults.
 */
inline void print(std::ostream& os, const Info& f, const Balance& b, double seconds) {
    const bool stored_split = f.stored != orbit::NO_SPLITTING;
    const double L = stored_split && f.stored.length_capping ? *f.stored.length_capping : *orbit::DEFAULT_LENGTH_CAPPING;
    const ulint bf = stored_split && f.stored.balancing ? *f.stored.balancing : *orbit::DEFAULT_BALANCING;
    const ulint cap = implied_cap(f.n, f.runs_split, L);
    const ulint bound = 2 * bf - 1;
    os << "split: structure=" << f.name << " stored_params=";
    if (stored_split) os << "{" << params_text(f.stored) << "}";
    else os << "NO_SPLITTING";
    os << " intervals=" << f.intervals << " n=" << f.n << " runs_stored=" << f.runs_stored
       << " runs_split=" << f.runs_split << " (" << f.runs_split_what << ")"
       << " n_over_r=" << double(f.n) / double(f.runs_split) << " cap_factor=" << L << " implied_cap=" << cap
       << " max_len=" << b.max_len << " within_cap=" << (b.max_len <= cap ? "yes" : "no") << "\n";
    if (!f.note.empty()) os << "split: structure=" << f.name << " note: " << f.note << "\n";
    const double m = double(b.intervals);
    os << "balance: structure=" << f.name << " intervals=" << b.intervals << " mean=" << double(b.total) / m
       << " p50=" << b.quantile(0.5) << " p99=" << b.quantile(0.99) << " p9999=" << b.quantile(0.9999)
       << " max=" << b.max << " max_at=" << b.max_at << " balancing_factor=" << bf << " guarantee_max=" << bound
       << " over=" << [&] { ulint s = 0; for (size_t k = bound + 1; k <= Balance::B; ++k) s += b.by_count[k]; return s; }()
       << " result=" << (b.bad ? "ERROR(invalid pointer or offset)" : b.max <= bound ? "PASS" : "FAIL")
       << " seconds=" << seconds << "\n";
    os << "balance-hist: structure=" << f.name;
    for (size_t k = 0; k <= Balance::B; ++k)
        if (b.by_count[k]) os << ' ' << k << (k == Balance::B ? "+" : "") << ':' << b.by_count[k];
    os << "\n";
    os.flush();
}

template <typename F>
inline std::pair<Balance, double> timed(F&& f) {
    const auto t0 = std::chrono::steady_clock::now();
    Balance b = f();
    return {b, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count()};
}

/** Maximal runs of equal values among chr(0), ..., chr(m - 1). */
template <typename Chr>
inline ulint count_runs(ulint m, Chr&& chr) {
    ulint runs = 0;
    for (ulint i = 0, prev = ~ulint(0); i < m; ++i) {
        const ulint c = chr(i);
        runs += c != prev;
        prev = c;
    }
    return runs;
}

/** Report an MSIX (Orbit-ms) index: LF only. */
inline int report_ms(const std::string& path, size_t threads, std::ostream& os) {
    auto t0 = std::chrono::steady_clock::now();
    auto opt = ms_serialize::read_index(path);
    if (!opt) { std::cerr << "Failed to load ms index: " << path << "\n"; return 1; }
    os << "index: " << path << " format=MSIX load_s="
       << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << "\n";
    const auto& lf = opt->impl();
    using Impl = std::remove_cv_t<std::remove_reference_t<decltype(lf)>>;
    constexpr size_t LEN = Impl::length_column(), PTR = Impl::pointer_column(), OFF = Impl::offset_column(),
                     CHR = Impl::character_column();
    const auto rd = lf.get_reader();
    const ulint m = lf.intervals();
    Info f;
    f.name = "lf";
    f.stored = lf.get_split_params();
    f.n = lf.domain();
    f.intervals = m;
    f.runs_stored = lf.runs();
    f.runs_split = count_runs(m, [&](ulint i) { return rd.template get<CHR>(i); });
    f.runs_split_what = "BWT runs, maximal runs of one character in the rows; LF splitting ran on these runs as "
                        "given to the build, which LCP splitting may have cut further";
    if (f.stored == orbit::NO_SPLITTING)
        f.note = "the format stores no LF split parameters: LF is built with NO_SPLITTING on rows that "
                 "apply_lf_splitting already cut with the build's --lf-split (default length capping 8, "
                 "balancing 16), so the stored params are NO_SPLITTING and runs_stored is the row count";
    auto [b, s] = timed([&] {
        return balance_lengths(
            m, [&](ulint i) { return rd.template get<LEN>(i); }, [&](ulint i) { return rd.template get<PTR>(i); },
            [&](ulint i) { return rd.template get<OFF>(i); }, threads);
    });
    print(os, f, b, s);
    return 0;
}

/** Load the parts in need of the tms index at path. */
inline std::optional<TmsIndex> load_tms(const std::string& path, const TmsParts& need, std::ostream& os,
                                        const char* what) {
    auto t0 = std::chrono::steady_clock::now();
    std::ifstream in(path, std::ios::binary);
    TmsIndex idx;
    try {
        idx.load(in, need);
    } catch (const std::exception& e) {
        std::cerr << "Failed to load tms index " << path << ": " << e.what() << "\n";
        return std::nullopt;
    }
    os << "load: " << what << " load_s=" << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count()
       << "\n";
    return idx;
}

/** Report a TMSX (Orbit-tms) index, one group of parts at a time. */
inline int report_tms(const std::string& path, size_t threads, std::ostream& os) {
    TmsParts held;
    try {
        held = TmsIndex::held_parts(path);
    } catch (const std::exception& e) {
        std::cerr << path << ": " << e.what() << "\n";
        return 1;
    }
    os << "index: " << path << " format=TMSX layout=" << TmsIndex::layout_name(held)
       << " phi_inv=" << (held.phi_inv ? "yes" : "no") << "\n";
    {
        // LF, with the PSI columns in the psi and full layouts.
        TmsParts need;
        need.psi = held.psi;
        auto opt = load_tms(path, need, os, "lf");
        if (!opt) return 1;
        const auto& lf = opt->lf();
        using LF = TmsIndex::LF;
        constexpr size_t LEN = LF::length_column(), PTR = LF::pointer_column(), OFF = LF::offset_column();
        constexpr size_t PSI_INT = LF::template data_column<TmsLFCols::PSI_INT>(),
                         PSI_OFF = LF::template data_column<TmsLFCols::PSI_OFF>();
        const auto rd = lf.get_reader();
        const ulint m = lf.intervals();
        const ulint bwt_runs = opt->bwt_runs();
        Info f;
        f.name = "lf";
        f.stored = lf.get_split_params();
        f.n = lf.domain();
        f.intervals = m;
        f.runs_stored = lf.runs();
        f.runs_split = held.psi ? bwt_runs : lf.runs();
        f.runs_split_what = held.psi ? "BWT runs, maximal runs of one character in the rows"
                                     : "runs stored with LF, those the build split";
        if (f.stored == orbit::NO_SPLITTING)
            f.note = held.psi ? "the format stores no LF split parameters: the rows are the union of LF's "
                                "intervals split with tms-build's --lf-split (default length capping 8, "
                                "balancing 16) and their images, built with NO_SPLITTING and not rebalanced, so "
                                "the stored params are NO_SPLITTING and runs_stored is the row count"
                              : "stored as NO_SPLITTING";
        auto [b, s] = timed([&] {
            return balance_lengths(
                m, [&](ulint i) { return rd.template get<LEN>(i); }, [&](ulint i) { return rd.template get<PTR>(i); },
                [&](ulint i) { return rd.template get<OFF>(i); }, threads);
        });
        print(os, f, b, s);
        if (held.psi) {
            // psi reads the same rows, in F order, with its own pointer and offset.
            Info g = f;
            g.name = "psi";
            g.note = "psi shares LF's rows and stored params; its pointer and offset are the PSI columns";
            auto [bp, sp] = timed([&] {
                return balance_lengths(
                    m, [&](ulint i) { return rd.template get<LEN>(i); },
                    [&](ulint i) { return rd.template get<PSI_INT>(i); },
                    [&](ulint i) { return rd.template get<PSI_OFF>(i); }, threads);
            });
            print(os, g, bp, sp);
        }
    }
    if (!held.phi) return 0;
    {
        // phi's lengths-based form.
        TmsParts need;
        need.phi = true;
        auto opt = load_tms(path, need, os, "phi_len");
        if (!opt) return 1;
        const auto& phi = opt->phi_len_perm();
        using P = TmsIndex::PhiLen;
        constexpr size_t LEN = P::length_column(), PTR = P::pointer_column(), OFF = P::offset_column();
        const auto rd = phi.get_reader();
        Info f;
        f.name = "phi_len";
        f.stored = phi.get_split_params();
        f.n = phi.domain();
        f.intervals = phi.intervals();
        f.runs_stored = phi.runs();
        f.runs_split = phi.runs();
        f.runs_split_what = "runs stored with phi: its sample records, those the build split";
        f.note = "phi's lengths-based form copies the intervals and params of the starts-based form";
        auto [b, s] = timed([&] {
            return balance_lengths(
                phi.intervals(), [&](ulint i) { return rd.template get<LEN>(i); },
                [&](ulint i) { return rd.template get<PTR>(i); }, [&](ulint i) { return rd.template get<OFF>(i); },
                threads);
        });
        print(os, f, b, s);
    }
    {
        // phi's starts-based form, and phi_inv.
        TmsParts need;
        need.phi = need.phi_starts = true;
        need.phi_inv = held.phi_inv;
        auto opt = load_tms(path, need, os, held.phi_inv ? "phi_starts+phi_inv" : "phi_starts");
        if (!opt) return 1;
        auto one = [&](const auto& perm, const char* name, const char* runs_what) {
            using P = std::remove_cv_t<std::remove_reference_t<decltype(perm)>>;
            constexpr size_t START = P::start_column(), PTR = P::pointer_column(), OFF = P::offset_column();
            const auto rd = perm.get_reader();
            Info f;
            f.name = name;
            f.stored = perm.get_split_params();
            f.n = perm.domain();
            f.intervals = perm.intervals();
            f.runs_stored = perm.runs();
            f.runs_split = perm.runs();
            f.runs_split_what = runs_what;
            auto [b, s] = timed([&] {
                return balance_starts(
                    perm.intervals(), perm.domain(), [&](ulint i) { return rd.template get<START>(i); },
                    [&](ulint i) { return rd.template get<PTR>(i); }, [&](ulint i) { return rd.template get<OFF>(i); },
                    threads);
            });
            print(os, f, b, s);
        };
        one(opt->phi_perm(), "phi_starts", "runs stored with phi: its sample records, those the build split");
        if (held.phi_inv)
            one(opt->phi_inv_perm(), "phi_inv",
                "runs stored with phi_inv: its sample records, split on their own with phi's params");
    }
    return 0;
}

/** Report the index at path, an MSIX or a TMSX file. */
inline int run(const std::string& path, size_t threads, std::ostream& os) {
    std::ifstream in(path, std::ios::binary);
    char magic[4] = {};
    in.read(magic, 4);
    if (!in.good()) { std::cerr << "Cannot read " << path << "\n"; return 1; }
    if (std::memcmp(magic, "MSIX", 4) == 0) return report_ms(path, threads, os);
    if (std::memcmp(magic, "TMSX", 4) == 0) return report_tms(path, threads, os);
    std::cerr << path << " is neither an MSIX nor a TMSX index\n";
    return 1;
}

}  // namespace split_report

#endif /* _SPLIT_REPORT_HPP */
