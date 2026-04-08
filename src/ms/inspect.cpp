/**
 * Bit-usage inspection and spillover TSV dump for MS index.
 *
 * Author: Ben Langmead (ben.langmead@gmail.com)
 * Date: Feb 17, 2026
 */

#include "inspect.hpp"
#include "ms_rlbwt.hpp"
#include <iostream>
#include <fstream>

namespace ms_inspect {

/**
 * Compute bit usage for different components of the index and print a brief
 * report to stdout.
 */
void run_inspect(MSIndexSpillLCP<false>& index) {
    const ulint r = index.move_runs(), n = index.domain();
    std::cout << "MS Index Inspection\n";
    std::cout << "  r (runs): " << r << ", n (input len): " << n << "\n";
    if (index.spill_align() > 1)
        std::cout << "  spill_align: " << index.spill_align() << " bytes\n";
    if (index.spill_split_bits() > 0)
        std::cout << "  spill_split_bits: " << static_cast<unsigned>(index.spill_split_bits()) << " (" << (1u << index.spill_split_bits()) << " arrays)\n";

    const auto& widths = index.get_widths();
    const uchar w_length = widths[0];
    const uchar w_pointer = widths[1];
    const uchar w_offset = widths[2];
    const uchar w_char = widths[3];
    const uchar w_top = widths[4];
    const uchar w_sub = widths[5];
    const uchar w_spill = widths[6];

    size_t base_bits = r * (w_length + w_pointer + w_offset + w_char);
    size_t lcp_bits = r * (w_top + w_sub + w_spill);
    size_t spill_bytes = index.spillover_total_bytes();
    size_t spill_bits = spill_bytes * 8;

    std::cout << "  base RLBWT bits: " << base_bits << "\n";
    std::cout << "    length: " << static_cast<unsigned>(w_length) << " bits x " << r << " runs\n";
    std::cout << "    pointer: " << static_cast<unsigned>(w_pointer) << " bits x " << r << " runs\n";
    std::cout << "    offset: " << static_cast<unsigned>(w_offset) << " bits x " << r << " runs\n";
    std::cout << "    character: " << static_cast<unsigned>(w_char) << " bits x " << r << " runs\n";
    std::cout << "  LCP columns bits: " << lcp_bits << "\n";
    std::cout << "    lcp_top: " << static_cast<unsigned>(w_top) << " bits x " << r << " runs\n";
    std::cout << "    lcp_min_sub: " << static_cast<unsigned>(w_sub) << " bits x " << r << " runs\n";
    std::cout << "    lcp_spill: " << static_cast<unsigned>(w_spill) << " bits x " << r << " runs (offset into spillover)\n";
    std::cout << "    spillover payload: variable (ULEB128), total in spill vector below\n";
    std::cout << "  spillover vector bits: " << spill_bits << "\n";
    std::cout << "  fraction for LCPs: "
        << (double)(lcp_bits + spill_bits) / (base_bits + lcp_bits + spill_bits) * 100.0
        << "%\n";

    // Now compute the bytes per row for the main data structure (7 fields)
    // Fields: length, pointer, offset, character, lcp_top, lcp_min_sub, lcp_spill (all as bits per row)
    size_t bits_per_row = static_cast<size_t>(w_length)
                        + static_cast<size_t>(w_pointer)
                        + static_cast<size_t>(w_offset)
                        + static_cast<size_t>(w_char)
                        + static_cast<size_t>(w_top)
                        + static_cast<size_t>(w_sub)
                        + static_cast<size_t>(w_spill);
    size_t bytes_per_row = (bits_per_row + 7) / 8; // round up to next byte

    std::cout << "  bytes per row (7 main-structure fields, not counting spillover): "
              << bytes_per_row << std::endl;
}

/**
 * Dump spillover table to TSV with the following columns:
 *   - run:          The run index (0-based)
 *   - jumbo:        Whether the run is a "jumbo" run (T/F); a jumbo run is one
 *                   whose LCP_TOP and LCP_MIN_SUB values are both at their
 *                   respective maximums
 *   - coalesced:    T if the spillover table contains any entries for this run
 *                   (i.e., at least 1 interior value), F otherwise. This means
 *                   the LCPs were too large to fit directly and were
 *                   "coalesced" (stored out of line)
 *   - num_interior: The number of interior (spillover) LCP values for the run
 *   - rest:         A comma-separated list encoding the spillover values as
 *                   offset,value pairs (offset0,value0,offset1,value1,...);
 *                   set to "-" if num_interior==0
 */
bool run_spillover_tsv(MSIndexSpillLCP<false>& index, const std::string& path) {
    std::ofstream out(path);
    if (!out.good()) return false;

    out << "run\tjumbo\tcoalesced\tnum_interior\trest\n";
    const ulint r = index.move_runs();
    const ulint max_top = index.max_lcp_top();
    const ulint max_sub = index.max_lcp_min_sub();

    for (ulint i = 0; i < r; ++i) {
        ulint top = index.template get<LCPSpillRunCols::LCP_TOP>(i);
        ulint sub = index.template get<LCPSpillRunCols::LCP_MIN_SUB>(i);
        ulint so = index.template get<LCPSpillRunCols::LCP_SPILL>(i);
        bool jumbo = (top == max_top && sub == max_sub);
        if (so == NO_SPILL) continue;

        out << i << "\t" << (jumbo ? "T" : "F") << "\t";

        const auto& spill = index.spillover_for_row(i);
        size_t p = index.spill_offset_bytes(so);
        if (jumbo) {
            p = skip_uleb128(spill, p);
            p = skip_uleb128(spill, p);
        }
        auto [n_int, np] = decode_uleb128(spill, p);
        p = np;

        out << (n_int > 0 ? "T" : "F") << "\t" << n_int << "\t";

        if (n_int == 0) {
            out << "-\n";
            continue;
        }

        std::string rest;
        for (ulint k = 0; k < n_int; ++k) {
            auto [o, no] = decode_uleb128(spill, p);
            auto [vv, nv] = decode_uleb128(spill, no);
            p = nv;
            if (k > 0) rest += ",";
            rest += std::to_string(o) + "," + std::to_string(vv);
        }
        out << rest << "\n";
    }
    return out.good();
}

}  // namespace ms_inspect
