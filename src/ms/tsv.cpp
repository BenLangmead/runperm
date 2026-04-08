/**
 * TSV loading and parsing implementations.  Right now, this follows the TSV
 * format from the movify.py script.  It is not an "official" format in any way.
 *
 * Author: Ben Langmead (ben.langmead@gmail.com)
 * Date: Feb 17, 2026
 */

#include "tsv.hpp"
#include "orbit/common.hpp"
#include "ms_constants.hpp"
#include <fstream>
#include <sstream>

namespace tsv {

/** Map movify '#' to SEP (1), keep 'A','C','G','T' as-is */
uchar map_movify_char(char c) {
    if (c == '#') return 1;
    return static_cast<uchar>(c);
}

/**
 * Parse a single LCP row: "0,0,2465,..." or "0,-,2465" (gap as -) ->
 * vector<ulint>
 */
std::vector<ulint> parse_lcp_row(const std::string& s) {
    std::vector<ulint> out;
    if (s.empty()) return out;
    std::istringstream iss(s);
    std::string tok;
    while (std::getline(iss, tok, ',')) {
        if (tok.empty() || tok == "-") {
            out.push_back(LCP_GAP);
        } else {
            out.push_back(static_cast<ulint>(std::stoull(tok)));
        }
    }
    return out;
}

/**
 * Load a TSV file and populate the given vectors.
 * TSV should have the following columns: id, off, len, c, sa, lcp.
 * The id and off columns are ignored as we can recalculate them.
 * The sa column is also ignored for now.
 */
bool load_tsv(
    const std::string& path,
    std::vector<uchar>& bwt_heads,
    std::vector<ulint>& bwt_run_lengths,
    std::vector<std::vector<ulint>>& lcps_per_run)
{
    std::ifstream f(path);
    if (!f.good()) return false;

    std::string line;
    while (std::getline(f, line) && line.size() >= 2 && line[0] == '#' && line[1] == '\t')
        ;

    if (line.empty() || line.find("lcp") == std::string::npos) return false;

    bwt_heads.clear();
    bwt_run_lengths.clear();
    lcps_per_run.clear();

    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::string id_s, off_s, len_s, c_s, sa_s, lcp_s;
        if (!(iss >> id_s >> off_s >> len_s >> c_s >> sa_s)) continue;
        std::getline(iss, lcp_s);
        if (!lcp_s.empty() && lcp_s[0] == '\t') lcp_s = lcp_s.substr(1);

        ulint len = std::stoull(len_s);
        if (c_s.empty()) return false;
        uchar ch = map_movify_char(c_s[0]);

        bwt_heads.push_back(ch);
        bwt_run_lengths.push_back(len);
        lcps_per_run.push_back(parse_lcp_row(lcp_s));
    }
    return !bwt_heads.empty();
}

}  // namespace tsv
