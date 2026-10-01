/**
 * Memory reports at the phase boundaries of an index build, on stderr: for
 * each phase, the resident set size at its end and the peak since the
 * previous report (Linux; elsewhere nothing is printed).  The peak is reset
 * at each report by writing 5 to /proc/self/clear_refs, so every report's
 * peak covers its own phase only.
 *
 * Author: Ben Langmead (ben.langmead@gmail.com)
 */

#ifndef _MS_MEM_REPORT_HPP
#define _MS_MEM_REPORT_HPP

#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

namespace mem_report {

/** VmRSS and VmHWM of this process in bytes, or 0 where unknown. */
inline void read_status(size_t& rss, size_t& hwm) {
    rss = hwm = 0;
    std::ifstream in("/proc/self/status");
    std::string line;
    while (std::getline(in, line)) {
        size_t kib = 0;
        if (std::sscanf(line.c_str(), "VmRSS: %zu kB", &kib) == 1) rss = kib * 1024;
        else if (std::sscanf(line.c_str(), "VmHWM: %zu kB", &kib) == 1) hwm = kib * 1024;
    }
}

/**
 * Print the phase's name, the seconds since the first report, the RSS now
 * and the peak since the last report, as a tab-separated line starting with
 * "mem", and reset the peak.  per gives a count (such as runs) to also report
 * the peak in bytes per unit.
 */
inline void phase(const char* name, size_t per = 0) {
    static const auto t0 = std::chrono::steady_clock::now();
    size_t rss, hwm;
    read_status(rss, hwm);
    if (hwm == 0) return;
    const double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    const double gib = 1024.0 * 1024.0 * 1024.0;
    char buf[256];
    std::snprintf(buf, sizeof(buf), "mem\t%s\tt=%.1fs\trss=%.3fGiB\tpeak=%.3fGiB", name, s, rss / gib, hwm / gib);
    std::cerr << buf;
    if (per > 0) {
        std::snprintf(buf, sizeof(buf), "\tpeak_per=%.1fB", static_cast<double>(hwm) / static_cast<double>(per));
        std::cerr << buf;
    }
    std::cerr << std::endl;
    std::ofstream clear("/proc/self/clear_refs");
    if (clear.good()) clear << "5";
}

}  // namespace mem_report

#endif /* _MS_MEM_REPORT_HPP */
