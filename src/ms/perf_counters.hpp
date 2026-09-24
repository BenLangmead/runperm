/**
 * Hardware event counts around the query phase of ms and tms batches, for
 * diagnosing memory stalls.  The environment variable MS_PERF lists events
 * as name=config pairs separated by commas, where config is a raw
 * perf_event config in hex (event | umask << 8 | cmask << 24 on Intel), or
 * one of the names cycles, instructions and branch-misses with no config.
 * Counting is for this thread in user mode only.  On systems without
 * perf_event_open, or with MS_PERF unset, nothing is counted.
 */

#ifndef _PERF_COUNTERS_HPP
#define _PERF_COUNTERS_HPP

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#ifdef __linux__
#include <linux/perf_event.h>
#include <sys/ioctl.h>
#include <sys/syscall.h>
#include <unistd.h>
#endif

class PerfCounters {
public:
    PerfCounters() {
#ifdef __linux__
        const char* spec = std::getenv("MS_PERF");
        if (!spec) return;
        std::string s(spec);
        size_t at = 0;
        while (at < s.size()) {
            size_t end = s.find(',', at);
            if (end == std::string::npos) end = s.size();
            const std::string item = s.substr(at, end - at);
            at = end + 1;
            if (item.empty()) continue;
            perf_event_attr attr{};
            attr.size = sizeof(attr);
            attr.disabled = 1;
            attr.exclude_kernel = 1;
            attr.exclude_hv = 1;
            const size_t eq = item.find('=');
            const std::string name = item.substr(0, eq);
            if (eq == std::string::npos) {
                attr.type = PERF_TYPE_HARDWARE;
                if (name == "cycles") attr.config = PERF_COUNT_HW_CPU_CYCLES;
                else if (name == "instructions") attr.config = PERF_COUNT_HW_INSTRUCTIONS;
                else if (name == "branch-misses") attr.config = PERF_COUNT_HW_BRANCH_MISSES;
                else { std::cerr << "perf: unknown event " << name << "\n"; continue; }
            } else {
                attr.type = PERF_TYPE_RAW;
                attr.config = std::stoull(item.substr(eq + 1), nullptr, 16);
            }
            const int fd = static_cast<int>(syscall(SYS_perf_event_open, &attr, 0, -1, -1, 0));
            if (fd < 0) { std::cerr << "perf: cannot open " << name << "\n"; continue; }
            names_.push_back(name);
            fds_.push_back(fd);
            totals_.push_back(0);
        }
#endif
    }
    ~PerfCounters() {
#ifdef __linux__
        for (int fd : fds_) close(fd);
#endif
    }
    void start() {
#ifdef __linux__
        for (int fd : fds_) { ioctl(fd, PERF_EVENT_IOC_RESET, 0); ioctl(fd, PERF_EVENT_IOC_ENABLE, 0); }
#endif
    }
    void stop() {
#ifdef __linux__
        for (size_t i = 0; i < fds_.size(); ++i) {
            ioctl(fds_[i], PERF_EVENT_IOC_DISABLE, 0);
            uint64_t v = 0;
            if (read(fds_[i], &v, sizeof(v)) == sizeof(v)) totals_[i] += v;
        }
#endif
    }
    /** One line with each event's count per base. */
    void report(std::ostream& os, size_t bases) const {
        if (names_.empty()) return;
        os << "perf:";
        for (size_t i = 0; i < names_.size(); ++i) os << " " << names_[i] << "/base=" << double(totals_[i]) / double(bases);
        os << "\n";
    }

private:
    std::vector<std::string> names_;
    std::vector<int> fds_;
    std::vector<uint64_t> totals_;
};

#endif /* _PERF_COUNTERS_HPP */
