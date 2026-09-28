/**
 * Block-parallel driver for batch queries: several threads each take the
 * next block of reads, query it, and format its output, and the formatted
 * blocks are written in input order.  Uses only the C++ standard library
 * (std::thread, std::mutex, std::condition_variable).
 *
 * Reading is serialized under one lock, since it parses a single stream.
 * Writing is serialized too, and done by whichever thread finds the next
 * block in input order ready, outside the lock that guards the queue of
 * finished blocks, so a thread depositing a block never waits for a write.
 * A thread waits before reading a new block only when that block would be
 * more than max_ahead blocks past the next one to write, which bounds the
 * memory held by finished blocks when one block is slow.
 *
 * Threads are not pinned to cores.  On machines with several memory nodes,
 * InterleavedAllocation spreads the index over all of them (see below).
 */

#ifndef _PARALLEL_BLOCKS_HPP
#define _PARALLEL_BLOCKS_HPP

#include <condition_variable>
#include <cstddef>
#include <exception>
#include <map>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#if defined(__linux__)
#include <sched.h>
#include <sys/syscall.h>
#include <unistd.h>
#endif

/**
 * Number of threads for a --threads value: n itself when positive, else the
 * number of CPUs this process may run on.  On Linux that is the affinity
 * mask, which schedulers such as Slurm and tools such as taskset and
 * numactl narrow; elsewhere it is the hardware's thread count (1 if it is
 * unknown).
 */
inline size_t resolve_thread_count(size_t n) {
    if (n > 0) return n;
#if defined(__linux__)
    cpu_set_t set;
    CPU_ZERO(&set);
    if (sched_getaffinity(0, sizeof(set), &set) == 0) {
        const int c = CPU_COUNT(&set);
        if (c > 0) return static_cast<size_t>(c);
    }
#endif
    const unsigned hw = std::thread::hardware_concurrency();
    return hw > 0 ? hw : 1;
}

/**
 * While it exists, the pages that the constructing thread allocates are
 * spread round-robin over the memory nodes the process may use, so that
 * threads on every node draw on the bandwidth of all nodes' memory rather
 * than of one node's.  By default Linux places each page on the node of the
 * thread that first touches it, which puts an index loaded by one thread on
 * one node.  It acts only on Linux with more than one allowed node, and only
 * when the thread has the default memory policy, so that placement chosen
 * with numactl (for example one process per node) is kept.  Elsewhere it
 * does nothing.
 */
class InterleavedAllocation {
public:
    InterleavedAllocation() {
#if defined(__linux__) && defined(SYS_get_mempolicy) && defined(SYS_set_mempolicy)
        constexpr int mpol_default = 0, mpol_interleave = 3, mpol_f_mems_allowed = 4;
        constexpr unsigned long bits = 1024;
        int mode = -1;
        if (syscall(SYS_get_mempolicy, &mode, nullptr, 0UL, nullptr, 0UL) != 0 || mode != mpol_default) return;
        unsigned long mask[bits / (8 * sizeof(unsigned long))] = {};
        if (syscall(SYS_get_mempolicy, nullptr, mask, bits, nullptr, static_cast<unsigned long>(mpol_f_mems_allowed)) != 0)
            return;
        int nodes = 0;
        for (unsigned long w : mask) nodes += __builtin_popcountl(w);
        if (nodes < 2) return;
        // The kernel reads maxnode - 1 bits of the mask.
        active_ = syscall(SYS_set_mempolicy, mpol_interleave, mask, bits + 1) == 0;
#endif
    }
    ~InterleavedAllocation() {
#if defined(__linux__) && defined(SYS_set_mempolicy)
        if (active_) syscall(SYS_set_mempolicy, 0, nullptr, 0UL);
#endif
    }
    InterleavedAllocation(const InterleavedAllocation&) = delete;
    InterleavedAllocation& operator=(const InterleavedAllocation&) = delete;

private:
    bool active_ = false;
};

/**
 * Runs the pipeline on `threads` threads (the calling thread alone when it
 * is 1).  Each thread owns a Block and a Worker:
 *
 * - make_worker(tid) returns thread tid's Worker, and is called on that
 *   thread, so a worker may hold thread-bound resources such as hardware
 *   counters;
 * - read(block) fills a block under the input lock, and returns false when
 *   the input is exhausted and the block is empty;
 * - worker.process(block, text) queries a block and appends its formatted
 *   output to text;
 * - write(text) writes one block's output, called for one block at a time,
 *   in input order;
 * - finish(worker) is called on each thread once its work is done, under a
 *   lock, to merge per-thread totals.
 *
 * An exception on any thread stops every thread from taking new blocks and
 * is rethrown on the calling thread.
 */
template <class Block, class MakeWorker, class Read, class Write, class Finish>
void run_parallel_blocks(size_t threads, size_t max_ahead, MakeWorker make_worker, Read read, Write write,
                         Finish finish) {
    if (threads < 1) threads = 1;
    if (max_ahead < threads) max_ahead = threads;

    std::mutex in_mu, out_mu, fin_mu;
    std::condition_variable room;
    size_t next_read = 0;    // id of the next block to read (under in_mu)
    bool input_done = false; // under in_mu
    size_t next_write = 0;   // id of the next block to write (under out_mu)
    bool writing = false;    // a thread is writing (under out_mu)
    bool failed = false;     // under out_mu
    std::map<size_t, std::string> ready;  // finished blocks not yet written
    std::exception_ptr error;

    // Write every ready block that is next in order.  Called with out_mu
    // held via lk; releases it around each write.
    auto drain = [&](std::unique_lock<std::mutex>& lk) {
        if (writing) return;
        writing = true;
        while (!failed && !ready.empty() && ready.begin()->first == next_write) {
            std::string text = std::move(ready.begin()->second);
            ready.erase(ready.begin());
            lk.unlock();
            write(text);
            lk.lock();
            ++next_write;
            room.notify_all();
        }
        writing = false;
    };

    auto run = [&](size_t tid) {
        try {
            auto worker = make_worker(tid);
            Block block;
            std::string text;
            for (;;) {
                size_t id;
                {
                    std::unique_lock<std::mutex> lk(in_mu);
                    if (input_done) break;
                    id = next_read;
                    {
                        // Wait for room while holding in_mu: the thread that
                        // would read block id has to wait anyway, and the
                        // others would only queue behind it.
                        std::unique_lock<std::mutex> olk(out_mu);
                        room.wait(olk, [&] { return failed || id < next_write + max_ahead; });
                        if (failed) break;
                    }
                    if (!read(block)) { input_done = true; break; }
                    ++next_read;
                }
                text.clear();
                worker.process(block, text);
                std::unique_lock<std::mutex> lk(out_mu);
                ready.emplace(id, std::move(text));
                text = std::string();
                drain(lk);
            }
            std::lock_guard<std::mutex> lk(fin_mu);
            finish(worker);
        } catch (...) {
            std::lock_guard<std::mutex> lk(out_mu);
            if (!failed) error = std::current_exception();
            failed = true;
            room.notify_all();
        }
    };

    if (threads == 1) {
        run(0);
    } else {
        std::vector<std::thread> pool;
        pool.reserve(threads);
        for (size_t t = 0; t < threads; ++t) pool.emplace_back(run, t);
        for (auto& th : pool) th.join();
    }
    if (error) std::rethrow_exception(error);
}

#endif /* _PARALLEL_BLOCKS_HPP */
