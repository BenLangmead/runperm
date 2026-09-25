// Microbenchmark: K independent chains of LF or psi steps over a TmsIndex,
// round robin, with each step split into start, prefetch and finish.
// Usage: tms_chain_bench INDEX [total_steps]
#include "tms_index.hpp"
#include <chrono>
#include <cstdio>
#include <fstream>
#include <random>

template <typename Pos, typename Start, typename Finish, typename Pf>
double run(std::vector<Pos> pos, size_t steps, Start start, Finish finish, Pf pf, bool prefetch, ulint& sink) {
    const size_t k = pos.size();
    for (auto& p : pos) { p = start(p); if (prefetch) pf(p.interval); }
    auto t0 = std::chrono::steady_clock::now();
    for (size_t s = 0; s < steps; ++s)
        for (size_t j = 0; j < k; ++j) {
            Pos p = finish(pos[j]);
            sink += p.interval;
            pos[j] = start(p);
            if (prefetch) pf(pos[j].interval);
        }
    auto t1 = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::nano>(t1 - t0).count() / double(steps * k);
}

int main(int argc, char** argv) {
    TmsIndex idx;
    std::ifstream in(argv[1], std::ios::binary);
    idx.load(in);
    const size_t total = argc > 2 ? std::stoul(argv[2]) : 2000000;
    std::mt19937_64 rng(1);
    ulint sink = 0;
    auto fresh_fl = [&](size_t k) {
        std::vector<TmsIndex::FLPos> v(k);
        for (auto& p : v) { p = TmsIndex::FLPos{}; p.interval = rng() % idx.psi_intervals(); }
        return v;
    };
    auto fresh_lf = [&](size_t k) {
        std::vector<TmsIndex::LFPos> v(k);
        for (auto& p : v) { p = TmsIndex::LFPos{}; p.interval = rng() % idx.move_runs(); }
        return v;
    };
    auto fs = [&](TmsIndex::FLPos p) { return idx.start_psi(p); };
    auto ff = [&](TmsIndex::FLPos p) { return idx.finish_psi(p); };
    auto fpf = [&](ulint i) { idx.prefetch_psi(i); };
    auto fpf2 = [&](ulint i) { idx.prefetch_psi(i); if (i + 1 < idx.psi_intervals()) idx.prefetch_psi(i + 1); };
    auto fpf_twice = [&](ulint i) { idx.prefetch_psi(i); idx.prefetch_psi(i); };
    auto fpf_next = [&](ulint i) { idx.prefetch_psi(i); idx.prefetch_psi(i + 1); };
    auto fpf_lf = [&](ulint i) { idx.prefetch_psi(i); idx.prefetch(i % idx.move_runs()); };
    auto ls = [&](TmsIndex::LFPos p) { return idx.start_LF(p); };
    auto lf = [&](TmsIndex::LFPos p) { return idx.finish_LF(p); };
    auto lpf = [&](ulint i) { idx.prefetch(i); };
    // Each run starts its chains at fresh random positions, so no run finds
    // another's lines in cache.
    for (size_t k : {1, 4, 16, 32, 64, 256}) {
        const size_t steps = total / k;
        double a = run(fresh_fl(k), steps, fs, ff, fpf, true, sink);
        double a2 = run(fresh_fl(k), steps, fs, ff, fpf2, true, sink);
        double b = run(fresh_fl(k), steps, fs, ff, fpf, false, sink);
        double c = run(fresh_lf(k), steps, ls, lf, lpf, true, sink);
        double d = run(fresh_lf(k), steps, ls, lf, lpf, false, sink);
        double e1 = run(fresh_fl(k), steps, fs, ff, fpf_twice, true, sink);
        double e2 = run(fresh_fl(k), steps, fs, ff, fpf_next, true, sink);
        double e3 = run(fresh_fl(k), steps, fs, ff, fpf_lf, true, sink);
        std::printf("K=%zu twice=%.1f next=%.1f plus_lf_row=%.1f ", k, e1, e2, e3);
        std::printf("K=%zu psi_pf=%.1f psi_pf2=%.1f psi_nopf=%.1f lf_pf=%.1f lf_nopf=%.1f ns/step\n", k, a, a2, b, c, d);
    }
    {
        // Independent random row reads, prefetched d reads ahead.
        const size_t N = 4000000;
        std::vector<ulint> fi(N), li(N);
        for (auto& x : fi) x = rng() % idx.psi_intervals();
        for (auto& x : li) x = rng() % idx.move_runs();
        for (size_t d : {0, 8, 32}) {
            auto t0 = std::chrono::steady_clock::now();
            for (size_t j = 0; j < N; ++j) { if (d && j + d < N) idx.prefetch_psi(fi[j + d]); sink += idx.lf().get_length(fi[j]); }
            auto t1 = std::chrono::steady_clock::now();
            for (size_t j = 0; j < N; ++j) { if (d && j + d < N) idx.prefetch(li[j + d]); sink += idx.lf().get_length(li[j]); }
            auto t2 = std::chrono::steady_clock::now();
            std::printf("random reads d=%zu fl=%.1f lf=%.1f ns\n", d,
                        std::chrono::duration<double, std::nano>(t1 - t0).count() / N,
                        std::chrono::duration<double, std::nano>(t2 - t1).count() / N);
        }
    }
    {
        // Fast-forward statistics of a single psi chain.
        TmsIndex::FLPos p{};
        p.interval = 12345;
        ulint ff = 0, far = 0, maxff = 0, maxoff = 0;
        for (size_t s = 0; s < 1000000; ++s) {
            auto u = idx.start_psi(p);
            maxoff = std::max<ulint>(maxoff, u.offset);
            p = idx.finish_psi(u);
            ff += p.interval - u.interval;
            far += p.interval - u.interval > 8;
            maxff = std::max<ulint>(maxff, p.interval - u.interval);
        }
        std::printf("psi ff/step=%.3f far_frac=%.4f max_ff=%lu max_unresolved_offset=%lu\n", ff / 1e6, far / 1e6,
                    (unsigned long)maxff, (unsigned long)maxoff);
    }
    std::printf("sink=%lu\n", (unsigned long)sink);
}
