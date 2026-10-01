#!/usr/bin/env python3
"""
Round-trip tests for fmd2rlbwt.

Random plain BWT strings over $ACGTN are turned into FMDs with
`ropebwt3 plain2fmd`; fmd2rlbwt's heads and lengths must equal the strings'
maximal runs.  Random per-row LCP arrays give a minima file in TeraLCP's
-f fmd numbering (sentinel runs split into length-1 runs), and
fmd2rlbwt --minima-in/--minima-out must turn it into the file the same
definitions give for the maximal runs.

Usage: test_fmd2rlbwt.py ROPEBWT3 FMD2RLBWT [TRIALS]
"""
import os
import random
import struct
import subprocess
import sys
import tempfile

HEAD = {'$': b'$', 'A': b'A', 'C': b'C', 'G': b'G', 'T': b'T', 'N': b'\x7f'}


def uleb(v):
    out = bytearray()
    while True:
        b = v & 0x7f
        v >>= 7
        out.append(b | (0x80 if v else 0))
        if not v:
            return bytes(out)


def minima_file(runs, lcp, n, flags, input_runs):
    """TLMSM v1 bytes for runs [(start, length)] over per-row lcp."""
    recs = bytearray()
    total = 0
    inf = float('inf')
    for i, (s, ln) in enumerate(runs):
        top = lcp[s]
        nxt = lcp[runs[i + 1][0]] if i + 1 < len(runs) else inf
        keep = set()
        m = top
        for o in range(1, ln):
            if lcp[s + o] < m:
                keep.add(o)
                m = lcp[s + o]
        m = nxt
        for o in range(ln - 1, 0, -1):
            if lcp[s + o] < m:
                keep.add(o)
                m = lcp[s + o]
        recs += uleb(top) + uleb(len(keep))
        prev = 0
        for o in sorted(keep):
            recs += uleb(o - prev) + uleb(lcp[s + o])
            prev = o
        total += len(keep)
    hdr = bytes([0x93]) + b'TLMSM\x00\x01' + bytes([flags]) + bytes(7)
    hdr += struct.pack('<QQQQ', len(runs), n, total, input_runs)
    return hdr + bytes(recs)


def maximal_runs(s):
    runs = []
    for i, ch in enumerate(s):
        if i and ch == s[i - 1]:
            runs[-1][1] += 1
        else:
            runs.append([i, 1])
    return [tuple(r) for r in runs]


def teralcp_runs(s):
    out = []
    for st, ln in maximal_runs(s):
        if s[st] == '$':
            out += [(st + k, 1) for k in range(ln)]
        else:
            out.append((st, ln))
    return out


def main():
    rb3, tool = sys.argv[1], sys.argv[2]
    trials = int(sys.argv[3]) if len(sys.argv) > 3 else 200
    rng = random.Random(7)
    checked = 0
    with tempfile.TemporaryDirectory() as d:
        for t in range(trials):
            alpha = rng.choice(['$ACGTN', '$ACGT', '$$$ACGTN', 'ACGT$NNNN'])
            parts = []
            for _ in range(rng.randint(1, 60)):
                parts.append(rng.choice(alpha) * rng.choice([1, 1, 2, 3, rng.randint(1, 300)]))
            s = ''.join(parts)
            if t % 10 == 0:
                s += 'A' * rng.randint(1000, 100000)
            plain = os.path.join(d, 'p.txt')
            with open(plain, 'w') as f:
                f.write(s)
            fmd = os.path.join(d, 'p.fmd')
            subprocess.run([rb3, 'plain2fmd', '-o', fmd, plain], check=True)
            n = len(s)
            lcp = [0] + [rng.choice([0, 1, 2, 5, rng.randint(0, 1 << rng.choice([3, 20, 45]))]) for _ in range(n - 1)]
            tr = teralcp_runs(s)
            mins_in = os.path.join(d, 'in.min')
            with open(mins_in, 'wb') as f:
                f.write(minima_file(tr, lcp, n, 3 | (4 if len(tr) != len(maximal_runs(s)) else 0), 0))
            base = os.path.join(d, 'out')
            res = subprocess.run([tool, fmd, base], check=True, capture_output=True, text=True)
            stats = dict(line.split('\t')[:2] for line in res.stdout.splitlines() if not line.startswith('symbol'))
            mr = maximal_runs(s)
            assert int(stats['maximal_runs']) == len(mr), (t, stats, len(mr))
            assert int(stats['n']) == n
            heads = open(base + '.bwt.heads', 'rb').read()
            lens = open(base + '.bwt.len', 'rb').read()
            assert heads == b''.join(HEAD[s[st]] for st, _ in mr), t
            assert lens == b''.join(ln.to_bytes(5, 'little') for _, ln in mr), t
            # The expected TeraLCP numbering assumes rld kept every run whole,
            # which small FMDs do; otherwise skip the minima check.
            if int(stats['teralcp_runs']) != len(tr):
                print('trial', t, 'rld split a run; minima check skipped', file=sys.stderr)
                continue
            mins_out = os.path.join(d, 'out.min')
            subprocess.run([tool, '--minima-in', mins_in, '--minima-out', mins_out, fmd, base],
                           check=True, capture_output=True, text=True)
            assert open(mins_out, 'rb').read() == minima_file(mr, lcp, n, 3, len(mr)), (t, s[:80])
            checked += 1
    print('ok', trials, 'trials,', checked, 'with minima')


if __name__ == '__main__':
    main()
