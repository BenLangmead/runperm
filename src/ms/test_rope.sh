#!/bin/bash
# Check tms-batch --report smem-all against ropebwt3 mem on a small FMD.
#
# A multi-contig reference with repeats (some reverse-complemented) is
# indexed by ropebwt3, converted with fmd2rlbwt and TeraLCP -ophi and built
# into a tms index.  For each cap N and minimum length, ropebwt3 mem -l X
# -p N and tms-batch --report smem-all --max-smem-positions N --min-smem-len
# X-1 --seq-ids FILE must give the same lines in the same order apart from
# the positions: the same read, start, end, count and n_pos.  For an SMEM
# with count <= N the positions must be the same set; above the cap each
# position tms-batch lists must be distinct and hold the SMEM's sequence
# (reverse-complemented on the - strand).  Without --seq-ids, the output must
# be the same but for text positions in place of the contig positions, and
# those must give the contig positions when mapped by hand.  A table that
# does not describe the index, or that splits it at other points so that a
# match crosses a contig's end, must make tms-batch fail, as must --seq-ids
# without smem-all and smem-all with --positions.
#
# Usage: test_rope.sh MS FMD2RLBWT TERALCP [ROPEBWT3] [WORK_DIR]
# TERALCP must support -ophi.
set -euo pipefail
MS=$1
FMD2RLBWT=$2
TERALCP=$3
RB3=${4:-ropebwt3}
W=${5:-$(mktemp -d)}
mkdir -p "$W"

python3 - "$W" <<'EOF'
import gzip, random, sys
w = sys.argv[1]
random.seed(3)
def rc(s): return s[::-1].translate(str.maketrans('ACGT', 'TGCA'))
rep = ''.join(random.choice('ACGT') for _ in range(60))
cs = []
for k in range(9):
    s = ''.join(random.choice('ACGT') for _ in range(random.randint(150, 2500)))
    for _ in range(random.randint(1, 5)):
        x = random.randint(0, len(s) - 60)
        s = s[:x] + (rep if random.random() < 0.6 else rc(rep)) + s[x + 60:]
    cs.append(s)
with open(f'{w}/ref.fa', 'w') as f:
    for k, s in enumerate(cs):
        f.write(f'>chr{k + 1} contig {k + 1}\n')
        for i in range(0, len(s), 70): f.write(s[i:i + 70] + '\n')
# ropebwt3's sequence-ID table: name and length per contig.
with gzip.open(f'{w}/ref.fmd.len.gz', 'wt') as f:
    for k, s in enumerate(cs): f.write(f'chr{k + 1}\t{len(s)}\n')
# The same total length split at other points, so that matches cross ends.
with gzip.open(f'{w}/shifted.len.gz', 'wt') as f:
    for k, s in enumerate(cs): f.write(f'chr{k + 1}\t{len(s) + (7 if k == 0 else -7 if k == 1 else 0)}\n')
with gzip.open(f'{w}/short.len.gz', 'wt') as f:
    for k, s in enumerate(cs[:-1]): f.write(f'chr{k + 1}\t{len(s)}\n')
random.seed(11)
with open(f'{w}/reads.fa', 'w') as f:
    for i in range(400):
        c = random.choice(cs)
        L = random.randint(30, 150)
        x = random.randint(0, max(0, len(c) - L))
        s = c[x:x + L]
        if random.random() < 0.5: s = rc(s)
        s = list(s)
        for _ in range(random.randint(0, 4)): s[random.randrange(len(s))] = random.choice('ACGT')
        if i % 50 == 7: s = [random.choice('ACGT') for _ in range(L)]
        if i % 50 == 13: s[len(s) // 2] = 'N'
        f.write(f'>r{i} extra\n{"".join(s)}\n')
    f.write('>allN\nNNNNNNNNNNNN\n')
EOF

"$RB3" build -d -o "$W/ref.fmd" "$W/ref.fa" 2>/dev/null
"$RB3" ssa -o "$W/ref.fmd.ssa" "$W/ref.fmd" 2>/dev/null
"$FMD2RLBWT" "$W/ref.fmd" "$W/ref" >/dev/null
"$TERALCP" -f fmd -i "$W/ref.fmd" -ophi "$W/ref.phi" -t "$W/teralcp.tmp" -v quiet
"$MS" tms-build "$W/ref.bwt.heads" "$W/ref.bwt.len" "$W/ref.idx" --phi "$W/ref.phi" >/dev/null 2>&1

for cap in 1 3 100 0; do
    for minlen in 1 20; do
        "$RB3" mem -l "$minlen" -p "$cap" "$W/ref.fmd" "$W/reads.fa" 2>/dev/null > "$W/rope.$cap.$minlen"
        "$MS" tms-batch "$W/ref.idx" "$W/reads.fa" --report smem-all --max-smem-positions "$cap" \
            --min-smem-len $((minlen - 1)) --seq-ids "$W/ref.fmd.len.gz" \
            -o "$W/orbit.$cap.$minlen" 2>/dev/null
        "$MS" tms-batch "$W/ref.idx" "$W/reads.fa" --report smem-all --max-smem-positions "$cap" \
            --min-smem-len $((minlen - 1)) -o "$W/raw.$cap.$minlen" 2>/dev/null
        python3 - "$W" "$cap" "$minlen" <<'EOF'
import gzip, sys
w, cap, minlen = sys.argv[1], int(sys.argv[2]), sys.argv[3]
def rc(s): return s[::-1].translate(str.maketrans('ACGT', 'TGCA'))
ref, name = {}, None
for l in open(f'{w}/ref.fa'):
    if l[0] == '>': name = l[1:].split()[0]; ref[name] = []
    else: ref[name].append(l.strip())
ref = {k: ''.join(v) for k, v in ref.items()}
reads, name = {}, None
for l in open(f'{w}/reads.fa'):
    if l[0] == '>': name = l[1:].split()[0]
    else: reads[name] = l.strip()
names = [l.split('\t')[0] for l in gzip.open(f'{w}/ref.fmd.len.gz', 'rt')]
def by_hand(p, L):
    at = 0
    for n in names:
        for strand in '+-':
            m = len(ref[n])
            if at <= p < at + m + 1:
                local = p - at
                assert local + L <= m, 'raw position crosses a contig end'
                return (n, strand, local if strand == '+' else m - (local + L))
            at += m + 1
    raise AssertionError('raw position past the text')
def parse(path):
    out = []
    for l in open(path):
        f = l.rstrip('\n').split('\t')
        st, en, cnt = int(f[1]), int(f[2]), int(f[3])
        npos = int(f[4]) if len(f) > 4 else 0
        pos = [tuple(x.rsplit(':', 2)) for x in f[5:]]
        pos = [(c, s, int(o)) for c, s, o in pos]
        assert len(pos) == npos, l
        out.append((f[0], st, en, cnt, npos, pos))
    return out
rope, orbit = parse(f'{w}/rope.{cap}.{minlen}'), parse(f'{w}/orbit.{cap}.{minlen}')
assert [x[:5] for x in rope] == [x[:5] for x in orbit], 'intervals, counts or n_pos differ'
sets = verified = 0
for a, b in zip(rope, orbit):
    nm, st, en, cnt, npos, _ = a
    assert npos == min(cap, cnt)
    assert len(set(b[5])) == len(b[5]), 'duplicate position'
    if cnt <= cap:
        assert set(a[5]) == set(b[5]), f'position sets differ: {a} {b}'
        sets += 1
    q = reads[nm][st:en]
    for c, s, o in b[5]:
        t = ref[c][o:o + en - st]
        assert (t if s == '+' else rc(t)) == q, f'{c}:{s}:{o} does not hold {q}'
        verified += 1
# Text positions mapped by hand give the contig positions.
raw = [l.rstrip('\n').split('\t') for l in open(f'{w}/raw.{cap}.{minlen}')]
assert len(raw) == len(orbit)
for f, o in zip(raw, orbit):
    head = (f[0], int(f[1]), int(f[2]), int(f[3]), int(f[4]) if len(f) > 4 else 0)
    assert head == o[:5], f'raw and contig lines differ at {f[0]}'
    st, en = int(f[1]), int(f[2])
    assert [by_hand(int(p), en - st) for p in f[5:]] == o[5], f'raw and contig positions differ at {f[0]}'
smems = [x for x in orbit if x[3] > 0]
zeros = len(orbit) - len(smems)
print(f'cap {cap} -l {minlen}: {len(smems)} SMEMs and {zeros} unmatched-character lines agree, {sets} complete position sets, {verified} positions verified')
EOF
    done
done

# A table for another text, and one that crosses contig ends, are refused.
for bad in short shifted; do
    if "$MS" tms-batch "$W/ref.idx" "$W/reads.fa" --report smem-all --max-smem-positions 3 \
        --seq-ids "$W/$bad.len.gz" -o "$W/bad.out" 2> "$W/bad.err"; then
        echo "tms-batch accepted $bad.len.gz"; exit 1
    fi
    echo "$bad.len.gz refused: $(tail -1 "$W/bad.err")"
done
grep -q crosses "$W/bad.err" || { echo "the shifted table was not refused for crossing a contig end"; exit 1; }
for args in "--report smem-one --seq-ids $W/ref.fmd.len.gz" "--seq-ids $W/ref.fmd.len.gz" \
    "--report smem-all --positions"; do
    # shellcheck disable=SC2086
    if "$MS" tms-batch "$W/ref.idx" "$W/reads.fa" $args -o "$W/bad.out" 2>/dev/null; then
        echo "tms-batch accepted '$args'"; exit 1
    fi
done
echo "test_rope: all passed"
