#!/bin/bash
# Check query --report smem-all against ropebwt3 mem on a small FMD.
#
# A multi-contig reference with repeats (some reverse-complemented) is
# indexed by ropebwt3, converted with fmd2rlbwt and TeraLCP -ophi and built
# into a tms index twice: with tms-build --seq-ids and its .len.gz, and
# without, the table then added by index-add-seq-ids.  The two must be equal
# apart from their provenance, and the second must keep the first index's
# sections byte for byte.  For each cap N and minimum length, ropebwt3 mem
# -l X -p N and query --report smem-all --max-smem-positions N
# --min-smem-len X must give the same lines in the same order apart from
# the positions: the same read, start, end, count and n_pos.  For an SMEM
# with count <= N the positions must be the same set; above the cap each
# position query lists must be distinct and hold the SMEM's sequence
# (reverse-complemented on the - strand).  On the index without a table, the
# output must be the same but for text positions in place of the contig
# positions, and those must give the contig positions when mapped by hand.
# A table whose text length is not the index's must be refused by both
# tms-build and index-add-seq-ids, and one whose string count is not the
# text's by tms-build, which counts the strings; one that fits both but
# splits the text at other points must make query fail when a match
# crosses a contig's end.
# smem-all with --positions must be refused.  Every counting method (--count
# walk, backward and adaptive) must give the default's output byte for byte.
# A psi-locate index, built with tms-build --layout psi-locate and derived
# from the full one by index-psi-locate, must be the same file apart from
# provenance, and give the full index's output byte for byte for psi
# matching statistics, --positions, smem-one and smem-all with every counting
# method; matching statistics in the modes that read phi's lengths-based
# form (phi, phi-walk, phiskip, dual, dual-walk) must be refused.
#
# Usage: test_rope.sh MS FMD2RLBWT TERALCP [ROPEBWT3] [WORK_DIR]
# TERALCP must support -ophi and -ominima.
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
# The first two contigs as one, with the same text length but a string fewer
# on each strand.
with gzip.open(f'{w}/merged.len.gz', 'wt') as f:
    f.write(f'chr1\t{len(cs[0]) + len(cs[1]) + 1}\n')
    for k, s in enumerate(cs[2:]): f.write(f'chr{k + 3}\t{len(s)}\n')
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
"$TERALCP" -f fmd -i "$W/ref.fmd" -ophi "$W/ref.phi" -ominima "$W/ref.min" -t "$W/teralcp.tmp" -v quiet
"$MS" tms-build "$W/ref.bwt.heads" "$W/ref.bwt.len" "$W/plain.idx" --phi "$W/ref.phi" --minima "$W/ref.min" >/dev/null 2>&1
"$MS" tms-build "$W/ref.bwt.heads" "$W/ref.bwt.len" "$W/ref.idx" --phi "$W/ref.phi" --minima "$W/ref.min" \
    --seq-ids "$W/ref.fmd.len.gz" >/dev/null 2>&1
"$MS" index-add-seq-ids "$W/plain.idx" "$W/ref.fmd.len.gz" "$W/added.idx" >/dev/null
"$MS" index-verify "$W/added.idx" >/dev/null
python3 - "$W" <<'PYEOF'
import struct, sys
w = sys.argv[1]
def parse(p):
    b = open(p, 'rb').read()
    ns, ni = struct.unpack_from('<H', b, 40)[0], struct.unpack_from('<H', b, 44)[0]
    secs = [b[128 + 80 * i:128 + 80 * (i + 1)] for i in range(ns)]
    return b, b[:44] + b[46:128], secs, b[128 + 80 * ns:128 + 80 * ns + ni].decode()
plain, built, added = (parse(f'{w}/{x}.idx') for x in ('plain', 'ref', 'added'))
# Equal but for the info text and its size.
assert built[1] == added[1] and built[2] == added[2] and built[0][4096:] == added[0][4096:], \
    'tms-build --seq-ids and index-add-seq-ids differ'
# The copy keeps the sections of the index it was given, at their offsets.
assert added[2][:-1] == plain[2] and added[0][4096:len(plain[0])] == plain[0][4096:], \
    'index-add-seq-ids changed a section'
assert added[3].startswith(plain[3]) and 'seq_ids_from=' in added[3]
print('tms-build --seq-ids and index-add-seq-ids agree; the copied sections are unchanged')
PYEOF
"$MS" tms-build "$W/ref.bwt.heads" "$W/ref.bwt.len" "$W/loc.idx" --phi "$W/ref.phi" --minima "$W/ref.min" \
    --seq-ids "$W/ref.fmd.len.gz" --layout psi-locate >/dev/null 2>&1
"$MS" index-psi-locate "$W/ref.idx" "$W/derived.idx" >/dev/null
"$MS" index-verify "$W/derived.idx" >/dev/null
python3 - "$W" <<'PYEOF'
import struct, sys
w = sys.argv[1]
def parse(p):
    b = open(p, 'rb').read()
    ns, ni = struct.unpack_from('<H', b, 40)[0], struct.unpack_from('<H', b, 44)[0]
    secs = [b[128 + 80 * i:128 + 80 * (i + 1)] for i in range(ns)]
    return b, b[:44] + b[46:128], secs, b[128 + 80 * ns:128 + 80 * ns + ni].decode()
full, loc, derived = (parse(f'{w}/{x}.idx') for x in ('ref', 'loc', 'derived'))
assert struct.unpack_from('<H', loc[0], 22)[0] == 4, 'tms-build --layout psi-locate wrote another layout'
assert loc[1] == derived[1] and loc[2] == derived[2] and loc[0][4096:] == derived[0][4096:], \
    'tms-build --layout psi-locate and index-psi-locate differ'
ids = lambda secs: [struct.unpack_from('<H', e, 0)[0] for e in secs]
assert ids(full[2]) == [1, 2, 3, 4, 5, 7] and ids(loc[2]) == [1, 2, 4, 5, 7], 'psi-locate must drop PHI_LEN alone'
# The kept sections' entries match the full index's but for their offsets.
strip = lambda e: e[:8] + e[16:]
assert [strip(e) for e in loc[2]] == [strip(e) for e in full[2] if struct.unpack_from('<H', e, 0)[0] != 3]
assert 'psi_locate_by=' in derived[3]
print('tms-build --layout psi-locate and index-psi-locate agree; they hold the full sections but PHI_LEN')
PYEOF
for q in "" "--positions" "--report smem-one"; do
    "$MS" query "$W/ref.idx" "$W/reads.fa" $q -o "$W/full.q" 2>/dev/null
    for x in loc derived; do
        "$MS" query "$W/$x.idx" "$W/reads.fa" $q -o "$W/$x.q" 2>/dev/null
        cmp -s "$W/full.q" "$W/$x.q" || { echo "psi-locate ($x) output differs for '$q'"; exit 1; }
    done
done
for mode in phi phi-walk phiskip dual dual-walk; do
    if "$MS" query "$W/loc.idx" "$W/reads.fa" --mode $mode -o "$W/bad.out" 2> "$W/bad.err"; then
        echo "psi-locate accepted --mode $mode"; exit 1
    fi
    grep -q "lengths-based" "$W/bad.err" || { echo "psi-locate refused --mode $mode for another reason"; exit 1; }
done
"$MS" index-plan "$W/loc.idx" --mode dual | grep -q refused || { echo "index-plan did not refuse dual"; exit 1; }
"$MS" index-plan "$W/loc.idx" --mode psi --report smem-all | grep -q "LF F_STARTS PHI_START PHI_INV SEQ_IDS" ||
    { echo "index-plan does not load psi smem-all's sections"; exit 1; }
if "$MS" index-psi-locate "$W/loc.idx" "$W/bad.idx" 2>/dev/null; then echo "index-psi-locate accepted a psi-locate index"; exit 1; fi
echo "psi-locate: psi MS, --positions and smem-one equal the full index's; phi, phi-walk, phiskip, dual and dual-walk refused"

for cap in 1 3 100 0; do
    for minlen in 1 20; do
        "$RB3" mem -l "$minlen" -p "$cap" "$W/ref.fmd" "$W/reads.fa" 2>/dev/null > "$W/rope.$cap.$minlen"
        "$MS" query "$W/ref.idx" "$W/reads.fa" --report smem-all --max-smem-positions "$cap" \
            --min-smem-len "$minlen" -o "$W/orbit.$cap.$minlen" 2>/dev/null
        "$MS" query "$W/added.idx" "$W/reads.fa" --report smem-all --max-smem-positions "$cap" \
            --min-smem-len "$minlen" -o "$W/added.$cap.$minlen" 2>/dev/null
        cmp -s "$W/orbit.$cap.$minlen" "$W/added.$cap.$minlen" || { echo "the two indexes' outputs differ"; exit 1; }
        # Every counting method, the walk, backward search alone or after a
        # short walk, gives the default's output byte for byte.
        for count in walk backward adaptive:0 adaptive:3 adaptive:3:12 adaptive:64; do
            IFS=: read -r method budget short <<< "$count"
            "$MS" query "$W/ref.idx" "$W/reads.fa" --report smem-all --max-smem-positions "$cap" \
                --min-smem-len "$minlen" --count "$method" ${budget:+--count-budget "$budget"} \
                ${short:+--count-short-len "$short"} --count-row-shift 1 -o "$W/count.$cap.$minlen" 2>/dev/null
            cmp -s "$W/orbit.$cap.$minlen" "$W/count.$cap.$minlen" || { echo "--count $count differs from the default"; exit 1; }
            for x in loc derived; do
                "$MS" query "$W/$x.idx" "$W/reads.fa" --report smem-all --max-smem-positions "$cap" \
                    --min-smem-len "$minlen" --count "$method" ${budget:+--count-budget "$budget"} \
                    ${short:+--count-short-len "$short"} --count-row-shift 1 -o "$W/count.$cap.$minlen" 2>/dev/null
                cmp -s "$W/orbit.$cap.$minlen" "$W/count.$cap.$minlen" ||
                    { echo "psi-locate ($x) --count $count differs from the full index"; exit 1; }
            done
        done
        # Every LCE mode gives psi's output, phi and dual reading the LCP
        # summaries built from the FMD's minima.
        for mode in phi dual phi-walk dual-walk; do
            "$MS" query "$W/ref.idx" "$W/reads.fa" --report smem-all --max-smem-positions "$cap" \
                --min-smem-len "$minlen" --mode "$mode" -o "$W/mode.$cap.$minlen" 2>/dev/null
            cmp -s "$W/orbit.$cap.$minlen" "$W/mode.$cap.$minlen" || { echo "--mode $mode differs from psi"; exit 1; }
        done
        "$MS" query "$W/plain.idx" "$W/reads.fa" --report smem-all --max-smem-positions "$cap" \
            --min-smem-len "$minlen" -o "$W/raw.$cap.$minlen" 2>/dev/null
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

# A table for another text is refused when it is stored.
for cmd in build add; do
    if [ $cmd = build ]; then
        run=("$MS" tms-build "$W/ref.bwt.heads" "$W/ref.bwt.len" "$W/bad.idx" --phi "$W/ref.phi" --minima "$W/ref.min" --seq-ids "$W/short.len.gz")
    else
        run=("$MS" index-add-seq-ids "$W/plain.idx" "$W/short.len.gz" "$W/bad.idx")
    fi
    if "${run[@]}" > "$W/bad.err" 2>&1; then echo "$cmd accepted short.len.gz"; exit 1; fi
    echo "$cmd refused short.len.gz: $(tail -1 "$W/bad.err" | sed "s|$W/||g")"
done
# tms-build counts the text's strings, which the header does not give.
if "$MS" tms-build "$W/ref.bwt.heads" "$W/ref.bwt.len" "$W/bad.idx" --phi "$W/ref.phi" --minima "$W/ref.min" \
    --seq-ids "$W/merged.len.gz" > "$W/bad.err" 2>&1; then
    echo "tms-build accepted a table with a string count other than the text's"; exit 1
fi
grep -q strings "$W/bad.err" || { echo "the merged table was not refused for its string count"; exit 1; }
echo "build refused merged.len.gz: $(tail -1 "$W/bad.err" | sed "s|$W/||g")"
# A table that fits the text's length and string count but splits it
# elsewhere is stored, and query fails when a match crosses an end.
"$MS" tms-build "$W/ref.bwt.heads" "$W/ref.bwt.len" "$W/shifted.idx" --phi "$W/ref.phi" --minima "$W/ref.min" \
    --seq-ids "$W/shifted.len.gz" >/dev/null 2>&1
if "$MS" query "$W/shifted.idx" "$W/reads.fa" --report smem-all --max-smem-positions 3 \
    -o "$W/bad.out" 2> "$W/bad.err"; then
    echo "query accepted a table that splits the text elsewhere"; exit 1
fi
grep -q crosses "$W/bad.err" || { echo "the shifted table did not fail for crossing a contig end"; exit 1; }
echo "shifted table: $(tail -1 "$W/bad.err")"
if "$MS" query "$W/ref.idx" "$W/reads.fa" --report smem-all --positions -o "$W/bad.out" 2>/dev/null; then
    echo "query accepted smem-all with --positions"; exit 1
fi
echo "test_rope: all passed"
