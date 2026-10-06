#!/bin/bash
# Checks that batch and tms-batch output is byte-identical for every
# combination of interleave K, thread count and block size, on indexes built
# from DATA_DIR/minishred1_20_002_lcp.tsv and reads drawn from its text with
# substitutions, plus random reads, empty reads and reads with N.
# Usage: test_threads.sh MS_BINARY DATA_DIR [WORK_DIR]
set -u
M=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
DATA=$2
W=${3:-$(mktemp -d)}
mkdir -p "$W"
$M build "$DATA/minishred1_20_002_lcp.tsv" "$W/mini.idx" > /dev/null 2>&1 || { echo "build failed"; exit 1; }
$M tms-build-tsv "$DATA/minishred1_20_002_lcp.tsv" "$W/mini.tms" --phi > /dev/null 2>&1 || { echo "tms-build failed"; exit 1; }
$M tms-text "$W/mini.tms" 2> /dev/null > "$W/text.txt"
python3 - "$W" <<'PY'
import random, sys
d = sys.argv[1]
t = open(d + '/text.txt').read().strip().replace('%', '').replace('$', '')
random.seed(7)
with open(d + '/reads.fa', 'w') as out:
    for i in range(6000):
        L = random.choice([0, 1, 5, 50, 150, 150, 150, 400])
        if L and random.random() < 0.8 and L < len(t):
            p = random.randrange(len(t) - L)
            s = list(t[p:p + L])
            er = random.choice([0, 0.002, 0.01, 0.05])
            for j in range(L):
                if random.random() < er: s[j] = random.choice('ACGTN')
            s = ''.join(s)
        else:
            s = ''.join(random.choice('ACGT') for _ in range(L))
        out.write('>r%d\n%s\n' % (i, s))
PY
cd "$W"
fails=0
while read -r cmd idx opts; do
    $M $cmd $idx reads.fa $opts -o ref.out --interleave 0 2> err.txt || { echo "failed: $cmd $opts"; cat err.txt; exit 1; }
    for k in 0 1 7 32; do for t in 1 2 3 8 17; do for b in 1 13 1000; do
        $M $cmd $idx reads.fa $opts -o o.out --interleave $k --threads $t --block-reads $b 2> err.txt \
            || { echo "failed: $cmd $opts k=$k t=$t b=$b"; cat err.txt; fails=$((fails + 1)); continue; }
        cmp -s ref.out o.out || { echo "differs: $cmd $opts k=$k t=$t b=$b"; fails=$((fails + 1)); }
    done; done; done
    echo "checked $cmd $opts"
done <<'CMDS'
batch mini.idx
tms-batch mini.tms
tms-batch mini.tms --mode phiskip --positions
tms-batch mini.tms --mode phi --positions
tms-batch mini.tms --report smem-all
tms-batch mini.tms --mode dual --report smem-one
batch mini.idx --orient rc
tms-batch mini.tms --mode dual --orient rc
CMDS
# The same reads as multi-line FASTA, FASTQ and one per line, with CRLF
# line endings, blank lines, lower case, and odd records at the end: output
# must not depend on thread count or block size, and the sequences' values
# must agree across formats.
python3 - "$W" <<'PY'
import random, sys
d = sys.argv[1]
t = open(d + '/text.txt').read().strip().replace('%', '').replace('$', '')
r = random.Random(11)
reads = []
for i in range(3000):
    L = r.choice([0, 1, 20, 150, 150, 300, 1000])
    p = r.randrange(len(t) - L)
    s = ''.join(c if r.random() > 0.01 else r.choice('acgtnACGTN') for c in t[p:p + L])
    reads.append(s.lower() if r.random() < 0.2 else s)
def nl(): return '\r\n' if r.random() < 0.3 else '\n'
with open(d + '/m.fa', 'w', newline='') as f:
    f.write(nl() * 2)
    for i, s in enumerate(reads):
        f.write('>r%d%s desc%s' % (i, r.choice(['', ' x', '\ty']), nl()))
        w = r.choice([60, 77, 5000])
        for j in range(0, len(s), w): f.write(s[j:j + w] + nl())
        if r.random() < 0.1: f.write(nl())
with open(d + '/m.fq', 'w', newline='') as f:
    for i, s in enumerate(reads):
        if r.random() < 0.05: f.write(nl())
        f.write('@q%d extra%s%s%s+%s%s%s' % (i, nl(), s, nl(), nl(), '@' * len(s), nl()))
    f.write('@trunc\nACGT\n@dangling\n')
with open(d + '/m.txt', 'w', newline='') as f:
    for s in reads:
        f.write(s + nl())
        if r.random() < 0.05: f.write(nl())
    f.write('ACGTACGTNNacgt')
PY
for f in m.fa m.fq m.txt; do
    $M batch mini.idx $f -o ref.out 2> err.txt || { echo "failed: $f"; cat err.txt; exit 1; }
    cp ref.out ref.$f
    for t in 1 3 8; do for b in 1 7 1000; do
        $M batch mini.idx $f -o o.out --threads $t --block-reads $b 2> /dev/null
        cmp -s ref.out o.out || { echo "differs: $f t=$t b=$b"; fails=$((fails + 1)); }
    done; done
    for t in 1 3 8; do for bb in 1 7 100 5000 0; do
        $M batch mini.idx $f -o o.out --threads $t --block-bytes $bb 2> /dev/null
        cmp -s ref.out o.out || { echo "differs: $f t=$t block-bytes=$bb"; fails=$((fails + 1)); }
    done; done
done
# A record much longer than the reader's 4 MB segments.
python3 -c "
import random; r = random.Random(5)
with open('long.fa', 'w') as f:
    f.write('>a\\nACGTTGCA\\n>big\\n')
    s = ''.join(r.choice('ACGT') for _ in range(10000000))
    for i in range(0, len(s), 80): f.write(s[i:i + 80] + '\\n')
    f.write('>c\\nGATTACA\\n')
"
$M batch mini.idx long.fa -o ref.out --threads 1 2> /dev/null
[ $(wc -l < ref.out) -eq 3 ] || { echo "long record: wrong line count"; fails=$((fails + 1)); }
$M batch mini.idx long.fa -o o.out --threads 3 --block-reads 1 2> /dev/null
cmp -s ref.out o.out || { echo "differs: long record"; fails=$((fails + 1)); }
for bb in 1 1000; do
    $M batch mini.idx long.fa -o o.out --threads 3 --block-bytes $bb 2> /dev/null
    cmp -s ref.out o.out || { echo "differs: long record, block-bytes=$bb"; fails=$((fails + 1)); }
done
# Output named *.gz is one gzip member per block, and decompresses to the
# plain output.
for t in 1 4; do
    $M batch mini.idx m.fq -o o.gz --threads $t --block-reads 50 2> /dev/null
    gzip -dc o.gz | cmp -s - ref.m.fq || { echo "differs: gzip output t=$t"; fails=$((fails + 1)); }
done
# Empty reads are skipped only by the one-per-line format.
cut -f2 ref.m.fa | grep -v '^$' > a.txt
head -n 3000 ref.m.fq | cut -f2 | grep -v '^$' > b.txt
head -n $(wc -l < a.txt) ref.m.txt | cut -f2 > c.txt
cmp -s a.txt b.txt && cmp -s a.txt c.txt || { echo "formats disagree"; fails=$((fails + 1)); }
echo "checked read formats"

# --orient qual on FASTQ with varied qualities: output must not depend on
# thread count, block size or K, and the flipped count must follow the rule
# (mean quality of the last min(20, len) bases above that of the first).
python3 - "$W" <<'PY'
import random, sys
d = sys.argv[1]
t = open(d + '/text.txt').read().strip().replace('%', '').replace('$', '')
r = random.Random(13)
flips = 0
with open(d + '/q.fq', 'w') as f:
    for i in range(2000):
        L = r.choice([0, 1, 15, 20, 21, 39, 40, 150, 150])
        p = r.randrange(len(t) - L)
        s = ''.join(c if r.random() > 0.02 else r.choice('ACGTN') for c in t[p:p + L])
        q = ''.join(chr(33 + r.randrange(42)) for _ in range(L))
        m = min(20, L)
        flips += sum(map(ord, q[L - m:])) > sum(map(ord, q[:m]))
        f.write('@q%d\n%s\n+\n%s\n' % (i, s, q))
open(d + '/q.flips', 'w').write('%d\n' % flips)
PY
$M batch mini.idx q.fq -o ref.out --orient qual --interleave 0 2> err.txt || { echo "failed: orient qual"; cat err.txt; exit 1; }
grep -q "flipped=$(cat q.flips) " err.txt || { echo "orient qual: wrong flipped count"; cat err.txt; fails=$((fails + 1)); }
for k in 1 32; do for t in 1 3 8; do for b in 1 13 1000; do
    $M batch mini.idx q.fq -o o.out --orient qual --interleave $k --threads $t --block-reads $b 2> /dev/null
    cmp -s ref.out o.out || { echo "differs: orient qual k=$k t=$t b=$b"; fails=$((fails + 1)); }
done; done; done
# Refused: qual without qualities, and positions or SMEM reports with rc.
$M batch mini.idx m.fa -o o.out --orient qual 2> /dev/null && { echo "orient qual accepted FASTA"; fails=$((fails + 1)); }
$M tms-batch mini.tms m.fa -o o.out --orient rc --mode phi --positions 2> /dev/null && { echo "orient rc accepted --positions"; fails=$((fails + 1)); }
$M tms-batch mini.tms m.fa -o o.out --orient rc --report smem-one 2> /dev/null && { echo "orient rc accepted an SMEM report"; fails=$((fails + 1)); }
echo "checked orientation"

# The oracle mode: outcomes recorded in any mode, with any K, T and B, are
# the same file, and replaying it with any K, T and B gives the recording
# modes' output, with and without positions, as does recording and replaying
# the reverse complements with --orient rc.
$M tms-batch mini.tms reads.fa -o ref.out --oracle-write ref.orc 2> err.txt || { echo "failed: oracle-write"; cat err.txt; exit 1; }
$M tms-batch mini.tms reads.fa -o ref.pos --mode phi --positions 2> /dev/null
for rec in "psi 0 1 1000" "phi 7 3 13" "dual 32 8 1" "phiskip 1 17 13"; do
    set -- $rec
    $M tms-batch mini.tms reads.fa -o o.out --mode $1 --interleave $2 --threads $3 --block-reads $4 --oracle-write o.orc 2> /dev/null
    cmp -s ref.orc o.orc || { echo "oracle differs: recorded with $rec"; fails=$((fails + 1)); }
    cmp -s ref.out o.out || { echo "differs: recording with $rec"; fails=$((fails + 1)); }
done
for k in 0 1 7 32; do for t in 1 3 8; do for b in 1 13 1000; do
    $M tms-batch mini.tms reads.fa -o o.out --mode oracle --oracle ref.orc --interleave $k --threads $t --block-reads $b 2> err.txt \
        || { echo "failed: oracle k=$k t=$t b=$b"; cat err.txt; fails=$((fails + 1)); continue; }
    cmp -s ref.out o.out || { echo "differs: oracle k=$k t=$t b=$b"; fails=$((fails + 1)); }
    $M tms-batch mini.tms reads.fa -o o.out --mode oracle --positions --oracle ref.orc --interleave $k --threads $t --block-reads $b 2> /dev/null
    cmp -s ref.pos o.out || { echo "differs: oracle --positions k=$k t=$t b=$b"; fails=$((fails + 1)); }
done; done; done
$M tms-batch mini.tms reads.fa -o ref.rc --mode dual --orient rc --oracle-write rc.orc 2> /dev/null
$M tms-batch mini.tms reads.fa -o o.out --mode oracle --orient rc --oracle rc.orc --threads 3 --block-reads 13 2> /dev/null
cmp -s ref.rc o.out || { echo "differs: oracle --orient rc"; fails=$((fails + 1)); }
# Refused: the oracle mode without an oracle, an oracle for other reads.
$M tms-batch mini.tms reads.fa -o o.out --mode oracle 2> /dev/null && { echo "oracle mode accepted no oracle"; fails=$((fails + 1)); }
$M tms-batch mini.tms m.fa -o o.out --mode oracle --oracle ref.orc 2> /dev/null && { echo "oracle accepted other reads"; fails=$((fails + 1)); }
echo "checked oracle"

if [ $fails -eq 0 ]; then echo "All thread checks PASSED"; else echo "$fails thread checks FAILED"; exit 1; fi
