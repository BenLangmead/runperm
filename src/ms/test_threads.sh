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
CMDS
if [ $fails -eq 0 ]; then echo "All thread checks PASSED"; else echo "$fails thread checks FAILED"; exit 1; fi
