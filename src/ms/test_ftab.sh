#!/bin/bash
# Checks stored ftabs end to end.  With DATA_DIR, builds ms and tms (full,
# psi and phi layout) indexes from DATA_DIR/minishred1_20_002_lcp.tsv and
# reads drawn from its text; with --indexes, takes existing ORBITIDX
# indexes and a reads file instead.  For each index it runs ftab-build
# (--auto and a given k), ftab-info and ftab-verify, and checks that batch
# or tms-batch output with --ftab-file is byte-identical to the output
# without it, in every mode the layout serves, with and without --positions
# and with each SMEM report the layout serves.  ftab-verify and --ftab-file
# must refuse a table built from another index.
# Usage: test_ftab.sh MS_BINARY DATA_DIR [WORK_DIR]
#        test_ftab.sh MS_BINARY --indexes READS WORK_DIR INDEX...
set -u
M=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
shift
fails=0
fail() { echo "FAILED: $*"; fails=$((fails + 1)); }

if [ "$1" = "--indexes" ]; then
    READS=$(cd "$(dirname "$2")" && pwd)/$(basename "$2")
    W=$3
    shift 3
    INDEXES=("$@")
    mkdir -p "$W"
    THREADS=${THREADS:-8}
else
    DATA=$1
    W=${2:-$(mktemp -d)}
    mkdir -p "$W"
    TSV="$DATA/minishred1_20_002_lcp.tsv"
    $M build "$TSV" "$W/mini.idx" > /dev/null 2>&1 || { echo "build failed"; exit 1; }
    $M build "$TSV" "$W/mini-split.idx" --lf-split 2 > /dev/null 2>&1 || { echo "build failed"; exit 1; }
    for l in full psi phi; do
        phi=--phi
        [ $l = psi ] && phi=
        $M tms-build-tsv "$TSV" "$W/mini-$l.tms" $phi --layout $l > /dev/null 2>&1 || { echo "tms-build failed"; exit 1; }
    done
    $M tms-build-tsv "$TSV" "$W/mini-full-split.tms" --phi --lf-split 2 > /dev/null 2>&1 || { echo "tms-build failed"; exit 1; }
    $M tms-text "$W/mini-full.tms" 2> /dev/null > "$W/text.txt"
    python3 - "$W" <<'PY'
import random, sys
d = sys.argv[1]
t = open(d + '/text.txt').read().strip().replace('%', '').replace('$', '')
random.seed(11)
with open(d + '/reads.fa', 'w') as out:
    for i in range(3000):
        L = random.choice([0, 1, 5, 9, 10, 11, 14, 50, 150, 150, 400])
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
    READS="$W/reads.fa"
    INDEXES=("$W/mini.idx" "$W/mini-full.tms" "$W/mini-psi.tms" "$W/mini-phi.tms")
    THREADS=${THREADS:-3}
    # Tables of other indexes are refused, naming the field.
    $M ftab-build "$W/mini-split.idx" "$W/mini-split.ftab" --k 6 > /dev/null || fail "ftab-build mini-split.idx"
    $M ftab-build "$W/mini-full-split.tms" "$W/mini-full-split.ftab" --k 6 > /dev/null || fail "ftab-build mini-full-split.tms"
fi

tool_of() { "$M" index-plan "$1" 2> /dev/null | grep -q 'tool: tms' && echo tms || echo ms; }

cd "$W"
for idx in "${INDEXES[@]}"; do
    base=$(basename "$idx")
    tool=$(tool_of "$idx")
    for kopt in "--auto" "--k 8"; do
        ft="$W/$base.$(echo $kopt | tr -d ' -').ftab"
        $M ftab-build "$idx" "$ft" $kopt --threads $THREADS > build.txt 2>&1 || { fail "ftab-build $base $kopt"; cat build.txt; continue; }
        echo "$base $kopt: $(grep -o 'k=[0-9]*, entry_bits=[0-9]*, bytes=[0-9]*, occupancy=[0-9.e-]*' build.txt)"
        $M ftab-info "$ft" > info.txt || fail "ftab-info $base"
        $M ftab-verify "$idx" "$ft" --crc > /dev/null || fail "ftab-verify $base $kopt"
        if [ "$tool" = ms ]; then
            cmds=("batch|")
            [ -f "$W/mini-split.ftab" ] && { $M ftab-verify "$idx" "$W/mini-split.ftab" 2> /dev/null && fail "ftab-verify accepted another index's table"; }
        else
            layout=$(grep -o 'layout=[a-z]*' info.txt | cut -d= -f2)
            cmds=()
            modes=()
            case $layout in
                full) modes=(psi phi phiskip dual) ;;
                psi) modes=(psi) ;;
                phi) modes=(phi) ;;
            esac
            for m in "${modes[@]}"; do
                cmds+=("tms-batch|--mode $m")
                if [ "$layout" != psi ]; then
                    cmds+=("tms-batch|--mode $m --positions" "tms-batch|--mode $m --report smem-one" "tms-batch|--mode $m --report smem-all --max-smem-positions 20")
                    cmds+=("tms-batch|--mode $m --phi-starts")
                fi
            done
            for other in "$W/mini-full-split.ftab" "$W/mini.idx.auto.ftab"; do
                [ -f "$other" ] || continue
                $M ftab-verify "$idx" "$other" 2> err.txt && fail "ftab-verify accepted $(basename $other) for $base"
                $M tms-batch "$idx" "$READS" --ftab-file "$other" --no-output 2> /dev/null && fail "tms-batch accepted $(basename $other) for $base"
            done
        fi
        for c in "${cmds[@]}"; do
            cmd=${c%%|*}
            opts=${c#*|}
            [ -f "ref.$base.$cmd.$(echo $opts | tr -d ' -')" ] || $M $cmd "$idx" "$READS" $opts --threads $THREADS -o "ref.$base.$cmd.$(echo $opts | tr -d ' -')" 2> err.txt \
                || { fail "$cmd $base $opts"; cat err.txt; continue; }
            for k in 0 32; do
                $M $cmd "$idx" "$READS" $opts --threads $THREADS --interleave $k --ftab-file "$ft" -o o.out 2> err.txt \
                    || { fail "$cmd $base $opts --ftab-file $kopt k=$k"; cat err.txt; continue; }
                cmp -s "ref.$base.$cmd.$(echo $opts | tr -d ' -')" o.out || fail "differs: $cmd $base $opts --ftab-file ($kopt) k=$k"
            done
            echo "  checked $cmd $opts"
        done
    done
done
if [ $fails -eq 0 ]; then echo "test_ftab: all checks PASSED"; else echo "test_ftab: $fails FAILED"; exit 1; fi
