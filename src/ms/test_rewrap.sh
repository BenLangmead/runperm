#!/bin/bash
# Check index-rewrap: each legacy index (MSIX v4 or TMSX v8) built by
# OLD_MS, rewrapped by NEW_MS, must equal the index NEW_MS builds from the
# same input with the same options, byte for byte apart from the header's
# provenance text, with one exception: a legacy tms index with phi stores
# offsets in LF's phi_off column, which NEW_MS builds with width 0 (every
# phi point it reads has offset 0), so its LF section may differ in that
# column's width alone.  Each new index must also pass index-verify, and
# queries on the rewrapped index must print what they print on the new one.
# An ftab of either holds the same table, and that of the new index links to
# the rewrapped one where their sections are equal.
#
# Usage: test_rewrap.sh OLD_MS NEW_MS DATA_DIR [WORK_DIR]
set -euo pipefail
OLD=$1
NEW=$2
DATA=$3
W=${4:-$(mktemp -d)}
H="$DATA/minishred1_20_002.bwt.heads"
L="$DATA/minishred1_20_002.bwt.len"
LCP="$DATA/minishred1_20_002.lcp"
mkdir -p "$W"

# Compare a new ORBITIDX file ($1) with a rewrapped one ($2): equal fixed
# fields, and per section equal entries and payloads, except that LF's
# phi_off column (code 19) may have width 0 in the new file and more in the
# rewrapped one, when LF's entry and payload are not compared beyond rows,
# domain and the other columns.
same() {
    python3 - "$1" "$2" <<'PYEOF'
import struct, sys
def parse(p):
    b = open(p, 'rb').read()
    ns = struct.unpack_from('<H', b, 40)[0]
    secs = []
    for i in range(ns):
        e = b[128 + 80 * i:128 + 80 * (i + 1)]
        sid, kind, nc = struct.unpack_from('<HHB', e)
        off, size, rows, dom = struct.unpack_from('<QQQQ', e, 8)
        cols = [tuple(e[48 + 2 * c:50 + 2 * c]) for c in range(nc)]
        secs.append(dict(id=sid, kind=kind, rows=rows, dom=dom, cols=cols, payload=b[off:off + size]))
    return b[:40] + b[46:128], secs
fa, sa = parse(sys.argv[1])
fb, sb = parse(sys.argv[2])
ok = fa == fb and [x['id'] for x in sa] == [x['id'] for x in sb]
for x, y in zip(sa, sb):
    if x == y:
        continue
    ok = ok and (x['id'] == 1 and x['rows'] == y['rows'] and x['dom'] == y['dom'] and
                 [c for c in x['cols'] if c[0] != 19] == [c for c in y['cols'] if c[0] != 19] and
                 (19, 0) in x['cols'] and (19, 0) not in y['cols'])
sys.exit(0 if ok else 1)
PYEOF
}

reads="$W/reads.txt"
python3 - > "$reads" <<'EOF'
import random, sys
random.seed(7)
for _ in range(200):
    print(''.join(random.choice('ACGT') for _ in range(random.randint(20, 120))))
EOF

check() {  # name, build command, options...
    local name=$1 cmd=$2
    shift 2
    "$OLD" "$cmd" "$H" "$L" "$W/$name.legacy" "$@" > /dev/null 2>&1
    "$NEW" "$cmd" "$H" "$L" "$W/$name.new" "$@" > /dev/null 2>&1
    "$NEW" index-rewrap "$W/$name.legacy" "$W/$name.rewrapped" > /dev/null
    same "$W/$name.new" "$W/$name.rewrapped" || { echo "$name: rewrapped index differs from a new build"; exit 1; }
    "$NEW" index-verify "$W/$name.rewrapped" > /dev/null || { echo "$name: index-verify failed"; exit 1; }
    # ftabs of the new build and of the rewrap hold the same table.  The
    # new build's links to the rewrap, whose info text differs, unless the
    # rewrap's LF keeps a phi_off column, which makes its LF another section.
    "$NEW" ftab-build "$W/$name.new" "$W/$name.ftab" --kmer-len 4 > /dev/null || { echo "$name: ftab-build failed"; exit 1; }
    "$NEW" ftab-build "$W/$name.rewrapped" "$W/$name.rw.ftab" --kmer-len 4 > /dev/null || { echo "$name: ftab-build failed"; exit 1; }
    cmp -s -i 4096 "$W/$name.ftab" "$W/$name.rw.ftab" || { echo "$name: the rewrap's ftab holds another table"; exit 1; }
    case "$name" in
        tms_full* | tms_phi*)
            if "$NEW" ftab-verify "$W/$name.rewrapped" "$W/$name.ftab" 2> "$W/err" > /dev/null || ! grep -q "section LF" "$W/err"; then
                echo "$name: ftab-verify must refuse a rewrap whose LF keeps phi_off"; exit 1
            fi ;;
        *)
            "$NEW" ftab-verify "$W/$name.rewrapped" "$W/$name.ftab" > /dev/null || { echo "$name: ftab-verify refused the rewrap"; exit 1; } ;;
    esac
    echo "$name: rewrap equals a new build"
}

check ms build-rlbwt --lcp-bin "$LCP"
check ms_split_bits build-rlbwt --lcp-bin "$LCP" --spill-split-bits 2
check ms_coalesce build-rlbwt --lcp-bin "$LCP" --coalesce-spillover --spill-align 4
check ms_minima build-rlbwt --lcp-bin "$LCP" --minima-only --lf-split 4
check tms_full tms-build --lcp-bin "$LCP"
check tms_full_noinv tms-build --lcp-bin "$LCP" --no-phi-inv
check tms_psi tms-build --layout psi
check tms_phi tms-build --lcp-bin "$LCP" --layout phi
check tms_phi_unsplit tms-build --lcp-bin "$LCP" --layout phi --phi-split 0 --lf-split 4

# Queries on rewrapped indexes print what they print on new builds.
"$NEW" query "$W/ms.new" "$reads" -o "$W/q.new" > /dev/null 2>&1
"$NEW" query "$W/ms.rewrapped" "$reads" -o "$W/q.rw" > /dev/null 2>&1
cmp -s "$W/q.new" "$W/q.rw" || { echo "ms query output differs"; exit 1; }
for name in tms_full tms_full_noinv tms_phi tms_phi_unsplit; do
    for q in "--mode psi" "--mode phi" "--mode dual --positions" "--mode phi --report smem-all"; do
        # The phi layout serves phi mode alone, and smem-all needs phi_inv.
        case "$name:$q" in
            tms_phi*:"--mode psi" | tms_phi*:"--mode dual --positions" | tms_full_noinv:*smem-all) continue ;;
        esac
        # shellcheck disable=SC2086
        "$NEW" query "$W/$name.new" "$reads" -o "$W/q.new" $q > /dev/null 2>&1
        # shellcheck disable=SC2086
        "$NEW" query "$W/$name.rewrapped" "$reads" -o "$W/q.rw" $q > /dev/null 2>&1
        cmp -s "$W/q.new" "$W/q.rw" || { echo "$name: tms query $q output differs"; exit 1; }
    done
done
echo "queries agree"

# A damaged section fails index-verify, and a truncated file fails to load.
cp "$W/tms_full.new" "$W/damaged"
printf '\xff' | dd of="$W/damaged" bs=1 seek=5000 conv=notrunc 2> /dev/null
if "$NEW" index-verify "$W/damaged" > /dev/null 2>&1; then echo "index-verify missed a damaged byte"; exit 1; fi
head -c 300000 "$W/tms_full.new" > "$W/truncated"
if "$NEW" tms-inspect "$W/truncated" > /dev/null 2>&1; then echo "a truncated index loaded"; exit 1; fi
echo "damage detected"
echo "ALL REWRAP TESTS PASSED"
