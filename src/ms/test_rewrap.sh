#!/bin/bash
# Check index-rewrap: each legacy index (MSIX v4 or TMSX v8) built by
# OLD_MS, rewrapped by NEW_MS, must equal the index NEW_MS builds from the
# same input with the same options, byte for byte apart from the header's
# provenance text, with one exception: a legacy tms index with phi stores
# offsets in LF's phi_off column, which NEW_MS builds with width 0 (every
# phi point it reads has offset 0), so its LF section may differ in that
# column's width alone.  Each new index must also pass index-verify, and
# queries on the rewrapped index must print what they print on the new one.
#
# Usage: test_rewrap.sh OLD_MS NEW_MS DATA_DIR [WORK_DIR]
set -euo pipefail
OLD=$1
NEW=$2
DATA=$3
W=${4:-$(mktemp -d)}
TSV="$DATA/minishred1_20_002_lcp.tsv"
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
python3 - "$TSV" > "$reads" <<'EOF'
import random, sys
random.seed(7)
for _ in range(200):
    print(''.join(random.choice('ACGT') for _ in range(random.randint(20, 120))))
EOF

check() {  # name, build command, options...
    local name=$1 cmd=$2
    shift 2
    "$OLD" "$cmd" "$TSV" "$W/$name.legacy" "$@" > /dev/null 2>&1
    "$NEW" "$cmd" "$TSV" "$W/$name.new" "$@" > /dev/null 2>&1
    "$NEW" index-rewrap "$W/$name.legacy" "$W/$name.rewrapped" > /dev/null
    same "$W/$name.new" "$W/$name.rewrapped" || { echo "$name: rewrapped index differs from a new build"; exit 1; }
    "$NEW" index-verify "$W/$name.rewrapped" > /dev/null || { echo "$name: index-verify failed"; exit 1; }
    echo "$name: rewrap equals a new build"
}

check ms build
check ms_split_bits build --spill-split-bits 2
check ms_coalesce build --coalesce-spillover --spill-align 4
check ms_minima build --minima-only --lf-split 4
check tms_full tms-build-tsv --phi
check tms_full_noinv tms-build-tsv --phi --no-phi-inv
check tms_psi tms-build-tsv --layout psi
check tms_phi tms-build-tsv --phi --layout phi
check tms_phi_unsplit tms-build-tsv --phi --layout phi --phi-split 0 --lf-split 4

# Queries on rewrapped indexes print what they print on new builds.
"$NEW" batch "$W/ms.new" "$reads" -o "$W/q.new" > /dev/null 2>&1
"$NEW" batch "$W/ms.rewrapped" "$reads" -o "$W/q.rw" > /dev/null 2>&1
cmp -s "$W/q.new" "$W/q.rw" || { echo "ms batch output differs"; exit 1; }
for name in tms_full tms_full_noinv tms_phi tms_phi_unsplit; do
    for q in "--mode psi" "--mode phi" "--mode dual --positions" "--mode phi --report smem-all"; do
        # The phi layout serves phi mode alone, and smem-all needs phi_inv.
        case "$name:$q" in
            tms_phi*:"--mode psi" | tms_phi*:"--mode dual --positions" | tms_full_noinv:*smem-all) continue ;;
        esac
        # shellcheck disable=SC2086
        "$NEW" tms-batch "$W/$name.new" "$reads" -o "$W/q.new" $q > /dev/null 2>&1
        # shellcheck disable=SC2086
        "$NEW" tms-batch "$W/$name.rewrapped" "$reads" -o "$W/q.rw" $q > /dev/null 2>&1
        cmp -s "$W/q.new" "$W/q.rw" || { echo "$name: tms-batch $q output differs"; exit 1; }
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
