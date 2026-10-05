#!/bin/bash
# Check index-rewrap: each legacy index (MSIX v4 or TMSX v8) built by
# OLD_MS, rewrapped by NEW_MS, must equal the index NEW_MS builds from the
# same input with the same options, byte for byte apart from the header's
# provenance text.  Each new index must also pass index-verify, and a query
# on the rewrapped index must print what it prints on the new one.
#
# Usage: test_rewrap.sh OLD_MS NEW_MS DATA_DIR [WORK_DIR]
set -euo pipefail
OLD=$1
NEW=$2
DATA=$3
W=${4:-$(mktemp -d)}
TSV="$DATA/minishred1_20_002_lcp.tsv"
mkdir -p "$W"

# Compare two ORBITIDX files with their provenance text blanked.
same() {
    python3 - "$1" "$2" <<'EOF'
import struct, sys
def norm(p):
    b = bytearray(open(p, 'rb').read())
    ns, ni = struct.unpack_from('<H', b, 40)[0], struct.unpack_from('<H', b, 44)[0]
    at = 128 + 80 * ns
    b[44:46] = b'\0\0'
    b[at:at + ni] = b'\0' * ni
    return bytes(b)
sys.exit(0 if norm(sys.argv[1]) == norm(sys.argv[2]) else 1)
EOF
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
for q in "--mode psi" "--mode phi" "--mode dual --positions" "--mode phi --report smem-all"; do
    # shellcheck disable=SC2086
    "$NEW" tms-batch "$W/tms_full.new" "$reads" -o "$W/q.new" $q > /dev/null 2>&1
    # shellcheck disable=SC2086
    "$NEW" tms-batch "$W/tms_full.rewrapped" "$reads" -o "$W/q.rw" $q > /dev/null 2>&1
    cmp -s "$W/q.new" "$W/q.rw" || { echo "tms-batch $q output differs"; exit 1; }
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
