#!/bin/bash

set -eox pipefail

BASE=data/minishred1_20_002_lcp
CONFIGS=()

NM=coal_split
CONFIGS+=("$NM")
./src/ms/ms build \
    ${BASE}.tsv \
    ${BASE}.tsv.${NM}.idx \
    --coalesce-spillover --split-threshold 1
./src/ms/ms inspect ${BASE}.tsv.${NM}.idx | tee ${BASE}.tsv.${NM}.idx.inspect.txt

for ALIGN in 2 3 4 5; do
    NM=align${ALIGN}_coal_split
    CONFIGS+=("$NM")
    ./src/ms/ms build \
        ${BASE}.tsv \
        ${BASE}.tsv.${NM}.idx \
        --coalesce-spillover --split-threshold 1 --spill-align ${ALIGN}
    ./src/ms/ms inspect ${BASE}.tsv.${NM}.idx | tee ${BASE}.tsv.${NM}.idx.inspect.txt
done

for MULTI in 2 3 4 5 6 7 8 9 10; do
    NM=multispill${MULTI}_coal_split
    CONFIGS+=("$NM")
    ./src/ms/ms build \
        ${BASE}.tsv \
        ${BASE}.tsv.${NM}.idx \
        --coalesce-spillover --split-threshold 1 --spill-split-bits ${MULTI}
    ./src/ms/ms inspect ${BASE}.tsv.${NM}.idx | tee ${BASE}.tsv.${NM}.idx.inspect.txt
done

NM=coal
CONFIGS+=("$NM")
./src/ms/ms build \
    ${BASE}.tsv \
    ${BASE}.tsv.${NM}.idx \
    --coalesce-spillover
./src/ms/ms inspect ${BASE}.tsv.${NM}.idx | tee ${BASE}.tsv.${NM}.idx.inspect.txt

NM=split
CONFIGS+=("$NM")
./src/ms/ms build \
    ${BASE}.tsv \
    ${BASE}.tsv.${NM}.idx \
    --split-threshold 1
./src/ms/ms inspect ${BASE}.tsv.${NM}.idx | tee ${BASE}.tsv.${NM}.idx.inspect.txt

NM=base
CONFIGS+=("$NM")
./src/ms/ms build \
    ${BASE}.tsv \
    ${BASE}.tsv.${NM}.idx
./src/ms/ms inspect ${BASE}.tsv.${NM}.idx | tee ${BASE}.tsv.${NM}.idx.inspect.txt

# Print summary table
set +x
echo ""
echo "=== Summary Table ==="

# Compute reference total from "base" config
REF_FILE=${BASE}.tsv.base.idx.inspect.txt
if [ -f "$REF_FILE" ]; then
    ref_base=$(awk -F: '/base RLBWT bits:/{gsub(/ /,"",$2); print $2}' "$REF_FILE")
    ref_lcp=$(awk -F: '/LCP columns bits:/{gsub(/ /,"",$2); print $2}' "$REF_FILE")
    ref_spill=$(awk -F: '/spillover vector bits:/{gsub(/ /,"",$2); print $2}' "$REF_FILE")
    REF_TOTAL=$((ref_base + ref_lcp + ref_spill))
else
    REF_TOTAL=0
fi

printf "%-25s  %10s  %11s  %10s  %10s  %10s  %20s  %8s\n" \
    "Config" "Base bits" "Bits per" "LCP bits" "Bits per" "Spill bits" "Total bits" "Pct LCP"
printf "%-25s  %10s  %11s  %10s  %10s  %10s  %20s  %8s\n" \
    "-------------------------" "----------" "-----------" "----------" "----------" "----------" "--------------------" "--------"

for NM in "${CONFIGS[@]}"; do
    FILE=${BASE}.tsv.${NM}.idx.inspect.txt
    [ -f "$FILE" ] || continue

    base_tot=$(awk -F: '/base RLBWT bits:/{gsub(/ /,"",$2); print $2}' "$FILE")
    len_b=$(awk '/length:/{print $2}' "$FILE")
    ptr_b=$(awk '/pointer:/{print $2}' "$FILE")
    off_b=$(awk '/offset:/{print $2}' "$FILE")
    chr_b=$(awk '/character:/{print $2}' "$FILE")

    lcp_tot=$(awk -F: '/LCP columns bits:/{gsub(/ /,"",$2); print $2}' "$FILE")
    lcp_t=$(awk '/lcp_top:/{print $2}' "$FILE")
    lcp_m=$(awk '/lcp_min_sub:/{print $2}' "$FILE")
    lcp_s=$(awk '/lcp_spill:/{print $2}' "$FILE")

    spill=$(awk -F: '/spillover vector bits:/{gsub(/ /,"",$2); print $2}' "$FILE")

    total=$((base_tot + lcp_tot + spill))
    if [ "$NM" = "base" ] || [ "$REF_TOTAL" -eq 0 ]; then
        total_str="$total (ref)"
    else
        pct=$(awk "BEGIN {printf \"%.1f\", ($total - $REF_TOTAL) * 100.0 / $REF_TOTAL}")
        total_str="$total (${pct}%)"
    fi

    pct_lcp=$(awk "BEGIN {printf \"%.1f%%\", ($lcp_tot + $spill) * 100.0 / $total}")

    printf "%-25s  %10s  %11s  %10s  %10s  %10s  %20s  %8s\n" \
        "$NM" "$base_tot" "${len_b}+${ptr_b}+${off_b}+${chr_b}" \
        "$lcp_tot" "${lcp_t}+${lcp_m}+${lcp_s}" "$spill" "$total_str" "$pct_lcp"
done
