#!/usr/bin/env bash

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

INPUT="results/functional_mutations.maf.gz"

TOTAL_CASES=583

OUTPUT="results/gene_case_frequency.tsv"
OUTPUT_PERCENT="results/gene_case_frequency_percent.tsv"

echo "============================================================"
echo "05 — GENE-LEVEL MUTATION FREQUENCY"
echo "============================================================"

# ------------------------------------------------------------
# 1. Count unique cases mutated in each gene
# ------------------------------------------------------------

echo "[1/2] Calculating gene-level mutation frequency..."

zcat "$INPUT" |
awk -F '\t' '
NR==1 {
    for(i=1;i<=NF;i++) {
        if($i=="Hugo_Symbol") gene=i
        if($i=="case_id") caseid=i
    }
    next
}

$gene!="" && $caseid!="" {
    key=$gene "\t" $caseid
    seen[key]=1
}

END {
    for(key in seen) {
        split(key,a,"\t")
        gene_cases[a[1]]++
    }

    for(g in gene_cases)
        print g "\t" gene_cases[g]
}' |
sort -k2,2nr > "$OUTPUT"

echo "Created: $OUTPUT"

# ------------------------------------------------------------
# 2. Calculate percentage of cohort cases
# ------------------------------------------------------------

echo
echo "[2/2] Calculating mutation prevalence percentages..."

awk -v total="$TOTAL_CASES" '{
    printf "%s\t%d\t%.2f%%\n", $1, $2, ($2/total)*100
}' "$OUTPUT" > "$OUTPUT_PERCENT"

echo "Created: $OUTPUT_PERCENT"

echo
echo "Top 20 mutated genes:"
echo "------------------------------------------------------------"

head -20 "$OUTPUT"

echo
echo "============================================================"
echo "05 — GENE-LEVEL FREQUENCY COMPLETE"
echo "============================================================"
