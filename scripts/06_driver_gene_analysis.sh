#!/usr/bin/env bash

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

INPUT="results/functional_mutations.maf.gz"
OUTPUT="results/driver_gene_functional_mutations.tsv"

echo "============================================================"
echo "06 — DRIVER GENE ANALYSIS"
echo "============================================================"

# Selected colorectal cancer driver genes
DRIVER_GENES="APC|TP53|KRAS|PIK3CA|FBXW7|SMAD4|BRAF|NRAS"

echo "[1/2] Extracting functional mutations in selected driver genes..."

zcat "$INPUT" |
awk -F '\t' -v genes="$DRIVER_GENES" '
BEGIN {
    OFS="\t"
    split(genes, g, "|")
    for (i in g)
        driver[g[i]]=1
}

NR==1 {
    for(i=1;i<=NF;i++) {
        if($i=="Hugo_Symbol") gene=i
    }
    print
    next
}

$gene in driver {
    print
}
' > "$OUTPUT"

echo "Created: $OUTPUT"

echo
echo "[2/2] Summary of selected driver genes..."

awk -F '\t' '
NR==1 {
    for(i=1;i<=NF;i++) {
        if($i=="Hugo_Symbol") gene=i
    }
    next
}
{
    count[$gene]++
}
END {
    for(g in count)
        print g "\t" count[g]
}
' "$OUTPUT" | sort -k2,2nr

echo
echo "============================================================"
echo "06 — DRIVER GENE ANALYSIS COMPLETE"
echo "============================================================"
