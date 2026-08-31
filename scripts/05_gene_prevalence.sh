#!/usr/bin/env bash

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

INPUT="results/case_gene_project.tsv"
OUTPUT="results/gene_prevalence_coad_read.tsv"

COAD_TOTAL=428
READ_TOTAL=155

echo "============================================================"
echo "05 — GENE MUTATION PREVALENCE"
echo "============================================================"

echo "[1/2] Calculating gene prevalence in COAD and READ..."

awk -F '\t' -v coad_total="$COAD_TOTAL" -v read_total="$READ_TOTAL" '
BEGIN {
    OFS="\t"
}

{
    caseid=$1
    gene=$2
    project=$3

    # Count each case-gene combination only once
    key=caseid SUBSEP gene

    if (!(key in seen)) {
        seen[key]=1

        if (project=="TCGA-COAD")
            coad[gene]++

        else if (project=="TCGA-READ")
            read[gene]++
    }
}

END {
    print "Gene","COAD_cases","COAD_percent","READ_cases","READ_percent"

    for (gene in coad) {
        c=coad[gene]
        r=read[gene]+0

        print gene,c,(c/coad_total)*100,r,(r/read_total)*100
    }

    for (gene in read) {
        if (!(gene in coad)) {
            c=0
            r=read[gene]

            print gene,c,(c/coad_total)*100,r,(r/read_total)*100
        }
    }
}
' "$INPUT" |
{
    head -n 1
    tail -n +2 | sort -k2,2nr
} > "$OUTPUT"

echo "Created: $OUTPUT"

echo
echo "[2/2] Checking final table..."

echo
echo "Top 20 genes:"
head -20 "$OUTPUT" | column -t

echo
echo "Total genes:"
tail -n +2 "$OUTPUT" | wc -l

echo
echo "============================================================"
echo "05 — GENE PREVALENCE COMPLETE"
echo "============================================================"
