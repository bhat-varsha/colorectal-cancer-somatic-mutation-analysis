#!/usr/bin/env bash

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

INPUT="results/functional_mutations.maf.gz"
MAPPING="results/case_project_mapping.tsv"
OUTPUT="results/mutation_burden_per_case.tsv"

echo "============================================================"
echo "03 — MUTATION BURDEN ANALYSIS"
echo "============================================================"

# ------------------------------------------------------------
# Calculate functional mutation burden per case
# ------------------------------------------------------------

echo "[1/3] Calculating mutation burden per case..."

zcat "$INPUT" |
awk -F '\t' '
NR==1 {
    for(i=1;i<=NF;i++)
        if($i=="case_id") caseid=i
    next
}

$caseid != "" {
    count[$caseid]++
}

END {
    for(x in count)
        print x "\t" count[x]
}' |
sort -k2,2nr > "$OUTPUT"

echo "Created: $OUTPUT"

# ------------------------------------------------------------
# Overall mutation burden summary
# ------------------------------------------------------------

echo
echo "[2/3] Overall mutation burden summary..."

awk '
{
    values[NR]=$2
    sum+=$2
}
END {
    n=NR

    asort(values)

    if(n%2==1)
        median=values[(n+1)/2]
    else
        median=(values[n/2]+values[n/2+1])/2

    print "Cases:", n
    print "Minimum:", values[1]
    print "Maximum:", values[n]
    print "Mean:", sum/n
    print "Median:", median
}' "$OUTPUT"

# ------------------------------------------------------------
# COAD vs READ summary
# ------------------------------------------------------------

echo
echo "[3/3] COAD vs READ mutation burden summary..."

awk -F '\t' '
BEGIN {
    while((getline < "'"$MAPPING"'") > 0)
        project[$1]=$2
}

{
    caseid=$1
    burden=$2
    p=project[caseid]

    if(p=="TCGA-COAD")
        coad[++n_coad]=burden

    else if(p=="TCGA-READ")
        read[++n_read]=burden
}

END {

    # ---------------- COAD ----------------

    asort(coad)

    sum=0
    for(i=1;i<=n_coad;i++)
        sum+=coad[i]

    print "TCGA-COAD"
    print "Cases:", n_coad
    print "Minimum:", coad[1]
    print "Maximum:", coad[n_coad]
    print "Mean:", sum/n_coad

    if(n_coad%2)
        print "Median:", coad[(n_coad+1)/2]
    else
        print "Median:", (coad[n_coad/2]+coad[n_coad/2+1])/2

    print ""

    # ---------------- READ ----------------

    asort(read)

    sum=0
    for(i=1;i<=n_read;i++)
        sum+=read[i]

    print "TCGA-READ"
    print "Cases:", n_read
    print "Minimum:", read[1]
    print "Maximum:", read[n_read]
    print "Mean:", sum/n_read

    if(n_read%2)
        print "Median:", read[(n_read+1)/2]
    else
        print "Median:", (read[n_read/2]+read[n_read/2+1])/2
}' "$OUTPUT"

echo
echo "============================================================"
echo "03 — MUTATION BURDEN COMPLETE"
echo "============================================================"
