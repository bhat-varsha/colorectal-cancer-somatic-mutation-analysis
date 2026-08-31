#!/usr/bin/env bash

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

INPUT="results/functional_mutations.maf.gz"
OUTPUT="results/mutation_spectrum.tsv"

echo "============================================================"
echo "10 — MUTATION SPECTRUM ANALYSIS"
echo "============================================================"

echo "[1/2] Calculating SNV mutation spectrum..."

zcat "$INPUT" |
awk -F '\t' '
BEGIN {
    OFS="\t"
}

NR==1 {
    for(i=1;i<=NF;i++) {
        if($i=="Reference_Allele") ref=i
        if($i=="Tumor_Seq_Allele2") alt=i
        if($i=="Variant_Type") type=i
    }
    next
}

toupper($type)=="SNP" {

    r=toupper($ref)
    a=toupper($alt)

    # Only retain simple single-base substitutions
    if(length(r)==1 && length(a)==1 &&
       (r=="A" || r=="C" || r=="G" || r=="T") &&
       (a=="A" || a=="C" || a=="G" || a=="T") &&
       r != a) {

        change=r ">" a
        count[change]++
        total++
    }
}

END {
    print "Mutation_Type", "Count"

    order[1]="A>G"
    order[2]="G>A"
    order[3]="C>T"
    order[4]="T>C"
    order[5]="A>C"
    order[6]="A>T"
    order[7]="C>A"
    order[8]="C>G"
    order[9]="G>C"
    order[10]="G>T"
    order[11]="T>A"
    order[12]="T>G"

    for(i=1;i<=12;i++) {
        x=order[i]
        print x, count[x]+0
    }

    print "TOTAL_SNV", total
}
' > "$OUTPUT"

echo "Created: $OUTPUT"

echo
echo "[2/2] Calculating transition/transversion ratio..."

awk -F '\t' '
NR==1 {
    next
}

$1=="A>G" || $1=="G>A" ||
$1=="C>T" || $1=="T>C" {
    transitions += $2
}

$1=="A>C" || $1=="A>T" ||
$1=="C>A" || $1=="C>G" ||
$1=="G>C" || $1=="G>T" ||
$1=="T>A" || $1=="T>G" {
    transversions += $2
}

END {
    print "Transitions:", transitions
    print "Transversions:", transversions

    if(transversions > 0)
        print "Ti/Tv ratio:", transitions / transversions
    else
        print "Ti/Tv ratio: undefined"
}
' "$OUTPUT"

echo
echo "============================================================"
echo "10 — MUTATION SPECTRUM COMPLETE"
echo "============================================================"
