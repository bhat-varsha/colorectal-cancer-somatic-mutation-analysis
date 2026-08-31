#!/usr/bin/env bash

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

MAF="data/raw/cohortMAF.2026-08-15.maf.gz"
OUT="results/functional_mutations.maf.gz"

echo "============================================================"
echo "02 — FUNCTIONAL MUTATION FILTERING"
echo "============================================================"

zcat "$MAF" |
awk -F '\t' '
NR==1 {
    for(i=1;i<=NF;i++) {
        if($i=="Variant_Classification") classification=i
    }

    print
    next
}

$classification=="Missense_Mutation" ||
$classification=="Nonsense_Mutation" ||
$classification=="Frame_Shift_Del" ||
$classification=="Frame_Shift_Ins" ||
$classification=="Splice_Site" ||
$classification=="In_Frame_Del" ||
$classification=="In_Frame_Ins" ||
$classification=="Translation_Start_Site" ||
$classification=="Nonstop_Mutation" {
    print
}' |
gzip > "$OUT"

echo
echo "Output:"
echo "$OUT"

echo
echo "Functional mutation records:"
zcat "$OUT" | tail -n +2 | wc -l

echo
echo "Functional mutation classifications:"

zcat "$OUT" |
awk -F '\t' '
NR==1 {
    for(i=1;i<=NF;i++)
        if($i=="Variant_Classification") classification=i
    next
}
{
    count[$classification]++
}
END {
    for(x in count)
        print x "\t" count[x]
}' |
sort -k2,2nr

echo
echo "Functional variant types:"

zcat "$OUT" |
awk -F '\t' '
NR==1 {
    for(i=1;i<=NF;i++)
        if($i=="Variant_Type") type=i
    next
}
{
    count[$type]++
}
END {
    for(x in count)
        print x "\t" count[x]
}' |
sort -k2,2nr

echo
echo "============================================================"
echo "02 — FUNCTIONAL FILTER COMPLETE"
echo "============================================================"
