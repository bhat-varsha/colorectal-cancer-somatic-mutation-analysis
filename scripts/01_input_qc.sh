#!/usr/bin/env bash

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

MAF="data/raw/cohortMAF.2026-08-15.maf.gz"
SAMPLE="data/raw/sample.tsv"

mkdir -p results

echo "============================================================"
echo "01 — INPUT QC AND COHORT VALIDATION"
echo "============================================================"

# ------------------------------------------------------------
# Check input files
# ------------------------------------------------------------

if [[ ! -f "$MAF" ]]; then
    echo "ERROR: MAF file not found: $MAF"
    exit 1
fi

if [[ ! -f "$SAMPLE" ]]; then
    echo "ERROR: sample.tsv not found: $SAMPLE"
    exit 1
fi

echo "MAF:    $MAF"
echo "Sample: $SAMPLE"
echo

# ------------------------------------------------------------
# Raw MAF overview
# ------------------------------------------------------------

echo "Raw MAF summary:"
echo

zcat "$MAF" | awk -F '\t' '
NR==1 {
    for(i=1;i<=NF;i++) {
        if($i=="Hugo_Symbol") gene=i
        if($i=="Variant_Classification") classification=i
        if($i=="Variant_Type") type=i
        if($i=="case_id") caseid=i
    }
    next
}
{
    mutations++
    genes[$gene]=1
    cases[$caseid]=1
    classifications[$classification]++
    types[$type]++
}
END {
    print "Total mutation records:", mutations
    print "Unique genes:", length(genes)
    print "Unique cases:", length(cases)

    print ""
    print "Variant Types:"
    for(x in types)
        print x, types[x]

    print ""
    print "Variant Classifications:"
    for(x in classifications)
        print x, classifications[x]
}'

echo

# ------------------------------------------------------------
# MAF case IDs
# ------------------------------------------------------------

echo "Creating MAF case ID list..."

zcat "$MAF" |
awk -F '\t' '
NR==1 {
    for(i=1;i<=NF;i++)
        if($i=="case_id") caseid=i
    next
}
$caseid != "" {
    print $caseid
}' |
sort -u > results/maf_case_ids.txt

# ------------------------------------------------------------
# Clinical cohort case IDs
# ------------------------------------------------------------

echo "Creating clinical cohort case ID list..."

awk -F '\t' '
NR==1 {
    for(i=1;i<=NF;i++)
        if($i=="cases.case_id") caseid=i
    next
}
$caseid != "" {
    print $caseid
}' "$SAMPLE" |
sort -u > results/cohort_case_ids.txt

# ------------------------------------------------------------
# Compare clinical cohort vs MAF
# ------------------------------------------------------------

comm -23 \
    results/cohort_case_ids.txt \
    results/maf_case_ids.txt \
    > results/cases_without_mutations.txt

echo
echo "Cohort / MAF reconciliation:"
echo "Clinical cohort cases: $(wc -l < results/cohort_case_ids.txt)"
echo "MAF cases:             $(wc -l < results/maf_case_ids.txt)"
echo "Cases absent from MAF: $(wc -l < results/cases_without_mutations.txt)"

# ------------------------------------------------------------
# Case -> project mapping
# ------------------------------------------------------------

echo
echo "Creating case -> project mapping..."

awk -F '\t' '
NR==1 {
    for(i=1;i<=NF;i++) {
        if($i=="cases.case_id") caseid=i
        if($i=="project.project_id") project=i
    }
    next
}
$caseid != "" && $project != "" {
    print $caseid "\t" $project
}' "$SAMPLE" |
sort -u > results/case_project_mapping.tsv

cp results/case_project_mapping.tsv results/case_project_lookup.tsv

# ------------------------------------------------------------
# Cohort counts
# ------------------------------------------------------------

echo
echo "Cohort project counts:"

awk -F '\t' '
{
    count[$2]++
}
END {
    for(x in count)
        print x "\t" count[x]
}' results/case_project_mapping.tsv | sort

# ------------------------------------------------------------
# Sample type distribution
# ------------------------------------------------------------

echo
echo "Sample type distribution:"

awk -F '\t' '
NR==1 {
    for(i=1;i<=NF;i++) {
        if($i=="project.project_id") project=i
        if($i=="samples.sample_type") type=i
    }
    next
}
$project != "" && $type != "" {
    count[$project,$type]++
}
END {
    print "Project\tSample_Type\tCount"
    for(x in count) {
        split(x,a,SUBSEP)
        print a[1] "\t" a[2] "\t" count[x]
    }
}' "$SAMPLE" | sort

# ------------------------------------------------------------
# Tumor + normal sample QC
# ------------------------------------------------------------

echo
echo "Checking cases with both tumor and normal samples..."

awk -F '\t' '
NR==1 {
    for(i=1;i<=NF;i++) {
        if($i=="cases.submitter_id") caseid=i
        if($i=="samples.sample_type") sampletype=i
    }
    next
}
{
    if($sampletype=="Primary Tumor")
        tumor[$caseid]=1

    if($sampletype=="Blood Derived Normal" ||
       $sampletype=="Solid Tissue Normal")
        normal[$caseid]=1
}
END {
    both=0

    for(x in tumor)
        if(x in normal)
            both++

    print "Cases with Primary Tumor + Normal sample:", both
}' "$SAMPLE"

# ------------------------------------------------------------
# Matched normal barcode check
# ------------------------------------------------------------

echo
echo "Checking matched-normal barcodes in MAF..."

zcat "$MAF" | awk -F '\t' '
NR==1 {
    for(i=1;i<=NF;i++)
        if($i=="Matched_Norm_Sample_Barcode") normal=i
    next
}
{
    if(normal && $normal != "")
        filled++
    else
        missing++
}
END {
    print "Records with matched normal:", filled
    print "Records without matched normal:", missing
}'

echo
echo "============================================================"
echo "01 — INPUT QC COMPLETE"
echo "============================================================"
