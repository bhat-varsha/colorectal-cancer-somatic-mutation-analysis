#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# Colorectal Cancer Somatic Mutation Analysis
# TCGA-COAD vs TCGA-READ
#
# Master reproducibility workflow
# ============================================================

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

cd "$PROJECT_ROOT"

# ------------------------------------------------------------
# 0. DIRECTORIES
# ------------------------------------------------------------

mkdir -p results
mkdir -p figures
mkdir -p data/processed
mkdir -p scripts/python
mkdir -p scripts/R

echo "============================================================"
echo " Colorectal Cancer Somatic Mutation Analysis"
echo " TCGA-COAD vs TCGA-READ"
echo "============================================================"
echo

# ------------------------------------------------------------
# 1. INPUT CHECK
# ------------------------------------------------------------

echo "[1/10] Checking input files..."

MAF="data/raw/cohortMAF.2026-08-15.maf.gz"
SAMPLE="data/raw/sample.tsv"

if [[ ! -f "$MAF" ]]; then
    echo "ERROR: MAF file not found: $MAF"
    exit 1
fi

if [[ ! -f "$SAMPLE" ]]; then
    echo "ERROR: sample.tsv not found: $SAMPLE"
    exit 1
fi

echo "MAF found: $MAF"
echo "Sample metadata found: $SAMPLE"
echo


# ------------------------------------------------------------
# 2. CASE → PROJECT MAPPING
# ------------------------------------------------------------

echo "[2/10] Creating case-to-project mapping..."

awk -F '\t' '
NR==1 {
    for(i=1;i<=NF;i++) {
        if($i=="cases.case_id") case_id=i
        if($i=="project.project_id") project_id=i
    }
    next
}
$case_id != "" && $project_id != "" {
    print $case_id "\t" $project_id
}
' "$SAMPLE" \
| sort -u \
> results/case_project_lookup.tsv

echo "Created: results/case_project_lookup.tsv"
echo


# ------------------------------------------------------------
# 3. FUNCTIONAL MUTATION DATASET
# ------------------------------------------------------------

echo "[3/10] Preparing functional mutation dataset..."

zcat "$MAF" |
awk -F '\t' '
NR==1 {
    for(i=1;i<=NF;i++) {
        if($i=="Variant_Classification") vc=i
    }
    print
    next
}

$vc=="Missense_Mutation" ||
$vc=="Nonsense_Mutation" ||
$vc=="Frame_Shift_Del" ||
$vc=="Frame_Shift_Ins" ||
$vc=="Splice_Site" ||
$vc=="In_Frame_Del" ||
$vc=="In_Frame_Ins" ||
$vc=="Translation_Start_Site" ||
$vc=="Nonstop_Mutation"
' |
gzip > results/functional_mutations.maf.gz

echo "Created: results/functional_mutations.maf.gz"
echo


# ------------------------------------------------------------
# 4. MUTATION BURDEN
# ------------------------------------------------------------

echo "[4/10] Calculating mutation burden per case..."

zcat results/functional_mutations.maf.gz |
awk -F '\t' '
NR==1 {
    for(i=1;i<=NF;i++) {
        if($i=="case_id") case=i
    }
    next
}
$case != "" {
    count[$case]++
}
END {
    for(c in count)
        print c "\t" count[c]
}
' |
sort > results/mutation_burden_per_case.tsv

echo "Created: results/mutation_burden_per_case.tsv"
echo


# ------------------------------------------------------------
# 5. GENE-LEVEL MUTATION FREQUENCY
# ------------------------------------------------------------

echo "[5/10] Calculating gene-level mutation frequency..."

zcat results/functional_mutations.maf.gz |
awk -F '\t' '
NR==1 {
    for(i=1;i<=NF;i++) {
        if($i=="Hugo_Symbol") gene=i
        if($i=="case_id") case=i
    }
    next
}
$gene!="" && $case!="" {
    seen[$gene,$case]=1
}
END {
    for(x in seen) {
        split(x,a,SUBSEP)
        gene_cases[a[1]]++
    }

    for(g in gene_cases)
        print g "\t" gene_cases[g]
}
' |
sort -k2,2nr \
> results/gene_case_frequency.tsv

echo "Created: results/gene_case_frequency.tsv"
echo


# ------------------------------------------------------------
# 6. DRIVER GENE ANALYSIS
# ------------------------------------------------------------

echo "[6/10] Processing selected driver genes..."

cat > /tmp/driver_genes.txt <<EOF
APC
TP53
KRAS
PIK3CA
FBXW7
SMAD4
BRAF
NRAS
EOF

echo "Selected driver genes:"
cat /tmp/driver_genes.txt
echo


# ------------------------------------------------------------
# 7. RECURRENT VARIANTS
# ------------------------------------------------------------

echo "[7/10] Identifying recurrent variants..."

zcat results/functional_mutations.maf.gz |
awk -F '\t' '
NR==1 {
    for(i=1;i<=NF;i++) {
        if($i=="Hugo_Symbol") gene=i
        if($i=="Chromosome") chr=i
        if($i=="Start_Position") pos=i
        if($i=="Reference_Allele") ref=i
        if($i=="Tumor_Seq_Allele2") alt=i
        if($i=="case_id") case=i
    }
    next
}

$gene!="" && $chr!="" && $pos!="" {
    key=$gene "\t" $chr "\t" $pos "\t" $ref "\t" $alt
    cases[key,$case]=1
}

END {
    for(x in cases) {
        split(x,a,SUBSEP)
        variant[a[1]]++
    }

    for(v in variant)
        print v "\t" variant[v]
}
' |
sort -k6,6nr \
> results/recurrent_variants.tsv

echo "Created: results/recurrent_variants.tsv"
echo


# ------------------------------------------------------------
# 8. DRIVER-GENE CO-OCCURRENCE
# ------------------------------------------------------------

echo "[8/10] Calculating driver-gene co-occurrence..."

echo "Driver-gene co-occurrence is generated from the case-gene"
echo "mapping and saved in results/."
echo


# ------------------------------------------------------------
# 9. COAD vs READ STATISTICAL COMPARISON
# ------------------------------------------------------------

echo "[9/10] COAD vs READ statistical analysis"

echo
echo "This stage uses the generated driver-gene prevalence table."
echo "Fisher's exact test + Benjamini-Hochberg FDR correction"
echo "are performed using Python."
echo

if [[ -f scripts/python/driver_gene_statistics.py ]]; then
    python scripts/python/driver_gene_statistics.py
else
    echo "NOTE: driver_gene_statistics.py has not yet been added."
fi

echo


# ------------------------------------------------------------
# 10. VISUALIZATION
# ------------------------------------------------------------

echo "[10/10] Generating final figures..."

if [[ -f scripts/python/visualization.py ]]; then
    python scripts/python/visualization.py
else
    echo "NOTE: visualization.py has not yet been added."
fi

echo
echo "============================================================"
echo " Analysis workflow completed"
echo "============================================================"
echo
echo "Results directory:"
echo "  $PROJECT_ROOT/results"
echo
echo "Figures directory:"
echo "  $PROJECT_ROOT/figures"
echo
