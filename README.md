# Colorectal Cancer Somatic Mutation Analysis

## TCGA-COAD vs TCGA-READ

This project performs cohort-level somatic mutation analysis of colorectal cancer using TCGA/GDC-derived mutation data.

The analysis focuses on:

- Input mutation-data quality control
- Functional mutation filtering
- Mutation burden per case
- Gene-level mutation frequency and prevalence
- Selected colorectal cancer driver genes
- Recurrent variants
- Driver-gene co-occurrence
- COAD vs READ comparison
- Fisher's exact test
- Benjamini-Hochberg FDR correction
- SNV mutation spectrum
- Transition/transversion analysis
- Visualization of major findings

---

## Project Objective

The primary objective is to characterize somatic mutation patterns in colorectal cancer and compare:

- TCGA-COAD — Colon Adenocarcinoma
- TCGA-READ — Rectum Adenocarcinoma

The analysis identifies frequently mutated genes, recurrent variants, driver-gene patterns, co-occurring driver genes, and differences in driver-gene prevalence between COAD and READ.

---

## Data

Primary mutation dataset:

data/raw/cohortMAF.2026-08-15.maf.gz

Sample/project metadata:

data/raw/sample.tsv

The analysis uses case-level identifiers to associate mutations with TCGA projects.

---

## Analysis Workflow

### 1. Input QC

Script:

scripts/01_input_qc.sh

Checks the input MAF and sample metadata and summarizes important mutation-data fields.

### 2. Functional Mutation Filtering

Script:

scripts/02_functional_filter.sh

Extracts functional mutation classes for downstream analysis.

Output:

results/functional_mutations.maf.gz

### 3. Mutation Burden

Script:

scripts/03_mutation_burden.sh

Calculates the number of functional mutations observed per case.

Output:

results/mutation_burden_per_case.tsv

Additional statistics:

scripts/python/mutation_burden_statistics.py

### 4. Gene-Level Mutation Frequency

Script:

scripts/05_gene_level_frequency.sh

Calculates the number of unique cases containing mutations in each gene.

Outputs:

results/gene_case_frequency.tsv
results/gene_case_frequency_percent.tsv
results/gene_frequency_coad_read.tsv

### 5. Gene Prevalence

Script:

scripts/05_gene_prevalence.sh

Calculates mutation prevalence separately in COAD and READ.

Output:

results/gene_prevalence_coad_read.tsv

### 6. Driver Gene Analysis

Script:

scripts/06_driver_gene_analysis.sh

Selected colorectal cancer driver genes:

APC
TP53
KRAS
PIK3CA
FBXW7
SMAD4
BRAF
NRAS

Functional mutations in these genes are extracted.

Output:

results/driver_gene_functional_mutations.tsv

Driver-gene prevalence:

results/driver_gene_prevalence_coad_read.tsv

### 7. Driver Gene Statistical Comparison

Script:

scripts/07_driver_gene_statistics.py

COAD and READ driver-gene prevalence are compared using:

- Fisher's exact test
- Benjamini-Hochberg FDR correction

Output:

results/driver_gene_coad_read_statistics.tsv

Genes significant after FDR correction at FDR < 0.05:

- BRAF
- PIK3CA
- TP53

### 8. Driver Gene Co-occurrence

Driver-gene co-occurrence is calculated at the case level.

Outputs:

results/driver_gene_cooccurrence.tsv
results/driver_gene_cooccurrence_matrix.tsv
results/driver_gene_cooccurrence_heatmap.png
results/driver_gene_cooccurrence_heatmap.pdf

The strongest observed co-occurrences include:

- APC + TP53
- APC + KRAS
- APC + PIK3CA
- TP53 + KRAS
- KRAS + PIK3CA

### 9. Recurrent Variant Analysis

Recurrent variants are identified using:

- Gene
- Chromosome
- Position
- Reference allele
- Alternate allele

Outputs:

results/recurrent_variants.tsv
results/driver_recurrent_variants.tsv
results/top20_driver_recurrent_variants.tsv
results/top_recurrent_driver_variants.png

The most recurrent driver variants include variants in:

- KRAS
- BRAF
- TP53
- APC
- PIK3CA

### 10. Mutation Spectrum

Script:

scripts/10_mutation_spectrum.sh

Calculates SNV substitution frequencies and the transition/transversion ratio.

Output:

results/mutation_spectrum.tsv

Mutation-spectrum visualization:

scripts/11_mutation_spectrum_plot.py

Current mutation-spectrum result:

Total SNVs: 205841
Transitions: 131609
Transversions: 74232
Ti/Tv ratio: 1.77294

---

## Selected Driver Gene Prevalence

| Gene | COAD % | READ % |
|------|-------:|-------:|
| APC | 75.99 | 90.58 |
| TP53 | 57.92 | 80.43 |
| KRAS | 45.30 | 45.65 |
| PIK3CA | 31.44 | 15.94 |
| FBXW7 | 17.57 | 18.84 |
| SMAD4 | 11.63 | 16.67 |
| BRAF | 15.35 | 3.62 |
| NRAS | 5.20 | 10.87 |

---

## COAD vs READ Statistical Results

| Gene | Odds Ratio | P-value | FDR | Significant |
|------|-----------:|--------:|----:|-------------|
| BRAF | 5.08 | 0.000054 | 0.000413 | Yes |
| PIK3CA | 2.55 | 0.000103 | 0.000413 | Yes |
| TP53 | 0.48 | 0.000274 | 0.000732 | Yes |
| APC | 0.61 | 0.032355 | 0.064709 | No |
| NRAS | 0.48 | 0.049453 | 0.079125 | No |
| SMAD4 | 0.71 | 0.247983 | 0.330644 | No |
| KRAS | 1.09 | 0.704380 | 0.805006 | No |
| FBXW7 | 0.99 | 1.000000 | 1.000000 | No |

---

## Mutation Spectrum

Total SNVs:

205841

Transitions:

131609

Transversions:

74232

Transition/Transversion ratio:

1.77294

The most frequent substitution classes are C>T and G>A.

---

## Project Structure

colorectal-cancer-somatic-mutation-analysis/
│
├── data/
│   ├── raw/
│   │   ├── cohortMAF.2026-08-15.maf.gz
│   │   └── sample.tsv
│   └── processed/
│
├── results/
│
├── figures/
│
├── scripts/
│   ├── 01_input_qc.sh
│   ├── 02_functional_filter.sh
│   ├── 03_mutation_burden.sh
│   ├── 05_gene_level_frequency.sh
│   ├── 05_gene_prevalence.sh
│   ├── 06_driver_gene_analysis.sh
│   ├── 07_driver_gene_statistics.py
│   ├── 09_driver_gene_pair_statistics.py
│   ├── 10_mutation_spectrum.sh
│   ├── 11_mutation_spectrum_plot.py
│   ├── python/
│   └── R/
│
├── README.md
├── all_analysis_commands.txt
└── full_analysis_history.txt

---

## Reproducibility

The analysis was performed using standard command-line tools and Python/R.

Main tools include:

- Bash
- awk
- grep
- sort
- gzip/zcat
- Python
- pandas
- SciPy
- Matplotlib
- R

No additional statsmodels dependency is required.

Fisher's exact test is performed using SciPy.

Benjamini-Hochberg FDR correction is implemented directly in the statistical analysis script.

---

## Reproducing the Analysis

From the project root:

bash scripts/01_input_qc.sh
bash scripts/02_functional_filter.sh
bash scripts/03_mutation_burden.sh
bash scripts/05_gene_level_frequency.sh
bash scripts/05_gene_prevalence.sh
bash scripts/06_driver_gene_analysis.sh
python scripts/07_driver_gene_statistics.py
python scripts/09_driver_gene_pair_statistics.py
bash scripts/10_mutation_spectrum.sh
python scripts/11_mutation_spectrum_plot.py

---

## Final Validation

The project contains:

- Raw input data
- Functional mutation dataset
- Mutation burden results
- Gene-level mutation frequency results
- Gene prevalence results
- Driver-gene analysis
- Driver-gene statistical comparison
- Driver-gene co-occurrence analysis
- Recurrent variant analysis
- Mutation spectrum analysis
- Visualization files
- Analysis scripts
- Reproducibility documentation
- Analysis command/history records

---

## Project Status

ANALYSIS COMPLETED

The project provides a reproducible cohort-level somatic mutation analysis of colorectal cancer, including comparison of TCGA-COAD and TCGA-READ and focused analysis of colorectal cancer driver genes.
