import pandas as pd
from scipy.stats import fisher_exact


INPUT = "results/case_gene_pairs.tsv"
OUTPUT = "results/driver_gene_pair_statistics.tsv"

TOTAL_CASES = 583

DRIVER_GENES = [
    "APC",
    "TP53",
    "KRAS",
    "PIK3CA",
    "FBXW7",
    "SMAD4",
    "BRAF",
    "NRAS"
]


# ------------------------------------------------------------
# Benjamini-Hochberg FDR correction
# ------------------------------------------------------------

def benjamini_hochberg(p_values):

    n = len(p_values)

    sorted_indices = sorted(
        range(n),
        key=lambda i: p_values[i]
    )

    fdr = [0.0] * n

    for rank, index in enumerate(sorted_indices, start=1):
        fdr[index] = p_values[index] * n / rank

    for i in range(n - 2, -1, -1):

        current_index = sorted_indices[i]
        next_index = sorted_indices[i + 1]

        fdr[current_index] = min(
            fdr[current_index],
            fdr[next_index]
        )

    fdr = [min(value, 1.0) for value in fdr]

    return fdr


# ------------------------------------------------------------
# Read case-gene pairs
# ------------------------------------------------------------

df = pd.read_csv(
    INPUT,
    sep="\t",
    header=None,
    names=["case_id", "gene"]
)


# ------------------------------------------------------------
# Keep only selected driver genes
# ------------------------------------------------------------

df = df[df["gene"].isin(DRIVER_GENES)]


# Remove duplicate case-gene combinations
df = df.drop_duplicates(
    subset=["case_id", "gene"]
)


# ------------------------------------------------------------
# Create case × driver-gene mutation matrix
# ------------------------------------------------------------

matrix = pd.crosstab(
    df["case_id"],
    df["gene"]
)

matrix = matrix.reindex(
    columns=DRIVER_GENES,
    fill_value=0
)

matrix = (matrix > 0).astype(int)


# ------------------------------------------------------------
# Pairwise Fisher's exact tests
# ------------------------------------------------------------

results = []

for i in range(len(DRIVER_GENES)):

    gene1 = DRIVER_GENES[i]

    for j in range(i + 1, len(DRIVER_GENES)):

        gene2 = DRIVER_GENES[j]

        gene1_mut = matrix[gene1]
        gene2_mut = matrix[gene2]

        both = int(
            ((gene1_mut == 1) & (gene2_mut == 1)).sum()
        )

        gene1_only = int(
            ((gene1_mut == 1) & (gene2_mut == 0)).sum()
        )

        gene2_only = int(
            ((gene1_mut == 0) & (gene2_mut == 1)).sum()
        )

        neither = TOTAL_CASES - (
            both +
            gene1_only +
            gene2_only
        )

        table = [
            [both, gene1_only],
            [gene2_only, neither]
        ]

        odds_ratio, p_value = fisher_exact(table)

        results.append([
            gene1,
            gene2,
            both,
            gene1_only,
            gene2_only,
            neither,
            odds_ratio,
            p_value
        ])


# ------------------------------------------------------------
# Create results table
# ------------------------------------------------------------

result = pd.DataFrame(
    results,
    columns=[
        "Gene1",
        "Gene2",
        "Cooccurring_Cases",
        "Gene1_only",
        "Gene2_only",
        "Neither",
        "Odds_Ratio",
        "P_value"
    ]
)


# ------------------------------------------------------------
# FDR correction
# ------------------------------------------------------------

result["FDR"] = benjamini_hochberg(
    result["P_value"].tolist()
)

result["Significant_FDR_0.05"] = (
    result["FDR"] < 0.05
)


# ------------------------------------------------------------
# Sort by p-value
# ------------------------------------------------------------

result = result.sort_values(
    "P_value"
)


# ------------------------------------------------------------
# Save results
# ------------------------------------------------------------

result.to_csv(
    OUTPUT,
    sep="\t",
    index=False
)


# ------------------------------------------------------------
# Print results
# ------------------------------------------------------------

print("Driver gene pairwise association analysis")
print("------------------------------------------")
print(result.to_string(index=False))
print()
print(f"Created: {OUTPUT}")
