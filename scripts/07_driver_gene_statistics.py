import pandas as pd
from scipy.stats import fisher_exact

INPUT = "results/driver_gene_prevalence_coad_read.tsv"
OUTPUT = "results/driver_gene_coad_read_statistics.tsv"

COAD_TOTAL = 428
READ_TOTAL = 155


# ------------------------------------------------------------
# Benjamini-Hochberg FDR correction
# ------------------------------------------------------------

def benjamini_hochberg(p_values):

    n = len(p_values)

    # Sort p-values while remembering their original positions
    sorted_indices = sorted(
        range(n),
        key=lambda i: p_values[i]
    )

    fdr = [0.0] * n

    # Calculate adjusted p-values
    for rank, index in enumerate(sorted_indices, start=1):
        fdr[index] = p_values[index] * n / rank

    # Ensure monotonicity from largest rank to smallest
    for i in range(n - 2, -1, -1):

        current_index = sorted_indices[i]
        next_index = sorted_indices[i + 1]

        fdr[current_index] = min(
            fdr[current_index],
            fdr[next_index]
        )

    # FDR cannot exceed 1
    fdr = [min(value, 1.0) for value in fdr]

    return fdr


# ------------------------------------------------------------
# Read input
# ------------------------------------------------------------

df = pd.read_csv(INPUT, sep="\t")

results = []


# ------------------------------------------------------------
# Fisher's exact test for each driver gene
# ------------------------------------------------------------

for _, row in df.iterrows():

    gene = row["Gene"]

    coad_mut = int(row["COAD_cases"])
    read_mut = int(row["READ_cases"])

    coad_not_mut = COAD_TOTAL - coad_mut
    read_not_mut = READ_TOTAL - read_mut

    table = [
        [coad_mut, coad_not_mut],
        [read_mut, read_not_mut]
    ]

    odds_ratio, p_value = fisher_exact(table)

    results.append([
        gene,
        coad_mut,
        row["COAD_percent"],
        read_mut,
        row["READ_percent"],
        odds_ratio,
        p_value
    ])


# ------------------------------------------------------------
# Create results table
# ------------------------------------------------------------

result = pd.DataFrame(
    results,
    columns=[
        "Gene",
        "COAD_cases",
        "COAD_percent",
        "READ_cases",
        "READ_percent",
        "Odds_Ratio",
        "P_value"
    ]
)


# ------------------------------------------------------------
# Benjamini-Hochberg FDR correction
# ------------------------------------------------------------

result["FDR"] = benjamini_hochberg(
    result["P_value"].tolist()
)

result["Significant_FDR_0.05"] = result["FDR"] < 0.05


# ------------------------------------------------------------
# Sort by p-value
# ------------------------------------------------------------

result = result.sort_values("P_value")


# ------------------------------------------------------------
# Save results
# ------------------------------------------------------------

result.to_csv(
    OUTPUT,
    sep="\t",
    index=False
)


# ------------------------------------------------------------
# Print summary
# ------------------------------------------------------------

print("Driver gene statistical comparison: COAD vs READ")
print("------------------------------------------------")
print(result.to_string(index=False))
print()
print(f"Created: {OUTPUT}")
