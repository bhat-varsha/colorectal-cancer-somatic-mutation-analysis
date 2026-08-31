import pandas as pd
from scipy.stats import mannwhitneyu

BURDEN_FILE = "results/mutation_burden_per_case.tsv"
PROJECT_FILE = "results/case_project_lookup.tsv"

burden = pd.read_csv(
    BURDEN_FILE,
    sep="\t",
    header=None,
    names=["case_id", "mutation_burden"]
)

project = pd.read_csv(
    PROJECT_FILE,
    sep="\t",
    header=None,
    names=["case_id", "project"]
)

df = burden.merge(
    project,
    on="case_id",
    how="inner"
)

coad = df.loc[
    df["project"] == "TCGA-COAD",
    "mutation_burden"
]

read = df.loc[
    df["project"] == "TCGA-READ",
    "mutation_burden"
]

u_stat, p_value = mannwhitneyu(
    coad,
    read,
    alternative="two-sided"
)

print("Mutation burden comparison: COAD vs READ")
print("------------------------------------------")
print(f"COAD cases: {len(coad)}")
print(f"READ cases: {len(read)}")
print(f"COAD median: {coad.median()}")
print(f"READ median: {read.median()}")
print(f"Mann-Whitney U statistic: {u_stat}")
print(f"P-value: {p_value}")

if p_value < 0.05:
    print("Result: Statistically significant difference (p < 0.05)")
else:
    print("Result: No statistically significant difference (p >= 0.05)")
