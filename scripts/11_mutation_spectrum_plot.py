import pandas as pd
import matplotlib.pyplot as plt

INPUT = "results/mutation_spectrum.tsv"
OUTPUT_PNG = "results/mutation_spectrum.png"
OUTPUT_PDF = "results/mutation_spectrum.pdf"

df = pd.read_csv(INPUT, sep="\t")

# Remove TOTAL_SNV row
df = df[df["Mutation_Type"] != "TOTAL_SNV"]

# Plot
plt.figure(figsize=(10, 6))

plt.bar(
    df["Mutation_Type"],
    df["Count"]
)

plt.xlabel("Mutation Type")
plt.ylabel("Number of SNVs")
plt.title("Somatic Mutation Spectrum")

plt.xticks(rotation=45)
plt.tight_layout()

plt.savefig(OUTPUT_PNG, dpi=300)
plt.savefig(OUTPUT_PDF)

plt.close()

print("Mutation spectrum visualization complete.")
print(f"Created: {OUTPUT_PNG}")
print(f"Created: {OUTPUT_PDF}")
