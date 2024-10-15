# Title: Matching gene and TE loci
# Author: Dr. Alice M. Godden

import pandas as pd

# Read the TE data from the file
te_data = pd.read_csv("matched_TEs_output_TEMP_ASM.txt", sep="\t", header=None, names=["Chrom", "Start", "End", "Name"])

# Read the gene data from the file
gene_data = pd.read_csv("genes_gender_group_temp_ASM", sep="\t", header=0)

# Initialize a list to store overlapping genes
overlapping_genes = []

# Check for overlaps
for _, te_row in te_data.iterrows():
    te_chrom = te_row["Chrom"]
    te_start = te_row["Start"]
    te_end = te_row["End"]

    # Filter genes for the current chromosome
    overlapping_genes_on_chrom = gene_data[(gene_data["Chrom"] == te_chrom) &
                                           (gene_data["Start"] < te_end) &
                                           (gene_data["End"] > te_start)]

    # Add overlapping genes to the list
    for _, gene_row in overlapping_genes_on_chrom.iterrows():
        overlapping_genes.append({
            "TE_Chrom": te_chrom,
            "TE_Start": te_start,
            "TE_End": te_end,
            "TE_Name": te_row["Name"],
            "Gene_Chrom": gene_row["Chrom"],
            "Gene_Start": gene_row["Start"],
            "Gene_End": gene_row["End"],
            "Gene_Name": gene_row["Name"]
        })

# Convert the list of overlapping genes to a DataFrame
overlapping_df = pd.DataFrame(overlapping_genes)

# Output the results
if not overlapping_df.empty:
    print("Overlapping Gene Positions Found:")
    print(overlapping_df)
else:
    print("No overlapping gene positions found.")
