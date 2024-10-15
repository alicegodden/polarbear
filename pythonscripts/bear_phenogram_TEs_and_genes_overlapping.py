# Title: Chromosomal plot of TEs and genes
# Author: Dr. Alice M. Godden

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Read chromosome end data
chrom_end_df = pd.read_csv("chrom_end_bear.txt", header=None, names=["Chromosome", "Length"], delim_whitespace=True)

# Read points data (TE data)
points_df = pd.read_csv("matched_TEs_output_TEMP_ASM.txt", names=["Chromosome", "Start", "End", "Name"],
                        delim_whitespace=True)

# Read gene data
gene_data = pd.read_csv("genes_gender_group_temp_ASM", sep="\t", header=0)

# Get unique chromosomes with hits
chromosomes_with_hits = points_df["Chromosome"].unique()

# Filter chromosome end data for chromosomes with hits and sort by chromosome number
chrom_end_df_filtered = chrom_end_df[chrom_end_df["Chromosome"].isin(chromosomes_with_hits)]
chrom_end_df_filtered = chrom_end_df_filtered.sort_values(by="Chromosome")

# Define colors for different element types
try:
    flare_palette = sns.color_palette("rocket")
    element_colors = {
        "DNA": flare_palette[0],
        "LINE": flare_palette[2],
        "SINE": flare_palette[4],
        "LTR": flare_palette[5]
    }
except IndexError:
    print("Error: Flare palette does not contain enough colors.")
    element_colors = {
        "DNA": "blue",
        "LINE": "green",
        "SINE": "red",
        "LTR": "purple"
    }

# Initialize a list to store overlapping genes
overlapping_genes = []

# Check for overlaps
for _, te_row in points_df.iterrows():
    te_chrom = te_row["Chromosome"]
    te_start = te_row["Start"]
    te_end = te_row["End"]

    # Filter genes for the current chromosome and check for overlaps
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

# Plot chromosomes with hits and corresponding points
plt.figure(figsize=(10, 6))
legend_handles = []
legend_labels = set()  # Set to store legend labels

for index, row in chrom_end_df_filtered.iterrows():
    chrom = row["Chromosome"]
    length = row["Length"]
    plt.plot([0, length], [chrom, chrom], color="black")  # Draw chromosome line
    points_on_chrom = points_df[points_df["Chromosome"] == chrom]

    for _, point_row in points_on_chrom.iterrows():
        start = point_row["Start"]
        end = point_row["End"]
        element_type = point_row["Name"].split("#")[1].split("/")[0]  # Extract element type from the Name column

        # Check if element_type is in element_colors
        if element_type in element_colors:
            color = element_colors[element_type]  # Get color for the element type
            plt.scatter([start, end], [chrom, chrom], color=color, alpha=0.6, zorder=10, s=20)  # Smaller points

            # Add to legend only if it hasn't been added before
            if element_type not in legend_labels:
                legend_handles.append(plt.scatter([], [], color=color, label=element_type, alpha=0.6))
                legend_labels.add(element_type)

# Plot overlapping genes
for gene in overlapping_genes:
    plt.scatter(gene["TE_Start"], gene["TE_Chrom"], color='navy', marker='D', s=50, label='Overlapping Genes')

# Add overlapping genes to the legend
legend_handles.append(plt.scatter([], [], color='navy', marker='D', label='Overlapping Genes', s=50))

plt.xlabel("Genomic Position", fontweight='bold')
plt.ylabel("Chromosome/Scaffold", fontweight='bold')
plt.title("Genomic Positions of Significantly Differentially Expressed TEs (Temperature)", fontweight='bold')
plt.yticks(chrom_end_df_filtered["Chromosome"], fontweight='bold', fontsize=8)  # Smaller font size for y-ticks
plt.xticks(fontweight='bold')
plt.legend(handles=legend_handles, title="TE Class")
plt.gca().invert_yaxis()  # Invert y-axis to plot chromosomes from top to bottom
plt.grid(False)  # Remove gridlines

# Save and show the plot
plt.savefig('polarbear_newgen_Adult_telescope_sig_deTEs_TEMPERATUREgenes.png', dpi=600, bbox_inches='tight')
plt.show()  # Optional: Show the plot
