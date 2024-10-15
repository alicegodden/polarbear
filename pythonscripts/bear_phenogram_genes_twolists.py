# Title: Bear phenogram for plotting genes
# Author: Dr. Alice M. Godden

# Import libraries
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Set the style to flare
sns.set(style='whitegrid')

# Read chromosome end data
chrom_end_df = pd.read_csv("chrom_end_bear.txt", header=None, names=["Chromosome", "Length"], delim_whitespace=True)

# Read primary gene data
genes_df = pd.read_csv("genes_gender_group_ASM", header=None, names=["Chrom", "Start", "End", "Name"],
                       delim_whitespace=True)

# Read additional gene data
additional_genes_df = pd.read_csv("genes_gender_group_temp_ASM", header=None, names=["Chrom", "Start", "End", "Name"],
                                  delim_whitespace=True)

# Convert Start and End to integers
for df in [genes_df, additional_genes_df]:
    df["Start"] = pd.to_numeric(df["Start"], errors='coerce')
    df["End"] = pd.to_numeric(df["End"], errors='coerce')

chrom_end_df["Length"] = pd.to_numeric(chrom_end_df["Length"], errors='coerce')

# Combine both gene dataframes
combined_genes_df = pd.concat([genes_df, additional_genes_df], ignore_index=True)

# Get unique chromosomes with genes
unique_chromosomes = combined_genes_df["Chrom"].unique()

# Filter chromosome end data for unique chromosomes with genes and sort by chromosome number
chrom_end_df_filtered = chrom_end_df[chrom_end_df["Chromosome"].isin(unique_chromosomes)]
chrom_end_df_filtered = chrom_end_df_filtered.sort_values(by="Chromosome")

# Print filtered chromosome end data for debugging
print("Filtered Chromosome End Data:")
print(chrom_end_df_filtered)

# Print gene data for debugging
print("Combined Gene Data:")
print(combined_genes_df)

# Plot chromosomes with corresponding gene points
plt.figure(figsize=(12, 8))  # Increased size for better clarity
unplotted_genes_count = 0  # Counter for genes that don't get plotted
plotted_genes_count = 0  # Counter for genes that do get plotted

# Define colors for different gene sets
flare_palette = sns.color_palette("flare", n_colors=2)
primary_color = flare_palette[0]  # First color from flare palette
additional_color = flare_palette[1]  # Second color from flare palette

# Determine y positions for chromosomes, increasing spacing for better readability
y_positions = np.arange(len(chrom_end_df_filtered)) * 3  # Increased spacing between chromosomes

# Iterate through filtered chromosome end data
for index, (y_pos, row) in enumerate(zip(y_positions, chrom_end_df_filtered.iterrows())):
    chrom = row[1]["Chromosome"]
    length = row[1]["Length"]

    # Draw chromosome line
    plt.plot([0, length], [y_pos, y_pos], color="black")  # Use index for y position

    # Filter genes for the current chromosome
    genes_on_chrom = combined_genes_df[combined_genes_df["Chrom"] == chrom]

    # Create a list to keep track of y offsets for labels
    used_offsets = []  # This will keep track of label positions

    # Plot genes
    for idx, (_, gene_row) in enumerate(genes_on_chrom.iterrows()):
        gene_start = gene_row["Start"]
        gene_name = gene_row["Name"]

        if gene_start <= length:
            # Determine color based on whether it's a primary or additional gene
            color = primary_color if gene_row["Chrom"] in genes_df["Chrom"].values else additional_color

            # Plot a point for each gene
            plt.scatter(gene_start, y_pos, color=color, s=30, zorder=10, alpha=0.8)

            # Alternate the offset based on the index and manage overlaps
            base_offset = 1.5  # Base offset
            label_offset = base_offset + (len(used_offsets) % 2) * 0.75  # Alternate above/below

            # Check for overlaps
            while label_offset in used_offsets:
                label_offset += 1.25  # Increase the offset if it's already used

            # Add the offset to the list of used offsets
            used_offsets.append(label_offset)

            # Annotate the gene using an arrow
            plt.annotate(gene_name, xy=(gene_start, y_pos),
                         xytext=(gene_start, y_pos + label_offset),  # Position with calculated offset
                         #arrowprops=dict(arrowstyle='-', color='black'),
                         fontsize=6, color='black', ha='center')
            plotted_genes_count += 1  # Increment plotted gene count
        else:
            print(
                f"Skipping gene {gene_name} on {chrom} due to start position {gene_start} exceeding chromosome length.")
            unplotted_genes_count += 1  # Increment unplotted gene count

# Plot customization
plt.xlabel("Genomic Position", fontweight='bold')
plt.ylabel("Chromosome/Scaffold", fontweight='bold')
plt.title("Genomic Positions of signigicantly differentially expressed genes from ASM1731132v1", fontweight='bold')

# Customize y-ticks and apply bold font
plt.yticks(y_positions, chrom_end_df_filtered["Chromosome"], fontweight='bold')
plt.xticks(fontweight='bold')
plt.tick_params(axis='y', which='both', length=5)  # Add tick marks on the y-axis
plt.gca().invert_yaxis()  # Invert y-axis to plot chromosomes from top to bottom
plt.grid(False)  # Remove gridlines

# Create a legend
plt.scatter([], [], color=primary_color, label='Sex + Group', alpha=0.8)
plt.scatter([], [], color=additional_color, label='Sex + Group + Temperature', alpha=0.8)
plt.legend(title='Gene Types')

plt.savefig('polarbear_newgen_Adult_telescope_sig_deGENES.png', dpi=600, bbox_inches='tight')

# Print out the number of plotted and unplotted genes
print(f"Number of genes that do get plotted: {plotted_genes_count}")
print(f"Number of genes that do not get plotted: {unplotted_genes_count}")

plt.show()  # Optional: Show the plot
