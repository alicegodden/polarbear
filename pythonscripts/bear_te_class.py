# Title: TE class plots
# Author: Dr. Alice M. Godden
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import hypergeom
import pandas as pd
import numpy as np

# Data
te_classes = ['LINE', 'LTR', 'DNA', 'SINE']
sig_de_counts = np.array([152, 4, 13, 10])  # temp, Significantly differentially expressed counts # countsASM temp 1249, 125, 94, 66
non_sig_de_counts = np.array([92636, 3213, 5320, 3528])  # Non-significant DE counts gender + group
# non_sig_de_counts = np.array([207395, 13228, 17427, 14428])  for temp+gender+group

# Calculate total counts
total_sig_de = sum(sig_de_counts)  # Total sig DE across all TEs
total_non_sig_de = sum(non_sig_de_counts)  # Total non-sig DE across all TEs
total_tes = total_sig_de + total_non_sig_de  # Total TEs

# Expected counts assuming equal probability
total_te_counts = sig_de_counts + non_sig_de_counts  # Total TE counts per class
expected_sig_counts = total_te_counts * (total_sig_de / total_tes)  # Expected sig DEs per class

# Perform enrichment test (hypergeometric)
p_values = []
for i, te_class in enumerate(te_classes):
    observed = sig_de_counts[i]
    expected = expected_sig_counts[i]
    total_in_class = total_te_counts[i]

    # Hypergeometric test: is sig DE count higher than expected?
    p_val = hypergeom.sf(observed - 1, total_tes, total_sig_de, total_in_class)
    p_values.append(p_val)

# Print p-values to console
print("\n### P-values for Each TE Class ###")
for te_class, p_val in zip(te_classes, p_values):
    print(f"{te_class}: p = {p_val:.4e}")

# Save p-values to CSV
p_value_df = pd.DataFrame({'TE Class': te_classes, 'Observed Sig DE': sig_de_counts,
                           'Expected Sig DE': expected_sig_counts, 'p-value': p_values})
p_value_df.to_csv("teclass_individual_hypergeom_p_values.csv", index=False)

# Define color palette
flare_palette = sns.color_palette("rocket")
element_colors = {"DNA": flare_palette[0], "LINE": flare_palette[2], "SINE": flare_palette[4], "LTR": flare_palette[5]}

# Create plot
fig, ax = plt.subplots(figsize=(10, 6))

x = np.arange(len(te_classes))
bar_width = 0.6

# Plot observed sig DE counts
bars = ax.bar(x, sig_de_counts, bar_width, label='Observed Sig DE', color=[element_colors[te] for te in te_classes])

# Plot expected sig DE counts as horizontal lines
for i, (expected, obs) in enumerate(zip(expected_sig_counts, sig_de_counts)):
    ax.hlines(expected, x[i] - bar_width / 2, x[i] + bar_width / 2, color='black', linewidth=2, linestyle='dashed', label="Expected" if i == 0 else "")

# Set labels and title
ax.set_xticks(x)
ax.set_xticklabels(te_classes, rotation=45, ha='right', fontweight='bold')
ax.set_xlabel("TE Class", fontweight='bold')
ax.set_ylabel("Number of Significant DE TEs", fontweight='bold')
ax.set_title("Observed vs. Expected Significant DE TEs per Class", fontweight='bold')

# Set the ticks fontweight
plt.xticks(fontweight='bold')
plt.yticks(fontweight='bold')

# Define significance thresholds and symbols
significance_levels = [(0.001, "***"), (0.01, "**"), (0.05, "*")]

# Add significance stars above bars
y_offset = max(sig_de_counts) * 0.05  # Space above bars for stars
for i, p_val in enumerate(p_values):
    star_label = next((s for p, s in significance_levels if p_val < p), "")
    if star_label:
        ax.text(x[i], sig_de_counts[i] + y_offset, star_label, ha='center', va='bottom', fontsize=14, fontweight='bold', color='black')

# Add legend
ax.legend()

# Adjust layout
plt.tight_layout()

# Save and show plot
plt.savefig('teclass_individual_hypergeom_enrichment_gender+group.png', dpi=600)
plt.show()
