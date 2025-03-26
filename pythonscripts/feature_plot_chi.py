# Title: Features bar chart
# Subtitle: Chi-square analysis
# Author: Dr. Alice M. Godden

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from scipy.stats import chi2_contingency

# Set Rocket color palette from Seaborn
sns.set_palette("rocket")

# Data from the table
data = {
    "Category": ["CDS", "exon", "gene", "lnc_RNA", "mRNA", "pseudogene", "snRNA", "transcript"],
    "LTR": [17, 18, 54, 5, 110, 2, 0, 2],
    "LINE": [724, 855, 600, 93, 1062, 11, 4, 56],
    "DNA": [4, 9, 59, 12, 155, 0, 0, 1],
    "SINE": [1, 1, 46, 2, 111, 0, 0, 8]
}

# Convert to DataFrame
df = pd.DataFrame(data)

# Remove columns with all zeros
numerical_cols = df.columns[1:]
non_zero_cols = []
for col in numerical_cols:
    if (df[col] != 0).any():
        non_zero_cols.append(col)

df_filtered = df[["Category"] + non_zero_cols]

# Remove rows with all zeros (excluding 'Category' column)
df_filtered = df_filtered[df_filtered.iloc[:, 1:].sum(axis=1) > 0]

# Add a small pseudocount to prevent zero expected values
df_filtered.iloc[:, 1:] += 1e-6

# Calculate p-values for each feature-TE family combination
p_values = {}
te_types = df_filtered.columns[1:]
categories = df_filtered["Category"]

for te in te_types:
    p_values[te] = {}
    for category in categories:
        observed = df_filtered[df_filtered["Category"] == category][te].values.reshape(1, 1)
        expected_row_sum = df_filtered[te].sum()
        expected_col_sum = df_filtered[df_filtered["Category"] == category][te].sum()
        if expected_row_sum == 0 or expected_col_sum == 0:
            p_values[te][category] = np.nan
            continue;
        observed = np.array([[df_filtered.loc[df_filtered['Category'] == category, te].values[0], expected_row_sum - df_filtered.loc[df_filtered['Category'] == category, te].values[0]], [expected_col_sum - df_filtered.loc[df_filtered['Category'] == category, te].values[0], (df_filtered.iloc[:, 1:].sum().sum() - expected_row_sum - expected_col_sum + df_filtered.loc[df_filtered['Category'] == category, te].values[0]) ]])

        chi2_stat, p_value, dof, expected_values = chi2_contingency(observed)
        p_values[te][category] = p_value

# Print p-value table and save to CSV
p_value_df = pd.DataFrame(p_values)
print("P-value Table:")
print(p_value_df)
p_value_df.to_csv("chi_square_p_values_pbear25_individual.csv")

# Plotting
fig, ax = plt.subplots(figsize=(12, 8))
x = np.arange(len(categories))

bar_width = 0.2
offsets = np.linspace(-1.5, 1.5, df_filtered.shape[1] - 1) * bar_width

# Get the default colors from the Seaborn palette
colors = sns.color_palette("rocket")

for i, te in enumerate(te_types):
    observed_values = df_filtered[te].values
    expected = df_filtered.iloc[:, 1:].sum(axis = 0)[te]/df_filtered.iloc[:, 1:].sum().sum()*df_filtered.iloc[:, 1:].sum(axis = 1)
    for j, observed_val in enumerate(observed_values):
        p_val = p_values[te][categories[j]]
        star_label = ""
        if p_val < 0.001:
            star_label = "***"
        elif p_val < 0.01:
            star_label = "**"
        elif p_val < 0.05:
            star_label = "*"

        te_color = colors[i]

        ax.bar(x[j] + offsets[i], observed_val, width=bar_width, label=te if j == 0 else None, color=te_color)
        if observed_val > expected[j] and star_label != "":
            ax.text(x[j] + offsets[i], observed_val + (max(observed_values) / 50), star_label, ha='center', va='bottom', fontsize=15, color='black')

# Labels and formatting
ax.set_xticks(x)
ax.set_xticklabels(categories, rotation=45, ha="right", fontweight="bold", fontsize=16)
ax.set_ylabel("Count", fontweight="bold", fontsize=16)
ax.set_xlabel("Feature Type", fontweight="bold", fontsize=16)
ax.set_title("TEs overlaping genomic features", fontweight="bold", fontsize=16)
ax.legend(title="TE Family", title_fontsize=12, prop={'weight': 'bold', 'size':10})

# Make tick labels bold
plt.xticks(fontweight='bold', fontsize=16)
plt.yticks(fontweight='bold', fontsize=16)

# Set log scale for Y-axis
ax.set_yscale("log")

# Display main plot
plt.tight_layout()
plt.savefig("te_family_counts_pbear25_features_CHISQ.png", dpi=600)
plt.show()
