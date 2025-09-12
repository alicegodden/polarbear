# Title : PCA plotting from multibam summary output
# Author : Dr. Alice M. Godden

# Load libraries
library(ggplot2)
library(pheatmap)
library(matrixStats)

# ---- Load data ----
mat <- read.csv("results_10kb_counts.csv", check.names=FALSE)

# Sample names and groups
sample_names <- c("SAMN16454142","SAMN16454143","SAMN16454144","SAMN16454145",
                  "SAMN16454146","SAMN16454147","SAMN16454148","SAMN16454149",
                  "SAMN16454150","SAMN16454151","SAMN16454152","SAMN16454153",
                  "SAMN16454154","SAMN16454157","SAMN16454158","SAMN16454159",
                  "SAMN16454160")
group_labels <- c("NEG","NEG","SEG","NEG","NEG","SEG","SEG","SEG","SEG","SEG",
                  "SEG","SEG","SEG","SEG","NEG","NEG","NEG")

colnames(mat) <- sample_names

# Ensure numeric matrix
mat_numeric <- as.matrix(mat)
mode(mat_numeric) <- "numeric"

# Remove zero-variance rows
mat_numeric <- mat_numeric[apply(mat_numeric, 1, var) != 0, ]

# ---- PCA ----
pca <- prcomp(t(mat_numeric), scale.=TRUE)  # transpose: rows = samples
pca_var_perc <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)

pca_df <- data.frame(pca$x, Group = factor(group_labels), Sample = rownames(pca$x))

# PCA plot
ggplot(pca_df, aes(x=PC1, y=PC2, color=Group)) +
  geom_point(size=4) +
  stat_ellipse(level=0.95, aes(fill=Group), geom="polygon", alpha=0.2, color=NA) +
  xlab(paste0("PC1 (", pca_var_perc[1], "%)")) +
  ylab(paste0("PC2 (", pca_var_perc[2], "%)")) +
  theme_classic() +
  theme(text=element_text(face="bold", size=12),
        legend.title=element_text(face="bold"),
        legend.text=element_text(face="bold")) +
  scale_color_manual(values=c("NEG"="#1f78b4","SEG"="orange")) +
  scale_fill_manual(values=c("NEG"="#1f78b4","SEG"="orange"))

