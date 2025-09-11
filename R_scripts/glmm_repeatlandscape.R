# Title: Negative binomial GLMM Repeatlandscape
# Subtitle: Is there a significant shift in TE age in NEG v SEG bears?
# Author: Dr. Alice M. Godden

# Title: Plotting repeat landscape from RepeatMasker divsum files
# Author: Gemini, based on Dr. Alice M. Godden's original script

# Load necessary libraries
# This script requires ggplot2, dplyr, readr, tidyr, and ggsci
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(ggsci)

# ---
# Step 1: Generate trimmed divsum files
# This section takes the raw divsum files and removes the header.

# Define the directory containing the original divsum files
input_dir <- "~/Desktop/IMMLER_PAPERS/polarbear/newgenome/denovo_assembly_ASMv1lib_repmask/divs"
output_dir <- "~/Desktop/IMMLER_PAPERS/polarbear/newgenome/denovo_assembly_ASMv1lib_repmask/divs/cleaned_divsum"

# Create the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Get a list of all divsum files in the directory
divsum_files <- list.files(input_dir, pattern = "\\.divsum$", full.names = TRUE)

# Function to clean a divsum file by removing the header
clean_divsum <- function(file_path, output_dir) {
  # Read the file
  lines <- readLines(file_path)
  
  # Find the line index where "Coverage for each repeat class and divergence (Kimura)" occurs
  start_index <- grep("Coverage for each repeat class and divergence \\(Kimura\\)", lines)
  
  # If the line is found, keep only the lines after it
  if (length(start_index) > 0) {
    cleaned_lines <- lines[(start_index + 1):length(lines)]
    
    # Write to a new file in the cleaned divsum directory
    output_file <- file.path(output_dir, basename(file_path))
    writeLines(cleaned_lines, output_file)
    
    cat(paste("Cleaned file saved:", output_file, "\n"))
  } else {
    cat(paste("Skipping file (header not found):", file_path, "\n"))
  }
}

# Apply the cleaning function to all divsum files
for (file in divsum_files) {
  clean_divsum(file, output_dir)
}

cat("Cleaning process completed!\n")

# ---
# Step 2: Read, process, and combine cleaned data
# This section reads the cleaned files, combines them into one dataframe,
# and assigns the TE families and sample conditions.

# Define directories for cleaned files
dir_neg <- "~/Desktop/IMMLER_PAPERS/polarbear/newgenome/denovo_assembly_ASMv1lib_repmask/divs/cleaned_divsum/cleaned_NEG"
dir_seg <- "~/Desktop/IMMLER_PAPERS/polarbear/newgenome/denovo_assembly_ASMv1lib_repmask/divs/cleaned_divsum/cleaned_SEG"

# List all cleaned divsum files
files_neg <- list.files(dir_neg, pattern = "\\.divsum$", full.names = TRUE)
files_seg <- list.files(dir_seg, pattern = "\\.divsum$", full.names = TRUE)

# Stop if no files found
if (length(files_neg) == 0 | length(files_seg) == 0) {
  stop("No .divsum files found in one of the directories!")
}

# Define genome size (adjust as needed for your analysis)
genome_size <- 146016386 # bear ASM1731132v1 assembly length in bp

# Function to process a single file and assign condition
process_divsum <- function(file, condition_label) {
  df <- read_delim(file, delim = " ", col_types = cols(), trim_ws = TRUE)
  sample_name <- gsub(".*/|\\.divsum$", "", file)
  df$Sample <- sample_name
  df$Condition <- condition_label
  
  df_long <- pivot_longer(df, cols = -c(Div, Sample, Condition), names_to = "Repeat_Type", values_to = "Count")
  df_long <- df_long %>% mutate(Normalized = (Count / genome_size) * 100)
  return(df_long)
}

# Read all files and assign condition
df_neg <- bind_rows(lapply(files_neg, process_divsum, condition_label = "NEG"))
df_seg <- bind_rows(lapply(files_seg, process_divsum, condition_label = "SEG"))

# Combine the data frames
plot_df <- bind_rows(df_neg, df_seg)

# Assign TE family groups
plot_df <- plot_df %>% 
  mutate(TE_Family = case_when(
    grepl("DNA", Repeat_Type) ~ "DNA",
    grepl("LINE", Repeat_Type) ~ "LINE",
    grepl("LTR", Repeat_Type) ~ "LTR",
    grepl("SINE", Repeat_Type) ~ "SINE",
    grepl("RC", Repeat_Type) ~ "RC",
    grepl("Satellite", Repeat_Type) ~ "Satellite",
    grepl("Retrotransposon", Repeat_Type) ~ "Retrotransposon",
    TRUE ~ "Other" # Label remaining types as "Other"
  )) %>% 
  filter(TE_Family != "Other") # Remove "Other" as per your plot code

# ---
# Step 3: Generate the final plot
# This section uses the prepared data to create the plot.

# Generate the ggplot visualization
final_plot <- ggplot(plot_df, aes(x = Div, y = Count, color = Condition)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "NULL") +
  facet_wrap(~ TE_Family, scales = "free_y") +
  theme_classic() +
  scale_color_manual(values = c("NEG" = "blue", "SEG" = "orange")) +
  coord_cartesian(xlim = c(0, 50)) +
  labs(title = "Kimura substitution - NEG v SEG",
       x = "Kimura Substitution level",
       y = "TE abundance",
       color = "Condition") +
  theme(
    text = element_text(face = "bold", color = "black"),
    strip.text = element_text(face = "bold", color = "black"),
    axis.title = element_text(face = "bold", color = "black"),
    axis.text = element_text(face = "bold", color = "black", size = 14)
  )

# Print the plot
print(final_plot)

# To save the plot to a file, you can uncomment and adjust the following line:
ggsave("kimura_plot_noloess.png", final_plot, width = 10, height = 8, dpi = 300)

# stats analysis 

# GLMMM
# ------------------------------
# GLMM Pipeline for TE Age Shift
# ------------------------------

library(dplyr)
library(lme4)
library(broom.mixed)
library(AER)  # for overdispersion check

# ------------------------------
# Step 1: Prepare Data
# ------------------------------

stats_df <- plot_df %>%
  dplyr::select(Div, Sample, Condition, Repeat_Type, TE_Family, Count, Normalized) %>%
  filter(!is.na(TE_Family))

stats_df$Condition <- factor(stats_df$Condition, levels = c("NEG", "SEG"))
stats_df$TE_Family <- factor(stats_df$TE_Family)
stats_df$Sample <- factor(stats_df$Sample)

response_var <- "Count"  # can switch to "Normalized" if needed

# ------------------------------
# Step 2: Fit GLMMs per TE Family
# ------------------------------

summary_table <- tibble()
model_list <- list()

for (te in unique(stats_df$TE_Family)) {
  
  cat("\n📌 TE Family:", te, "\n")
  
  te_data <- stats_df %>% filter(TE_Family == te)
  formula_glmm <- as.formula(paste(response_var, "~ Condition * Div + (1|Sample)"))
  
  # Poisson GLMM
  poisson_glmm <- tryCatch({
    glmer(formula_glmm, data=te_data, family="poisson",
          control=glmerControl(optimizer="bobyqa"))
  }, error=function(e) NULL)
  
  # Negative Binomial GLMM
  nb_glmm <- tryCatch({
    glmer.nb(formula_glmm, data=te_data,
             control=glmerControl(optimizer="bobyqa"))
  }, error=function(e) NULL)
  
  # Overdispersion for Poisson
  overdisp <- if(!is.null(poisson_glmm)){
    rdf <- df.residual(poisson_glmm)
    res <- residuals(poisson_glmm, type="pearson")
    chi_sq <- sum(res^2)
    disp_ratio <- chi_sq / rdf
    p_val <- pchisq(chi_sq, df=rdf, lower.tail=FALSE)
    list(dispersion=disp_ratio, p_value=p_val)
  } else list(dispersion=NA_real_, p_value=NA_real_)
  
  # Extract NB fixed effects
  if(!is.null(nb_glmm)){
    coef_sum <- tidy(nb_glmm, effects="fixed")
    
    # Extract p-values and estimates
    p_cond <- coef_sum$p.value[coef_sum$term=="ConditionSEG"]
    p_div <- coef_sum$p.value[coef_sum$term=="Div"]
    p_inter <- coef_sum$p.value[coef_sum$term=="ConditionSEG:Div"]
    
    est_cond <- coef_sum$estimate[coef_sum$term=="ConditionSEG"]
    est_div <- coef_sum$estimate[coef_sum$term=="Div"]
    est_inter <- coef_sum$estimate[coef_sum$term=="ConditionSEG:Div"]
    
    # Significance function
    star <- function(p){
      if(is.na(p)) return(NA_character_)
      else if(p<0.001) return("***")
      else if(p<0.01) return("**")
      else if(p<0.05) return("*")
      else return("ns")
    }
    
    summary_table <- bind_rows(summary_table, tibble(
      TE_Family = te,
      Poisson_AIC = if(!is.null(poisson_glmm)) AIC(poisson_glmm) else NA_real_,
      NB_AIC = AIC(nb_glmm),
      Overdispersion_ratio = overdisp$dispersion,
      Overdispersion_p = overdisp$p_value,
      Est_Condition = est_cond,
      P_Condition = p_cond,
      Sig_Condition = star(p_cond),
      Est_Div = est_div,
      P_Div = p_div,
      Sig_Div = star(p_div),
      Est_Interaction = est_inter,
      P_Interaction = p_inter,
      Sig_Interaction = star(p_inter)
    ))
    
    model_list[[te]] <- list(poisson_glmm=poisson_glmm, nb_glmm=nb_glmm)
    cat("✅ GLMMs fit successfully\n")
    
  } else {
    summary_table <- bind_rows(summary_table, tibble(
      TE_Family = te,
      Poisson_AIC = NA_real_,
      NB_AIC = NA_real_,
      Overdispersion_ratio = NA_real_,
      Overdispersion_p = NA_real_,
      Est_Condition = NA_real_,
      P_Condition = NA_real_,
      Sig_Condition = NA_character_,
      Est_Div = NA_real_,
      P_Div = NA_real_,
      Sig_Div = NA_character_,
      Est_Interaction = NA_real_,
      P_Interaction = NA_real_,
      Sig_Interaction = NA_character_
    ))
    cat("❌ Failed to fit GLMMs for", te, "\n")
  }
}

# ------------------------------
# Step 3: Save and Print
# ------------------------------

summary_table <- summary_table %>% arrange(P_Interaction)
print(summary_table)

write.csv(summary_table, "TE_Family_GLMM_Interaction_Results.csv", row.names = FALSE)
cat("\n✅ GLMM analysis with interaction complete! Results saved to 'TE_Family_GLMM_Interaction_Results.csv'.\n")



library(ggplot2)
library(dplyr)

# ------------------------------
# Step 1: Prepare significance labels
# ------------------------------

# Replace with your GLMM results
glmm_results <- data.frame(
  TE_Family = c("LINE", "LTR"),
  Sig_Interaction = c("***", "*")
)

# ------------------------------
# Step 2: Summarize counts for plotting
# ------------------------------

plot_df_summary <- plot_df %>%
  filter(TE_Family %in% glmm_results$TE_Family) %>%
  group_by(TE_Family, Condition, Div) %>%
  summarise(
    Mean_Count = mean(Count, na.rm=TRUE),
    SD_Count = sd(Count, na.rm=TRUE),
    .groups="drop"
  )

# Merge significance labels
plot_df_summary <- plot_df_summary %>%
  left_join(glmm_results, by="TE_Family")

# ------------------------------
# Step 3: Plot
# ------------------------------

ggplot(plot_df_summary, aes(x=Div, y=Mean_Count, color=Condition)) +
  geom_line(size=1) +
  geom_ribbon(aes(ymin=Mean_Count-SD_Count, ymax=Mean_Count+SD_Count, fill=Condition),
              alpha=0.2, color=NA) +
  facet_wrap(~TE_Family, scales="free_y") +
  scale_color_manual(values=c("NEG"="blue", "SEG"="orange")) +
  scale_fill_manual(values=c("NEG"="blue", "SEG"="orange")) +
  labs(
    title="TE Abundance vs Kimura Divergence",
    subtitle="Interaction indicates TEs are younger in SEG",
    x="Kimura Divergence (%)",
    y="Mean TE Count",
    color="Condition",
    fill="Condition"
  ) +
  theme_classic() +
  theme(
    text=element_text(face="bold", color="black"),
    axis.text=element_text(face="bold", color="black", size=12),
    strip.text=element_text(face="bold", size=14)
  ) +
  # Add interaction significance labels at top of panel
  geom_text(
    data = glmm_results,
    aes(x=Inf, y=Inf, label=Sig_Interaction),
    inherit.aes=FALSE,
    hjust=1.1, vjust=1.5,
    size=6,
    fontface="bold",
    color="red"
  )

