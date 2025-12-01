#!/usr/bin/env Rscript

# R code to create a nice forrest plot from the csv file above
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(gridExtra)
library(metafor)

# Read the data from the CSV file
df <- read.csv("data/sumtable_forestplot_CC_Varcorr.csv", sep = ",")
# ST14 PGS performance metrics across DCM, HCM & CC MTAG

# Replace "\n" in the Outcome column with actual line breaks
df$Outcome <- gsub("\\\\n", "\n", df$Outcome)

# Define colors for each outcome
colors <- c("DCM-MTAG" = "#7030A0", "HCM-MTAG" = "#00B050", "CC-MTAG" = "royalblue1")

# Prepare data for forest plot
df_forest <- df %>%
  mutate(
         outcome = factor(Outcome, levels = unique(Outcome)),
         prs = factor(PRS, levels = unique(PRS)),
         OR_CI = sprintf("%.2f [%.2f, %.2f]", OR.full, CI.L.OR.full, CI.H.OR.full), # OR [95% CI] with 2 digits
         P_value = toupper(format(P.full, scientific = TRUE, digits = 3)),        # P-value with 2 digits before the 'E'
         AUC = sprintf("%.3f", AUC.prs),                                         # AUC with 3 digits
         AUPRC = sprintf("%.3f", AUPRC.prs),                                     # AUPRC with 3 digits
         Ngk_R2 = sprintf("%.3f", ngk_R2_delta_o_resid),                         # Ngk-R2 with 3 digits
         Liab_R2 = ifelse(is.na(as.numeric(liab_R2_delta_o_resid)),"-",          # If value non-numeric, replace with "-"
         sprintf("%.3f", as.numeric(liab_R2_delta_o_resid))),                    # If numeric, with 3 digits
         color = colors[PRS]) %>%
  arrange(outcome, prs)

# Create a forest plot using metafor
forest_plot <- metafor::forest(
  x = df_forest$OR.full,                   # The point estimates (odds ratios) to be plotted.
  ci.lb = df_forest$CI.L.OR.full,          # The lower bounds of the confidence intervals.
  ci.ub = df_forest$CI.H.OR.full,          # The upper bounds of the confidence intervals.
  slab = NA,                               # Remove the default slab text.
  xlab = "Odds Ratio",                     # Label for the x-axis.
  alim = c(0.5, 3.22),                      # Limits for the x-axis (range of the odds ratios to be displayed).
  pch = 16,                                # Plotting character (symbol) for the points (16 is a filled circle).
  xlim = c(-0.9, 8.6),                    # Adjusted limits for the x-axis (including space for labels).
  at = c(1, 1.5, 2, 2.5, 3),          # Positions of the tick marks on the x-axis.
  refline = 1,                             # Position of the reference line (typically 1 for odds ratios).
  header = FALSE,                             # Header for the plot.
  ilab = cbind(df_forest$OR_CI, df_forest$P_value, df_forest$AUC, df_forest$AUPRC, df_forest$Ngk_R2, df_forest$Liab_R2),  # Additional labels to be displayed in the plot.
  ilab.lab = c("OR [95% CI]", "P-value", "AUC", "AUPRC", "Ngk-R2", "Liab-R2"),  # Labels for the additional labels.
  ilab.xpos = c(3.45, 4.8, 5.7, 6.4, 7.155, 7.9),            # Positions of the additional labels on the x-axis.
  ilab.pos = 4,                            # Position of the labels relative to the points (4 = on the right).
  annotate = FALSE,                        # Disable default annotation of OR[CI]
  col = df_forest$color,                   # Colors for the observed outcomes
  psize = 1                                # Set a constant size for the points
)

# Add dotted vertical lines at the tick marks, constrained to the plot area
tick_positions <- c(1, 1.5, 2, 2.5, 3)  # Positions of the tick marks

for (pos in tick_positions) {
  segments(
    x0 = pos, y0 = 0.01,  # Bottom (start) of the line
    x1 = pos, y1 = nrow(df_forest) + 0.99,  # Top (end) of the line
    col = "gray", lty = "dotted"  # Line color and type
  )
}

# Add column names "cohort" and "PGS" above the Outcome and PRS text
text(
  x = -0.9,
  y = nrow(df_forest) + 2,
  labels = "Cohort",
  pos = 4,
  cex = 1,
  col = "black",
  font = 2
)

text(
  x = 0,
  y = nrow(df_forest) + 2,
  labels = "PGS",
  pos = 4,
  cex = 1,
  col = "black",
  font = 2
)

# Add Outcome and PRS labels using text function
unique_outcomes <- unique(df_forest$outcome)
for (outcome in unique_outcomes) {
  prs_indices <- which(df_forest$outcome == outcome)
  mid_index <- median(prs_indices)
  
  # Add Outcome label
  text(
    x = -0.9,
    y = nrow(df_forest) - mid_index + 1,
    labels = as.character(outcome),
    pos = 4,
    cex = 0.95,
    col = "black"
  )
  
  # Add number of cases below the outcome text
  if (outcome == "DCM") {
    text(
      x = -0.9,
      y = nrow(df_forest) - mid_index + 0.4,
      labels = "(n=1053)",
      pos = 4,
      cex = 0.65,
      col = "black"
    )
  } else if (outcome == "HCM") {
    text(
      x = -0.9,
      y = nrow(df_forest) - mid_index + 0.4,
      labels = "(n=562)",
      pos = 4,
      cex = 0.65,
      col = "black"
    )
  }
  
  # Add PRS labels
  for (i in prs_indices) {
    text(
      x = 0,
      y = nrow(df_forest) - i + 1,
      labels = as.character(df_forest$prs[i]),
      pos = 4,
      cex = 0.95,
      col = "black"
    )
  }

  x_position = -0.1

  # Calculate the start and end positions for the vertical line
  y_start <- nrow(df_forest) - max(prs_indices) + 1  # Top of the block
  y_end <- nrow(df_forest) - min(prs_indices) + 1    # Bottom of the block
  
  # Add a single vertical line for the current outcome
  segments(
    x0 = x_position, y0 = y_start,  # Start of the line
    x1 = x_position, y1 = y_end,    # End of the line
    col = "black", lwd = 1  # Line color and width
  )
}

# Save as PDF with 9.19 x 6 inches
# OR PDF with 11.50 x 5.88 inches