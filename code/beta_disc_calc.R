############################################################
# Compare discovery and validation effect sizes
# Discovery  = CC-GWAS results
# Validation = external HCM/DCM sumstats
#
# Steps:
# 1. Read discovery and validation datasets
# 2. Merge by rsID
# 3. Align validation beta to discovery effect allele
# 4. Flip all validation betas
# 5. For SNPs with negative discovery beta, flip both axes
#    so all discovery betas are positive on x-axis
# 6. Run correlation and regression
# 7. Plot and save PNG
############################################################

library(data.table)
library(dplyr)
library(ggplot2)

#-----------------------------------------------------------
# 0. File paths (edit to point to your own files)
#-----------------------------------------------------------
discovery_file  <- "path/to/discovery_ccgwas_results_with_scaled_betas.tsv"
validation_file <- "path/to/validation_dcm_hcm_sumstats.tsv.gz"
top_snps_file   <- "path/to/top_snps.txt"

#-----------------------------------------------------------
# 1. Read discovery results
#-----------------------------------------------------------
discovery <- fread(discovery_file)

# Keep relevant columns
discovery_sub <- discovery %>%
  select(
    SNP,
    EA,
    NEA,
    OLS_beta,
    OLS_se,
    BETA_DCMvHCM_OLSbackcomputed,
    SE_DCMvHCM_OLSbackcomputed,
    BETA_DCMvHCM_exactbackcomputed,
    SE_DCMvHCM_exactbackcomputed
  )

#-----------------------------------------------------------
# 2. Read validation sumstats
#-----------------------------------------------------------
validation <- fread(validation_file)

# Remove rows without rsID
validation <- validation %>%
  filter(rsID != ".")

# Keep relevant columns
validation_sub <- validation %>%
  select(
    rsID,
    A1,
    A2,
    BETA,
    P  )

#-----------------------------------------------------------
# 3. Read top SNP list
#-----------------------------------------------------------
top_snps <- fread(top_snps_file)

# Use first column as SNP list
top_snp_vector <- top_snps[[1]]

#-----------------------------------------------------------
# 4. Merge discovery and validation by SNP / rsID
#-----------------------------------------------------------
df <- discovery_sub %>%
  inner_join(validation_sub, by = c("SNP" = "rsID"))

cat("Number of overlapping SNPs before filtering:", nrow(df), "\n")

#-----------------------------------------------------------
# 5. Harmonize alleles and prepare comparison variables
#-----------------------------------------------------------
df <- df %>%
  mutate(
    # Align validation beta to discovery effect allele
    allele_flip = case_when(
      EA == A1 ~ 1,
      EA == A2 ~ -1,
      TRUE ~ NA_real_
    ),
    validation_beta_aligned = BETA * allele_flip
  ) %>%

  # Keep only harmonized SNPs, remove extreme validation beta values,
  # and restrict to top discovery SNPs
  filter(
    !is.na(allele_flip),
    SNP %in% top_snp_vector
  ) %>%

  mutate(
    # Use OLS backcomputed beta as discovery effect
    discovery_beta_raw = BETA_DCMvHCM_OLSbackcomputed,

    # Flip all validation betas globally
    validation_beta_globally_flipped = -validation_beta_aligned,

    # If discovery beta is negative, flip both axes
    discovery_beta_plot = if_else(discovery_beta_raw < 0,
                                  -discovery_beta_raw,
                                  discovery_beta_raw),

    validation_beta_plot = if_else(discovery_beta_raw < 0,
                                   -validation_beta_globally_flipped,
                                   validation_beta_globally_flipped)
  )

cat("Number of SNPs after harmonization/filtering:", nrow(df), "\n")

#-----------------------------------------------------------
# 6. Inspect prepared data
#-----------------------------------------------------------
print(
  df %>%
    select(
      SNP, EA, NEA, A1, A2,
      BETA_DCMvHCM_OLSbackcomputed,
      BETA,
      validation_beta_aligned,
      validation_beta_globally_flipped,
      discovery_beta_plot,
      validation_beta_plot
    ) %>%
    head()
)

#-----------------------------------------------------------
# 7. Correlation
#-----------------------------------------------------------
cor_val <- cor(
  df$discovery_beta_plot,
  df$validation_beta_plot,
  use = "complete.obs"
)

cat("Correlation:", cor_val, "\n")

#-----------------------------------------------------------
# 8. Linear regression
#-----------------------------------------------------------
model <- lm(validation_beta_plot ~ discovery_beta_plot, data = df)
model_summary <- summary(model)

print(model_summary)

# Call:
# lm(formula = validation_beta_plot ~ discovery_beta_plot, data = df)

# Residuals:
#       Min        1Q    Median        3Q       Max
# -0.277078 -0.074677  0.008256  0.060232  0.295518

# Coefficients:
#                     Estimate Std. Error t value Pr(>|t|)
# (Intercept)         0.008146   0.019580   0.416    0.678
# discovery_beta_plot 0.962936   0.083774  11.495   <2e-16 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Residual standard error: 0.09654 on 122 degrees of freedom
# Multiple R-squared:  0.5199,    Adjusted R-squared:  0.516
# F-statistic: 132.1 on 1 and 122 DF,  p-value: < 2.2e-16

#-----------------------------------------------------------
# 9. Create plot
#-----------------------------------------------------------
p <- ggplot(df, aes(
  x = discovery_beta_plot,
  y = validation_beta_plot
)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Discovery beta (CC-GWAS OLS backcomputed, sign-normalized)",
    y = "Validation beta (allele-aligned and sign-normalized)",
    title = "Discovery vs validation effect sizes"
  )

print(p)

#-----------------------------------------------------------
# 10. Save plot
#-----------------------------------------------------------
ggsave(
  filename = "discovery_vs_validation_OLS_backcomputed.png",
  plot = p,
  width = 6,
  height = 6,
  dpi = 300
)

#-----------------------------------------------------------
# 11. Optional: save merged results
#-----------------------------------------------------------
fwrite(df, "discovery_validation_merged_topSNPs.tsv", sep = "\t")

#-----------------------------------------------------------
# 12. Print compact summary for reporting
#-----------------------------------------------------------
cat("\n================ MODEL SUMMARY ================\n")
cat("Number of SNPs:", nrow(df), "\n")
cat("Correlation:", round(cor_val, 4), "\n")
cat("Intercept:", round(coef(model)[1], 6), "\n")
cat("Slope:", round(coef(model)[2], 6), "\n")
cat("R-squared:", round(model_summary$r.squared, 4), "\n")
cat("Adj. R-squared:", round(model_summary$adj.r.squared, 4), "\n")
cat("P-value (slope):", signif(coef(model_summary)[2, 4], 4), "\n")
cat("================================================\n")

# ================ MODEL SUMMARY ================
# Number of SNPs: 124
# Correlation: 0.7211
# Intercept: 0.008146
# Slope: 0.962936
# R-squared: 0.5199
# Adj. R-squared: 0.516
# P-value (slope): 3.605e-21
# ================================================