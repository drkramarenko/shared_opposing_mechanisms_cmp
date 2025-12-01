#!/usr/bin/env Rscript

## shared_effect_meta_analysis_DCM_HCM.r
##
## Post-process MTAG fixed-effects meta-analysis across DCM and HCM:
##  - Attach DCM and HCM effect sizes
##  - Filter on extreme SE ratios
##  - Run random-effects meta-analysis for suggestive variants
##  - Choose the more conservative P-value (fixed vs random)
##  - Map back to build 37 coordinates and optionally exclude MYBPC3 region

suppressPackageStartupMessages({
  library(data.table)
  library(tidyr)
  library(meta)
})

project_dir <- "/path/to/project"

## ---- 1. Load MTAG fixed-effects meta results ----
mtag_file <- file.path(
  project_dir,
  "results",
  "shared_effect_meta_analysis_DCM_HCM_mtag_meta.txt"
)

dat <- fread(mtag_file, data.table = FALSE)
dat <- dat[order(dat$mtag_pval), ]

## ---- 2. Attach DCM and HCM summary statistics ----
# DCM: harmonized_dcm.tsv.gz
dcm <- fread(
  file.path(project_dir, "sum_stats", "harmonized_dcm.tsv.gz"),
  data.table = FALSE,
  select = c("rsID", "EA", "NEA", "BETA", "SE", "P")
)
colnames(dcm) <- c("SNP", "EA", "NEA", "BETA_DCM", "SE_DCM", "P_DCM")

dcm1 <- dcm2 <- dcm
dcm1$SNP_EA <- paste0(dcm1$SNP, "_", dcm1$EA)
dcm2$SNP_EA <- paste0(dcm2$SNP, "_", dcm2$NEA)
dcm2$BETA_DCM <- -1 * dcm2$BETA_DCM
dcm_long <- rbind(dcm1, dcm2)[, c("SNP_EA", "BETA_DCM", "SE_DCM", "P_DCM")]

# HCM: harmonized_hcm.tsv.gz
hcm <- fread(
  file.path(project_dir, "sum_stats", "harmonized_hcm.tsv.gz"),
  data.table = FALSE,
  select = c("rsid", "effect_allele", "noneffect_allele", "beta", "se", "pvalue")
)
colnames(hcm) <- c("SNP", "EA", "NEA", "BETA_HCM", "SE_HCM", "P_HCM")

hcm1 <- hcm2 <- hcm
hcm1$SNP_EA <- paste0(hcm1$SNP, "_", hcm1$EA)
hcm2$SNP_EA <- paste0(hcm2$SNP, "_", hcm2$NEA)
hcm2$BETA_HCM <- -1 * hcm2$BETA_HCM
hcm_long <- rbind(hcm1, hcm2)[, c("SNP_EA", "BETA_HCM", "SE_HCM", "P_HCM")]

# Attach to MTAG results
dat$SNP_EA <- paste0(dat$SNP, "_", dat$A1)
dat <- merge(dat, dcm_long, by = "SNP_EA", all.x = TRUE, all.y = FALSE)
dat <- merge(dat, hcm_long, by = "SNP_EA", all.x = TRUE, all.y = FALSE)

message("N with missing DCM betas: ", sum(is.na(dat$BETA_DCM)))
message("N with missing HCM betas: ", sum(is.na(dat$BETA_HCM)))
message("Total rows after merge    : ", nrow(dat))

## ---- 3. Filter on SE ratios (precision outliers) ----
dat$SE_ratio <- dat$SE_DCM / dat$SE_HCM

# Remove outliers where one trait dominates due to extreme SE differences
idx_rm <- which(dat$SE_ratio < 0.65 | dat$SE_ratio > 1.2)
if (length(idx_rm) > 0) {
  dat <- dat[-idx_rm, ]
}

message("Rows retained after SE-ratio filter: ", nrow(dat))

dat <- dat[order(dat$mtag_pval), ]

## ---- 4. (Optional) sanity checks on direction of effects ----
dat$Directions <- paste0(
  ifelse(dat$BETA_DCM >= 0, "+", "-"),
  ifelse(dat$BETA_HCM >= 0, "+", "-")
)

## ---- 5. Save QC’d fixed-effects meta results ----
out_full  <- file.path(
  project_dir,
  "results",
  "shared_effect_meta_analysis_DCM_HCM_mtag_meta_QCfull.txt"
)
out_short <- file.path(
  project_dir,
  "results",
  "shared_effect_meta_analysis_DCM_HCM_mtag_meta_QCshort.txt"
)

write.table(
  dat,
  file = out_full,
  col.names = TRUE,
  row.names = FALSE,
  quote     = FALSE,
  sep       = "\t"
)

write.table(
  dat[, c("SNP", "CHR", "BP", "A1", "A2",
          "meta_freq", "mtag_beta", "mtag_se", "mtag_pval")],
  file = out_short,
  col.names = TRUE,
  row.names = FALSE,
  quote     = FALSE,
  sep       = "\t"
)

system(paste("gzip -f", out_full))
system(paste("gzip -f", out_short))

## ---- 6. Random-effects meta-analysis for suggestive variants ----
# Reload QCfull (gzipped)
dat <- fread(paste0(out_full, ".gz"), data.table = FALSE)
message("Rows in QCfull for random-effects step: ", nrow(dat))

# Run random-effects meta for variants with P < 1e-4 in MTAG
idx_re <- which(dat$mtag_pval < 1e-4)

if (length(idx_re) > 0) {
  dat$metagen_res <- NA_character_

  dat[idx_re, "metagen_res"] <- apply(
    X = dat[idx_re, c("BETA_DCM", "SE_DCM", "BETA_HCM", "SE_HCM")],
    MARGIN = 1,
    FUN = function(line) {
      TEs   <- c(line[1], line[3])
      seTEs <- c(line[2], line[4])
      met   <- metagen(TE = TEs, seTE = seTEs, comb.fixed = FALSE, comb.random = TRUE)
      Q     <- met$Q
      p.Q   <- met$pval.Q
      beta  <- met$TE.random
      se    <- met$seTE.random
      p     <- met$pval.random
      paste0(c(Q, p.Q, beta, se, p), collapse = ";")
    }
  )

  # Split random-effects results into columns
  dat_re  <- dat[idx_re, ]
  dat_re  <- tidyr::separate(
    dat_re,
    col  = "metagen_res",
    into = c("Q", "Q.pvalue", "random_beta", "random_se", "random_pval"),
    sep  = ";"
  )

  dat_no_re <- dat[-idx_re, ]
  dat_no_re$metagen_res <- NULL

  dat_no_re[, c("Q", "Q.pvalue", "random_beta", "random_se", "random_pval")] <- NA

  dat <- rbind(dat_re, dat_no_re)
  dat <- dat[order(dat$CHR, dat$BP), ]

  # Convert to numeric
  numeric_cols <- c("Q", "Q.pvalue", "random_beta", "random_se", "random_pval")
  dat[numeric_cols] <- lapply(dat[numeric_cols], as.numeric)
} else {
  # If no variants pass threshold, add empty random-effects columns
  dat[, c("Q", "Q.pvalue", "random_beta", "random_se", "random_pval")] <- NA
}

message("Rows after random-effects step: ", nrow(dat))

## ---- 7. Choose harmonized (more conservative) P-values ----
dat$harmonized_beta  <- dat$mtag_beta
dat$harmonized_se    <- dat$mtag_se
dat$harmonized_pval  <- dat$mtag_pval

idx_use_random <- which(!is.na(dat$random_pval) & dat$random_pval > dat$mtag_pval)

if (length(idx_use_random) > 0) {
  dat[idx_use_random, c("harmonized_beta", "harmonized_se", "harmonized_pval")] <-
    dat[idx_use_random, c("random_beta", "random_se", "random_pval")]
}

## ---- 8. Map back to build 37 coordinates using harmonized DCM ----
dcm_coord <- fread(
  file.path(project_dir, "sum_stats", "harmonized_dcm.tsv.gz"),
  data.table = FALSE,
  select = c("rsID", "CHRBP_B37")
)

dcm_coord <- dcm_coord[!duplicated(dcm_coord$rsID), ]

dat_new <- merge(dcm_coord, dat, by.x = "rsID", by.y = "SNP", all = FALSE)
message("Rows after merging with build 37 coords: ", nrow(dat_new))

dat_new <- tidyr::separate(
  dat_new,
  col  = "CHRBP_B37",
  into = c("CHR_B37", "BP_B37"),
  sep  = ":"
)

dat_new$CHR_B37 <- as.integer(dat_new$CHR_B37)
dat_new$BP_B37  <- as.integer(dat_new$BP_B37)
dat_new <- dat_new[order(dat_new$CHR_B37, dat_new$BP_B37), ]

## ---- 9. Write final outputs ----
out_full2  <- file.path(
  project_dir,
  "results",
  "shared_effect_meta_analysis_DCM_HCM_mtagandrandom_meta_QCfull.txt"
)
out_short2 <- file.path(
  project_dir,
  "results",
  "shared_effect_meta_analysis_DCM_HCM_mtag_meta_QCshort.txt"
)
out_harm   <- file.path(
  project_dir,
  "results",
  "shared_effect_meta_analysis_DCM_HCM_harmonized_meta_QCshort.txt"
)

write.table(
  dat_new,
  file = out_full2,
  col.names = TRUE,
  row.names = FALSE,
  quote     = FALSE,
  sep       = "\t"
)

write.table(
  dat_new[, c("rsID", "CHR_B37", "BP_B37", "A1", "A2",
              "meta_freq", "mtag_beta", "mtag_se", "mtag_pval")],
  file = out_short2,
  col.names = TRUE,
  row.names = FALSE,
  quote     = FALSE,
  sep       = "\t"
)

write.table(
  dat_new[, c("rsID", "CHR_B37", "BP_B37", "A1", "A2",
              "meta_freq", "harmonized_beta", "harmonized_se", "harmonized_pval")],
  file = out_harm,
  col.names = TRUE,
  row.names = FALSE,
  quote     = FALSE,
  sep       = "\t"
)

system(paste("gzip -f", out_full2))
system(paste("gzip -f", out_short2))
system(paste("gzip -f", out_harm))

## ---- 10. (Optional) exclude MYBPC3 region (chr11:30–80 Mb) ----
dat_new <- fread(paste0(out_full2, ".gz"), data.table = FALSE)

idx_keep <- !(dat_new$CHR_B37 == 11 &
              dat_new$BP_B37 >= 30000000 &
              dat_new$BP_B37 <= 80000000)

dat_new <- dat_new[idx_keep, ]
message("Rows after excluding MYBPC3 region: ", nrow(dat_new))

out_full2_excl  <- file.path(
  project_dir,
  "results",
  "shared_effect_meta_analysis_DCM_HCM_mtagandrandom_meta_QCfull_exclMYBPC3reg.txt"
)
out_short2_excl <- file.path(
  project_dir,
  "results",
  "shared_effect_meta_analysis_DCM_HCM_mtag_meta_QCshort_exclMYBPC3reg.txt"
)
out_harm_excl   <- file.path(
  project_dir,
  "results",
  "shared_effect_meta_analysis_DCM_HCM_harmonized_meta_QCshort_exclMYBPC3reg.txt"
)

write.table(
  dat_new,
  file = out_full2_excl,
  col.names = TRUE,
  row.names = FALSE,
  quote     = FALSE,
  sep       = "\t"
)

write.table(
  dat_new[, c("rsID", "CHR_B37", "BP_B37", "A1", "A2",
              "meta_freq", "mtag_beta", "mtag_se", "mtag_pval")],
  file = out_short2_excl,
  col.names = TRUE,
  row.names = FALSE,
  quote     = FALSE,
  sep       = "\t"
)

write.table(
  dat_new[, c("rsID", "CHR_B37", "BP_B37", "A1", "A2",
              "meta_freq", "harmonized_beta", "harmonized_se", "harmonized_pval")],
  file = out_harm_excl,
  col.names = TRUE,
  row.names = FALSE,
  quote     = FALSE,
  sep       = "\t"
)

system(paste("gzip -f", out_full2_excl))
system(paste("gzip -f", out_short2_excl))
system(paste("gzip -f", out_harm_excl))

message("Shared-effects meta-analysis completed.")
