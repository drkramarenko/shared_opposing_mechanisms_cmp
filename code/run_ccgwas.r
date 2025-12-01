#!/usr/bin/env Rscript

##
## CC-GWAS for DCM vs HCM
## Clean, reproducible script using harmonized summary statistics
##

suppressPackageStartupMessages({
  library(data.table)
  library(R.utils)
  library(CCGWAS)
})

# ----------------------------------------------------------
# 1. User-defined paths (EDIT for your environment)
# ----------------------------------------------------------

sumstats_dcm <- "sum_stats/harmonized_dcm.tsv.gz"
sumstats_hcm <- "sum_stats/harmonized_hcm.tsv.gz"
out_prefix   <- "CC_GWAS__DCM__HCM"

# ----------------------------------------------------------
# 2. CC-GWAS parameters (pre-computed from LDSC)
# ----------------------------------------------------------

CCGWAS(
  outcome_file = paste0(out_prefix, ".out"),
  A_name = "DCM",
  B_name = "HCM",

  sumstats_fileA1A0 = sumstats_dcm,
  sumstats_fileB1B0 = sumstats_hcm,

  ## Population prevalences / liability-scale parameters
  K_A1A0      = 0.004,
  K_A1A0_high = 0.01,
  K_A1A0_low  = 0.002,

  K_B1B0      = 0.002,
  K_B1B0_high = 0.005,
  K_B1B0_low  = 0.001,

  ## LDSC-derived parameters
  h2l_A1A0            = 0.142,
  h2l_B1B0            = 0.1798,
  rg_A1A0_B1B0        = -0.5642,
  intercept_A1A0_B1B0 = 0.0107,

  ## Number of causal SNPs (recommended default)
  m = 1200,

  ## Effective sample sizes
  N_A1 = 9365 * 0.9,
  N_B1 = 5900 * 0.85,
  N_A0 = 946368 * 0.9,
  N_B0 = 68359 * 0.85,
  N_overlap_A0B0 = 40283 * 0.85
)

# ----------------------------------------------------------
# 3. Load CC-GWAS output
# ----------------------------------------------------------

cc <- fread(paste0(out_prefix, ".out.results.gz"), data.table = FALSE)

message("Genome-wide significant CC-GWAS loci (OLS_pval < 5e-8 & CCGWAS_signif = 1): ",
        sum(cc$OLS_pval < 5e-8 & cc$CCGWAS_signif == 1))

# ----------------------------------------------------------
# 4. Sort by chromosome and base-pair position
# ----------------------------------------------------------

cc <- cc[order(cc$CHR, cc$BP), ]

# ----------------------------------------------------------
# 5. Remove MYBPC3 region
# ----------------------------------------------------------

mybpc3_chr   <- 11
mybpc3_start <- 30000000
mybpc3_stop  <- 80000000

message("Original number of SNPs: ", nrow(cc))

cc_filtered <- cc[
  !(CHR == mybpc3_chr & BP >= mybpc3_start & BP <= mybpc3_stop),
]

message("After removing MYBPC3 region: ", nrow(cc_filtered))

# ----------------------------------------------------------
# 6. Write filtered output
# ----------------------------------------------------------

outfile_txt <- paste0(out_prefix, "_exclMYBPC3reg.txt")

write.table(
  cc_filtered,
  file = outfile_txt,
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

R.utils::gzip(outfile_txt, overwrite = TRUE)

message("Filtered CC-GWAS output written to: ", outfile_txt, ".gz")
