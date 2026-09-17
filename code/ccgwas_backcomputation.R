# Approximate back-conversion of CC-GWAS Exact_beta to direct case-case logistic scale
# Inputs:
#   exact_beta : CC-GWAS Exact_beta (standardized observed 50/50 scale)
#   or_a       : odds ratio from disorder A case-control GWAS (A1 vs A0)
#   frq_a0     : effect-allele frequency in A controls (A0)
#   or_b       : odds ratio from disorder B case-control GWAS (B1 vs B0)
#   frq_b0     : effect-allele frequency in B controls (B0)
#
# Returns:
#   A list with:
#     logOR_A1B1 : approximate direct case-case log(OR) for A1 vs B1
#     OR_A1B1    : approximate direct case-case OR
#     p_A1       : inferred effect-allele frequency in A cases
#     p_B1       : inferred effect-allele frequency in B cases
#     p_casecase : pooled case-case allele frequency used in inversion
#     b_unscaled : unscaled linear-regression coefficient before re-conversion

ccgwas_exact_to_casecase_logor <- function(exact_beta, or_a, frq_a0, or_b, frq_b0) {
  # Basic checks
  stopifnot(length(exact_beta) == 1L || is.vector(exact_beta))
  stopifnot(all(or_a > 0, na.rm = TRUE))
  stopifnot(all(or_b > 0, na.rm = TRUE))
  stopifnot(all(frq_a0 > 0 & frq_a0 < 1, na.rm = TRUE))
  stopifnot(all(frq_b0 > 0 & frq_b0 < 1, na.rm = TRUE))

  # Reconstruct case allele frequencies from control AF and case-control OR:
  # p_case = OR * p_control / (1 - p_control + OR * p_control)
  p_A1 <- (or_a * frq_a0) / (1 - frq_a0 + or_a * frq_a0)
  p_B1 <- (or_b * frq_b0) / (1 - frq_b0 + or_b * frq_b0)

  # Pooled AF in the direct A1 vs B1 sample (assuming 50/50 A1:B1 sampling)
  p_casecase <- 0.5 * (p_A1 + p_B1)

  # Guard against numerical edge cases
  eps <- .Machine$double.eps^0.5
  p_casecase <- pmin(pmax(p_casecase, eps), 1 - eps)

  # Undo CC-GWAS scaling:
  # beta_scaled = 2 * b * sqrt(2 p (1-p))
  b_unscaled <- exact_beta / (2 * sqrt(2 * p_casecase * (1 - p_casecase)))

  # Invert the quadratic relation used in OR_to_scaled_linreg()
  A <- p_casecase * (1 - p_casecase)
  OR_A1B1 <- (A * b_unscaled^2 + 0.5 * b_unscaled + 0.25) /
    (A * b_unscaled^2 - 0.5 * b_unscaled + 0.25)

  logOR_A1B1 <- log(OR_A1B1)

  #list(
  #  logOR_A1B1 = logOR_A1B1,
  #  OR_A1B1    = OR_A1B1,
  #  p_A1       = p_A1,
  #  p_B1       = p_B1,
  #  p_casecase = p_casecase,
  #  b_unscaled = b_unscaled
  #)
  return(logOR_A1B1)
}

cc$BETA_DCMvHCM_OLSbackcomputed <- ccgwas_exact_to_casecase_logor(cc$OLS_beta, cc$OR_DCM, cc$FRQ_DCM, cc$OR_HCM, cc$FRQ_HCM)
cc$SE_DCMvHCM_OLSbackcomputed <- abs(cc$BETA_DCMvHCM_OLSbackcomputed) / (qnorm(cc$OLS_pval/2, lower.tail=F))
cc$OR_DCMvHCM_OLSbackcomputed <- exp(cc$BETA_DCMvHCM_OLSbackcomputed)
cc$BETA_DCMvHCM_exactbackcomputed <- ccgwas_exact_to_casecase_logor(cc$Exact_beta, cc$OR_DCM, cc$FRQ_DCM, cc$OR_HCM, cc$FRQ_HCM)
cc$SE_DCMvHCM_exactbackcomputed <- abs(cc$BETA_DCMvHCM_exactbackcomputed) / (qnorm(cc$OLS_pval/2, lower.tail=F))
cc$OR_DCMvHCM_exactbackcomputed <- exp(cc$BETA_DCMvHCM_exactbackcomputed)