#!/usr/bin/env Rscript

## make_ldsc_binary_annot.R
##
## Create a binary LDSC annotation file from:
##  - a BED file with regions of interest
##  - a PLINK .bim file (SNP positions)
##
## Usage:
##   Rscript make_ldsc_binary_annot.R \
##     <regions.bed> \
##     <reference.bim> \
##     <out.annot.gz> \
##     <mode>
##
## Arguments:
##   1) regions.bed   : BED file with at least 3 columns: chr, start, end
##                      (chr may be "chr1" or "1"; start is 0-based)
##   2) reference.bim : PLINK BIM file (e.g. 1000G.EUR.QC.1.bim)
##   3) out.annot.gz  : Output path for LDSC .annot.gz file
##   4) mode          : "full-annot" (default) → LDSC format: CHR, BP, SNP, CM, ANNOT
##                      "thin-annot"          → only ANNOT column (0/1)
##
## Notes:
##   - The script uses GenomicRanges to intersect SNP positions with BED regions.
##   - BED start positions are converted from 0-based to 1-based for R/GRanges.

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3 || length(args) > 4) {
  stop(
    "Usage:\n",
    "  Rscript make_ldsc_binary_annot.R <regions.bed> <reference.bim> <out.annot.gz> [mode]\n",
    "where mode is either 'full-annot' (default) or 'thin-annot'.\n",
    call. = FALSE
  )
}

bed_file   <- args[1]
bim_file   <- args[2]
annot_file <- args[3]
mode       <- ifelse(length(args) == 4, args[4], "full-annot")

suppressPackageStartupMessages({
  library(GenomicRanges)
})

## ----------------------------
## 1. Read and process BED file
## ----------------------------

annot_bed <- read.table(bed_file, header = FALSE, stringsAsFactors = FALSE)[, 1:3]
colnames(annot_bed) <- c("chr", "start", "end")

## Remove "chr" prefix if present
annot_bed$chr <- gsub("^chr", "", annot_bed$chr)

## Convert 0-based start to 1-based (BED → genomic coordinates)
annot_bed$start <- annot_bed$start + 1

## ----------------------------
## 2. Prepare output directory
## ----------------------------

dir_annot_ldsc <- dirname(annot_file)
dir.create(dir_annot_ldsc, showWarnings = FALSE, recursive = TRUE)

## ----------------------------
## 3. Read BIM and build GRanges
## ----------------------------

bim_data <- read.table(bim_file, header = FALSE, stringsAsFactors = FALSE)
colnames(bim_data) <- c("CHR", "SNP", "CM", "BP", "A1", "A2")

## LDSC expects: CHR, BP, SNP, CM, ANNOT (or just ANNOT for thin mode)
annot_ldsc <- bim_data[, c("CHR", "BP", "SNP", "CM")]

annot_ldsc_gr <- makeGRangesFromDataFrame(
  annot_ldsc,
  seqnames.field   = "CHR",
  start.field      = "BP",
  end.field        = "BP",
  keep.extra.columns = TRUE
)

annot_bed_gr <- makeGRangesFromDataFrame(
  annot_bed,
  seqnames.field   = "chr",
  start.field      = "start",
  end.field        = "end",
  keep.extra.columns = TRUE
)

## ----------------------------
## 4. Find overlaps and set ANNOT
## ----------------------------

overlaps <- findOverlaps(
  query   = annot_ldsc_gr,
  subject = annot_bed_gr
)

annot_ldsc$ANNOT <- 0L
if (length(overlaps) > 0) {
  idx_overlaps <- as.data.frame(overlaps)
  annot_ldsc[idx_overlaps$queryHits, "ANNOT"] <- 1L
}

## ----------------------------
## 5. Write output
## ----------------------------

if (mode == "thin-annot") {
  ## Only ANNOT column (0/1), as some workflows prefer this minimal format
  out_dat <- data.frame(ANNOT = annot_ldsc[["ANNOT"]])
} else {
  ## Full LDSC format: CHR, BP, SNP, CM, ANNOT
  out_dat <- annot_ldsc
}

con <- gzfile(annot_file, open = "wt")
write.table(
  out_dat,
  con,
  col.names = TRUE,
  row.names = FALSE,
  quote     = FALSE,
  sep       = "\t"
)
close(con)

cat("LDSC annotation file saved to:", annot_file, "\n")
cat("Mode:", mode, "\n")
