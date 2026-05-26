INPUT_DIR=/${WD_PROJECT}/data/sumst_processed/UK_DCM_HCM
OUTPUT_DIR=/${WD_PROJECT}/data/sumst_processed/ldsc/MUNGE
MERGE_ALLELES=/gpfs/work5/0/gusr0607/dominicz/GWAS_SCA/HeritabilityCorrelation/ldsc/w_hm3.snplist
export WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM"
export WD_TOOLS="/home/dkramarenk/projects/tools/"

cd $WD_PROJECT
conda activate LAVA_2024 
conda activate r_env


# Preprocess in R

library(data.table)
library(stringr)

input_file <- "/home/dkramarenk/projects/LAVA/DCM_HCM/data/sumst_processed/UK_DCM_HCM/HCM_DCM_sumstats_SEP25.tsv.gz"
vcf_file <- "/home/dkramarenk/projects/tools/hg38/00-common_all.vcf.gz"
output_file <- file.path(Sys.getenv("OUTPUT_DIR"), "HCM_DCM_sumstats_SEP25.clean_rsid.tsv")

# read sumstats
dt <- fread(input_file)

# parse alleles column: ['C', 'T'] -> A1, A2
dt[, alleles_clean := gsub("\\[|\\]|'", "", alleles)]
dt[, c("A1", "A2") := tstrsplit(alleles_clean, ",\\s*")]

# keep relevant columns and drop missing values
dt_clean <- dt[
  !is.na(locus) & locus != "" &
  !is.na(A1) & A1 != "" &
  !is.na(A2) & A2 != "" &
  !is.na(beta) &
  !is.na(p_value),
  .(
    SNP = locus,
    A1 = A1,
    A2 = A2,
    BETA = beta,
    P = p_value
  )
]

# read VCF mapping: skip header lines beginning with ##
vcf <- fread(
  cmd = paste(
    "zcat", vcf_file,
    "| grep -v '^##'"
  ),
  sep = "\t",
  header = TRUE,
  select = c("#CHROM", "POS", "ID")
)

# rename columns
setnames(vcf, c("#CHROM", "POS", "ID"), c("CHR", "BP", "rsID"))

# remove missing IDs
vcf <- vcf[!is.na(rsID) & rsID != "." & rsID != "NA"]

# make chr:pos key to match dt_clean
vcf[, SNP := paste0("chr", CHR, ":", BP)]

# keep unique mapping
vcf_map <- unique(vcf[, .(SNP, rsID)])

# merge rsID into dt_clean
dt_clean <- merge(dt_clean, vcf_map, by = "SNP", all.x = TRUE)

# optional checks
cat("Total variants:", nrow(dt_clean), "\n")
cat("Mapped to rsID:", sum(!is.na(dt_clean$rsID)), "\n")
cat("Unmapped:", sum(is.na(dt_clean$rsID)), "\n")

# Total variants: 7916007
# Mapped to rsID: 7504689
# Unmapped: 411318

# keep only mapped rows and replace SNP with rsID
dt_clean_mapped <- dt_clean[
  !is.na(rsID)]

# write cleaned file
fwrite(dt_clean_mapped, "HCM_DCM_sumstats_SEP25.clean_rsid_chr_pos.tsv", sep = "\t", quote = FALSE, na = "NA")

library(data.table)
library(stringr)

# 1. Load targets
target_snps <- fread("/gpfs/work5/0/gusr0607/dkramarenko/LAVA/DCM_HCM/data/sumst_processed/ldsc/PART/top_snps/Top_SNP_CC.txt", header = FALSE, col.names = "rsID")

top_novel_snps <- fread("/gpfs/work5/0/gusr0607/dkramarenko/LAVA/DCM_HCM/data/sumst_processed/CC_power/top_snps_novel.txt", header = FALSE)

# 2. Filter mapped data
dt_targets <- dt_clean_mapped[rsID %in% target_snps$rsID]

# 3. Extract and format coordinates
dt_targets[, c("CHR", "BP") := tstrsplit(SNP, ":", fixed = TRUE)]
dt_targets[, BP := as.numeric(BP)]

# Ensure CHR has "chr" prefix (prevents "chrchr1" if already present)
dt_targets[, CHR := ifelse(startsWith(CHR, "chr"), CHR, paste0("chr", CHR))]

# 4. Calculate Windows
dt_targets[, `:=`(
  START = pmax(0, BP - 250000),
  END   = BP + 250000
)]

# 5. Create temp numeric CHR for proper sorting (e.g., chr2 before chr10)
dt_targets[, chr_num := as.numeric(gsub("chr", "", CHR, ignore.case = TRUE))]
setorder(dt_targets, chr_num, START)

# 6. Final BED selection and write
bed_output <- dt_targets[, .(CHR, START, END, rsID)]

# 5. Write to file
name <- "CC_repl_rsID"
fwrite(bed_output, 
       paste0("gws_loci_", tolower(name), "_merged.bed"), 
       sep = "\t", 
       col.names = FALSE, 
       quote = FALSE)

bed_output <- dt_targets[, .(CHR, START, END)]

# 5. Write to file
name <- "CC_repl"
fwrite(bed_output, 
       paste0("gws_loci_", tolower(name), "_merged.bed"), 
       sep = "\t", 
       col.names = FALSE, 
       quote = FALSE)



conda activate ldsc

export WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM"
export WD_TOOLS="/home/dkramarenk/projects/tools/"
INPUT_DIR=/${WD_PROJECT}/data/sumst_processed/UK_DCM_HCM
OUTPUT_DIR=/${WD_PROJECT}/data/sumst_processed/ldsc/MUNGE
MERGE_ALLELES=/gpfs/work5/0/gusr0607/dominicz/GWAS_SCA/HeritabilityCorrelation/ldsc/w_hm3.snplist



# Munge
${WD_TOOLS}/ldsc/munge_sumstats.py \
  --merge-alleles ${MERGE_ALLELES} \
  --sumstats HCM_DCM_sumstats_SEP25.clean_rsid.tsv \
  --snp SNP \
  --a1 A1 \
  --a2 A2 \
  --p P \
  --N 2683 \
  --signed-sumstats BETA,0 \
  --chunksize 500000 \
  --out ${OUTPUT_DIR}/HCM_DCM_sumstats_SEP25_repl

# Read 1217311 SNPs for allele merge.
# Reading sumstats from HCM_DCM_sumstats_SEP25.clean_rsid.tsv into memory 500000 SNPs at a time.
# ................ done
# Read 7504689 SNPs from --sumstats file.
# Removed 6349338 SNPs not in --merge-alleles.
# Removed 0 SNPs with missing values.
# Removed 0 SNPs with INFO <= 0.9.
# Removed 0 SNPs with MAF <= 0.01.
# Removed 0 SNPs with out-of-bounds p-values.
# Removed 0 variants that were not SNPs or were strand-ambiguous.
# 1155351 SNPs remain.
# Removed 86 SNPs with duplicated rs numbers (1155265 SNPs remain).
# Using N = 2683.0
# Median value of SIGNED_SUMSTATS was 4.15588059233e-05, which seems sensible.
# Removed 43 SNPs whose alleles did not match --merge-alleles (1155222 SNPs remain).
# Writing summary statistics for 1217311 SNPs (1155222 with nonmissing beta) to //home/dkramarenk/projects/LAVA/DCM_HCM/data/sumst_processed/ldsc/MUNGE/HCM_DCM_sumstats_SEP25_repl.sumstats.gz.

# CC loci

name="CC_repl"

awk 'BEGIN{FS=OFS="\t"}
NR==FNR{rs[$1]; next}
FNR==1{
  for(i=1;i<=NF;i++){
    if($i=="CHR") chr=i
    if($i=="BP") bp=i
    if($i=="rsID") rsid=i
  }
  next
}
($rsid in rs){
  s=$bp-250000; if(s<0)s=0
  e=$bp+250000
  print "chr"$chr, s, e, $rsid
}' top_snps/Top_SNP_CC.txt \
   /${WD_PROJECT}/data/sumst_processed/UK_DCM_HCM/HCM_DCM_sumstats_SEP25.ldsc.rs.tsv > gws_loci_${name,,}_merged.bed

########################################################
library(data.table)
library(stringr)

# Define file paths and name
name <- "CC_repl"
wd_project <- Sys.getenv("WD_PROJECT") # Ensure this env var is set

# 1. Load the target list of SNPs (assuming 1st col is rsID)
target_snps <- fread(paste0("/", wd_project, "/data/sumst_processed/ldsc/PART/top_snps/Top_SNP_CC.txt"), header = FALSE, col.names = "rsID")

# 2. Load the rsID mapping file
rs_map <- fread(paste0("/", wd_project, "/data/sumst_processed/UK_DCM_HCM/HCM_DCM_sumstats_SEP25.ldsc.rs.tsv"), 
                select = c("rsID"))

# 3. Load the file with chr:pos (assuming it matches rs_map by row order)
pos_map <- fread(paste0("/", wd_project, "/data/sumst_processed/UK_DCM_HCM/HCM_DCM_sumstats_SEP25.ldsc.tsv"), select = c("rsID"))
setnames(pos_map, "rsID", "chr_pos")

# 4. Combine and Filter
# Bind columns since they represent the same data rows
mapping <- cbind(rs_map, pos_map)
# Filter for only your target SNPs
final_data <- mapping[rsID %in% target_snps$rsID]

# 5. Extract Chromosome and Position from "chr1:822354"
# Split by ":" and remove "chr" prefix if needed
coords <- str_split_fixed(final_data$chr_pos, ":", 2)
final_data[, chr := coords[, 1]]
final_data[, bp := as.numeric(coords[, 2])]

# 6. Calculate Windows (±250kb)
final_data[, s := pmax(0, bp - 250000)]
final_data[, e := bp + 250000]

# 7. Select BED columns and write output
bed_output <- final_data[, .(chr, s, e, rsID)]

output_file <- paste0("gws_loci_", tolower(name), "_merged.bed")
fwrite(bed_output, output_file, sep = "\t", col.names = FALSE, quote = FALSE)
########################################################

cd ${WD_PROJECT}//data/sumst_processed/ldsc/PART

conda activate ldsc

for chr in {1..22}; do
Rscript ${WD_PROJECT}/scripts//make_ldsc_binary_annot.R \
  gws_loci_cc_repl_merged.bed \
  ${WD_TOOLS}/LDSC_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr}.bim \
  annot_shared_ld/gws_loci_${name,,}.${chr}.annot.gz "full-annot"
done

import pandas as pd

for chr in range(1, 23):
    path = f"annot_shared_ld/gws_loci_cc_repl.{chr}.l2.ldscore.gz"
    try:
        df = pd.read_csv(path, sep=r"\s+")
        print(f"OK: chr {chr}, shape={df.shape}")
    except Exception as e:
        print(f"FAIL: chr {chr}")
        print(e)
name="cc_repl"

for chr in {1..1}; do
python ${WD_TOOLS}/ldsc/ldsc.py \
  --l2 \
  --bfile ${WD_TOOLS}/LDSC_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} \
  --print-snps ${WD_TOOLS}/LDSC_files/hm3_no_MHC.list.txt \
  --ld-wind-cm 1 \
  --annot annot_shared_ld/gws_loci_${name,,}.${chr}.annot.gz \
  --out annot_shared_ld/gws_loci_${name,,}.${chr}
done
#!/bin/bash
#SBATCH --job-name=cc_loci_annot_process       # Job name
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --output=cc_loci_annot_process_%j.out

source /sw/arch/RHEL8/EB_production/2022/software/Anaconda3/2022.05/etc/profile.d/conda.sh

conda activate ldsc

for chr in {1..22}; do
Rscript ${WD_PROJECT}/scripts//make_ldsc_binary_annot.R \
  gws_loci_cc_repl_merged.bed \
  ${WD_TOOLS}/LDSC_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr}.bim \
  annot_shared_ld/gws_loci_${name,,}.${chr}.annot.gz "full-annot"
done


name="cc_repl"

for chr in {1..1}; do
python ${WD_TOOLS}/ldsc/ldsc.py \
  --l2 \
  --bfile ${WD_TOOLS}/LDSC_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} \
  --print-snps ${WD_TOOLS}/LDSC_files/hm3_no_MHC.list.txt \
  --ld-wind-cm 1 \
  --annot annot_shared_ld/gws_loci_${name,,}.${chr}.annot.gz \
  --out annot_shared_ld/gws_loci_${name,,}.${chr}
done

#!/bin/bash
#SBATCH --job-name=cc_loci_annot_process       # Job name
#SBATCH --time=03:00:00
#SBATCH --mem=32G
#SBATCH --output=cc_loci_part_hert_%j.out

source /sw/arch/RHEL8/EB_production/2022/software/Anaconda3/2022.05/etc/profile.d/conda.sh
conda activate LAVA_2024 
INPUT_DIR=/${WD_PROJECT}/data/sumst_processed/UK_DCM_HCM
OUTPUT_DIR=/${WD_PROJECT}/data/sumst_processed/ldsc/MUNGE
MERGE_ALLELES=/gpfs/work5/0/gusr0607/dominicz/GWAS_SCA/HeritabilityCorrelation/ldsc/w_hm3.snplist
export WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM"
export WD_TOOLS="/home/dkramarenk/projects/tools/"

for name in CC_GWAS CC_MTAG; do
# /gpfs/work5/0/gusr0607/dkramarenko/tools/bedtools2/bin/bedtools = bedtools
/gpfs/work5/0/gusr0607/dkramarenko/tools/bedtools2/bin/bedtools sort -i "gws_loci_${name,,}.bed" | /gpfs/work5/0/gusr0607/dkramarenko/tools/bedtools2/bin/bedtools merge -i - > "gws_loci_${name,,}_merged.bed"
done

cat gws_loci_cc_gwas.bed gws_loci_cc_mtag.bed \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - \
  > gws_loci_cc_combined_merged.bed


# 2. Generate .annot.gz  files using make_annot.py


for chr in {1..22}; do
Rscript ${WD_PROJECT}/scripts//make_ldsc_binary_annot.R \
  gws_loci_cc_combined_merged.bed \
  ${WD_TOOLS}/LDSC_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr}.bim \
  annot_shared_ld/gws_loci_cc_combined_merged.${chr}.annot.gz "full-annot"
done

for chr in {1..22}; do
python ${WD_TOOLS}/ldsc/ldsc.py \
  --l2 \
  --bfile ${WD_TOOLS}/LDSC_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} \
  --print-snps ${WD_TOOLS}/LDSC_files/hm3_no_MHC.list.txt \
  --ld-wind-cm 1 \
  --annot annot_shared_ld/gws_loci_cc_combined_merged.${chr}.annot.gz \
  --out annot_shared_ld/gws_loci_cc_combined_merged.${chr}
done

# 3. Making LD scores for annot file
for name in CC_GWAS CC_MTAG; do
for chr in {1..22}; do
python ${WD_TOOLS}/ldsc/ldsc.py \
  --l2 \
  --bfile ${WD_TOOLS}/LDSC_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} \
  --print-snps ${WD_TOOLS}/LDSC_files/hm3_no_MHC.list.txt \
  --ld-wind-cm 1 \
  --annot annot_shared_ld/gws_loci_${name,,}.${chr}.annot.gz \
  --out annot_shared_ld/gws_loci_${name,,}.${chr}
done
done 


conda activate ldsc
INPUT_DIR=/${WD_PROJECT}/data/sumst_processed/UK_DCM_HCM
OUTPUT_DIR=/${WD_PROJECT}/data/sumst_processed/ldsc/MUNGE
MERGE_ALLELES=/gpfs/work5/0/gusr0607/dominicz/GWAS_SCA/HeritabilityCorrelation/ldsc/w_hm3.snplist
export WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM"
export WD_TOOLS="/home/dkramarenk/projects/tools/"


python ${WD_TOOLS}/ldsc/ldsc.py \
--h2 /${WD_PROJECT}/data/sumst_processed/ldsc/MUNGE/DCM_GWAS_37_exclMYBPC3reg.sumstats.gz \
--ref-ld-chr /${WD_TOOLS}/LDSC_files/baselineLD.,/${WD_PROJECT}/data/sumst_processed/ldsc/PART/annot_shared_ld/gws_loci_cc_combined_merged. \
--frqfile-chr ${WD_TOOLS}/LDSC_files/1000G_Phase3_frq/1000G.EUR.QC. \
--w-ld-chr ${WD_TOOLS}/LDSC_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--overlap-annot --print-cov --print-coefficients --print-delete-vals \
--out part_out_cc_march_2026/march_2026_cc_loci_DCM_sign_baselineLD

python ${WD_TOOLS}/ldsc/ldsc.py \
--h2 /${WD_PROJECT}/data/sumst_processed/ldsc/MUNGE/HCM_GWAS_37_exclMYBPC3reg.sumstats.gz \
--ref-ld-chr /${WD_TOOLS}/LDSC_files/baselineLD.,/${WD_PROJECT}/data/sumst_processed/ldsc/PART/annot_shared_ld/gws_loci_cc_combined_merged. \
--frqfile-chr ${WD_TOOLS}/LDSC_files/1000G_Phase3_frq/1000G.EUR.QC. \
--w-ld-chr ${WD_TOOLS}/LDSC_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--overlap-annot --print-cov --print-coefficients --print-delete-vals \
--out part_out_cc_march_2026/march_2026_cc_loci_HCM_sign_baselineLD

--h2 /${WD_PROJECT}/data/sumst_processed/ldsc/MUNGE/CC_GWAS_37_exclMYBPC3reg.sumstats.gz \
--ref-ld-chr /${WD_TOOLS}/LDSC_files/baselineLD.,/${WD_PROJECT}/data/sumst_processed/ldsc/PART/annot_shared_ld/gws_loci_cc_combined_merged. \
--frqfile-chr ${WD_TOOLS}/LDSC_files/1000G_Phase3_frq/1000G.EUR.QC. \
--w-ld-chr ${WD_TOOLS}/LDSC_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--overlap-annot --print-cov --print-coefficients --print-delete-vals \
--out part_out_cc_march_2026/march_2026_cc_loci_CCGWAS_sign_baselineLD_overlap

python ${WD_TOOLS}/ldsc/ldsc.py \
--h2 /${WD_PROJECT}/data/sumst_processed/ldsc/MUNGE/HCM_DCM_sumstats_SEP25.ldsc.rs.sumstats.gz \
--ref-ld-chr /${WD_TOOLS}/LDSC_files/baselineLD.,/${WD_PROJECT}/data/sumst_processed/ldsc/PART/annot_shared_ld/gws_loci_cc_combined_merged. \
--frqfile-chr ${WD_TOOLS}/LDSC_files/1000G_Phase3_frq/1000G.EUR.QC. \
--w-ld-chr ${WD_TOOLS}/LDSC_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--overlap-annot --print-cov --print-coefficients --print-delete-vals \
--out part_out_cc_march_2026/march_2026_cc_loci_CCrepl_sept_sign_baselineLD_overlap

python ${WD_TOOLS}/ldsc/ldsc.py \
--h2 /${WD_PROJECT}/data/sumst_processed/ldsc/MUNGE/OCT25_results_NoOverlaps.tsv.gz.tsv.gz.sumstats.gz \
--ref-ld-chr /${WD_TOOLS}/LDSC_files/baselineLD.,/${WD_PROJECT}/data/sumst_processed/ldsc/PART/annot_shared_ld/gws_loci_cc_combined_merged. \
--frqfile-chr ${WD_TOOLS}/LDSC_files/1000G_Phase3_frq/1000G.EUR.QC. \
--w-ld-chr ${WD_TOOLS}/LDSC_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--overlap-annot --print-cov --print-coefficients --print-delete-vals \
--out part_out_cc_march_2026/march_2026_cc_loci_CCrepl_oct_sign_baselineLD_overlap

conda activate ldsc

python ${WD_TOOLS}/ldsc/ldsc.py \
--h2 /${WD_PROJECT}/data/sumst_processed/ldsc/MUNGE/OCT25_results_NoOverlaps.tsv.gz.tsv.gz.sumstats.gz \
--ref-ld-chr /${WD_TOOLS}/LDSC_files/baselineLD. \
--frqfile-chr ${WD_TOOLS}/LDSC_files/1000G_Phase3_frq/1000G.EUR.QC. \
--w-ld-chr ${WD_TOOLS}/LDSC_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--print-cov --print-coefficients --print-delete-vals \
--out part_out_cc_gwas_repl/cc_loci_sign_baselineLD



