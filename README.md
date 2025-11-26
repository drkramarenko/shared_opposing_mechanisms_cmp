# Leveraging shared and opposing genetic mechanisms in heritable cardiomyopathies
Internal code and data layout for the project **shared_opposing_mechanisms_cmp**.

This repository documents the analysis pipeline and scripts used in the manuscript
"Leveraging the shared and opposing genetic mechanisms in the heritable cardiomyopathies".

All figures are generated using R (v4.3.1) and the tidyverse / Bioconductor ecosystem.  
A reproducible description of the R environment is provided in the **Code Availability** section.

Detailed documentation for specific figure inputs can be found in: [`README_figures.md.sh`](code/example_munging.sh) 

## Directory layout (high level)

- `code/` – all analysis scripts (R, Bash, etc.)
- `data/` –  input data

> Note: Only `code/` and small helper files will be shared publicly with the paper.
---

## Step 1 – Obtain DCM and HCM GWAS summary statistics

1. Identify the two primary GWAS datasets:
   - DCM GWAS summary statistics from: https://www.nature.com/articles/s41588-024-01975-5
    Summary statistics for our GWAS meta-analyses have been made available for download through the Cardiovascular Disease Knowledge Portal (https://cvd.hugeamp.org/downloads.html)

   - HCM GWAS summary statistics from: https://www.nature.com/articles/s41588-025-02087-4
     Full GWAS summary statistics of HCM, HCMSARC−, HCMSARC+, MTAG and ten LV traits are available on the GWAS catalog (accession IDs GCST90435254 (https://www.ebi.ac.uk/gwas/studies/GCST90435254) –GCST90435267(https://www.ebi.ac.uk/gwas/studies/GCST90435267))


## Step 2 — Summary of QC for processed summary statistics (DCM/HCM)

QC filters:

1. Effect allele frequency filter:  
   - `EAFREQ >= 0.005 & EAFREQ <= 0.995`

2. Sample-size filters:  
   - For DCM GWAS: `N_cases >= 0.7 * max(N_cases)` 
   - For DCM MTAG: `N_Neff >= 0.7 * max(N_Neff)`
   - For HCM GWAS: `N_cases >= 0.96 * max(N_cases)`  
   - For HCM MTAG: `N_Neff >= 0.96 * max(N_Neff)` 

3. Exclusion of MHC-like extended region on chromosome 11:  
   - DCM: exclude variants with `CHR == 11 & POS ∈ [29,978,453; 80,288,956]`  
   - HCM: exclude variants with `CHR == 11 & POS ∈ [30,000,000; 80,000,000]`

Output: 
- harmonized_dcm.tsv.gz
- harmonized_hcm.tsv.gz
<!-- #### ===== packages and versions 
# conda activate LAVA_2024
# conda list
# conda env export > environment_LAVA_2024.yml
# R --version -->
<!-- Rscript -e "packageVersion('optparse')"
Rscript -e "packageVersion('data.table')"
Rscript -e "packageVersion('dplyr')"
Rscript -e "packageVersion('readr')"
Rscript -e "packageVersion('tidyr')" -->

---
## Step 3 - Genetic correlations

### 3.1 Global genetic correlation (rg)

Estimation of univariate SNP-heritability and bivariate genetic correlation using LD Score Regression (LDSC).

#### 3.1.1 LDSC software
We use LD Score Regression from Bulik-Sullivan et al.
Repository: https://github.com/bulik/ldsc

##### 3.1.2 Munging 

Example script: [`code/example_munging.sh`](code/example_munging.sh)

##### 3.1.3 LDSC Genetic Correlation 
Example script: [`code/example_ldsc_rg.sh`](code/example_ldsc_rg.sh)

##### TO ADD === version ==== source activate ldsc  # activate your conda env
Version recorded automatically in:
##### code/environment/ldsc_version.txt### 4.3 CC-GWAS: DCM vs HCM

##### 3.1.4 Figure 2a: Genetic correlation heatmap (DCM–HCM + MRI traits):
Example script:  [`code/figure2a.r`](code/figure2a.R)

##### 3.1.4 Figure 2b: SNP effect concordance between DCM and HCM:
Example script: [`code/figure2b.r`](code/figure2a.R)

### 3.2 Local genetic correlation (LAVA)

LAVA was executed in a dedicated conda environment (`LAVA_2024`) to ensure consistent versions of PLINK, PLINK2, bcftools, and supporting R packages. Details and additional input files and script examples here: https://github.com/josefin-werme/LAVA 

Version used for our paper: https://github.com/josefin-werme/LAVA/releases/tag/v0.1.0 

#### 3.2.1 Activate environment
```bash
conda activate LAVA_2024
# R version (used for LAVA)
conda install -c bioconda r-base=4.4
```

```R
# Install LAVA from GitHub
install.packages("remotes")
remotes::install_github("josefin-werme/LAVA")
```

#### 3.2.2 Example script: [`code/run_lava_local_rg.r`](code/run_lava_local_rg.R)

#### 3.2.3 Figures 2d and 2f (LAVA locus annotation and Manhattan-type plots)


#### 3.3.2 Input Configuration File for LAVA

LAVA requires a simple tab-delimited configuration file specifying, for each phenotype: [`data/LAVA/input.info.txt`](data/input.info.txt)

