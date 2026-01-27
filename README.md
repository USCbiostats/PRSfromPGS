# PRSfromPGS

Calculate Polygenic Risk Scores (PRS) from genetic data using [PGS Catalog](https://www.pgscatalog.org/)-format scoring models and bgzipped VCF files.

## Installation

PRSfromPGS depends on Bioconductor packages. Install them first, then install PRSfromPGS from GitHub:

```r
# Install Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("VariantAnnotation", "Rsamtools"))

# Install PRSfromPGS
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("USCbiostats/PRSfromPGS")
```

For parallel processing support, also install:

```r
# Parallel chromosome processing
install.packages(c("future", "future.apply", "parallelly"))

# Parallel VCF indexing
BiocManager::install("BiocParallel")
```

## Quick Start

```r
library(PRSfromPGS)

# Step 1: Prepare VCF files (creates per-chromosome RDS summaries)
vcf_files <- list.files("data/", pattern = "\\.vcf\\.gz$", full.names = TRUE)
prepareVCFfiles(vcf_files, genome = "hg38", tabix = TRUE)

# Step 2: Read a PGS Catalog scoring file
pgs_model <- readPGSmodel("PGS000001.txt")

# Step 3: Calculate PRS across all autosomes
result <- calculatePRS(pgs_model)
prs_values <- result$prs  # Named numeric vector (one score per sample)
```

## Usage

### Preparing VCF Data

VCF files must be bgzipped (`.vcf.gz`). Each file should contain variants for a single chromosome.

```r
# Index and summarize VCF files
vcf_files <- c("chr1.vcf.gz", "chr2.vcf.gz", "chr3.vcf.gz")
prepareVCFfiles(vcf_files, genome = "hg38", tabix = TRUE)

# This creates chromosome1VCFsummary.rds, chromosome2VCFsummary.rds, etc.
# in the current working directory
```

### Reading PGS Models

Reads tab-delimited scoring files from the [PGS Catalog](https://www.pgscatalog.org/). Handles harmonized columns (`hm_chr`, `hm_pos`, `hm_inferOtherAllele`) automatically.

```r
pgs_model <- readPGSmodel("PGS000001.txt")
# Returns a data frame with columns:
#   rsID, chr_name, chr_position, effect_allele, other_allele, effect_weight
```

### Calculating PRS

#### Single model

```r
result <- calculatePRS(pgs_model)
result$prs                # Named numeric vector of PRS values
result$unmatched_by_chr   # Count of unmatched variants per chromosome
result$unmatched_rsIDs    # rsIDs of unmatched variants per chromosome
result$excluded_chr_counts # Variants on non-processed chromosomes
```

#### Multiple models (efficient)

`fitPRSmodels()` reads each chromosome's VCF data only once across all models:

```r
model_files <- c("PGS000001.txt", "PGS000002.txt", "PGS000003.txt")
results <- fitPRSmodels(model_files)
results$prs  # Matrix: rows = samples, columns = models
```

#### Parallel processing

```r
# Single model with parallel chromosome processing
result <- calculatePRS(pgs_model, parallel = TRUE, n_cores = 4)

# Multiple models in parallel
results <- fitPRSmodels(model_files, parallel = TRUE, n_cores = 4)
```

### Lower-Level Functions

For more control, you can use the underlying functions directly:

```r
# Summarize a single VCF file
vcf_summary <- summarizevcf("chr21.vcf.gz", genome = "hg38")

# Read dosage values for specific SNP indices
dosages <- getvcfsnp(vcf_summary, snp_indices = c(1, 5, 10))

# Match PGS model variants to VCF variants
matches <- findSNPs(pgs_model, vcf_summary$snp_summary$snp_df)

# Calculate PRS separately for REF and ALT matches
prs_alt <- calculatePRSalt(matches$alt_matches, pgs_model, vcf_summary)
prs_ref <- calculatePRSref(matches$ref_matches, pgs_model, vcf_summary)
```

## How It Works

PRSfromPGS matches PGS model variants to VCF variants by genomic position and allele, then computes weighted dosage sums:

- **ALT matches**: When the effect allele matches the ALT allele, the VCF dosage (DS field) is used directly
- **REF matches**: When the effect allele matches the REF allele, the effect allele dosage is `2 - DS`

The final PRS for each sample is the sum of `effect_weight * effect_allele_dosage` across all matched variants and all chromosomes.

## Dependencies

| Package | Source | Required |
|---------|--------|----------|
| VariantAnnotation | Bioconductor | Yes |
| Rsamtools | Bioconductor | Yes |
| data.table | CRAN | Yes |
| BiocParallel | Bioconductor | No (parallel VCF indexing) |
| future | CRAN | No (parallel chromosome processing) |
| future.apply | CRAN | No (parallel chromosome processing) |
| parallelly | CRAN | No (parallel chromosome processing) |

## License

MIT
