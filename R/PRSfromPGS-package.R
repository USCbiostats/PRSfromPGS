#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import VariantAnnotation
#' @import Rsamtools
#' @importFrom data.table as.data.table fread .I
#' @importFrom GenomicRanges seqnames
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges ranges
#' @importFrom Biostrings DNAStringSet
#' @importFrom methods is
#' @importFrom stats sd start
## usethis namespace: end
NULL

# Suppress R CMD check notes for data.table NSE symbols
utils::globalVariables(c(
  "prs_idx",
  "vcf_idx",
  ".I",
  ".",
  ":=",
  "position",
  "REF",
  "ALT",
  "i.prs_idx",
  "x.vcf_idx",
  "chr_position",
  "effect_allele"
))
