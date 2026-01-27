# Helpers to build test data frames matching expected schemas

make_prs <- function(positions, effect_alleles, rsIDs = NULL, weights = rep(1, length(positions))) {
  if (is.null(rsIDs)) rsIDs <- paste0("rs", seq_along(positions))
  data.frame(
    rsID = rsIDs,
    chr_name = rep("1", length(positions)),
    chr_position = positions,
    effect_allele = effect_alleles,
    other_allele = rep("N", length(positions)),
    effect_weight = weights,
    stringsAsFactors = FALSE
  )
}

make_vcf <- function(positions, refs, alts, ids = NULL) {
  if (is.null(ids)) ids <- paste0("snp", seq_along(positions))
  data.frame(
    SNP_ID = ids,
    seqnames = rep("1", length(positions)),
    position = positions,
    REF = refs,
    ALT = alts,
    stringsAsFactors = FALSE
  )
}

test_that("ALT-only matches: effect_allele matches ALT at same position", {
  prs <- make_prs(positions = c(100, 200), effect_alleles = c("A", "T"))
  vcf <- make_vcf(positions = c(100, 200), refs = c("G", "C"), alts = c("A", "T"))

  result <- findSNPs(prs, vcf)

  expect_length(result$alt_matches, 2)
  expect_length(result$ref_matches, 0)
  expect_length(result$no_matches, 0)

  # Each alt_match element is c(prs_idx, vcf_idx, ...)
  alt_prs_indices <- vapply(result$alt_matches, `[`, integer(1), 1)
  expect_setequal(alt_prs_indices, c(1L, 2L))
})

test_that("REF-only matches: effect_allele matches REF at same position", {
  prs <- make_prs(positions = c(100, 200), effect_alleles = c("G", "C"))
  vcf <- make_vcf(positions = c(100, 200), refs = c("G", "C"), alts = c("A", "T"))

  result <- findSNPs(prs, vcf)

  expect_length(result$ref_matches, 2)
  expect_length(result$alt_matches, 0)
  expect_length(result$no_matches, 0)

  ref_prs_indices <- vapply(result$ref_matches, `[`, integer(1), 1)
  expect_setequal(ref_prs_indices, c(1L, 2L))
})

test_that("mixed REF and ALT matches in one call", {
  prs <- make_prs(positions = c(100, 200), effect_alleles = c("G", "T"))
  vcf <- make_vcf(positions = c(100, 200), refs = c("G", "C"), alts = c("A", "T"))

  result <- findSNPs(prs, vcf)

  expect_length(result$ref_matches, 1)
  expect_length(result$alt_matches, 1)
  expect_length(result$no_matches, 0)

  # PRS row 1 (effect_allele=G) matches VCF row 1 REF=G
  expect_equal(result$ref_matches[[1]][1], 1L)
  # PRS row 2 (effect_allele=T) matches VCF row 2 ALT=T
  expect_equal(result$alt_matches[[1]][1], 2L)
})

test_that("no matches when position or allele mismatch", {
  prs <- make_prs(positions = c(100, 200), effect_alleles = c("A", "T"))
  vcf <- make_vcf(positions = c(300, 400), refs = c("G", "C"), alts = c("A", "T"))

  result <- findSNPs(prs, vcf)

  expect_length(result$ref_matches, 0)
  expect_length(result$alt_matches, 0)
  expect_equal(result$no_matches, c(1L, 2L))
})

test_that("all variants unmatched returns full index vector in no_matches", {
  prs <- make_prs(positions = c(100, 200, 300), effect_alleles = c("X", "Y", "Z"))
  vcf <- make_vcf(positions = c(100, 200, 300), refs = c("A", "B", "C"), alts = c("D", "E", "F"))

  result <- findSNPs(prs, vcf)

  expect_equal(sort(result$no_matches), 1:3)
  expect_length(result$ref_matches, 0)
  expect_length(result$alt_matches, 0)
})

test_that("multiple VCF rows matching one PRS variant (cartesian)", {
  prs <- make_prs(positions = 100, effect_alleles = "A")
  vcf <- make_vcf(
    positions = c(100, 100, 100),
    refs = c("G", "G", "G"),
    alts = c("A", "A", "T"),
    ids = c("snp1", "snp2", "snp3")
  )

  result <- findSNPs(prs, vcf)

  # PRS row 1 should match VCF rows 1 and 2 on ALT
  expect_length(result$alt_matches, 1)
  match_vec <- result$alt_matches[[1]]
  expect_equal(match_vec[1], 1L)  # prs_idx
  expect_true(length(match_vec) >= 3)  # prs_idx + at least 2 vcf indices
  expect_setequal(match_vec[-1], c(1L, 2L))
  expect_length(result$no_matches, 0)
})

test_that("empty PRS input returns empty match lists", {
  prs <- data.frame(
    rsID = character(0), chr_name = character(0), chr_position = integer(0),
    effect_allele = character(0), other_allele = character(0),
    effect_weight = numeric(0), stringsAsFactors = FALSE
  )
  vcf <- make_vcf(positions = c(100, 200), refs = c("A", "G"), alts = c("T", "C"))

  result <- findSNPs(prs, vcf)

  expect_length(result$ref_matches, 0)
  expect_length(result$alt_matches, 0)
  expect_length(result$no_matches, 0)
})

test_that("output structure: list with ref_matches, alt_matches, no_matches", {
  prs <- make_prs(positions = 100, effect_alleles = "A")
  vcf <- make_vcf(positions = 100, refs = c("G"), alts = c("A"))

  result <- findSNPs(prs, vcf)

  expect_type(result, "list")
  expect_named(result, c("ref_matches", "alt_matches", "no_matches"))
  expect_type(result$ref_matches, "list")
  expect_type(result$alt_matches, "list")
  expect_type(result$no_matches, "integer")
})
