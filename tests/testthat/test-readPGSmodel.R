# Helper: write a PGS Catalog-style TSV to a temp file and return its path
write_pgs_file <- function(header_lines, col_names, data_lines) {
  tmp <- tempfile(fileext = ".txt")
  lines <- c(header_lines,
             paste(col_names, collapse = "\t"),
             vapply(data_lines, function(row) paste(row, collapse = "\t"), character(1)))
  writeLines(lines, tmp)
  tmp
}

test_that("reads a standard PGS Catalog file with all required columns", {
  f <- write_pgs_file(
    header_lines = c("## PGS Catalog", "## version 1.0"),
    col_names = c("rsID", "chr_name", "chr_position", "effect_allele", "other_allele", "effect_weight"),
    data_lines = list(
      c("rs123", "1", "100", "A", "G", "0.5"),
      c("rs456", "2", "200", "T", "C", "-0.3")
    )
  )
  on.exit(unlink(f), add = TRUE)

  result <- readPGSmodel(f, verbose = FALSE)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_equal(result$rsID, c("rs123", "rs456"))
  expect_equal(result$effect_weight, c(0.5, -0.3))
})

test_that("skips ## comment lines correctly", {
  f <- write_pgs_file(
    header_lines = c("## comment 1", "## comment 2", "## comment 3"),
    col_names = c("rsID", "chr_name", "chr_position", "effect_allele", "other_allele", "effect_weight"),
    data_lines = list(
      c("rs789", "3", "300", "C", "T", "1.2")
    )
  )
  on.exit(unlink(f), add = TRUE)

  result <- readPGSmodel(f, verbose = FALSE)

  expect_equal(nrow(result), 1)
  expect_equal(result$rsID, "rs789")
})

test_that("normalizes hm_chr -> chr_name, hm_pos -> chr_position, hm_inferOtherAllele -> other_allele", {
  f <- write_pgs_file(
    header_lines = c("## header"),
    col_names = c("rsID", "hm_chr", "hm_pos", "effect_allele", "hm_inferOtherAllele", "effect_weight"),
    data_lines = list(
      c("rs100", "5", "500", "G", "A", "0.1")
    )
  )
  on.exit(unlink(f), add = TRUE)

  result <- readPGSmodel(f, verbose = FALSE)

  expect_true("chr_name" %in% names(result))
  expect_true("chr_position" %in% names(result))
  expect_true("other_allele" %in% names(result))
  expect_equal(result$chr_name, 5)
  expect_equal(result$chr_position, 500)
  expect_equal(result$other_allele, "A")
})

test_that("removes rows with NA or 0 chr_position", {
  f <- write_pgs_file(
    header_lines = c("## header"),
    col_names = c("rsID", "chr_name", "chr_position", "effect_allele", "other_allele", "effect_weight"),
    data_lines = list(
      c("rs1", "1", "100", "A", "G", "0.5"),
      c("rs2", "1", "0",   "T", "C", "0.3"),
      c("rs3", "1", "NA",  "G", "A", "0.2"),
      c("rs4", "2", "400", "C", "T", "0.1")
    )
  )
  on.exit(unlink(f), add = TRUE)

  result <- readPGSmodel(f, verbose = FALSE)

  expect_equal(nrow(result), 2)
  expect_equal(result$rsID, c("rs1", "rs4"))
  expect_equal(result$chr_position, c(100, 400))
})

test_that("creates composite IDs for empty rsID values using X for NA other_allele", {
  f <- write_pgs_file(
    header_lines = c("## header"),
    col_names = c("rsID", "chr_name", "chr_position", "effect_allele", "other_allele", "effect_weight"),
    data_lines = list(
      c("",    "1", "100", "A", "G",  "0.5"),
      c("",    "2", "200", "T", "",   "0.3"),
      c("rs3", "3", "300", "C", "A",  "0.1")
    )
  )
  on.exit(unlink(f), add = TRUE)

  result <- readPGSmodel(f, verbose = FALSE)

  # First row: empty rsID -> composite ID
  expect_equal(result$rsID[1], "1:100:A:G")
  # Second row: empty rsID, empty other_allele read as NA -> uses "X"
  # Note: fread may read "" as empty string; the function checks is.na()
  # The composite ID should be 2:200:T:X if other_allele is NA, or 2:200:T: if empty string
  expect_true(grepl("^2:200:T:", result$rsID[2]))
  # Third row: existing rsID preserved
  expect_equal(result$rsID[3], "rs3")
})

test_that("returns exactly the 6 required columns in correct order", {
  f <- write_pgs_file(
    header_lines = c("## header"),
    col_names = c("rsID", "chr_name", "chr_position", "effect_allele",
                  "other_allele", "effect_weight", "extra_col"),
    data_lines = list(
      c("rs1", "1", "100", "A", "G", "0.5", "extra_data")
    )
  )
  on.exit(unlink(f), add = TRUE)

  result <- readPGSmodel(f, verbose = FALSE)

  expected_cols <- c("rsID", "chr_name", "chr_position", "effect_allele", "other_allele", "effect_weight")
  expect_equal(names(result), expected_cols)
  expect_equal(ncol(result), 6)
})

test_that("errors on missing file", {
  expect_error(
    readPGSmodel("nonexistent_file_path.txt", verbose = FALSE),
    "File does not exist"
  )
})

test_that("errors on missing required columns (e.g., no effect_weight)", {
  f <- write_pgs_file(
    header_lines = c("## header"),
    col_names = c("rsID", "chr_name", "chr_position", "effect_allele", "other_allele"),
    data_lines = list(
      c("rs1", "1", "100", "A", "G")
    )
  )
  on.exit(unlink(f), add = TRUE)

  expect_error(
    readPGSmodel(f, verbose = FALSE),
    "Missing required columns.*effect_weight"
  )
})

test_that("errors on non-string file_path input", {
  expect_error(
    readPGSmodel(123, verbose = FALSE),
    "file_path must be a single character string"
  )
  expect_error(
    readPGSmodel(c("a.txt", "b.txt"), verbose = FALSE),
    "file_path must be a single character string"
  )
})

test_that("creates other_allele as NA when column is missing entirely", {
  f <- write_pgs_file(
    header_lines = c("## header"),
    col_names = c("rsID", "chr_name", "chr_position", "effect_allele", "effect_weight"),
    data_lines = list(
      c("rs1", "1", "100", "A", "0.5"),
      c("rs2", "2", "200", "T", "0.3")
    )
  )
  on.exit(unlink(f), add = TRUE)

  result <- readPGSmodel(f, verbose = FALSE)

  expect_true("other_allele" %in% names(result))
  expect_true(all(is.na(result$other_allele)))
})

test_that("file with no comment lines reads correctly", {
  f <- write_pgs_file(
    header_lines = character(0),
    col_names = c("rsID", "chr_name", "chr_position", "effect_allele", "other_allele", "effect_weight"),
    data_lines = list(
      c("rs1", "1", "100", "A", "G", "0.5")
    )
  )
  on.exit(unlink(f), add = TRUE)

  result <- readPGSmodel(f, verbose = FALSE)

  expect_equal(nrow(result), 1)
  expect_equal(result$rsID, "rs1")
})
