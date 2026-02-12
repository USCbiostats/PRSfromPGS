test_that("calculatePRS produces expected results for PGS000001 test model", {
  # Skip if test files not available
  mdlfile <- system.file("extdata", "PGS000001_test.txt", package = "PRSfromPGS")
  skip_if(mdlfile == "", "Test model file PGS000001_test.txt not available")

  res2file <- system.file("extdata", "PGS000001_test.rds", package = "PRSfromPGS")
  skip_if(res2file == "", "Reference result file PGS000001_test.rds not available")

  # Check for VCF files - need at least chr3
  vcffile3 <- system.file("extdata", "chr3.vcf.gz", package = "PRSfromPGS")
  skip_if(vcffile3 == "", "Test VCF file chr3.vcf.gz not available")

  # Create temp directory for VCF summaries
  temp_dir_path <- tempdir()

  # Find all available chromosome VCF files
  datapath <- system.file("extdata", package = "PRSfromPGS")
  all_chrs <- 1:22
  vcffiles <- paste0(datapath, "/chr", all_chrs, ".vcf.gz")
  available_vcfs <- vcffiles[file.exists(vcffiles)]

  skip_if(length(available_vcfs) == 0, "No VCF files available for testing")

  # Prepare VCF files
  prepareVCFfiles(available_vcfs, output_dir = temp_dir_path, verbose = FALSE)

  # Determine which chromosomes are available
  available_chrs <- all_chrs[file.exists(vcffiles)]

  # Read PGS model and calculate PRS
  mdl <- readPGSmodel(mdlfile, verbose = FALSE)
  res <- calculatePRS(mdl, available_chrs, temp_dir_path, verbose = FALSE)

  # Read expected result
  res2 <- readRDS(res2file)

  # Verify structure - calculatePRS returns a list
  expect_type(res, "list")
  expect_true("prs" %in% names(res), "Result should contain 'prs' element")
  expect_type(res$prs, "double")
  expect_true(!is.null(names(res$prs)), "PRS result should be a named numeric vector")

  # res2 is also a list with $prs element
  expect_type(res2, "list")
  expect_type(res2$prs, "double")

  # Main assertion: res$prs and res2$prs should be identical
  expect_equal(res$prs, res2$prs,
               info = "Calculated PRS should match reference result")

  # Additional checks
  expect_true(length(res$prs) > 0, "PRS result should contain sample scores")
  expect_true(all(!is.na(res$prs)), "PRS scores should not contain NA values")

  # Clean up VCF summary files
  for (chr in available_chrs) {
    summary_file <- file.path(temp_dir_path, paste0("chromosome", chr, "VCFsummary.rds"))
    if (file.exists(summary_file)) {
      unlink(summary_file)
    }
  }
})

test_that("calculatePRS handles single chromosome correctly", {
  # Skip if test files not available
  mdlfile <- system.file("extdata", "PGS000001_test.txt", package = "PRSfromPGS")
  skip_if(mdlfile == "", "Test model file PGS000001_test.txt not available")

  vcffile3 <- system.file("extdata", "chr3.vcf.gz", package = "PRSfromPGS")
  skip_if(vcffile3 == "", "Test VCF file chr3.vcf.gz not available")

  # Create temp directory
  temp_dir_path <- tempdir()

  # Prepare single VCF file
  prepareVCFfiles(vcffile3, output_dir = temp_dir_path, verbose = FALSE)

  # Calculate PRS for only chromosome 3
  mdl <- readPGSmodel(mdlfile, verbose = FALSE)
  res <- calculatePRS(mdl, 3, temp_dir_path, verbose = FALSE)

  # Verify result
  expect_type(res, "list")
  expect_true("prs" %in% names(res))
  if (!is.null(res$prs)) {
    expect_type(res$prs, "double")
    expect_true(!is.null(names(res$prs)))
    expect_true(length(res$prs) > 0)
  }

  # Clean up
  summary_file <- file.path(temp_dir_path, "chromosome3VCFsummary.rds")
  unlink(summary_file)
})

test_that("calculatePRS returns zero scores when no variants match", {
  # Create a PGS model with variants that won't match any real data
  tmp_model <- tempfile(fileext = ".txt")
  model_content <- paste(
    "## Test PGS model with non-matching variants",
    "rsID\tchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight",
    "rs_fake1\t3\t99999999\tA\tG\t1.0",
    "rs_fake2\t3\t99999998\tT\tC\t0.5",
    sep = "\n"
  )
  writeLines(model_content, tmp_model)
  on.exit(unlink(tmp_model), add = TRUE)

  # Skip if VCF not available
  vcffile3 <- system.file("extdata", "chr3.vcf.gz", package = "PRSfromPGS")
  skip_if(vcffile3 == "", "Test VCF file chr3.vcf.gz not available")

  temp_dir_path <- tempdir()
  prepareVCFfiles(vcffile3, output_dir = temp_dir_path, verbose = FALSE)

  mdl <- readPGSmodel(tmp_model, verbose = FALSE)
  res <- suppressWarnings(calculatePRS(mdl, 3, temp_dir_path, verbose = FALSE))

  # Should return NULL prs when no variants match
  expect_type(res, "list")
  expect_true("prs" %in% names(res))
  expect_null(res$prs, "PRS should be NULL when no variants match")

  # Clean up
  summary_file <- file.path(temp_dir_path, "chromosome3VCFsummary.rds")
  unlink(summary_file)
})

test_that("calculatePRS works with default directory", {
  # Skip if test files not available
  mdlfile <- system.file("extdata", "PGS000001_test.txt", package = "PRSfromPGS")
  skip_if(mdlfile == "", "Test model file PGS000001_test.txt not available")

  vcffile3 <- system.file("extdata", "chr3.vcf.gz", package = "PRSfromPGS")
  skip_if(vcffile3 == "", "Test VCF file chr3.vcf.gz not available")

  # Save current directory
  old_wd <- getwd()

  # Create and use a temp directory
  temp_dir <- tempfile()
  dir.create(temp_dir)
  setwd(temp_dir)

  on.exit({
    setwd(old_wd)
    unlink(temp_dir, recursive = TRUE)
  }, add = TRUE)

  # Prepare VCF file in current directory (default)
  prepareVCFfiles(vcffile3, verbose = FALSE)

  # Calculate PRS using default directory (current working directory)
  mdl <- readPGSmodel(mdlfile, verbose = FALSE)
  res <- suppressWarnings(calculatePRS(mdl, 3, verbose = FALSE))

  # Verify result
  expect_type(res, "list")
  expect_true("prs" %in% names(res))
  if (!is.null(res$prs)) {
    expect_type(res$prs, "double")
    expect_true(length(res$prs) > 0)
  }
})

test_that("fitPRSmodels processes multiple models efficiently", {
  # Skip if test files not available
  vcffile3 <- system.file("extdata", "chr3.vcf.gz", package = "PRSfromPGS")
  skip_if(vcffile3 == "", "Test VCF file chr3.vcf.gz not available")

  # Find available test model files
  datapath <- system.file("extdata", package = "PRSfromPGS")
  model_files <- c(
    system.file("extdata", "PGS000001_test.txt", package = "PRSfromPGS"),
    system.file("extdata", "PGS000002_test.txt", package = "PRSfromPGS"),
    system.file("extdata", "PGS000003_test.txt", package = "PRSfromPGS")
  )

  # Filter to only existing files
  model_files <- model_files[file.exists(model_files) & model_files != ""]
  skip_if(length(model_files) < 2, "Need at least 2 test model files for fitPRSmodels test")

  # Create temp directory
  temp_dir_path <- tempdir()

  # Prepare VCF file
  prepareVCFfiles(vcffile3, output_dir = temp_dir_path, verbose = FALSE)

  # Fit multiple models using file paths
  res <- suppressWarnings(fitPRSmodels(model_files, 3, temp_dir_path, verbose = FALSE))

  # Verify result structure - fitPRSmodels returns a list with prs matrix
  expect_type(res, "list")
  expect_true("prs" %in% names(res), "Result should contain 'prs' element")
  expect_true("unmatched_by_chr" %in% names(res))
  expect_true("unmatched_rsIDs" %in% names(res))
  expect_true("excluded_chr_counts" %in% names(res))

  # PRS should be a matrix with samples as rows and models as columns
  if (!is.null(res$prs)) {
    expect_true(is.matrix(res$prs), "PRS should be a matrix")
    expect_equal(ncol(res$prs), length(model_files),
                 info = "PRS matrix should have one column per model")
    expect_true(nrow(res$prs) > 0, "PRS matrix should have sample data")
  }

  # Clean up
  summary_file <- file.path(temp_dir_path, "chromosome3VCFsummary.rds")
  unlink(summary_file)
})

test_that("calculatePRS handles multiple chromosomes", {
  # Skip if test files not available
  mdlfile <- system.file("extdata", "PGS000001_test.txt", package = "PRSfromPGS")
  skip_if(mdlfile == "", "Test model file PGS000001_test.txt not available")

  vcffile3 <- system.file("extdata", "chr3.vcf.gz", package = "PRSfromPGS")
  vcffile4 <- system.file("extdata", "chr4.vcf.gz", package = "PRSfromPGS")
  skip_if(vcffile3 == "" || vcffile4 == "",
          "Test VCF files chr3.vcf.gz and chr4.vcf.gz not available")

  # Create temp directory
  temp_dir_path <- tempdir()

  # Prepare multiple VCF files
  prepareVCFfiles(c(vcffile3, vcffile4), output_dir = temp_dir_path, verbose = FALSE)

  # Calculate PRS across chromosomes 3 and 4
  mdl <- readPGSmodel(mdlfile, verbose = FALSE)
  res <- suppressWarnings(calculatePRS(mdl, c(3, 4), temp_dir_path, verbose = FALSE))

  # Verify result
  expect_type(res, "list")
  expect_true("prs" %in% names(res))
  if (!is.null(res$prs)) {
    expect_type(res$prs, "double")
    expect_true(length(res$prs) > 0)
    expect_true(all(!is.na(res$prs)))
  }

  # Clean up
  unlink(file.path(temp_dir_path, "chromosome3VCFsummary.rds"))
  unlink(file.path(temp_dir_path, "chromosome4VCFsummary.rds"))
})

test_that("calculatePRS works with parallel processing", {
  # Skip if required packages not available
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  skip_if_not_installed("parallelly")

  # Skip if test files not available
  mdlfile <- system.file("extdata", "PGS000001_test.txt", package = "PRSfromPGS")
  skip_if(mdlfile == "", "Test model file PGS000001_test.txt not available")

  vcffile3 <- system.file("extdata", "chr3.vcf.gz", package = "PRSfromPGS")
  vcffile4 <- system.file("extdata", "chr4.vcf.gz", package = "PRSfromPGS")
  skip_if(vcffile3 == "" || vcffile4 == "",
          "Test VCF files chr3.vcf.gz and chr4.vcf.gz not available")

  # Create temp directory
  temp_dir_path <- tempdir()

  # Prepare multiple VCF files (need multiple chromosomes for parallel to be meaningful)
  prepareVCFfiles(c(vcffile3, vcffile4), output_dir = temp_dir_path, verbose = FALSE)

  # Calculate PRS with parallel processing
  mdl <- readPGSmodel(mdlfile, verbose = FALSE)
  res_parallel <- suppressWarnings(
    calculatePRS(mdl, c(3, 4), temp_dir_path, parallel = TRUE, n_cores = 2, verbose = FALSE)
  )

  # Calculate PRS without parallel processing for comparison
  res_sequential <- suppressWarnings(
    calculatePRS(mdl, c(3, 4), temp_dir_path, parallel = FALSE, verbose = FALSE)
  )

  # Verify parallel result structure
  expect_type(res_parallel, "list")
  expect_true("prs" %in% names(res_parallel))

  # Results should be identical whether parallel or sequential
  if (!is.null(res_parallel$prs) && !is.null(res_sequential$prs)) {
    expect_equal(res_parallel$prs, res_sequential$prs,
                 info = "Parallel and sequential PRS should produce identical results")
  }

  # Clean up
  unlink(file.path(temp_dir_path, "chromosome3VCFsummary.rds"))
  unlink(file.path(temp_dir_path, "chromosome4VCFsummary.rds"))
})

test_that("fitPRSmodels works with parallel processing", {
  # Skip if required packages not available
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  skip_if_not_installed("parallelly")

  # Skip if test files not available
  vcffile3 <- system.file("extdata", "chr3.vcf.gz", package = "PRSfromPGS")
  vcffile4 <- system.file("extdata", "chr4.vcf.gz", package = "PRSfromPGS")
  skip_if(vcffile3 == "" || vcffile4 == "",
          "Test VCF files chr3.vcf.gz and chr4.vcf.gz not available")

  # Find available test model files
  model_files <- c(
    system.file("extdata", "PGS000001_test.txt", package = "PRSfromPGS"),
    system.file("extdata", "PGS000002_test.txt", package = "PRSfromPGS")
  )
  model_files <- model_files[file.exists(model_files) & model_files != ""]
  skip_if(length(model_files) < 2, "Need at least 2 test model files for parallel test")

  # Create temp directory
  temp_dir_path <- tempdir()

  # Prepare VCF files
  prepareVCFfiles(c(vcffile3, vcffile4), output_dir = temp_dir_path, verbose = FALSE)

  # Fit models with parallel processing
  res_parallel <- suppressWarnings(
    fitPRSmodels(model_files, c(3, 4), temp_dir_path, parallel = TRUE, n_cores = 2, verbose = FALSE)
  )

  # Fit models without parallel processing for comparison
  res_sequential <- suppressWarnings(
    fitPRSmodels(model_files, c(3, 4), temp_dir_path, parallel = FALSE, verbose = FALSE)
  )

  # Verify parallel result structure
  expect_type(res_parallel, "list")
  expect_true("prs" %in% names(res_parallel))

  # Results should be identical whether parallel or sequential
  if (!is.null(res_parallel$prs) && !is.null(res_sequential$prs)) {
    expect_equal(res_parallel$prs, res_sequential$prs,
                 info = "Parallel and sequential fitPRSmodels should produce identical results")
  }

  # Clean up
  unlink(file.path(temp_dir_path, "chromosome3VCFsummary.rds"))
  unlink(file.path(temp_dir_path, "chromosome4VCFsummary.rds"))
})
