test_that("prepareVCFfiles creates correct VCF summary file", {
  # Skip if VCF file not available
  vcffile3 <- system.file("extdata", "chr3.vcf.gz", package = "PRSfromPGS")
  skip_if(vcffile3 == "", "Test VCF file chr3.vcf.gz not available")

  # Create temp directory for output
  temp_dir_path <- tempdir()

  # Run prepareVCFfiles
  prepareVCFfiles(vcffile3, output_dir = temp_dir_path, verbose = FALSE)

  # Read the generated summary
  generated_summary_path <- file.path(temp_dir_path, "chromosome3VCFsummary.rds")
  expect_true(file.exists(generated_summary_path),
              "prepareVCFfiles should create chromosome3VCFsummary.rds")

  res <- readRDS(generated_summary_path)

  # Read the reference summary
  summaryfile <- system.file("extdata", "chromosome3VCFsummary.rds",
                              package = "PRSfromPGS")
  skip_if(summaryfile == "", "Reference chromosome3VCFsummary.rds not available")

  res2 <- readRDS(summaryfile)

  # Verify structure of the VCF summary
  expect_type(res, "list")
  expect_true("vcf_file" %in% names(res),
              "VCF summary should contain 'vcf_file' element")
  expect_true("sample_ids" %in% names(res),
              "VCF summary should contain 'sample_ids' element")
  expect_true("snp_summary" %in% names(res),
              "VCF summary should contain 'snp_summary' element")
  expect_true("genome" %in% names(res),
              "VCF summary should contain 'genome' element")

  # Check snp_summary structure (it's a list with snp_df and rr)
  expect_type(res$snp_summary, "list")
  expect_true("snp_df" %in% names(res$snp_summary),
              "snp_summary should contain 'snp_df' element")
  expect_s3_class(res$snp_summary$snp_df, "data.frame")
  expect_true(nrow(res$snp_summary$snp_df) > 0,
              "snp_df should contain variant data")

  # Check that key elements match reference
  expect_equal(names(res), names(res2),
               info = "Generated and reference summaries should have same structure")
  expect_equal(nrow(res$snp_summary$snp_df), nrow(res2$snp_summary$snp_df),
               info = "Generated and reference should have same number of variants")
  expect_equal(length(res$sample_ids), length(res2$sample_ids),
               info = "Generated and reference should have same number of samples")

  # Clean up
  unlink(generated_summary_path)
})

test_that("prepareVCFfiles works with default output_dir", {
  # Skip if VCF file not available
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

  # Run prepareVCFfiles with default output_dir (current directory)
  prepareVCFfiles(vcffile3, verbose = FALSE)

  # Check that output file was created in current directory
  expected_file <- file.path(temp_dir, "chromosome3VCFsummary.rds")
  expect_true(file.exists(expected_file),
              "prepareVCFfiles should create output in current directory when output_dir not specified")

  res <- readRDS(expected_file)
  expect_type(res, "list")
  expect_true("vcf_file" %in% names(res))
})

test_that("prepareVCFfiles handles multiple VCF files", {
  # Skip if VCF files not available
  vcffile3 <- system.file("extdata", "chr3.vcf.gz", package = "PRSfromPGS")
  vcffile4 <- system.file("extdata", "chr4.vcf.gz", package = "PRSfromPGS")
  skip_if(vcffile3 == "" || vcffile4 == "",
          "Test VCF files chr3.vcf.gz and chr4.vcf.gz not available")

  # Create temp directory for output
  temp_dir_path <- tempdir()

  # Run prepareVCFfiles with multiple files
  prepareVCFfiles(c(vcffile3, vcffile4), output_dir = temp_dir_path, verbose = FALSE)

  # Check that both summary files were created
  chr3_summary <- file.path(temp_dir_path, "chromosome3VCFsummary.rds")
  chr4_summary <- file.path(temp_dir_path, "chromosome4VCFsummary.rds")

  expect_true(file.exists(chr3_summary),
              "prepareVCFfiles should create chromosome3VCFsummary.rds")
  expect_true(file.exists(chr4_summary),
              "prepareVCFfiles should create chromosome4VCFsummary.rds")

  # Verify both summaries are valid
  res3 <- readRDS(chr3_summary)
  res4 <- readRDS(chr4_summary)

  expect_type(res3, "list")
  expect_type(res4, "list")
  expect_true("snp_summary" %in% names(res3))
  expect_true("snp_summary" %in% names(res4))

  # Clean up
  unlink(chr3_summary)
  unlink(chr4_summary)
})
