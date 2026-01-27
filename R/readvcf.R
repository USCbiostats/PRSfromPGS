#' Create Tabix Index Files for VCF Files
#'
#' @param vcf_files Character vector of paths to VCF files (must be bgzipped)
#' @param ncores Integer; number of cores for parallel processing. Default is 1 (sequential).
#' @param verbose Logical; if TRUE, display informational messages. Default is TRUE.
#' @return A data frame with file paths and indexing status
#' @export
#' @examples
#' \dontrun{
#' vcf_files <- c("data/file1.vcf.gz", "data/file2.vcf.gz")
#' create_tabix_files(vcf_files)
#' # Use parallel processing with 4 cores
#' create_tabix_files(vcf_files, ncores = 4)
#' }
create_tabix_files <- function(vcf_files, ncores = 1, verbose = TRUE) {
  # Helper function to index a single file
  index_single_file <- function(vcf_file) {
    if (!file.exists(vcf_file)) {
      return(list(indexed = FALSE, error = "File does not exist"))
    }

    tryCatch({
      indexTabix(vcf_file, format = "vcf")
      list(indexed = TRUE, error = NA_character_)
    }, error = function(e) {
      list(indexed = FALSE, error = as.character(e$message))
    })
  }

  # Process files (parallel or sequential)
  if (ncores > 1 && length(vcf_files) > 1) {
    if (!requireNamespace("BiocParallel", quietly = TRUE)) {
      warning("BiocParallel not available, falling back to sequential processing")
      ncores <- 1
    } else {
      if (verbose) message(sprintf("Creating tabix indices for %d files using %d cores...", length(vcf_files), ncores))

      bp_param <- BiocParallel::MulticoreParam(workers = ncores)
      results_list <- BiocParallel::bplapply(vcf_files, index_single_file, BPPARAM = bp_param)

      results <- data.frame(
        file = vcf_files,
        indexed = vapply(results_list, function(x) x$indexed, logical(1)),
        error = vapply(results_list, function(x) x$error, character(1)),
        stringsAsFactors = FALSE
      )

      n_success <- sum(results$indexed)
      if (verbose) message(sprintf("\nIndexing complete: %d/%d files successfully indexed", n_success, length(vcf_files)))

      return(results)
    }
  }

  # Sequential processing (fallback or ncores = 1)
  results <- data.frame(
    file = vcf_files,
    indexed = FALSE,
    error = NA_character_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(vcf_files)) {
    vcf_file <- vcf_files[i]
    if (verbose) message(sprintf("Creating tabix index for: %s", vcf_file))

    res <- index_single_file(vcf_file)
    results$indexed[i] <- res$indexed
    results$error[i] <- res$error

    if (verbose) {
      if (res$indexed) {
        message(sprintf("Successfully created index for: %s", vcf_file))
      } else {
        message(sprintf("Error indexing %s: %s", vcf_file, res$error))
      }
    }
  }

  n_success <- sum(results$indexed)
  if (verbose) message(sprintf("\nIndexing complete: %d/%d files successfully indexed", n_success, length(vcf_files)))

  return(results)
}

#' Read VCF Header Information
#'
#' @param vcf_file Character string with path to VCF file (must be bgzipped)
#' @param check_tabix Logical; if TRUE, checks for tabix index file
#' @param verbose Logical; if TRUE, display informational messages. Default is TRUE.
#' @return A VCFHeader object containing header information
#' @export
#' @examples
#' \dontrun{
#' vcf_file <- "data/axiom_acs_aus_nf_chr21.vcf.gz"
#' header <- getvcfheader(vcf_file)
#' # Access header components
#' samples(header)
#' meta(header)
#' geno(header)
#' }
getvcfheader <- function(vcf_file, check_tabix = TRUE, verbose = TRUE) {
  # Check if VCF file exists
  if (!file.exists(vcf_file)) {
    stop(sprintf("VCF file does not exist: %s", vcf_file))
  }

  # Check for tabix index
  tabix_file <- paste0(vcf_file, ".tbi")
  if (check_tabix) {
    if (!file.exists(tabix_file)) {
      warning(sprintf("Tabix index file not found: %s\nConsider running create_tabix_files() first.", tabix_file))
    } else {
      if (verbose) message(sprintf("Found tabix index: %s", tabix_file))
    }
  }

  # Read VCF header
  tryCatch({
    if (verbose) message(sprintf("Reading VCF header from: %s", vcf_file))
    header <- scanVcfHeader(vcf_file)
    if (verbose) message("Header successfully read")

    # Print summary information
    n_samples <- length(samples(header))
    if (verbose) message(sprintf("Number of samples: %d", n_samples))

    return(header)
  }, error = function(e) {
    stop(sprintf("Error reading VCF header from %s: %s", vcf_file, e$message))
  })
}

#' Read SNP Information Without Genotype Data
#'
#' @param vcf_file Character string with path to VCF file (must be bgzipped)
#' @param genome Character string specifying genome build (e.g., "hg19", "hg38")
#' @param verbose Logical; if TRUE, display informational messages. Default is TRUE.
#' @return A VCF object containing SNP information (position, alleles, INFO fields) without genotype data
#' @export
#' @examples
#' \dontrun{
#' vcf_file <- "data/axiom_acs_aus_nf_chr21.vcf.gz"
#' snps <- getsnpinfo(vcf_file, genome = "hg19")
#' # Access SNP information
#' rowRanges(snps)  # Genomic positions and alleles
#' info(snps)       # INFO fields (AF, MAF, R2, etc.)
#' }
getsnpinfo <- function(vcf_file, genome = "hg19", verbose = TRUE) {
  # Check if VCF file exists
  if (!file.exists(vcf_file)) {
    stop(sprintf("VCF file does not exist: %s", vcf_file))
  }

  # Check for tabix index
  tabix_file <- paste0(vcf_file, ".tbi")
  if (!file.exists(tabix_file)) {
    warning(sprintf("Tabix index file not found: %s\nConsider running create_tabix_files() first.", tabix_file))
  }

  # Create ScanVcfParam to exclude genotype data
  param <- ScanVcfParam(geno = NA)

  # Read VCF without genotype data
  tryCatch({
    if (verbose) {
      message(sprintf("Reading SNP information from: %s", vcf_file))
      message("Excluding genotype data to reduce memory usage")
    }

    vcf <- readVcf(vcf_file, genome = genome, param = param)

    n_snps <- nrow(vcf)
    if (verbose) message(sprintf("Successfully read %d SNPs", n_snps))

    return(vcf)
  }, error = function(e) {
    stop(sprintf("Error reading SNPs from %s: %s", vcf_file, e$message))
  })
}

#' Create SNP Summary Data Frame
#'
#' @param vcf_snps VCF object from getsnpinfo() function
#' @param verbose Logical; if TRUE, display informational messages. Default is TRUE.
#' @return A list containing snp_df (data frame with SNP ID, chromosome, position, reference allele, and alternate allele) and rr (rowRanges object)
#' @export
#' @examples
#' \dontrun{
#' vcf_file <- "data/axiom_acs_aus_nf_chr21.vcf.gz"
#' snps <- getsnpinfo(vcf_file)
#' result <- summarizesnpinfo(snps)
#' snp_df <- result$snp_df
#' rr <- result$rr
#' }
summarizesnpinfo <- function(vcf_snps, verbose = TRUE) {
  # Extract rowRanges
  rr <- rowRanges(vcf_snps)

  # Create data frame
  snp_df <- data.frame(
    SNP_ID = names(rr),
    seqnames = as.character(seqnames(rr)),
    position = start(ranges(rr)),
    REF = as.character(ref(vcf_snps)),
    ALT = vapply(alt(vcf_snps), paste, character(1), collapse = ","),
    stringsAsFactors = FALSE
  )

  if (verbose) message(sprintf("Created summary data frame with %d SNPs", nrow(snp_df)))

  return(list(snp_df = snp_df, rr = rr))
}

#' Read Dosage (DS) Genotype Information for SNPs
#'
#' @param vcf_summary List returned from summarizevcf() function containing vcf_file, genome, snp_summary, and sample_ids
#' @param snp_indices Integer vector of SNP indices to retrieve (1-based). Can be a single index or multiple indices.
#' @param verbose Logical; if TRUE, display informational messages during execution. Default is FALSE.
#' @return If single SNP requested, a numeric vector of DS values. If multiple SNPs, a matrix (rows = SNPs, columns = samples).
#' @export
#' @examples
#' \dontrun{
#' vcf_file <- "data/axiom_acs_aus_nf_chr21.vcf.gz"
#' vcf_summary <- summarizevcf(vcf_file)
#'
#' # Read a single SNP
#' dosage <- getvcfsnp(vcf_summary, snp_indices = 1)
#'
#' # Read multiple SNPs efficiently
#' dosages <- getvcfsnp(vcf_summary, snp_indices = c(1, 2, 5, 10))
#'
#' # Read a range of SNPs with verbose output
#' dosages <- getvcfsnp(vcf_summary, snp_indices = 1:100, verbose = TRUE)
#' }
getvcfsnp <- function(vcf_summary, snp_indices, verbose = FALSE) {
  # Validate input
  if (!is.list(vcf_summary) || !all(c("vcf_file", "genome", "snp_summary") %in% names(vcf_summary))) {
    stop("vcf_summary must be a list returned from summarizevcf() containing vcf_file, genome, and snp_summary")
  }

  # Extract components from vcf_summary
  vcf_file <- vcf_summary$vcf_file
  genome <- vcf_summary$genome
  row_ranges <- vcf_summary$snp_summary$rr

  # Check if VCF file exists
  if (!file.exists(vcf_file)) {
    stop(sprintf("VCF file does not exist: %s", vcf_file))
  }

  # Check for tabix index
  tabix_file <- paste0(vcf_file, ".tbi")
  if (!file.exists(tabix_file)) {
    stop(sprintf("Tabix index file not found: %s\nRun create_tabix_files() first.", tabix_file))
  }

  tryCatch({
    # Validate all indices are in range
    if (any(snp_indices < 1 | snp_indices > length(row_ranges))) {
      stop(sprintf("Some SNP indices are out of range. File contains %d SNPs.", length(row_ranges)))
    }

    # Check for duplicate indices
    if (any(duplicated(snp_indices))) {
      duplicates <- unique(snp_indices[duplicated(snp_indices)])
      stop(sprintf("Duplicate SNP indices found: %s. Each SNP index must be unique.",
                   paste(duplicates, collapse = ", ")))
    }

    if (verbose) {
      message(sprintf("Reading DS genotype data for %d SNP(s)", length(snp_indices)))
    }

    # Get the genomic ranges for the specified SNPs
    selected_ranges <- row_ranges[snp_indices]
    original_n_indices <- length(snp_indices)

    # Create data frame of requested SNPs from snp_summary$snp_df
    snp_df <- vcf_summary$snp_summary$snp_df
    requested_snps_df <- snp_df[snp_indices, ]

    # Check for duplicate start positions and remove them
    start_positions <- start(ranges(selected_ranges))
    if (any(duplicated(start_positions))) {
      unique_indices <- !duplicated(start_positions)
      n_duplicates <- sum(!unique_indices)
      duplicate_positions <- unique(start_positions[duplicated(start_positions)])

      warning(sprintf("Found %d duplicate start position(s): %s. Keeping only the first occurrence of each position.",
                      n_duplicates,
                      paste(duplicate_positions, collapse = ", ")))

      # Filter to keep only unique start positions
      # Note: requested_snps_df is NOT filtered here - it stays based on unique snp_indices
      selected_ranges <- selected_ranges[unique_indices]
      snp_indices <- snp_indices[unique_indices]

      if (verbose) {
        message(sprintf("After removing duplicates: %d unique SNP(s) remaining", length(snp_indices)))
      }
    }

    # Create ScanVcfParam to read only these SNPs' DS genotype field
    param <- ScanVcfParam(
      which = selected_ranges,
      geno = "DS"
    )

    # Read the VCF with only DS genotype data for these SNPs
    # Note: readVcf may return multiple SNPs if they fall within the same ranges
    vcf <- readVcf(vcf_file, genome = genome, param = param)

    # Get SNP information from the VCF
    vcf_rowranges <- rowRanges(vcf)
    actual_snp_ids <- names(vcf_rowranges)

    # Extract DS values (matrix: SNPs x Samples)
    ds_matrix <- geno(vcf)$DS

    # Create actual SNPs data frame for filtering
    actual_snps_df <- data.frame(
      SNP_ID = actual_snp_ids,
      seqnames = as.character(seqnames(vcf_rowranges)),
      position = start(ranges(vcf_rowranges)),
      REF = as.character(ref(vcf)),
      ALT = vapply(alt(vcf), paste, character(1), collapse = ","),
      stringsAsFactors = FALSE
    )

    if (verbose) {
      message(sprintf("VCF returned DS genotype data: %d SNPs x %d samples", nrow(ds_matrix), ncol(ds_matrix)))
    }

    # Create composite SNP IDs from seqnames:position:REF:ALT
    composite_snp_ids <- paste(actual_snps_df$seqnames,
                                actual_snps_df$position,
                                actual_snps_df$REF,
                                actual_snps_df$ALT,
                                sep = ":")

    # Set row names to the composite SNP IDs
    rownames(ds_matrix) <- composite_snp_ids

    # Filter dosage matrix to keep only the requested SNPs
    # Match based on SNP_ID, position, REF, and ALT using vectorized matching
    requested_keys <- paste(requested_snps_df$SNP_ID,
                            requested_snps_df$position,
                            requested_snps_df$REF,
                            requested_snps_df$ALT, sep = ":")
    actual_keys <- paste(actual_snps_df$SNP_ID,
                         actual_snps_df$position,
                         actual_snps_df$REF,
                         actual_snps_df$ALT, sep = ":")
    keep_rows <- actual_keys %in% requested_keys

    ds_matrix <- ds_matrix[keep_rows, , drop = FALSE]

    if (verbose) {
      message(sprintf("Filtered to requested SNPs: %d SNPs x %d samples", nrow(ds_matrix), ncol(ds_matrix)))
    }

    # Return dosage data
    if (original_n_indices == 1 && nrow(ds_matrix) == 1) {
      return(ds_matrix[1, ])
    } else {
      return(ds_matrix)
    }
  }, error = function(e) {
    stop(sprintf("Error reading SNPs from %s: %s", vcf_file, e$message))
  })
}

#' Summarize VCF File
#'
#' @param vcf_file Character string with path to VCF file (must be bgzipped)
#' @param genome Character string specifying genome build (e.g., "hg19", "hg38")
#' @param verbose Logical; if TRUE, display informational messages. Default is TRUE.
#' @return A list containing:
#'   - vcf_file: The VCF file name
#'   - genome: The genome build used
#'   - snp_summary: The list returned from summarizesnpinfo (contains snp_df and rr)
#'   - sample_ids: Vector of sample IDs from the VCF file
#' @export
#' @examples
#' \dontrun{
#' vcf_file <- "data/axiom_acs_aus_nf_chr21.vcf.gz"
#' summary <- summarizevcf(vcf_file)
#' summary$vcf_file
#' summary$genome
#' summary$snp_summary$snp_df
#' summary$sample_ids
#' }
summarizevcf <- function(vcf_file, genome = "hg19", verbose = TRUE) {
  # Check if VCF file exists
  if (!file.exists(vcf_file)) {
    stop(sprintf("VCF file does not exist: %s", vcf_file))
  }

  # Read SNP information
  if (verbose) message(sprintf("Summarizing VCF file: %s", vcf_file))
  vcf_snps <- getsnpinfo(vcf_file, genome = genome, verbose = verbose)

  # Get SNP summary
  snp_summary <- summarizesnpinfo(vcf_snps, verbose = verbose)

  # Extract sample IDs from VCF object header (avoids re-reading the file)
  sample_ids <- samples(header(vcf_snps))
  if (verbose) message(sprintf("Found %d samples", length(sample_ids)))

  # Return combined information
  return(list(
    vcf_file = vcf_file,
    genome = genome,
    snp_summary = snp_summary,
    sample_ids = sample_ids
  ))
}

#' Prepare VCF Files for Analysis
#'
#' @param vcf_files Character vector of paths to VCF files (must be bgzipped)
#' @param genome Character string specifying genome build (e.g., "hg19", "hg38"). Applied to all files.
#' @param tabix Logical; if TRUE, create tabix index files for the VCF files. Default is FALSE.
#' @param ncores Integer; number of cores for parallel processing. Default is 1 (sequential).
#' @param verbose Logical; if TRUE, display informational messages. Default is TRUE.
#' @return Nothing (NULL invisibly)
#' @export
#' @examples
#' \dontrun{
#' vcf_files <- c("data/chr21.vcf.gz", "data/chr22.vcf.gz")
#' prepareVCFfiles(vcf_files, genome = "hg19", tabix = TRUE)
#' # Use parallel processing with 4 cores
#' prepareVCFfiles(vcf_files, genome = "hg19", ncores = 4)
#' # Suppress messages
#' prepareVCFfiles(vcf_files, verbose = FALSE)
#' }
prepareVCFfiles <- function(vcf_files, genome = "hg19", tabix = FALSE, ncores = 1, verbose = TRUE) {
  # Validate inputs
  if (length(vcf_files) == 0) {
    stop("vcf_files vector is empty")
  }

  # Create tabix index files if requested
  if (tabix) {
    if (verbose) message("Creating tabix index files...")
    create_tabix_files(vcf_files, ncores = ncores, verbose = verbose)
  }

  # Helper function to process a single VCF file
  process_single_vcf <- function(vcf_file, genome, verbose) {
    # Summarize VCF
    vcf_summary <- summarizevcf(vcf_file, genome = genome, verbose = verbose)

    # Save as RDS file using chromosome name from the VCF data
    chr_name <- vcf_summary$snp_summary$snp_df$seqnames[1]
    rds_filename <- sprintf("chromosome%sVCFsummary.rds", chr_name)
    saveRDS(vcf_summary, file = rds_filename)

    return(rds_filename)
  }

  # Process VCF files (parallel or sequential)
  if (ncores > 1 && length(vcf_files) > 1) {
    # Parallel processing using BiocParallel
    if (!requireNamespace("BiocParallel", quietly = TRUE)) {
      warning("BiocParallel not available, falling back to sequential processing")
      ncores <- 1
    } else {
      if (verbose) message(sprintf("Processing %d VCF files in parallel using %d cores...", length(vcf_files), ncores))

      # Set up parallel backend
      bp_param <- BiocParallel::MulticoreParam(workers = ncores)

      # Process files in parallel (verbose=FALSE in workers to avoid interleaved output)
      results <- BiocParallel::bplapply(vcf_files, function(vcf_file) {
        tryCatch({
          rds_file <- process_single_vcf(vcf_file, genome, verbose = FALSE)
          list(success = TRUE, file = vcf_file, rds = rds_file, error = NA)
        }, error = function(e) {
          list(success = FALSE, file = vcf_file, rds = NA, error = e$message)
        })
      }, BPPARAM = bp_param)

      # Report results
      n_success <- sum(vapply(results, function(x) x$success, logical(1)))
      if (verbose) message(sprintf("\n=== Successfully processed %d/%d VCF files ===", n_success, length(vcf_files)))

      # Report any errors
      for (res in results) {
        if (!res$success) {
          warning(sprintf("Failed to process %s: %s", res$file, res$error))
        } else {
          if (verbose) message(sprintf("Saved summary to: %s", res$rds))
        }
      }

      return(invisible(NULL))
    }
  }

  # Sequential processing (fallback or ncores = 1)
  for (i in seq_along(vcf_files)) {
    vcf_file <- vcf_files[i]

    if (verbose) message(sprintf("\n=== Processing VCF file %d/%d: %s ===", i, length(vcf_files), vcf_file))

    rds_filename <- process_single_vcf(vcf_file, genome, verbose)
    if (verbose) message(sprintf("Saved summary to: %s", rds_filename))
  }

  if (verbose) message(sprintf("\n=== Successfully processed %d VCF files ===", length(vcf_files)))

  invisible(NULL)
}
