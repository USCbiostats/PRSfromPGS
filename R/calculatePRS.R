#' Read PGS Model from Local File
#'
#' Reads a Polygenic Score (PGS) model from a local file in PGS Catalog format.
#' The file format should be:
#'   - Comment lines starting with # (skipped)
#'   - Header line with column names
#'   - Tab-delimited data rows
#'
#' @param file_path Character string specifying the path to the PGS model file
#' @param verbose Logical; if TRUE (default), display informational messages
#' @return A data frame containing the PGS model data
#' @export
#' @examples
#' \dontrun{
#' # Read a PGS model file
#' pgs_model <- readPGSmodel("pgsmodel.txt")
#'
#' # Read without messages
#' pgs_model <- readPGSmodel("pgsmodel.txt", verbose = FALSE)
#' }
readPGSmodel <- function(file_path, verbose = TRUE) {

  # Validate input
  if (!is.character(file_path) || length(file_path) != 1) {
    stop("file_path must be a single character string")
  }

  # Check if file exists
  if (!file.exists(file_path)) {
    stop(sprintf("File does not exist: %s", file_path))
  }

  if (verbose) {
    message(sprintf("Reading PGS model from: %s", file_path))
  }

  # Read file using data.table::fread for speed
  # Count comment lines to skip
  pgs_data <- tryCatch({
    # Read first lines to count comment lines starting with ##
    con <- file(file_path, "r")
    n_skip <- 0
    while (TRUE) {
      line <- readLines(con, n = 1)
      if (length(line) == 0 || !grepl("^##", line)) break
      n_skip <- n_skip + 1
    }
    close(con)

    # Read data using fread (much faster than read.table)
    dt <- data.table::fread(file_path,
                            header = TRUE,
                            sep = "\t",
                            skip = n_skip,
                            quote = "",
                            check.names = FALSE)
    as.data.frame(dt)
  }, error = function(e) {
    stop(sprintf("Error reading file %s: %s", file_path, e$message))
  })

  # Ensure result is a data frame
  if (!is.data.frame(pgs_data)) {
    stop("File did not produce a data frame. Check file format.")
  }

  if (verbose) {
    message(sprintf("Successfully read %d rows and %d columns",
                   nrow(pgs_data), ncol(pgs_data)))
    message(sprintf("Columns: %s", paste(names(pgs_data), collapse = ", ")))
  }

  # Create chr_name column if it doesn't exist
  if (!"chr_name" %in% names(pgs_data)) {
    if ("hm_chr" %in% names(pgs_data)) {
      pgs_data$chr_name <- pgs_data$hm_chr
      if (verbose) {
        message("Column 'chr_name' not found - using hm_chr values")
      }
    }
  }

  # Create chr_position column if it doesn't exist
  if (!"chr_position" %in% names(pgs_data)) {
    if ("hm_pos" %in% names(pgs_data)) {
      pgs_data$chr_position <- pgs_data$hm_pos
      if (verbose) {
        message("Column 'chr_position' not found - using hm_pos values")
      }
    }
  }

  # Remove rows where chr_position is NA or 0
  if ("chr_position" %in% names(pgs_data)) {
    invalid_pos <- is.na(pgs_data$chr_position) | pgs_data$chr_position == 0
    n_invalid <- sum(invalid_pos, na.rm = TRUE)
    if (n_invalid > 0) {
      pgs_data <- pgs_data[!invalid_pos, ]
      if (verbose) {
        message(sprintf("Removed %d rows with NA or 0 chr_position", n_invalid))
        message(sprintf("Remaining rows: %d", nrow(pgs_data)))
      }
    }
  }

  # Create other_allele column if it doesn't exist
  if (!"other_allele" %in% names(pgs_data)) {
    if ("hm_inferOtherAllele" %in% names(pgs_data)) {
      pgs_data$other_allele <- pgs_data$hm_inferOtherAllele
      if (verbose) {
        message("Column 'other_allele' not found - using hm_inferOtherAllele values")
      }
    } else {
      pgs_data$other_allele <- NA_character_
      if (verbose) {
        message("Column 'other_allele' not found - created with NA values")
      }
    }
  }

  # For rows where rsID is empty, create composite ID
  if ("rsID" %in% names(pgs_data)) {
    empty_rsID <- pgs_data$rsID == "" | is.na(pgs_data$rsID)
    n_empty <- sum(empty_rsID, na.rm = TRUE)

    if (n_empty > 0) {
      # Check required columns for composite ID
      required_cols <- c("chr_name", "chr_position", "effect_allele")
      missing_cols <- setdiff(required_cols, names(pgs_data))

      if (length(missing_cols) > 0) {
        warning(sprintf("Cannot create composite IDs: missing columns %s",
                       paste(missing_cols, collapse = ", ")))
      } else {
        # Create composite IDs (use 'X' for NA other_allele)
        pgs_data$rsID[empty_rsID] <- paste(
          pgs_data$chr_name[empty_rsID],
          pgs_data$chr_position[empty_rsID],
          pgs_data$effect_allele[empty_rsID],
          ifelse(is.na(pgs_data$other_allele[empty_rsID]), "X", pgs_data$other_allele[empty_rsID]),
          sep = ":"
        )

        if (verbose) {
          message(sprintf("Created composite IDs for %d rows with empty rsID", n_empty))
        }
      }
    }
  } else {
    # Create rsID column with composite IDs if it doesn't exist
    required_cols <- c("chr_name", "chr_position", "effect_allele")
    missing_cols <- setdiff(required_cols, names(pgs_data))

    if (length(missing_cols) > 0) {
      warning(sprintf("Cannot create composite IDs: missing columns %s",
                     paste(missing_cols, collapse = ", ")))
      pgs_data$rsID <- ""
    } else {
      # Create composite IDs (use 'X' for NA other_allele)
      pgs_data$rsID <- paste(
        pgs_data$chr_name,
        pgs_data$chr_position,
        pgs_data$effect_allele,
        ifelse(is.na(pgs_data$other_allele), "X", pgs_data$other_allele),
        sep = ":"
      )
      if (verbose) {
        message(sprintf("Column 'rsID' not found - created composite IDs for %d rows", nrow(pgs_data)))
      }
    }
  }

  # Return only the required columns in specified order
  required_cols <- c("rsID", "chr_name", "chr_position", "effect_allele", "other_allele", "effect_weight")
  missing_cols <- setdiff(required_cols, names(pgs_data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  return(pgs_data[, required_cols])
}

#' Process a Single Chromosome for PRS Calculation
#'
#' Internal helper function that processes one chromosome for PRS calculation.
#' This function is designed to be called in parallel across chromosomes.
#'
#' @param chr Chromosome identifier
#' @param pgsmodel Data frame with PGS model data
#' @param batch_size Integer; maximum number of SNPs to read from VCF at once
#' @param vcf_dir Character string specifying directory containing VCF summary files
#' @return A list containing chr_prs, unmatched_count, unmatched_rsIDs, and status info
#' @keywords internal
processSingleChromosome <- function(chr, pgsmodel, batch_size, vcf_dir) {
  # Filter pgsmodel for current chromosome
  chr_pgsmodel <- pgsmodel[pgsmodel$chr_name == chr, ]

  if (nrow(chr_pgsmodel) == 0) {
    return(list(
      chr = chr,
      chr_prs = NULL,
      unmatched_count = 0,
      unmatched_rsIDs = character(0),
      ref_matches = 0,
      alt_matches = 0,
      status = "no_variants"
    ))
  }

  # Load VCF summary for this chromosome
  vcf_summary_file <- file.path(vcf_dir, sprintf("chromosome%sVCFsummary.rds", chr))

  if (!file.exists(vcf_summary_file)) {
    return(list(
      chr = chr,
      chr_prs = NULL,
      unmatched_count = 0,
      unmatched_rsIDs = character(0),
      ref_matches = 0,
      alt_matches = 0,
      status = "file_not_found"
    ))
  }

  vcf_summary <- tryCatch({
    readRDS(vcf_summary_file)
  }, error = function(e) {
    return(NULL)
  })

  if (is.null(vcf_summary)) {
    return(list(
      chr = chr,
      chr_prs = NULL,
      unmatched_count = 0,
      unmatched_rsIDs = character(0),
      ref_matches = 0,
      alt_matches = 0,
      status = "read_error"
    ))
  }

  # Find matching SNPs

snp_matches <- findSNPs(chr_pgsmodel, vcf_summary$snp_summary$snp_df)

  # Track unmatched SNPs
  unmatched_count <- length(snp_matches$no_matches)
  if (unmatched_count > 0) {
    unmatched_rsIDs <- chr_pgsmodel$rsID[snp_matches$no_matches]
  } else {
    unmatched_rsIDs <- character(0)
  }

  # Calculate PRS for ALT matches
  prs_alt <- NULL
  if (length(snp_matches$alt_matches) > 0) {
    prs_alt <- calculatePRSalt(snp_matches$alt_matches, chr_pgsmodel, vcf_summary, batch_size = batch_size, verbose = FALSE)
  }

  # Calculate PRS for REF matches
  prs_ref <- NULL
  if (length(snp_matches$ref_matches) > 0) {
    prs_ref <- calculatePRSref(snp_matches$ref_matches, chr_pgsmodel, vcf_summary, batch_size = batch_size, verbose = FALSE)
  }

  # Combine PRS from ALT and REF matches
  chr_prs <- NULL

  if (!is.null(prs_alt) && !is.null(prs_ref)) {
    if (!all(names(prs_alt) == names(prs_ref))) {
      prs_ref <- prs_ref[names(prs_alt)]
    }
    chr_prs <- prs_alt + prs_ref
  } else if (!is.null(prs_alt)) {
    chr_prs <- prs_alt
  } else if (!is.null(prs_ref)) {
    chr_prs <- prs_ref
  }

  return(list(
    chr = chr,
    chr_prs = chr_prs,
    unmatched_count = unmatched_count,
    unmatched_rsIDs = unmatched_rsIDs,
    ref_matches = length(snp_matches$ref_matches),
    alt_matches = length(snp_matches$alt_matches),
    status = "success"
  ))
}

#' Calculate Polygenic Risk Score Across Chromosomes
#'
#' Calculates PRS by processing variants from multiple chromosomes. For each chromosome,
#' the function loads VCF summary data, matches variants, and calculates PRS contributions
#' from both REF and ALT allele matches. The final PRS is the sum across all chromosomes.
#'
#' @param pgsmodel Data frame returned from readPGSmodel() containing PGS variant information.
#'   Expected columns: rsID, chr_name, chr_position, effect_allele, other_allele, effect_weight
#' @param chromosomes Integer vector of chromosome numbers to process. Default is 1:22 (autosomes).
#' @param vcf_dir Character string specifying directory containing VCF summary files. Default is "." (current directory).
#' @param batch_size Integer; maximum number of SNPs to read from VCF at once. Default is 1000.
#' @param parallel Logical; if TRUE, process chromosomes in parallel using future. Default is FALSE.
#' @param n_cores Integer; number of cores to use for parallel processing. Default is NULL (use all available).
#' @param verbose Logical; if TRUE (default), display informational messages during execution
#' @return A list containing:
#'   - prs: A numeric vector of PRS values for all samples
#'   - unmatched_by_chr: A named vector with counts of unmatched SNPs by chromosome
#'   - unmatched_rsIDs: A named list with rsIDs of unmatched SNPs by chromosome
#'   - excluded_chr_counts: A named vector with counts of SNPs not on processed chromosomes
#' @export
#' @examples
#' \dontrun{
#' # Read PGS model
#' pgs_model <- readPGSmodel("pgsmodel.txt")
#'
#' # Calculate PRS for all autosomes
#' result <- calculatePRS(pgs_model)
#' prs_values <- result$prs
#' unmatched_counts <- result$unmatched_by_chr
#' unmatched_ids <- result$unmatched_rsIDs
#' excluded_counts <- result$excluded_chr_counts
#'
#' # View unmatched SNPs for chromosome 1
#' result$unmatched_rsIDs[["1"]]
#'
#' # View SNPs on non-autosomal chromosomes (X, Y, MT, etc.)
#' result$excluded_chr_counts
#'
#' # Calculate PRS for specific chromosomes
#' result <- calculatePRS(pgs_model, chromosomes = c(1, 2, 3))
#'
#' # Calculate PRS quietly
#' result <- calculatePRS(pgs_model, verbose = FALSE)
#'
#' # Calculate PRS with custom batch size
#' result <- calculatePRS(pgs_model, batch_size = 500)
#'
#' # Calculate PRS in parallel
#' result <- calculatePRS(pgs_model, parallel = TRUE)
#'
#' # Calculate PRS in parallel with specific number of cores
#' result <- calculatePRS(pgs_model, parallel = TRUE, n_cores = 4)
#'
#' # Calculate PRS using VCF summaries from a specific directory
#' result <- calculatePRS(pgs_model, vcf_dir = "vcf_summaries")
#' }
calculatePRS <- function(pgsmodel, chromosomes = 1:22, vcf_dir = ".", batch_size = 1000, parallel = FALSE, n_cores = NULL, verbose = TRUE) {
  # Validate inputs
  if (!is.data.frame(pgsmodel)) {
    stop("pgsmodel must be a data frame returned from readPGSmodel()")
  }

  # Check required columns
  required_cols <- c("rsID", "chr_name", "chr_position", "effect_allele", "effect_weight")
  missing_cols <- setdiff(required_cols, names(pgsmodel))
  if (length(missing_cols) > 0) {
    stop(sprintf("pgsmodel is missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  # Identify SNPs not on chromosomes being processed
  excluded_snps <- pgsmodel[!pgsmodel$chr_name %in% chromosomes, ]
  excluded_chr_counts <- table(excluded_snps$chr_name)

  if (verbose) {
    message(sprintf("Calculating PRS using %d variants across %d chromosomes",
                   nrow(pgsmodel), length(chromosomes)))

    if (nrow(excluded_snps) > 0) {
      message(sprintf("Note: %d variants not on processed chromosomes will be excluded:", nrow(excluded_snps)))
      for (chr_name in names(excluded_chr_counts)) {
        message(sprintf("  Chromosome %s: %d variants", chr_name, excluded_chr_counts[chr_name]))
      }
    }
  }

  # Process chromosomes (parallel or sequential)
  if (parallel) {
    # Check for required packages
    required_pkgs <- c("future", "future.apply", "parallelly")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
      stop(sprintf("Packages '%s' are required for parallel processing. Install with: install.packages(c('future', 'future.apply', 'parallelly'))",
                   paste(missing_pkgs, collapse = "', '")))
    }

    # Set up parallel backend
    if (is.null(n_cores)) {
      n_cores <- parallelly::availableCores()
    }

    # Use multicore on Linux/macOS (forking), multisession on Windows
    if (.Platform$OS.type == "unix") {
      parallel_type <- "multicore"
      if (verbose) {
        message(sprintf("Processing chromosomes in parallel using %d cores (multicore/forking)", n_cores))
      }
      old_plan <- future::plan()
      on.exit(future::plan(old_plan), add = TRUE)
      future::plan(future::multicore, workers = n_cores)
    } else {
      parallel_type <- "multisession"
      if (verbose) {
        message(sprintf("Processing chromosomes in parallel using %d cores (multisession)", n_cores))
      }
      old_plan <- future::plan()
      on.exit(future::plan(old_plan), add = TRUE)
      future::plan(future::multisession, workers = n_cores)
    }

    # Process chromosomes in parallel
    chr_results <- future.apply::future_lapply(chromosomes, function(chr) {
      processSingleChromosome(chr, pgsmodel, batch_size, vcf_dir)
    }, future.seed = TRUE)

  } else {
    # Sequential processing
    chr_results <- lapply(chromosomes, function(chr) {
      if (verbose) {
        message(sprintf("\n=== Processing chromosome %s ===", chr))
      }

      result <- processSingleChromosome(chr, pgsmodel, batch_size, vcf_dir)

      if (verbose) {
        if (result$status == "no_variants") {
          message(sprintf("No variants found for chromosome %s in PGS model - skipping", chr))
        } else if (result$status == "file_not_found") {
          message(sprintf("VCF summary file not found for chromosome %s - skipping", chr))
        } else if (result$status == "read_error") {
          message(sprintf("Error reading VCF summary for chromosome %s - skipping", chr))
        } else {
          message(sprintf("  REF matches: %d", result$ref_matches))
          message(sprintf("  ALT matches: %d", result$alt_matches))
          message(sprintf("  No matches: %d", result$unmatched_count))
          if (!is.null(result$chr_prs)) {
            message(sprintf("Chromosome %s PRS contribution: mean=%.4f, sd=%.4f",
                           chr, mean(result$chr_prs), sd(result$chr_prs)))
          }
        }
      }

      result
    })
  }

  # Combine results from all chromosomes
  total_prs <- NULL
  sample_names <- NULL
  unmatched_counts <- integer()
  unmatched_rsIDs <- list()

  for (result in chr_results) {
    chr_key <- as.character(result$chr)
    unmatched_counts[chr_key] <- result$unmatched_count
    unmatched_rsIDs[[chr_key]] <- result$unmatched_rsIDs

    if (!is.null(result$chr_prs)) {
      if (is.null(total_prs)) {
        total_prs <- result$chr_prs
        sample_names <- names(result$chr_prs)
      } else {
        # Ensure sample names match before adding
        if (!all(names(result$chr_prs) == sample_names)) {
          result$chr_prs <- result$chr_prs[sample_names]
        }
        total_prs <- total_prs + result$chr_prs
      }
    }
  }

  # Final summary
  if (is.null(total_prs)) {
    warning("No PRS calculated - no matching variants found in any chromosome")
    return(list(
      prs = NULL,
      unmatched_by_chr = unmatched_counts,
      unmatched_rsIDs = unmatched_rsIDs,
      excluded_chr_counts = excluded_chr_counts
    ))
  }

  if (verbose) {
    message(sprintf("\n=== PRS Calculation Complete ==="))
    message(sprintf("Total samples: %d", length(total_prs)))
    message(sprintf("Final PRS: mean=%.4f, sd=%.4f, min=%.4f, max=%.4f",
                   mean(total_prs), sd(total_prs), min(total_prs), max(total_prs)))

    total_unmatched <- sum(unmatched_counts)
    if (total_unmatched > 0) {
      message(sprintf("Total unmatched SNPs across all chromosomes: %d", total_unmatched))
    }

    if (length(excluded_chr_counts) > 0) {
      message(sprintf("Total excluded SNPs (not on processed chromosomes): %d", sum(excluded_chr_counts)))
    }
  }

  return(list(
    prs = total_prs,
    unmatched_by_chr = unmatched_counts,
    unmatched_rsIDs = unmatched_rsIDs,
    excluded_chr_counts = excluded_chr_counts
  ))
}

#' Find SNPs in VCF Data
#'
#' @param prssnps Data frame with PRS SNP information containing columns:
#'   - rsID: SNP identifier
#'   - chr_position: Chromosome position
#'   - effect_weight: Effect weight for PRS calculation
#'   - effect_allele: Effect allele
#' @param vcfsnpinfo Data frame with VCF SNP information (format from summarizesnpinfo) containing columns:
#'   - SNP_ID: SNP identifier
#'   - seqnames: Chromosome name
#'   - position: Position on chromosome
#'   - REF: Reference allele
#'   - ALT: Alternate allele(s)
#' @return A list containing three elements:
#'   - ref_matches: List where each element is a vector containing the prssnps index followed by matching vcfsnpinfo indices (position and REF match)
#'   - alt_matches: List where each element is a vector containing the prssnps index followed by matching vcfsnpinfo indices (position and ALT match)
#'   - no_matches: Vector of prssnps indices that had no matching values in vcfsnpinfo
#' @export
#' @examples
#' \dontrun{
#' result <- findSNPs(prssnps, vcfsnpinfo)
#' ref_matches <- result$ref_matches
#' alt_matches <- result$alt_matches
#' no_matches <- result$no_matches
#' }
findSNPs <- function(prssnps, vcfsnpinfo) {
  # Use data.table for fast joins
  prs_dt <- data.table::as.data.table(prssnps)
  prs_dt[, prs_idx := .I]

  vcf_dt <- data.table::as.data.table(vcfsnpinfo)
  vcf_dt[, vcf_idx := .I]

  # REF matches - vectorized join on position and REF allele
  ref_joined <- vcf_dt[prs_dt,
                       on = .(position = chr_position, REF = effect_allele),
                       nomatch = NULL,
                       allow.cartesian = TRUE,
                       .(prs_idx = i.prs_idx, vcf_idx = x.vcf_idx)]

  # ALT matches - vectorized join on position and ALT allele
  alt_joined <- vcf_dt[prs_dt,
                       on = .(position = chr_position, ALT = effect_allele),
                       nomatch = NULL,
                       allow.cartesian = TRUE,
                       .(prs_idx = i.prs_idx, vcf_idx = x.vcf_idx)]

  # Convert to list format (pre-allocate lists)
  # Group by prs_idx and create vectors of [prs_idx, vcf_idx1, vcf_idx2, ...]
  if (nrow(ref_joined) > 0) {
    ref_split <- split(ref_joined$vcf_idx, ref_joined$prs_idx)
    ref_prs_indices <- as.integer(names(ref_split))
    matched_list <- vector("list", length(ref_split))
    for (i in seq_along(ref_split)) {
      matched_list[[i]] <- c(ref_prs_indices[i], ref_split[[i]])
    }
  } else {
    matched_list <- list()
  }

  if (nrow(alt_joined) > 0) {
    alt_split <- split(alt_joined$vcf_idx, alt_joined$prs_idx)
    alt_prs_indices <- as.integer(names(alt_split))
    alt_matched_list <- vector("list", length(alt_split))
    for (i in seq_along(alt_split)) {
      alt_matched_list[[i]] <- c(alt_prs_indices[i], alt_split[[i]])
    }
  } else {
    alt_matched_list <- list()
  }

  # Track indices that had at least one match
  matched_prs_indices <- unique(c(ref_joined$prs_idx, alt_joined$prs_idx))
  all_prs_indices <- seq_len(nrow(prssnps))
  no_match_indices <- setdiff(all_prs_indices, matched_prs_indices)

  return(list(ref_matches = matched_list, alt_matches = alt_matched_list, no_matches = no_match_indices))
}

#' Calculate PRS for ALT Allele Matches
#'
#' @param alt_matches List where each element is a vector containing the prssnps index followed by matching vcfsnpinfo indices (position and ALT match)
#' @param prssnps Data frame with PRS SNP information containing columns:
#'   - rsID: SNP identifier
#'   - chr_position: Chromosome position
#'   - effect_weight: Effect weight for PRS calculation
#'   - effect_allele: Effect allele
#' @param vcf_summary List returned from summarizevcf() function containing vcf_file, genome, snp_summary, and sample_ids
#' @param batch_size Integer; maximum number of SNPs to read from VCF at once. Default is 1000.
#' @param verbose Logical; if TRUE, display informational messages during execution. Default is FALSE.
#' @return A numeric vector of PRS values for all samples
#' @export
#' @examples
#' \dontrun{
#' vcf_file <- "data/axiom_acs_aus_nf_chr21.vcf.gz"
#' vcf_summary <- summarizevcf(vcf_file)
#' snp_matches <- findSNPs(prssnps, vcf_summary$snp_summary$snp_df)
#' prs_alt <- calculatePRSalt(snp_matches$alt_matches, prssnps, vcf_summary)
#' }
calculatePRSalt <- function(alt_matches, prssnps, vcf_summary, batch_size = 1000, verbose = FALSE) {
  # Validate inputs
  if (!is.list(alt_matches)) {
    stop("alt_matches must be a list")
  }

  if (length(alt_matches) == 0) {
    warning("alt_matches is empty, returning NULL")
    return(NULL)
  }

  if (!is.data.frame(prssnps) || !"effect_weight" %in% names(prssnps)) {
    stop("prssnps must be a data frame containing 'effect_weight' column")
  }

  if (!is.list(vcf_summary) || !all(c("vcf_file", "genome", "snp_summary") %in% names(vcf_summary))) {
    stop("vcf_summary must be a list returned from summarizevcf()")
  }

  if (verbose) {
    message(sprintf("Calculating PRS for %d ALT allele matches", length(alt_matches)))
  }

  # Collect all unique VCF indices needed
  all_vcf_indices <- unique(unlist(lapply(alt_matches, function(x) x[-1])))

  if (verbose) {
    message(sprintf("Loading dosages for %d unique VCF SNPs in batches of %d", length(all_vcf_indices), batch_size))
  }

  # Pre-load all dosages in batches
  dosage_list <- list()
  n_batches <- ceiling(length(all_vcf_indices) / batch_size)

  for (batch_num in seq_len(n_batches)) {
    start_idx <- (batch_num - 1) * batch_size + 1
    end_idx <- min(batch_num * batch_size, length(all_vcf_indices))
    batch_indices <- all_vcf_indices[start_idx:end_idx]

    if (verbose) {
      message(sprintf("  Loading batch %d/%d (%d SNPs)", batch_num, n_batches, length(batch_indices)))
    }

    batch_dosage <- getvcfsnp(vcf_summary, snp_indices = batch_indices, verbose = FALSE)

    # Store dosages indexed by VCF index
    if (is.matrix(batch_dosage)) {
      for (i in seq_along(batch_indices)) {
        dosage_list[[as.character(batch_indices[i])]] <- batch_dosage[i, ]
      }
    } else {
      # Single SNP returned as vector
      dosage_list[[as.character(batch_indices[1])]] <- batch_dosage
    }
  }

  # Initialize PRS vector
  sample_names <- names(dosage_list[[1]])
  n_samples <- length(sample_names)
  prs <- numeric(n_samples)
  names(prs) <- sample_names

  # Process each match using pre-loaded dosages
  for (match_idx in seq_along(alt_matches)) {
    match_vec <- alt_matches[[match_idx]]

    # First element is prssnps index
    prssnps_idx <- match_vec[1]

    # Remaining elements are vcfsnpinfo indices
    vcfsnpinfo_indices <- match_vec[-1]

    # Get effect weight from prssnps
    effect_weight <- prssnps$effect_weight[prssnps_idx]

    # Get dosages from pre-loaded data and sum if multiple
    if (length(vcfsnpinfo_indices) == 1) {
      dosage_sum <- dosage_list[[as.character(vcfsnpinfo_indices)]]
    } else {
      dosage_matrix <- do.call(rbind, lapply(vcfsnpinfo_indices, function(idx) dosage_list[[as.character(idx)]]))
      dosage_sum <- colSums(dosage_matrix)
    }

    # Add weighted dosage to PRS
    prs <- prs + dosage_sum * effect_weight
  }

  if (verbose) {
    message(sprintf("PRS calculation complete for %d samples using %d ALT matches", n_samples, length(alt_matches)))
  }

  return(prs)
}

#' Calculate PRS for REF Allele Matches
#'
#' @param ref_matches List where each element is a vector containing the prssnps index followed by matching vcfsnpinfo indices (position and REF match)
#' @param prssnps Data frame with PRS SNP information containing columns:
#'   - rsID: SNP identifier
#'   - chr_position: Chromosome position
#'   - effect_weight: Effect weight for PRS calculation
#'   - effect_allele: Effect allele
#' @param vcf_summary List returned from summarizevcf() function containing vcf_file, genome, snp_summary, and sample_ids
#' @param batch_size Integer; maximum number of SNPs to read from VCF at once. Default is 1000.
#' @param verbose Logical; if TRUE, display informational messages during execution. Default is FALSE.
#' @return A numeric vector of PRS values for all samples
#' @details For REF matches, the VCF dosage represents the non-effect (alternate) allele dosage.
#'   The effect allele dosage is calculated as (2 - dosage_sum), then multiplied by the effect weight.
#' @export
#' @examples
#' \dontrun{
#' vcf_file <- "data/axiom_acs_aus_nf_chr21.vcf.gz"
#' vcf_summary <- summarizevcf(vcf_file)
#' snp_matches <- findSNPs(prssnps, vcf_summary$snp_summary$snp_df)
#' prs_ref <- calculatePRSref(snp_matches$ref_matches, prssnps, vcf_summary)
#' }
calculatePRSref <- function(ref_matches, prssnps, vcf_summary, batch_size = 1000, verbose = FALSE) {
  # Validate inputs
  if (!is.list(ref_matches)) {
    stop("ref_matches must be a list")
  }

  if (length(ref_matches) == 0) {
    warning("ref_matches is empty, returning NULL")
    return(NULL)
  }

  if (!is.data.frame(prssnps) || !"effect_weight" %in% names(prssnps)) {
    stop("prssnps must be a data frame containing 'effect_weight' column")
  }

  if (!is.list(vcf_summary) || !all(c("vcf_file", "genome", "snp_summary") %in% names(vcf_summary))) {
    stop("vcf_summary must be a list returned from summarizevcf()")
  }

  if (verbose) {
    message(sprintf("Calculating PRS for %d REF allele matches", length(ref_matches)))
  }

  # Collect all unique VCF indices needed
  all_vcf_indices <- unique(unlist(lapply(ref_matches, function(x) x[-1])))

  if (verbose) {
    message(sprintf("Loading dosages for %d unique VCF SNPs in batches of %d", length(all_vcf_indices), batch_size))
  }

  # Pre-load all dosages in batches
  dosage_list <- list()
  n_batches <- ceiling(length(all_vcf_indices) / batch_size)

  for (batch_num in seq_len(n_batches)) {
    start_idx <- (batch_num - 1) * batch_size + 1
    end_idx <- min(batch_num * batch_size, length(all_vcf_indices))
    batch_indices <- all_vcf_indices[start_idx:end_idx]

    if (verbose) {
      message(sprintf("  Loading batch %d/%d (%d SNPs)", batch_num, n_batches, length(batch_indices)))
    }

    batch_dosage <- getvcfsnp(vcf_summary, snp_indices = batch_indices, verbose = FALSE)

    # Store dosages indexed by VCF index
    if (is.matrix(batch_dosage)) {
      for (i in seq_along(batch_indices)) {
        dosage_list[[as.character(batch_indices[i])]] <- batch_dosage[i, ]
      }
    } else {
      # Single SNP returned as vector
      dosage_list[[as.character(batch_indices[1])]] <- batch_dosage
    }
  }

  # Initialize PRS vector
  sample_names <- names(dosage_list[[1]])
  n_samples <- length(sample_names)
  prs <- numeric(n_samples)
  names(prs) <- sample_names

  # Process each match using pre-loaded dosages
  for (match_idx in seq_along(ref_matches)) {
    match_vec <- ref_matches[[match_idx]]

    # First element is prssnps index
    prssnps_idx <- match_vec[1]

    # Remaining elements are vcfsnpinfo indices
    vcfsnpinfo_indices <- match_vec[-1]

    # Get effect weight from prssnps
    effect_weight <- prssnps$effect_weight[prssnps_idx]

    # Get dosages from pre-loaded data and sum if multiple
    if (length(vcfsnpinfo_indices) == 1) {
      dosage_sum <- dosage_list[[as.character(vcfsnpinfo_indices)]]
    } else {
      dosage_matrix <- do.call(rbind, lapply(vcfsnpinfo_indices, function(idx) dosage_list[[as.character(idx)]]))
      dosage_sum <- colSums(dosage_matrix)
    }

    # For REF matches, the effect allele is the reference allele
    # VCF dosage represents alternate allele, so effect allele dosage = 2 - dosage_sum
    # Add weighted effect allele dosage to PRS
    prs <- prs + (2 - dosage_sum) * effect_weight
  }

  if (verbose) {
    message(sprintf("PRS calculation complete for %d samples using %d REF matches", n_samples, length(ref_matches)))
  }

  return(prs)
}

#' Process a Single Chromosome for Multiple PRS Models
#'
#' Internal helper function that processes one chromosome for multiple PRS models.
#' This function is designed to be called in parallel across chromosomes.
#'
#' @param chr Chromosome identifier
#' @param models Named list of PGS model data frames
#' @param model_names Character vector of model names
#' @param batch_size Integer; maximum number of SNPs to read from VCF at once
#' @param vcf_dir Character string specifying directory containing VCF summary files
#' @return A list containing results for all models for this chromosome
#' @keywords internal
processChromosomeForModels <- function(chr, models, model_names, batch_size, vcf_dir) {
  # Load VCF summary for this chromosome
  vcf_summary_file <- file.path(vcf_dir, sprintf("chromosome%sVCFsummary.rds", chr))

  if (!file.exists(vcf_summary_file)) {
    # Return empty results for all models
    empty_results <- lapply(model_names, function(mn) {
      list(
        model_name = mn,
        chr = chr,
        chr_prs = NULL,
        unmatched_count = 0,
        unmatched_rsIDs = character(0),
        status = "file_not_found"
      )
    })
    names(empty_results) <- model_names
    return(list(chr = chr, model_results = empty_results, status = "file_not_found"))
  }

  vcf_summary <- tryCatch({
    readRDS(vcf_summary_file)
  }, error = function(e) {
    return(NULL)
  })

  if (is.null(vcf_summary)) {
    empty_results <- lapply(model_names, function(mn) {
      list(
        model_name = mn,
        chr = chr,
        chr_prs = NULL,
        unmatched_count = 0,
        unmatched_rsIDs = character(0),
        status = "read_error"
      )
    })
    names(empty_results) <- model_names
    return(list(chr = chr, model_results = empty_results, status = "read_error"))
  }

  # Process each model for this chromosome
  model_results <- list()

  for (model_name in model_names) {
    pgsmodel <- models[[model_name]]

    # Filter model for current chromosome
    chr_pgsmodel <- pgsmodel[pgsmodel$chr_name == chr, ]

    if (nrow(chr_pgsmodel) == 0) {
      model_results[[model_name]] <- list(
        model_name = model_name,
        chr = chr,
        chr_prs = NULL,
        unmatched_count = 0,
        unmatched_rsIDs = character(0),
        status = "no_variants"
      )
      next
    }

    # Find matching SNPs
    snp_matches <- findSNPs(chr_pgsmodel, vcf_summary$snp_summary$snp_df)

    # Track unmatched SNPs
    unmatched_count <- length(snp_matches$no_matches)
    if (unmatched_count > 0) {
      unmatched_rsIDs <- chr_pgsmodel$rsID[snp_matches$no_matches]
    } else {
      unmatched_rsIDs <- character(0)
    }

    # Calculate PRS for ALT matches
    prs_alt <- NULL
    if (length(snp_matches$alt_matches) > 0) {
      prs_alt <- calculatePRSalt(snp_matches$alt_matches, chr_pgsmodel, vcf_summary, batch_size = batch_size, verbose = FALSE)
    }

    # Calculate PRS for REF matches
    prs_ref <- NULL
    if (length(snp_matches$ref_matches) > 0) {
      prs_ref <- calculatePRSref(snp_matches$ref_matches, chr_pgsmodel, vcf_summary, batch_size = batch_size, verbose = FALSE)
    }

    # Combine PRS from ALT and REF matches
    chr_prs <- NULL

    if (!is.null(prs_alt) && !is.null(prs_ref)) {
      if (!all(names(prs_alt) == names(prs_ref))) {
        prs_ref <- prs_ref[names(prs_alt)]
      }
      chr_prs <- prs_alt + prs_ref
    } else if (!is.null(prs_alt)) {
      chr_prs <- prs_alt
    } else if (!is.null(prs_ref)) {
      chr_prs <- prs_ref
    }

    model_results[[model_name]] <- list(
      model_name = model_name,
      chr = chr,
      chr_prs = chr_prs,
      unmatched_count = unmatched_count,
      unmatched_rsIDs = unmatched_rsIDs,
      ref_matches = length(snp_matches$ref_matches),
      alt_matches = length(snp_matches$alt_matches),
      status = "success"
    )
  }

  return(list(chr = chr, model_results = model_results, status = "success"))
}

#' Fit Multiple PRS Models Efficiently
#'
#' Calculates PRS for multiple models by reading each chromosome's VCF data only once.
#' This is much more efficient than calling calculatePRS() separately for each model.
#'
#' @param model_files Character vector of file paths to PGS model files
#' @param chromosomes Integer vector of chromosome numbers to process. Default is 1:22 (autosomes).
#' @param vcf_dir Character string specifying directory containing VCF summary files. Default is "." (current directory).
#' @param batch_size Integer; maximum number of SNPs to read from VCF at once. Default is 1000.
#' @param parallel Logical; if TRUE, process chromosomes in parallel using future. Default is FALSE.
#' @param n_cores Integer; number of cores to use for parallel processing. Default is NULL (use all available).
#' @param verbose Logical; if TRUE (default), display informational messages during execution
#' @return A list containing:
#'   - prs: A matrix where rows are samples and columns are models (named after model files)
#'   - unmatched_by_chr: A named list where each entry corresponds to a model and contains
#'       a named vector with counts of unmatched SNPs by chromosome
#'   - unmatched_rsIDs: A named list where each entry corresponds to a model and contains
#'       a named list with rsIDs of unmatched SNPs by chromosome
#'   - excluded_chr_counts: A named list where each entry corresponds to a model and contains
#'       a named vector with counts of SNPs not on processed chromosomes
#' @export
#' @examples
#' \dontrun{
#' # Fit multiple PRS models
#' model_files <- c("model1.txt", "model2.txt", "model3.txt")
#' results <- fitPRSmodels(model_files)
#' prs_matrix <- results$prs  # Rows = samples, columns = models
#'
#' # Access results for a specific model
#' model1_prs <- results$prs[, "model1.txt"]
#' model1_unmatched <- results$unmatched_by_chr[["model1.txt"]]
#'
#' # Fit models for specific chromosomes
#' results <- fitPRSmodels(model_files, chromosomes = c(1, 2, 3))
#'
#' # Fit models quietly
#' results <- fitPRSmodels(model_files, verbose = FALSE)
#'
#' # Fit models with custom batch size
#' results <- fitPRSmodels(model_files, batch_size = 500)
#'
#' # Fit models in parallel
#' results <- fitPRSmodels(model_files, parallel = TRUE)
#'
#' # Fit models in parallel with specific number of cores
#' results <- fitPRSmodels(model_files, parallel = TRUE, n_cores = 4)
#'
#' # Fit models using VCF summaries from a specific directory
#' results <- fitPRSmodels(model_files, vcf_dir = "vcf_summaries")
#' }
fitPRSmodels <- function(model_files, chromosomes = 1:22, vcf_dir = ".", batch_size = 1000, parallel = FALSE, n_cores = NULL, verbose = TRUE) {
  # Validate inputs
  if (!is.character(model_files) || length(model_files) == 0) {
    stop("model_files must be a non-empty character vector")
  }

  # Read all models
  if (verbose) {
    message(sprintf("Reading %d PGS model files...", length(model_files)))
  }

  models <- list()
  model_names <- character(length(model_files))

  for (i in seq_along(model_files)) {
    file_path <- model_files[i]

    # Use part of basename before first period or underscore as model name
    model_names[i] <- sub("[._].*", "", basename(file_path))

    if (verbose) {
      message(sprintf("  [%d/%d] Reading: %s", i, length(model_files), file_path))
    }

    models[[model_names[i]]] <- tryCatch({
      readPGSmodel(file_path, verbose = FALSE)
    }, error = function(e) {
      stop(sprintf("Error reading model file %s: %s", file_path, e$message))
    })
  }

  # Initialize output structures for each model
  excluded_chr_counts_list <- vector("list", length(models))
  names(excluded_chr_counts_list) <- model_names

  # Calculate excluded SNPs for each model
  for (model_name in model_names) {
    pgsmodel <- models[[model_name]]
    excluded_snps <- pgsmodel[!pgsmodel$chr_name %in% chromosomes, ]
    excluded_chr_counts_list[[model_name]] <- table(excluded_snps$chr_name)

    if (verbose && nrow(excluded_snps) > 0) {
      message(sprintf("Model '%s': %d variants not on processed chromosomes will be excluded",
                     model_name, nrow(excluded_snps)))
    }
  }

  # Process chromosomes (parallel or sequential)
  if (parallel) {
    # Check for required packages
    required_pkgs <- c("future", "future.apply", "parallelly")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
      stop(sprintf("Packages '%s' are required for parallel processing. Install with: install.packages(c('future', 'future.apply', 'parallelly'))",
                   paste(missing_pkgs, collapse = "', '")))
    }

    # Set up parallel backend
    if (is.null(n_cores)) {
      n_cores <- parallelly::availableCores()
    }

    # Use multicore on Linux/macOS (forking), multisession on Windows
    if (.Platform$OS.type == "unix") {
      parallel_type <- "multicore"
      if (verbose) {
        message(sprintf("Processing chromosomes in parallel using %d cores (multicore/forking)", n_cores))
      }
      old_plan <- future::plan()
      on.exit(future::plan(old_plan), add = TRUE)
      future::plan(future::multicore, workers = n_cores)
    } else {
      parallel_type <- "multisession"
      if (verbose) {
        message(sprintf("Processing chromosomes in parallel using %d cores (multisession)", n_cores))
      }
      old_plan <- future::plan()
      on.exit(future::plan(old_plan), add = TRUE)
      future::plan(future::multisession, workers = n_cores)
    }

    # Process chromosomes in parallel
    chr_results <- future.apply::future_lapply(chromosomes, function(chr) {
      processChromosomeForModels(chr, models, model_names, batch_size, vcf_dir)
    }, future.seed = TRUE)

  } else {
    # Sequential processing
    chr_results <- lapply(chromosomes, function(chr) {
      if (verbose) {
        message(sprintf("\n=== Processing chromosome %s ===", chr))
      }

      result <- processChromosomeForModels(chr, models, model_names, batch_size, vcf_dir)

      if (verbose) {
        if (result$status == "file_not_found") {
          message(sprintf("VCF summary file not found for chromosome %s - skipping", chr))
        } else if (result$status == "read_error") {
          message(sprintf("Error reading VCF summary for chromosome %s - skipping", chr))
        } else {
          for (model_name in model_names) {
            mr <- result$model_results[[model_name]]
            if (mr$status == "no_variants") {
              message(sprintf("  Model '%s': No variants for chromosome %s - skipping", model_name, chr))
            } else if (mr$status == "success") {
              message(sprintf("  Model '%s': REF matches: %d, ALT matches: %d, No matches: %d",
                             model_name, mr$ref_matches, mr$alt_matches, mr$unmatched_count))
            }
          }
        }
      }

      result
    })
  }

  # Combine results from all chromosomes
  total_prs_list <- vector("list", length(models))
  names(total_prs_list) <- model_names

  unmatched_counts_list <- vector("list", length(models))
  names(unmatched_counts_list) <- model_names

  unmatched_rsIDs_list <- vector("list", length(models))
  names(unmatched_rsIDs_list) <- model_names

  # Initialize each model's tracking structures
  for (model_name in model_names) {
    total_prs_list[[model_name]] <- NULL
    unmatched_counts_list[[model_name]] <- integer()
    unmatched_rsIDs_list[[model_name]] <- list()
  }

  sample_names <- NULL

  for (chr_result in chr_results) {
    chr_key <- as.character(chr_result$chr)

    for (model_name in model_names) {
      mr <- chr_result$model_results[[model_name]]

      # Track unmatched
      unmatched_counts_list[[model_name]][chr_key] <- mr$unmatched_count
      unmatched_rsIDs_list[[model_name]][[chr_key]] <- mr$unmatched_rsIDs

      if (!is.null(mr$chr_prs)) {
        # Set sample names from first result with data
        if (is.null(sample_names)) {
          sample_names <- names(mr$chr_prs)
        }

        if (is.null(total_prs_list[[model_name]])) {
          total_prs_list[[model_name]] <- mr$chr_prs
        } else {
          # Ensure sample names match before adding
          if (!all(names(mr$chr_prs) == sample_names)) {
            mr$chr_prs <- mr$chr_prs[sample_names]
          }
          total_prs_list[[model_name]] <- total_prs_list[[model_name]] + mr$chr_prs
        }
      }
    }
  }

  # Combine PRS vectors into a matrix
  if (is.null(sample_names)) {
    warning("No PRS calculated - no matching variants found in any chromosome for any model")
    return(list(
      prs = NULL,
      unmatched_by_chr = unmatched_counts_list,
      unmatched_rsIDs = unmatched_rsIDs_list,
      excluded_chr_counts = excluded_chr_counts_list
    ))
  }

  # Create PRS matrix
  prs_matrix <- matrix(NA, nrow = length(sample_names), ncol = length(model_names))
  rownames(prs_matrix) <- sample_names
  colnames(prs_matrix) <- model_names

  for (i in seq_along(model_names)) {
    model_name <- model_names[i]
    if (!is.null(total_prs_list[[model_name]])) {
      # Ensure correct order
      prs_matrix[, i] <- total_prs_list[[model_name]][sample_names]
    } else {
      warning(sprintf("No PRS calculated for model '%s' - no matching variants found", model_name))
    }
  }

  # Final summary
  if (verbose) {
    message(sprintf("\n=== PRS Calculation Complete ==="))
    message(sprintf("Total samples: %d", length(sample_names)))
    message(sprintf("Total models: %d", length(model_names)))

    for (model_name in model_names) {
      if (!is.null(total_prs_list[[model_name]])) {
        prs_vec <- total_prs_list[[model_name]]
        message(sprintf("\nModel '%s':", model_name))
        message(sprintf("  PRS: mean=%.4f, sd=%.4f, min=%.4f, max=%.4f",
                       mean(prs_vec), sd(prs_vec), min(prs_vec), max(prs_vec)))

        total_unmatched <- sum(unmatched_counts_list[[model_name]])
        if (total_unmatched > 0) {
          message(sprintf("  Unmatched SNPs: %d", total_unmatched))
        }

        if (length(excluded_chr_counts_list[[model_name]]) > 0) {
          message(sprintf("  Excluded SNPs: %d", sum(excluded_chr_counts_list[[model_name]])))
        }
      }
    }
  }

  return(list(
    prs = prs_matrix,
    unmatched_by_chr = unmatched_counts_list,
    unmatched_rsIDs = unmatched_rsIDs_list,
    excluded_chr_counts = excluded_chr_counts_list
  ))
}
