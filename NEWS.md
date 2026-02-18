# PRSfromPGS 0.99.0

## New Features

* Initial Bioconductor submission.
* `prepareVCFfiles()` — batch-indexes and summarizes bgzipped VCF files,
  saving per-chromosome RDS summaries for fast downstream access.
* `readPGSmodel()` — parses PGS Catalog-format scoring files with automatic
  column normalisation and composite variant ID creation.
* `calculatePRS()` — computes per-sample Polygenic Risk Scores across all
  autosomes, with separate REF/ALT dosage paths and optional parallel
  processing via `future`/`future.apply`.
* `fitPRSmodels()` — efficiently scores multiple PGS models by reading each
  chromosome's VCF data once.
* Optional parallelism via `BiocParallel` (VCF indexing) and
  `future`/`future.apply` (chromosome-level scoring).
