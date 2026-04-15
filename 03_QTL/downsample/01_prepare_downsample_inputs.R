#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

read_table_auto <- function(path) {
  if (!grepl("\\.gz$", path)) {
    return(fread(path))
  }

  direct <- tryCatch(
    fread(path),
    error = function(e) e
  )
  if (!inherits(direct, "error")) {
    return(direct)
  }

  stream_cmd <- sprintf("gzip -dc %s", shQuote(path))
  fread(cmd = stream_cmd)
}

normalize_chr_values <- function(chr_values, genotype_has_chr_prefix) {
  chr_values <- as.character(chr_values)
  has_chr <- startsWith(chr_values, "chr")
  if (genotype_has_chr_prefix) {
    chr_values[!has_chr] <- paste0("chr", chr_values[!has_chr])
  } else {
    chr_values[has_chr] <- sub("^chr", "", chr_values[has_chr])
  }
  chr_values
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  set.seed(2026)

  base_input <- "../01_pre"
  out_dir <- "downsample/pre_input"
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  stage_map <- c(stage1 = "Stage1", stage2 = "Stage2", stage3 = "Stage3")
  sample_lists <- lapply(stage_map, function(stage_label) {
    readLines(file.path(base_input, sprintf("sample.%s.final", stage_label)))
  })

  target_n <- min(lengths(sample_lists))
  message(sprintf("[INFO] Downsampling all stages to n = %d", target_n))

  summary_rows <- list()
  downsample_samples <- list()

  for (stage_key in names(stage_map)) {
    stage_label <- stage_map[[stage_key]]
    samples <- sample_lists[[stage_key]]
    if (length(samples) < target_n) {
      stop(sprintf("Stage %s has only %d samples, below target_n=%d", stage_label, length(samples), target_n))
    }

    selected <- sort(sample(samples, target_n, replace = FALSE))
    downsample_samples[[stage_label]] <- selected
    writeLines(selected, file.path(out_dir, sprintf("sample.%s.downsample", stage_label)))

    cov <- read_table_auto(file.path(base_input, sprintf("covariates_for_qtl.%s.txt", stage_label)))
    keep_cov_cols <- c(names(cov)[1], selected)
    fwrite(cov[, ..keep_cov_cols], file.path(out_dir, sprintf("covariates_for_qtl.%s.txt", stage_label)),
           sep = "\t", quote = FALSE)

    pheno <- read_table_auto(file.path(base_input, sprintf("phenotype.%s.bed.gz", stage_label)))
    # Downsampled genotypes are regenerated from VCF and use chr-prefixed chromosomes.
    pheno[, `#Chr` := normalize_chr_values(`#Chr`, TRUE)]
    keep_pheno_cols <- c("#Chr", "start", "end", "pid", selected)
    fwrite(pheno[, ..keep_pheno_cols], file.path(out_dir, sprintf("phenotype.%s.bed", stage_label)),
           sep = "\t", quote = FALSE)

    summary_rows[[stage_label]] <- data.table(
      stage = stage_label,
      original_n = length(samples),
      downsample_n = length(selected)
    )
  }

  write_json(downsample_samples, file.path(out_dir, "downsample_samples.json"), auto_unbox = TRUE, pretty = TRUE)
  summary_dt <- rbindlist(summary_rows)
  fwrite(summary_dt, file.path(out_dir, "downsample_summary.tsv"), sep = "\t", quote = FALSE)
  message("[INFO] Downsample input preparation complete")
}

if (sys.nframe() == 0) {
  main()
}
