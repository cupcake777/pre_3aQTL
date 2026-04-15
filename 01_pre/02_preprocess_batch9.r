#!/usr/bin/env Rscript
# Batch9 preprocessing script
# - Imputation only (without batch correction)
# - log2 transformation
# - Output for comparison with other batches

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
  library(readr)
  library(missForest)
  library(flashier)
  library(irlba)
  library(ebnm)
  library(ggplot2)
  library(ggrepel)
})

check_system_tool <- function(tool_name) {
  result <- system(paste("which", tool_name), intern = TRUE, ignore.stderr = TRUE)
  return(length(result) > 0 && nzchar(result[1]))
}

cap_pc_count <- function(n_samples, n_features, requested_pcs) {
  max_pcs <- min(n_samples, n_features, requested_pcs)
  as.integer(max(max_pcs, 0))
}

main <- function() {
  # Check for required system tools
  missing_tools <- c()
  if (!check_system_tool("bgzip")) missing_tools <- c(missing_tools, "bgzip")
  if (!check_system_tool("tabix")) missing_tools <- c(missing_tools, "tabix")

  if (length(missing_tools) > 0) {
    warning("Missing system tools: ", paste(missing_tools, collapse = ", "))
  }

  set.seed(2026)

  # File paths
  apa_file <- "../../raw_data/batch9.apa.count"
  out_dir <- "./"

  if (!file.exists(apa_file)) {
    stop("Missing file: ", apa_file)
  }

  message("Loading batch9 data...")
  data <- fread(apa_file, check.names = FALSE)

  # Parse Gene column
  if (grepl("\\|", data$Gene[1])) {
    data[, c("transcript", "gene", "chr", "strand") := tstrsplit(Gene, "|", fixed = TRUE)]
  } else {
    stop("Gene column does not contain '|'")
  }

  data[, phenotype_id := paste(transcript, gene, Loci, strand, sep = "|")]
  saving_names <- data$phenotype_id

  # Extract phenotype matrix
  info_cols_indices <- 1:4
  info_data <- cbind(data[, ..info_cols_indices], data[, .(phenotype_id)])
  fwrite(info_data, file.path(out_dir, "batch9.apa.info"), sep = '\t', quote = FALSE)

  cols_to_exclude <- c(names(data)[1:4], "transcript", "gene", "chr", "strand", "phenotype_id")
  expr_data <- data[, !cols_to_exclude, with = FALSE]
  expr_data <- as.data.frame(expr_data)
  rownames(expr_data) <- saving_names

  # Sample QC: remove samples with >80% missing
  sample_na_rate <- colMeans(is.na(expr_data))
  keep_samples <- sample_na_rate <= 0.8
  expr_data_sample_qc <- expr_data[, keep_samples]
  message(sprintf("Samples after QC: %d / %d", sum(keep_samples), ncol(expr_data)))

  # Phenotype QC: remove rows with >50% missing or zero variance
  row_na_rate <- rowMeans(is.na(expr_data_sample_qc))
  row_sd <- apply(expr_data_sample_qc, 1, sd, na.rm = TRUE)
  keep_rows <- (row_na_rate <= 0.5) & (row_sd >= 0.001) & (!is.na(row_sd))
  expr_data_qc <- expr_data_sample_qc[keep_rows, ]
  message(sprintf("Phenotypes after QC: %d / %d", sum(keep_rows), nrow(expr_data_sample_qc)))

  # Imputation using softImpute + flash
  message("Running imputation...")
  lambda <- if (nrow(expr_data_qc) < 1500) 5 else 30
  rank <- 50
  X_mis_C <- as(as.matrix(expr_data_qc), "Incomplete")
  fit1 <- softImpute::softImpute(X_mis_C, rank = rank, lambda = lambda, type = "svd")
  pheno_soft <- softImpute::complete(as.matrix(expr_data_qc), fit1)

  # Flash for better imputation
  min_sd <- min(apply(pheno_soft, 1, sd))
  pca_res <- irlba::irlba(pheno_soft, nv = 10)
  pca_res_list <- list(d = pca_res$d, u = pca_res$u, v = pca_res$v)

  fl_pca <- flash_init(as.matrix(expr_data_qc), S = min_sd, var_type = 1) |>
    flash_factors_init(pca_res_list, ebnm_fn = ebnm::ebnm_point_laplace) |>
    flash_backfit(maxiter = 300)

  expr_data_imputed <- ifelse(is.na(as.matrix(expr_data_qc)), fitted(fl_pca), as.matrix(expr_data_qc))

  # Fill any remaining NAs with row means
  if (any(is.na(expr_data_imputed))) {
    message("Warning: Post-imputation NAs found. Filling with row means.")
    row_means <- rowMeans(expr_data_imputed, na.rm = TRUE)
    for (i in seq_len(nrow(expr_data_imputed))) {
      expr_data_imputed[i, is.na(expr_data_imputed[i, ])] <- row_means[i]
    }
  }

  # Save raw imputed data
  fwrite(as.data.frame(expr_data_imputed), file.path(out_dir, "batch9.after.impute.txt"),
         sep = '\t', quote = FALSE, row.names = TRUE)

  # Log2 transformation (no batch correction for batch9)
  message("Running log2 transformation...")
  epsilon <- 0.001
  dat_bounded <- pmin(pmax(expr_data_imputed, epsilon), 1 - epsilon)
  expr_data_scaled <- log2(dat_bounded / (1 - dat_bounded))

  # Save final output
  fwrite(as.data.frame(expr_data_scaled), file.path(out_dir, "batch9.after.combat.txt"),
         sep = '\t', quote = FALSE, row.names = TRUE)

  message("Batch9 preprocessing complete!")
  message("Output files:")
  message("  - batch9.after.impute.txt (raw imputed)")
  message("  - batch9.after.combat.txt (log2 transformed)")
  message("  - batch9.apa.info (phenotype info)")
}

main()