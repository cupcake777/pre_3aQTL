#!/usr/bin/env Rscript
# Check for required system tools
check_system_tool <- function(tool_name) {
  result <- system(paste("which", tool_name), intern = TRUE, ignore.stderr = TRUE)
  return(length(result) > 0 && nzchar(result[1]))
}

get_psam_iid <- function(psam_data) {
  iid_candidates <- c("IID", "#IID", "iid")
  iid_col <- iid_candidates[iid_candidates %in% colnames(psam_data)][1]
  if (!is.na(iid_col) && nzchar(iid_col)) {
    return(psam_data[[iid_col]])
  }
  if (ncol(psam_data) >= 2) {
    return(psam_data[[2]])
  }
  stop("Unable to determine IID column from .psam file")
}

read_eigenvec_pcs <- function(eigenvec_file, max_pcs = 5) {
  eigenvec_data <- data.table::fread(eigenvec_file, data.table = FALSE)
  iid_candidates <- c("IID", "#IID", "iid")
  iid_col <- iid_candidates[iid_candidates %in% colnames(eigenvec_data)][1]

  if (is.na(iid_col) || !nzchar(iid_col)) {
    if (ncol(eigenvec_data) < 2) {
      stop(sprintf("Unable to determine sample ID column from eigenvec file: %s", eigenvec_file))
    }
    iid_col <- colnames(eigenvec_data)[1]
  }

  id_cols <- unique(c("FID", "#FID", "fid", iid_candidates, iid_col))
  pc_cols <- setdiff(colnames(eigenvec_data), id_cols)
  if (length(pc_cols) == 0) {
    stop(sprintf("No PC columns found in eigenvec file: %s", eigenvec_file))
  }

  pc_cols <- pc_cols[seq_len(min(max_pcs, length(pc_cols)))]
  geno_pcs <- eigenvec_data[, pc_cols, drop = FALSE]
  rownames(geno_pcs) <- eigenvec_data[[iid_col]]
  colnames(geno_pcs) <- paste0("genoPC", seq_len(ncol(geno_pcs)))
  geno_pcs
}

align_covariate_samples <- function(sample_ids, known_covs, geno_pcs, expr_pcs) {
  common_samples <- sample_ids
  if (nrow(known_covs) > 0) {
    common_samples <- intersect(common_samples, rownames(known_covs))
  }
  if (nrow(geno_pcs) > 0) {
    common_samples <- intersect(common_samples, rownames(geno_pcs))
  }
  common_samples <- intersect(common_samples, rownames(expr_pcs))

  list(
    sample_ids = common_samples,
    known_covs = known_covs[common_samples, , drop = FALSE],
    geno_pcs = geno_pcs[common_samples, , drop = FALSE],
    expr_pcs = expr_pcs[common_samples, , drop = FALSE]
  )
}

detect_chr_prefix <- function(pvar_file) {
  con <- file(pvar_file, open = "r")
  on.exit(close(con), add = TRUE)

  repeat {
    pvar_lines <- readLines(con, warn = FALSE, n = 1000)
    if (length(pvar_lines) == 0) {
      stop("Unable to detect chromosome format: no variant records found in .pvar file")
    }
    pvar_lines <- pvar_lines[nzchar(pvar_lines)]
    data_lines <- pvar_lines[!grepl("^#", pvar_lines)]
    if (length(data_lines) > 0) {
      first_chr <- strsplit(data_lines[1], "\t", fixed = TRUE)[[1]][1]
      return(grepl("^chr", first_chr))
    }
  }
}

cap_pc_count <- function(n_samples, n_features, requested_pcs) {
  max_pcs <- min(n_samples, n_features, requested_pcs)
  as.integer(max(max_pcs, 0))
}

run_command <- function(command, failure_message, success_message = NULL) {
  exit_code <- system(command)
  if (exit_code != 0) {
    stop(sprintf("%s (exit code %d)", failure_message, exit_code))
  }
  if (!is.null(success_message)) {
    message(success_message)
  }
  invisible(exit_code)
}

main <- function() {
  # Check for required system tools and warn if missing
  missing_tools <- c()
  if (!check_system_tool("bgzip")) missing_tools <- c(missing_tools, "bgzip")
  if (!check_system_tool("tabix")) missing_tools <- c(missing_tools, "tabix")
  if (!check_system_tool("plink2")) missing_tools <- c(missing_tools, "plink2")

  if (length(missing_tools) > 0) {
    warning("Missing system tools: ", paste(missing_tools, collapse = ", "))
    warning("Some parts of the pipeline may fail without these tools.")
  }

  suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(tibble)
    library(readr)
    library(missForest)
    library(flashier)
    library(PCAForQTL)
    library(sva)
    library(softImpute)
    library(irlba)
    library(ebnm)
    library(ggplot2)
    library(ggrepel)
  })

  # Set random seed for reproducibility
  set.seed(2026)

  # Check for required input files
  apa_file  <- "../../raw_data/all.apa.count"
  meta_file <- "../../raw_data/CHB.all.meta"
  geno_prefix <- "../../raw_data/genotype"
  out_dir   <- "./"

  # Check for plink2 genotype files (.pgen, .psam, .pvar)
  geno_files <- c(
    paste0(geno_prefix, ".pgen"),
    paste0(geno_prefix, ".psam"),
    paste0(geno_prefix, ".pvar")
  )
  input_files <- c(apa_file, meta_file, geno_files)
  missing_files <- input_files[!file.exists(input_files)]

  if (length(missing_files) > 0) {
    stop("Missing required input files: ", paste(missing_files, collapse = ", "))
  }

  stage_col   <- c(Stage1="#56B4E9", Stage2="#FFB482", Stage3="#8DE5A1")
  batch_col   <- c(batch1="#9FB1A3", batch2="#AFC3B3", batch3="#C0D7C5", batch4="#D3ECD9",
                   batch5="#E3EFFC", batch6="#BDD0D6", batch7="#B4C8CD", batch8="#BBCDD2",
                   batch9="#FFB482")

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  data <- fread(apa_file, check.names=FALSE)

  if (grepl("\\|", data$Gene[1])) {
    data[, c("transcript", "gene", "chr", "strand") := tstrsplit(Gene, "|", fixed=TRUE)]
  } else {
    stop("Gene column does not contain '|', cannot construct phenotype_id safely.")
  }
  data[, phenotype_id := paste(transcript, gene, Loci, strand, sep="|")]
  saving_names <- data$phenotype_id

  info_cols_indices <- 1:4
  info_data <- cbind(data[, ..info_cols_indices], data[, .(phenotype_id)])
  fwrite(info_data, file.path(out_dir, "proximal.apa.info"), sep='\t', quote=FALSE)
  cols_to_exclude <- c(names(data)[1:4], "transcript", "gene", "chr", "strand", "phenotype_id")
  expr_data <- data[, !cols_to_exclude, with=FALSE]
  expr_data <- as.data.frame(expr_data)
  rownames(expr_data) <- saving_names

  meta <- fread(meta_file, data.table=FALSE)
  rownames(meta) <- meta[,1]
  meta <- meta[,-1]

  common_samples <- intersect(colnames(expr_data), rownames(meta))
  expr_data <- expr_data[, common_samples]
  meta <- meta[common_samples, , drop=FALSE]

  sample_na_rate <- colMeans(is.na(expr_data))
  keep_samples <- sample_na_rate <= 0.8
  expr_data_sample_qc <- expr_data[, keep_samples]
  meta_clean <- meta[keep_samples, , drop=FALSE]

  row_na_rate <- rowMeans(is.na(expr_data_sample_qc))
  row_sd <- apply(expr_data_sample_qc, 1, sd, na.rm = TRUE)
  keep_rows <- (row_na_rate <= 0.5) & (row_sd >= 0.001) & (!is.na(row_sd))
  expr_data_qc <- expr_data_sample_qc[keep_rows, ]

  message("Imputation running...")
  lambda <- if (nrow(expr_data_qc) < 1500) 5 else 30
  rank <- 50
  X_mis_C <- as(as.matrix(expr_data_qc), "Incomplete")
  fit1 <- softImpute::softImpute(X_mis_C, rank = rank, lambda = lambda, type = "svd")
  pheno_soft <- softImpute::complete(as.matrix(expr_data_qc), fit1)

  min_sd <- min(apply(pheno_soft, 1, sd))
  pca_res <- irlba::irlba(pheno_soft, nv = 10)
  pca_res_list <- list(d = pca_res$d, u = pca_res$u, v = pca_res$v)

  fl_pca <- flash_init(as.matrix(expr_data_qc), S = min_sd, var_type = 1) |>
    flash_factors_init(pca_res_list, ebnm_fn = ebnm::ebnm_point_laplace) |>
    flash_backfit(maxiter = 300)

  expr_data_imputed <- ifelse(is.na(as.matrix(expr_data_qc)), fitted(fl_pca), as.matrix(expr_data_qc))

  if (any(is.na(expr_data_imputed))) {
    message("Warning: Post-imputation NAs found. Filling with row means.")
    row_means <- rowMeans(expr_data_imputed, na.rm = TRUE)
    for (i in seq_len(nrow(expr_data_imputed))) {
      expr_data_imputed[i, is.na(expr_data_imputed[i, ])] <- row_means[i]
    }
  }
  run_mahalanobis_outlier <- function(expr_matrix, metadata, n_pcs = 5, p_value_cutoff = 0.001) {
    vars <- apply(expr_matrix, 1, var)
    expr_for_pca <- expr_matrix[vars > 1e-10, , drop = FALSE]

    if (nrow(expr_for_pca) < 2 || ncol(expr_for_pca) < 3) {
      warning("Skipping Mahalanobis outlier detection: insufficient matrix dimensions after variance filter.")
      pca_data <- data.frame(
        sample = rownames(metadata),
        PC1 = 0,
        PC2 = 0,
        M_Dist = NA_real_,
        P_Value = NA_real_,
        LifeStage = metadata$LifeStage,
        IsOutlier = FALSE
      )
      return(list(data = pca_data, plot = ggplot(pca_data, aes(x = PC1, y = PC2, color = LifeStage)) +
                    geom_point(alpha = 0.6) + scale_color_manual(values = stage_col) + theme_bw()))
    }

    pca_obj <- prcomp(t(expr_for_pca), center = TRUE, scale. = TRUE)
    available_pcs <- cap_pc_count(nrow(pca_obj$x), ncol(pca_obj$x), n_pcs)
    if (available_pcs < 2) {
      warning("Skipping Mahalanobis outlier detection: fewer than 2 usable PCs available.")
      pca_data <- data.frame(
        sample = rownames(metadata),
        PC1 = 0,
        PC2 = 0,
        M_Dist = NA_real_,
        P_Value = NA_real_,
        LifeStage = metadata$LifeStage,
        IsOutlier = FALSE
      )
      return(list(data = pca_data, plot = ggplot(pca_data, aes(x = PC1, y = PC2, color = LifeStage)) +
                    geom_point(alpha = 0.6) + scale_color_manual(values = stage_col) + theme_bw()))
    }

    message(sprintf("Running Outlier Detection with Mahalanobis distance (PCs=1:%d)...", available_pcs))
    pcs <- pca_obj$x[, seq_len(available_pcs), drop = FALSE]

    pca_data_list <- list()
    for (stage in unique(metadata$LifeStage)) {
      stage_samples <- rownames(metadata[metadata$LifeStage == stage, , drop = FALSE])
      stage_samples <- intersect(stage_samples, rownames(pcs))

      if (length(stage_samples) <= available_pcs) {
        warning(sprintf("Sample size for stage %s is too small for df=%d. Skipping.", stage, available_pcs))
        next
      }

      group_pcs <- pcs[stage_samples, , drop = FALSE]
      group_center <- colMeans(group_pcs)
      group_cov <- cov(group_pcs)
      if (anyNA(group_cov) || qr(group_cov)$rank < ncol(group_pcs)) {
        warning(sprintf("Covariance matrix is singular for stage %s. Skipping outlier detection for this stage.", stage))
        next
      }

      m_dist <- mahalanobis(group_pcs, center = group_center, cov = group_cov)
      p_vals <- 1 - pchisq(m_dist, df = available_pcs)

      stage_df <- data.frame(
        sample = stage_samples,
        PC1 = group_pcs[, 1],
        PC2 = group_pcs[, 2],
        M_Dist = m_dist,
        P_Value = p_vals,
        LifeStage = stage,
        IsOutlier = p_vals < p_value_cutoff
      )
      pca_data_list[[stage]] <- stage_df
    }

    if (length(pca_data_list) == 0) {
      pca_data <- data.frame(
        sample = rownames(metadata),
        PC1 = 0,
        PC2 = 0,
        M_Dist = NA_real_,
        P_Value = NA_real_,
        LifeStage = metadata$LifeStage,
        IsOutlier = FALSE
      )
    } else {
      pca_data <- do.call(rbind, pca_data_list)
    }

    p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = LifeStage)) +
      geom_point(aes(size = M_Dist, alpha = IsOutlier), shape = 16) +
      scale_color_manual(values = stage_col) +
      scale_alpha_manual(values = c("FALSE" = 0.4, "TRUE" = 1), guide = "none") +
      geom_text_repel(data = subset(pca_data, IsOutlier),
                      aes(label = sample), color = "red", size = 3, fontface = "bold") +
      labs(title = "Outlier Detection: Mahalanobis Distance in PC Space",
           subtitle = sprintf("Red labels: p < %s (Chi-square dist, df=%d)", p_value_cutoff, available_pcs),
           x = "PC1", y = "PC2", size = "Mahalanobis Dist") +
      theme_bw() + theme(legend.position = "right", plot.title = element_text(face = "bold"))

    list(data = pca_data, plot = p)
  }

  outlier_res <- run_mahalanobis_outlier(expr_data_imputed, meta_clean, n_pcs = 5)
  ggsave(file.path(out_dir, "outlier_detect_mahalanobis.png"), outlier_res$plot, width = 8, height = 6)

  outliers <- outlier_res$data$sample[outlier_res$data$IsOutlier]
  if (length(outliers) > 0) {
    message(sprintf("Removing %d outlier samples based on Mahalanobis distance.", length(outliers)))
    samples_keep <- setdiff(colnames(expr_data_imputed), outliers)
    expr_data_imputed <- expr_data_imputed[, samples_keep]
    meta_clean <- meta_clean[samples_keep, ]
  }

  fwrite(as.data.frame(expr_data_imputed), file.path(out_dir, "after.impute.all.raw"), sep='\t', quote=FALSE, row.names=TRUE)

  apply_int <- function(expr_matrix) {
    apply(as.matrix(expr_matrix), 1, function(x) {
      n_val <- sum(!is.na(x))
      qnorm(rank(x, na.last = "keep", ties.method = "average") / (n_val + 1))
    }) %>% t()
  }

  epsilon <- 0.001
  dat_bounded <- pmin(pmax(expr_data_imputed, epsilon), 1 - epsilon)
  expr_data_scaled <- log2(dat_bounded / (1 - dat_bounded))

  if ("batch" %in% colnames(meta_clean)) {
    batch_vec <- as.factor(meta_clean$batch)
    mod <- model.matrix(~ LifeStage, data = meta_clean)
    ref_batch <- "batch9"
    message("Running ComBat with ", ref_batch, " as reference batch...")
    expr_clean <- ComBat(dat = as.matrix(expr_data_scaled), batch = batch_vec, mod = mod,
                         ref.batch = if (ref_batch %in% levels(batch_vec)) ref_batch else NULL)
  } else {
    expr_clean <- as.matrix(expr_data_scaled)
  }
  fwrite(as.data.frame(expr_clean), file.path(out_dir, "after.combat.txt"), sep='\t', quote=FALSE, row.names=TRUE)

  message("Performing global genotype alignment...")
  psam_file <- paste0(geno_prefix, ".psam")
  psam_data <- fread(psam_file, data.table = FALSE)
  geno_samples <- get_psam_iid(psam_data)
  common_samples_global <- intersect(colnames(expr_clean), geno_samples)
  message(sprintf("Global alignment: %d expression samples -> %d samples with genotype",
                  ncol(expr_clean), length(common_samples_global)))

  if (length(common_samples_global) < 10) {
    stop("Insufficient overlapping samples between expression and genotype data!")
  }

  expr_clean <- expr_clean[, common_samples_global, drop = FALSE]
  meta_clean <- meta_clean[common_samples_global, , drop = FALSE]

  stages <- unique(meta_clean$LifeStage)
  message(sprintf("After genotype alignment, processing %d stages: %s",
                  length(stages), paste(stages, collapse = ", ")))

  for (stage in stages) {
    message(sprintf("\n========== Processing %s ==========", stage))

    stage_samples <- rownames(meta_clean)[meta_clean$LifeStage == stage]
    if (length(stage_samples) < 5) {
      message(sprintf("[%s] Skipping: only %d samples available (< 5 minimum required)", stage, length(stage_samples)))
      next
    }

    sub_expr_clean <- expr_clean[, stage_samples, drop = FALSE]
    sub_meta <- meta_clean[stage_samples, , drop = FALSE]

    sub_expr_norm <- apply_int(sub_expr_clean)
    sub_expr_norm <- sub_expr_norm[complete.cases(sub_expr_norm), , drop = FALSE]
    out_filename <- sprintf("pre.phenotype.combat.int.%s.txt", stage)
    fwrite(as.data.frame(sub_expr_norm), file.path(out_dir, out_filename), sep = '\t', quote = FALSE, row.names = TRUE)

    stage_sample_file <- file.path(out_dir, sprintf("samples.%s.txt", tolower(stage)))
    writeLines(stage_samples, stage_sample_file)

    pgen_prefix <- file.path(out_dir, paste0("genotype.", tolower(stage)))
    plink_cmd <- sprintf("plink2 --pfile %s --keep %s --make-pgen --out %s",
                         geno_prefix, stage_sample_file, pgen_prefix)
    run_command(plink_cmd, sprintf("[%s] plink2 genotype extraction failed", stage))

    stage_pca_prefix <- file.path(out_dir, paste0("genotype_pca.", tolower(stage)))
    pca_cmd <- sprintf("plink2 --pfile %s --pca 20 --output-chr chrM --out %s",
                       pgen_prefix, stage_pca_prefix)
    run_command(pca_cmd, sprintf("[%s] plink2 PCA failed", stage))

    stage_eigenvec <- paste0(stage_pca_prefix, ".eigenvec")
    if (!file.exists(stage_eigenvec)) {
      stop(sprintf("[%s] Eigenvec file not found: %s", stage, stage_eigenvec))
    }

    vars_norm <- apply(sub_expr_norm, 1, var)
    sub_expr_for_pca <- sub_expr_norm[vars_norm > 1e-10, , drop = FALSE]

    if (ncol(sub_expr_for_pca) >= 2 && nrow(sub_expr_for_pca) >= 2) {
      sub_pca_final <- prcomp(t(sub_expr_for_pca), center = TRUE, scale. = TRUE)
      max_available_pcs <- ncol(sub_pca_final$x)
      if (max_available_pcs < 1) {
        warning(sprintf("[%s] No usable expression PCs available. Skipping.", stage))
        next
      }

      K_BE <- suppressWarnings(as.integer(PCAForQTL::runBE(t(sub_expr_for_pca), B = 20)$numOfPCsChosen))
      if (is.na(K_BE) || K_BE < 1) {
        K_BE <- 10L
      }
      K_final <- min(max(K_BE, 10L), 50L, max_available_pcs)

      remove_col <- c("batch", "sample_id", "LifeStage")
      known_covs <- sub_meta %>%
        mutate(across(where(is.character), as.factor)) %>%
        mutate(across(where(is.factor), as.numeric)) %>%
        select(where(is.numeric)) %>%
        select(-any_of(remove_col))

      pcs_all <- as.data.frame(sub_pca_final$x[, seq_len(K_final), drop = FALSE])
      geno_sub <- read_eigenvec_pcs(stage_eigenvec, max_pcs = 5)

      if (ncol(known_covs) == 0) {
        known_covs <- data.frame(row.names = rownames(sub_meta))
      }

      aligned_covs <- align_covariate_samples(
        sample_ids = rownames(sub_meta),
        known_covs = known_covs,
        geno_pcs = geno_sub,
        expr_pcs = pcs_all
      )
      common_cov_samples <- aligned_covs$sample_ids
      known_covs <- aligned_covs$known_covs
      geno_sub <- aligned_covs$geno_pcs
      PCsTop <- aligned_covs$expr_pcs

      if (ncol(geno_sub) > 0) {
        known_covs <- cbind(known_covs, geno_sub)
      }

      if (nrow(PCsTop) == 0) {
        warning(sprintf("[%s] No overlapping samples remained for covariate construction. Skipping.", stage))
        next
      }

      if (ncol(known_covs) > 0) {
        covsFiltered <- PCAForQTL::filterKnownCovariates(known_covs, PCsTop, unadjustedR2_cutoff = 0.9)
      } else {
        covsFiltered <- known_covs
      }
      covsToUse <- cbind(covsFiltered, PCsTop)

      out_filename <- sprintf("covariates_for_qtl.%s.txt", stage)
      fwrite(as.data.frame(t(covsToUse)), file.path(out_dir, out_filename), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

      sub_expr_final <- as.data.frame(sub_expr_norm)
      sub_expr_final$pid <- rownames(sub_expr_norm)
      setDT(sub_expr_final)
      sub_expr_final[, c("transcript", "gene", "Loci", "strand") := tstrsplit(pid, "|", fixed = TRUE)]
      sub_expr_final[, c("chr", "region") := tstrsplit(Loci, ":", fixed = TRUE)]
      sub_expr_final[, c("start", "end") := tstrsplit(region, "-", fixed = TRUE)]
      sub_expr_final[, `:=`(start = as.numeric(start), end = as.numeric(end), `#Chr` = chr)]
      final_samples <- rownames(covsToUse)
      if (!all(final_samples %in% colnames(sub_expr_final))) {
        stop("Critical Error: Mismatch between covariate samples and expression samples!")
      }
      format_cols <- c("#Chr", "start", "end", "pid", final_samples)

      format_data <- sub_expr_final[, ..format_cols]
      row_count_before <- nrow(format_data)
      format_data <- unique(format_data, by = "pid")
      row_count_after <- nrow(format_data)
      if (row_count_before > row_count_after) {
        message(sprintf("[%s] Removed %d duplicated phenotype rows.", stage, row_count_before - row_count_after))
      }
      setorder(format_data, `#Chr`, start, end)

      bed_file <- file.path(out_dir, sprintf("phenotype.%s.bed", stage))
      fwrite(format_data, bed_file, sep = '\t', quote = FALSE)

      has_chr_prefix <- detect_chr_prefix(paste0(pgen_prefix, ".pvar"))
      pheno_chr_sample <- as.character(format_data$`#Chr`[1])
      pheno_has_prefix <- grepl("^chr", pheno_chr_sample)
      if (has_chr_prefix != pheno_has_prefix) {
        if (has_chr_prefix) {
          format_data[, `#Chr` := paste0("chr", `#Chr`)]
        } else {
          format_data[, `#Chr` := sub("^chr", "", `#Chr`)]
        }
        fwrite(format_data, bed_file, sep = '\t', quote = FALSE)
      }

      run_command(sprintf("bgzip -f %s", bed_file), sprintf("[%s] bgzip failed for phenotype BED", stage))
      run_command(sprintf("tabix -p bed -f %s", paste0(bed_file, ".gz")),
                  sprintf("[%s] tabix indexing failed for phenotype BED", stage))

      sample_file <- file.path(out_dir, sprintf("sample.%s.final", stage))
      writeLines(final_samples, sample_file)
      message(sprintf("[%s] Outputs ready: %s, %s, %s.gz, %s",
                      stage, out_filename, bed_file, bed_file, sample_file))
    } else {
      warning(sprintf("[%s] Insufficient samples or phenotypes for covariate analysis. Skipping.", stage))
    }

    message(sprintf("[%s] Processing complete.", stage))
  }
}

if (sys.nframe() == 0) {
  main()
}
