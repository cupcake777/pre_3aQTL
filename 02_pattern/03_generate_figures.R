suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(ggrepel)
  library(limma)
  library(splines)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

stage_col  <- c(Stage1 = "#56B4E9", Stage2 = "#FFB482", Stage3 = "#8DE5A1")
cl_col     <- c("1" = "#D55E00", "2" = "#E69F00", "3" = "#4DBBD5", "4" = "#0072B2")
stage_labs <- c("Stage 1", "Stage 2", "Stage 3")
stage_x_labels <- c("Stage 1", "Stage 2", "Stage 3")
stage_levels   <- c("Stage1", "Stage2", "Stage3")

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE))
} else {
  normalizePath(getwd(), mustWork = FALSE)
}
figure_dir <- if (length(commandArgs(trailingOnly = TRUE)) >= 1) {
  normalizePath(commandArgs(trailingOnly = TRUE)[1], mustWork = FALSE)
} else {
  script_dir
}
if (!dir.exists(figure_dir)) dir.create(figure_dir, recursive = TRUE)

resolve_existing_path <- function(label, candidates) {
  candidates <- unique(candidates[nzchar(candidates)])
  for (candidate in candidates) {
    if (file.exists(candidate)) {
      return(normalizePath(candidate, mustWork = FALSE))
    }
  }
  stop(sprintf("Missing %s. Checked:\n%s", label, paste(candidates, collapse = "\n")))
}

resolve_latest_gsea_file <- function(script_dir, prefix, fallback_k = 3L) {
  candidates <- Sys.glob(file.path(script_dir, sprintf("%s_K*_all.csv", prefix)))
  if (length(candidates) == 0) {
    return(file.path(script_dir, sprintf("%s_K%d_all.csv", prefix, fallback_k)))
  }
  candidate_k <- suppressWarnings(as.integer(sub(sprintf("^%s_K([0-9]+)_all\\.csv$", prefix), "\\1", basename(candidates))))
  candidates[[order(candidate_k, decreasing = TRUE, na.last = TRUE)[1]]]
}

input_paths <- list(
  apa_matrix   = file.path(script_dir, "../01_pre/after.combat.txt"),
  iso_matrix   = file.path(script_dir, "../../raw_data/isoform.after.combat.txt"),
  metadata     = file.path(script_dir, "../../raw_data/CHB.all.meta"),
  gene2ensembl = resolve_existing_path(
    "gene2ensembl reference",
    c(
      Sys.getenv("GENE2ENSEMBL_FILE", unset = ""),
      file.path(script_dir, "../../ref/gene2ensembl.gz"),
      file.path(script_dir, "../../../ref/gene2ensembl.gz"),
      "/mnt/share_group_folder/ref/gene2ensembl.gz"
    )
  ),
  kegg_gsea    = resolve_latest_gsea_file(script_dir, "05b_GSEA_KEGG"),
  go_gsea      = resolve_latest_gsea_file(script_dir, "05b_GSEA_GO_BP"),
  cluster_core = file.path(script_dir, "06_final_CORE_cluster_assignments.csv"),
  cluster_summary = file.path(script_dir, "05_core_gene_filter_summary.csv"),
  representative = file.path(script_dir, "06_representative_core_genes.csv"),
  smr_eur_dir  = file.path(script_dir, "../smr/results/EUR"),
  smr_eas_dir  = file.path(script_dir, "../smr/results/EAS")
)

required_inputs <- c(
  input_paths$apa_matrix,
  input_paths$iso_matrix,
  input_paths$metadata,
  input_paths$gene2ensembl,
  input_paths$kegg_gsea,
  input_paths$go_gsea,
  input_paths$cluster_core,
  input_paths$cluster_summary,
  input_paths$representative
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) {
  stop(sprintf("03_generate_figures.R missing inputs:\n%s", paste(missing_inputs, collapse = "\n")))
}

figure_output_specs <- list(
  fig1 = list(file = "figure_01_trajectory_overview.pdf", width = 16, height = 6),
  fig2a = list(file = "figure_02a_peak_valley_enrichment.pdf", width = 14, height = 7),
  fig2b = list(file = "figure_02b_cluster2_enrichment.pdf", width = 11, height = 5),
  fig3 = list(file = "figure_03_cluster4_energy_program.pdf", width = 9, height = 4),
  fig4 = list(file = "figure_04_lifespan_apa_model.pdf", width = 13, height = 6.75),
  fig5 = list(file = "figure_05_cluster1_translation_program.pdf", width = 14, height = 6.5),
  fig6 = list(file = "figure_06_cluster2_metabolic_program.pdf", width = 14, height = 6.5),
  fig7 = list(file = "figure_07_cluster_trait_associations.pdf", width = 14, height = 7.5)
)

theme_ppt <- function(base_size = 13) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      axis.text        = element_text(size = base_size),
      axis.title       = element_text(size = base_size + 1, face = "bold"),
      plot.title       = element_text(size = base_size + 2, face = "bold", hjust = 0.5),
      plot.subtitle    = element_text(size = base_size - 1, hjust = 0.5, color = "gray40"),
      legend.text      = element_text(size = base_size - 1),
      legend.title     = element_text(size = base_size, face = "bold"),
      legend.key.size  = unit(0.7, "lines"),
      strip.text       = element_text(size = base_size, face = "bold"),
      plot.margin      = margin(8, 12, 8, 12)
    )
}

message("Loading shared data matrices...")

apa_raw  <- fread(input_paths$apa_matrix,  data.table = FALSE)
meta_raw <- fread(input_paths$metadata, data.table = FALSE)
rownames(apa_raw)  <- apa_raw[, 1];  apa_raw  <- apa_raw[, -1]
rownames(meta_raw) <- meta_raw[, 1]; meta_raw <- meta_raw[, -1]
common_apa <- intersect(colnames(apa_raw), rownames(meta_raw))
apa_mat    <- as.matrix(apa_raw[, common_apa, drop = FALSE])
meta_shared <- meta_raw[common_apa, , drop = FALSE]

sex_vec_sh    <- ifelse(meta_shared$sex == "Male", 1, 0)
rin_vec_sh    <- ns(meta_shared$RIN, df = 2)
design_prot_sh <- model.matrix(~ LifeStage, data = meta_shared)
apa_mat        <- removeBatchEffect(apa_mat,
                                    covariates = cbind(sex = sex_vec_sh, rin_vec_sh),
                                    design     = design_prot_sh)
stages_apa <- meta_shared$LifeStage
apa_rowids <- rownames(apa_mat)
apa_enst   <- sub("\\|.*", "", apa_rowids)
message(sprintf("  APA matrix: %d features x %d samples", nrow(apa_mat), ncol(apa_mat)))

iso_raw  <- fread(input_paths$iso_matrix, data.table = FALSE)
rownames(iso_raw) <- iso_raw[, 1]; iso_raw <- iso_raw[, -1]
common_iso <- intersect(colnames(iso_raw), common_apa)
iso_mat    <- as.matrix(iso_raw[, common_iso, drop = FALSE])
rownames(iso_mat) <- sub("\\..*", "", rownames(iso_mat))   # strip ENST version
stages_iso <- meta_shared[common_iso, "LifeStage"]
message(sprintf("  Isoform matrix: %d transcripts x %d samples", nrow(iso_mat), ncol(iso_mat)))

g2e_raw   <- fread(input_paths$gene2ensembl, data.table = FALSE, sep = "\t",
                   col.names = c("tax_id","GeneID","Ensembl_gene",
                                 "RNA_accession","Ensembl_rna",
                                 "Protein_accession","Ensembl_protein"))
g2e_human_all <- g2e_raw[g2e_raw$tax_id == 9606, ]
message("  gene2ensembl loaded.")

cluster_core_df <- read.csv(input_paths$cluster_core, stringsAsFactors = FALSE)
cluster_summary_df <- read.csv(input_paths$cluster_summary, stringsAsFactors = FALSE)
representative_df <- read.csv(input_paths$representative, stringsAsFactors = FALSE)

cluster_levels <- sort(unique(cluster_core_df$Cluster))
cluster_labels <- c("Peak", "Up", "Down", "Valley", "Flat")[seq_along(cluster_levels)]
cluster_name_map <- setNames(cluster_labels, cluster_levels)
cluster_display_map <- setNames(
  sprintf("Cluster %s\n(%s)", cluster_levels, cluster_name_map[as.character(cluster_levels)]),
  cluster_levels
)

stage_mean_mat <- t(sapply(stage_levels, function(st) {
  s_samples <- rownames(meta_shared)[meta_shared$LifeStage == st]
  if (length(s_samples) > 1) return(rowMeans(apa_mat[, s_samples, drop = FALSE], na.rm = TRUE))
  if (length(s_samples) == 1) return(apa_mat[, s_samples, drop = TRUE])
  rep(NA_real_, nrow(apa_mat))
}))
stage_mean_mat <- t(stage_mean_mat)
rownames(stage_mean_mat) <- rownames(apa_mat)
colnames(stage_mean_mat) <- stage_levels
stage_mean_z <- t(scale(t(stage_mean_mat)))
stage_mean_z <- stage_mean_z[complete.cases(stage_mean_z), , drop = FALSE]

centroids <- do.call(rbind, lapply(cluster_levels, function(cl) {
  genes <- intersect(cluster_core_df$Gene[cluster_core_df$Cluster == cl], rownames(stage_mean_z))
  ctr <- colMeans(stage_mean_z[genes, , drop = FALSE], na.rm = TRUE)
  data.frame(
    ClusterID = as.character(cl),
    Cluster = cluster_name_map[as.character(cl)],
    Stage = factor(stage_labs, levels = stage_labs),
    Zscore = ctr,
    N = sum(cluster_summary_df$Total_Genes[cluster_summary_df$Cluster == cl], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

save_plot <- function(plot_obj, spec, label) {
  output_file <- file.path(figure_dir, spec$file)
  ggsave(output_file, plot_obj, width = spec$width, height = spec$height, device = "pdf")
  message(sprintf("✓ %s saved: %s", label, output_file))
}

load_enrichment_tables <- function() {
  list(
    kegg = read.csv(input_paths$kegg_gsea, stringsAsFactors = FALSE),
    go   = read.csv(input_paths$go_gsea, stringsAsFactors = FALSE)
  )
}

parse_core_enrichment_ids <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) {
    return(integer(0))
  }
  unique(as.integer(unlist(strsplit(x, "/"))))
}

parse_core_symbol_ids <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) {
    return(character(0))
  }
  unique(unlist(strsplit(x, "/")))
}

symbols_to_entrez <- function(symbols) {
  symbols <- unique(symbols[!is.na(symbols) & nzchar(symbols)])
  if (length(symbols) == 0) {
    return(integer(0))
  }
  converted <- tryCatch(
    bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
    error = function(e) NULL
  )
  if (is.null(converted) || nrow(converted) == 0) {
    return(integer(0))
  }
  unique(as.integer(converted$ENTREZID))
}

build_program_data <- function(entrez_ids) {
  g2e_sub  <- g2e_human_all[g2e_human_all$GeneID %in% entrez_ids, ]
  g2e_sub  <- g2e_sub[!duplicated(paste(g2e_sub$GeneID, g2e_sub$Ensembl_rna)), ]

  g2e_sub$enst_id <- sub("\\..*", "", g2e_sub$Ensembl_rna)
  g2e_sub$apa_row_idx <- match(g2e_sub$enst_id, sub("\\..*", "", apa_enst))
  g2e_apa_sub <- g2e_sub[!is.na(g2e_sub$apa_row_idx), ]

  uniq_gids <- unique(g2e_apa_sub$GeneID)
  if (length(uniq_gids) == 0) {
    return(list(
      combined_long = data.frame(),
      median_line = data.frame(),
      n_genes = 0L,
      keep_genes = character(0),
      cor_df = data.frame()
    ))
  }
  apa_stage_df <- t(vapply(uniq_gids, function(gid) {
    rows <- g2e_apa_sub$apa_row_idx[g2e_apa_sub$GeneID == gid]
    mat  <- apa_mat[rows, , drop = FALSE]
    vapply(stage_levels, function(st)
      median(mat[, which(stages_apa == st)], na.rm = TRUE), numeric(1))
  }, numeric(length(stage_levels))))
  rownames(apa_stage_df) <- as.character(uniq_gids)
  colnames(apa_stage_df) <- stage_levels

  iso_rowids_nov <- rownames(iso_mat)
  iso_stage_df <- t(vapply(uniq_gids, function(gid) {
    enst_ids <- g2e_apa_sub$enst_id[g2e_apa_sub$GeneID == gid]
    rows     <- which(iso_rowids_nov %in% enst_ids)
    if (length(rows) == 0) return(rep(NA_real_, length(stage_levels)))
    mat <- iso_mat[rows, , drop = FALSE]
    vapply(stage_levels, function(st) {
      cols <- which(stages_iso == st)
      if (length(cols) == 0) return(NA_real_)
      median(mat[, cols], na.rm = TRUE)
    }, numeric(1))
  }, numeric(length(stage_levels))))
  rownames(iso_stage_df) <- as.character(uniq_gids)
  colnames(iso_stage_df) <- stage_levels

  ## Keep genes with complete data in both matrices
  keep <- rownames(apa_stage_df)[
    rownames(apa_stage_df) %in% rownames(iso_stage_df) &
    rowSums(!is.na(iso_stage_df[rownames(apa_stage_df), ])) == length(stage_levels)
  ]
  n_genes <- length(keep)
  if (n_genes == 0) {
    return(list(
      combined_long = data.frame(),
      median_line = data.frame(),
      n_genes = 0L,
      keep_genes = character(0),
      cor_df = data.frame()
    ))
  }

  ## Long data frame for per-gene lines
  combined_long <- data.frame(
    Stage      = factor(rep(stage_levels, each = n_genes), levels = stage_levels),
    GeneID     = rep(keep, times = length(stage_levels)),
    UTR_length = as.vector(apa_stage_df[keep, ]),
    Expression = as.vector(iso_stage_df[keep, ]),
    stringsAsFactors = FALSE
  )

  ## Population median
  median_line <- data.frame(
    Stage      = factor(stage_levels, levels = stage_levels),
    UTR_length = apply(apa_stage_df[keep, ], 2, median, na.rm = TRUE),
    Expression = apply(iso_stage_df[keep, ], 2, median, na.rm = TRUE)
  )

  iso_rowids_all <- rownames(iso_mat)
  uniq_gids_cor  <- unique(g2e_apa_sub$GeneID)
  cor_vals <- sapply(uniq_gids_cor, function(gid) {
    sub_rows <- g2e_apa_sub[g2e_apa_sub$GeneID == gid, ]
    for (i in seq_len(nrow(sub_rows))) {
      row_apa <- sub_rows$apa_row_idx[i]
      row_iso <- match(sub_rows$enst_id[i], iso_rowids_all)
      if (!is.na(row_apa) && !is.na(row_iso)) {
        apa_vals <- as.numeric(apa_mat[row_apa, common_iso])  # restrict to shared samples
        iso_vals <- as.numeric(iso_mat[row_iso, ])
        valid <- !is.na(apa_vals) & !is.na(iso_vals)
        if (sum(valid) >= 10)
          return(cor(apa_vals[valid], iso_vals[valid], method = "spearman"))
      }
    }
    NA_real_
  })
  cor_df <- data.frame(
    GeneID      = as.character(uniq_gids_cor),
    SpearmanRho = cor_vals,
    stringsAsFactors = FALSE
  )
  cor_df <- cor_df[!is.na(cor_df$SpearmanRho), ]

  list(combined_long = combined_long, median_line = median_line,
       n_genes = n_genes, keep_genes = keep, cor_df = cor_df)
}

summarize_rho_sign_test <- function(cor_df, direction = c("neg", "pos", "two.sided")) {
  direction <- match.arg(direction)
  n_neg <- sum(cor_df$SpearmanRho < 0, na.rm = TRUE)
  n_pos <- sum(cor_df$SpearmanRho > 0, na.rm = TRUE)
  n_tied <- sum(cor_df$SpearmanRho == 0, na.rm = TRUE)
  n_test <- n_neg + n_pos

  if (n_test == 0) {
    return(list(
      n_neg = n_neg,
      n_pos = n_pos,
      n_tied = n_tied,
      n_test = n_test,
      p_value = NA_real_,
      label = "no non-zero rho",
      x = NA_integer_
    ))
  }

  if (direction == "neg") {
    bt <- binom.test(x = n_neg, n = n_test, p = 0.5, alternative = "greater")
    label <- "neg > 50%"
    x <- n_neg
  } else if (direction == "pos") {
    bt <- binom.test(x = n_pos, n = n_test, p = 0.5, alternative = "greater")
    label <- "pos > 50%"
    x <- n_pos
  } else {
    bt <- binom.test(x = n_neg, n = n_test, p = 0.5, alternative = "two.sided")
    label <- "two-sided"
    x <- n_neg
  }

  list(
    n_neg = n_neg,
    n_pos = n_pos,
    n_tied = n_tied,
    n_test = n_test,
    p_value = bt$p.value,
    label = label,
    x = x
  )
}

build_dual_axis_program_figure <- function(vd, spec) {
  utr_med  <- setNames(vd$median_line$UTR_length, as.character(vd$median_line$Stage))
  expr_med <- setNames(vd$median_line$Expression, as.character(vd$median_line$Stage))

  if (spec$transform_mode == "inverted") {
    a <- (utr_med["Stage3"] - utr_med["Stage1"]) / (expr_med["Stage1"] - expr_med["Stage3"])
  } else {
    a <- (utr_med["Stage3"] - utr_med["Stage1"]) / (expr_med["Stage3"] - expr_med["Stage1"])
  }
  b <- utr_med["Stage1"] - a * expr_med["Stage1"]
  expr_to_utr <- function(x) a * x + b
  utr_to_expr <- function(y) (y - b) / a

  med <- vd$median_line
  med$Expr_scaled <- expr_to_utr(med$Expression)

  gene_long <- vd$combined_long
  gene_long$Expr_scaled <- expr_to_utr(gene_long$Expression)

  ribbon_utr <- gene_long %>%
    group_by(Stage) %>%
    summarise(
      q1 = quantile(UTR_length, 0.25, na.rm = TRUE),
      q3 = quantile(UTR_length, 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  yvals <- c(ribbon_utr$q1, ribbon_utr$q3, med$UTR_length, med$Expr_scaled)
  ypad <- diff(range(yvals, na.rm = TRUE)) * 0.14

  utr_end <- med$UTR_length[3]
  expr_end <- med$Expr_scaled[3]
  gap <- abs(utr_end - expr_end)
  nudge <- if (gap < 0.3) 0.20 else 0
  utr_label_y <- utr_end + spec$utr_label_nudge_sign * nudge
  expr_label_y <- expr_end + spec$expr_label_nudge_sign * nudge

  p_main <- ggplot() +
    geom_ribbon(
      data = ribbon_utr,
      aes(x = as.integer(Stage), ymin = q1, ymax = q3),
      fill = spec$col_utr, alpha = 0.20
    ) +
    geom_line(
      data = med,
      aes(x = as.integer(Stage), y = UTR_length, group = 1),
      color = spec$col_utr, linewidth = 2.2
    ) +
    geom_point(
      data = med,
      aes(x = as.integer(Stage), y = UTR_length),
      color = spec$col_utr, size = 5.5, shape = 16
    ) +
    geom_line(
      data = med,
      aes(x = as.integer(Stage), y = Expr_scaled, group = 1),
      color = spec$col_expr, linewidth = 2.2, linetype = "dashed"
    ) +
    geom_point(
      data = med,
      aes(x = as.integer(Stage), y = Expr_scaled),
      color = spec$col_expr, size = 5.5, shape = 17
    ) +
    annotate(
      "segment",
      x = 3.08, xend = 3.20,
      y = utr_label_y, yend = utr_label_y,
      color = spec$col_utr, linewidth = 0.9
    ) +
    annotate(
      "text",
      x = 3.22, y = utr_label_y,
      label = spec$utr_label_text,
      color = spec$col_utr, size = spec$label_size,
      fontface = "bold", hjust = 0, vjust = 0.5
    ) +
    annotate(
      "segment",
      x = 3.08, xend = 3.20,
      y = expr_label_y, yend = expr_label_y,
      color = spec$col_expr, linewidth = 0.9
    ) +
    annotate(
      "text",
      x = 3.22, y = expr_label_y,
      label = spec$expr_label_text,
      color = spec$col_expr, size = spec$label_size,
      fontface = "bold", hjust = 0, vjust = 0.5
    ) +
    annotate(
      "text",
      x = 0.68,
      y = min(yvals, na.rm = TRUE) + ypad * 0.15,
      label = "Shaded = 3'UTR IQR",
      color = spec$col_utr, size = 3.5,
      hjust = 0, fontface = "italic", alpha = 0.7
    ) +
    scale_x_continuous(
      breaks = 1:3,
      labels = stage_x_labels,
      limits = spec$x_limits,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      name = "3'UTR Length (corrected)",
      limits = c(min(yvals, na.rm = TRUE) - ypad, max(yvals, na.rm = TRUE) + ypad),
      sec.axis = sec_axis(~ utr_to_expr(.), name = "Isoform Expression (log2)")
    ) +
    labs(title = spec$panel_title, subtitle = spec$panel_subtitle, x = NULL) +
    theme_ppt() +
    theme(
      axis.title.y.left  = element_text(color = spec$col_utr, face = "bold"),
      axis.text.y.left   = element_text(color = spec$col_utr),
      axis.title.y.right = element_text(color = spec$col_expr, face = "bold"),
      axis.text.y.right  = element_text(color = spec$col_expr),
      panel.grid.major.y = element_line(color = "gray92", linewidth = 0.4)
    )

  rho_stats <- summarize_rho_sign_test(vd$cor_df, direction = spec$rho_direction)
  rho_med <- median(vd$cor_df$SpearmanRho)
  cor_sorted <- vd$cor_df[order(vd$cor_df$SpearmanRho), ]
  cor_sorted$Rank <- seq_len(nrow(cor_sorted))
  cor_sorted$Direction <- ifelse(
    cor_sorted$SpearmanRho < 0,
    "Negative (rho < 0)",
    "Positive (rho > 0)"
  )

  if (is.na(rho_stats$p_value)) {
    zero_count <- if (spec$rho_direction == "neg") rho_stats$n_neg else rho_stats$n_pos
    rho_label <- sprintf("%d / %d genes\nno non-zero rho", zero_count, rho_stats$n_test)
  } else if (spec$rho_direction == "neg") {
    rho_label <- sprintf("%d / %d genes\nrho < 0\np = %.3f",
                         rho_stats$n_neg, rho_stats$n_test, rho_stats$p_value)
  } else {
    rho_label <- sprintf("%d / %d genes\nrho > 0\np = %.3f",
                         rho_stats$n_pos, rho_stats$n_test, rho_stats$p_value)
  }

  fill_rect <- if (spec$rho_direction == "neg") {
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
             fill = spec$col_utr, alpha = 0.06)
  } else {
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = spec$col_utr, alpha = 0.06)
  }

  rho_median_hjust <- if (spec$rho_text_side == "left") 1 else 0
  rho_median_x <- if (spec$rho_text_side == "left") rho_med - 0.03 else rho_med + 0.03
  rho_count_hjust <- if (spec$rho_text_side == "left") 0 else 1
  rho_count_x <- if (spec$rho_text_side == "left") {
    min(vd$cor_df$SpearmanRho) + 0.02
  } else {
    max(vd$cor_df$SpearmanRho) - 0.02
  }
  positive_color <- if (!is.null(spec$positive_color)) spec$positive_color else "#E69F00"

  p_rho <- ggplot(cor_sorted, aes(x = SpearmanRho, y = Rank, color = Direction)) +
    fill_rect +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
    geom_vline(xintercept = rho_med, linetype = "solid", color = "gray20", linewidth = 1.2, alpha = 0.7) +
    geom_point(size = 4, alpha = 0.9) +
    annotate(
      "text",
      x = rho_median_x, y = nrow(cor_sorted) * 0.93,
      label = sprintf("median rho = %.2f", rho_med),
      size = 3.6, hjust = rho_median_hjust, color = "gray20", fontface = "italic"
    ) +
    annotate(
      "text",
      x = rho_count_x, y = 3.5,
      label = rho_label,
      size = 4, hjust = rho_count_hjust, color = spec$col_utr,
      fontface = "bold", lineheight = 1.3
    ) +
    scale_color_manual(
      values = c("Negative (rho < 0)" = spec$negative_color, "Positive (rho > 0)" = positive_color),
      name = NULL
    ) +
    scale_y_continuous(expand = expansion(add = c(0.5, 1))) +
    labs(
      title = "B  Per-gene Spearman rho",
      subtitle = "3'UTR length vs. Expression",
      x = "Spearman rho", y = "Gene (ranked)"
    ) +
    theme_ppt() +
    theme(
      legend.position = spec$legend_position,
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.text = element_text(size = 11)
    )

  (p_main + p_rho) +
    plot_layout(widths = spec$layout_widths) +
    plot_annotation(
      title = spec$plot_title,
      subtitle = spec$plot_subtitle,
      theme = theme(
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40")
      )
    )
}

## =============================================================================
## FIGURE 1 — Trajectory overview
## Panel A: centroid line plot  |  Panel B: cluster size bars  |  Panel C: representative genes
## =============================================================================

## ── Fig 1A: centroid trends ───────────────────────────────────────────────────
cent_long <- centroids %>%
  mutate(
    ClusterID = ClusterID,
    Label = sprintf("C%s: %s (n=%s)", ClusterID, Cluster, format(N, big.mark = ",")),
    Stage = factor(Stage, levels = stage_labs)
  )

p1a <- ggplot(cent_long, aes(x = Stage, y = Zscore,
                              color = ClusterID, group = ClusterID)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70", linewidth = 0.5) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 4) +
  scale_color_manual(values = cl_col[names(cl_col) %in% unique(cent_long$ClusterID)],
                     labels = unique(cent_long$Label),
                     name = NULL) +
  labs(title = "A  3'UTR Length Trajectories",
       x = NULL, y = "Centroid Z-score") +
  theme_ppt() +
  theme(legend.position = "bottom")

## ── Fig 1B: cluster sizes ─────────────────────────────────────────────────────
size_df <- cluster_summary_df %>%
  transmute(
    Cluster = factor(cluster_display_map[as.character(Cluster)],
                     levels = cluster_display_map[as.character(cluster_levels)]),
    ClusterID = as.character(Cluster),
    N = Total_Genes
  )

p1b <- ggplot(size_df, aes(x = Cluster, y = N, fill = ClusterID)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.4) +
  geom_text(aes(label = format(N, big.mark = ",")),
            vjust = -0.4, size = 4.5, fontface = "bold") +
  scale_fill_manual(values = cl_col, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "B  Cluster Size", x = NULL, y = "Number of APA events") +
  theme_ppt() +
  theme(axis.text.x = element_text(size = 12))

## ── Fig 1C: representative gene APA profiles ─────────────────────────────────
rep_genes <- representative_df %>%
  mutate(
    ClusterID = as.character(Cluster),
    TranscriptShort = sub("\\|.*", "", Gene),
    GeneLabel = ifelse(!is.na(GeneName) & nzchar(GeneName), GeneName, TranscriptShort),
    Panel = sprintf("%s\n(C%s %s)", GeneLabel, ClusterID, cluster_name_map[ClusterID])
  ) %>%
  rowwise() %>%
  do({
    gene_id <- .$Gene
    zvals <- stage_mean_z[gene_id, stage_levels]
    data.frame(
      Stage = factor(stage_labs, levels = stage_labs),
      Panel = .$Panel,
      ClusterID = .$ClusterID,
      Zscore = as.numeric(zvals),
      stringsAsFactors = FALSE
    )
  }) %>%
  ungroup()

p1c <- ggplot(rep_genes, aes(x = Stage, y = Zscore,
                              color = ClusterID, group = Panel)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70", linewidth = 0.5) +
  geom_line(linewidth = 1.4) +
  geom_point(size = 3.5) +
  facet_wrap(~ Panel, nrow = 1) +
  scale_color_manual(values = cl_col, guide = "none") +
  labs(title = "C  Representative Gene Profiles",
       x = NULL, y = "3'UTR Length Z-score") +
  theme_ppt() +
  theme(strip.background = element_rect(fill = "gray93", color = NA))

## ── Combine Fig 1 ────────────────────────────────────────────────────────────
fig1 <- (p1a | p1b | p1c) + plot_layout(widths = c(3, 1.5, 3))

save_plot(fig1, figure_output_specs$fig1, "Fig 1")

## =============================================================================
## FIGURE 2A — Peak vs Valley mirror enrichment
## =============================================================================

enrichment_tables <- load_enrichment_tables()
kegg <- enrichment_tables$kegg
go   <- enrichment_tables$go

shared_kegg_desc <- intersect(
  kegg$Description[kegg$Cluster == 1],
  kegg$Description[kegg$Cluster == 4]
)
shared_go_desc <- intersect(
  go$Description[go$Cluster == 1],
  go$Description[go$Cluster == 4]
)

kegg_select <- intersect(shared_kegg_desc, c(
  "Oxidative phosphorylation",
  "Pathways of neurodegeneration - multiple diseases",
  "Alzheimer disease",
  "Huntington disease"
))

go_select <- intersect(shared_go_desc, c(
  "ATP metabolic process",
  "macroautophagy",
  "process utilizing autophagic mechanism",
  "purine nucleoside triphosphate metabolic process"
))

butterfly_data <- bind_rows(
  kegg %>% filter(Description %in% kegg_select, Cluster %in% c(1, 4)) %>%
    mutate(DB = "KEGG"),
  go %>% filter(Description %in% go_select, Cluster %in% c(1, 4)) %>%
    mutate(DB = "GO:BP")
) %>%
  dplyr::select(Description, NES, p.adjust, Cluster, DB) %>%
  mutate(
    ClusterLabel = ifelse(Cluster == 1, "C1: Peak", "C4: Valley"),
    Significant  = p.adjust < 0.05
  )

butterfly_data <- butterfly_data %>%
  mutate(Description = recode(Description,
    "Pathways of neurodegeneration - multiple diseases" = "Neurodegeneration pathways",
    "ATP metabolic process" = "ATP metabolic process",
    "process utilizing autophagic mechanism" = "Autophagic mechanism",
    "purine nucleoside triphosphate metabolic process" = "Purine NTP metabolism"
  ))

term_order <- butterfly_data %>%
  filter(Cluster == 4) %>%
  arrange(NES) %>%
  pull(Description)

butterfly_data <- butterfly_data %>%
  mutate(Description = factor(Description, levels = term_order))

p2a <- ggplot(butterfly_data,
              aes(x = NES, y = Description,
                  fill = NES > 0)) +
  geom_col(aes(linewidth = ifelse(Significant, 0.8, 0)),
           color = "gray30", width = 0.7) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.8) +
  geom_text(data = butterfly_data %>% filter(Cluster == 4 & Significant),
            aes(x = NES + 0.08, label = "*"),
            hjust = 0, size = 5, color = "black") +
  geom_text(data = butterfly_data %>% filter(Cluster == 1 & Significant),
            aes(x = NES - 0.08, label = "*"),
            hjust = 1, size = 5, color = "black") +
  scale_fill_manual(values = c("TRUE" = "#E69F00", "FALSE" = "#4DBBD5"),
                    labels = c("TRUE" = "Enriched (NES > 0)",
                               "FALSE" = "Depleted (NES < 0)"),
                    name = NULL) +
  scale_linewidth_identity() +
  facet_wrap(~ ClusterLabel, ncol = 2) +
  labs(title = "Mirror-Image Enrichment: Peak vs Valley Programs",
       subtitle = "* p.adj < 0.05  |  All bars: p.adj < 0.20",
       x = "Normalized Enrichment Score (NES)", y = NULL) +
  theme_ppt() +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "gray93", color = NA),
        panel.grid.major.x = element_line(color = "gray90"))

save_plot(p2a, figure_output_specs$fig2a, "Fig 2A")

## =============================================================================
## FIGURE 2B — Cluster 2 functional enrichment (horizontal bar)
## =============================================================================

c2_kegg <- kegg %>% filter(Cluster == 2) %>% dplyr::select(Description, NES, p.adjust)
c2_go   <- go   %>% filter(Cluster == 2) %>% dplyr::select(Description, NES, p.adjust)

## Curated terms for C2
c2_terms <- bind_rows(
  c2_kegg %>% filter(Description %in% c(
    "Oxidative phosphorylation",
    "Glycolysis / Gluconeogenesis",
    "Chemical carcinogenesis - reactive oxygen species"
  )) %>% mutate(DB = "KEGG"),
  c2_go %>% filter(Description %in% c(
    "regulation of chromatin organization",
    "regulation of heterochromatin formation",
    "axon guidance",
    "neuron projection guidance"
  )) %>% mutate(DB = "GO:BP")
) %>%
  distinct(Description, .keep_all = TRUE) %>%
  mutate(
    Significant = p.adjust < 0.05,
    Direction   = ifelse(NES > 0, "Enriched in late-life elongation", "Depleted in late-life elongation"),
    Description = recode(Description,
      "regulation of chromatin organization" = "Chromatin organization",
      "regulation of heterochromatin formation" = "Heterochromatin formation",
      "neuron projection guidance" = "Neuron projection guidance",
      "Chemical carcinogenesis - reactive oxygen species" = "Reactive oxygen species"
    ),
    Description = factor(Description,
                         levels = Description[order(NES)])
  )

p2b <- ggplot(c2_terms, aes(x = NES, y = Description, fill = NES > 0)) +
  geom_col(aes(color = Significant),
           linewidth = 0.7, width = 0.65) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.8) +
  scale_fill_manual(values = c("TRUE" = "#D55E00", "FALSE" = "#999999"),
                    labels = c("TRUE" = "Chromatin / guidance axis",
                               "FALSE" = "Energy metabolism axis"),
                    name = "3'UTR Direction") +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "gray70"),
                     guide = "none") +
  annotate("text", x = 2.7, y = 0.8,
           label = "Positive: chromatin organization and axon guidance\nNegative: oxidative and glycolytic metabolism",
           hjust = 1, size = 3.8, color = "gray30", fontface = "italic") +
  labs(title = "Cluster 2 (Up): Functional Enrichment",
       subtitle = "Bold border = p.adj < 0.05  |  Mechanism-first summary for report panels",
       x = "Normalized Enrichment Score (NES)", y = NULL) +
  theme_ppt() +
  theme(legend.position = "right",
        panel.grid.major.x = element_line(color = "gray90"))

save_plot(p2b, figure_output_specs$fig2b, "Fig 2B")

message("Fig 3: building Cluster 4 energy-program panel...")
c4_energy_ids <- unique(c(
  parse_core_enrichment_ids(
    kegg$core_enrichment[kegg$Cluster == 4 &
      kegg$Description %in% c("Oxidative phosphorylation", "Glycolysis / Gluconeogenesis")]
  ),
  symbols_to_entrez(parse_core_symbol_ids(
    go$core_enrichment[go$Cluster == 4 &
      go$Description %in% c(
        "generation of precursor metabolites and energy",
        "ATP metabolic process",
        "macroautophagy",
        "process utilizing autophagic mechanism"
      )]
  ))
))

vd3 <- build_program_data(c4_energy_ids)
message(sprintf("  Cluster 4 energy program: %d genes with APA + isoform data", vd3$n_genes))

fig3 <- build_dual_axis_program_figure(vd3, list(
  transform_mode = "inverted",
  col_utr = "#2980b9",
  col_expr = "#c0392b",
  negative_color = "#2980b9",
  positive_color = "#E69F00",
  x_limits = c(0.6, 4),
  utr_label_nudge_sign = -1,
  expr_label_nudge_sign = 1,
  utr_label_text = "3'UTR shortens",
  expr_label_text = "Expression rises",
  label_size = 3.5,
  panel_title = "A  3'UTR Shortening Marks the Valley Energy-homeostasis Program",
  panel_subtitle = NULL,
  rho_direction = "neg",
  rho_text_side = "left",
  legend_position = "top",
  layout_widths = c(1.2, 1),
  plot_title = "Energy Program",
  plot_subtitle = "Cluster 4 (Valley) genes linking mitochondrial ATP programs and autophagy-related processes"
))

save_plot(fig3, figure_output_specs$fig3, "Fig 3")

xmin_s1 <- 0.0;  xmax_s1 <- 1.0
xmin_s2 <- 1.0;  xmax_s2 <- 2.5
xmin_s3 <- 2.5;  xmax_s3 <- 4.0

x_pts <- c(0.5, 1.75, 3.25)

curve_df <- data.frame(
  x        = rep(x_pts, 4),
  y        = c(-0.99,  0.87,  0.12,
               -1.10,  0.42,  0.68,
                1.01, -0.35, -0.66,
                1.15, -0.62, -0.53),
  Cluster  = rep(c("C1: Peak", "C2: Up", "C3: Down", "C4: Valley"), each = 3),
  ClusterID = rep(c("1", "2", "3", "4"), each = 3)
)

rect_df <- data.frame(
  xmin  = c(xmin_s1, xmin_s2, xmin_s3),
  xmax  = c(xmax_s1, xmax_s2, xmax_s3),
  Stage = c("Stage1", "Stage2", "Stage3")
)

annot_cl1 <- data.frame(x = 1.75, y = 1.35,
  label = "Translation / ribosome\n(transient postnatal peak)")
annot_cl2 <- data.frame(x = 3.25, y = 1.10,
  label = "Chromatin organization\nAxon guidance\n(progressive elongation)")
annot_cl3 <- data.frame(x = 0.5, y = 1.55,
  label = "Energy-homeostasis\nmitochondrial ATP + autophagy\n(valley trajectory)")
annot_cl4 <- data.frame(x = 3.25, y = -0.95,
  label = "Down trajectory\nrespiratory recovery\nwithout a dominant report example")

p4 <- ggplot() +
  ## stage background bands
  geom_rect(data = rect_df,
            aes(xmin = xmin, xmax = xmax, ymin = -1.6, ymax = 1.7, fill = Stage),
            alpha = 0.18, inherit.aes = FALSE) +
  ## stage boundary lines
  geom_vline(xintercept = c(xmin_s2, xmin_s3),
             linetype = "dashed", color = "gray50", linewidth = 0.7) +
  geom_line(data = curve_df,
            aes(x = x, y = y, color = ClusterID, group = Cluster),
            linewidth = 2.0) +
  geom_point(data = curve_df,
             aes(x = x, y = y, color = ClusterID),
             size = 3.5) +
  ## annotation boxes
  geom_label(data = annot_cl1,
             aes(x = x, y = y, label = label),
             fill = scales::alpha(cl_col["1"], 0.15),
             color = cl_col["1"], size = 3.8, fontface = "italic",
             label.padding = unit(0.4, "lines"), linewidth = 0.4,
             inherit.aes = FALSE) +
  geom_label(data = annot_cl2,
             aes(x = x, y = y, label = label),
             fill = scales::alpha(cl_col["2"], 0.15),
             color = cl_col["2"], size = 3.8, fontface = "italic",
             label.padding = unit(0.4, "lines"), linewidth = 0.4,
             inherit.aes = FALSE) +
  geom_label(data = annot_cl3,
             aes(x = x, y = y, label = label),
             fill = scales::alpha(cl_col["3"], 0.15),
             color = cl_col["3"], size = 3.8, fontface = "italic",
             label.padding = unit(0.4, "lines"), linewidth = 0.4,
             inherit.aes = FALSE) +
  geom_label(data = annot_cl4,
             aes(x = x, y = y, label = label),
             fill = scales::alpha(cl_col["4"], 0.15),
             color = cl_col["4"], size = 3.5, fontface = "italic",
             label.padding = unit(0.35, "lines"), linewidth = 0.4,
             inherit.aes = FALSE) +
  ## stage labels at top
  annotate("text", x = 0.5,  y = -1.55, label = "Stage 1\n(Prenatal)",      size = 4, fontface = "bold", color = "gray30") +
  annotate("text", x = 1.75, y = -1.55, label = "Stage 2\n(0-65 yr)",        size = 4, fontface = "bold", color = "gray30") +
  annotate("text", x = 3.25, y = -1.55, label = "Stage 3\n(>65 yr)",         size = 4, fontface = "bold", color = "gray30") +
  ## timeline arrow
  annotate("segment", x = 0.05, xend = 3.95, y = -1.9, yend = -1.9,
           arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
           color = "gray20", linewidth = 1.0) +
  annotate("text", x = 2.0, y = -2.05,
           label = "Human Lifespan", size = 4.2, color = "gray20", fontface = "bold") +
  ## cluster legend labels
  scale_color_manual(values = cl_col,
                     labels = c("1" = "C1: Peak", "2" = "C2: Up", "3" = "C3: Down", "4" = "C4: Valley"),
                     name = "APA Cluster") +
  scale_fill_manual(values = stage_col, guide = "none") +
  scale_x_continuous(limits = c(0, 4), expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(-2.15, 1.9)) +
  labs(title = "Model of APA Regulation across the Human Lifespan",
       subtitle = "Each cluster represents a distinct 3'UTR length trajectory linked to specific biological programs",
       x = NULL, y = "3'UTR Length (Z-score)") +
  theme_ppt() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    panel.grid   = element_blank()
  )

save_plot(p4, figure_output_specs$fig4, "Fig 4")

message("Fig 5: building Cluster 1 translation-program panel...")

c1_ids <- unique(c(
  parse_core_enrichment_ids(
    kegg$core_enrichment[kegg$Cluster == 1 & kegg$Description %in% c("Ribosome")]
  ),
  symbols_to_entrez(parse_core_symbol_ids(
    go$core_enrichment[go$Cluster == 1 & go$Description %in% c("cytoplasmic translation")]
  ))
))

vd5 <- build_program_data(c1_ids)
message(sprintf("  Cluster 1 translation program: %d genes with APA + isoform data", vd5$n_genes))
fig5 <- build_dual_axis_program_figure(vd5, list(
  transform_mode = "forward",
  col_utr = "#D55E00",
  col_expr = "#8B4000",
  negative_color = "#999999",
  positive_color = "#D55E00",
  x_limits = c(0.6, 3.75),
  utr_label_nudge_sign = 1,
  expr_label_nudge_sign = -1,
  utr_label_text = "3'UTR",
  expr_label_text = "Expression",
  label_size = 4.5,
  panel_title = "A  Translation-associated Genes Peak in the Peak Program",
  panel_subtitle = paste0(
    "n = ", length(vd5$keep_genes),
    " Cluster 1 genes (ribosome + cytoplasmic translation)  |  median 3'UTR (+ IQR band)  |  median Expression"
  ),
  rho_direction = "pos",
  rho_text_side = "right",
  legend_position = "bottom",
  layout_widths = c(3, 1.4),
  plot_title = "Translation Program",
  plot_subtitle = "Cluster 1 (Peak) genes - transient 3'UTR lengthening accompanies a translation-linked program"
))

save_plot(fig5, figure_output_specs$fig5, "Fig 5")

message("Fig 6: building Cluster 2 metabolic-program panel...")

c2_ids <- unique(parse_core_enrichment_ids(
  kegg$core_enrichment[kegg$Cluster == 2 &
    kegg$Description %in% c(
      "Oxidative phosphorylation",
      "Glycolysis / Gluconeogenesis",
      "Chemical carcinogenesis - reactive oxygen species"
    )]
))

vd6 <- build_program_data(c2_ids)
message(sprintf("  Cluster 2 metabolic program: %d genes with APA + isoform data", vd6$n_genes))
fig6 <- build_dual_axis_program_figure(vd6, list(
  transform_mode = "inverted",
  col_utr = "#E69F00",
  col_expr = "#B07000",
  negative_color = "#E69F00",
  positive_color = "#999999",
  x_limits = c(0.6, 3.75),
  utr_label_nudge_sign = -1,
  expr_label_nudge_sign = 1,
  utr_label_text = "3'UTR down",
  expr_label_text = "Expression up",
  label_size = 4.5,
  panel_title = "A  Metabolic Remodeling Genes Shorten in the Up Program",
  panel_subtitle = paste0(
    "n = ", length(vd6$keep_genes),
    " Cluster 2 genes (oxidative phosphorylation + glycolysis)  |  median 3'UTR (+ IQR band)  |  median Expression"
  ),
  rho_direction = "neg",
  rho_text_side = "left",
  legend_position = "bottom",
  layout_widths = c(3, 1.4),
  plot_title = "Metabolic Program",
  plot_subtitle = "Cluster 2 (Up) genes - 3'UTR shortening highlights a residual metabolic axis beneath the dominant elongation trajectory"
))

save_plot(fig6, figure_output_specs$fig6, "Fig 6")

## =============================================================================
## FIG 7 — Systematic APA cluster–disease association via SMR
## Strategy: load ALL SMR results (EUR + EAS, 13 traits, 4 stages),
##   apply standard filters (p_SMR < 0.05, p_HEIDI > 0.05),
##   then test enrichment of each trait's significant genes in each cluster
##   using Fisher's exact test (background = all 12,261 core cluster genes).
##   No disease pre-selection — all 13 GWAS traits are tested.
## =============================================================================
message("\n--- Fig 7: Systematic SMR-based disease–cluster enrichment ---")

## ── Helper: load all .smr files from one directory ───────────────────────────
load_smr_dir <- function(smr_dir, pop) {
  files <- list.files(smr_dir, pattern = "[.]smr$", full.names = TRUE)
  lst <- lapply(files, function(f) {
    bn    <- basename(f)
    parts <- strsplit(sub("[.]smr$", "", bn), "_")[[1]]
    ## filename: {POP}_{TRAIT}_{stage}.smr  (POP already stripped from dirname)
    trait <- paste(parts[2:(length(parts) - 1)], collapse = "_")
    stage <- parts[length(parts)]
    d <- tryCatch(read.table(f, header = TRUE, sep = "\t"),
                  error = function(e) NULL)
    if (is.null(d) || nrow(d) == 0) return(NULL)
    d$Trait <- trait
    d$Stage <- stage
    d$Pop   <- pop
    d
  })
  do.call(rbind, Filter(Negate(is.null), lst))
}

## ── Load SMR data ─────────────────────────────────────────────────────────────
smr_eur <- load_smr_dir(input_paths$smr_eur_dir, "EUR")
smr_eas <- load_smr_dir(input_paths$smr_eas_dir, "EAS")
all_smr <- rbind(smr_eur, smr_eas)

## Standard SMR significance filter (p_SMR < 0.05 AND p_HEIDI > 0.05)
smr_sig <- all_smr[
  all_smr$p_SMR < 0.05 &
    !is.na(all_smr$p_HEIDI) & all_smr$p_HEIDI > 0.05, ]

## ── Cluster gene names ────────────────────────────────────────────────────────
cluster_df <- read.csv(input_paths$cluster_core,
                       stringsAsFactors = FALSE)
cluster_df$gene_name <- sub("^[^|]+[|]([^|]+)[|].*", "\\1",
                             cluster_df$Gene)
total_genes <- nrow(cluster_df)

## ── Per-trait, per-cluster Fisher enrichment ─────────────────────────────────
traits <- sort(unique(smr_sig$Trait))

enrich_lst <- lapply(traits, function(tr) {
  sig_genes <- unique(smr_sig[smr_sig$Trait == tr, "Gene"])
  n_sig_total <- length(sig_genes)
  lapply(1:3, function(cl) {
    cl_genes   <- cluster_df[cluster_df$Cluster == cl, "gene_name"]
    not_cl     <- cluster_df[cluster_df$Cluster != cl, "gene_name"]
    a <- sum(cl_genes %in% sig_genes)              # cluster ∩ sig
    b <- length(cl_genes) - a                      # cluster ∩ not-sig
    c <- sum(not_cl %in% sig_genes)                # ¬cluster ∩ sig
    d <- length(not_cl) - c                        # ¬cluster ∩ not-sig
    ft <- fisher.test(matrix(c(a, b, c, d), 2, 2), alternative = "greater")
    data.frame(
      Trait            = tr,
      Cluster          = cl,
      OR               = as.numeric(ft$estimate),
      P                = ft$p.value,
      n_sig_in_cluster = a,
      n_cluster        = length(cl_genes),
      n_sig_total      = n_sig_total
    )
  })
})
fig7_res <- do.call(rbind, do.call(c, enrich_lst))
fig7_res$FDR <- p.adjust(fig7_res$P, "BH")

## ── Trait labels (display-friendly names) ────────────────────────────────────
trait_labels <- c(
  AD   = "Alzheimer's disease",
  ADHD = "ADHD",
  AIS  = "Ischemic stroke",
  ALS  = "ALS",
  ANX  = "Anxiety",
  ASD  = "Autism",
  BIP  = "Bipolar disorder",
  INS  = "Insomnia",
  INT  = "Intelligence",
  MDD  = "Depression",
  PD   = "Parkinson's disease",
  PTSD = "PTSD",
  SCZ  = "Schizophrenia"
)

fig7_res$Trait_label <- trait_labels[fig7_res$Trait]
fig7_res$Trait_label[is.na(fig7_res$Trait_label)] <- fig7_res$Trait[is.na(fig7_res$Trait_label)]

## ── Cluster factor labels ────────────────────────────────────────────────────
fig7_res$Cluster_label <- factor(
  c("1" = "C1\n(Peak)", "2" = "C2\n(Up)", "3" = "C3\n(Down)")[
    as.character(fig7_res$Cluster)],
  levels = c("C1\n(Peak)", "C2\n(Up)", "C3\n(Down)")
)

## ── Order traits by max OR across all clusters ───────────────────────────────
trait_order <- fig7_res %>%
  group_by(Trait_label) %>%
  summarise(max_or = max(OR, na.rm = TRUE), .groups = "drop") %>%
  arrange(max_or) %>%
  pull(Trait_label)
fig7_res$Trait_label <- factor(fig7_res$Trait_label, levels = trait_order)

## ── Significance labels ───────────────────────────────────────────────────────
fig7_res$sig_label <- ifelse(fig7_res$FDR < 0.001, "***",
                      ifelse(fig7_res$FDR < 0.01,  "**",
                      ifelse(fig7_res$FDR < 0.05,  "*",
                      ifelse(fig7_res$FDR < 0.10,  "·", ""))))

## ── Cap OR for display (colour scale) ────────────────────────────────────────
fig7_res$log2OR <- log2(pmax(fig7_res$OR, 0.01))
fig7_res$log2OR_capped <- pmin(pmax(fig7_res$log2OR, -2), 4)

## ── Panel A: Heatmap (trait × cluster, colour = log2OR, dot = sig) ───────────
cl_x_col <- c(
  "C1\n(Peak)"  = "#D55E00",
  "C2\n(Up)"    = "#E69F00",
  "C3\n(Down)"  = "#4DBBD5"
)

pA <- ggplot(fig7_res,
             aes(x = Cluster_label, y = Trait_label)) +
  geom_tile(aes(fill = log2OR_capped), color = "white", linewidth = 0.8) +
  ## significance star inside tile
  geom_text(aes(label = sig_label),
            size = 5, color = "white", fontface = "bold", vjust = 0.5) +
  scale_fill_gradient2(
    low      = "#2166AC",
    mid      = "#F7F7F7",
    high     = "#B2182B",
    midpoint = 0,
    limits   = c(-2, 4),
    oob      = scales::squish,
    name     = "log2(OR)",
    guide    = guide_colorbar(
      barheight = unit(4.5, "cm"),
      barwidth  = unit(0.5, "cm"),
      title.position = "top"
    )
  ) +
  scale_x_discrete(expand = expansion(mult = 0.05)) +
  scale_y_discrete(expand = expansion(mult = 0.04)) +
  labs(
    title    = "APA Trajectory Clusters Show Broad Brain Disease Relevance",
    subtitle = paste0(
      "SMR (p < 0.05, p_HEIDI > 0.05) | Fisher's exact test (one-sided) | ",
      "BH-FDR correction\n",
      "Background: all ", total_genes, " core APA genes | ",
      "13 GWAS traits, EUR + EAS populations"
    ),
    x       = NULL,
    y       = NULL,
    caption = "  * FDR<0.05   ** FDR<0.01   *** FDR<0.001   · FDR<0.10"
  ) +
  theme_ppt(base_size = 13) +
  theme(
    panel.grid      = element_blank(),
    axis.line       = element_blank(),
    axis.ticks      = element_blank(),
    axis.text.x     = element_text(size = 13, face = "bold",
                                   color = cl_x_col[
                                     c("C1\n(Peak)", "C2\n(Up)", "C3\n(Down)")]),
    axis.text.y     = element_text(size = 12),
    legend.position = "right",
    plot.title      = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(size = 9.5, color = "gray40", hjust = 0.5,
                                   lineheight = 1.3),
    plot.caption    = element_text(size = 9, color = "gray45", hjust = 0),
    plot.margin     = margin(10, 10, 8, 10)
  )

## ── Panel B: dot plot — OR (x) × trait (y), facet by cluster ─────────────────
## Keep only OR > 1 rows for clarity; colour by -log10(FDR)
fig7_dot <- fig7_res %>%
  filter(OR > 1) %>%
  mutate(
    neg_log10_fdr = pmin(-log10(FDR + 1e-10), 6),
    sig           = FDR < 0.05
  )

pB <- ggplot(fig7_dot,
             aes(x = OR, y = Trait_label)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "grey60", linewidth = 0.5) +
  geom_point(aes(color = neg_log10_fdr, size = n_sig_in_cluster),
             alpha = 0.85) +
  ## star for sig hits
  geom_text(data = filter(fig7_dot, FDR < 0.05),
            aes(x = OR * 1.15, label = sig_label),
            size = 4.5, fontface = "bold", color = "#8B0000", hjust = 0) +
  scale_color_gradient(
    low    = "grey80",
    high   = "#B2182B",
    name   = expression(-log[10](FDR)),
    limits = c(0, 6), oob = scales::squish,
    guide  = guide_colorbar(barheight = unit(3, "cm"),
                            barwidth  = unit(0.4, "cm"),
                            title.position = "top")
  ) +
  scale_size_continuous(
    range  = c(2, 8),
    name   = "SMR hits\nin cluster",
    breaks = c(10, 30, 60)
  ) +
  scale_x_log10(
    breaks = c(1, 2, 5, 10, 25),
    labels = c("1", "2", "5", "10", "25"),
    expand = expansion(mult = c(0.05, 0.18))
  ) +
  facet_wrap(~ Cluster_label, nrow = 1, scales = "free_x") +
  labs(
    title    = "Odds Ratios (enrichment only)",
    subtitle = "Colour = significance   Size = SMR gene count in cluster",
    x        = "Odds Ratio",
    y        = NULL
  ) +
  theme_ppt(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(size = 12, face = "bold"),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    axis.line.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    axis.text.y      = element_text(size = 11),
    legend.position  = "right",
    plot.title       = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(size = 9, color = "gray40", hjust = 0.5),
    plot.margin      = margin(10, 10, 8, 4)
  )

## ── Assemble 16:9 ─────────────────────────────────────────────────────────────
fig7 <- (pA | pB) +
  plot_layout(widths = c(1, 1.8)) +
  plot_annotation(
      title    = "APA Trajectory Clusters Are Systematically Associated with Brain Diseases",
    subtitle = paste0(
      "Genetic evidence from SMR across 13 neuropsychiatric/neurodegenerative GWAS traits\n",
      "C3 (Down: 3'UTR shortening in aging) shows the broadest disease enrichment"
    ),
    theme = theme(
      plot.title    = element_text(size = 15, face = "bold", hjust = 0.5,
                                   margin = margin(b = 3)),
      plot.subtitle = element_text(size = 10.5, color = "gray35", hjust = 0.5,
                                   lineheight = 1.3, margin = margin(b = 6))
    )
  )

save_plot(fig7, figure_output_specs$fig7, "Fig 7")

## =============================================================================
message("\n=== Figure export completed ===")
message("  figure_01_trajectory_overview.pdf")
message("  figure_02a_peak_valley_enrichment.pdf")
message("  figure_02b_cluster2_enrichment.pdf")
message("  figure_03_cluster4_energy_program.pdf")
message("  figure_04_lifespan_apa_model.pdf")
message("  figure_05_cluster1_translation_program.pdf")
message("  figure_06_cluster2_metabolic_program.pdf")
message("  figure_07_cluster_trait_associations.pdf")
