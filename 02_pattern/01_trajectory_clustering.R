#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(reshape2)
  library(patchwork)
  library(tibble)
  library(tidyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(stringr)
  library(factoextra)
  library(mclust)
  library(ComplexHeatmap)
  library(circlize)
  library(ggrepel)
  library(limma)
  library(splines)
  library(cluster)
  library(Mfuzz)
  library(Biobase)
  library(scales)
  library(fgsea)
})
args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE))
} else {
  normalizePath(getwd(), mustWork = FALSE)
}

OUTPUT_DIR <- if (length(args) >= 1) {
  normalizePath(args[1], mustWork = FALSE)
} else {
  script_dir
}
output_dir <- OUTPUT_DIR
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

resolve_existing_path <- function(label, candidates) {
  candidates <- unique(candidates[nzchar(candidates)])
  for (candidate in candidates) {
    if (file.exists(candidate)) {
      return(normalizePath(candidate, mustWork = FALSE))
    }
  }
  stop(sprintf("Missing %s. Checked:\n%s", label, paste(candidates, collapse = "\n")))
}

expr_file <- resolve_existing_path(
  "APA matrix",
  c(file.path(script_dir, "../01_pre/after.combat.txt"))
)
meta_file <- resolve_existing_path(
  "metadata",
  c(file.path(script_dir, "../../raw_data/CHB.all.meta"))
)
ISO_MATRIX_FILE <- resolve_existing_path(
  "isoform matrix",
  c(
    Sys.getenv("ISO_MATRIX_FILE", unset = ""),
    file.path(script_dir, "../../raw_data/isoform.after.combat.txt")
  )
)
GENE2ENSEMBL_FILE <- resolve_existing_path(
  "gene2ensembl reference",
  c(
    Sys.getenv("GENE2ENSEMBL_FILE", unset = ""),
    file.path(script_dir, "../../ref/gene2ensembl.gz"),
    file.path(script_dir, "../../../ref/gene2ensembl.gz"),
    "/mnt/share_group_folder/ref/gene2ensembl.gz"
  )
)

K_MAX                   <- 10
CORE_GENE_COR_THRESHOLD <- 0.8
CORE_GENE_MIN_SIZE      <- 10
STAGE_FDR_CUTOFF        <- 0.05

cluster_palette <- c(
  "1" = "#56B4E9", "2" = "#FFB482", "3" = "#8DE5A1", "4" = "#FEC44F",
  "5" = "#B2ABD2", "6" = "#80B1D3", "7" = "#9F79EE", "8" = "#F4A582"
)

batch_col <- c(
  batch1 = "#284F56", batch2 = "#2494BA", batch3 = "#28AFB0", batch4 = "#83BFBF",
  batch5 = "#F38F2C", batch6 = "#F0BA1A", batch7 = "#BDA055", batch8 = "#EAD291",
  batch9 = "#FFB482"
)

stage_col <- c(Stage1 = "#56B4E9", Stage2 = "#FFB482", Stage3 = "#8DE5A1")

get_cluster_colors <- function(k) {
  if (k <= length(cluster_palette)) {
    cluster_palette[seq_len(k)]
  } else {
    setNames(colorRampPalette(unname(cluster_palette))(k), as.character(seq_len(k)))
  }
}

if (!file.exists(expr_file) || !file.exists(meta_file)) stop("Check your file path")
expr_raw <- fread(expr_file, data.table = FALSE)
meta_raw <- fread(meta_file, data.table = FALSE)
rownames(expr_raw) <- expr_raw[, 1]; expr_raw <- expr_raw[, -1]
rownames(meta_raw) <- meta_raw[, 1]; meta_raw <- meta_raw[, -1]
common_samples <- intersect(colnames(expr_raw), rownames(meta_raw))
if (length(common_samples) == 0) stop("No shared sample, check file input!")
expr_final  <- expr_raw[, common_samples, drop = FALSE]
meta_final  <- meta_raw[common_samples, , drop = FALSE]
normalize_sex <- function(x) {
  x_chr <- trimws(as.character(x))
  ifelse(
    x_chr %in% c("Male", "M", "1"),
    "Male",
    ifelse(x_chr %in% c("Female", "F", "0"), "Female", x_chr)
  )
}
meta_final$sex_label <- normalize_sex(meta_final$sex)
data_matrix     <- as.matrix(expr_final)
data_matrix_raw <- data_matrix
cat("Matrix dimension", nrow(data_matrix), "features,", ncol(data_matrix), "samples\n")

sex_vec    <- ifelse(meta_final$sex_label == "Male", 1, 0)
rin_vec    <- ns(meta_final$RIN, df = 2)
covariates <- cbind(sex = sex_vec, rin_vec)
design_protect <- model.matrix(~ LifeStage, data = meta_final)
data_matrix <- removeBatchEffect(data_matrix,
                                 covariates = covariates,
                                 design     = design_protect)

make_pca_plot <- function(pca_df, pct_var, color_var, color_vals,
                          shape_var = NULL, title, subtitle) {
  aes_base <- if (!is.null(shape_var)) {
    aes(x = PC1, y = PC2,
        color = .data[[color_var]], shape = .data[[shape_var]])
  } else {
    aes(x = PC1, y = PC2, color = .data[[color_var]])
  }
  p <- ggplot(pca_df, aes_base) +
    geom_point(size = 2.2, alpha = 0.85) +
    scale_color_manual(values = color_vals, name = color_var, na.value = "grey70") +
    labs(title    = title,
         subtitle = subtitle,
         x = sprintf("PC1 (%.1f%%)", pct_var[1]),
         y = sprintf("PC2 (%.1f%%)", pct_var[2])) +
    theme_bw(base_size = 11) +
    theme(legend.position  = "right",
          legend.key.size  = unit(0.45, "cm"),
          legend.text      = element_text(size = 8),
          panel.grid.minor = element_blank(),
          plot.title       = element_text(face = "bold", size = 11),
          plot.subtitle    = element_text(size = 8, color = "grey40"))
  if (!is.null(shape_var))
    p <- p + scale_shape_manual(values = c(16, 17, 15, 18, 8),
                                 name = shape_var)
  p
}

pca_before    <- prcomp(t(data_matrix_raw), center = TRUE, scale. = FALSE)
pct_before    <- round(100 * pca_before$sdev^2 / sum(pca_before$sdev^2), 1)
pca_before_df <- as.data.frame(pca_before$x[, 1:2])
pca_before_df$sample    <- rownames(pca_before_df)
pca_before_df$LifeStage <- meta_final[pca_before_df$sample, "LifeStage"]
pca_before_df$batch     <- if ("batch" %in% colnames(meta_final))
  meta_final[pca_before_df$sample, "batch"] else "unknown"
pca_before_df$sex <- meta_final[pca_before_df$sample, "sex_label"]
pca_before_df$RIN <- as.numeric(meta_final[pca_before_df$sample, "RIN"])

p_before_stage <- make_pca_plot(
  pca_before_df, pct_before,
  color_var  = "LifeStage", color_vals = stage_col,
  title      = "Before correction — LifeStage",
  subtitle   = "Color: LifeStage"
)
p_before_batch <- make_pca_plot(
  pca_before_df, pct_before,
  color_var  = "batch", color_vals = batch_col,
  title      = "Before correction — Batch",
  subtitle   = "Color: Batch"
)
p_before_RIN <- ggplot(pca_before_df, aes(x = PC1, y = PC2, color = RIN)) +
  geom_point(size = 2.2, alpha = 0.85) +
  scale_color_gradient2(low = "#4575B4", mid = "#FFFFBF", high = "#D73027",
                        midpoint = median(pca_before_df$RIN), name = "RIN") +
  labs(title    = "Before correction — RIN",
       subtitle = "Color: RIN (continuous)",
       x = sprintf("PC1 (%.1f%%)", pct_before[1]),
       y = sprintf("PC2 (%.1f%%)", pct_before[2])) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8, color = "grey40"))

sex_col_map <- c(Male = "#2166AC", Female = "#D6604D")
p_before_sex <- make_pca_plot(
  pca_before_df, pct_before,
  color_var  = "sex", color_vals = sex_col_map,
  title      = "Before correction — Sex",
  subtitle   = "Color: Sex"
)

pca_after    <- prcomp(t(data_matrix), center = TRUE, scale. = FALSE)
pct_after    <- round(100 * pca_after$sdev^2 / sum(pca_after$sdev^2), 1)
pca_after_df <- as.data.frame(pca_after$x[, 1:2])
pca_after_df$sample    <- rownames(pca_after_df)
pca_after_df$LifeStage <- meta_final[pca_after_df$sample, "LifeStage"]
pca_after_df$batch     <- if ("batch" %in% colnames(meta_final))
  meta_final[pca_after_df$sample, "batch"] else "unknown"
pca_after_df$sex <- meta_final[pca_after_df$sample, "sex_label"]
pca_after_df$RIN <- as.numeric(meta_final[pca_after_df$sample, "RIN"])

p_after_stage <- make_pca_plot(
  pca_after_df, pct_after,
  color_var  = "LifeStage", color_vals = stage_col,
  title      = "After correction — LifeStage",
  subtitle   = "Stage structure retained"
)
p_after_batch <- make_pca_plot(
  pca_after_df, pct_after,
  color_var  = "batch", color_vals = batch_col,
  title      = "After correction — Batch",
  subtitle   = "Batch structure reduced"
)
p_after_RIN <- ggplot(pca_after_df, aes(x = PC1, y = PC2, color = RIN)) +
  geom_point(size = 2.2, alpha = 0.85) +
  scale_color_gradient2(low = "#4575B4", mid = "#FFFFBF", high = "#D73027",
                        midpoint = median(pca_after_df$RIN), name = "RIN") +
  labs(title    = "After correction — RIN",
       subtitle = "RIN gradient reduced",
       x = sprintf("PC1 (%.1f%%)", pct_after[1]),
       y = sprintf("PC2 (%.1f%%)", pct_after[2])) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8, color = "grey40"))
p_after_sex <- make_pca_plot(
  pca_after_df, pct_after,
  color_var  = "sex", color_vals = sex_col_map,
  title      = "After correction — Sex",
  subtitle   = "Sex structure reduced"
)

p_pca_combined <- (p_before_stage | p_after_stage) /
                  (p_before_batch  | p_after_batch)  /
                  (p_before_RIN    | p_after_RIN)    /
                  (p_before_sex    | p_after_sex)    +
  plot_annotation(
    title    = "PCA: Covariate Correction QC (Before vs After)",
    subtitle = "sex + RIN (ns df=2) regressed out via limma::removeBatchEffect (design protects LifeStage)",
    theme    = theme(plot.title    = element_text(face = "bold", size = 13),
                     plot.subtitle = element_text(size = 9, color = "grey40"))
  )

ggsave(file.path(output_dir, "00_pca_residual_check.pdf"),
       p_pca_combined, width = 14, height = 16)

sex_for_test <- ifelse(meta_final$sex_label == "Male", 1, 0)
rin_for_test    <- ns(meta_final$RIN, df = 2)
design_test  <- model.matrix(
  ~ LifeStage + sex_for_test + rin_for_test,
  data = data.frame(
    LifeStage    = meta_final$LifeStage,
    sex_for_test = sex_for_test,
    rin_for_test = rin_for_test
  )
)

fit_stage  <- lmFit(data_matrix, design_test)
fit_stage  <- eBayes(fit_stage)
stage_coef <- grep("^LifeStage", colnames(design_test), value = TRUE)
top_stage  <- topTable(fit_stage, coef = stage_coef,
                       number = Inf, sort.by = "none", adjust.method = "BH")

stage_sig_df <- data.frame(
  Gene      = rownames(top_stage),
  F_stat    = top_stage$F,
  P_val     = top_stage$P.Value,
  FDR       = top_stage$adj.P.Val,
  Sig       = top_stage$adj.P.Val < STAGE_FDR_CUTOFF,
  stringsAsFactors = FALSE
)
write.csv(stage_sig_df,
          file.path(output_dir, "00_stage_sig_fdr0.05.csv"),
          row.names = FALSE)

n_sig   <- sum(stage_sig_df$Sig, na.rm = TRUE)
n_total <- nrow(stage_sig_df)

p_sig_hist <- ggplot(stage_sig_df, aes(x = FDR)) +
  geom_histogram(bins = 50, fill = "#56B4E9", color = "white", alpha = 0.8) +
  geom_vline(xintercept = STAGE_FDR_CUTOFF, linetype = "dashed", color = "red", linewidth = 0.8) +
  annotate("text", x = STAGE_FDR_CUTOFF + 0.02, y = Inf, vjust = 2,
           label = sprintf("FDR=%.2f\n(%d dynamic)", STAGE_FDR_CUTOFF, n_sig),
           color = "red", hjust = 0, size = 3.5) +
  labs(title = "APA Stage Significance Distribution (All Events)",
       x = "BH-adjusted p-value (FDR)", y = "Count") +
  theme_bw()

stage_sig_df_flag <- stage_sig_df %>%
  mutate(Group = factor(
    ifelse(Sig, "Dynamic\n(FDR<cutoff)", "Stable\n(FDR≥cutoff)"),
    levels = c("Dynamic\n(FDR<cutoff)", "Stable\n(FDR≥cutoff)")
  ))

p_boxplot_filter <- ggplot(stage_sig_df_flag,
                           aes(x = Group, y = log2(F_stat + 1), fill = Group)) +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.25, linewidth = 0.5,
               width = 0.55, notch = TRUE) +
  scale_fill_manual(values = c("Dynamic\n(FDR<cutoff)" = "#56B4E9",
                                "Stable\n(FDR≥cutoff)"  = "#D3D3D3"),
                    guide = "none") +
  labs(title    = "Stage-effect Strength: Dynamic vs Stable APA Events",
       subtitle = sprintf(
         "limma F-statistic|Dynamic: %d  Stable: %d  (FDR cutoff=%.2f)",
         n_sig, n_total - n_sig, STAGE_FDR_CUTOFF),
       x = NULL,
       y = "log2(F-statistic + 1)") +
  theme_bw(base_size = 11) +
  theme(plot.subtitle = element_text(size = 8, color = "grey40"),
        axis.text.x   = element_text(size = 10))

ggsave(file.path(output_dir, "00_stage_significance_histogram.pdf"),
       p_sig_hist / p_boxplot_filter, width = 4, height = 8)

sig_genes          <- stage_sig_df$Gene[stage_sig_df$Sig]
stage_dynamic_flag <- setNames(
  ifelse(stage_sig_df$Sig, "Dynamic", "Stable"),
  stage_sig_df$Gene
)
cat(sprintf("  Stage-dynamic (FDR<%.2f): %d / %d APA events (%.1f%%) retained for clustering\n",
            STAGE_FDR_CUTOFF, length(sig_genes), nrow(stage_sig_df),
            100 * length(sig_genes) / nrow(stage_sig_df)))
stage_names <- c("Stage1", "Stage2", "Stage3")
data_scaled <- t(scale(t(sapply(stage_names, function(s) {
  s_samples <- rownames(meta_final)[meta_final$LifeStage == s]
  cat("Stage:", s, "→", length(s_samples), "samples\n")
  if (length(s_samples) > 1)  return(rowMeans(data_matrix[, s_samples], na.rm = TRUE))
  if (length(s_samples) == 1) return(data_matrix[, s_samples])
  return(rep(NA_real_, nrow(data_matrix)))
}))))
data_scaled <- data_scaled[complete.cases(data_scaled), ]
data_scaled <- data_scaled[rownames(data_scaled) %in% sig_genes, , drop = FALSE]
cat(sprintf("  Clustering input after FDR<%.2f filtering: %d APA events\n", STAGE_FDR_CUTOFF, nrow(data_scaled)))
all_gene_parts <- str_split_fixed(rownames(data_scaled), "\\|", 5)
gene_meta_df <- data.frame(
  Gene         = rownames(data_scaled),
  TranscriptID = all_gene_parts[, 1],
  GeneName     = all_gene_parts[, 2],
  stringsAsFactors = FALSE
)

cat("\n=== Data-driven K selection (Mfuzz Dmin + Silhouette) ===\n")
set.seed(2025)
eset <- new("ExpressionSet",
            exprs = data_scaled)

m_val <- mestimate(eset)
cat(sprintf("  Mfuzz fuzzifier m = %.4f\n", m_val))

dist_mat <- dist(t(data_scaled))  # sample-space dist not needed; use feature-space below
feat_dist <- dist(data_scaled)    # distance matrix for silhouette (genes × stages)

Dmin_vals <- numeric(K_MAX)
Sil_vals  <- numeric(K_MAX)
Dmin_vals[1] <- NA
Sil_vals[1]  <- NA
for (k in 2:K_MAX) {
  cl_obj <- mfuzz(eset, c = k, m = m_val)
  centers <- cl_obj$centers
  dists   <- dist(centers)
  Dmin_vals[k] <- min(as.numeric(dists))
  # Hard-assign each gene to its highest-membership cluster for silhouette
  hard_labels <- apply(cl_obj$membership, 1, which.max)
  if (length(unique(hard_labels)) >= 2) {
    sil_obj <- silhouette(hard_labels, feat_dist)
    Sil_vals[k] <- mean(sil_obj[, 3])
  } else {
    Sil_vals[k] <- NA_real_
  }
  cat(sprintf("  [Mfuzz]  K=%d: Dmin = %.6f  avg_sil = %.4f\n",
              k, Dmin_vals[k], Sil_vals[k]))
}
Dmin_df <- data.frame(K = 2:K_MAX, Dmin = Dmin_vals[2:K_MAX], AvgSil = Sil_vals[2:K_MAX])

# Dmin elbow: K immediately after the largest relative drop
k_seq    <- Dmin_df$K
d_seq    <- Dmin_df$Dmin
rel_drop <- diff(d_seq) / d_seq[-length(d_seq)]
i_max    <- which.min(rel_drop)
K_dmin   <- max(2L, as.integer(k_seq[i_max + 1]))

# Silhouette: K that maximises average silhouette width
K_sil    <- max(2L, as.integer(k_seq[which.max(Dmin_df$AvgSil)]))

# Conservative selection: take the smaller of the two to avoid over-segmentation
FINAL_K  <- min(K_dmin, K_sil)
cat(sprintf("\n>>> K_dmin = %d  K_sil = %d  →  FINAL_K = %d (conservative min) <<<\n",
            K_dmin, K_sil, FINAL_K))
cat(sprintf("    Dmin values (K=2..%d): %s\n", K_MAX,
            paste(round(Dmin_df$Dmin, 4), collapse = ", ")))
cat(sprintf("    AvgSil values:         %s\n",
            paste(round(Dmin_df$AvgSil, 4), collapse = ", ")))

set.seed(2025)
cl_final <- mfuzz(eset, c = FINAL_K, m = m_val)
cat(sprintf("  Mfuzz final run: K=%d, m=%.4f\n", FINAL_K, m_val))

p_dmin <- ggplot(Dmin_df, aes(x = K, y = Dmin)) +
  geom_line(color = "#E74C3C", linewidth = 1.2) +
  geom_point(size = 3.5, color = "#E74C3C") +
  geom_point(data = filter(Dmin_df, K == FINAL_K),
             size = 6, shape = 23, fill = "#E41A1C", color = "#E41A1C") +
  geom_vline(xintercept = FINAL_K, linetype = "dashed",
             color = "#2166AC", linewidth = 0.8, alpha = 0.7) +
  annotate("text", x = FINAL_K + 0.3, y = max(Dmin_df$Dmin) * 0.95,
           label = paste0("FINAL_K=", FINAL_K), color = "#2166AC", hjust = 0, size = 3.5) +
  geom_text(aes(label = round(Dmin, 4)), vjust = -1, size = 3, color = "grey30") +
  scale_x_continuous(breaks = 2:K_MAX) +
  labs(title    = "Mfuzz Dmin (minimum centroid distance)",
       subtitle = sprintf("Dmin elbow → K=%d; used together with silhouette (see below)", K_dmin),
       x = "K", y = "Dmin") +
  theme_bw() + theme(panel.grid.minor = element_blank())

p_sil <- ggplot(Dmin_df, aes(x = K, y = AvgSil)) +
  geom_line(color = "#2CA02C", linewidth = 1.2) +
  geom_point(size = 3.5, color = "#2CA02C") +
  geom_point(data = filter(Dmin_df, K == FINAL_K),
             size = 6, shape = 23, fill = "#1A7F1A", color = "#1A7F1A") +
  geom_vline(xintercept = FINAL_K, linetype = "dashed",
             color = "#2166AC", linewidth = 0.8, alpha = 0.7) +
  annotate("text", x = FINAL_K + 0.3, y = max(Dmin_df$AvgSil, na.rm = TRUE) * 0.95,
           label = paste0("FINAL_K=", FINAL_K), color = "#2166AC", hjust = 0, size = 3.5) +
  geom_text(aes(label = round(AvgSil, 4)), vjust = -1, size = 3, color = "grey30") +
  scale_x_continuous(breaks = 2:K_MAX) +
  labs(title    = "Average silhouette width",
       subtitle = sprintf("Silhouette max → K=%d; FINAL_K = min(K_dmin, K_sil) = %d", K_sil, FINAL_K),
       x = "K", y = "Avg silhouette") +
  theme_bw() + theme(panel.grid.minor = element_blank())

pdf(file.path(output_dir, "01_K_selection_evidence.pdf"), width = 8, height = 9)
print(p_dmin / p_sil)
dev.off()
cat("Saved 01_K_selection_evidence.pdf\n")

acore_list <- acore(eset, cl = cl_final, min.acore = 0.5)

n_genes_total <- nrow(data_scaled)
final_raw_labels <- integer(n_genes_total)
names(final_raw_labels) <- rownames(data_scaled)

mem_mat <- cl_final$membership
final_raw_labels <- apply(mem_mat, 1, which.max)
names(final_raw_labels) <- rownames(data_scaled)

cat(sprintf("  acore genes (mem>=0.5) per cluster:\n"))
for (cl_i in seq_len(FINAL_K)) {
  n_core <- nrow(acore_list[[cl_i]])
  n_total_cl <- sum(final_raw_labels == cl_i)
  cat(sprintf("    Cluster %d: %d core genes (mem>=0.5) / %d total (max-mem assigned)\n",
              cl_i, n_core, n_total_cl))
}

name_pattern <- function(centroid) {
  s1 <- centroid[1]; s2 <- centroid[2]; s3 <- centroid[3]
  if (s1 <= s2 && s2 >= s3 && s2 > s1) return("Peak")
  if (s1 <= s2 && s2 <= s3 && s3 > s1) return("Up")
  if (s1 >= s2 && s2 >= s3 && s1 > s3) return("Down")
  if (s1 >= s2 && s2 <= s3 && s2 < s1) return("Valley")
  return("Flat")
}

pattern_priority <- c("Peak" = 1, "Up" = 2, "Down" = 3, "Valley" = 4, "Flat" = 5)

centroid_info <- do.call(rbind, lapply(sort(unique(final_raw_labels)), function(cl) {
  genes <- names(final_raw_labels)[final_raw_labels == cl]
  ctr <- colMeans(data_scaled[genes, , drop = FALSE], na.rm = TRUE)
  pname <- name_pattern(ctr)
  data.frame(raw_cl = cl, pattern = pname, priority = pattern_priority[pname],
             n_genes = length(genes), s1 = ctr[1], s2 = ctr[2], s3 = ctr[3],
             stringsAsFactors = FALSE)
}))

centroid_info <- centroid_info[order(centroid_info$priority, -abs(centroid_info$s2 - (centroid_info$s1 + centroid_info$s3)/2)), ]

pattern_counts <- table(centroid_info$pattern)
dup_patterns <- names(pattern_counts[pattern_counts > 1])
if (length(dup_patterns) > 0) {
  for (dp in dup_patterns) {
    idx <- which(centroid_info$pattern == dp)
    centroid_info$pattern[idx] <- paste0(dp, "_", seq_along(idx))
  }
}

centroid_info$new_cl <- seq_len(nrow(centroid_info))
old_to_new <- setNames(centroid_info$new_cl, centroid_info$raw_cl)
pattern_labels <- setNames(centroid_info$pattern, as.character(centroid_info$new_cl))

final_labels <- as.integer(old_to_new[as.character(final_raw_labels)])
names(final_labels) <- names(final_raw_labels)

cat("\n=== Final Clustering (K =", FINAL_K, ") ===\n")
for (i in seq_len(nrow(centroid_info))) {
  cat(sprintf("  Cluster %d (%s): %d genes  centroid=[%.3f, %.3f, %.3f]\n",
              centroid_info$new_cl[i], centroid_info$pattern[i], centroid_info$n_genes[i],
              centroid_info$s1[i], centroid_info$s2[i], centroid_info$s3[i]))
}

classify_trajectory <- function(zvec) {
  z1 <- zvec[1]; z2 <- zvec[2]; z3 <- zvec[3]
  pname <- name_pattern(zvec)
  paste0("Pattern_", pattern_priority[pname], "_", pname)
}

ss_final <- silhouette(final_labels, dist(data_scaled))
pdf(file.path(output_dir, "02_silhouette_final.pdf"), width = 10, height = 8)
plot(ss_final, col = get_cluster_colors(FINAL_K)[sort(unique(final_labels))],
     main = sprintf("Silhouette Plot (K=%d, avg=%.3f)", FINAL_K, mean(ss_final[, 3])),
     border = NA)
dev.off()
cat("Saved 02_silhouette_final.pdf\n")

centroid_plot_df <- do.call(rbind, lapply(seq_len(FINAL_K), function(cl) {
  genes <- names(final_labels)[final_labels == cl]
  ctr <- colMeans(data_scaled[genes, , drop = FALSE], na.rm = TRUE)
  data.frame(
    Cluster = factor(cl),
    Stage   = factor(stage_names, levels = stage_names),
    value   = ctr,
    stringsAsFactors = FALSE
  )
}))
centroid_plot_df$ClusterLabel <- paste0("C", centroid_plot_df$Cluster, " ",
                                        pattern_labels[as.character(centroid_plot_df$Cluster)],
                                        " (n=", centroid_info$n_genes[as.integer(as.character(centroid_plot_df$Cluster))], ")")

color_map <- get_cluster_colors(FINAL_K)

p_centroid <- ggplot(centroid_plot_df, aes(x = as.numeric(Stage), y = value,
                                           group = Cluster, color = Cluster)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) +
  scale_color_manual(values = color_map, labels = unique(centroid_plot_df$ClusterLabel)) +
  scale_x_continuous(breaks = seq_along(stage_names), labels = stage_names) +
  labs(title = paste0("Cluster Centroid Trends (K=", FINAL_K, ", data-driven)"),
       subtitle = "Pattern names auto-assigned by centroid shape",
       x = "Stage", y = "Z-score", color = "Cluster") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position   = "right",
        panel.grid.minor  = element_blank())

ggsave(file.path(output_dir, "03_centroid_trends.pdf"), p_centroid, width = 10, height = 6)
cat("Saved 03_centroid_trends.pdf\n")

nrow_plot <- ceiling(sqrt(FINAL_K))
ncol_plot <- ceiling(FINAL_K / nrow_plot)
pdf(file.path(output_dir, "03b_mfuzz_plot2.pdf"),
    width = ncol_plot * 4, height = nrow_plot * 4)
mfuzz_colo <- colorRampPalette(c("#CCCCCC", "#4DBBD5"))(100)

.mfuzz_plot2_patched <- Mfuzz::mfuzz.plot2
body(.mfuzz_plot2_patched)[[6]] <- quote(
  if (missing(colo)) {
    colo <- c("#FF0000","#FF1800","#FF3000","#FF4800","#FF6000","#FF7800",
              "#FF8F00","#FFA700","#FFBF00","#FFD700","#FFEF00","#F7FF00",
              "#DFFF00","#C7FF00","#AFFF00","#97FF00","#80FF00","#68FF00",
              "#50FF00","#38FF00","#20FF00","#08FF00","#00FF10","#00FF28",
              "#00FF40","#00FF58","#00FF70","#00FF87","#00FF9F","#00FFB7",
              "#00FFCF","#00FFE7","#00FFFF","#00E7FF","#00CFFF","#00B7FF",
              "#009FFF","#0087FF","#0070FF","#0058FF","#0040FF","#0028FF",
              "#0010FF","#0800FF","#2000FF","#3800FF","#5000FF","#6800FF",
              "#8000FF","#9700FF","#AF00FF","#C700FF","#DF00FF","#F700FF",
              "#FF00EF","#FF00D7","#FF00BF","#FF00A7","#FF008F","#FF0078",
              "#FF0060","#FF0048","#FF0030","#FF0018")
  }
)

.mfuzz_plot2_patched(eset, cl = cl_final,
            mfrow   = c(nrow_plot, ncol_plot),
            centre  = TRUE,
            centre.col = "black",
            centre.lwd = 2,
            colo    = mfuzz_colo,
            min.mem = 0.7,
            time.labels = stage_names,
            x11     = FALSE)
dev.off()
cat("Saved 03b_mfuzz_plot2.pdf\n")
cat("Saved 03_centroid_trends.pdf\n")


if (any(grepl("\\|", rownames(data_matrix)))) {
  all_genes_symbol <- str_split(rownames(data_matrix), "\\|", simplify = TRUE)[, 2]
} else {
  all_genes_symbol <- rownames(data_matrix)
  message("No pipe-delimited rownames detected; using rownames as symbols")
}
universe_genes <- unique(all_genes_symbol[all_genes_symbol != "" & !is.na(all_genes_symbol)])
universe_entrez <- tryCatch({
  bitr(universe_genes, fromType = "SYMBOL", toType = "ENTREZID",
       OrgDb = org.Hs.eg.db)$ENTREZID
}, error = function(e) { warning("Failed to map universe symbols to Entrez IDs"); NULL })
cat(sprintf("Background universe from full data_matrix: %d unique gene symbols, %d Entrez IDs\n",
            length(universe_genes), length(universe_entrez)))


extract_symbols <- function(gene_ids) {
  if (any(grepl("\\|", gene_ids))) {
    syms <- str_split(gene_ids, "\\|", simplify = TRUE)[, 2]
    return(unique(syms[syms != "" & !is.na(syms)]))
  } else {
    return(unique(gene_ids))
  }
}


enrich_percluster_GO <- function(gene_list,          # named list: cluster -> SYMBOL vector
                                  universe,
                                  orgdb  = org.Hs.eg.db,
                                  ont    = "BP",
                                  padj_method = "BH") {
  results <- list()
  for (cl in names(gene_list)) {
    syms <- gene_list[[cl]]
    if (length(syms) < 5) next
     res <- tryCatch(
       enrichGO(gene          = syms,
                universe      = universe,
                OrgDb         = orgdb,
                keyType       = "SYMBOL",
                ont           = ont,
                pAdjustMethod = padj_method,
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2,
                readable      = FALSE),
       error = function(e) NULL
     )
    if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
      df <- as.data.frame(res)
      df$Cluster <- cl
      results[[cl]] <- df
    }
  }
  if (length(results) == 0) return(NULL)
  bind_rows(results)
}

enrich_percluster_KEGG <- function(entrez_list,      # named list: cluster -> ENTREZID vector
                                    universe_entrez,
                                    organism    = "hsa",
                                    padj_method = "BH") {
  results <- list()
  for (cl in names(entrez_list)) {
    ids <- entrez_list[[cl]]
    if (length(ids) < 5) next
     res <- tryCatch(
       enrichKEGG(gene          = ids,
                  universe      = universe_entrez,
                  organism      = organism,
                  keyType       = "kegg",
                  pAdjustMethod = padj_method,
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2),
       error = function(e) NULL
     )
    if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
      df <- as.data.frame(res)
      df$Cluster <- cl
      results[[cl]] <- df
    }
  }
  if (length(results) == 0) return(NULL)
  bind_rows(results)
}

make_compare_cluster_result <- function(df) {
  obj <- methods::new("compareClusterResult",
                      compareClusterResult = df,
                      fun                  = "enrichGO",
                      .call                = call("enrichGO"))
  obj
}

select_terms_per_cluster <- function(res_df, n_per_cluster = 5, min_count = 3) {
  res_df <- res_df %>%
    mutate(
      GeneRatio_numeric = sapply(strsplit(as.character(GeneRatio), "/"),
                                 function(x) as.numeric(x[1]) / as.numeric(x[2])),
      Count_numeric = as.numeric(Count)
    ) %>%
    filter(Count_numeric >= min_count)

  if (nrow(res_df) == 0) return(list(plot_data = data.frame(), cluster_strategies = data.frame()))

  cluster_plot_data  <- list()
  cluster_strategies <- data.frame()

  for (cl in unique(res_df$Cluster)) {
    cl_data <- filter(res_df, Cluster == cl)
    sig_05  <- filter(cl_data, p.adjust < 0.05)
    sig_02  <- filter(cl_data, p.adjust < 0.20)

    if (nrow(sig_05) >= n_per_cluster) {
      selected <- slice_head(arrange(sig_05, p.adjust), n = n_per_cluster)
      strategy <- "Significant (p.adj<0.05)"
    } else if (nrow(sig_02) >= n_per_cluster) {
      selected <- slice_head(arrange(sig_02, p.adjust), n = n_per_cluster)
      strategy <- "Marginally sig (p.adj<0.2)"
    } else {
      selected <- slice_head(arrange(cl_data, desc(GeneRatio_numeric), desc(Count_numeric)),
                             n = n_per_cluster)
      strategy <- sprintf("Exploratory (GeneRatio, %d sig)", sum(cl_data$p.adjust < 0.05, na.rm = TRUE))
    }

    cluster_plot_data[[as.character(cl)]] <- selected
    cluster_strategies <- rbind(cluster_strategies, data.frame(
      Cluster        = cl,
      Strategy       = strategy,
      Total_Terms    = nrow(cl_data),
      Selected_Terms = nrow(selected),
      N_Sig_0.05     = sum(cl_data$p.adjust < 0.05, na.rm = TRUE),
      N_Sig_0.2      = sum(cl_data$p.adjust < 0.20, na.rm = TRUE)
    ))
  }
  list(plot_data = bind_rows(cluster_plot_data), cluster_strategies = cluster_strategies)
}

filter_core_genes <- function(data_scaled, cluster_labels,
                              cor_threshold = CORE_GENE_COR_THRESHOLD,
                              min_size = CORE_GENE_MIN_SIZE) {
  core_genes_list <- list()
  filter_stats    <- data.frame()
  all_gene_cors   <- data.frame()

  for (cl in sort(unique(cluster_labels))) {
    cluster_genes  <- names(cluster_labels)[cluster_labels == cl]
    cluster_expr   <- data_scaled[cluster_genes, , drop = FALSE]
    cluster_center <- colMeans(cluster_expr, na.rm = TRUE)
    gene_cors <- apply(cluster_expr, 1, function(g)
      cor(g, cluster_center, use = "pairwise.complete.obs"))

    all_gene_cors <- rbind(all_gene_cors, data.frame(
      Gene = cluster_genes, Cluster = cl, Correlation = gene_cors, stringsAsFactors = FALSE))

    core_genes <- cluster_genes[gene_cors >= cor_threshold]
    if (length(core_genes) < min_size) {
      n_take     <- min(min_size, length(cluster_genes))
      core_genes <- cluster_genes[order(gene_cors, decreasing = TRUE)[seq_len(n_take)]]
      warning(sprintf("Cluster %d: only %d genes met threshold %.2f; using top %d",
                      cl, sum(gene_cors >= cor_threshold), cor_threshold, n_take))
    }

    core_genes_list[[as.character(cl)]] <- core_genes
    filter_stats <- rbind(filter_stats, data.frame(
      Cluster          = cl,
      Total_Genes      = length(cluster_genes),
      Core_Genes       = length(core_genes),
      Filter_Rate      = round(length(core_genes) / length(cluster_genes), 3),
      Min_Correlation  = round(min(gene_cors[core_genes]), 3),
      Mean_Correlation = round(mean(gene_cors[core_genes]), 3)
    ))
  }
  list(core_genes = core_genes_list, filter_stats = filter_stats, gene_correlations = all_gene_cors)
}

cat(sprintf("\n=== Core gene filtering K = %d ===\n", FINAL_K))
final_core_result <- filter_core_genes(data_scaled, final_labels)
cat("\n=== Core Gene Stats ===\n")
print(final_core_result$filter_stats)
write.csv(final_core_result$filter_stats,
          file.path(output_dir, "05_core_gene_filter_summary.csv"),
          row.names = FALSE)

cat(sprintf("\n=== GSEA K = %d (per-cluster Mfuzz correlation ranking) ===\n", FINAL_K))

build_gsea_rank <- function(cluster_id, data_sc, labels) {
  ctr <- colMeans(data_sc[names(labels)[labels == cluster_id], , drop = FALSE], na.rm = TRUE)
  cors <- apply(data_sc, 1, function(g) cor(g, ctr, use = "pairwise.complete.obs"))
  syms <- if (any(grepl("\\|", names(cors)))) {
    str_split(names(cors), "\\|", simplify = TRUE)[, 2]
  } else {
    names(cors)
  }
  df <- data.frame(Symbol = syms, Correlation = cors, stringsAsFactors = FALSE) %>%
    filter(Symbol != "" & !is.na(Symbol) & !is.na(Correlation))
  df <- df %>%
    group_by(Symbol) %>%
    slice_max(order_by = abs(Correlation), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(Correlation))
  rank_vec <- setNames(df$Correlation, df$Symbol)
  rank_vec
}

cluster_label_names <- setNames(
  c("Peak", "Up", "Down", "Valley", "Flat")[seq_len(FINAL_K)],
  as.character(seq_len(FINAL_K))
)

gsea_go_results   <- list()
gsea_kegg_results <- list()

set.seed(42)
for (cl in sort(unique(final_labels))) {
  cl_chr  <- as.character(cl)
  cl_name <- cluster_label_names[cl_chr]
  cat(sprintf("  Cluster %d (%s): building ranked list...\n", cl, cl_name))

  rank_vec <- build_gsea_rank(cl, data_scaled, final_labels)
  cat(sprintf("    Ranked genes: %d\n", length(rank_vec)))

  # ---------- GSEA GO BP ----------
  gsea_go <- tryCatch({
    gseGO(
      geneList     = rank_vec,
      OrgDb        = org.Hs.eg.db,
      keyType      = "SYMBOL",
      ont          = "BP",
      minGSSize    = 15,
      maxGSSize    = 500,
      pvalueCutoff = 0.2,
      pAdjustMethod = "BH",
      verbose      = FALSE,
      eps          = 0
    )
  }, error = function(e) { cat("  gseGO failed:", e$message, "\n"); NULL })

  if (!is.null(gsea_go) && nrow(gsea_go@result) > 0) {
    gsea_go_simp <- tryCatch(
      simplify(gsea_go, cutoff = 0.7, by = "p.adjust", select_fun = min),
      error = function(e) gsea_go
    )
    go_df <- as.data.frame(gsea_go_simp@result) %>%
      mutate(Cluster = cl, ClusterName = cl_name)
    gsea_go_results[[cl_chr]] <- go_df
    cat(sprintf("    GO BP GSEA: %d terms (after simplify)\n", nrow(go_df)))
  } else {
    cat("    GO BP GSEA: no significant results (p.adj < 0.05)\n")
  }

  entrez_map <- tryCatch(
    bitr(names(rank_vec), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
    error = function(e) NULL
  )
  if (!is.null(entrez_map) && nrow(entrez_map) > 0) {
    rank_entrez <- rank_vec[entrez_map$SYMBOL]
    names(rank_entrez) <- entrez_map$ENTREZID
    rank_entrez <- rank_entrez[!duplicated(names(rank_entrez))]
    rank_entrez <- sort(rank_entrez, decreasing = TRUE)

    gsea_kegg <- tryCatch({
      gseKEGG(
        geneList     = rank_entrez,
        organism     = "hsa",
        minGSSize    = 15,
        maxGSSize    = 500,
        pvalueCutoff = 0.2,
        pAdjustMethod = "BH",
        verbose      = FALSE,
        eps          = 0,
        use_internal_data = FALSE
      )
    }, error = function(e) { cat("  gseKEGG failed:", e$message, "\n"); NULL })

    if (!is.null(gsea_kegg) && nrow(gsea_kegg@result) > 0) {
      kegg_df <- as.data.frame(gsea_kegg@result) %>%
        mutate(Cluster = cl, ClusterName = cl_name)
      gsea_kegg_results[[cl_chr]] <- kegg_df
      cat(sprintf("    KEGG GSEA: %d terms\n", nrow(kegg_df)))
    } else {
      cat("    KEGG GSEA: no significant results (p.adj < 0.05)\n")
    }
  }
}

if (length(gsea_go_results) > 0) {
  gsea_go_all <- bind_rows(gsea_go_results)
  write.csv(gsea_go_all,
            file.path(output_dir, sprintf("05b_GSEA_GO_BP_K%d_all.csv", FINAL_K)),
            row.names = FALSE)
  cat(sprintf("Saved GSEA GO BP results: %d rows\n", nrow(gsea_go_all)))
} else {
  cat("GSEA GO BP: no significant results across clusters\n")
  gsea_go_all <- NULL
}

if (length(gsea_kegg_results) > 0) {
  gsea_kegg_all <- bind_rows(gsea_kegg_results)
  write.csv(gsea_kegg_all,
            file.path(output_dir, sprintf("05b_GSEA_KEGG_K%d_all.csv", FINAL_K)),
            row.names = FALSE)
  cat(sprintf("Saved GSEA KEGG results: %d rows\n", nrow(gsea_kegg_all)))
} else {
  cat("GSEA KEGG: no significant results across clusters\n")
  gsea_kegg_all <- NULL
}

make_gsea_dotplot <- function(gsea_df, top_n = 5, title = "GSEA Results",
                               sig_cutoff = 0.05) {
  if (is.null(gsea_df) || nrow(gsea_df) == 0) return(NULL)

  pattern_names <- c("1" = "Peak", "2" = "Up", "3" = "Down",
                      "4" = "Valley", "5" = "Flat")

  top_df <- gsea_df %>%
    mutate(ClusterName = pattern_names[as.character(Cluster)]) %>%
    group_by(Cluster) %>%
    slice_min(order_by = p.adjust, n = top_n, with_ties = FALSE) %>%
    ungroup()

  shown_terms <- unique(top_df$Description)

  all_clusters <- sort(unique(gsea_df$Cluster))
  full_grid <- expand.grid(
    Cluster     = all_clusters,
    Description = shown_terms,
    stringsAsFactors = FALSE
  )
  plot_df <- full_grid %>%
    left_join(gsea_df %>% select(Cluster, Description, NES, p.adjust, setSize),
              by = c("Cluster", "Description")) %>%
    mutate(
      ClusterName  = pattern_names[as.character(Cluster)],
      ClusterLabel = paste0("Cluster ", Cluster, "\n(", ClusterName, ")"),
      neg_log10p   = ifelse(is.na(p.adjust), NA_real_, -log10(p.adjust)),
      Significant  = !is.na(p.adjust) & p.adjust < sig_cutoff,
      Description  = str_wrap(Description, width = 45)
    )

  term_order <- plot_df %>%
    filter(!is.na(NES)) %>%
    group_by(Description) %>%
    summarise(n_cl = n(), mean_abs_nes = mean(abs(NES), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(n_cl), desc(mean_abs_nes)) %>%
    pull(Description)
  plot_df$Description <- factor(plot_df$Description, levels = rev(term_order))

  nes_lim <- max(abs(plot_df$NES), na.rm = TRUE)
  nes_lim <- ceiling(nes_lim * 10) / 10

  p <- ggplot(plot_df, aes(x = ClusterLabel, y = Description)) +
    geom_point(data = filter(plot_df, is.na(NES)),
               size = 3, color = "grey88", shape = 16) +
    geom_point(data = filter(plot_df, !is.na(NES)),
               aes(size = neg_log10p, color = NES), alpha = 0.92) +
    geom_point(data = filter(plot_df, Significant),
               aes(size = neg_log10p), shape = 21,
               fill = NA, color = "black", stroke = 0.7) +
    scale_color_gradient2(
      low      = "#3B82C4",
      mid      = "white",
      high     = "#C0392B",
      midpoint = 0,
      limits   = c(-nes_lim, nes_lim),
      name     = "NES"
    ) +
    scale_size_continuous(
      range = c(2.5, 9),
      name  = expression(-log[10](p.adj)),
      breaks = function(x) pretty(x, n = 4)
    ) +
    labs(title    = title,
         subtitle = sprintf("Circle size = -log10(p.adjust) | Black border = p.adj < %.2f | Blue = suppressed, Red = activated",
                            sig_cutoff),
         x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x        = element_text(size = 9, face = "bold"),
      axis.text.y        = element_text(size = 8),
      plot.title         = element_text(face = "bold", size = 11),
      plot.subtitle      = element_text(size = 7.5, color = "grey40"),
      panel.grid.major   = element_line(colour = "grey92", linewidth = 0.3),
      legend.key.height  = unit(0.5, "cm"),
      legend.title       = element_text(size = 8),
      legend.text        = element_text(size = 7.5)
    )
  p
}

if (!is.null(gsea_go_all)) {
  p_gsea_go <- make_gsea_dotplot(
    gsea_go_all, top_n = 5, sig_cutoff = 0.05,
    title = sprintf("K=%d GSEA GO BP (Mfuzz APA correlation, top 5 per cluster)", FINAL_K)
  )
  if (!is.null(p_gsea_go)) {
    n_terms_go <- length(unique(gsea_go_all$Description[
      gsea_go_all$Description %in% (gsea_go_all %>%
        group_by(Cluster) %>% slice_min(p.adjust, n=5, with_ties=FALSE) %>%
        ungroup() %>% pull(Description))
    ]))
    h_go <- max(4, n_terms_go * 0.45 + 2)
    ggsave(file.path(output_dir, sprintf("05b_GSEA_GO_BP_K%d.pdf", FINAL_K)),
           p_gsea_go, width = 3 * FINAL_K + 2, height = h_go)
    cat("Saved GSEA GO BP dot plot\n")
  }
}

if (!is.null(gsea_kegg_all)) {
  p_gsea_kegg <- make_gsea_dotplot(
    gsea_kegg_all, top_n = 5, sig_cutoff = 0.05,
    title = sprintf("K=%d GSEA KEGG (Mfuzz APA correlation, top 5 per cluster)", FINAL_K)
  )
  if (!is.null(p_gsea_kegg)) {
    n_terms_kegg <- length(unique(gsea_kegg_all$Description[
      gsea_kegg_all$Description %in% (gsea_kegg_all %>%
        group_by(Cluster) %>% slice_min(p.adjust, n=5, with_ties=FALSE) %>%
        ungroup() %>% pull(Description))
    ]))
    h_kegg <- max(4, n_terms_kegg * 0.45 + 2)
    ggsave(file.path(output_dir, sprintf("05b_GSEA_KEGG_K%d.pdf", FINAL_K)),
           p_gsea_kegg, width = 3 * FINAL_K + 2, height = h_kegg)
    cat("Saved GSEA KEGG dot plot\n")
  }
}

cat("GSEA analysis completed\n")

final_assignments_full <- final_core_result$gene_correlations %>%
  rename(Cluster_Correlation = Correlation) %>%
  mutate(Is_Core_Gene = Gene %in% unlist(final_core_result$core_genes))
write.csv(final_assignments_full,
          file.path(output_dir, "06_final_cluster_assignments_with_correlation.csv"),
          row.names = FALSE)

final_core_assignments <- filter(final_assignments_full, Is_Core_Gene) %>%
  select(Gene, Cluster, Cluster_Correlation)
write.csv(final_core_assignments,
          file.path(output_dir, "06_final_CORE_cluster_assignments.csv"),
          row.names = FALSE)

traj_shape_core <- apply(data_scaled[unlist(final_core_result$core_genes), ], 1,
                         classify_trajectory)
ordered_core_genes <- final_core_assignments %>%
  mutate(
    slope      = data_scaled[Gene, ncol(data_scaled)] - data_scaled[Gene, 1],
    traj_shape = traj_shape_core[Gene]
  ) %>%
  arrange(Cluster, traj_shape, desc(slope)) %>%
  pull(Gene)

plot_matrix_core <- data_scaled[ordered_core_genes, ]
core_labels      <- final_labels[ordered_core_genes]
traj_shape_ordered <- traj_shape_core[ordered_core_genes]

all_possible_shapes <- c(
  "Pattern_1_Peak"   = "#D55E00",
  "Pattern_2_Up"     = "#E69F00",
  "Pattern_3_Down"   = "#4DBBD5",
  "Pattern_4_Valley" = "#0072B2",
  "Pattern_5_Flat"   = "#999999"
)
shapes_present <- unique(traj_shape_ordered)
shape_col <- all_possible_shapes[names(all_possible_shapes) %in% shapes_present]


row_ha_core <- rowAnnotation(
  Cluster    = as.factor(core_labels),
  TrajectoryShape = traj_shape_ordered,
  col        = list(Cluster         = cluster_palette,
                    TrajectoryShape = shape_col),
  show_annotation_name = TRUE,
  annotation_name_gp   = gpar(fontsize = 8)
)

pdf(file.path(output_dir, "06_heatmap_core_genes.pdf"), width = 15, height = 16)
draw(Heatmap(plot_matrix_core, name = "Z-score",
             cluster_rows = FALSE, cluster_columns = FALSE,
             left_annotation = row_ha_core,
             show_row_names = FALSE, show_column_names = TRUE,
             column_names_rot = 0, column_names_centered = TRUE,
             col = colorRamp2(c(-2, 0, 2), c("#7EC0EE", "white", "#FFD700")),
             row_split = core_labels, row_gap = unit(2, "mm"), border = TRUE,
             column_title = paste0("Core Gene Expression (K=", FINAL_K,
                                   ", cor >= ", CORE_GENE_COR_THRESHOLD, ")"),
             row_title = "Core Genes"))
dev.off()


centroid_lookup <- setNames(
  lapply(sort(unique(final_labels)), function(cl)
    colMeans(data_scaled[final_labels == cl, , drop = FALSE])),
  as.character(sort(unique(final_labels)))
)

N_BG_GENES <- 20

representative_genes_top <- final_core_assignments %>%
  mutate(Cluster_chr = as.character(Cluster)) %>%
  rowwise() %>%
  mutate(Distance = sqrt(sum((data_scaled[Gene, ] - centroid_lookup[[Cluster_chr]])^2))) %>%
  ungroup() %>%
  group_by(Cluster) %>%
  arrange(Distance) %>%
  slice_head(n = N_BG_GENES) %>%
  ungroup() %>%
  left_join(gene_meta_df, by = "Gene")

representative_genes <- representative_genes_top %>%
  group_by(Cluster) %>%
  slice_head(n = 1) %>%
  ungroup()

write.csv(representative_genes,
          file.path(output_dir, "06_representative_core_genes.csv"), row.names = FALSE)

rep_bg_df <- data_scaled[representative_genes_top$Gene, , drop = FALSE] %>%
  as.data.frame() %>% rownames_to_column("Gene") %>%
  pivot_longer(cols = all_of(stage_names), names_to = "Stage", values_to = "Zscore") %>%
  left_join(representative_genes_top %>% select(Gene, Cluster), by = "Gene") %>%
  mutate(
    Cluster  = factor(Cluster, levels = sort(unique(Cluster))),
    Stage    = factor(Stage, levels = stage_names),
    StageNum = as.numeric(Stage)
  )

centroid_df <- do.call(rbind, lapply(sort(unique(final_labels)), function(cl) {
  data.frame(
    Stage    = factor(stage_names, levels = stage_names),
    StageNum = seq_along(stage_names),
    Zscore   = centroid_lookup[[as.character(cl)]],
    Cluster  = factor(cl, levels = sort(unique(final_labels)))
  )
}))

centroid_summary <- centroid_df %>%
  group_by(Cluster) %>%
  summarise(
    S1 = Zscore[StageNum == 1],
    S2 = Zscore[StageNum == 2],
    S3 = Zscore[StageNum == max(StageNum)],
    .groups = "drop"
  ) %>%
  mutate(Pattern = pattern_labels[as.character(Cluster)])

centroid_df <- centroid_df %>%
  left_join(centroid_summary %>% select(Cluster, Pattern), by = "Cluster") %>%
  mutate(ClusterLabel = paste0("Cluster ", Cluster, "\n(", Pattern, ", n=",
    table(final_core_assignments$Cluster)[as.character(Cluster)], ")"))

rep_bg_df <- rep_bg_df %>%
  left_join(centroid_df %>% distinct(Cluster, ClusterLabel), by = "Cluster")

rep_gene_label_df <- data_scaled[representative_genes$Gene, , drop = FALSE] %>%
  as.data.frame() %>% rownames_to_column("Gene") %>%
  pivot_longer(cols = all_of(stage_names), names_to = "Stage", values_to = "Zscore") %>%
  left_join(representative_genes, by = "Gene") %>%
  mutate(
    Cluster      = factor(Cluster, levels = sort(unique(Cluster))),
    GeneName     = ifelse(is.na(GeneName)     | GeneName     == "", TranscriptID, GeneName),
    TranscriptID = ifelse(is.na(TranscriptID) | TranscriptID == "", Gene, TranscriptID),
    Label        = str_wrap(paste(TranscriptID, GeneName, sep = "|"), width = 18),
    Stage        = factor(Stage, levels = stage_names),
    StageNum     = as.numeric(Stage)
  ) %>%
  left_join(centroid_df %>% distinct(Cluster, ClusterLabel), by = "Cluster")

write.csv(rep_gene_label_df,
          file.path(output_dir, "06_representative_core_gene_trends.csv"), row.names = FALSE)

stage_shading <- tibble(
  xmin = seq_along(stage_names) - 0.5,
  xmax = seq_along(stage_names) + 0.5,
  fill = rep(c("#F7F9FB", "#FFFFFF"), length.out = length(stage_names))
)
label_positions <- filter(rep_gene_label_df, Stage == tail(stage_names, 1))
cluster_count   <- length(unique(centroid_df$Cluster))

p_rep_gene <- ggplot() +
  geom_rect(data = stage_shading, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
            alpha = 0.6, color = NA) +
  geom_hline(yintercept = 0, color = "#CBD5E1", linewidth = 0.4) +
  geom_line(data = rep_bg_df,
            aes(x = StageNum, y = Zscore, group = Gene),
            color = "#C0C8D8", linewidth = 0.4, alpha = 0.6) +
  geom_line(data = centroid_df,
            aes(x = StageNum, y = Zscore, color = Cluster, group = Cluster),
            linewidth = 2.0) +
  geom_point(data = centroid_df,
             aes(x = StageNum, y = Zscore, color = Cluster, group = Cluster),
             size = 4, shape = 21, fill = "white", stroke = 1.5) +
  geom_text_repel(data = label_positions,
                  aes(x = StageNum, y = Zscore, label = Label),
                  nudge_x = 0.35, nudge_y = -0.3, direction = "both", hjust = 0,
                  size = 3.0, color = "#64748B",
                  segment.color = "#94A3B8", segment.size = 0.3,
                  min.segment.length = 0, max.overlaps = Inf, box.padding = 0.4) +
  facet_wrap(~ ClusterLabel, ncol = cluster_count) +
  scale_color_manual(values = cluster_palette, drop = FALSE, guide = "none") +
  scale_fill_identity() +
  scale_x_continuous(breaks = seq_along(stage_names), labels = stage_names,
                     expand = expansion(mult = c(0.02, 0.25))) +
  labs(title    = paste0("Cluster Centroid Trajectories (K=", FINAL_K, ", data-driven)"),
       subtitle = paste0(
         "Bold colored line = cluster mean (centroid); ",
         "Gray lines = top ", N_BG_GENES, " core genes closest to centroid\n",
         "Core genes: cor >= ", CORE_GENE_COR_THRESHOLD,
         " | Stage filter: FDR < ", STAGE_FDR_CUTOFF),
       x = "Stage", y = "Z-score (mean across samples)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor    = element_blank(),
        panel.grid.major.x  = element_blank(),
        panel.grid.major.y  = element_line(color = "#E2E8F0", linewidth = 0.3),
        strip.text          = element_text(face = "bold"),
        plot.title          = element_text(face = "bold", size = 14),
        panel.spacing       = unit(1.2, "lines"),
        plot.margin         = margin(20, 50, 20, 20)) +
  coord_cartesian(clip = "off")

ggsave(file.path(output_dir, "06_representative_core_gene_trends.pdf"),
       p_rep_gene, width = 14, height = 5)

rep_sample_df <- lapply(seq_len(nrow(representative_genes)), function(i) {
  gid   <- representative_genes$Gene[i]
  cl    <- representative_genes$Cluster[i]
  label <- if (!is.na(representative_genes$GeneName[i]) && representative_genes$GeneName[i] != "")
    paste0(representative_genes$GeneName[i], "\n(Cluster ", cl, ")")
  else paste0(representative_genes$TranscriptID[i], "\n(Cluster ", cl, ")")
  expr_vals <- data_matrix[gid, , drop = TRUE]
  data.frame(
    Sample    = names(expr_vals),
    Expr      = as.numeric(expr_vals),
    LifeStage = meta_final[names(expr_vals), "LifeStage"],
    Gene      = gid,
    Label     = label,
    Cluster   = as.character(cl),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows() %>%
  mutate(LifeStage = factor(LifeStage, levels = stage_names),
         StageNum  = as.numeric(LifeStage),
         Label     = factor(Label, levels = unique(Label)))

rep_stage_mean_df <- rep_sample_df %>%
  group_by(Label, Cluster, LifeStage, StageNum) %>%
  summarise(Mean = mean(Expr, na.rm = TRUE), .groups = "drop")

n_rep <- length(unique(rep_sample_df$Label))

p_rep_sample <- ggplot(rep_sample_df,
                       aes(x = StageNum, y = Expr, color = LifeStage)) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.5, show.legend = TRUE) +
  geom_line(data = rep_stage_mean_df,
            aes(x = StageNum, y = Mean, group = Label),
            color = "black", linewidth = 1.1, inherit.aes = FALSE) +
  geom_point(data = rep_stage_mean_df,
             aes(x = StageNum, y = Mean, fill = LifeStage),
             shape = 21, size = 3, color = "black", stroke = 0.5,
             inherit.aes = FALSE) +
  facet_wrap(~ Label, ncol = 4, scales = "free_y") +
  scale_color_manual(values = stage_col, name = "LifeStage") +
  scale_fill_manual(values  = stage_col, name = "LifeStage") +
  scale_x_continuous(breaks = seq_along(stage_names), labels = stage_names) +
  labs(title    = paste0("Representative Gene: Sample-level Expression (K=", FINAL_K, ")"),
       subtitle = "Dots = individual samples; Black line = stage mean",
       x = "Stage", y = "Corrected Expression") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor   = element_blank(),
        strip.text         = element_text(face = "bold", size = 9),
        plot.title         = element_text(face = "bold", size = 13),
        legend.position    = "bottom",
        panel.spacing      = unit(1, "lines"))

ggsave(file.path(output_dir, "06_representative_sample_scatter.pdf"),
       p_rep_sample, width = 16, height = 4 * ceiling(n_rep / 4))
cat("Saved 06_representative_sample_scatter.pdf\n")

cat("\n=== sample heatmap ===\n")

all_core_genes_vec <- unlist(final_core_result$core_genes)
sample_expr        <- data_matrix[all_core_genes_vec, , drop = FALSE]

traj_shape_sample <- apply(data_scaled[all_core_genes_vec, ], 1, classify_trajectory)
slope_vec         <- data_scaled[all_core_genes_vec, ncol(data_scaled)] -
                     data_scaled[all_core_genes_vec, 1]
order_df_sample <- data.frame(
  Gene       = all_core_genes_vec,
  Cluster    = final_labels[all_core_genes_vec],
  traj_shape = traj_shape_sample[all_core_genes_vec],
  slope      = slope_vec,
  stringsAsFactors = FALSE
) %>% arrange(Cluster, traj_shape, desc(slope))

ordered_core_sample <- order_df_sample$Gene
sample_expr_ordered <- sample_expr[ordered_core_sample, ]

sample_expr_ordered <- t(scale(t(sample_expr_ordered)))

stage_labels <- meta_final[colnames(sample_expr_ordered), "LifeStage"]
col_split    <- factor(stage_labels, levels = c("Stage1", "Stage2", "Stage3"))

cluster_annot_sample  <- final_labels[ordered_core_sample]
traj_shape_annot      <- traj_shape_sample[ordered_core_sample]

shapes_present_sample <- unique(traj_shape_annot)
shape_col_sample      <- all_possible_shapes[names(all_possible_shapes) %in% shapes_present_sample]

row_ha_sample <- rowAnnotation(
  Cluster         = as.factor(cluster_annot_sample),
  TrajectoryShape = traj_shape_annot,
  col             = list(Cluster         = cluster_palette,
                         TrajectoryShape = shape_col_sample),
  show_annotation_name = TRUE,
  annotation_name_gp   = gpar(fontsize = 8)
)

pdf(file.path(output_dir, "07_heatmap_by_stage.pdf"), width = 16, height = 14)
draw(Heatmap(sample_expr_ordered, name = "Expression",
             cluster_rows = FALSE, cluster_columns = FALSE,
             left_annotation = row_ha_sample,
             show_row_names = FALSE, show_column_names = FALSE,
             col = colorRamp2(c(-2, 0, 2), c("#7EC0EE", "white", "#FFD700")),
             row_split = cluster_annot_sample,
             column_split = col_split,
             row_gap = unit(2, "mm"), column_gap = unit(3, "mm"),
             border = TRUE,
             column_title_gp = gpar(fontsize = 12, fontface = "bold"),
             row_title = "Core Genes", row_title_gp = gpar(fontsize = 10)))
dev.off()
cat("Saved 07_heatmap_by_stage.pdf\n")
cat("Additional isoform/OXPHOS validation was removed from the main workflow; generate it separately in 03_generate_figures.R if needed.\n")
message("\n====== Finish ALL ======")
message("FINAL_K (data-driven)    : ", FINAL_K)
message("  Mfuzz Dmin-driven K    : ", FINAL_K, " (max-curvature method)")
message("  Dmin values (K=2..MAX) : ", paste(round(Dmin_df$Dmin, 4), collapse = ", "))
message("STAGE_FDR_CUTOFF         : ", STAGE_FDR_CUTOFF)
message("CORE_GENE_COR_THRESHOLD  : ", CORE_GENE_COR_THRESHOLD)
message("CORE_GENE_MIN_SIZE       : ", CORE_GENE_MIN_SIZE)
message("Pattern assignments      : ")
for (i in seq_len(nrow(centroid_info))) {
  message(sprintf("  Cluster %d: %s (%d genes)", centroid_info$new_cl[i],
                  centroid_info$pattern[i], centroid_info$n_genes[i]))
}
quit(save = "no")
