#!/usr/bin/env Rscript
# Batch effect recalibration evaluation script

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(gridExtra)
  library(pheatmap)
  library(RColorBrewer)
  library(reshape2)
  library(dplyr)
  library(cowplot)
  library(ggrepel)
  library(patchwork)
  library(cluster)
  library(tibble)
  library(grid)
  library(scales)
})

cat("Loading data...\n")

# Load data with correct dimensions
expr_before <- fread("after.impute.all.raw", data.table = FALSE) %>% 
  column_to_rownames("V1") %>% as.matrix()

expr_after <- fread("after.combat.txt", data.table = FALSE) %>% 
  column_to_rownames("V1") %>% as.matrix()

meta <- fread("../../raw_data/CHB.all.meta", data.table = FALSE) %>% 
  column_to_rownames("sample")

# Batch9 data (transpose so samples are rows)
standar_data <- fread("batch9.after.combat.txt", data.table = FALSE) %>% 
  tibble::column_to_rownames("V1") %>% t() %>% as.matrix()
common_standar <- intersect(rownames(standar_data), rownames(meta))
standar_data <- standar_data[common_standar, , drop = FALSE]

recal_data <- fread("pre.phenotype.combat.int.Stage3.txt", data.table = FALSE) %>% 
  column_to_rownames("V1") %>% as.matrix()

batch_col <- c(batch1="#284F56", batch2="#2494BA", batch3="#28AFB0", batch4="#83BFBF",
  batch5="#F0BA1A", batch6="#BDA055", batch7="#DCBC9B", batch8="#EAD291", batch9="#FFB482")
process_col <- c(Before="#E69F00", After="#1874CD")
LifeStage_col <- c(Stage1="#56B4E9", Stage2="#FFB482", Stage3="#8DE5A1")

# Replicate pairs
rep_samples <- grep("_rep$", colnames(expr_after), value = TRUE)
replicate_pairs <- data.frame(
  Original = sub("_rep$", "", rep_samples),
  Replicate = rep_samples,
  stringsAsFactors = FALSE
) %>% filter(Original %in% colnames(expr_after))

cat("Found", nrow(replicate_pairs), "replicate pairs\n")

# PCA function with error handling
run_pca_df <- function(mat, process_name) {
  mat[is.na(mat)] <- 0
  col_vars <- apply(mat, 2, sd, na.rm = TRUE)
  col_vars[is.na(col_vars)] <- 0
  mat <- mat[, col_vars > 1e-10, drop = FALSE]
  if (ncol(mat) < 2) return(NULL)
  pca <- tryCatch(prcomp(t(mat), center = TRUE, scale. = TRUE), error = function(e) NULL)
  if (is.null(pca)) return(NULL)
  df <- as.data.frame(pca$x[, 1:2]) %>%
    rownames_to_column("Sample") %>%
    mutate(Process = process_name)
  df$Batch <- meta[df$Sample, "batch"]
  df$LifeStage <- meta[df$Sample, "LifeStage"]
  list(df = df, pca = pca)
}

cat("Running PCA...\n")
pca_res_b <- run_pca_df(expr_before, "Before")
pca_res_a <- run_pca_df(expr_after, "After")
pca_res_stand <- run_pca_df(standar_data, "Standard")
pca_res_recal <- run_pca_df(recal_data, "Recalibrated")

get_pca_var <- function(pca_obj) {
  if (is.null(pca_obj)) return(c(0, 0))
  var_expl <- pca_obj$sdev^2 / sum(pca_obj$sdev^2)
  if (length(var_expl) < 2) var_expl <- c(0, 0)
  var_expl[1:2]
}

# Plot helper
plot_pca <- function(res, title) {
  if (is.null(res)) return(NULL)
  ggplot(res$df, aes(PC1, PC2)) +
    geom_point(aes(color = Batch, shape = LifeStage), size = 3, alpha = 0.7) +
    scale_color_manual(values = batch_col) +
    labs(title = title, 
         x = paste0("PC1 (", round(get_pca_var(res$pca)[1]*100, 1), "%)"),
         y = paste0("PC2 (", round(get_pca_var(res$pca)[2]*100, 1), "%)")) +
    theme_bw()
}

cat("Creating plots...\n")

p1 <- plot_pca(pca_res_b, "Before Recalibration")
p2 <- plot_pca(pca_res_a, "After Recalibration")
p3 <- plot_pca(pca_res_stand, "Batch9 (Standard)")
p4 <- plot_pca(pca_res_recal, "Recalibrated")

# Filter NULL plots
plots <- list(p1, p2, p3, p4)
plots <- plots[!sapply(plots, is.null)]

cat("Saving plots...\n")

pdf("evaluation_plots.pdf", width = 14, height = 10)
do.call(grid.arrange, c(plots, ncol = min(2, length(plots))))
dev.off()

cat("Success! Check evaluation_plots.pdf\n")