#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(stringr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE))
} else {
  normalizePath(getwd(), mustWork = FALSE)
}

resolve_existing_path <- function(label, candidates) {
  candidates <- unique(candidates[nzchar(candidates)])
  for (candidate in candidates) {
    if (file.exists(candidate)) {
      return(normalizePath(candidate, mustWork = FALSE))
    }
  }
  stop(sprintf("Missing %s. Checked:\n%s", label, paste(candidates, collapse = "\n")))
}

resolve_latest_file <- function(dir_path, prefix, fallback_k = 3L) {
  candidates <- Sys.glob(file.path(dir_path, sprintf("%s_K*_all.csv", prefix)))
  if (length(candidates) == 0) {
    return(file.path(dir_path, sprintf("%s_K%d_all.csv", prefix, fallback_k)))
  }
  candidate_k <- suppressWarnings(as.integer(sub(sprintf("^%s_K([0-9]+)_all\\.csv$", prefix), "\\1", basename(candidates))))
  candidates[[order(candidate_k, decreasing = TRUE, na.last = TRUE)[1]]]
}

gsea_input_file <- if (length(args) >= 1) {
  normalizePath(args[1], mustWork = FALSE)
} else {
  resolve_latest_file(script_dir, "05b_GSEA_KEGG")
}

gsea_output_dir <- if (length(args) >= 2) {
  normalizePath(args[2], mustWork = FALSE)
} else {
  file.path(script_dir, "GSEA")
}
if (!dir.exists(gsea_output_dir)) dir.create(gsea_output_dir, recursive = TRUE)

core_file <- resolve_existing_path(
  "core cluster assignments",
  c(file.path(script_dir, "06_final_CORE_cluster_assignments.csv"))
)
full_file <- resolve_existing_path(
  "full cluster assignments",
  c(file.path(script_dir, "06_final_cluster_assignments_with_correlation.csv"))
)

if (!file.exists(gsea_input_file)) {
  stop(sprintf("Missing GSEA input: %s\nRun 01_trajectory_clustering_and_ranked_enrichment.R first.", gsea_input_file))
}

extract_token2 <- function(gene_ids) {
  parts <- str_split_fixed(as.character(gene_ids), fixed("|"), 4)
  if (ncol(parts) >= 2) parts[, 2] else rep("", nrow(parts))
}

build_id_map <- function(token2_vec) {
  token2_vec <- unique(token2_vec[nzchar(token2_vec) & !is.na(token2_vec)])
  symbol_ids <- unique(token2_vec[!grepl("^ENSG", token2_vec)])
  ensembl_ids <- unique(token2_vec[grepl("^ENSG", token2_vec)])

  symbol_map <- data.frame(token2 = character(), SYMBOL = character(), ENTREZID = character(), stringsAsFactors = FALSE)
  if (length(symbol_ids) > 0) {
    symbol_res <- suppressMessages(suppressWarnings(
      bitr(symbol_ids, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    ))
    if (!is.null(symbol_res) && nrow(symbol_res) > 0) {
      symbol_map <- data.frame(
        token2 = symbol_res$SYMBOL,
        SYMBOL = symbol_res$SYMBOL,
        ENTREZID = as.character(symbol_res$ENTREZID),
        stringsAsFactors = FALSE
      )
    }
  }

  ensembl_map <- data.frame(token2 = character(), SYMBOL = character(), ENTREZID = character(), stringsAsFactors = FALSE)
  if (length(ensembl_ids) > 0) {
    chunks <- split(ensembl_ids, ceiling(seq_along(ensembl_ids) / 500))
    ensembl_list <- lapply(chunks, function(chunk_ids) {
      res <- tryCatch(
        suppressMessages(suppressWarnings(
          bitr(chunk_ids, fromType = "ENSEMBL", toType = c("SYMBOL", "ENTREZID"), OrgDb = org.Hs.eg.db)
        )),
        error = function(e) NULL
      )
      if (is.null(res) || nrow(res) == 0) return(NULL)
      data.frame(
        token2 = res$ENSEMBL,
        SYMBOL = res$SYMBOL,
        ENTREZID = as.character(res$ENTREZID),
        stringsAsFactors = FALSE
      )
    })
    ensembl_list <- ensembl_list[!vapply(ensembl_list, is.null, logical(1))]
    if (length(ensembl_list) > 0) {
      ensembl_map <- bind_rows(ensembl_list)
    }
  }

  bind_rows(symbol_map, ensembl_map) %>%
    distinct(token2, .keep_all = TRUE)
}

run_enrich_go <- function(gene_symbols, universe_symbols) {
  if (length(gene_symbols) < 5 || length(universe_symbols) < 10) return(NULL)
  tryCatch(
    enrichGO(
      gene = gene_symbols,
      universe = universe_symbols,
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = FALSE
    ),
    error = function(e) NULL
  )
}

run_enrich_kegg <- function(gene_entrez, universe_entrez) {
  if (length(gene_entrez) < 5 || length(universe_entrez) < 10) return(NULL)
  tryCatch(
    enrichKEGG(
      gene = gene_entrez,
      universe = universe_entrez,
      organism = "hsa",
      keyType = "kegg",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
    ),
    error = function(e) NULL
  )
}

collect_ora <- function(cluster_df, universe_symbols, universe_entrez, type = c("GO", "KEGG")) {
  type <- match.arg(type)
  results <- list()

  for (cl in sort(unique(cluster_df$Cluster))) {
    cl_df <- cluster_df %>% filter(Cluster == cl)
    res <- if (type == "GO") {
      run_enrich_go(unique(na.omit(cl_df$SYMBOL)), universe_symbols)
    } else {
      run_enrich_kegg(unique(na.omit(cl_df$ENTREZID)), universe_entrez)
    }
    if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
      results[[as.character(cl)]] <- as.data.frame(res) %>%
        mutate(Cluster = cl)
    }
  }

  if (length(results) == 0) return(NULL)
  bind_rows(results)
}

make_ora_dotplot <- function(res_df, title) {
  if (is.null(res_df) || nrow(res_df) == 0) return(NULL)

  plot_df <- res_df %>%
    group_by(Cluster) %>%
    slice_min(order_by = p.adjust, n = 5, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      Description = str_wrap(Description, width = 42),
      neglog10 = -log10(p.adjust + 1e-300),
      Cluster = factor(as.character(Cluster), levels = sort(unique(as.character(Cluster))))
    )

  ggplot(plot_df, aes(x = Cluster, y = Description)) +
    geom_point(aes(size = Count, color = neglog10), alpha = 0.9) +
    scale_color_gradient(low = "#9ecae1", high = "#c0392b", name = expression(-log[10](p.adjust))) +
    scale_size_continuous(name = "Gene count") +
    labs(title = title, x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_text(size = 8), plot.title = element_text(face = "bold"))
}

write_result_or_empty <- function(df, file) {
  if (is.null(df)) {
    write.csv(data.frame(), file, row.names = FALSE)
  } else {
    write.csv(df, file, row.names = FALSE)
  }
}

gsea_raw <- read_csv(gsea_input_file, show_col_types = FALSE)
required_cols <- c("Cluster", "Description", "NES", "p.adjust", "pvalue")
missing_cols <- setdiff(required_cols, colnames(gsea_raw))
if (length(missing_cols) > 0) {
  stop(sprintf("GSEA table missing columns: %s", paste(missing_cols, collapse = ", ")))
}

# Annotate all results with significance label; keep full table for archival purposes.
gsea_annotated <- gsea_raw %>%
  mutate(
    Cluster = as.character(Cluster),
    sig_label = case_when(
      !is.na(p.adjust) & p.adjust < 0.05 ~ "FDR<0.05",
      !is.na(pvalue)   & pvalue   < 0.05 ~ "P<0.05",
      TRUE ~ "NS"
    )
  ) %>%
  arrange(Cluster, p.adjust, desc(abs(NES)))

# Full table including nominally significant results (P<0.05, not FDR-corrected) — for reference only.
write_csv(
  filter(gsea_annotated, sig_label != "NS"),
  file.path(gsea_output_dir, "all_clusters_gsea_results_nominal_included.csv")
)

# Primary analysis: restrict to FDR<0.05 for all figures and summary statistics.
sig_gsea <- filter(gsea_annotated, sig_label == "FDR<0.05")

if (nrow(sig_gsea) == 0) {
  stop(sprintf("No FDR<0.05 KEGG GSEA rows found in %s", basename(gsea_input_file)))
}

write_csv(sig_gsea, file.path(gsea_output_dir, "all_clusters_gsea_results.csv"))

summary_df <- sig_gsea %>%
  group_by(Cluster) %>%
  summarise(
    n_pathways = n(),
    n_fdr = sum(sig_label == "FDR<0.05", na.rm = TRUE),
    top_pathway = Description[which.min(p.adjust)],
    top_nes = NES[which.min(p.adjust)],
    .groups = "drop"
  )
write_csv(summary_df, file.path(gsea_output_dir, "cluster_summary_statistics.csv"))

top_terms <- sig_gsea %>%
  group_by(Cluster) %>%
  slice_min(order_by = p.adjust, n = 10, with_ties = FALSE) %>%
  ungroup()

pathway_matrix <- top_terms %>%
  dplyr::select(Cluster, Description, NES) %>%
  distinct() %>%
  pivot_wider(names_from = Cluster, values_from = NES) %>%
  arrange(Description)
write_csv(pathway_matrix, file.path(gsea_output_dir, "pathway_cluster_matrix.csv"))

pathway_sharing <- top_terms %>%
  count(Description, name = "n_clusters") %>%
  arrange(desc(n_clusters), Description)
write_csv(pathway_sharing, file.path(gsea_output_dir, "pathway_sharing_stats.csv"))

matrix_df <- as.data.frame(pathway_matrix)
rownames(matrix_df) <- matrix_df$Description
matrix_df$Description <- NULL
matrix_mat <- as.matrix(matrix_df)
matrix_mat[is.na(matrix_mat)] <- 0

if (nrow(matrix_mat) > 1 && ncol(matrix_mat) > 0) {
  pdf(file.path(gsea_output_dir, "cross_cluster_NES_heatmap.pdf"), width = 8, height = max(5, 0.28 * nrow(matrix_mat) + 2))
  pheatmap(
    matrix_mat,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("#3B82C4", "white", "#C0392B"))(100),
    breaks = seq(-max(abs(matrix_mat), na.rm = TRUE), max(abs(matrix_mat), na.rm = TRUE), length.out = 101),
    fontsize_row = 8,
    fontsize_col = 10,
    main = "Cross-cluster KEGG NES heatmap"
  )
  dev.off()
}

gsea_plot_df <- top_terms %>%
  mutate(
    Description = str_wrap(Description, width = 42),
    neglog10 = -log10(p.adjust),
    Cluster = factor(Cluster, levels = sort(unique(Cluster)))
  )

gsea_plot <- ggplot(gsea_plot_df, aes(x = Cluster, y = Description)) +
  geom_point(aes(size = neglog10, color = NES), alpha = 0.9) +
  scale_color_gradient2(low = "#3B82C4", mid = "white", high = "#C0392B", midpoint = 0, name = "NES") +
  scale_size_continuous(name = expression(-log[10](p.adjust))) +
  labs(
    title = "Cross-cluster ranked enrichment summary",
    subtitle = "KEGG GSEA top pathways from 01_trajectory_clustering_and_ranked_enrichment.R",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_text(size = 8), plot.title = element_text(face = "bold"))

ggsave(
  file.path(gsea_output_dir, "cross_cluster_dotplot.pdf"),
  gsea_plot,
  width = 9,
  height = max(5, 0.3 * length(unique(gsea_plot_df$Description)) + 2)
)

full_assignments <- fread(full_file, data.table = FALSE)
core_assignments <- fread(core_file, data.table = FALSE)
full_assignments$token2 <- extract_token2(full_assignments$Gene)
core_assignments$token2 <- extract_token2(core_assignments$Gene)
id_map <- build_id_map(full_assignments$token2)

full_mapped <- full_assignments %>% left_join(id_map, by = "token2")
core_mapped <- core_assignments %>% left_join(id_map, by = "token2")

universe_symbols <- full_mapped %>%
  filter(!is.na(SYMBOL), nzchar(SYMBOL)) %>%
  pull(SYMBOL) %>%
  unique()
universe_entrez <- full_mapped %>%
  filter(!is.na(ENTREZID), nzchar(ENTREZID)) %>%
  pull(ENTREZID) %>%
  unique()

k_value <- length(unique(core_mapped$Cluster))

mapping_summary <- core_mapped %>%
  group_by(Cluster) %>%
  summarise(
    core_gene_rows = n(),
    unique_gene_ids = n_distinct(Gene),
    mapped_symbols = n_distinct(SYMBOL[!is.na(SYMBOL) & nzchar(SYMBOL)]),
    mapped_entrez = n_distinct(ENTREZID[!is.na(ENTREZID) & nzchar(ENTREZID)]),
    .groups = "drop"
  ) %>%
  mutate(
    universe_symbols = length(universe_symbols),
    universe_entrez = length(universe_entrez)
  )

go_res <- collect_ora(core_mapped, universe_symbols, universe_entrez, type = "GO")
kegg_res <- collect_ora(core_mapped, universe_symbols, universe_entrez, type = "KEGG")

write_result_or_empty(go_res, file.path(script_dir, sprintf("05c_ORA_GO_BP_core_K%d_all.csv", k_value)))
write_result_or_empty(kegg_res, file.path(script_dir, sprintf("05c_ORA_KEGG_core_K%d_all.csv", k_value)))
write.csv(mapping_summary, file.path(script_dir, "05c_ORA_core_mapping_summary.csv"), row.names = FALSE)

go_plot <- make_ora_dotplot(go_res, sprintf("Core-gene ORA GO BP summary (K=%d)", k_value))
if (!is.null(go_plot)) {
  ggsave(
    file.path(script_dir, sprintf("05c_ORA_GO_BP_core_K%d.pdf", k_value)),
    go_plot,
    width = 9,
    height = max(5, 0.3 * length(unique(go_res$Description)) + 2)
  )
}

kegg_plot <- make_ora_dotplot(kegg_res, sprintf("Core-gene ORA KEGG summary (K=%d)", k_value))
if (!is.null(kegg_plot)) {
  ggsave(
    file.path(script_dir, sprintf("05c_ORA_KEGG_core_K%d.pdf", k_value)),
    kegg_plot,
    width = 9,
    height = max(5, 0.3 * length(unique(kegg_res$Description)) + 2)
  )
}

cat(sprintf("Ranked enrichment summary written to %s\n", gsea_output_dir))
cat(sprintf("Core-gene ORA written to %s\n", script_dir))
