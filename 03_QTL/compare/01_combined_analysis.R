#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(UpSetR)
  library(qvalue)
})

enable_pi1 <- identical(tolower(Sys.getenv("ENABLE_PI1", "0")), "1")

dir.create("compare", showWarnings = FALSE)
setwd("compare")

stage_files <- c(
  Stage1 = "../downsample/stage1_sig_QTL.txt.gz",
  Stage2 = "../downsample/stage2_sig_QTL.txt.gz",
  Stage3 = "../downsample/stage3_sig_QTL.txt.gz"
)

nominal_files <- c(
  Stage1 = "../downsample/nominal/stage1.cis_qtl_pairs.txt.gz",
  Stage2 = "../downsample/nominal/stage2.cis_qtl_pairs.txt.gz",
  Stage3 = "../downsample/nominal/stage3.cis_qtl_pairs.txt.gz"
)

reference_candidates <- list.files(".", pattern = "\\.3aQTL\\.txt\\.gz$", full.names = FALSE)

make_pheno <- function(phenotype_id) {
  parts <- strsplit(as.character(phenotype_id), "\\|")
  transcript <- sub("\\..*$", "", sapply(parts, `[`, 1))
  gene <- sapply(parts, function(x) if (length(x) >= 2) x[2] else x[1])
  paste(transcript, gene, sep = "@")
}

load_stage <- function(path) {
  dt <- fread(path)
  dt[, pheno := make_pheno(phenotype_id)]
  dt[, key_id := paste(variant_id, pheno, sep = "_")]
  dt
}

load_reference <- function(path) {
  dt <- fread(path)
  if (all(c("rs", "transcript") %in% names(dt))) {
    dt[, key_id := paste(rs, transcript, sep = "_")]
  } else if (all(c("SNP", "transcript") %in% names(dt))) {
    dt[, key_id := paste(SNP, transcript, sep = "_")]
  } else if (all(c("variant_id", "phenotype_id") %in% names(dt))) {
    dt[, key_id := paste(variant_id, make_pheno(phenotype_id), sep = "_")]
  } else {
    stop(sprintf("Unsupported reference format: %s", path))
  }
  dt
}

# NOTE on pi1 estimation: pi1 is computed on the subset of SNP-gene pairs that
# overlap with an external reference QTL set (matched by key_id). If the reference
# set was itself filtered by significance, the resulting p-value distribution may be
# enriched for small p-values, causing pi1 to be overestimated. Results should be
# interpreted as indicative rather than exact. A minimum of n_overlap >= 100 is
# required before estimation (enforced below).
load_nominal_for_pi1 <- function(path) {
  dt <- fread(cmd = sprintf("zcat %s | cut -f1,2,8", shQuote(path)))
  setnames(dt, c("phenotype_id", "variant_id", "pval_nominal"))
  dt[, pheno := make_pheno(phenotype_id)]
  dt[, key_id := paste(variant_id, pheno, sep = "_")]
  dt[, .(key_id, pval_nominal)]
}

stage_data <- lapply(stage_files, load_stage)
stage_sets <- lapply(stage_data, function(dt) unique(dt$key_id))

overlap_summary <- rbindlist(lapply(combn(names(stage_sets), 2, simplify = FALSE), function(pair) {
  s1 <- pair[1]
  s2 <- pair[2]
  shared <- intersect(stage_sets[[s1]], stage_sets[[s2]])
  data.table(
    comparison = paste0(s1, "_vs_", s2),
    stage1 = s1,
    stage2 = s2,
    stage1_n = length(stage_sets[[s1]]),
    stage2_n = length(stage_sets[[s2]]),
    shared_n = length(shared),
    jaccard = length(shared) / length(union(stage_sets[[s1]], stage_sets[[s2]]))
  )
}))
fwrite(overlap_summary, "table_stage_overlap_summary.txt", sep = "\t", quote = FALSE)

pdf("figure_stage_overlap_upset.pdf", width = 10, height = 6)
upset(fromList(stage_sets),
      nsets = length(stage_sets),
      keep.order = TRUE,
      order.by = "freq",
      sets.bar.color = c("#56B4E9", "#FFB482", "#8DE5A1"),
      main.bar.color = "#3C5488",
      matrix.color = "#3C5488")
dev.off()

# Effect-size correlation: lookup approach (avoids winner's curse).
# For each ordered pair (s_anchor, s_lookup):
#   - Take all significant pairs in s_anchor
#   - Look up their nominal slopes in s_lookup (no significance requirement)
# This is done in both directions for each unordered stage pair.
load_nominal_slopes <- function(path) {
  dt <- fread(cmd = sprintf("zcat %s | cut -f2,13", shQuote(path)))
  setnames(dt, c("variant_id", "slope_nominal"))
  # phenotype_id is not available from this cut; use variant_id only for lookup
  dt
}

nominal_slope_data <- lapply(nominal_files, function(path) {
  dt <- fread(cmd = sprintf("zcat %s", shQuote(path)))
  dt[, pheno := make_pheno(phenotype_id)]
  dt[, key_id := paste(variant_id, pheno, sep = "_")]
  dt[, .(key_id, slope_nominal = slope)]
})

cor_rows <- list()
cor_plots <- list()
for (pair in combn(names(stage_data), 2, simplify = FALSE)) {
  s1 <- pair[1]
  s2 <- pair[2]

  for (direction in list(c(s1, s2), c(s2, s1))) {
    s_anchor  <- direction[1]
    s_lookup  <- direction[2]
    dir_label <- paste0(s_anchor, "_sig_lookup_", s_lookup, "_nominal")

    anchor_dt <- stage_data[[s_anchor]][, .(key_id, slope_anchor = slope)]
    lookup_dt <- nominal_slope_data[[s_lookup]]

    merged <- merge(anchor_dt, lookup_dt, by = "key_id")
    if (nrow(merged) < 3) next

    cor_val <- cor(merged$slope_anchor, merged$slope_nominal)
    cor_rows[[dir_label]] <- data.table(
      comparison      = dir_label,
      anchor_stage    = s_anchor,
      lookup_stage    = s_lookup,
      n_anchor_sig    = nrow(anchor_dt),
      n_lookup_found  = nrow(merged),
      pearson_r       = cor_val,
      slope_sign_concordance = mean(sign(merged$slope_anchor) == sign(merged$slope_nominal))
    )

    p <- ggplot(merged, aes(slope_anchor, slope_nominal)) +
      geom_point(alpha = 0.4, size = 0.6, color = "#3C5488") +
      geom_smooth(method = "lm", se = FALSE, color = "#E64B35") +
      labs(
        title    = paste(s_anchor, "(sig) vs", s_lookup, "(nominal lookup)"),
        subtitle = sprintf("r = %.3f, n = %d", cor_val, nrow(merged)),
        x        = paste(s_anchor, "slope (significant pairs)"),
        y        = paste(s_lookup, "slope (nominal lookup)")
      ) +
      theme_bw(base_size = 12)
    cor_plots[[dir_label]] <- p
  }
}

if (length(cor_rows) > 0) {
  fwrite(rbindlist(cor_rows), "table_stage_correlation_summary.txt", sep = "\t", quote = FALSE)
}

if (length(cor_plots) > 0) {
  pdf("figure_stage_effect_size_correlation.pdf", width = 12, height = 4 * length(cor_plots))
  for (nm in names(cor_plots)) print(cor_plots[[nm]])
  dev.off()
}

if (enable_pi1) {
  pi1_rows <- list()
  nominal_pi1_data <- lapply(nominal_files, load_nominal_for_pi1)
  for (ref_name in reference_candidates) {
    ref <- tryCatch(load_reference(ref_name), error = function(e) NULL)
    if (is.null(ref)) next
    ref_keys <- unique(ref$key_id)

    for (stage in names(nominal_files)) {
      nom <- nominal_pi1_data[[stage]]
      sub <- nom[key_id %in% ref_keys]
      if (!"pval_nominal" %in% names(sub) || nrow(sub) < 100) next
      pvals <- sub$pval_nominal
      pvals <- pvals[is.finite(pvals) & !is.na(pvals) & pvals > 0 & pvals <= 1]
      if (length(pvals) < 100) next
      qobj <- qvalue(pvals, pi0.method = "bootstrap")
      pi1_rows[[paste(ref_name, stage, sep = "_")]] <- data.table(
        reference = ref_name,
        stage = stage,
        n_overlap = length(pvals),
        pi1 = 1 - qobj$pi0
      )
    }
  }

  if (length(pi1_rows) > 0) {
    pi1_dt <- rbindlist(pi1_rows)
    fwrite(pi1_dt, "table_reference_pi1_summary.txt", sep = "\t", quote = FALSE)

    p_pi1 <- ggplot(pi1_dt, aes(stage, reference, fill = pi1)) +
      geom_tile(color = "white") +
      geom_text(aes(label = sprintf("%.3f", pi1)), size = 3.5) +
      scale_fill_gradient(low = "#F7FBFF", high = "#08306B") +
      theme_minimal(base_size = 12) +
      labs(title = "pi1 replication against external 3'aQTL references", x = NULL, y = NULL)
    ggsave("figure_reference_pi1_heatmap.pdf", p_pi1, width = 8, height = 4)
  }
}
