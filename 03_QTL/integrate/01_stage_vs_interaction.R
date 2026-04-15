#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
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

  fread(cmd = sprintf("gzip -dc %s", shQuote(path)))
}

make_pheno <- function(phenotype_id) {
  parts <- strsplit(as.character(phenotype_id), "\\|")
  transcript <- sub("\\..*$", "", sapply(parts, `[`, 1))
  gene <- sapply(parts, function(x) if (length(x) >= 2) x[2] else x[1])
  paste(transcript, gene, sep = "@")
}

interaction_path <- "../interaction/interaction_sig_QTL.txt.gz"
if (!file.exists(interaction_path)) {
  stop("Missing ../interaction/interaction_sig_QTL.txt.gz")
}

interaction_dt <- read_table_auto(interaction_path)
interaction_dt[, pheno := make_pheno(phenotype_id)]
interaction_dt[, key_id := paste(variant_id, pheno, sep = "_")]

stage_files <- c(
  Stage1 = "../stage1_sig_QTL.txt.gz",
  Stage2 = "../stage2_sig_QTL.txt.gz",
  Stage3 = "../stage3_sig_QTL.txt.gz"
)

summary_rows <- list()
plot_rows <- list()
for (stage in names(stage_files)) {
  dt <- read_table_auto(stage_files[[stage]])
  dt[, pheno := make_pheno(phenotype_id)]
  dt[, key_id := paste(variant_id, pheno, sep = "_")]
  merged <- merge(
    dt[, .(key_id, slope_stage = slope)],
    interaction_dt[, .(key_id, slope_interaction = slope)],
    by = "key_id"
  )
  summary_rows[[stage]] <- data.table(
    stage = stage,
    stage_n = nrow(dt),
    interaction_n = nrow(interaction_dt),
    shared_n = nrow(merged),
    concordance = if (nrow(merged) == 0) NA_real_ else mean(sign(merged$slope_stage) == sign(merged$slope_interaction))
  )
  if (nrow(merged) > 0) {
    merged[, stage := stage]
    plot_rows[[stage]] <- merged
  }
}

summary_dt <- rbindlist(summary_rows)
fwrite(summary_dt, "table_stage_vs_interaction_summary.txt", sep = "\t", quote = FALSE)

if (length(plot_rows) > 0) {
  plot_dt <- rbindlist(plot_rows)
  p <- ggplot(plot_dt, aes(slope_stage, slope_interaction)) +
    geom_point(alpha = 0.5, size = 0.8, color = "#3C5488") +
    geom_smooth(method = "lm", se = FALSE, color = "#E64B35") +
    facet_wrap(~stage, scales = "free") +
    theme_bw(base_size = 12) +
    labs(
      title = "Stage-specific vs interaction effect-size comparison",
      x = "Stage-specific slope",
      y = "Interaction slope"
    )
  ggsave("figure_stage_vs_interaction_scatter.pdf", p, width = 11, height = 4)
}
