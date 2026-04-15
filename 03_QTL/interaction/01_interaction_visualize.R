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

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)])
script_dir <- dirname(normalizePath(script_path))

sig_path <- file.path(script_dir, "interaction_sig_QTL.txt.gz")
if (!file.exists(sig_path)) {
  stop("Missing interaction_sig_QTL.txt.gz. Run 04_interaction.sh first.")
}

dt <- read_table_auto(sig_path)
if (!"slope" %in% names(dt)) {
  stop("Expected a slope column in interaction_sig_QTL.txt.gz.")
}

p <- ggplot(dt, aes(x = slope)) +
  geom_histogram(bins = 50, fill = "#B2182B", color = "white") +
  theme_bw(base_size = 12) +
  labs(
    title = "Interaction QTL effect-size distribution",
    x = "Interaction slope",
    y = "Count"
  )

ggsave(file.path(script_dir, "figure_interaction_effect_size_histogram.pdf"), p, width = 7, height = 5)
