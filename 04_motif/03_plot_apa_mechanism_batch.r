#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else {
  getwd()
}
script_path <- file.path(script_dir, "02_plot_apa_mechanism.r")

forwarded_args <- c("--batch", args)
cmd <- sprintf("Rscript %s %s",
               shQuote(script_path),
               paste(shQuote(forwarded_args), collapse = " "))

status <- system(cmd)
quit(status = status)
