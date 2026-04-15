#!/usr/bin/env Rscript
# =============================================================================
# APA Mechanism Validation Figure — Publication Grade
# =============================================================================
# This script builds a single-page composite figure:
# - Left panel: PDUI distribution by genotype with pairwise Wilcoxon tests.
# - Right panel: per-sample coverage tracks around PAS/SNP plus motif annotation.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggpubr)
  library(tibble)
  library(Gviz)
  library(grid)
  library(scales)
  library(GenomicFeatures)
  library(dplyr)
})

# =============================================================================
# 1. CONFIGURATION
# =============================================================================
C <- list(
  snp         = "rs35683888",
  gene        = "RPL22L1",
  id          = "ENST00000295830.13|RPL22L1|chr3:170864875-170866524|-",
  chr         = "chr3",
  snp_pos     = 170866285,
  strand      = "-",

  utr         = c(170864875L, 170866524L),

  pas         = c(0L, 0L),

  pas_ref     = "AAAAAA",
  pas_alt     = "AAGAAA",
  pas_type    = "PAS_switch",   # PAS_creation / PAS_switch / PAS_disruption

  ref_allele  = "A",
  alt_allele  = "G",

  wig_dir     = "/home/lyc/share_group_folder/CHB/wig",
  vcf_file    = "/home/lyc/share_group_folder/raw_data/CHB_all_anno.vcf.gz",
  pheno_file  = "../01_pre/after.combat.txt",
  gtf_file    = "/home/lyc/share_group_folder/ref/gencode.v40.annotation.gtf.gz",
  motif_table = "Motif_Mechanism_Table.txt",

  flank       = 1500L,
  n_per_group = 2L,

  cols        = c("#56B4E9", "#FFB482", "#8DE5A1"),

  out_pdf     = "APA_mechanism_RPL22L1_rs35683888.pdf",
  out_width   = 14,
  out_height  = 5
)

# =============================================================================
# 1b. COMMAND-LINE OVERRIDES
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(args, flag) {
  idx <- which(args == flag)
  if (length(idx) > 0 && idx < length(args)) args[idx + 1] else NULL
}

snp_override      <- parse_arg(args, "--snp")
gene_override     <- parse_arg(args, "--gene")
chr_override      <- parse_arg(args, "--chr")
pos_override      <- parse_arg(args, "--pos")
strand_override   <- parse_arg(args, "--strand")
utr_override      <- parse_arg(args, "--utr")
pas_override      <- parse_arg(args, "--pas")
pas_ref_override  <- parse_arg(args, "--pas-ref")
pas_alt_override  <- parse_arg(args, "--pas-alt")
pas_type_override <- parse_arg(args, "--pas-type")
ref_override      <- parse_arg(args, "--ref")
alt_override      <- parse_arg(args, "--alt")
id_override       <- parse_arg(args, "--id")
out_pdf_override  <- parse_arg(args, "--out-pdf")
batch_mode        <- "--batch" %in% args
table_override    <- parse_arg(args, "--table")
output_base_arg   <- parse_arg(args, "--output-base")
n_override        <- parse_arg(args, "--n")
type_override     <- parse_arg(args, "--type")
pval_override     <- parse_arg(args, "--pval")

if (!is.null(snp_override))      C$snp        <- snp_override
if (!is.null(gene_override))     C$gene       <- gene_override
if (!is.null(chr_override))      C$chr        <- chr_override
if (!is.null(pos_override))      C$snp_pos    <- as.integer(pos_override)
if (!is.null(strand_override))   C$strand     <- strand_override
if (!is.null(ref_override))      C$ref_allele <- ref_override
if (!is.null(alt_override))      C$alt_allele <- alt_override
if (!is.null(pas_ref_override))  C$pas_ref    <- pas_ref_override
if (!is.null(pas_alt_override))  C$pas_alt    <- pas_alt_override
if (!is.null(pas_type_override)) C$pas_type   <- pas_type_override
if (!is.null(id_override))       C$id         <- id_override
if (!is.null(out_pdf_override))  C$out_pdf    <- out_pdf_override

if (!is.null(utr_override)) {
  utr_parts <- as.integer(strsplit(utr_override, ",")[[1]])
  if (length(utr_parts) == 2) C$utr <- utr_parts
}

if (!is.null(pas_override)) {
  pas_parts <- as.integer(strsplit(pas_override, ",")[[1]])
  if (length(pas_parts) >= 1) C$pas[1] <- pas_parts[1]
  if (length(pas_parts) >= 2) C$pas[2] <- pas_parts[2]
}

if (batch_mode) {
  message("=== Batch APA Mechanism Figure Generator ===")

  table_file <- if (!is.null(table_override)) table_override else C$motif_table
  output_base <- if (!is.null(output_base_arg)) output_base_arg else "."
  n_per_type <- if (!is.null(n_override)) as.integer(n_override) else 3L
  pas_type_filter <- if (!is.null(type_override)) type_override else "all"
  pval_thresh <- if (!is.null(pval_override)) as.numeric(pval_override) else 0.05
  pheno_file <- C$pheno_file
  vcf_file <- C$vcf_file

  message(sprintf("  Input table:    %s", table_file))
  message(sprintf("  Output base:    %s", output_base))
  message(sprintf("  Phenotype file: %s", pheno_file))
  message(sprintf("  VCF file:       %s", vcf_file))
  message(sprintf("  Genes per type: %s", n_per_type))
  message(sprintf("  PAS_type filter:%s", pas_type_filter))
  message(sprintf("  Wilcoxon p<=:   %s (all 3 pairs must pass)", pval_thresh))

  if (!file.exists(table_file)) stop("Batch table not found: ", table_file)
  if (!file.exists(pheno_file)) stop("Phenotype file not found: ", pheno_file)

  dt <- fread(table_file)
  setnames(dt, tolower(names(dt)))

  dt <- dt[pas_type %in% c("PAS_creation", "PAS_switch", "PAS_disruption") &
             !is.na(pas_type) & pas_type != "None"]
  if (!is.null(pas_type_filter) && !is.na(pas_type_filter) && pas_type_filter != "all") {
    dt <- dt[pas_type == pas_type_filter]
  }

  setkeyv(dt, c("gene", "pas_type"))
  dt_uniq <- unique(dt, by = c("gene", "pas_type"))
  if ("pval_nominal" %in% names(dt_uniq)) {
    dt_uniq <- dt_uniq[order(pval_nominal)]
  }
  message(sprintf("\n  Unique candidates: %d (after PAS_type filter)", nrow(dt_uniq)))

  message("\n[Pre-filter] Loading phenotype matrix once...")
  pheno_raw <- fread(pheno_file, data.table = FALSE, check.names = FALSE)
  rownames(pheno_raw) <- pheno_raw[[1]]
  pheno_raw <- pheno_raw[, -1, drop = FALSE]
  pheno_samples <- colnames(pheno_raw)
  message(sprintf("  Phenotype: %d features x %d samples", nrow(pheno_raw), ncol(pheno_raw)))

  .fix_samples_batch <- function(pheno_s, geno_df) {
    n0 <- length(intersect(pheno_s, geno_df$sample))
    if (n0 > 0) return(list(pheno_cols = pheno_s, geno = geno_df))

    g1 <- geno_df %>% mutate(sample = sub("^(.+)_\\1$", "\\1", sample))
    n1 <- length(intersect(pheno_s, g1$sample))

    g2 <- geno_df %>% mutate(sample = sub("[_\\.].*$", "", sample))
    n2 <- length(intersect(pheno_s, g2$sample))

    p3 <- sub("^(.+)_\\1$", "\\1", pheno_s)
    n3 <- length(intersect(p3, geno_df$sample))

    best <- which.max(c(n1, n2, n3))
    if (max(n1, n2, n3) == 0) return(NULL)

    if (best == 1) {
      list(pheno_cols = pheno_s, geno = g1)
    } else if (best == 2) {
      list(pheno_cols = pheno_s, geno = g2)
    } else {
      list(pheno_cols = p3, geno = geno_df)
    }
  }

  .detect_chrom_batch <- function(chrom, pos, vcf) {
    bare <- sub("^chr", "", chrom)
    ucsc <- paste0("chr", bare)
    for (ch in c(bare, ucsc)) {
      res <- tryCatch(trimws(system(
        sprintf("bcftools query -r %s:%d-%d -f '%%CHROM\\n' %s 2>/dev/null | head -1",
                ch, pos, pos, vcf), intern = TRUE)),
        error = function(e) ""
      )
      if (length(res) > 0 && nchar(res[1]) > 0) return(ch)
    }
    return(NA_character_)
  }

  .all_pairs_sig <- function(row_dt, pval_thresh) {
    pid <- as.character(row_dt$phenotype_id)
    chr <- as.character(row_dt$chr)
    pos <- as.integer(row_dt$pos)
    ref <- as.character(row_dt$ref)
    alt <- as.character(row_dt$alt)

    if (!pid %in% rownames(pheno_raw)) return(FALSE)

    chr_use <- .detect_chrom_batch(chr, pos, vcf_file)
    if (is.na(chr_use)) return(FALSE)

    geno_cmd <- sprintf(
      "bcftools query -r %s:%d -f '[%%SAMPLE\\t%%GT\\n]' %s | sed 's/0\\/0/0/g; s/0\\/1/1/g; s/1\\/0/1/g; s/1\\/1/2/g; s/0|0/0/g; s/0|1/1/g; s/1|0/1/g; s/1|1/2/g'",
      chr_use, pos, vcf_file
    )
    geno <- tryCatch(
      fread(cmd = geno_cmd, col.names = c("sample", "dosage"),
            colClasses = c("character", "integer"), showProgress = FALSE),
      error = function(e) NULL
    )
    if (is.null(geno) || nrow(geno) == 0) return(FALSE)
    geno <- geno %>% filter(!is.na(dosage), dosage %in% 0:2)

    fix <- .fix_samples_batch(pheno_samples, geno)
    if (is.null(fix)) return(FALSE)
    geno_fixed <- fix$geno
    pheno_cols <- fix$pheno_cols

    geno_labels <- c(
      "0" = paste0(ref, ref),
      "1" = paste0(ref, alt),
      "2" = paste0(alt, alt)
    )

    pdui_vec <- as.numeric(pheno_raw[pid, , drop = TRUE])
    names(pdui_vec) <- pheno_cols

    m_data <- data.frame(
      sample = geno_fixed$sample,
      dosage = geno_fixed$dosage,
      stringsAsFactors = FALSE
    ) %>%
      mutate(PDUI = pdui_vec[sample]) %>%
      filter(!is.na(PDUI), dosage %in% 0:2) %>%
      mutate(Genotype = geno_labels[as.character(dosage)])

    grp_counts <- table(m_data$dosage)
    if (any(c("0", "1", "2") %notin% names(grp_counts))) return(FALSE)
    if (any(grp_counts[c("0", "1", "2")] < 3)) return(FALSE)

    pairs <- list(c(0L, 1L), c(0L, 2L), c(1L, 2L))
    for (pair in pairs) {
      x <- m_data$PDUI[m_data$dosage == pair[1]]
      y <- m_data$PDUI[m_data$dosage == pair[2]]
      p <- tryCatch(wilcox.test(x, y)$p.value, error = function(e) NA_real_)
      if (is.na(p) || p >= pval_thresh) return(FALSE)
    }
    return(TRUE)
  }

  `%notin%` <- Negate(`%in%`)

  message("[Pre-filter] Testing all 3 pairwise Wilcoxon for each candidate...")
  keep_idx <- logical(nrow(dt_uniq))
  for (i in seq_len(nrow(dt_uniq))) {
    row_dt <- dt_uniq[i]
    keep_idx[i] <- .all_pairs_sig(row_dt, pval_thresh)
    if (i %% 10 == 0) {
      message(sprintf("  Tested %d / %d  (kept %d so far)",
                      i, nrow(dt_uniq), sum(keep_idx[1:i])))
    }
  }

  dt_sig <- dt_uniq[keep_idx]
  message(sprintf("\n[Pre-filter] Kept %d / %d candidates (all 3 pairs p < %g)",
                  nrow(dt_sig), nrow(dt_uniq), pval_thresh))
  if (nrow(dt_sig) == 0) quit(status = 0)

  if ("pval_nominal" %in% names(dt_sig)) {
    dt_sig <- dt_sig[order(pval_nominal)]
  }
  selected <- dt_sig[, if (.N > n_per_type) .SD[1:n_per_type] else .SD, by = pas_type]

  .parse_utr <- function(phenotype_id) {
    pid <- as.character(phenotype_id)
    parts <- strsplit(pid, "[|]")[[1]]
    if (length(parts) >= 4) {
      coord_str <- parts[3]
      coords <- strsplit(coord_str, "[-:]")[[1]]
      if (length(coords) >= 3) {
        return(list(
          utr_start = as.integer(coords[2]),
          utr_end = as.integer(coords[3]),
          strand = parts[4]
        ))
      }
    }
    return(NULL)
  }

  script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
  if (length(script_arg) == 0) stop("Cannot resolve current script path for batch mode.")
  script_path <- normalizePath(sub("^--file=", "", script_arg[1]))

  message(sprintf("\n=== Final selection: %d genes to plot ===", nrow(selected)))
  for (i in seq_len(nrow(selected))) {
    row <- selected[i]

    gene <- as.character(row$gene)
    snp <- as.character(row$variant_id)
    chr <- as.character(row$chr)
    pos <- as.integer(row$pos)
    strand <- as.character(row$strand)
    ref <- as.character(row$ref)
    alt <- as.character(row$alt)
    pas_type <- as.character(row$pas_type)
    pas_ref <- as.character(row$pas_ref)
    pas_alt <- as.character(row$pas_alt)
    phenotype_id <- as.character(row$phenotype_id)
    pA_site <- as.integer(row$pa_site)

    utr_info <- .parse_utr(phenotype_id)
    if (is.null(utr_info)) {
      warning(sprintf("  Skip %s: Cannot parse UTR from phenotype_id", gene))
      next
    }

    out_dir <- file.path(output_base, paste0(gene, "_", snp))
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    out_pdf <- file.path(out_dir, sprintf("APA_mechanism_%s_%s.pdf", gene, snp))

    pas_type_label <- switch(
      pas_type,
      "PAS_creation" = "PAS Creation",
      "PAS_switch" = "PAS Switch",
      "PAS_disruption" = "PAS Disruption",
      pas_type
    )

    cmd_args <- sprintf(
      "--snp %s --gene %s --chr %s --pos %d --strand %s --ref %s --alt %s --id '%s' --utr %d,%d --pas %d,%d --pas-ref '%s' --pas-alt '%s' --pas-type '%s' --out-pdf '%s'",
      snp, gene, chr, pos, utr_info$strand, ref, alt,
      phenotype_id,
      utr_info$utr_start, utr_info$utr_end,
      pA_site, utr_info$utr_start,
      pas_ref, pas_alt, pas_type_label, out_pdf
    )

    run_cmd <- sprintf("Rscript %s %s", shQuote(script_path), cmd_args)
    message(sprintf("\n--- Processing: %s | %s (%s) ---", gene, snp, pas_type))
    result <- tryCatch(system(run_cmd), error = function(e) 1)
    if (result == 0) {
      message(sprintf("  Saved: %s", out_pdf))
    } else {
      warning(sprintf("  Failed (exit code: %d)", result))
    }
  }

  message("\n=== Batch processing complete ===")
  quit(status = 0)
}

# Parse motif strings and resolve display/action logic by PAS mechanism type.
.parse_motifs <- function(x) {
  if (is.null(x) || is.na(x) || x == "" || tolower(trimws(x)) == "none") return(character(0))
  ms <- trimws(strsplit(as.character(x), ",")[[1]])
  ms[nchar(ms) == 6 & toupper(ms) != "NONE" & ms != ""]
}

.first_motif <- function(x) {
  ms <- .parse_motifs(x)
  if (length(ms) == 0) NA_character_ else ms[1]
}

.resolve_motif <- function(pas_ref, pas_alt, pas_type) {
  r_all <- .parse_motifs(pas_ref)
  a_all <- .parse_motifs(pas_alt)
  type_low <- tolower(pas_type)

  if (grepl("creation", type_low)) {
    a <- if (length(a_all) > 0) a_all[1] else NA_character_
    r <- if (length(r_all) > 0) r_all[1] else a
    list(show_ref = r,
         show_alt = a,
         label    = "PAS Creation",
         action   = if (!is.na(a)) sprintf("creates %s", a) else "")

  } else if (grepl("disruption", type_low)) {
    r <- if (length(r_all) > 0) r_all[1] else NA_character_
    a <- if (length(a_all) > 0) a_all[1] else r
    list(show_ref = r,
         show_alt = a,
         label    = "PAS Disruption",
         action   = if (!is.na(r)) sprintf("disrupts %s", r) else "")

  } else if (grepl("switch", type_low)) {
    ref_only <- setdiff(r_all, a_all)
    alt_only <- setdiff(a_all, r_all)

    if (length(ref_only) > 0 && length(alt_only) > 0) {
      r <- ref_only[1]; a <- alt_only[1]
    } else if (length(ref_only) > 0) {
      r <- ref_only[1]; a <- if (length(a_all) > 0) a_all[1] else r
    } else if (length(alt_only) > 0) {
      r <- if (length(r_all) > 0) r_all[1] else alt_only[1]; a <- alt_only[1]
    } else {
      r <- if (length(r_all) > 0) r_all[1] else NA_character_
      a <- if (length(a_all) > 0) a_all[1] else NA_character_
    }

    action <- if (!is.na(r) && !is.na(a) && r != a) sprintf("%s -> %s", r, a) else
              if (!is.na(r)) r else ""
    list(show_ref = r,
         show_alt = a,
         label    = "PAS Switch",
         action   = action)

  } else {
    r <- if (length(r_all) > 0) r_all[1] else NA_character_
    a <- if (length(a_all) > 0) a_all[1] else NA_character_
    list(show_ref = r,
         show_alt = a,
         label    = pas_type,
         action   = "")
  }
}

motif_info <- .resolve_motif(C$pas_ref, C$pas_alt, C$pas_type)

# =============================================================================
# 1c. OPTIONAL LOOKUP FROM MOTIF TABLE
# =============================================================================
.lookup_motif_table <- function(snp_id, motif_file) {
  if (!file.exists(motif_file)) {
    message(sprintf("  Motif table not found: %s (using manual config)", motif_file))
    return(NULL)
  }
  
  dt <- tryCatch(
    fread(motif_file, select = c("variant_id", "pA_site", "PAS_Type", "PAS_Ref", "PAS_Alt", 
                                  "URich_Type", "URich_Ref", "URich_Alt", "dist_to_pA"),
          stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  
  # Match by rsID (variant_id column)
  row <- dt[variant_id == snp_id]
  if (nrow(row) == 0) return(NULL)
  
  # Use the first match
  row <- row[1]
  
  message(sprintf("  [Auto] Found SNP %s in motif table:", snp_id))
  message(sprintf("    pA_site=%d, dist_to_pA=%d", row$pA_site, row$dist_to_pA))
  message(sprintf("    PAS: %s (%s -> %s)", row$PAS_Type, row$PAS_Ref, row$PAS_Alt))
  
  list(
    pA_site  = as.integer(row$pA_site),
    PAS_Type = row$PAS_Type,
    PAS_Ref  = row$PAS_Ref,
    PAS_Alt  = row$PAS_Alt,
    dist_to_pA = as.integer(row$dist_to_pA)
  )
}

motif_lookup <- .lookup_motif_table(C$snp, C$motif_table)

if (!is.null(motif_lookup)) {
  C$pas[1] <- motif_lookup$pA_site
  C$pas_type <- motif_lookup$PAS_Type
  C$pas_ref <- motif_lookup$PAS_Ref
  C$pas_alt <- motif_lookup$PAS_Alt
  
  motif_info <- .resolve_motif(C$pas_ref, C$pas_alt, C$pas_type)
  
  message(sprintf("  [Auto] Updated pas[1]=%d from motif table", C$pas[1]))
} else {
  message(sprintf("  [Warn] SNP %s not found in motif table, using config values", C$snp))
  
  id_parts <- strsplit(C$id, "[|]")[[1]]
  if (length(id_parts) >= 4) {
    region <- id_parts[3]  # "chr:start-end"
    coords <- strsplit(region, "[:-]")[[1]]
    if (length(coords) >= 3) {
      utr_start <- as.integer(coords[2])
      utr_end <- as.integer(coords[3])
      if (utr_start > 0 && utr_end > utr_start) {
        C$utr <- c(utr_start, utr_end)
        message(sprintf("  [Auto] Updated UTR from id: %d-%d", utr_start, utr_end))
      }
    }
  }
  if (C$pas[1] == 0) {
    C$pas[1] <- C$snp_pos - 30
    message(sprintf("  [Warn] Using estimated pas[1]=%d (snp_pos - 30)", C$pas[1]))
  }
}

# Derive the plotting window around PAS and ensure SNP/UTR context is included.
geno_labels <- setNames(
  c(paste0(C$ref_allele, C$ref_allele),
    paste0(C$ref_allele, C$alt_allele),
    paste0(C$alt_allele, C$alt_allele)),
  c("0", "1", "2")
)
col_map <- setNames(C$cols, geno_labels)
is_neg  <- (C$strand == "-")


win <- c(
  C$pas[1] - C$flank,
  C$pas[1] + C$flank
)

if (C$snp_pos < win[1]) win[1] <- C$snp_pos - C$flank
if (C$snp_pos > win[2]) win[2] <- C$snp_pos + C$flank

dist_utr_to_pas <- if (C$strand == "-") C$utr[1] - C$pas[1] else C$utr[2] - C$pas[1]
if (dist_utr_to_pas > 0 && dist_utr_to_pas < 20000) {
  if (C$utr[1] < win[1]) win[1] <- C$utr[1] - C$flank
  if (C$utr[2] > win[2]) win[2] <- C$utr[2] + C$flank
  message(sprintf("  [Info] UTR within %.1f kb of PAS - included in window", dist_utr_to_pas/1000))
} else {
  message(sprintf("  [Info] UTR is %.1f kb from PAS - focusing on PAS region only", dist_utr_to_pas/1000))
}

message("=== APA Mechanism Figure ===")
message(sprintf("  Gene: %s  SNP: %s (%s>%s)  Mechanism: %s",
                C$gene, C$snp, C$ref_allele, C$alt_allele, motif_info$label))
message(sprintf("  Window: %s:%d-%d  Strand: %s", C$chr, win[1], win[2], C$strand))
message(sprintf("  SNP_pos=%d  pA_site=%d  window_size=%d bp",
                C$snp_pos, C$pas[1], win[2] - win[1]))

# =============================================================================
# 2. GENOTYPE LOADING
# =============================================================================
message("[1/5] Loading genotypes from VCF...")

.detect_chrom <- function(chrom, pos, vcf) {
  bare <- sub("^chr", "", chrom)
  ucsc <- paste0("chr", bare)
  for (ch in c(bare, ucsc)) {
    res <- tryCatch(trimws(system(
      sprintf("bcftools query -r %s:%d-%d -f '%%CHROM\\n' %s 2>/dev/null | head -1",
              ch, pos, pos, vcf), intern = TRUE)),
      error = function(e) "")
    if (length(res) > 0 && nchar(res[1]) > 0) return(ch)
  }
  stop(sprintf("Cannot find %s:%d in %s", chrom, pos, vcf))
}

chrom_use <- .detect_chrom(C$chr, C$snp_pos, C$vcf_file)

geno_cmd <- sprintf(
  "bcftools query -r %s:%d -f '[%%SAMPLE\\t%%GT\\n]' %s | sed 's/0\\/0/0/g; s/0\\/1/1/g; s/1\\/0/1/g; s/1\\/1/2/g; s/0|0/0/g; s/0|1/1/g; s/1|0/1/g; s/1|1/2/g'",
  chrom_use, C$snp_pos, C$vcf_file
)
geno <- tryCatch(
  fread(cmd = geno_cmd, col.names = c("sample", "dosage"),
        colClasses = c("character", "integer"), showProgress = FALSE),
  error = function(e) stop("bcftools failed: ", e$message)
) %>% filter(!is.na(dosage), dosage %in% 0:2)

if (nrow(geno) == 0) stop("No valid genotypes extracted.")
message(sprintf("  Samples: %s=%d  %s=%d  %s=%d",
                geno_labels["0"], sum(geno$dosage == 0),
                geno_labels["1"], sum(geno$dosage == 1),
                geno_labels["2"], sum(geno$dosage == 2)))

# =============================================================================
# 3. PHENOTYPE LOADING & MERGE
# =============================================================================
message("[2/5] Loading PDUI phenotype...")

pheno_raw <- fread(C$pheno_file, data.table = FALSE, check.names = FALSE)
rownames(pheno_raw) <- pheno_raw[[1]]
pheno_raw <- pheno_raw[, -1, drop = FALSE]

if (!C$id %in% rownames(pheno_raw))
  stop("Phenotype ID not found: ", C$id)

pheno_samples <- colnames(pheno_raw)
geno_samples  <- geno$sample

.fix_samples <- function(pheno_s, geno_df) {
  n0 <- length(intersect(pheno_s, geno_df$sample))
  if (n0 > 0) return(list(pheno_cols = pheno_s, geno = geno_df, n = n0))

  g1 <- geno_df %>% mutate(sample = sub("^(.+)_\\1$", "\\1", sample))
  n1 <- length(intersect(pheno_s, g1$sample))

  g2 <- geno_df %>% mutate(sample = sub("[_\\.].*$", "", sample))
  n2 <- length(intersect(pheno_s, g2$sample))

  p3 <- sub("^(.+)_\\1$", "\\1", pheno_s)
  n3 <- length(intersect(p3, geno_df$sample))

  best <- which.max(c(n1, n2, n3))
  if (max(n1, n2, n3) == 0)
    stop("Cannot match sample names between phenotype and genotype.")

  if (best == 1) {
    message(sprintf("  Auto-fix: de-duplicate VCF IDs (%d matches)", n1))
    list(pheno_cols = pheno_s, geno = g1, n = n1)
  } else if (best == 2) {
    message(sprintf("  Auto-fix: strip VCF ID suffix (%d matches)", n2))
    list(pheno_cols = pheno_s, geno = g2, n = n2)
  } else {
    message(sprintf("  Auto-fix: de-duplicate pheno IDs (%d matches)", n3))
    list(pheno_cols = p3, geno = geno_df, n = n3)
  }
}

fix  <- .fix_samples(pheno_samples, geno)
geno <- fix$geno
colnames(pheno_raw) <- fix$pheno_cols

m_data <- as.data.frame(t(pheno_raw[C$id, , drop = FALSE])) %>%
  rownames_to_column("sample") %>%
  setNames(c("sample", "PDUI")) %>%
  mutate(PDUI = as.numeric(PDUI)) %>%
  inner_join(dplyr::select(geno, sample, dosage), by = "sample") %>%
  filter(!is.na(dosage), !is.na(PDUI)) %>%
  mutate(Genotype = factor(geno_labels[as.character(dosage)],
                           levels = unname(geno_labels)))

message(sprintf("  Merged samples: %d", nrow(m_data)))
if (nrow(m_data) == 0) stop("No samples after merge.")

# =============================================================================
# 4. COVERAGE EXTRACTION
# =============================================================================
message("[3/5] Extracting per-sample coverage tracks...")
options(ucscChromosomeNames = FALSE)

.detect_ext <- function(wig_dir, sample_id) {
  for (ext in c(".bw", ".bigWig", ".bigwig", ".wig", ".bg", ".bedGraph")) {
    f <- file.path(wig_dir, paste0(sample_id, ext))
    if (file.exists(f)) return(ext)
  }
  return(NA_character_)
}

.read_bw <- function(bw_path, chrom, start, end) {
  cmd <- sprintf("bigWigToBedGraph -chrom=%s -start=%d -end=%d %s /dev/stdout 2>/dev/null",
                 chrom, start - 1L, end, bw_path)
  dt <- tryCatch(
    fread(cmd = cmd, col.names = c("chr", "s0", "e1", "value"),
          sep = "\t", showProgress = FALSE),
    error = function(e) NULL
  )
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  data.frame(start = dt$s0 + 1L, end = dt$e1, value = dt$value)
}

.read_wig <- function(wig_path, chrom, start, end) {
  hdr <- tryCatch(readLines(wig_path, n = 5), error = function(e) character(0))
  is_bedgraph  <- any(grepl("^(chr|[0-9])", hdr) & sapply(hdr, function(l) length(strsplit(l, "\t")[[1]]) == 4))
  is_fixedStep <- any(grepl("^fixedStep", hdr))
  is_varStep   <- any(grepl("^variableStep", hdr))

  if (is_bedgraph || (!is_fixedStep && !is_varStep)) {
    awk_script <- "
$1 == c && NF >= 4 {
  st = $2 + 1; en = $3; v = $4 + 0
  if (en < s || st > e) next
  if (st < s) st = s
  if (en > e) en = e
  for (i = st; i <= en; i++) { sum[i] += v; cnt[i]++ }
}
END {
  for (p in sum) printf \"%d\\t%d\\t%.4f\\n\", p, p+1, sum[p]/cnt[p]
}
"
  } else {
    awk_script <- "
BEGIN { chrom_ok=0; step=1; pos=0; past_end=0 }
past_end { next }
/^variableStep/ {
  chrom_ok = (index($0, c) > 0) ? 1 : 0; step=1; pos=0; next
}
/^fixedStep/ {
  chrom_ok = (index($0, c) > 0) ? 1 : 0
  if (chrom_ok) {
    match($0, /start=([0-9]+)/, a); pos = a[1]+0
    match($0, /step=([0-9]+)/,  b); step = b[1]+0
    if (pos > e) past_end = 1
  }
  next
}
/^track|^browser/ { next }
chrom_ok && NF == 2 {
  p = $1+0
  if (p > e) { past_end=1; next }
  if (p >= s) { sum[p] += $2+0; cnt[p]++ }
  next
}
chrom_ok && NF == 1 {
  if (pos >= s && pos <= e) { sum[pos] += $1+0; cnt[pos]++ }
  if (pos > e) past_end = 1
  pos += step; next
}
END {
  for (p in sum) printf \"%d\\t%d\\t%.4f\\n\", p, p+1, sum[p]/cnt[p]
}
"
  }

  awk_f <- tempfile(fileext = ".awk")
  writeLines(awk_script, awk_f)
  on.exit(unlink(awk_f), add = TRUE)

  cmd <- sprintf("awk -v c=%s -v s=%d -v e=%d -f %s %s | sort -k1,1n",
                 chrom, start, end, awk_f, shQuote(wig_path))
  dt <- tryCatch(
    fread(cmd = cmd, col.names = c("start", "end", "value"),
          sep = "\t", showProgress = FALSE),
    error = function(e) { warning("AWK error: ", e$message); NULL }
  )
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  as.data.frame(dt)
}

# .normalize_cov <- function(df) {
#   mx <- max(df$value, na.rm = TRUE)
#   if (is.na(mx) || mx == 0) mx <- 1
#   df$value <- df$value / mx * 100
#   df
# }

.extract_sample_cov <- function(sample_id, wig_dir, chrom, start, end) {
  ext <- .detect_ext(wig_dir, sample_id)
  if (is.na(ext)) {
    warning(sprintf("  No coverage file for sample: %s", sample_id))
    return(NULL)
  }
  fpath <- file.path(wig_dir, paste0(sample_id, ext))
  dt <- if (ext %in% c(".bw", ".bigWig", ".bigwig")) {
    .read_bw(fpath, chrom, start, end)
  } else {
    .read_wig(fpath, chrom, start, end)
  }
  if (is.null(dt) || nrow(dt) == 0) {
    warning(sprintf("  Empty coverage for sample: %s", sample_id))
    return(NULL)
  }
  attr(dt, "max_raw") <- max(dt$value, na.rm = TRUE)
  dt
}

sample_tracks <- list()

for (d in 0:2) {
  label <- geno_labels[as.character(d)]
  color <- col_map[label]

  grp <- m_data %>% filter(dosage == d)
  if (nrow(grp) == 0) next

  med_pdui <- median(grp$PDUI, na.rm = TRUE)
  reps     <- grp %>%
    mutate(.dist = abs(PDUI - med_pdui)) %>%
    arrange(.dist) %>%
    head(C$n_per_group) %>%
    pull(sample)

  message(sprintf("  [%s] Representative samples (PDUI median=%.3f): %s",
                  label, med_pdui, paste(reps, collapse = ", ")))

  for (s_id in reps) {
    cov <- .extract_sample_cov(s_id, C$wig_dir, C$chr, win[1], win[2])
    if (!is.null(cov)) {
      max_raw <- attr(cov, "max_raw")
      if (is.null(max_raw) || is.na(max_raw) || max_raw == 0) max_raw <- 1
      # Store RAW coverage (no normalization)
      sample_tracks[[length(sample_tracks) + 1]] <- list(
        label = label, color = color, sample_id = s_id, cov = cov,
        max_raw = max_raw
      )
      message(sprintf("    %s: %d bins, max depth=%.0f", s_id, nrow(cov), max_raw))
    }
  }
}

if (length(sample_tracks) == 0)
  stop("No coverage data extracted. Check WIG/bigWig paths and format.")

# ---- Calculate Global Y-limit ----
all_max <- 0
for (tk in sample_tracks) {
  if (!is.null(tk$cov)) {
    mx <- max(tk$cov$value, na.rm = TRUE)
    if (!is.na(mx) && mx > all_max) all_max <- mx
  }
}
ylim_global <- c(0, all_max * 1.1)
message(sprintf("  Global Y-axis limit: %.1f", ylim_global[2]))

# =============================================================================
# 5. BUILD FIGURE
# =============================================================================
message("[4/5] Building composite figure...")

# ---- 5a. Left panel: PDUI boxplot ----
present_genos   <- levels(droplevels(m_data$Genotype))
comparisons_all <- list(c(geno_labels["0"], geno_labels["2"]),
                        c(geno_labels["0"], geno_labels["1"]),
                        c(geno_labels["1"], geno_labels["2"]))
comparisons_use <- Filter(function(x) all(x %in% present_genos), comparisons_all)

box_title    <- sprintf("%s  %s", C$gene, C$snp)
box_subtitle <- sprintf("%s  (n=%d)", motif_info$label, nrow(m_data))

p_box <- ggplot(m_data, aes(x = Genotype, y = PDUI, fill = Genotype)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85,
               color = "grey20", linewidth = 0.5) +
  geom_jitter(width = 0.12, size = 1.4, alpha = 0.45, color = "grey30") +
  stat_compare_means(
    comparisons   = comparisons_use,
    method        = "wilcox.test",
    label         = "p.signif",
    tip.length    = 0.012,
    step.increase = 0.12,
    size          = 3.8,
    vjust         = 0.5
  ) +
  scale_fill_manual(values = col_map, guide = "none") +
  scale_x_discrete(labels = present_genos) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.28))) +
  theme_classic(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 11, hjust = 0.5),
    plot.subtitle   = element_text(size = 9, hjust = 0.5, color = "#E63946",
                                   face = "bold"),
    axis.title      = element_text(size = 10),
    axis.text       = element_text(size = 9),
    axis.line       = element_line(color = "grey20", linewidth = 0.4),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
    plot.margin     = margin(8, 6, 8, 6)
  ) +
  labs(
    title    = box_title,
    subtitle = box_subtitle,
    x        = "Genotype",
    y        = "normalized PDUI"
  )

# ---- 5b. Right panel: Gviz tracks ----

axis_tk <- GenomeAxisTrack(
  fontsize = 9, cex = 0.8, fontcolor = "black", col = "black", lwd = 1.0,
  genome = "hg38", labelPos = "below", littleTicks = FALSE,
  ticks = 3,  # only 3 major ticks
  showTitle = FALSE
)

.make_data_track <- function(tk) {
  DataTrack(
    start          = tk$cov$start,
    end            = tk$cov$end,
    data           = matrix(tk$cov$value, nrow = 1),
    chr            = C$chr,
    genome         = "hg38",
    name           = tk$label,
    type           = "histogram",
    fill.histogram = tk$color,
    col.histogram  = NA,
    ylim           = ylim_global,
    background.title  = "white",
    fontcolor.title   = tk$color,
    fontface.title    = 2,
    fontsize.title    = 7,
    col.axis          = "grey40",
    cex.axis          = 0.45,
    background.panel  = "white"
  )
}

data_tks <- lapply(sample_tracks, .make_data_track)

# Gene model
gene_tk <- tryCatch({
  message("  Building GeneRegionTrack from GTF...")
  txdb <- makeTxDbFromGFF(
    C$gtf_file,
    format    = "gtf",
    organism  = "Homo sapiens",
    circ_seqs = character(0)
  )
  GeneRegionTrack(
    txdb,
    chromosome       = C$chr,
    start            = win[1],
    end              = win[2],
    genome           = "hg38",
    strand           = C$strand,
    name             = C$gene,
    showId           = TRUE,
    geneSymbol       = TRUE,
    fill             = "black",
    col              = "black",
    col.line         = "black",
    utr3             = "grey60",
    utr5             = "grey60",
    protein_coding   = "black",
    collapseTranscripts = "longest",
    fontsize.group   = 6,
    fontcolor.group  = "black",
    background.title = "white",
    fontcolor.title  = "grey20",
    fontface.title   = 2,
    fontsize.title   = 7
  )
}, error = function(e) {
  warning("GeneRegionTrack failed; using fallback.")
  AnnotationTrack(
    start = min(C$utr[1], C$pas[1]), end = max(C$utr[2], C$pas[1]),
    chr = C$chr, genome = "hg38", strand = C$strand, id = C$gene,
    name = C$gene, shape = "box",
    fill = "grey55", col = "grey35",
    fontcolor.item = "white", fontface.item = 2, fontsize.item = 7,
    background.title = "white", fontcolor.title = "grey20",
    fontface.title = 2, fontsize.title = 7,
    showFeatureId = TRUE
  )
})

n_data  <- length(data_tks)
all_tks <- c(list(axis_tk), data_tks, list(gene_tk))
tk_sizes <- c(
  0.50,           # axis
  rep(0.80, n_data),  # coverage tracks
  0.60            # gene model
)
tk_sizes_norm <- tk_sizes / sum(tk_sizes) * 10

# =============================================================================
# 6. MOTIF ANNOTATION HELPER
# =============================================================================
.draw_motif <- function(show_ref, show_alt, mechanism_label,
                        snp_id, ref_allele, alt_allele, action,
                        x_center = 0.5, y_pos = 0.06, snp_x_npc = NULL) {

  has_ref <- !is.null(show_ref) && !is.na(show_ref) && nchar(show_ref) == 6
  has_alt <- !is.null(show_alt) && !is.na(show_alt) && nchar(show_alt) == 6

  if (!has_ref && !has_alt) return(invisible(NULL))

  type_low <- tolower(mechanism_label)

  if (grepl("creation", type_low)) {
    seq_show <- if (has_alt) show_alt else show_ref
    r_vec <- if (has_ref) strsplit(show_ref, "")[[1]] else rep("?", 6)
    a_vec <- strsplit(seq_show, "")[[1]]
    snp_idx <- integer(0)
    label_prefix <- "creates: "
  } else if (grepl("disruption", type_low)) {
    seq_show <- if (has_ref) show_ref else show_alt
    r_vec <- strsplit(seq_show, "")[[1]]
    a_vec <- if (has_alt) strsplit(show_alt, "")[[1]] else r_vec
    snp_idx <- if (length(r_vec) == length(a_vec)) which(r_vec != a_vec) else integer(0)
    # Fallback: if pas_alt=None, ref_allele position is inferred by matching ref_allele
    # in the motif sequence (pick the central-most match)
    if (length(snp_idx) == 0 && !is.null(ref_allele) && !is.na(ref_allele) && nchar(ref_allele) == 1) {
      matches <- which(toupper(r_vec) == toupper(ref_allele))
      if (length(matches) > 0) {
        # pick the one closest to the center of the 6bp motif
        center_pos <- 3.5
        snp_idx <- matches[which.min(abs(matches - center_pos))]
      }
    }
    label_prefix <- "disrupts: "
  } else {
    seq_show <- if (has_ref) show_ref else show_alt
    r_vec <- if (has_ref) strsplit(show_ref, "")[[1]] else rep("?", 6)
    a_vec <- if (has_alt) strsplit(show_alt, "")[[1]] else r_vec
    snp_idx <- if (length(r_vec) == length(a_vec)) which(r_vec != a_vec) else integer(0)
    label_prefix <- "ref motif: "
  }

  n       <- length(r_vec)
  unit_w  <- 0.018
  start_x <- x_center - (n * unit_w) / 2

  grid.text(
    label_prefix,
    x    = unit(start_x - 0.004, "npc"),
    y    = unit(y_pos, "npc"),
    just = "right",
    gp   = gpar(fontsize = 8, fontface = "italic", col = "grey40")
  )

  for (i in seq_along(r_vec)) {
    is_mut <- i %in% snp_idx
    grid.text(
      r_vec[i],
      x  = unit(start_x + (i - 0.5) * unit_w, "npc"),
      y  = unit(y_pos, "npc"),
      gp = gpar(
        col        = if (is_mut) "#E63946" else "black",
        fontface   = if (is_mut) 2L else 1L,
        fontsize   = 14,
        fontfamily = "mono"
      )
    )
  }

  # alt_row_offset = 0.07 npc (>font size 14 ~0.019 npc) to avoid overlap with ref row
  alt_row_offset <- 0.07
  if (grepl("switch", type_low) && has_alt) {
    a_vec2 <- strsplit(show_alt, "")[[1]]
    grid.text(
      "alt motif: ",
      x    = unit(start_x - 0.004, "npc"),
      y    = unit(y_pos - alt_row_offset, "npc"),
      just = "right",
      gp   = gpar(fontsize = 8, fontface = "italic", col = "grey40")
    )
    for (i in seq_along(a_vec2)) {
      is_mut <- i %in% snp_idx
      grid.text(
        a_vec2[i],
        x  = unit(start_x + (i - 0.5) * unit_w, "npc"),
        y  = unit(y_pos - alt_row_offset, "npc"),
        gp = gpar(
          col        = if (is_mut) "#E63946" else "black",
          fontface   = if (is_mut) 2L else 1L,
          fontsize   = 14,
          fontfamily = "mono"
        )
      )
    }
  }

  sx <- if (!is.null(snp_x_npc) && !is.na(snp_x_npc)) snp_x_npc else x_center

  # For switch mode: labels go below the alt motif row; otherwise below the ref row
  label_y_base <- if (grepl("switch", type_low) && has_alt) y_pos - alt_row_offset else y_pos

  grid.segments(
    x0 = unit(sx, "npc"), y0 = unit(label_y_base - 0.010, "npc"),
    x1 = unit(sx, "npc"), y1 = unit(y_pos + 0.028, "npc"),
    gp = gpar(col = "#E63946", lty = 2, lwd = 1.2)
  )

  # rsid + allele on one compact line below the connector
  snp_label    <- if (!is.null(snp_id) && !is.na(snp_id)) snp_id else ""
  allele_label <- sprintf("(%s>%s)", ref_allele, alt_allele)

  grid.text(
    snp_label,
    x = unit(sx, "npc"), y = unit(label_y_base - 0.022, "npc"),
    just = "centre",
    gp = gpar(fontsize = 7.5, fontface = "bold", col = "#E63946")
  )
  grid.text(
    allele_label,
    x = unit(sx, "npc"), y = unit(label_y_base - 0.034, "npc"),
    just = "centre",
    gp = gpar(fontsize = 7, col = "#E63946")
  )

}

# =============================================================================
# 7. RENDER & EXPORT
# =============================================================================
message(sprintf("[5/5] Rendering to %s ...", C$out_pdf))

# Pre-compute SNP x position in npc within right viewport
{
  plot_xmin_npc <- 0.12
  plot_xmax_npc <- 1.00
  win_bp        <- win[2] - win[1]
  snp_frac      <- (C$snp_pos - win[1]) / win_bp
  if (is_neg) snp_frac <- 1 - snp_frac
  snp_x_npc_right <- plot_xmin_npc + snp_frac * (plot_xmax_npc - plot_xmin_npc)
}

left_frac  <- 2.2 / 10   # narrower left panel (boxplot)
right_frac <- 7.8 / 10

pdf(C$out_pdf, width = C$out_width, height = C$out_height, useDingbats = FALSE)
grid.newpage()

# ── Left panel ──
vp_left <- viewport(
  x      = unit(left_frac / 2, "npc"),
  y      = unit(0.5, "npc"),
  width  = unit(left_frac, "npc"),
  height = unit(1, "npc"),
  just   = c("centre", "centre")
)
print(p_box, vp = vp_left)

# ── Right panel ──
vp_right <- viewport(
  x      = unit(left_frac + right_frac / 2, "npc"),
  y      = unit(0.48, "npc"),  # slightly lower to give top margin
  width  = unit(right_frac, "npc"),
  height = unit(0.88, "npc"),  # reduced height for axis labels
  just   = c("centre", "centre")
)
pushViewport(vp_right)

# Pre-compute track layout fractions (used for SNP line and motif position)
tk_total      <- 0.50 + n_data * 0.80 + 0.60   # axis + cov*n + gene
gene_frac_npc <- 0.60 / tk_total        # npc fraction for gene track height
motif_frac_npc <- (0.60 + 0.10) / tk_total  # above gene track to avoid overlap

# Highlight ranges - SNP (lightblue) + Proximal PAS (mistyrose)
hl_ranges <- tryCatch(
  GenomicRanges::GRanges(
    seqnames = C$chr,
    ranges   = IRanges::IRanges(
      start = c(C$snp_pos - 3L, C$pas[1] - 5L),
      end   = c(C$snp_pos + 3L, C$pas[1] + 5L)
    ),
    col      = c("lightblue", "mistyrose")
  ),
  error = function(e) NULL
)

plotTracks(
  trackList        = all_tks,
  from             = win[1],
  to               = win[2],
  chromosome       = C$chr,
  sizes            = tk_sizes_norm,
  reverseStrand    = is_neg,
  highlight        = hl_ranges,
  innerMargin      = 8,       # increased for axis labels
  background.panel = "white",
  background.title = "white",
  col.border.title = "white",
  showTitle        = TRUE,
  add              = TRUE
)

# ── Red dashed SNP line ──
snp_line_y0 <- motif_frac_npc   # just above gene track
snp_line_y1 <- 0.94
grid.segments(
  x0 = unit(snp_x_npc_right, "npc"), y0 = unit(snp_line_y0, "npc"),
  x1 = unit(snp_x_npc_right, "npc"), y1 = unit(snp_line_y1, "npc"),
  gp = gpar(col = "#E63946", lty = 2, lwd = 1.5)
)

motif_track_y <- motif_frac_npc  # just above gene track

motif_seq <- if (!is.na(motif_info$show_ref)) motif_info$show_ref else 
             if (!is.na(motif_info$show_alt)) motif_info$show_alt else "N/A"
motif_vec <- strsplit(motif_seq, "")[[1]]

snp_idx <- integer(0)
if (grepl("disruption", tolower(motif_info$label))) {
  if (length(motif_vec) == 6 && nchar(C$ref_allele) == 1) {
    matches <- which(toupper(motif_vec) == toupper(C$ref_allele))
    if (length(matches) > 0) {
      center_pos <- 3.5
      snp_idx <- matches[which.min(abs(matches - center_pos))]
    }
  }
} else if (grepl("switch", tolower(motif_info$label))) {
  ref_seq <- motif_info$show_ref
  alt_seq <- motif_info$show_alt
  if (!is.na(ref_seq) && !is.na(alt_seq) && nchar(ref_seq) == 6 && nchar(alt_seq) == 6) {
    r_vec <- strsplit(ref_seq, "")[[1]]
    a_vec <- strsplit(alt_seq, "")[[1]]
    snp_idx <- which(r_vec != a_vec)
  }
} else if (grepl("creation", tolower(motif_info$label))) {
  ref_seq <- motif_info$show_ref
  alt_seq <- motif_info$show_alt
  if (!is.na(ref_seq) && !is.na(alt_seq) && nchar(ref_seq) == 6 && nchar(alt_seq) == 6) {
    r_vec <- strsplit(ref_seq, "")[[1]]
    a_vec <- strsplit(alt_seq, "")[[1]]
    snp_idx <- which(r_vec != a_vec)
  } else {
    snp_idx <- integer(0)
  }
}

is_switch_mode <- grepl("switch", tolower(motif_info$label))
unit_w <- 0.012

action_prefix <- switch(
  motif_info$label,
  "PAS Disruption" = "disrupts:",
  "PAS Creation"   = "creates:",
  "PAS Switch"      = "ref:",
  tolower(motif_info$label)
)

if (is_switch_mode) {
  motif_center_x <- snp_x_npc_right+0.01
  motif_start_x <- motif_center_x - (6 * unit_w) / 2
  ref_y <- motif_track_y + 0.025  # higher position
  alt_y <- motif_track_y - 0.010  # lower position
  
  grid.text(
    "ref:",
    x = unit(motif_start_x - 0.025, "npc"),
    y = unit(ref_y, "npc"),
    just = "right",
    gp = gpar(fontsize = 8, fontface = "italic", col = "grey40")
  )
  for (i in seq_along(motif_vec)) {
    is_mut <- i %in% snp_idx
    grid.text(
      motif_vec[i],
      x = unit(motif_start_x + (i - 1) * unit_w, "npc"),
      y = unit(ref_y, "npc"),
      gp = gpar(
        col = if (is_mut) "#E63946" else "black",
        fontface = if (is_mut) 2L else 1L,
        fontsize = 11,
        fontfamily = "mono"
      )
    )
  }
  
  # alt motif (lower row)
  alt_seq <- motif_info$show_alt
  if (!is.na(alt_seq) && nchar(alt_seq) == 6) {
    alt_vec <- strsplit(alt_seq, "")[[1]]
    grid.text(
      "alt:",
      x = unit(motif_start_x - 0.025, "npc"),
      y = unit(alt_y, "npc"),
      just = "right",
      gp = gpar(fontsize = 8, fontface = "italic", col = "grey40")
    )
    for (i in seq_along(alt_vec)) {
      is_mut <- i %in% snp_idx
      grid.text(
        alt_vec[i],
        x = unit(motif_start_x + (i - 1) * unit_w, "npc"),
        y = unit(alt_y, "npc"),
        gp = gpar(
          col = if (is_mut) "#E63946" else "black",
          fontface = if (is_mut) 2L else 1L,
          fontsize = 11,
          fontfamily = "mono"
        )
      )
    }
  }
  
  # rsid below alt row - position higher to avoid gene track overlap
  rsid_y <- alt_y - 0.020
} else {
  motif_center_x <- snp_x_npc_right+0.01
  motif_start_x <- motif_center_x - (6 * unit_w) / 2
  grid.text(
    action_prefix,
    x = unit(motif_start_x - 0.025, "npc"),
    y = unit(motif_track_y, "npc"),
    just = "right",
    gp = gpar(fontsize = 9, fontface = "italic", col = "grey40")
  )
  
  for (i in seq_along(motif_vec)) {
    is_mut <- i %in% snp_idx
    grid.text(
      motif_vec[i],
      x = unit(motif_start_x + (i - 1) * unit_w, "npc"),
      y = unit(motif_track_y, "npc"),
      gp = gpar(
        col = if (is_mut) "#E63946" else "black",
        fontface = if (is_mut) 2L else 1L,
        fontsize = 12,
        fontfamily = "mono"
      )
    )
  }
  rsid_y <- motif_track_y
}

rsid_x <- snp_x_npc_right + 0.04
grid.text(
  sprintf("%s (%s>%s)", C$snp, C$ref_allele, C$alt_allele),
  x = unit(rsid_x, "npc"),
  y = unit(rsid_y, "npc"),
  just = "left",
  gp = gpar(fontsize = 9, fontface = "bold", col = "#E63946")
)

# ── 1 kb scale bar (top-right) ──
win_total_bp   <- win[2] - win[1]
scale_bp       <- min(500L, as.integer(win_total_bp / 2))  # 500 bp or half window
scale_frac     <- scale_bp / win_total_bp
plot_width_npc <- plot_xmax_npc - plot_xmin_npc
scale_npc      <- scale_frac * plot_width_npc
bar_xr  <- plot_xmax_npc - 0.01
bar_xl  <- bar_xr - scale_npc
bar_y   <- 0.96

grid.segments(
  x0 = unit(bar_xl, "npc"), y0 = unit(bar_y, "npc"),
  x1 = unit(bar_xr, "npc"), y1 = unit(bar_y, "npc"),
  gp = gpar(col = "black", lwd = 1.5)
)
for (bx in c(bar_xl, bar_xr)) {
  grid.segments(
    x0 = unit(bx, "npc"), y0 = unit(bar_y - 0.006, "npc"),
    x1 = unit(bx, "npc"), y1 = unit(bar_y + 0.006, "npc"),
    gp = gpar(col = "black", lwd = 1.5)
  )
}
grid.text(
  sprintf("%d bp", scale_bp),
  x    = unit((bar_xl + bar_xr) / 2, "npc"),
  y    = unit(bar_y + 0.016, "npc"),
  just = "centre",
  gp   = gpar(fontsize = 7, fontface = "bold", col = "black")
)

# ── Motif annotation ──
snp_line_y0 <- motif_frac_npc - 0.02   # just below motif track
snp_line_y1 <- 0.94

popViewport()

dev.off()
message(sprintf("Done. Figure saved: %s", C$out_pdf))

# =============================================================================
# 8. SUMMARY STATISTICS
# =============================================================================
cat("\n=== PDUI Summary by Genotype ===\n")
m_data %>%
  group_by(Genotype) %>%
  summarise(
    n      = n(),
    median = round(median(PDUI, na.rm = TRUE), 3),
    mean   = round(mean(PDUI,   na.rm = TRUE), 3),
    sd     = round(sd(PDUI,     na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  as.data.frame() %>%
  print()

cat("\n=== Coverage tracks rendered ===\n")
for (tk in sample_tracks) {
  cat(sprintf("  [%s] %s — %d bins\n", tk$label, tk$sample_id, nrow(tk$cov)))
}
