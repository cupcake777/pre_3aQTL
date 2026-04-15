#!/usr/bin/env Rscript
#
# PAS Motif Alteration Analysis for 3'aQTLs
#
# Goal: Identify 3'aVariants that alter PAS motifs or uridylate-rich elements
#       near annotated poly(A) sites (PolyA_DB), following Li et al. methodology.
#
# Flow:
#   1. Load 3'aQTL results (Stage1-3)
#   2. Load annotated polyA sites from PolyA_DB
#   3. Filter SNPs within 50bp upstream of annotated polyA sites
#   4. Extract genomic context around polyA site, substitute alt allele
#   5. Scan for PAS motif and U-rich element changes (independently)
#   6. Classify: Disruption / Creation / Switch
#   7. Statistics and visualization

suppressPackageStartupMessages({
  library(data.table)
  library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
export_ppt <- "--ppt" %in% args

BSgenome <- BSgenome.Hsapiens.UCSC.hg38

# --- Motif definitions ---
# 15 PAS hexamer variants (Beaudoing et al. 2000 / Tian et al.)
PAS_MOTIFS <- c(
  "AATAAA", "ATTAAA", "AGTAAA", "TATAAA", "GATAAA", "AAGAAA",
  "AATATA", "AATACA", "CATAAA", "ACTAAA", "AATAGA", "AAAAAA",
  "TTTAAA", "AATGAA", "AAACAA"
)

# U-rich / DSE motifs (IUPAC codes: W=A/T, S=C/G, K=G/T, N=any)
URICH_MOTIFS <- c(
  "WTTSTTW", "ATTTA", "WTTTW", "AWTAAA",
  "WWWTTTWWW", "TTTTT", "WWWWTTTWWWW",
  "GTTTG", "WTTKTTW", "WWTTTWW",
  "WTTNTTW", "WWWWWTTTWWWWW"
)

# --- Parameters ---
# SNP must be within 50bp upstream of annotated polyA site
UPSTREAM_WINDOW <- 50
# Sequence context: extract [pA_site - 60, pA_site + 10] to cover
# PAS region (~10-30bp upstream) and downstream U-rich elements
SEQ_UPSTREAM  <- 60
SEQ_DOWNSTREAM <- 10

# ==============================
# 1. Load QTL results
# ==============================
message("[1/7] Loading QTL results...")

stage_files <- list(
  "Stage1" = "../03_QTL/stage1_sig_QTL.txt.gz",
  "Stage2" = "../03_QTL/stage2_sig_QTL.txt.gz",
  "Stage3" = "../03_QTL/stage3_sig_QTL.txt.gz"
)

pvar_files <- list(
  "Stage1" = "../01_pre/genotype.stage1.pvar",
  "Stage2" = "../01_pre/genotype.stage2.pvar",
  "Stage3" = "../01_pre/genotype.stage3.pvar"
)

bed_variant_map_file <- "ALL_sigQTL_snps.bed"

qtl_list <- list()
for (stg in names(stage_files)) {
  if (!file.exists(stage_files[[stg]])) {
    warning(paste("File not found:", stage_files[[stg]]))
    next
  }
  dt <- fread(stage_files[[stg]])
  dt[, Stage := stg]
  qtl_list[[stg]] <- dt
}
all_qtl <- rbindlist(qtl_list, fill = TRUE)

message(sprintf("   Loaded %d QTL records across %d stages.",
                nrow(all_qtl), length(qtl_list)))

# Parse phenotype_id: transcript|gene|chr:start-end|strand
all_qtl[, c("transcript", "gene", "region", "strand") :=
          tstrsplit(phenotype_id, "|", fixed = TRUE)]
all_qtl[, c("chr_pheno", "range") := tstrsplit(region, ":")]
all_qtl[, c("region_start", "region_end") :=
          tstrsplit(range, "-", fixed = TRUE)]
all_qtl[, `:=`(
  region_start = as.numeric(region_start),
  region_end   = as.numeric(region_end)
)]

# Input "#chr" column has bare numbers (1, 2, ...); add "chr" prefix
# Use setnames to handle the "#chr" column name from fread
if ("#chr" %in% names(all_qtl)) {
  setnames(all_qtl, "#chr", "chr_snp")
} else if ("chr" %in% names(all_qtl)) {
  # chr was already overwritten by phenotype parsing — recover from original
  setnames(all_qtl, "chr", "chr_snp")
}
# The SNP chr needs "chr" prefix for BSgenome
all_qtl[, chr := ifelse(grepl("^chr", chr_pheno), chr_pheno,
                        paste0("chr", chr_snp))]

# Deduplicate: unique SNP-phenotype combinations
all_qtl <- unique(all_qtl, by = c("variant_id", "phenotype_id", "Stage"))

message(sprintf("   After dedup: %d unique records.", nrow(all_qtl)))

if (!"pos" %in% names(all_qtl)) {
  message("   Input QTL files do not contain genomic coordinates; recovering from local BED/pvar maps...")

  bed_map <- NULL
  if (file.exists(bed_variant_map_file)) {
    bed_dt <- fread(
      bed_variant_map_file,
      col.names = c("chr_snp", "start0", "end1", "variant_tag", "score", "strand_bed")
    )
    bed_dt[, c("variant_id", "ref", "alt", "gene_bed", "Stage") :=
             tstrsplit(variant_tag, "|", fixed = TRUE)]
    bed_dt[, pos := as.numeric(start0) + 1]
    bed_map <- unique(
      bed_dt[, .(variant_id, gene = gene_bed, Stage, chr_snp, pos, ref, alt)],
      by = c("variant_id", "gene", "Stage")
    )
  }

  if (!is.null(bed_map)) {
    all_qtl <- merge(
      all_qtl,
      bed_map,
      by = c("variant_id", "gene", "Stage"),
      all.x = TRUE
    )

    # data.table::merge may suffix overlapping columns; normalize them here
    if (!"chr_snp" %in% names(all_qtl)) {
      chr_candidates <- intersect(c("chr_snp.x", "chr_snp.y"), names(all_qtl))
      if (length(chr_candidates) > 0) {
        all_qtl[, chr_snp := get(chr_candidates[1])]
        if (length(chr_candidates) == 2) {
          all_qtl[is.na(chr_snp), chr_snp := get(chr_candidates[2])]
        }
      }
    }
    if (!"pos" %in% names(all_qtl)) {
      pos_candidates <- intersect(c("pos.x", "pos.y"), names(all_qtl))
      if (length(pos_candidates) > 0) {
        all_qtl[, pos := get(pos_candidates[1])]
        if (length(pos_candidates) == 2) {
          all_qtl[is.na(pos), pos := get(pos_candidates[2])]
        }
      }
    }
    if (!"ref" %in% names(all_qtl)) {
      ref_candidates <- intersect(c("ref.x", "ref.y"), names(all_qtl))
      if (length(ref_candidates) > 0) {
        all_qtl[, ref := get(ref_candidates[1])]
        if (length(ref_candidates) == 2) {
          all_qtl[is.na(ref), ref := get(ref_candidates[2])]
        }
      }
    }
    if (!"alt" %in% names(all_qtl)) {
      alt_candidates <- intersect(c("alt.x", "alt.y"), names(all_qtl))
      if (length(alt_candidates) > 0) {
        all_qtl[, alt := get(alt_candidates[1])]
        if (length(alt_candidates) == 2) {
          all_qtl[is.na(alt), alt := get(alt_candidates[2])]
        }
      }
    }
    all_qtl[, intersect(c("chr_snp.x", "chr_snp.y", "pos.x", "pos.y",
                          "ref.x", "ref.y", "alt.x", "alt.y"), names(all_qtl)) := NULL]
  }

  if (!"chr_snp" %in% names(all_qtl)) all_qtl[, chr_snp := NA_character_]
  if (!"pos" %in% names(all_qtl)) all_qtl[, pos := NA_real_]
  if (!"ref" %in% names(all_qtl)) all_qtl[, ref := NA_character_]
  if (!"alt" %in% names(all_qtl)) all_qtl[, alt := NA_character_]

  if (!"ref" %in% names(all_qtl) || all_qtl[is.na(pos), .N] > 0) {
    pvar_list <- list()
    for (stg in names(pvar_files)) {
      pvar_path <- pvar_files[[stg]]
      if (!file.exists(pvar_path)) {
        warning(paste("pvar file not found:", pvar_path))
        next
      }

      pvar_dt <- fread(
        cmd = sprintf("grep -v '^##' %s", shQuote(pvar_path))
      )
      setnames(pvar_dt, "#CHROM", "chr_snp")
      setnames(pvar_dt, "POS", "pos")
      setnames(pvar_dt, "ID", "variant_key")
      pvar_dt[, Stage := stg]
      pvar_list[[stg]] <- pvar_dt[, .(variant_key, Stage, chr_snp_pvar = chr_snp, pos_pvar = pos,
                                      ref_pvar = REF, alt_pvar = ALT)]
    }

    if (length(pvar_list) > 0) {
      pvar_map <- rbindlist(pvar_list, fill = TRUE)
      all_qtl[, variant_key := variant_id]
      all_qtl <- merge(all_qtl, pvar_map, by = c("variant_key", "Stage"), all.x = TRUE)
      all_qtl[is.na(chr_snp) & !is.na(chr_snp_pvar), chr_snp := chr_snp_pvar]
      all_qtl[is.na(pos) & !is.na(pos_pvar), pos := pos_pvar]
      all_qtl[is.na(ref) & !is.na(ref_pvar), ref := ref_pvar]
      all_qtl[is.na(alt) & !is.na(alt_pvar), alt := alt_pvar]
      all_qtl[, c("chr_snp_pvar", "pos_pvar", "ref_pvar", "alt_pvar", "variant_key") := NULL]
    }
  }

  missing_pos <- all_qtl[is.na(pos) | is.na(ref) | is.na(alt), .N]
  if (missing_pos > 0) {
    warning(sprintf("   %d QTL rows could not be mapped to chr/pos/ref/alt and will be dropped.", missing_pos))
    all_qtl <- all_qtl[!is.na(pos) & !is.na(ref) & !is.na(alt)]
  }

  all_qtl[, pos := as.numeric(pos)]
  all_qtl[, chr := ifelse(grepl("^chr", chr_snp), chr_snp, paste0("chr", chr_snp))]
}

# ==============================
# 2. Load annotated polyA sites (PolyA_DB)
# ==============================
message("[2/7] Loading PolyA_DB annotations...")

polya_db_file <- "PolyA_DB_hg38.txt.gz"
if (!file.exists(polya_db_file)) stop("PolyA_DB file not found!")

polya_db <- fread(polya_db_file)
# Columns: chr, pA_site, strand, Gene_Symbol, PAS_ID, PAS_ID_hg19, PAS_Signal
message(sprintf("   Loaded %d annotated polyA sites.", nrow(polya_db)))

# ==============================
# 3. Match SNPs to nearby polyA sites (within 50bp upstream)
# ==============================
message("[3/7] Matching SNPs to annotated polyA sites...")

# For + strand gene: SNP upstream of pA means SNP_pos in [pA_site - 50, pA_site]
# For - strand gene: SNP upstream of pA (in transcript) means SNP_pos in [pA_site, pA_site + 50]
#
# Strategy: non-equi join by chr and position range

polya_plus  <- polya_db[strand == "+",
  .(chr, pA_site, pa_strand = strand, Gene_Symbol, PAS_ID, PAS_Signal,
    pa_start = pA_site - UPSTREAM_WINDOW,
    pa_end   = pA_site)]

polya_minus <- polya_db[strand == "-",
  .(chr, pA_site, pa_strand = strand, Gene_Symbol, PAS_ID, PAS_Signal,
    pa_start = pA_site,
    pa_end   = pA_site + UPSTREAM_WINDOW)]

polya_window <- rbind(polya_plus, polya_minus)

# Join: find SNPs that fall within any polyA site's upstream window
setkey(all_qtl, chr, pos)

# Non-equi join: SNP falls in [pa_start, pa_end] on the same chr
matched <- polya_window[all_qtl,
  on = .(chr = chr, pa_start <= pos, pa_end >= pos),
  nomatch = NULL,
  allow.cartesian = TRUE]

# Fix column names after non-equi join
setnames(matched, c("pa_start", "pa_end"), c("snp_pos", "snp_pos2"))
matched[, snp_pos2 := NULL]

# Rename to clarify
setnames(matched, "snp_pos", "pos")

message(sprintf("   Found %d SNP-polyA site pairs within %dbp upstream window.",
                nrow(matched), UPSTREAM_WINDOW))

if (nrow(matched) == 0) {
  stop("No SNPs found near annotated polyA sites. Check chromosome naming or input data.")
}

# Calculate distance: SNP to polyA site (signed, transcript-oriented)
# Negative = upstream of pA site in transcript
matched[, dist_to_pA := ifelse(pa_strand == "+",
                               pos - pA_site,
                               pA_site - pos)]

# ==============================
# 4. Extract context sequences and scan motifs
# ==============================
message("[4/7] Extracting sequences and scanning motifs...")

# Motif detection function (IUPAC-aware via fixed=FALSE)
find_motifs <- function(sequence, motif_vec) {
  seq_dna <- DNAString(sequence)
  hits <- sapply(motif_vec, function(m) {
    length(matchPattern(DNAString(m), seq_dna, fixed = FALSE)) > 0
  })
  motif_vec[hits]
}

# Main scanning function
scan_motif_change <- function(dt) {
  results_list <- vector("list", nrow(dt))
  n_ok <- 0

  for (i in seq_len(nrow(dt))) {
    row <- dt[i]

    tryCatch({
      # Build window around polyA site
      if (row$pa_strand == "+") {
        w_start <- row$pA_site - SEQ_UPSTREAM
        w_end   <- row$pA_site + SEQ_DOWNSTREAM
      } else {
        w_start <- row$pA_site - SEQ_DOWNSTREAM
        w_end   <- row$pA_site + SEQ_UPSTREAM
      }

      # Get reference sequence
      seq_ref <- getSeq(BSgenome, row$chr, start = w_start, end = w_end)

      # Position of SNP within the window
      rel_pos <- row$pos - w_start + 1

      if (rel_pos < 1 || rel_pos > length(seq_ref)) next

      # Create alt sequence
      seq_alt_char <- as.character(seq_ref)
      substr(seq_alt_char, rel_pos, rel_pos) <- row$alt

      # Reverse complement for - strand (so we scan in transcript orientation)
      if (row$pa_strand == "-") {
        final_ref <- as.character(reverseComplement(seq_ref))
        final_alt <- as.character(reverseComplement(DNAString(seq_alt_char)))
      } else {
        final_ref <- as.character(seq_ref)
        final_alt <- seq_alt_char
      }

      # --- PAS motif scanning (independent) ---
      pas_ref  <- find_motifs(final_ref, PAS_MOTIFS)
      pas_alt  <- find_motifs(final_alt, PAS_MOTIFS)

      pas_type <- "None"
      if (length(pas_ref) > 0 && length(pas_alt) == 0) {
        pas_type <- "PAS_disruption"
      } else if (length(pas_ref) == 0 && length(pas_alt) > 0) {
        pas_type <- "PAS_creation"
      } else if (!identical(sort(pas_ref), sort(pas_alt))) {
        pas_type <- "PAS_switch"
      }

      # --- URich motif scanning (independent) ---
      urich_ref <- find_motifs(final_ref, URICH_MOTIFS)
      urich_alt <- find_motifs(final_alt, URICH_MOTIFS)

      urich_type <- "None"
      if (length(urich_ref) > 0 && length(urich_alt) == 0) {
        urich_type <- "URich_disruption"
      } else if (length(urich_ref) == 0 && length(urich_alt) > 0) {
        urich_type <- "URich_creation"
      } else if (!identical(sort(urich_ref), sort(urich_alt))) {
        urich_type <- "URich_switch"
      }

      # Format motif strings
      fmt <- function(m) if (length(m) == 0) "None" else paste(sort(unique(m)), collapse = ",")

      n_ok <- n_ok + 1
      results_list[[n_ok]] <- data.table(
        variant_id   = row$variant_id,
        phenotype_id = row$phenotype_id,
        gene         = row$gene,
        chr          = row$chr,
        pos          = row$pos,
        ref          = row$ref,
        alt          = row$alt,
        strand       = row$pa_strand,
        Stage        = row$Stage,
        pA_site      = row$pA_site,
        PAS_ID       = row$PAS_ID,
        dist_to_pA   = row$dist_to_pA,
        slope        = row$slope,
        pval_nominal = row$pval_nominal,
        PAS_Type     = pas_type,
        PAS_Ref      = fmt(pas_ref),
        PAS_Alt      = fmt(pas_alt),
        URich_Type   = urich_type,
        URich_Ref    = fmt(urich_ref),
        URich_Alt    = fmt(urich_alt)
      )

    }, error = function(e) NULL)

    if (i %% 2000 == 0) message(sprintf("   Processed %d / %d ...", i, nrow(dt)))
  }

  rbindlist(results_list[seq_len(n_ok)])
}

motif_results <- scan_motif_change(matched)

message(sprintf("   Scanned %d SNP-polyA pairs. Results: %d rows.",
                nrow(matched), nrow(motif_results)))

if (nrow(motif_results) == 0) {
  stop("No results produced. Check input data and genome build.")
}

# ==============================
# 5. Classify and summarize
# ==============================
message("[5/7] Classifying results...")

# Any motif alteration (PAS or URich)
motif_results[, is_PAS_alter  := PAS_Type != "None"]
motif_results[, is_URich_alter := URich_Type != "None"]
motif_results[, is_any_alter  := is_PAS_alter | is_URich_alter]

# Combined motif type for summary
motif_results[, Motif_Type := fcase(
  is_PAS_alter & is_URich_alter, paste(PAS_Type, "+", URich_Type),
  is_PAS_alter,                  PAS_Type,
  is_URich_alter,                URich_Type,
  default = "None"
)]

# ==============================
# 6. Save results
# ==============================
message("[6/7] Saving results...")

fwrite(motif_results, "Motif_Mechanism_Table.txt", sep = "\t")

# --- Summary statistics ---

# 6a. Overall counts
alter_summary <- motif_results[, .(
  n_total       = .N,
  n_PAS_alter   = sum(is_PAS_alter),
  n_URich_alter = sum(is_URich_alter),
  n_any_alter   = sum(is_any_alter)
)]
message(sprintf("   Total: %d | PAS-altering: %d | URich-altering: %d | Any: %d",
                alter_summary$n_total, alter_summary$n_PAS_alter,
                alter_summary$n_URich_alter, alter_summary$n_any_alter))

# 6b. Effect size comparison (motif-altering vs non-altering)
if (alter_summary$n_any_alter > 0 && alter_summary$n_any_alter < alter_summary$n_total) {
  effect_test <- wilcox.test(abs(slope) ~ is_any_alter, data = motif_results)
  capture.output(effect_test, file = "EffectSize_Wilcox.txt")
}

# 6c. Stage enrichment
stage_summary <- motif_results[, .(
  n_total       = .N,
  n_PAS_alter   = sum(is_PAS_alter),
  n_URich_alter = sum(is_URich_alter),
  n_any_alter   = sum(is_any_alter)
), by = Stage]
stage_summary[, pct_PAS   := round(n_PAS_alter / n_total, 4)]
stage_summary[, pct_URich := round(n_URich_alter / n_total, 4)]
stage_summary[, pct_any   := round(n_any_alter / n_total, 4)]

fwrite(stage_summary, "Stage_Motif_Enrichment.txt", sep = "\t")

# 6d. PAS motif type breakdown (like Fig. 4b in the paper)
pas_detail <- motif_results[is_PAS_alter == TRUE, .N, by = .(PAS_Type, PAS_Ref, PAS_Alt)]
pas_detail <- pas_detail[order(-N)]
fwrite(pas_detail, "PAS_Detail_Breakdown.txt", sep = "\t")

# 6e. URich motif type breakdown
urich_detail <- motif_results[is_URich_alter == TRUE, .N, by = .(URich_Type, URich_Ref, URich_Alt)]
urich_detail <- urich_detail[order(-N)]
fwrite(urich_detail, "URich_Detail_Breakdown.txt", sep = "\t")

# ==============================
# 7. Visualization
# ==============================
message("[7/7] Generating figures...")

# Fixed stage colors
stage_colors <- c("Stage1" = "#56B4E9", "Stage2" = "#FFB482", "Stage3" = "#8DE5A1")

# Shared publication theme
theme_pub <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      legend.position    = "top",
      legend.title       = element_blank(),
      legend.key.size    = unit(0.4, "cm"),
      strip.background   = element_blank(),
      strip.text         = element_text(face = "bold", size = base_size),
      axis.line          = element_line(colour = "black"),
      panel.border       = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      plot.title         = element_text(face = "bold", size = base_size + 1),
      plot.subtitle      = element_text(size = base_size - 2, colour = "grey40"),
      axis.text          = element_text(size = base_size - 1),
      axis.title         = element_text(size = base_size)
    )
}

# --- Fig 1: PAS & URich alteration counts (faceted by motif class) ---
# Alteration category (disruption/creation/switch) mapped to unified colors
# so the same category reads identically in both PAS and URich panels.
pas_plot_dt   <- motif_results[is_PAS_alter   == TRUE, .N, by = .(Stage, Type = PAS_Type)]
pas_plot_dt[, Class := "PAS"]
urich_plot_dt <- motif_results[is_URich_alter == TRUE, .N, by = .(Stage, Type = URich_Type)]
urich_plot_dt[, Class := "URich"]
count_dt <- rbind(pas_plot_dt, urich_plot_dt)

if (nrow(count_dt) > 0) {
  # Strip the class prefix so fill legend is shared: Disruption / Creation / Switch
  count_dt[, AltType := sub("^(PAS|URich)_", "", Type)]
  count_dt[, AltType := factor(tools::toTitleCase(AltType),
                               levels = c("Disruption", "Creation", "Switch"))]
  count_dt[, Class   := factor(Class, levels = c("PAS", "URich"))]
  count_dt[, Stage   := factor(Stage, levels = c("Stage1", "Stage2", "Stage3"))]

  # Colorblind-friendly: red = disruption, blue = creation, green = switch
  alt_colors <- c("Disruption" = "#D62728", "Creation" = "#1F77B4", "Switch" = "#2CA02C")

  p1 <- ggplot(count_dt, aes(x = Stage, y = N, fill = AltType)) +
    geom_bar(stat = "identity", position = "stack", width = 0.65) +
    geom_text(aes(label = N), position = position_stack(vjust = 0.5),
              color = "white", size = 3, fontface = "bold") +
    scale_fill_manual(values = alt_colors) +
    facet_wrap(~Class, scales = "free_y") +
    theme_pub() +
    labs(title    = "Motif-altering 3'aQTLs per Stage",
         subtitle = "Stacked by alteration type (Disruption / Creation / Switch)",
         y = "Number of SNP-polyA pairs", x = NULL)

  ggsave("Fig1_Motif_Counts.pdf", p1, width = 9, height = 5)
}

# --- Fig 2: Proportion of motif-altering variants (faceted PAS vs URich) ---
if (nrow(stage_summary) > 0) {
  # Merge total N for annotation
  stage_long <- melt(stage_summary,
    id.vars       = c("Stage", "n_total"),
    measure.vars  = c("pct_PAS", "pct_URich"),
    variable.name = "Motif_Class",
    value.name    = "Fraction")
  stage_long[, Motif_Class := factor(gsub("pct_", "", Motif_Class),
                                     levels = c("PAS", "URich"))]
  stage_long[, Pct := Fraction * 100]

  # Build annotation: "Stage1\n(n=819)" on x-axis
  stage_long[, StageN := paste0(Stage, "\n(n=", n_total, ")")]
  stage_long[, StageN := factor(StageN,
    levels = unique(stage_long[order(Stage), StageN]))]

  p2 <- ggplot(stage_long, aes(x = StageN, y = Pct, fill = Stage)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = sprintf("%.1f%%", Pct)),
              vjust = -0.4, size = 3.5, fontface = "bold") +
    scale_fill_manual(values = stage_colors, guide = "none") +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       expand = expansion(mult = c(0, 0.15))) +
    facet_wrap(~Motif_Class) +
    theme_pub() +
    labs(title    = "Proportion of Motif-altering 3'aVariants",
         subtitle = "Among all SNP-polyA pairs within 50 bp upstream window",
         y = "Percentage (%)", x = NULL)

  ggsave("Fig2_Motif_Proportion.pdf", p2, width = 8, height = 5)
}

# --- Fig 3: Motif frequency overview with Stage breakdown (faceted PAS vs URich) ---
# IUPAC codes used in U-rich motifs: W=A/T, S=C/G, K=G/T, N=any
iupac_note <- "IUPAC: W=A/T, K=G/T, S=C/G, N=any"

all_motif_events <- rbind(
  motif_results[is_PAS_alter   == TRUE, .(Motif = PAS_Ref,   Class = "PAS",   Stage)],
  motif_results[is_URich_alter == TRUE, .(Motif = URich_Ref, Class = "URich",  Stage)]
)
# Filter out 'None' (creation events where Ref had no motif)
all_motif_events <- all_motif_events[Motif != "None"]

if (nrow(all_motif_events) > 0) {
  # Expand comma-separated motifs
  all_motif_events <- all_motif_events[,
    .(Motif = unlist(strsplit(Motif, ","))), by = .(Class, Stage)]
  motif_freq <- all_motif_events[, .N, by = .(Class, Stage, Motif)]

  # Order motifs by total frequency within each Class (ascending for coord_flip)
  motif_order <- motif_freq[, .(total = sum(N)), by = .(Class, Motif)]
  motif_order <- motif_order[order(Class, total)]
  motif_freq[, Motif := factor(Motif, levels = unique(motif_order$Motif))]
  motif_freq[, Class := factor(Class, levels = c("PAS", "URich"))]
  motif_freq[, Stage := factor(Stage, levels = c("Stage1", "Stage2", "Stage3"))]

  # Count motifs per panel to set relative heights
  n_pas   <- motif_freq[Class == "PAS",   uniqueN(Motif)]
  n_urich <- motif_freq[Class == "URich",  uniqueN(Motif)]
  row_h   <- 0.45   # cm per motif row
  fig3_h  <- max(6, (n_pas + n_urich) * row_h + 3)

  p3 <- ggplot(motif_freq, aes(x = Motif, y = N, fill = Stage)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    scale_fill_manual(values = stage_colors) +
    coord_flip() +
    facet_wrap(~Class, scales = "free", ncol = 1) +
    theme_pub(base_size = 12) +
    labs(title    = "Frequency of Altered Reference Motifs by Stage",
         subtitle = iupac_note,
         x = NULL, y = "Count")

  ggsave("Fig3_Motif_Frequency.pdf", p3,
         width = 7, height = fig3_h, limitsize = FALSE)
}

message("\nAnalysis completed successfully.")

if (export_ppt) {
  message("[7b/7] Generating PPT-ready figures...")

  library(scales)

  theme_science <- function(base_size = 14) {
    theme_classic(base_size = base_size) +
      theme(
        text               = element_text(family = "sans", color = "black"),
        plot.title         = element_text(face = "bold", size = rel(1.2), hjust = 0),
        plot.subtitle      = element_text(size = rel(0.9), color = "grey30"),
        plot.margin        = margin(15, 15, 15, 15),
        legend.position    = "top",
        legend.direction   = "horizontal",
        legend.background  = element_blank(),
        legend.title       = element_text(face = "bold", size = rel(0.8)),
        legend.text        = element_text(size = rel(0.8)),
        strip.background   = element_rect(fill = "grey95", color = NA),
        strip.text         = element_text(
          face = "bold",
          size = rel(0.9),
          margin = margin(5, 0, 5, 0)
        ),
        axis.line          = element_line(color = "black", linewidth = 0.6),
        axis.text          = element_text(color = "black", size = rel(0.85)),
        axis.title         = element_text(face = "bold", size = rel(0.95)),
        panel.grid.major.y = element_line(color = "grey95", linetype = "solid")
      )
  }

  if (nrow(count_dt) > 0) {
    count_dt_ppt <- copy(count_dt)
    count_dt_ppt[, AltType := factor(
      tools::toTitleCase(sub("^(PAS|URich)_", "", Type)),
      levels = c("Disruption", "Creation", "Switch")
    )]
    count_dt_ppt[, Class := ifelse(Class == "PAS", "PAS motif", "U-rich element")]

    p1_ppt <- ggplot(count_dt_ppt, aes(x = Stage, y = N, fill = AltType)) +
      geom_bar(
        stat = "identity",
        position = position_stack(reverse = TRUE),
        width = 0.7,
        color = "white",
        linewidth = 0.2
      ) +
      geom_text(
        aes(label = ifelse(N > 5, N, "")),
        position = position_stack(vjust = 0.5, reverse = TRUE),
        color = "black",
        size = 3.5,
        fontface = "bold"
      ) +
      scale_fill_manual(values = alt_colors, name = "Alteration Type") +
      facet_wrap(~Class, scales = "free_y") +
      theme_science() +
      labs(
        title = "Number of Motif-Altering 3'aVariants",
        subtitle = "Classified by mechanism: Disruption, Creation, or Motif Switching",
        y = "Counts (SNP-polyA pairs)",
        x = NULL
      )

    ggsave("Fig1_Motif_Counts_PPT.pdf", p1_ppt, width = 10, height = 6)
  }

  if (nrow(stage_summary) > 0) {
    stage_long_ppt <- melt(
      stage_summary,
      id.vars = c("Stage", "n_total"),
      measure.vars = c("pct_PAS", "pct_URich"),
      variable.name = "Class",
      value.name = "Fraction"
    )
    stage_long_ppt[, Class := factor(
      ifelse(Class == "pct_PAS", "PAS motif", "U-rich element")
    )]

    p2_ppt <- ggplot(stage_long_ppt, aes(x = Stage, y = Fraction, fill = Stage)) +
      geom_col(width = 0.65, color = "black", linewidth = 0.1) +
      geom_text(
        aes(label = percent(Fraction, accuracy = 0.1)),
        vjust = -0.5,
        size = 4,
        fontface = "bold",
        color = "black"
      ) +
      scale_fill_manual(values = stage_colors, guide = "none") +
      scale_y_continuous(
        labels = percent_format(),
        expand = expansion(mult = c(0, 0.2))
      ) +
      facet_wrap(~Class) +
      theme_science() +
      labs(
        title = "Prevalence of Motif-Altering Variants",
        subtitle = "Percentage of aVariants within 50bp of a polyA site that alter core motifs",
        y = "Percentage of 3'aVariants (%)",
        x = NULL
      )

    ggsave("Fig2_Motif_Proportion_PPT.pdf", p2_ppt, width = 8, height = 6)
  }

  all_motif_events_ppt <- rbind(
    motif_results[
      is_PAS_alter == TRUE & PAS_Ref != "None",
      .(Motif = PAS_Ref, Class = "PAS motif", Stage)
    ],
    motif_results[
      is_URich_alter == TRUE & URich_Ref != "None",
      .(Motif = URich_Ref, Class = "U-rich element", Stage)
    ]
  )

  if (nrow(all_motif_events_ppt) > 0) {
    motif_freq_ppt <- all_motif_events_ppt[
      ,
      .(Motif = unlist(strsplit(Motif, ","))),
      by = .(Class, Stage)
    ]
    motif_sum_ppt <- motif_freq_ppt[, .N, by = .(Class, Stage, Motif)]
    motif_order_ppt <- motif_sum_ppt[, .(total = sum(N)), by = Motif][order(total)]
    motif_sum_ppt[, Motif := factor(Motif, levels = motif_order_ppt$Motif)]

    p3_ppt <- ggplot(motif_sum_ppt, aes(x = Motif, y = N, fill = Stage)) +
      geom_bar(stat = "identity", position = "stack", width = 0.75, alpha = 0.7) +
      scale_fill_manual(values = stage_colors) +
      coord_flip() +
      facet_wrap(~Class, scales = "free", ncol = 2) +
      theme_science(base_size = 12) +
      theme(panel.grid.major.x = element_line(color = "grey90")) +
      labs(
        title = "Spectrum of Altered Reference Motifs",
        subtitle = "Distribution across developmental stages",
        x = NULL,
        y = "Number of Events"
      )

    ggsave("Fig3_Motif_Frequency_PPT.pdf", p3_ppt, width = 12, height = 8)
  }

  message("[Success] PPT-ready figures saved as PDF.")
}
