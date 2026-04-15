#!/usr/bin/env Rscript
library(data.table)
hg38 <- fread("hg38.pas.bed", col.names = c("chr", "start", "end", "PAS_ID", "score", "strand"))
anno <- fread("~/share_group_folder/ref/human.PAS.txt.gz", select = c("PAS_ID", "Gene Symbol", "PAS Signal", "Mean RPM"))
setnames(anno, "Gene Symbol", "Gene_Symbol")
setnames(anno, "PAS Signal", "PAS_Signal")
setnames(anno, "Mean RPM", "Mean_RPM")
merged <- merge(hg38, anno, by = "PAS_ID")
standard_chrs <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
merged <- merged[chr %in% standard_chrs]
merged[, pA_site := end]
merged[, PAS_ID_hg19 := PAS_ID]
merged[, PAS_ID := paste0(chr, ":", pA_site, ":", strand)]
final_ref <- merged[, .(chr, pA_site, strand, Gene_Symbol, PAS_ID, PAS_ID_hg19, PAS_Signal)]
fwrite(final_ref, "PolyA_DB_hg38.txt.gz", sep = "\t")
