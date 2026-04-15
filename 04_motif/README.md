# 04_motif

This directory contains the PAS motif analysis module for the Three_stages rerun. The module links significant 3'aQTLs to annotated polyadenylation sites, tests whether the associated variants alter canonical PAS or U-rich motifs, summarizes the alteration burden across stages, and generates representative APA mechanism plots for publication-facing review.

## Public scripts

- [`00_build_polya_db_reference.r`](/mnt/share_group_folder/tmp/work/Pre_3aQTL/Three_stages/04_motif/00_build_polya_db_reference.r)
  Prepares PolyA_DB-derived reference resources used by the downstream motif workflow.

- [`01_run_pas_motif_scan.r`](/mnt/share_group_folder/tmp/work/Pre_3aQTL/Three_stages/04_motif/01_run_pas_motif_scan.r)
  Main analysis entrypoint. It loads significant QTL tables, maps variants to annotated polyA sites, scans local sequence context for PAS and U-rich motif changes, writes the main result tables, and exports summary figures.

- [`02_plot_apa_mechanism.r`](/mnt/share_group_folder/tmp/work/Pre_3aQTL/Three_stages/04_motif/02_plot_apa_mechanism.r)
  Unified plotting entrypoint for APA mechanism figures. It supports single-locus plotting and internal batch mode, combining genotype-stratified PDUI distributions with representative coverage tracks and motif annotations.

- [`03_plot_apa_mechanism_batch.r`](/mnt/share_group_folder/tmp/work/Pre_3aQTL/Three_stages/04_motif/03_plot_apa_mechanism_batch.r)
  Compatibility wrapper for batch plotting. It forwards command-line arguments to `02_plot_apa_mechanism.r --batch`.

## Expected outputs

Running `01_run_pas_motif_scan.r` regenerates the core tables and summary figures:

- `Motif_Mechanism_Table.txt`
- `Stage_Motif_Enrichment.txt`
- `PAS_Detail_Breakdown.txt`
- `URich_Detail_Breakdown.txt`
- `EffectSize_Wilcox.txt`
- `Fig1_Motif_Counts.pdf`
- `Fig2_Motif_Proportion.pdf`
- `Fig3_Motif_Frequency.pdf`

Running the mechanism plotting step regenerates per-locus figure folders named as `{gene}_{variant}` and places one publication-style PDF inside each folder.

## Runtime requirements

The module depends on local references and cohort-specific resources that are not bundled into this directory, including the hg38 genome package, the local CHB VCF, and local WIG or bigWig coverage tracks. The scripts are therefore intended to be reproducible inside the project workspace rather than as a standalone portable package.
