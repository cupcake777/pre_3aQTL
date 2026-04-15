# 02_pattern

`02_pattern/` contains the released trajectory-clustering module derived from the current preprocessing workflow.

## Primary scripts

- [`run.sh`](/mnt/share_group_folder/tmp/work/Pre_3aQTL/Three_stages/02_pattern/run.sh)
  Convenience entry point for the module.

- [`01_trajectory_clustering.R`](/mnt/share_group_folder/tmp/work/Pre_3aQTL/Three_stages/02_pattern/01_trajectory_clustering.R)
  Core clustering workflow and model selection.

- [`02_functional_enrichment.R`](/mnt/share_group_folder/tmp/work/Pre_3aQTL/Three_stages/02_pattern/02_functional_enrichment.R)
  Functional interpretation of cluster-level signatures.

- [`03_generate_figures.R`](/mnt/share_group_folder/tmp/work/Pre_3aQTL/Three_stages/02_pattern/03_generate_figures.R)
  Figure generation for summary visualization.

## Key outputs

- `00_stage_sig_fdr0.05.csv`
- `05_core_gene_filter_summary.csv`
- `06_final_CORE_cluster_assignments.csv`
- `06_final_cluster_assignments_with_correlation.csv`
- `06_representative_core_gene_trends.csv`
- `GSEA/all_clusters_gsea_results.csv`
- `01_K_selection_evidence.pdf`
- `06_heatmap_core_genes.pdf`
- `07_heatmap_by_stage.pdf`

## Reproducibility note

The scripts depend on local inputs from `../01_pre/` and upstream metadata under `../../raw_data/`. This release provides a clear workflow contract and summary outputs, but not a standalone environment-independent package.
