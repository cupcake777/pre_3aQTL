# 06_smr

`06_smr/` contains the SMR follow-up layer used after the stage-QTL and
colocalization preparation steps. The module consumes stage-level BESD eQTL
summaries together with GWAS `.ma` files and produces compact per-trait SMR
outputs plus a final integrative summary table.

## Primary scripts

- [`scripts/05_run_from_05_coloc.sh`](/mnt/share_group_folder/tmp/work/Pre_3aQTL/Three_stages/06_smr/scripts/05_run_from_05_coloc.sh)
  Release-facing entry point that reads GWAS `.ma` files from [`05_coloc/`](/mnt/share_group_folder/tmp/work/Pre_3aQTL/Three_stages/05_coloc/README.md) and runs SMR for `stage1` to `stage3`.

- [`scripts/01_clean_eqtl.r`](/mnt/share_group_folder/tmp/work/Pre_3aQTL/Three_stages/06_smr/scripts/01_clean_eqtl.r)
  QC and deduplication helper for stage-level eQTL summary tables before BESD conversion.

- [`scripts/02_clean_gwas.r`](/mnt/share_group_folder/tmp/work/Pre_3aQTL/Three_stages/06_smr/scripts/02_clean_gwas.r)
  Legacy GWAS cleanup helper retained for reproducibility notes.

- [`scripts/03_combine_results.r`](/mnt/share_group_folder/tmp/work/Pre_3aQTL/Three_stages/06_smr/scripts/03_combine_results.r)
  Merges population- and stage-level SMR outputs into downstream summary artifacts.

- [`scripts/04_transform.py`](/mnt/share_group_folder/tmp/work/Pre_3aQTL/Three_stages/06_smr/scripts/04_transform.py)
  Legacy GWAS format conversion helper retained for historical reproducibility.

## Intended checked-in outputs

- `results/Final_Integrative_SMR_Coloc_Hits.csv`
- `results/SMR_Comprehensive_Bubble_Plot.pdf`
- `results/{EAS,EUR}/*.smr`

## Stage naming note

In this release snapshot, `stage1` to `stage3` correspond to the three primary developmental-stage QTL analyses. Files labeled `stage4` are retained because they represent the interaction-QTL follow-up layer rather than an extra developmental stage.

## Release boundary

The repository excludes regenerated heavy inputs and binary intermediates,
including:

- GWAS `.ma` files
- BESD sidecar files (`.besd`, `.esi`, `.epi`, `.summary`)
- frequency-check lists (`*.snp_failed_freq_ck.list`)
- ad hoc smoke-test directories

## Reproducibility note

This module depends on local SMR, PLINK, and reference panels, plus prepared
inputs from [`05_coloc/`](/mnt/share_group_folder/tmp/work/Pre_3aQTL/Three_stages/05_coloc/README.md).
It should be interpreted as a documented analysis layer with compact released
outputs rather than a fully self-contained portable package.
