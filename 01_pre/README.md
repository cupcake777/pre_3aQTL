# 01_pre

`01_pre/` contains preprocessing scripts and the derived files required by [`02_pattern/`](/mnt/share_group_folder/tmp/work/Pre_3aQTL/Three_stages/02_pattern/README.md) and [`03_QTL/`](/mnt/share_group_folder/tmp/work/Pre_3aQTL/Three_stages/03_QTL/README.md).

## Primary scripts

- [`01_preprocess_all.r`](/mnt/share_group_folder/tmp/work/Pre_3aQTL/Three_stages/01_pre/01_preprocess_all.r)
  Main preprocessing entry point from raw APA counts, metadata, and genotype input.

- [`02_preprocess_batch9.r`](/mnt/share_group_folder/tmp/work/Pre_3aQTL/Three_stages/01_pre/02_preprocess_batch9.r)
  Batch-specific preprocessing for batch 9 samples.

- [`03_evaluate.r`](/mnt/share_group_folder/tmp/work/Pre_3aQTL/Three_stages/01_pre/03_evaluate.r)
  Post-processing evaluation and QC.

- [`04_summary.r`](/mnt/share_group_folder/tmp/work/Pre_3aQTL/Three_stages/01_pre/04_summary.r)
  Summary reporting for the preprocessing stage.

## Key outputs

- `covariates_for_qtl.Stage{1,2,3}.txt`
- `sample.Stage{1,2,3}.final`
- `phenotype.Stage{1,2,3}.bed.gz`
- `outlier_detect_mahalanobis.png`

## Release boundary

Large regenerated intermediates, including expression matrices and PLINK binaries, are excluded from standard version control. The top-level `.gitignore` captures these files explicitly.

## Reproducibility note

The scripts assume local upstream resources under `../../raw_data/`. This module is therefore released as a documented workflow entry point rather than a fully self-contained package.
