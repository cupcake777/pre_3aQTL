# downsample

`downsample/` standardizes the three developmental stages to a more comparable sample size before rerunning the QTL workflow. The purpose is to reduce sample-size-driven bias in cross-stage comparisons.

## Primary script

- [`01_prepare_downsample_inputs.R`](/mnt/share_group_folder/tmp/work/Pre_3aQTL/Three_stages/03_QTL/downsample/01_prepare_downsample_inputs.R)
  Generates the downsampled sample lists, covariates, and input matrices.

## Scope

This directory is not an independent analysis theme. It is the normalization layer for `03_QTL` and should be interpreted together with `02_downsample.sh` in the parent module.
