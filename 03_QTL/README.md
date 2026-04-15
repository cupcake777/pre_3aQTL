# Three-stage 3'aQTL Analysis Pipeline

This directory contains a modular five-step workflow for stage-specific 3'aQTL mapping, downsampling-based harmonization, cross-stage and external validation, interaction QTL analysis, and final integrative comparison. The workflow is designed so that each analytical layer has a single entry script and a stable set of expected outputs.

## Workflow Overview

1. `01_stage_QTL.sh`
   Generate stage-specific `nominal`, `permutation`, `SuSiE`, and significant QTL-pair outputs for `Stage1`, `Stage2`, and `Stage3` using the processed inputs from `../01_pre`.

2. `02_downsample.sh`
   Rebuild harmonized downsampled inputs and rerun the same QTL workflow so that cross-stage comparisons are not driven by sample-size differences.

3. `03_compare.sh`
   Compare downsampled stage-specific results across developmental stages and against external reference 3'aQTL datasets.

4. `04_interaction.sh`
   Run interaction QTL analysis on interaction-ready inputs and generate basic visualization outputs.

5. `05_integrate.sh`
   Compare stage-specific and interaction-derived signals to identify shared and distinct regulatory patterns.

The intended execution order is:

```bash
bash 01_stage_QTL.sh
bash 02_downsample.sh
bash 03_compare.sh
bash 04_interaction.sh
bash 05_integrate.sh
```

Each script is also designed to be runnable independently once its upstream inputs already exist.

## Directory Structure

- `nominal/`
  Stage-specific nominal TensorQTL results from the original sample sets.

- `permutation/`
  Stage-specific permutation-based TensorQTL results and significant eGene summaries.

- `susie/`
  Stage-specific SuSiE fine-mapping summaries and optional visualizations.

- `downsample/`
  Downsampled inputs and outputs, including harmonized stage-specific genotype/phenotype/covariate files and rerun TensorQTL results.

- `compare/`
  Cross-stage and external validation analyses, including overlap tables, replication summaries, and figures.

- `interaction/`
  Interaction-QTL inputs, outputs, and downstream visualizations.

- `integrate/`
  Final comparison outputs between stage-specific and interaction-derived results.

- `scripts/`
  Reusable helper scripts shared across the five top-level workflow entry points.

## Script-level Contract

### `01_stage_QTL.sh`

Inputs:

- `../01_pre/genotype.stage1`, `../01_pre/genotype.stage2`, `../01_pre/genotype.stage3`
- `../01_pre/phenotype.Stage1.bed.gz`, `../01_pre/phenotype.Stage2.bed.gz`, `../01_pre/phenotype.Stage3.bed.gz`
- `../01_pre/covariates_for_qtl.Stage1.txt`, `../01_pre/covariates_for_qtl.Stage2.txt`, `../01_pre/covariates_for_qtl.Stage3.txt`

Outputs:

- `nominal/stage1.cis_qtl_pairs.txt.gz`, `nominal/stage2.cis_qtl_pairs.txt.gz`, `nominal/stage3.cis_qtl_pairs.txt.gz`
- `permutation/stage1.cis_qtl.txt.gz`, `permutation/stage2.cis_qtl.txt.gz`, `permutation/stage3.cis_qtl.txt.gz`
- `permutation/stage1.sig_egenes.fdr05.txt`, `permutation/stage2.sig_egenes.fdr05.txt`, `permutation/stage3.sig_egenes.fdr05.txt`
- `susie/stage1.SuSiE_summary.txt.gz`, `susie/stage2.SuSiE_summary.txt.gz`, `susie/stage3.SuSiE_summary.txt.gz`
- `stage1_sig_QTL.txt.gz`, `stage2_sig_QTL.txt.gz`, `stage3_sig_QTL.txt.gz`

### `02_downsample.sh`

Inputs:

- Original stage-specific inputs from `../01_pre`
- Downsampling targets inferred from stage sample counts

Outputs:

- `downsample/pre_input/sample.Stage1.priority`, `downsample/pre_input/sample.Stage2.priority`, `downsample/pre_input/sample.Stage3.priority`
- `downsample/pre_input/genotype.stage1.priority.*`, `downsample/pre_input/genotype.stage2.priority.*`, `downsample/pre_input/genotype.stage3.priority.*`
- `downsample/nominal/stage1.cis_qtl_pairs.txt.gz`, `downsample/nominal/stage2.cis_qtl_pairs.txt.gz`, `downsample/nominal/stage3.cis_qtl_pairs.txt.gz`
- `downsample/permutation/stage1.cis_qtl.txt.gz`, `downsample/permutation/stage2.cis_qtl.txt.gz`, `downsample/permutation/stage3.cis_qtl.txt.gz`
- `downsample/susie/stage1.SuSiE_summary.txt.gz`, `downsample/susie/stage2.SuSiE_summary.txt.gz`, `downsample/susie/stage3.SuSiE_summary.txt.gz`
- `downsample/stage1_sig_QTL.txt.gz`, `downsample/stage2_sig_QTL.txt.gz`, `downsample/stage3_sig_QTL.txt.gz`

### `03_compare.sh`

Inputs:

- `downsample/stage*_sig_QTL.txt.gz`
- `downsample/nominal/stage*.cis_qtl_pairs.txt.gz`
- External reference datasets under `compare/*.3aQTL.txt.gz`

Outputs:

- `compare/table_stage_overlap_summary.txt`
- `compare/table_stage_correlation_summary.txt`
- `compare/figure_stage_overlap_upset.pdf`
- `compare/figure_stage_effect_size_correlation.pdf`
- Optional: `compare/table_reference_pi1_summary.txt`
- Optional: `compare/figure_reference_pi1_heatmap.pdf`

### `04_interaction.sh`

Inputs:

- Interaction-ready genotype, phenotype, and covariate data under `interaction/pre_input/`
- An interaction term encoded in the covariate model

Outputs:

- `interaction/nominal/*`
- `interaction/interaction_sig_QTL.txt.gz`
- `interaction/figure_interaction_effect_size_histogram.pdf`
- Prepared inputs retained under `interaction/pre_input/`

### `05_integrate.sh`

Inputs:

- Stage-specific significant QTL outputs
- Interaction significant QTL outputs
- Optional nominal and SuSiE summaries for deeper comparison

Outputs:

- `integrate/table_stage_vs_interaction_summary.txt`
- `integrate/figure_stage_vs_interaction_scatter.pdf`

## Design Principles

- Stage-specific discovery and cross-stage comparison are intentionally separated.
- Cross-stage comparison is based on downsampled data to minimize sample-size-driven bias.
- Interaction analysis is treated as a distinct analytical layer rather than being merged into the primary stage-specific workflow.
- Final integration is reserved for interpretation after both stage-specific and interaction analyses are complete.

## Current Status

This directory is the most complete module in the current `Three_stages` snapshot and can be read as the public-facing contract for the stage-QTL workflow.

The workflow structure, entry scripts, and expected outputs are already stabilized enough for repository sharing. However, the module still has two practical limitations that should be stated explicitly:

- some scripts depend on local reference resources under `../../raw_data/` or site-specific software environments;
- run-time logs and caches may still be generated locally and should not be treated as publication artifacts.

Result files use a simple naming convention:

- `table_*.txt` for plain-text summary tables;
- `figure_*.pdf` for manuscript-facing plots;
- script files retain ordered numeric prefixes such as `01_*.R` when execution order matters.

For these reasons, the directory is suitable for manuscript-linked code release and method inspection, but not yet a fully self-contained reproducible package.
