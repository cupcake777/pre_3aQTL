#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path

import pandas as pd


# Linear encoding (0/1/2) assumes the genetic effect changes *linearly* across
# stages (Stage1→2 shift equals Stage2→3 shift). This is the default interaction
# term used by TensorQTL --interaction.
#
# Limitation: trajectory clustering (02_pattern) shows that many APA events
# follow non-linear patterns (Peak, Valley). Linear encoding has low power to
# detect interactions driven by those patterns. A contrast (dummy-variable)
# encoding is written alongside the linear term (see write_covariates_and_samples)
# and can be passed to TensorQTL via --interaction-mode contrast (see 04_interaction.sh).
STAGE_ENCODING = {"Stage1": 0, "Stage2": 1, "Stage3": 2}

# Contrast encoding: two binary columns capture Stage2-vs-Stage1 and Stage3-vs-Stage1
# effects independently, allowing non-linear (e.g. Peak/Valley) pattern detection.
CONTRAST_ENCODING: dict[str, list[int]] = {
    "Stage1": [0, 0],
    "Stage2": [1, 0],
    "Stage3": [0, 1],
}


def prefix_path(prefix: Path, suffix: str) -> Path:
    return Path(f"{prefix}{suffix}")


def read_stage_covariates(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", index_col=0)


def build_combined_covariates(stage_covariates: dict[str, pd.DataFrame], stage_encoding: dict[str, int]):
    stage_names = list(stage_encoding)
    missing = [stage for stage in stage_names if stage not in stage_covariates]
    if missing:
        raise ValueError(f"Missing stage covariates for: {', '.join(missing)}")

    common_rows = set(stage_covariates[stage_names[0]].index)
    for stage in stage_names[1:]:
        common_rows &= set(stage_covariates[stage].index)
    if not common_rows:
        raise ValueError("No common covariate rows shared across stages.")

    ordered_common_rows = [idx for idx in stage_covariates[stage_names[0]].index if idx in common_rows]
    sample_order = []
    combined_parts = []
    interaction_values = []

    for stage in stage_names:
        cov_df = stage_covariates[stage].loc[ordered_common_rows]
        stage_samples = cov_df.columns.tolist()
        duplicated = set(sample_order) & set(stage_samples)
        if duplicated:
            raise ValueError(f"Duplicated sample IDs across stages: {sorted(duplicated)[:5]}")
        sample_order.extend(stage_samples)
        combined_parts.append(cov_df)
        interaction_values.extend([stage_encoding[stage]] * len(stage_samples))

    combined_covariates = pd.concat(combined_parts, axis=1)
    interaction_df = pd.DataFrame({"LifeStage": interaction_values}, index=sample_order)
    return combined_covariates, interaction_df, sample_order


def detect_chr_prefix(pvar_path: Path) -> bool:
    with pvar_path.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            chrom = line.split("\t", 1)[0]
            return chrom.startswith("chr")
    raise ValueError(f"No variant records found in {pvar_path}")


def build_phenotype_bed(expression_df: pd.DataFrame, sample_order: list[str], genotype_has_chr_prefix: bool) -> pd.DataFrame:
    missing = [sample for sample in sample_order if sample not in expression_df.columns]
    if missing:
        raise ValueError(f"Expression matrix is missing {len(missing)} samples from covariates.")

    phenotype_df = expression_df.loc[:, sample_order].copy()
    phenotype_df["pid"] = phenotype_df.index
    parsed = phenotype_df["pid"].str.split("|", expand=True, regex=False)
    if parsed.shape[1] < 4:
        raise ValueError("Phenotype IDs do not follow transcript|gene|chr:start-end|strand format.")

    loci = parsed[2].str.split(":", expand=True, regex=False)
    coords = loci[1].str.split("-", expand=True, regex=False)
    bed_df = phenotype_df.reset_index(drop=True)
    bed_df.insert(0, "end", coords[1].astype(int).to_numpy())
    bed_df.insert(0, "start", coords[0].astype(int).to_numpy())
    chrom = loci[0].astype(str)
    if genotype_has_chr_prefix:
        chrom = chrom.map(lambda x: x if x.startswith("chr") else f"chr{x}")
    else:
        chrom = chrom.str.replace("^chr", "", regex=True)
    bed_df.insert(0, "#Chr", chrom.to_numpy())
    bed_df = bed_df[["#Chr", "start", "end", "pid", *sample_order]].drop_duplicates(subset="pid")
    return bed_df.sort_values(["#Chr", "start", "end", "pid"], kind="mergesort").reset_index(drop=True)


def choose_effective_genotype_prefix(genotype_prefix: Path, extracted_genotype_prefix: Path | None = None) -> Path:
    if extracted_genotype_prefix is not None and prefix_path(extracted_genotype_prefix, ".pvar").exists():
        return extracted_genotype_prefix
    return genotype_prefix


def write_covariates_and_samples(
    stage_cov_paths: dict[str, Path],
    out_dir: Path,
) -> dict[str, Path | list[str]]:
    stage_covariates = {stage: read_stage_covariates(path) for stage, path in stage_cov_paths.items()}
    combined_covariates, interaction_df, sample_order = build_combined_covariates(stage_covariates, STAGE_ENCODING)

    out_dir.mkdir(parents=True, exist_ok=True)
    covariate_path = out_dir / "covariates_for_qtl.interaction.txt"
    interaction_path = out_dir / "interaction_term.txt"
    sample_path = out_dir / "samples.interaction.txt"

    combined_covariates.to_csv(covariate_path, sep="\t")
    interaction_df.to_csv(interaction_path, sep="\t", index_label="sample_id")
    sample_path.write_text("".join(f"{sample}\n" for sample in sample_order))

    # Also write a contrast (dummy-variable) interaction term for optional non-linear analysis.
    # Columns: Stage2_vs_Stage1 (1 if Stage2 else 0), Stage3_vs_Stage1 (1 if Stage3 else 0).
    stage_lookup: dict[str, str] = {}
    for stage_name, cov_df in zip(STAGE_ENCODING, [
        stage_covariates[s] for s in STAGE_ENCODING
    ]):
        for sample in cov_df.columns:
            stage_lookup[sample] = stage_name

    contrast_values = [CONTRAST_ENCODING[stage_lookup[s]] for s in sample_order]
    contrast_df = pd.DataFrame(
        contrast_values,
        index=sample_order,
        columns=["Stage2_vs_Stage1", "Stage3_vs_Stage1"],
    )
    contrast_path = out_dir / "interaction_term_contrast.txt"
    contrast_df.to_csv(contrast_path, sep="\t", index_label="sample_id")

    return {
        "covariates": covariate_path,
        "interaction": interaction_path,
        "interaction_contrast": contrast_path,
        "samples": sample_path,
        "sample_order": sample_order,
    }


def write_phenotype_output(
    expression_path: Path,
    sample_order: list[str],
    genotype_prefix: Path,
    out_dir: Path,
) -> Path:
    expression_df = pd.read_csv(expression_path, sep="\t", index_col=0)
    genotype_has_chr_prefix = detect_chr_prefix(prefix_path(genotype_prefix, ".pvar"))
    phenotype_bed = build_phenotype_bed(expression_df, sample_order, genotype_has_chr_prefix)
    phenotype_path = out_dir / "phenotype.interaction.bed"
    phenotype_bed.to_csv(phenotype_path, sep="\t", index=False)
    phenotype_gz = phenotype_path.with_suffix(".bed.gz")
    phenotype_tbi = phenotype_path.with_suffix(".bed.gz.tbi")
    phenotype_gz.unlink(missing_ok=True)
    phenotype_tbi.unlink(missing_ok=True)
    return phenotype_path


def run_plink_extract(plink2_bin: str, genotype_prefix: Path, sample_path: Path, out_prefix: Path) -> None:
    command = [
        plink2_bin,
        "--pfile",
        str(genotype_prefix),
        "--keep",
        str(sample_path),
        "--make-pgen",
        "--out",
        str(out_prefix),
    ]
    subprocess.run(command, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare merged TensorQTL interaction inputs using LifeStage 0/1/2.")
    parser.add_argument("--expression", default="../01_pre/after.combat.txt")
    parser.add_argument("--cov-dir", default="../01_pre")
    parser.add_argument("--genotype-prefix", default="../../raw_data/genotype")
    parser.add_argument("--out-dir", default="interaction/pre_input")
    parser.add_argument("--plink2", default="plink2")
    parser.add_argument("--skip-genotype", action="store_true")
    args = parser.parse_args()

    cov_dir = Path(args.cov_dir)
    stage_cov_paths = {
        stage: cov_dir / f"covariates_for_qtl.{stage}.txt"
        for stage in STAGE_ENCODING
    }
    outputs = write_covariates_and_samples(
        stage_cov_paths=stage_cov_paths,
        out_dir=Path(args.out_dir),
    )
    extracted_prefix = Path(args.out_dir) / "genotype.interaction"

    if not args.skip_genotype:
        run_plink_extract(
            plink2_bin=args.plink2,
            genotype_prefix=Path(args.genotype_prefix),
            sample_path=outputs["samples"],
            out_prefix=extracted_prefix,
        )

    effective_genotype_prefix = choose_effective_genotype_prefix(Path(args.genotype_prefix), extracted_prefix)
    phenotype_path = write_phenotype_output(
        expression_path=Path(args.expression),
        sample_order=outputs["sample_order"],
        genotype_prefix=effective_genotype_prefix,
        out_dir=Path(args.out_dir),
    )
    outputs["phenotype_bed"] = phenotype_path


if __name__ == "__main__":
    main()
