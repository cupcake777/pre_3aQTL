#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd


def _resolve_interaction_columns(df: pd.DataFrame, interaction_name: str) -> tuple[str, str | None, str]:
    modern_beta_col = "b_gi"
    modern_beta_se_col = "b_gi_se"
    modern_pval_col = "pval_gi"
    if modern_beta_col in df.columns and modern_pval_col in df.columns:
        beta_se_col = modern_beta_se_col if modern_beta_se_col in df.columns else None
        return modern_beta_col, beta_se_col, modern_pval_col

    legacy_beta_col = f"b_g-{interaction_name}"
    legacy_beta_se_col = f"b_g-{interaction_name}_se"
    legacy_pval_col = f"pval_g-{interaction_name}"
    if legacy_beta_col in df.columns and legacy_pval_col in df.columns:
        beta_se_col = legacy_beta_se_col if legacy_beta_se_col in df.columns else None
        return legacy_beta_col, beta_se_col, legacy_pval_col

    raise ValueError(f"Expected interaction columns for {interaction_name} were not found.")


def summarize_interaction_top_assoc(input_path: str | Path, output_path: str | Path, fdr: float = 0.05, interaction_name: str = "LifeStage") -> None:
    df = pd.read_csv(input_path, sep="\t")
    pval_col = "pval_adj_bh" if "pval_adj_bh" in df.columns else "pval_emt"
    if pval_col not in df.columns:
        raise ValueError("TensorQTL interaction output is missing both pval_adj_bh and pval_emt.")

    beta_col, beta_se_col, pval_gi_col = _resolve_interaction_columns(df, interaction_name)

    sig_df = df.loc[df[pval_col] <= fdr].copy()
    sig_df["slope"] = sig_df[beta_col]
    if beta_se_col is not None:
        sig_df["slope_se"] = sig_df[beta_se_col]
    sig_df["interaction_pval"] = sig_df[pval_gi_col]
    sig_df.to_csv(output_path, sep="\t", index=False, compression="gzip")


def main() -> None:
    parser = argparse.ArgumentParser(description="Filter TensorQTL interaction top associations into workflow-compatible output.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--interaction-name", default="LifeStage")
    args = parser.parse_args()

    summarize_interaction_top_assoc(
        input_path=args.input,
        output_path=args.output,
        fdr=args.fdr,
        interaction_name=args.interaction_name,
    )


if __name__ == "__main__":
    main()
