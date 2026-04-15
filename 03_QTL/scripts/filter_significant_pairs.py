#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd


def detect_column(df: pd.DataFrame, candidates: list[str], label: str) -> str:
    for col in candidates:
        if col in df.columns:
            return col
    raise ValueError(f"Missing required column for {label}: {candidates}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Filter nominal QTL pairs using permutation-derived gene thresholds.")
    parser.add_argument("--permutation", required=True, help="Permutation cis result .txt.gz path.")
    parser.add_argument("--nominal", required=True, help="Nominal cis pair result .txt.gz path.")
    parser.add_argument("--output", required=True, help="Output .txt.gz path for significant SNP-gene pairs.")
    parser.add_argument("--qval-threshold", type=float, default=0.05, help="FDR threshold used to retain genes.")
    args = parser.parse_args()

    perm = pd.read_csv(args.permutation, sep="\t")
    perm_pheno_col = detect_column(perm, ["phenotype_id", "pid"], "permutation phenotype")
    qcol = detect_column(perm, ["qval", "qvalue", "q_value"], "q-value")
    thrcol = detect_column(perm, ["pval_nominal_threshold"], "nominal threshold")

    sig_gene_thresholds = (
        perm.loc[perm[qcol] <= args.qval_threshold, [perm_pheno_col, thrcol]]
        .drop_duplicates()
        .set_index(perm_pheno_col)[thrcol]
        .to_dict()
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    sig_pairs = 0
    sig_genes: set[str] = set()
    first_chunk = True
    nom_pheno_col = None

    for chunk in pd.read_csv(args.nominal, sep="\t", chunksize=500000):
        if nom_pheno_col is None:
            nom_pheno_col = detect_column(chunk, ["phenotype_id", "pid"], "nominal phenotype")
            pval_col = detect_column(chunk, ["pval_nominal", "npval"], "nominal p-value")

        chunk = chunk[chunk[nom_pheno_col].isin(sig_gene_thresholds)]
        if chunk.empty:
            continue

        thresholds = chunk[nom_pheno_col].map(sig_gene_thresholds)
        filtered = chunk.loc[chunk[pval_col] <= thresholds].copy()
        if filtered.empty:
            continue

        filtered.to_csv(
            out_path,
            sep="\t",
            index=False,
            compression="gzip",
            float_format="%.6g",
            mode="wt" if first_chunk else "at",
            header=first_chunk,
        )
        first_chunk = False
        sig_pairs += len(filtered)
        sig_genes.update(filtered[nom_pheno_col].astype(str).unique().tolist())

    if first_chunk:
        raise ValueError("No significant SNP-gene pairs passed the filtering criteria.")

    print(f"[INFO] Significant SNP-gene pairs: {sig_pairs:,}")
    print(f"[INFO] Significant phenotypes: {len(sig_genes):,}")
    print(f"[INFO] Saved filtered pairs to: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
