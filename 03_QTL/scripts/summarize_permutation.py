#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd


def detect_qval_col(df: pd.DataFrame) -> str:
    for col in ("qval", "qvalue", "q_value"):
        if col in df.columns:
            return col
    raise ValueError("No q-value column found in permutation result.")


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize TensorQTL permutation results.")
    parser.add_argument("--input", required=True, help="Permutation result .txt.gz path.")
    parser.add_argument("--sig-output", required=True, help="Output path for FDR<0.05 eGenes.")
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")
    qcol = detect_qval_col(df)

    total = len(df)
    sig05 = df[df[qcol] < 0.05].copy()
    sig01 = df[df[qcol] < 0.01].copy()

    out_path = Path(args.sig_output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    sig05.to_csv(out_path, sep="\t", index=False)

    pct05 = 100.0 * len(sig05) / total if total else 0.0
    pct01 = 100.0 * len(sig01) / total if total else 0.0

    print("-" * 40)
    print("PERMUTATION QTL REPORT")
    print("-" * 40)
    print(f"Total Genes Tested:      {total}")
    print(f"eGenes (FDR < 0.05):     {len(sig05)} ({pct05:.2f}%)")
    print(f"eGenes (FDR < 0.01):     {len(sig01)} ({pct01:.2f}%)")
    print(f"Saved significant eGenes to: {out_path}")
    print("-" * 40)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
