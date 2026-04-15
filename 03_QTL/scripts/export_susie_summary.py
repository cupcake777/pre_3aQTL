#!/usr/bin/env python3
import argparse
import pickle
from pathlib import Path

import numpy as np
import pandas as pd


def parse_phenotype(phenotype_id: str) -> dict[str, object]:
    parts = str(phenotype_id).split("|")
    region = parts[2] if len(parts) > 2 else ""
    chrom = ""
    start = np.nan
    end = np.nan
    if ":" in region and "-" in region:
        chrom, coords = region.split(":", 1)
        start_s, end_s = coords.split("-", 1)
        start = float(start_s)
        end = float(end_s)
    return {
        "transcript": parts[0] if len(parts) > 0 else "",
        "gene": parts[1] if len(parts) > 1 else "",
        "chrom": chrom,
        "start": start,
        "end": end,
        "strand": parts[3] if len(parts) > 3 else "",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Export a compact text summary from TensorQTL SuSiE pickle output.")
    parser.add_argument("--pickle", required=True, help="Input .SuSiE.pickle file.")
    parser.add_argument("--output", required=True, help="Output .txt.gz summary path.")
    args = parser.parse_args()

    with open(args.pickle, "rb") as handle:
        data = pickle.load(handle)

    rows: list[dict[str, object]] = []
    for phenotype_id, result in data.items():
        info = parse_phenotype(phenotype_id)
        pip_df = result.get("pip")
        sets = result.get("sets", {})
        cs_map = sets.get("cs", {})
        purity = sets.get("purity")
        if pip_df is None:
            continue

        for cs_label, indices in cs_map.items():
            cs_pips = pip_df["pip"].iloc[indices]
            top_idx = cs_pips.idxmax()
            pur_row = purity.loc[cs_label] if purity is not None and cs_label in purity.index else None
            rows.append(
                {
                    "phenotype_id": phenotype_id,
                    "gene": info["gene"],
                    "chrom": info["chrom"],
                    "start": info["start"],
                    "end": info["end"],
                    "cs_label": cs_label,
                    "cs_size": len(indices),
                    "top_variant": top_idx,
                    "top_pip": float(cs_pips.max()),
                    "sum_pip": float(cs_pips.sum()),
                    "min_abs_corr": np.nan if pur_row is None else pur_row.get("min_abs_corr", np.nan),
                    "mean_abs_corr": np.nan if pur_row is None else pur_row.get("mean_abs_corr", np.nan),
                    "converged": bool(result.get("converged", False)),
                    "n_cs": len(cs_map),
                }
            )

    df = pd.DataFrame(rows)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False, compression="gzip", float_format="%.6g")
    print(f"[INFO] Exported {len(df):,} credible-set summary rows to {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
