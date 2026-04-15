#!/usr/bin/env python3
import argparse
import glob
import gzip
import sys
from pathlib import Path

import pandas as pd


def main() -> int:
    parser = argparse.ArgumentParser(description="Merge TensorQTL cis_nominal parquet shards.")
    parser.add_argument("--prefix", required=True, help="TensorQTL output prefix before parquet shard suffixes.")
    parser.add_argument("--output", required=True, help="Output .txt.gz path.")
    args = parser.parse_args()

    pattern = f"{args.prefix}.cis_qtl_pairs.*.parquet"
    files = sorted(glob.glob(pattern))
    if not files:
        print(f"[ERROR] No parquet files matched {pattern}", file=sys.stderr)
        return 1

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    total_rows = 0
    unique_pheno: set[str] = set()
    append_mode = out_path.exists() and out_path.stat().st_size > 0
    open_mode = "at" if append_mode else "wt"

    with gzip.open(out_path, open_mode) as handle:
        for idx, path in enumerate(files):
            df = pd.read_parquet(path)
            total_rows += len(df)
            if "phenotype_id" in df.columns:
                unique_pheno.update(df["phenotype_id"].astype(str).unique().tolist())
            df.to_csv(
                handle,
                sep="\t",
                index=False,
                header=(idx == 0 and not append_mode),
                float_format="%.6g",
            )
            Path(path).unlink(missing_ok=True)
            print(f"[INFO] Appended shard {idx + 1}/{len(files)}: {path}", flush=True)

    print(f"[INFO] Merged {len(files)} shards into {out_path}", flush=True)
    print(f"[INFO] Total associations: {total_rows:,}", flush=True)
    print(f"[INFO] Unique phenotypes: {len(unique_pheno):,}", flush=True)
    print("[INFO] Removed temporary parquet shards", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
