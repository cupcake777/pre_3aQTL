#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge per-chromosome TensorQTL interaction top-association outputs.")
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    files = [Path(i) for i in args.inputs if Path(i).is_file()]
    if not files:
      raise FileNotFoundError("No interaction top-association chunk files were found.")

    merged = pd.concat([pd.read_csv(path, sep="\t", index_col=0, dtype=str) for path in files], axis=0)
    merged.to_csv(args.output, sep="\t", compression="gzip")


if __name__ == "__main__":
    main()
