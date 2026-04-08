#!/usr/bin/env python3
"""Merge per-sample Tn5 count TSVs into a wide count matrix."""

from __future__ import annotations

import argparse
import csv
from typing import Dict, List


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a wide Tn5 count matrix from per-sample count TSVs."
    )
    parser.add_argument(
        "--counts",
        nargs="+",
        required=True,
        help="Per-sample count TSVs",
    )
    parser.add_argument("--output", required=True, help="Output matrix TSV")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows_by_bin: Dict[str, Dict[str, str]] = {}
    samples: List[str] = []
    ordered_bins: List[str] = []
    metadata_columns: List[str] | None = None

    for path in args.counts:
        with open(path, "r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                continue
            current_metadata_columns = [
                field
                for field in reader.fieldnames
                if field not in {"sample", "tn5_site_count"}
            ]
            if metadata_columns is None:
                metadata_columns = current_metadata_columns
            elif metadata_columns != current_metadata_columns:
                raise ValueError(
                    f"Inconsistent count TSV schema in {path}: "
                    f"expected {metadata_columns}, found {current_metadata_columns}"
                )
            sample_name = None
            for row in reader:
                if sample_name is None:
                    sample_name = row["sample"]
                    if sample_name not in samples:
                        samples.append(sample_name)
                bin_id = row["bin_id"]
                if bin_id not in rows_by_bin:
                    rows_by_bin[bin_id] = {
                        column: row[column] for column in metadata_columns
                    }
                    ordered_bins.append(bin_id)
                rows_by_bin[bin_id][sample_name] = row["tn5_site_count"]

    if metadata_columns is None:
        raise ValueError("No count rows were found to build the matrix")

    with open(args.output, "w", encoding="utf-8", newline="") as handle:
        fieldnames = metadata_columns + samples
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for bin_id in ordered_bins:
            row = rows_by_bin[bin_id]
            for sample in samples:
                row.setdefault(sample, "0")
            writer.writerow(row)


if __name__ == "__main__":
    main()
