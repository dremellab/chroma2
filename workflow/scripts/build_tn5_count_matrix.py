#!/usr/bin/env python3
"""
Merge per-sample Tn5 count TSVs into a wide count matrix.

VERSION: 1.0 (2026-06-19)
  - Compatible with updated build_tn5_midpoint_bins.py that supports repeatMasker GTFs
  - Works with auto-generated feature IDs from GTFs without standard gene_id attributes
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
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


def validate_count_inputs(paths: List[str]) -> None:
    missing: List[str] = []
    not_files: List[str] = []
    empty: List[str] = []

    for raw_path in paths:
        path = Path(raw_path)
        if not path.exists():
            missing.append(raw_path)
        elif not path.is_file():
            not_files.append(raw_path)
        elif path.stat().st_size == 0:
            empty.append(raw_path)

    if not (missing or not_files or empty):
        print(
            f"Validated {len(paths)} Tn5 count input TSVs before matrix build",
            file=sys.stderr,
        )
        return

    print("Invalid Tn5 count matrix input TSVs:", file=sys.stderr)
    for label, bad_paths in (
        ("missing", missing),
        ("not regular files", not_files),
        ("empty", empty),
    ):
        if not bad_paths:
            continue
        print(f"  {label}:", file=sys.stderr)
        for bad_path in bad_paths:
            print(f"    {bad_path}", file=sys.stderr)
    raise RuntimeError("Cannot build Tn5 count matrix with invalid input TSVs")


def validate_output(path: str) -> None:
    output = Path(path)
    if not output.exists():
        raise RuntimeError(f"Output matrix TSV was not created: {path}")
    if not output.is_file():
        raise RuntimeError(f"Output matrix path is not a regular file: {path}")
    if output.stat().st_size == 0:
        raise RuntimeError(f"Output matrix TSV is empty: {path}")


def main() -> None:
    args = parse_args()
    validate_count_inputs(args.counts)
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
            row_count = 0
            for row in reader:
                row_count += 1
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
            if row_count == 0:
                # Handle empty files gracefully (e.g., when no bins overlap BAM records)
                pass

    if metadata_columns is None or not ordered_bins:
        # If no data rows found across all input files, create empty output with just headers
        # This can happen when feature regions (e.g., rRNA) are empty or don't overlap input BAM
        metadata_columns = [
            "chrom",
            "start",
            "end",
            "bin_id",
            "gene_id",
            "gene_name",
            "gene_type",
            "strand",
            "tss",
        ]
        # Write header-only output
        with open(args.output, "w", encoding="utf-8", newline="") as handle:
            fieldnames = metadata_columns + samples
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
        return

    with open(args.output, "w", encoding="utf-8", newline="") as handle:
        fieldnames = metadata_columns + samples
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for bin_id in ordered_bins:
            row = rows_by_bin[bin_id]
            for sample in samples:
                row.setdefault(sample, "0")
            writer.writerow(row)
    validate_output(args.output)
    print(
        f"Created Tn5 count matrix: {args.output} "
        f"({len(ordered_bins)} rows, {len(samples)} samples)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
