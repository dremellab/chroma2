#!/usr/bin/env python3
"""
Build fixed-size non-overlapping Tn5 counting bins across viral genomes.

Creates uniform bins across each viral genome/contig with a configurable bin size.
One bin per fixed-size interval, with the last bin potentially smaller.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create fixed-size non-overlapping Tn5 counting bins across viral genome(s). "
            "Useful for coverage analysis of complete viral genomes."
        )
    )
    parser.add_argument(
        "--chromsizes", required=True, help="Chromsizes TSV (viral genome)"
    )
    parser.add_argument(
        "--bin-size",
        type=int,
        default=200,
        help="Fixed bin size in bp (default: 200). Last bin may be smaller.",
    )
    parser.add_argument("--output", required=True, help="Output BED path")
    return parser.parse_args()


def load_chromsizes(path: str) -> Dict[str, int]:
    chromsizes: Dict[str, int] = {}
    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            chrom, size = line.split("\t", 1)
            chromsizes[chrom] = int(size)
    return chromsizes


def write_bins(path: str, chromsizes: Dict[str, int], bin_size: int) -> None:
    rows = []

    for chrom in sorted(chromsizes.keys()):
        chrom_size = chromsizes[chrom]
        start = 0
        bin_num = 1

        while start < chrom_size:
            end = min(start + bin_size, chrom_size)
            bin_id = f"{chrom}:{start}-{end}"

            rows.append(
                (
                    chrom,
                    str(start),
                    str(end),
                    bin_id,
                    str(bin_num),
                    str(end - start),  # actual bin length
                )
            )

            start = end
            bin_num += 1

    with open(path, "w", encoding="utf-8") as out:
        for row in rows:
            out.write("\t".join(row) + "\n")


def main() -> None:
    args = parse_args()

    if args.bin_size <= 0:
        raise ValueError("--bin-size must be positive")

    chromsizes = load_chromsizes(args.chromsizes)

    if not chromsizes:
        raise ValueError(f"No chromosomes found in {args.chromsizes}")

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    write_bins(args.output, chromsizes, args.bin_size)

    total_bins = sum(
        (size + args.bin_size - 1) // args.bin_size for size in chromsizes.values()
    )

    print(
        f"✓ Created {total_bins} bins ({args.bin_size}bp) across {len(chromsizes)} contig(s) in {args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
