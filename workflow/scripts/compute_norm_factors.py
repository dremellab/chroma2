#!/usr/bin/env python3
"""Compute read-depth-based normalization factors from idxstats_summary.tsv.

For each sample:
    non_mt_trimmed   = post_trimming_input_sequences - chrM_raw   (host chrM only)
    million_non_mt   = non_mt_trimmed / 1e6
    bw_scale_factor  = 1 / million_non_mt
    deseq2_size_factor = bw_scale_factor / mean(bw_scale_factor across all samples)

Writes a TSV with columns: sample, bw_scale_factor, deseq2_size_factor.
"""

from __future__ import annotations

import argparse
import csv
import statistics
import sys


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--idxstats-summary", required=True, help="Path to idxstats_summary.tsv"
    )
    parser.add_argument(
        "--output", required=True, help="Path to write norm_factors.tsv"
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    samples = []
    bw_scale_factors = {}
    with open(args.idxstats_summary, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sample = row["sample"]
            post_trimming = row["post_trimming_input_sequences"]
            chrm_raw = row["chrM_raw"]
            if post_trimming == "NA" or chrm_raw == "NA":
                sys.exit(
                    f"ERROR: sample '{sample}' has missing post_trimming_input_sequences "
                    f"or chrM_raw in {args.idxstats_summary}; cannot compute a read-depth "
                    "normalization factor for it."
                )
            non_mt_trimmed = int(post_trimming) - int(chrm_raw)
            if non_mt_trimmed <= 0:
                sys.exit(
                    f"ERROR: sample '{sample}' has non-positive non-MT trimmed read count "
                    f"({non_mt_trimmed}); cannot compute a normalization factor for it."
                )
            million_non_mt = non_mt_trimmed / 1e6
            samples.append(sample)
            bw_scale_factors[sample] = 1.0 / million_non_mt

    if not samples:
        sys.exit(f"ERROR: no samples found in {args.idxstats_summary}")

    mean_scale_factor = statistics.mean(bw_scale_factors.values())

    with open(args.output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "bw_scale_factor", "deseq2_size_factor"])
        for sample in samples:
            bw_scale_factor = bw_scale_factors[sample]
            deseq2_size_factor = bw_scale_factor / mean_scale_factor
            writer.writerow(
                [sample, f"{bw_scale_factor:.10g}", f"{deseq2_size_factor:.10g}"]
            )


if __name__ == "__main__":
    main()
