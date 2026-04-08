#!/usr/bin/env python3
"""
Count Tn5 cut sites from a BAM inside a supplied BED of bins.
"""

from __future__ import annotations

import argparse
import bisect
from pathlib import Path
from typing import Dict, List, Tuple

import pysam

from extract_tn5_motifs import cut_sites_from_record


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Count Tn5 cut sites from a BAM inside custom bins."
    )
    parser.add_argument("--bam", required=True, help="Input BAM")
    parser.add_argument("--bins-bed", required=True, help="BED of bins to count into")
    parser.add_argument("--output", required=True, help="Output TSV")
    parser.add_argument("--sample", required=True, help="Sample label")
    parser.add_argument("--threads", type=int, default=2, help="BAM read threads")
    parser.add_argument("--mapq-min", type=int, default=0, help="Minimum MAPQ")
    parser.add_argument(
        "--exclude-secondary",
        action="store_true",
        help="Exclude secondary alignments",
    )
    parser.add_argument(
        "--exclude-supplementary",
        action="store_true",
        help="Exclude supplementary alignments",
    )
    return parser.parse_args()


def load_bins(path: str):
    bins_by_chrom: Dict[str, List[Tuple[int, int, List[str]]]] = {}
    starts_by_chrom: Dict[str, List[int]] = {}
    counts: Dict[str, int] = {}
    max_width = 0

    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            bin_id = fields[3]
            bins_by_chrom.setdefault(chrom, []).append((start, end, fields))
            counts[bin_id] = 0
            max_width = max(max_width, end - start)

    for chrom, entries in bins_by_chrom.items():
        entries.sort(key=lambda item: (item[0], item[1], item[2][3]))
        starts_by_chrom[chrom] = [item[0] for item in entries]

    return bins_by_chrom, starts_by_chrom, counts, max_width


def overlapping_bin_ids(
    chrom: str,
    pos: int,
    bins_by_chrom,
    starts_by_chrom,
    max_width: int,
) -> List[str]:
    entries = bins_by_chrom.get(chrom)
    if not entries:
        return []
    starts = starts_by_chrom[chrom]
    idx = bisect.bisect_right(starts, pos) - 1
    hits: List[str] = []
    while idx >= 0:
        start, end, fields = entries[idx]
        if pos - start > max_width:
            break
        if start <= pos < end:
            hits.append(fields[3])
        idx -= 1
    return hits


def write_counts(
    path: str,
    sample: str,
    bins_bed: str,
    counts: Dict[str, int],
) -> None:
    with open(bins_bed, "r", encoding="utf-8") as inp, open(
        path, "w", encoding="utf-8"
    ) as out:
        out.write(
            "sample\tchrom\tstart\tend\tbin_id\tgene_id\tgene_name\tgene_type\tstrand\ttss\ttn5_site_count\n"
        )
        for raw_line in inp:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            out.write(
                "\t".join(
                    [
                        sample,
                        fields[0],
                        fields[1],
                        fields[2],
                        fields[3],
                        fields[4],
                        fields[5],
                        fields[6],
                        fields[7],
                        fields[8],
                        str(counts[fields[3]]),
                    ]
                )
                + "\n"
            )


def main() -> None:
    args = parse_args()
    bins_by_chrom, starts_by_chrom, counts, max_width = load_bins(args.bins_bed)

    with pysam.AlignmentFile(args.bam, "rb", threads=args.threads) as bam:
        for rec in bam.fetch(until_eof=True):
            for site in cut_sites_from_record(
                rec=rec,
                mapq_min=args.mapq_min,
                exclude_secondary=args.exclude_secondary,
                exclude_supplementary=args.exclude_supplementary,
            ):
                for bin_id in overlapping_bin_ids(
                    site.chrom, site.start, bins_by_chrom, starts_by_chrom, max_width
                ):
                    counts[bin_id] += 1

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    write_counts(args.output, args.sample, args.bins_bed, counts)


if __name__ == "__main__":
    main()
