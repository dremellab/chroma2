#!/usr/bin/env python3
"""
Count Tn5 cut sites from a BAM inside a supplied BED of bins.
"""

from __future__ import annotations

import argparse
import bisect
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pysam

from extract_tn5_motifs import cut_sites_from_record


def _non_negative_int(value: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError(f"must be >= 0, got {parsed}")
    return parsed


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Count Tn5 cut sites from a BAM inside custom bins."
    )
    parser.add_argument("--bam", required=True, help="Input BAM")
    parser.add_argument("--bins-bed", required=True, help="BED of bins to count into")
    parser.add_argument("--output", required=True, help="Output TSV")
    parser.add_argument("--sample", required=True, help="Sample label")
    parser.add_argument("--threads", type=int, default=2, help="BAM read threads")
    parser.add_argument(
        "--mapq-min",
        type=_non_negative_int,
        default=0,
        help="Minimum MAPQ (must be >= 0)",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=1_000_000,
        help="Report progress every N BAM records processed (default: 1000000)",
    )
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
    parser.add_argument(
        "--fractional-counting",
        action="store_true",
        help="Use NH-weighted fractional counting for multi-mapping reads (1/NH per alignment)",
    )
    return parser.parse_args(argv)


def log(message: str) -> None:
    timestamp = datetime.now().isoformat(timespec="seconds")
    print(f"[{timestamp}] {message}", file=sys.stderr, flush=True)


def load_bins(path: str):
    bins_by_chrom: Dict[str, List[Tuple[int, int, List[str]]]] = {}
    starts_by_chrom: Dict[str, List[int]] = {}
    counts: Dict[str, float] = {}
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


def get_nh_value(rec: pysam.AlignedSegment) -> int:
    """
    Safely extract NH tag (number of reported alignments) from BAM record.
    Returns 1 if NH tag is missing (single alignment).
    """
    try:
        return rec.get_tag("NH")
    except KeyError:
        return 1


def write_counts(
    path: str,
    sample: str,
    bins_bed: str,
    counts: Dict[str, float],
) -> int:
    rows_written = 0
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
            rounded_count = round(counts[fields[3]])
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
                        str(rounded_count),
                    ]
                )
                + "\n"
            )
            rows_written += 1
    return rows_written


def count_lines(path: Path) -> int:
    with open(path, "r", encoding="utf-8") as handle:
        return sum(1 for _ in handle)


def validate_output(path: str, rows_written: int, expected_rows: int) -> int:
    output_path = Path(path)
    if not output_path.exists():
        raise RuntimeError(f"Output TSV was not created: {path}")
    output_size = output_path.stat().st_size
    if output_size == 0:
        raise RuntimeError(f"Output TSV is empty: {path}")
    if rows_written != expected_rows:
        raise RuntimeError(
            f"Output TSV row count mismatch while writing {path}: "
            f"wrote {rows_written} data rows, expected {expected_rows}"
        )
    line_count = count_lines(output_path)
    expected_lines = expected_rows + 1
    if line_count != expected_lines:
        raise RuntimeError(
            f"Output TSV line count mismatch for {path}: "
            f"found {line_count} lines, expected {expected_lines} "
            f"({expected_rows} data rows plus header)"
        )
    return output_size


def main() -> None:
    args = parse_args()
    log(f"Starting Tn5 bin counting for sample {args.sample}")
    log(f"Input BAM: {args.bam}")
    log(f"Input bins: {args.bins_bed}")
    log(f"Output TSV: {args.output}")
    log(
        f"Counting mode: {'NH-weighted fractional' if args.fractional_counting else 'integer'}"
    )
    bins_by_chrom, starts_by_chrom, counts, max_width = load_bins(args.bins_bed)
    log(f"Loaded {len(counts)} bins across {len(bins_by_chrom)} contigs")

    records_processed = 0
    cut_sites_seen = 0
    overlapping_sites = 0
    with pysam.AlignmentFile(args.bam, "rb", threads=args.threads) as bam:
        for rec in bam.fetch(until_eof=True):
            records_processed += 1
            for site in cut_sites_from_record(
                rec=rec,
                mapq_min=args.mapq_min,
                exclude_secondary=args.exclude_secondary,
                exclude_supplementary=args.exclude_supplementary,
            ):
                cut_sites_seen += 1
                weight = 1.0 / get_nh_value(rec) if args.fractional_counting else 1.0
                for bin_id in overlapping_bin_ids(
                    site.chrom, site.start, bins_by_chrom, starts_by_chrom, max_width
                ):
                    counts[bin_id] += weight
                    overlapping_sites += 1
            if args.progress_every > 0 and records_processed % args.progress_every == 0:
                log(
                    "Processed "
                    f"{records_processed} BAM records, observed {cut_sites_seen} "
                    f"cut sites, assigned {overlapping_sites} bin overlaps"
                )

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    rows_written = write_counts(args.output, args.sample, args.bins_bed, counts)
    output_size = validate_output(args.output, rows_written, len(counts))
    log(
        "Finished Tn5 bin counting: "
        f"processed {records_processed} BAM records, observed {cut_sites_seen} "
        f"cut sites, assigned {overlapping_sites} bin overlaps"
    )
    log(
        "Created output TSV: "
        f"{args.output} ({rows_written} data rows, {output_size} bytes)"
    )


if __name__ == "__main__":
    main()
