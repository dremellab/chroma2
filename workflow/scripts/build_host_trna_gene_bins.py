#!/usr/bin/env python3
"""
Build fixed-width host Tn5 counting bins centered on tRNA gene body midpoints.

Extracts all genes/transcripts from a tRNA-specific GTF file and creates bins
centered on the gene body midpoint (not TSS). One bin per tRNA gene, filtered
to host contigs only.

Output columns:

+-----+-----------+-----------------------------------------------------+
| Col | Name      | Description                                         |
+=====+===========+=====================================================+
| 1   | chrom     | Host contig/chromosome name.                       |
| 2   | start     | 0-based bin start, clipped at 0.                   |
| 3   | end       | BED-style end, clipped at chromosome length.       |
| 4   | bin_id    | Synthetic ID: "{gene_id}|{gene_type}".             |
| 5   | gene_id   | GTF gene/transcript identifier.                    |
| 6   | gene_name | GTF gene name, or gene_id if no name is present.   |
| 7   | gene_type | GTF gene_type attribute, defaulting to "tRNA".     |
| 8   | strand    | GTF strand value.                                  |
| 9   | center    | 0-based midpoint used to create the bin.           |
+-----+-----------+-----------------------------------------------------+
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path
from typing import Dict, Iterable, Iterator, Optional


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create fixed-width host Tn5 counting bins centered on tRNA gene body midpoints."
        )
    )
    parser.add_argument("--trna-gtf", required=True, help="tRNA GTF file")
    parser.add_argument("--host-regions", required=True, help="Host .fa.regions file")
    parser.add_argument("--chromsizes", required=True, help="Host chromsizes TSV")
    parser.add_argument(
        "--flank-size",
        type=int,
        default=100,
        help="Flank size around the gene center in bp (default: 100)",
    )
    parser.add_argument("--output", required=True, help="Output BED path")
    return parser.parse_args()


def parse_attributes(attr_text: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for raw_field in attr_text.strip().strip(";").split(";"):
        field = raw_field.strip()
        if not field:
            continue
        if "=" in field and '"' not in field:
            key, value = field.split("=", 1)
        elif " " in field:
            key, value = field.split(" ", 1)
        else:
            continue
        attrs[key.strip()] = value.strip().strip('"')
    return attrs


def load_host_contigs(path: str) -> set[str]:
    contigs: set[str] = set()
    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                continue
            contigs.update(contig for contig in fields[1].split() if contig)
    return contigs


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


def gene_type(attrs: Dict[str, str]) -> str:
    return attrs.get("gene_type", "tRNA")


def feature_name(attrs: Dict[str, str], fallback: str) -> str:
    for key in ("gene_name", "Name", "gene", "transcript_name"):
        value = attrs.get(key)
        if value:
            return value
    return fallback


def feature_id(attrs: Dict[str, str]) -> Optional[str]:
    for key in ("gene_id", "ID", "transcript_id"):
        value = attrs.get(key)
        if value:
            return value
    return None


def gene_center_from_fields(fields: list[str]) -> int:
    start_0based = int(fields[3]) - 1
    end_0based = int(fields[4]) - 1  # GTF end is 1-based inclusive → 0-based inclusive
    return (start_0based + end_0based) // 2


def iter_gtf_records(path: str) -> Iterator[list[str]]:
    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) >= 9:
                yield fields


def collect_trna_genes(
    gtf_path: str,
    host_contigs: set[str],
    include_types: Optional[Iterable[str]] = None,
) -> Dict[str, Dict[str, object]]:
    include = {item.lower() for item in include_types} if include_types else None
    genes: Dict[str, Dict[str, object]] = {}

    for fields in iter_gtf_records(gtf_path):
        feature = fields[2]
        if feature not in {"transcript", "gene"}:
            continue
        chrom = fields[0]
        if chrom not in host_contigs:
            continue
        attrs = parse_attributes(fields[8])
        gid = feature_id(attrs)
        if not gid:
            continue
        gtype = gene_type(attrs)
        if include is not None and gtype.lower() not in include:
            continue

        strand = fields[6]
        center = gene_center_from_fields(fields)
        name = feature_name(attrs, gid)

        existing = genes.get(gid)
        if existing is None:
            genes[gid] = {
                "chrom": chrom,
                "strand": strand,
                "center": center,
                "gene_name": name,
                "gene_type": gtype,
            }
            continue

        if existing["chrom"] != chrom or existing["strand"] != strand:
            continue
        # keep the outermost span → recompute center from min start / max end
        existing_center = int(existing["center"])
        existing["center"] = (existing_center + center) // 2

    return genes


def main() -> None:
    args = parse_args()
    if args.flank_size < 0 or args.flank_size > 1000:
        raise ValueError("--flank-size must be between 0 and 1000 bp")

    host_contigs = load_host_contigs(args.host_regions)
    chromsizes = load_chromsizes(args.chromsizes)

    trna_genes = collect_trna_genes(args.trna_gtf, host_contigs)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    rows = list(trna_genes.items())
    rows.sort(key=lambda item: (str(item[1]["chrom"]), int(item[1]["center"]), item[0]))

    with open(args.output, "w", encoding="utf-8") as out:
        for gene_id, meta in rows:
            chrom = str(meta["chrom"])
            center = int(meta["center"])
            start = max(0, center - args.flank_size)
            end = min(chromsizes[chrom], center + args.flank_size + 1)
            bin_id = f'{gene_id}|{meta["gene_type"]}'
            out.write(
                "\t".join(
                    [
                        chrom,
                        str(start),
                        str(end),
                        bin_id,
                        gene_id,
                        str(meta["gene_name"]),
                        str(meta["gene_type"]),
                        str(meta["strand"]),
                        str(center),
                    ]
                )
                + "\n"
            )


if __name__ == "__main__":
    main()
