#!/usr/bin/env python3
"""
Build fixed-width Tn5 counting bins around gene midpoint.

Bins are created centered on the midpoint of the gene with customizable flank size:
  bin_start = midpoint - flank_size
  bin_end = midpoint + flank_size

Suitable for small genes like tRNAs, snRNAs, and other compact features where
the entire gene is relevant for analysis.
One bin is emitted per gene.
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
            "Create fixed-width Tn5 counting bins centered on gene midpoints. "
            "Suitable for small genes like tRNAs where the entire gene is relevant."
        )
    )
    parser.add_argument("--gtf", required=True, help="Input GTF file")
    parser.add_argument("--host-regions", required=True, help="Host .fa.regions file")
    parser.add_argument("--chromsizes", required=True, help="Chromsizes TSV")
    parser.add_argument(
        "--gene-types",
        nargs="+",
        default=None,
        help="Gene types to include (default: None = all gene types). If specified, only these types are included.",
    )
    parser.add_argument(
        "--flank-size",
        type=int,
        default=100,
        help="Flank size around midpoint in bp (default: 100). Bin will span midpoint±flank_size. "
        "Use -1 to use the actual feature start and end coordinates (full gene span).",
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
    for key in (
        "gene_type",
        "gene_biotype",
        "transcript_type",
        "transcript_biotype",
        "gene_biotype_ENSEMBL",
    ):
        value = attrs.get(key)
        if value:
            return value
    return ""


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
    """Extract gene center (0-based coordinate) from GTF fields.

    Returns the midpoint between gene start and end.
    """
    start_0based = int(fields[3]) - 1
    end_0based = int(fields[4]) - 1
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


def collect_gene_centers(
    gtf_path: str,
    host_contigs: set[str],
    include_types: Iterable[str] | None,
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
        start_0based = int(fields[3]) - 1
        end_0based = int(fields[4]) - 1
        center = (start_0based + end_0based) // 2
        name = feature_name(attrs, gid)

        existing = genes.get(gid)
        if existing is None:
            genes[gid] = {
                "chrom": chrom,
                "strand": strand,
                "center": center,
                "start": start_0based,
                "end": end_0based,
                "gene_name": name,
                "gene_type": gtype,
            }
            continue

        if existing["chrom"] != chrom or existing["strand"] != strand:
            continue
        # keep the outermost span → recompute center and expand start/end
        existing_center = int(existing["center"])
        existing["center"] = (existing_center + center) // 2
        existing["start"] = min(int(existing["start"]), start_0based)
        existing["end"] = max(int(existing["end"]), end_0based)

    return genes


def write_bins(
    path: str,
    chromsizes: Dict[str, int],
    flank_size: int,
    genes: Dict[str, Dict[str, object]],
) -> None:
    rows = []
    for gene_id, meta in genes.items():
        rows.append((gene_id, meta))

    rows.sort(key=lambda item: (str(item[1]["chrom"]), int(item[1]["center"]), item[0]))

    with open(path, "w", encoding="utf-8") as out:
        for gene_id, meta in rows:
            chrom = str(meta["chrom"])
            if flank_size == -1:
                # Use full feature span
                start = int(meta["start"])
                end = int(meta["end"]) + 1
                reference_pos = int(meta["center"])
            else:
                # Use midpoint ± flank_size
                center = int(meta["center"])
                start = max(0, center - flank_size)
                end = min(chromsizes[chrom], center + flank_size + 1)
                reference_pos = center
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
                        str(reference_pos),
                    ]
                )
                + "\n"
            )


def main() -> None:
    args = parse_args()
    if args.flank_size < -1 or (args.flank_size > 1000 and args.flank_size != -1):
        raise ValueError(
            "--flank-size must be between 0 and 1000 bp, or -1 to use full feature span"
        )

    host_contigs = load_host_contigs(args.host_regions)
    chromsizes = load_chromsizes(args.chromsizes)

    genes = collect_gene_centers(
        args.gtf,
        host_contigs,
        include_types=args.gene_types,
    )

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    write_bins(
        args.output,
        chromsizes,
        args.flank_size,
        genes,
    )
    print(f"Created {len(genes)} midpoint bins in {args.output}", flush=True)


if __name__ == "__main__":
    main()
