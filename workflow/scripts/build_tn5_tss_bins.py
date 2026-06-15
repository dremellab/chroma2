#!/usr/bin/env python3
"""
Build fixed-width Tn5 counting bins around TSS (Transcription Start Site).

Bins are created centered on TSS with customizable flank size:
  bin_start = TSS - flank_size
  bin_end = TSS + flank_size

Suitable for protein-coding genes, rRNA genes, and other genes where TSS is meaningful.
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
            "Create fixed-width Tn5 counting bins centered on gene TSS positions. "
            "Suitable for protein-coding and rRNA genes."
        )
    )
    parser.add_argument("--gtf", required=True, help="Input GTF file")
    parser.add_argument("--host-regions", required=True, help="Host .fa.regions file")
    parser.add_argument("--chromsizes", required=True, help="Chromsizes TSV")
    parser.add_argument(
        "--gene-types",
        nargs="+",
        default=["protein_coding", "protein-coding"],
        help="Gene types to include (default: protein_coding, protein-coding)",
    )
    parser.add_argument(
        "--flank-size",
        type=int,
        default=250,
        help="Flank size around TSS in bp (default: 250). Bin will span TSS±flank_size.",
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


def tss_from_fields(fields: list[str]) -> int:
    """Extract TSS (0-based coordinate) from GTF fields.

    For + strand: TSS is the start position (0-based)
    For - strand: TSS is the end position (0-based, which is GTF end - 1)
    """
    start_1based = int(fields[3])
    end_1based = int(fields[4])
    strand = fields[6]
    return start_1based - 1 if strand == "+" else end_1based - 1


def iter_gtf_records(path: str) -> Iterator[list[str]]:
    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) >= 9:
                yield fields


def collect_gene_tss(
    gtf_path: str,
    host_contigs: set[str],
    include_types: Iterable[str],
) -> Dict[str, Dict[str, object]]:
    include = {item.lower() for item in include_types}
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
        gtype = gene_type(attrs).lower()
        if gtype not in include:
            continue

        strand = fields[6]
        tss = tss_from_fields(fields)
        name = feature_name(attrs, gid)

        existing = genes.get(gid)
        if existing is None:
            genes[gid] = {
                "chrom": chrom,
                "strand": strand,
                "tss": tss,
                "gene_name": name,
                "gene_type": gtype,
            }
            continue

        if existing["chrom"] != chrom or existing["strand"] != strand:
            continue
        if strand == "+":
            existing["tss"] = min(int(existing["tss"]), tss)
        else:
            existing["tss"] = max(int(existing["tss"]), tss)

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

    rows.sort(key=lambda item: (str(item[1]["chrom"]), int(item[1]["tss"]), item[0]))

    with open(path, "w", encoding="utf-8") as out:
        for gene_id, meta in rows:
            chrom = str(meta["chrom"])
            tss = int(meta["tss"])
            start = max(0, tss - flank_size)
            end = min(chromsizes[chrom], tss + flank_size + 1)
            bin_id = f"{gene_id}|tss"
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
                        str(tss),
                    ]
                )
                + "\n"
            )


def main() -> None:
    args = parse_args()
    if args.flank_size < 0 or args.flank_size > 1000:
        raise ValueError("--flank-size must be between 0 and 1000 bp")

    host_contigs = load_host_contigs(args.host_regions)
    chromsizes = load_chromsizes(args.chromsizes)

    genes = collect_gene_tss(
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
    print(f"Created {len(genes)} TSS bins in {args.output}", flush=True)


if __name__ == "__main__":
    main()
