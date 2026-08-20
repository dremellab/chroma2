#!/usr/bin/env python3
"""
GTF-parsing helpers shared by build_tn5_midpoint_bins.py, build_tn5_tss_bins.py,
and build_viral_genome_bins.py.
"""

from __future__ import annotations

import gzip
from typing import Dict, Iterator, Optional


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


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


def iter_gtf_records(path: str) -> Iterator[list[str]]:
    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) >= 9:
                yield fields
