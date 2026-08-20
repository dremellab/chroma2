#!/usr/bin/env python3
"""
Build fixed-width Tn5 counting bins around gene midpoint with dual-mode GTF processing.

Bins are created centered on the midpoint of the gene with customizable flank size:
  bin_start = midpoint - flank_size
  bin_end = midpoint + flank_size

DUAL-MODE OPERATION:
  Mode A (--gene-types specified): Filter GTF by feature type (standard operation)
    - Groups GTF lines by gene_id, merges overlapping features
    - Output: 9-column BED
    - Use case: Pol3 GTFs, standard gene GTFs with type filtering

  Mode B (--gene-types NOT specified): Per-line feature processing (repeatMasker mode)
    - Treats each GTF line as individual feature, no merging
    - Accepts all feature types (SINE/Alu, LINE/L1, etc.)
    - Output: 13-column BED (gene_id/gene_name + repeat metadata: feature_type,
      repeat_name, repeat_family, sw_score)
    - Use case: repeatMasker GTF files without standard gene_id attributes

VERSION: 2.2 (2026-08-20)
  - Mode B now also writes gene_id/gene_name (previously grouped-mode only),
    giving both schemas a common gene_id/gene_name/gene_type/strand prefix
  - Both modes' column schemas live in gtf_common.py (GROUPED_BIN_METADATA_COLUMNS/
    PERLINE_BIN_METADATA_COLUMNS), shared with count_tn5_sites_in_bins.py so the
    two sides can't silently drift apart on column count/order again
  - Metadata flows through count_tn5_sites_in_bins.py to final count matrices
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, Iterable

from gtf_common import (
    GROUPED_BIN_METADATA_COLUMNS,
    PERLINE_BIN_METADATA_COLUMNS,
    feature_id,
    feature_name,
    gene_type,
    iter_gtf_records,
    load_chromsizes,
    load_host_contigs,
    parse_attributes,
)


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
    parser.add_argument(
        "--max-size",
        type=int,
        default=1000,
        help="Maximum feature size in bp (default: 1000). Only features smaller than this are included. "
        "Use 0 or negative to disable size filtering.",
    )
    parser.add_argument("--output", required=True, help="Output BED path")
    return parser.parse_args()


def gene_center_from_fields(fields: list[str]) -> int:
    """Extract gene center (0-based coordinate) from GTF fields.

    Returns the midpoint between gene start and end.
    """
    start_0based = int(fields[3]) - 1
    end_0based = int(fields[4]) - 1
    return (start_0based + end_0based) // 2


def collect_gene_centers(
    gtf_path: str,
    host_contigs: set[str],
    include_types: Iterable[str] | None,
) -> Dict[str, Dict[str, object]]:
    include = {item.lower() for item in include_types} if include_types else None
    genes: Dict[str, Dict[str, object]] = {}

    for fields in iter_gtf_records(gtf_path):
        feature = fields[2]
        chrom = fields[0]
        if chrom not in host_contigs:
            continue
        attrs = parse_attributes(fields[8])

        # Dual-mode processing:
        # - If include_types specified: filter by feature type (for targeted filtering)
        # - If NOT specified: accept all features, each line is unique (for repeatMasker GTFs)
        if include is not None:
            if feature.lower() not in include:
                continue

        gid = feature_id(attrs)
        gtype = gene_type(attrs)

        # For repeatMasker GTFs or other files without standard gene_id, use feature-based ID
        if not gid:
            # Create ID from feature type, coordinates, and repeat_name if available
            repeat_name = attrs.get("repeat_name", "")
            gid = f"{feature}_{fields[3]}_{fields[4]}_{repeat_name}".replace(
                " ", "_"
            ).replace("/", "_")

        strand = fields[6]
        start_0based = int(fields[3]) - 1
        end_0based = int(fields[4]) - 1
        center = (start_0based + end_0based) // 2
        name = feature_name(attrs, gid if gid else repeat_name or feature)

        # Build base metadata
        metadata = {
            "chrom": chrom,
            "strand": strand,
            "center": center,
            "start": start_0based,
            "end": end_0based,
            "gene_name": name,
            "gene_type": gtype if gtype else feature,
        }

        # If processing without --gene-types (no filtering), add repeat metadata
        if include is None:
            metadata["feature_type"] = feature
            metadata["repeat_name"] = attrs.get("repeat_name", "")
            metadata["repeat_family"] = attrs.get("repeat_family", "")
            metadata["sw_score"] = attrs.get("sw_score", "")

        existing = genes.get(gid)
        if existing is None:
            genes[gid] = metadata
            continue

        # Only merge if both have same chrom/strand
        if existing["chrom"] != chrom or existing["strand"] != strand:
            continue

        # Merge: keep outermost span, preserve metadata
        if include is not None:
            # In filtered mode, merge is expected
            existing["start"] = min(int(existing["start"]), start_0based)
            existing["end"] = max(int(existing["end"]), end_0based)
        # In per-line mode (include is None), don't merge - each line is unique due to ID

    # Recompute center from the final (possibly multi-line-merged) span, rather than
    # accumulating an order-dependent running average of each line's own center.
    for meta in genes.values():
        meta["center"] = (int(meta["start"]) + int(meta["end"])) // 2

    return genes


def write_bins(
    path: str,
    chromsizes: Dict[str, int],
    flank_size: int,
    max_size: int,
    genes: Dict[str, Dict[str, object]],
) -> None:
    rows = []
    for gene_id, meta in genes.items():
        rows.append((gene_id, meta))

    rows.sort(key=lambda item: (str(item[1]["chrom"]), int(item[1]["center"]), item[0]))

    # Check for chromosomes in GTF that aren't in chromsizes
    chroms_in_genes = set(str(meta["chrom"]) for _, meta in rows)
    missing_chroms = chroms_in_genes - set(chromsizes.keys())
    if missing_chroms:
        print(
            f"⚠️  WARNING: {len(missing_chroms)} chromosome(s) in GTF not found in chromsizes: {sorted(missing_chroms)}",
            file=sys.stderr,
        )

    # Track filtering and truncation statistics
    filtered_by_size = 0
    truncated_count = 0
    truncated_start = 0
    truncated_end = 0

    with open(path, "w", encoding="utf-8") as out:
        for gene_id, meta in rows:
            chrom = str(meta["chrom"])
            if chrom not in chromsizes:
                continue

            # Apply size filter. meta["start"]/meta["end"] are 0-based
            # *inclusive* coordinates, so true length is end - start + 1.
            if max_size > 0:
                feature_size = int(meta["end"]) - int(meta["start"]) + 1
                if feature_size > max_size:
                    filtered_by_size += 1
                    continue

            if flank_size == -1:
                # Use full feature span
                start = int(meta["start"])
                end = int(meta["end"]) + 1
                reference_pos = int(meta["center"])
            else:
                # Use midpoint ± flank_size
                center = int(meta["center"])
                desired_start = center - flank_size
                desired_end = center + flank_size + 1
                start = max(0, desired_start)
                end = min(chromsizes[chrom], desired_end)
                reference_pos = center

                if start != desired_start or end != desired_end:
                    truncated_count += 1
                    if start != desired_start:
                        truncated_start += 1
                    if end != desired_end:
                        truncated_end += 1

            # Generate bin_id: format differs for grouped vs per-line mode
            if "feature_type" in meta:
                # Per-line mode: feature_type_chrom_start_end_name
                feature_clean = meta["gene_type"].replace("/", "_")
                start_1based = int(meta["start"]) + 1
                end_1based = int(meta["end"]) + 1
                repeat_name = meta.get("repeat_name", "")
                bin_id = f"{feature_clean}_{chrom}_{start_1based}_{end_1based}"
                if repeat_name:
                    bin_id += f"_{repeat_name}"
            else:
                # Grouped mode: include gene_id and gene_name to distinguish overlapping liftover coordinates
                gene_name = meta.get("gene_name", gene_id)
                bin_id = f"{gene_id}_{gene_name}"

            # gene_id/gene_name/gene_type/strand are common to both schemas;
            # per-line mode adds reference_pos/repeat metadata instead of tss.
            # Building this as a name->value lookup (rather than positional
            # list-building) means a future schema change that isn't matched
            # here raises KeyError immediately instead of silently shifting
            # every column after it.
            metadata_values = {
                "gene_id": str(gene_id),
                "gene_name": str(meta.get("gene_name", gene_id)),
                "gene_type": str(meta["gene_type"]),
                "strand": str(meta["strand"]),
            }
            if "feature_type" in meta:
                schema = PERLINE_BIN_METADATA_COLUMNS
                metadata_values["reference_pos"] = str(reference_pos)
                metadata_values["feature_type"] = str(meta["feature_type"])
                metadata_values["repeat_name"] = str(meta["repeat_name"])
                metadata_values["repeat_family"] = str(meta["repeat_family"])
                metadata_values["sw_score"] = str(meta["sw_score"])
            else:
                schema = GROUPED_BIN_METADATA_COLUMNS
                metadata_values["tss"] = str(reference_pos)

            output_fields = [chrom, str(start), str(end), bin_id] + [
                metadata_values[column] for column in schema
            ]

            out.write("\t".join(output_fields) + "\n")

    # Report size filtering statistics
    if max_size > 0 and filtered_by_size > 0:
        print(
            f"ℹ️  SIZE FILTER REPORT: {filtered_by_size} features excluded (larger than {max_size} bp)",
            file=sys.stderr,
        )

    # Report truncation statistics (only for midpoint mode, not full feature span)
    if flank_size != -1 and truncated_count > 0:
        print(
            f"⚠️  TRUNCATION REPORT: {truncated_count} bins were truncated to fit chromosome boundaries",
            file=sys.stderr,
        )
        print(
            f"   - {truncated_start} truncated at start (midpoint near chromosome start)",
            file=sys.stderr,
        )
        print(
            f"   - {truncated_end} truncated at end (midpoint near chromosome end)",
            file=sys.stderr,
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
        args.max_size,
        genes,
    )

    # Count bins actually written (excluding those on missing chromosomes)
    chroms_in_genes = set(str(meta["chrom"]) for meta in genes.values())
    missing_chroms = chroms_in_genes - set(chromsizes.keys())
    skipped_genes = sum(
        1 for meta in genes.values() if str(meta["chrom"]) in missing_chroms
    )
    written_bins = len(genes) - skipped_genes

    mode = (
        "full feature span"
        if args.flank_size == -1
        else f"midpoint±{args.flank_size}bp"
    )
    print(
        f"✓ Created {written_bins} bins ({mode}) in {args.output}",
        flush=True,
    )
    if skipped_genes > 0:
        print(
            f"  ({skipped_genes} genes skipped due to missing chromosome in chromsizes)",
            flush=True,
        )


if __name__ == "__main__":
    main()
