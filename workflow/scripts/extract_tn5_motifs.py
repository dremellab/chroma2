#!/usr/bin/env python3
"""
Extract Tn5 nicking sites from a BAM, split them into host and per-virus groups,
write BED and strand-aware FASTA windows, build a position frequency matrix, and
render a sequence logo for a single caller-labeled output set. The same helper
can also emit virus-only fixed-width Tn5 bin counts from peak-caller input BAMs.

Outputs are written under:
  <outdir>/<group>/
"""

from __future__ import annotations

import argparse
import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, TextIO

import pysam

try:
    import pandas as pd
except ImportError:  # pragma: no cover
    pd = None

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:  # pragma: no cover
    plt = None

try:
    import logomaker
except ImportError:  # pragma: no cover
    logomaker = None


BASES = ("A", "C", "G", "T")
DEFAULT_FLANK_SIZE = 10


@dataclass(frozen=True)
class CutSite:
    chrom: str
    start: int
    end: int
    name: str
    score: int
    strand: str


class GroupWriter:
    def __init__(
        self,
        outdir: Path,
        sample: str,
        group: str,
        scenario_name: str,
        flank_size: int,
        dedup: bool,
        logo_format: str,
        skip_flank_output: bool = False,
    ) -> None:
        self.outdir = outdir / group
        self.outdir.mkdir(parents=True, exist_ok=True)

        window_size = 2 * flank_size + 1
        prefix = f"{sample}.{group}.{scenario_name}.tn5_sites.{window_size}bp"
        self.exact_bed_path = (
            self.outdir / f"{sample}.{group}.{scenario_name}.tn5_sites.1bp.bed.gz"
        )
        self.flank_bed_path = self.outdir / f"{prefix}.bed"
        self.fasta_path = self.outdir / f"{prefix}.fa"
        self.pfm_path = self.outdir / f"{prefix}.pfm.tsv.gz"
        self.logo_path = self.outdir / f"{prefix}.logo.{logo_format}"

        self.exact_bed_fh = gzip.open(self.exact_bed_path, "wt", encoding="utf-8")
        self.flank_bed_fh = (
            self.flank_bed_path.open("w", encoding="utf-8")
            if not skip_flank_output
            else None
        )
        self.fasta_fh = (
            self.fasta_path.open("w", encoding="utf-8")
            if not skip_flank_output
            else None
        )

        self.flank_size = flank_size
        self.exact_bed_count = 0
        self.flank_bed_count = 0
        self.fasta_count = 0
        self.pfm = [{base: 0 for base in BASES} for _ in range(window_size)]
        self.seen = set() if dedup else None

    def close(self) -> None:
        self.exact_bed_fh.close()
        if self.flank_bed_fh is not None:
            self.flank_bed_fh.close()
        if self.fasta_fh is not None:
            self.fasta_fh.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract Tn5 nicking sites from a BAM, write BED/FASTA/PFM outputs, "
            "render caller-specific sequence logos per host/virus group, or emit "
            "virus-only fixed-width Tn5 bin counts."
        )
    )
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument(
        "--fasta", required=True, help="Reference FASTA used for alignment"
    )
    parser.add_argument("--sample", required=True, help="Sample name for output naming")
    parser.add_argument(
        "--scenario-name", required=True, help="Caller label for output naming"
    )
    parser.add_argument(
        "--outdir", required=True, help="Caller-specific output directory"
    )
    parser.add_argument(
        "--host-regions",
        default=None,
        help="Host .fa.regions file used to define host contigs",
    )
    parser.add_argument(
        "--virus-regions",
        action="append",
        default=[],
        help="Virus mapping in the form NAME=/path/to/virus.fa.regions; repeatable",
    )
    parser.add_argument(
        "--flank-size",
        type=int,
        default=DEFAULT_FLANK_SIZE,
        help="Bases to include on each side of the nick site (default: 10)",
    )
    parser.add_argument(
        "--virus-bin-size",
        type=int,
        default=100,
        help="Fixed bin size used for viral Tn5 site count matrices (default: 100)",
    )
    parser.add_argument(
        "--bin-counts-only",
        action="store_true",
        help="Only write viral bin-count matrices; skip BED/FASTA/PFM/logo outputs",
    )
    parser.add_argument(
        "--bin-counts-output",
        default=None,
        help="Output TSV path used with --bin-counts-only",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Threads for BAM reading (default: 4)",
    )
    parser.add_argument(
        "--mapq-min",
        type=int,
        default=0,
        help="Minimum MAPQ to keep an alignment",
    )
    parser.add_argument(
        "--dedup",
        action="store_true",
        help="Drop duplicate cut sites within each group by chrom/start/end/strand",
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
        "--logo-format",
        choices=("png", "pdf", "svg"),
        default="png",
        help="Output format for sequence logo image (default: png)",
    )
    parser.add_argument(
        "--skip-logo",
        action="store_true",
        help="Skip sequence logo generation",
    )
    parser.add_argument(
        "--skip-flank-output",
        action="store_true",
        help="Skip writing flank BED/FASTA outputs (PFM still computed in-memory)",
    )
    return parser.parse_args()


def parse_regions_contigs(path: str) -> set[str]:
    contigs: set[str] = set()
    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                continue
            for contig in fields[1].split():
                if contig:
                    contigs.add(contig)
    return contigs


def load_group_contigs(args: argparse.Namespace) -> Dict[str, set[str]]:
    groups: Dict[str, set[str]] = {}
    if args.host_regions:
        groups["host"] = parse_regions_contigs(args.host_regions)
    for item in args.virus_regions:
        if "=" not in item:
            raise ValueError(
                f"--virus-regions must be NAME=/path/to/file.fa.regions, got: {item}"
            )
        name, path = item.split("=", 1)
        if not name:
            raise ValueError(f"Empty virus name in --virus-regions: {item}")
        groups[name] = parse_regions_contigs(path)
    return groups


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def alignment_passes_filters(
    rec: pysam.AlignedSegment,
    mapq_min: int,
    exclude_secondary: bool,
    exclude_supplementary: bool,
) -> bool:
    if rec.is_unmapped:
        return False
    if rec.mapping_quality < mapq_min:
        return False
    if rec.is_secondary and exclude_secondary:
        return False
    if rec.is_supplementary and exclude_supplementary:
        return False
    return True


def single_cut_site(rec: pysam.AlignedSegment, ordinal: int) -> List[CutSite]:
    if rec.reference_start is None or rec.reference_end is None:
        return []

    if rec.is_reverse:
        cut = rec.reference_end - 5
        strand = "-"
    else:
        cut = rec.reference_start + 4
        strand = "+"

    if cut < 0:
        return []

    return [
        CutSite(
            rec.reference_name,
            cut,
            cut + 1,
            f"{rec.query_name}|read{ordinal}",
            rec.mapping_quality,
            strand,
        )
    ]


def cut_sites_from_records(
    records: List[pysam.AlignedSegment],
    mapq_min: int,
    exclude_secondary: bool,
    exclude_supplementary: bool,
) -> List[CutSite]:
    usable = [
        rec
        for rec in records
        if alignment_passes_filters(
            rec,
            mapq_min=mapq_min,
            exclude_secondary=exclude_secondary,
            exclude_supplementary=exclude_supplementary,
        )
    ]
    if not usable:
        return []

    # Each mapped read independently reveals exactly one Tn5 cut site via its own
    # strand (single_cut_site handles both strands correctly); a proper pair simply
    # yields two sites, one from each mate, without needing to derive one mate's
    # site from the other's coordinates.
    sites: List[CutSite] = []
    for i, rec in enumerate(usable, start=1):
        sites.extend(single_cut_site(rec, i))
    return sites


def cut_sites_from_record(
    rec: pysam.AlignedSegment,
    mapq_min: int,
    exclude_secondary: bool,
    exclude_supplementary: bool,
) -> List[CutSite]:
    if not alignment_passes_filters(
        rec,
        mapq_min=mapq_min,
        exclude_secondary=exclude_secondary,
        exclude_supplementary=exclude_supplementary,
    ):
        return []
    return single_cut_site(rec, 1)


def assign_group(chrom: str, group_contigs: Dict[str, set[str]]) -> Optional[str]:
    hits = [group for group, contigs in group_contigs.items() if chrom in contigs]
    if not hits:
        return None
    if len(hits) > 1:
        raise ValueError(f"Contig {chrom} appears in multiple groups: {hits}")
    return hits[0]


def fetch_window_sequence(
    fasta: pysam.FastaFile,
    chrom: str,
    cut_start: int,
    strand: str,
    flank_size: int,
) -> Optional[str]:
    chrom_len = fasta.get_reference_length(chrom)
    left = cut_start - flank_size
    right = cut_start + flank_size + 1
    window_size = 2 * flank_size + 1

    if left < 0 or right > chrom_len:
        return None

    seq = fasta.fetch(chrom, left, right).upper()
    if len(seq) != window_size:
        return None
    if strand == "-":
        seq = reverse_complement(seq)
    return seq


def write_exact_bed(fh: TextIO, site: CutSite) -> None:
    fh.write(
        f"{site.chrom}\t{site.start}\t{site.end}\t"
        f"{site.name}\t{site.score}\t{site.strand}\n"
    )


def write_flank_bed(fh: TextIO, site: CutSite, flank_size: int) -> None:
    fh.write(
        f"{site.chrom}\t{site.start - flank_size}\t{site.end + flank_size}\t"
        f"{site.name}\t{site.score}\t{site.strand}\n"
    )


def update_writer(writer: GroupWriter, fasta: pysam.FastaFile, site: CutSite) -> None:
    dedup_key = (site.chrom, site.start, site.end, site.strand)
    if writer.seen is not None and dedup_key in writer.seen:
        return
    if writer.seen is not None:
        writer.seen.add(dedup_key)

    write_exact_bed(writer.exact_bed_fh, site)
    writer.exact_bed_count += 1
    if writer.flank_bed_fh is not None:
        write_flank_bed(writer.flank_bed_fh, site, writer.flank_size)
        writer.flank_bed_count += 1

    seq = fetch_window_sequence(
        fasta=fasta,
        chrom=site.chrom,
        cut_start=site.start,
        strand=site.strand,
        flank_size=writer.flank_size,
    )
    if seq is None or any(base not in BASES for base in seq):
        return

    if writer.fasta_fh is not None:
        header = f"{site.name}|{site.chrom}:{site.start + 1}-{site.end}|{site.strand}"
        writer.fasta_fh.write(f">{header}\n{seq}\n")
        writer.fasta_count += 1
    for i, base in enumerate(seq):
        writer.pfm[i][base] += 1


def write_pfm(path: Path, pfm: List[Dict[str, int]]) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as out:
        out.write("position\tA\tC\tG\tT\n")
        for i, counts in enumerate(pfm, start=1):
            out.write(
                f"{i}\t{counts['A']}\t{counts['C']}\t{counts['G']}\t{counts['T']}\n"
            )


def write_bin_counts(
    path: Path,
    sample: str,
    bin_counts: Dict[str, List[int]],
    chrom_lengths: Dict[str, int],
    bin_size: int,
) -> None:
    with path.open("w", encoding="utf-8") as out:
        out.write("sample\tchrom\tstart\tend\tbin_id\tbin_index\ttn5_site_count\n")
        for chrom, counts in bin_counts.items():
            chrom_len = chrom_lengths[chrom]
            for bin_idx, count in enumerate(counts):
                bin_start = bin_idx * bin_size
                bin_end = min(chrom_len, bin_start + bin_size)
                bin_id = f"{chrom}:{bin_start}-{bin_end}"
                out.write(
                    f"{sample}\t{chrom}\t{bin_start}\t{bin_end}\t{bin_id}\t{bin_idx}\t{count}\n"
                )


def pfm_to_probability_df(pfm: List[Dict[str, int]]) -> pd.DataFrame:
    rows = []
    for counts in pfm:
        total = sum(counts.values())
        if total == 0:
            rows.append({base: 0.0 for base in BASES})
        else:
            rows.append({base: counts[base] / total for base in BASES})
    df = pd.DataFrame(rows, columns=list(BASES))
    df.index = range(1, len(df) + 1)
    return df


def write_logo(path: Path, pfm: List[Dict[str, int]], title: str) -> bool:
    if logomaker is None or plt is None or pd is None:
        return False
    df = pfm_to_probability_df(pfm)
    fig, ax = plt.subplots(figsize=(max(6, len(df) * 0.35), 3.5))
    logomaker.Logo(df, ax=ax)
    ax.set_title(title)
    ax.set_xlabel("Position")
    ax.set_ylabel("Frequency")
    ax.set_xlim(-0.5, len(df) - 0.5)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return True


def main() -> None:
    args = parse_args()
    if args.flank_size < 0:
        raise ValueError("--flank-size must be >= 0")
    if args.virus_bin_size <= 0:
        raise ValueError("--virus-bin-size must be > 0")
    if args.mapq_min < 0:
        raise ValueError("--mapq-min must be >= 0")
    if args.bin_counts_only and not args.bin_counts_output:
        raise ValueError("--bin-counts-output is required with --bin-counts-only")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    fasta_fai = Path(f"{args.fasta}.fai")
    if not fasta_fai.exists():
        pysam.faidx(args.fasta)

    group_contigs = load_group_contigs(args)
    virus_groups = {item.split("=", 1)[0] for item in args.virus_regions}
    if args.bin_counts_only:
        if len(virus_groups) != 1:
            raise ValueError(
                "--bin-counts-only expects exactly one virus group/input per invocation"
            )
        virus_group = next(iter(virus_groups))
        bin_counts_output = Path(args.bin_counts_output)
        chrom_lengths: Dict[str, int] = {}
        bin_counts: Dict[str, List[int]] = {}

        with pysam.FastaFile(args.fasta) as fasta:
            chrom_lengths = {
                chrom: fasta.get_reference_length(chrom)
                for chrom in sorted(group_contigs[virus_group])
            }
            bin_counts = {
                chrom: [0]
                * ((chrom_len + args.virus_bin_size - 1) // args.virus_bin_size)
                for chrom, chrom_len in chrom_lengths.items()
            }

            with pysam.AlignmentFile(args.bam, "rb", threads=args.threads) as bam:
                for rec in bam.fetch(until_eof=True):
                    for site in cut_sites_from_record(
                        rec=rec,
                        mapq_min=args.mapq_min,
                        exclude_secondary=args.exclude_secondary,
                        exclude_supplementary=args.exclude_supplementary,
                    ):
                        group = assign_group(site.chrom, group_contigs)
                        if group != virus_group:
                            continue
                        bin_idx = site.start // args.virus_bin_size
                        bin_counts[site.chrom][bin_idx] += 1

        bin_counts_output.parent.mkdir(parents=True, exist_ok=True)
        write_bin_counts(
            bin_counts_output,
            args.sample,
            bin_counts,
            chrom_lengths,
            args.virus_bin_size,
        )
        total_sites = sum(sum(counts) for counts in bin_counts.values())
        print(
            f"[{args.scenario_name}][{virus_group}] bin_counts_only total_sites={total_sites} "
            f"bin_counts={bin_counts_output}"
        )
        return

    writers: Dict[str, GroupWriter] = {}

    with pysam.FastaFile(args.fasta) as fasta, pysam.AlignmentFile(
        args.bam, "rb", threads=args.threads
    ) as bam:
        current_qname: Optional[str] = None
        bucket: List[pysam.AlignedSegment] = []

        def flush(records: List[pysam.AlignedSegment]) -> None:
            if not records:
                return
            for site in cut_sites_from_records(
                records=records,
                mapq_min=args.mapq_min,
                exclude_secondary=args.exclude_secondary,
                exclude_supplementary=args.exclude_supplementary,
            ):
                group = assign_group(site.chrom, group_contigs)
                if group is None:
                    continue
                if group not in writers:
                    writers[group] = GroupWriter(
                        outdir=outdir,
                        sample=args.sample,
                        group=group,
                        scenario_name=args.scenario_name,
                        flank_size=args.flank_size,
                        dedup=args.dedup,
                        logo_format=args.logo_format,
                        skip_flank_output=args.skip_flank_output,
                    )
                update_writer(writers[group], fasta, site)

        for rec in bam.fetch(until_eof=True):
            qname = rec.query_name
            if current_qname is None:
                current_qname = qname
            elif qname != current_qname:
                flush(bucket)
                bucket = []
                current_qname = qname
            bucket.append(rec)
        flush(bucket)

    for group, writer in sorted(writers.items()):
        writer.close()
        write_pfm(writer.pfm_path, writer.pfm)
        if not args.skip_logo:
            logo_written = write_logo(
                writer.logo_path,
                writer.pfm,
                title=f"{args.sample} {group} Tn5 motif ({args.scenario_name})",
            )
        else:
            logo_written = False
        print(
            f"[{args.scenario_name}][{group}] exact_bed_sites={writer.exact_bed_count} "
            f"flank_bed_sites={writer.flank_bed_count} fasta_seqs={writer.fasta_count} "
            f"exact_bed={writer.exact_bed_path} flank_bed={writer.flank_bed_path} "
            f"fasta={writer.fasta_path} pfm={writer.pfm_path} "
            f"logo={writer.logo_path if logo_written else 'SKIPPED (install logomaker/matplotlib)'}"
        )


if __name__ == "__main__":
    main()
