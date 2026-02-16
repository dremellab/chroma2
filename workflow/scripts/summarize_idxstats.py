#!/usr/bin/env python3
"""Summarize idxstats files into host/virus/chrM totals per step.

Reads idxstats for steps (aligned, clean, dedup, final), and uses regions files
to sum host and per-virus alignments. chrM is reported explicitly.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple
import zipfile


DEFAULT_STEPS = ["aligned", "clean", "dedup", "final"]


def detect_step(filename: str, steps: List[str]) -> str | None:
    if filename.endswith(".aligned.idxstats.txt"):
        return "aligned"
    for step in steps:
        if step == "aligned":
            continue
        if f".aligned.{step}.idxstats.txt" in filename:
            return step
    return None


def detect_sample(filename: str) -> str | None:
    if ".aligned." in filename:
        return filename.split(".aligned.")[0]
    if filename.endswith(".aligned.idxstats.txt"):
        return filename[: -len(".aligned.idxstats.txt")]
    return None


def parse_idxstats(path: Path) -> Tuple[Dict[str, int], int, int]:
    mapped_by_contig: Dict[str, int] = {}
    total_mapped = 0
    unmapped = 0
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            rname = parts[0]
            mapped = int(parts[2])
            unmapped_count = int(parts[3])
            if rname == "*":
                unmapped += unmapped_count
                continue
            mapped_by_contig[rname] = mapped_by_contig.get(rname, 0) + mapped
            total_mapped += mapped
    return mapped_by_contig, total_mapped, unmapped


def read_regions_file(path: Path) -> Dict[str, Set[str]]:
    genomes: Dict[str, Set[str]] = {}
    if not path.exists():
        return genomes
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            genome = parts[0]
            contigs = [c for c in parts[1:] if c]
            genomes.setdefault(genome, set()).update(contigs)
    return genomes


def sum_contigs(mapped_by_contig: Dict[str, int], contigs: Iterable[str]) -> int:
    return sum(mapped_by_contig.get(c, 0) for c in contigs)


def default_ref_dir(results_dir: Path) -> Path | None:
    if results_dir.name == "results" and (results_dir.parent / "ref").exists():
        return results_dir.parent / "ref"
    if Path("ref").exists():
        return Path("ref")
    return None


def find_fastqc_zip(
    results_dir: Path, sample: str, kind: str, read: str
) -> Path | None:
    pattern = f"{sample}.{kind}_{read}_fastqc.zip"
    for path in results_dir.rglob(pattern):
        return path
    return None


def extract_total_sequences(fastqc_zip: Path) -> int | None:
    if fastqc_zip is None or not fastqc_zip.exists():
        return None
    try:
        with zipfile.ZipFile(fastqc_zip) as zf:
            for name in zf.namelist():
                if name.endswith("fastqc_data.txt"):
                    with zf.open(name) as handle:
                        for raw in handle:
                            try:
                                line = raw.decode("utf-8").strip()
                            except UnicodeDecodeError:
                                continue
                            if line.startswith("Total Sequences"):
                                parts = line.split()
                                if parts:
                                    return int(parts[-1])
                    break
    except (zipfile.BadZipFile, OSError, ValueError):
        return None
    return None


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Summarize idxstats files into host/virus/chrM totals per step."
    )
    parser.add_argument(
        "--results-dir",
        default="results",
        help="Results directory to scan (default: results)",
    )
    parser.add_argument(
        "--ref-dir",
        default="",
        help="Ref directory containing ref.fa.regions.* files (default: infer from results dir)",
    )
    parser.add_argument(
        "--regions-host",
        default="",
        help="Host regions file (default: <ref-dir>/ref.fa.regions.host)",
    )
    parser.add_argument(
        "--regions-viruses",
        default="",
        help="Viruses regions file (default: <ref-dir>/ref.fa.regions.viruses)",
    )
    parser.add_argument(
        "--chrM-names",
        default="chrM,MT",
        help="Comma-separated chrM contig names (default: chrM,MT)",
    )
    parser.add_argument(
        "--output",
        default="-",
        help="Output TSV file (default: stdout)",
    )
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    if not results_dir.exists():
        print(f"ERROR: results dir not found: {results_dir}", file=sys.stderr)
        return 2

    steps = DEFAULT_STEPS[:]

    ref_dir = Path(args.ref_dir) if args.ref_dir else default_ref_dir(results_dir)
    regions_host = Path(args.regions_host) if args.regions_host else None
    regions_viruses = Path(args.regions_viruses) if args.regions_viruses else None

    if regions_host is None or regions_viruses is None:
        if ref_dir is None:
            print(
                "ERROR: could not infer ref dir; provide --ref-dir or --regions-host/--regions-viruses",
                file=sys.stderr,
            )
            return 2
        if regions_host is None:
            regions_host = ref_dir / "ref.fa.regions.host"
        if regions_viruses is None:
            regions_viruses = ref_dir / "ref.fa.regions.viruses"

    host_genomes = read_regions_file(regions_host)
    virus_genomes = read_regions_file(regions_viruses)
    host_contigs: Set[str] = set()
    for contigs in host_genomes.values():
        host_contigs.update(contigs)

    chrM_names = {c.strip() for c in args.chrM_names.split(",") if c.strip()}

    sample_step: Dict[str, Dict[str, Dict[str, int]]] = {}
    samples: Set[str] = set()

    for path in results_dir.glob("*/align/*.idxstats.txt"):
        step = detect_step(path.name, steps)
        if step is None:
            continue
        sample = detect_sample(path.name)
        if not sample:
            continue
        mapped_by_contig, total_mapped, unmapped = parse_idxstats(path)
        samples.add(sample)

        host_mapped = sum_contigs(mapped_by_contig, host_contigs)
        chrM_mapped = sum_contigs(mapped_by_contig, chrM_names)
        host_no_chrM = (
            host_mapped - chrM_mapped if host_mapped >= chrM_mapped else host_mapped
        )

        data: Dict[str, int] = {
            "total_mapped": total_mapped,
            "unmapped": unmapped,
            "host_mapped": host_mapped,
            "host_no_chrM_mapped": host_no_chrM,
            "chrM_mapped": chrM_mapped,
        }
        for virus, contigs in virus_genomes.items():
            data[f"virus_{virus}_mapped"] = sum_contigs(mapped_by_contig, contigs)

        sample_step.setdefault(sample, {})[step] = data

    virus_names = sorted(virus_genomes.keys())
    samples_sorted = sorted(samples)

    header_main = [
        "sample",
        "raw_input_sequences",
        "post_trimming_input_sequences",
        "total_unmapped",
        "host_raw",
        "host_no_multimappers",
        "host_dedup",
        "host_mapqfiltered",
        "chrM_raw",
        "chrM_no_multimappers",
        "chrM_dedup",
        "chrM_mapqfiltered",
    ]
    for virus in virus_names:
        header_main.extend(
            [
                f"{virus}_raw",
                f"{virus}_no_multimappers",
                f"{virus}_dedup",
                f"{virus}_mapqfiltered",
            ]
        )

    if args.output == "-":
        out_main = sys.stdout
    else:
        out_main = open(args.output, "w")

    try:
        out_main.write("\t".join(header_main) + "\n")
        for sample in samples_sorted:
            aligned = sample_step.get(sample, {}).get("aligned", {})
            clean = sample_step.get(sample, {}).get("clean", {})
            dedup = sample_step.get(sample, {}).get("dedup", {})
            final = sample_step.get(sample, {}).get("final", {})

            raw_r1 = extract_total_sequences(
                find_fastqc_zip(results_dir, sample, "raw", "R1")
            )
            trimmed_r1 = extract_total_sequences(
                find_fastqc_zip(results_dir, sample, "trimmed", "R1")
            )
            has_r2 = (
                find_fastqc_zip(results_dir, sample, "raw", "R2") is not None
                or find_fastqc_zip(results_dir, sample, "trimmed", "R2") is not None
            )
            pe_factor = 2 if has_r2 else 1
            raw_input = raw_r1 * pe_factor if raw_r1 is not None else "NA"
            trimmed_input = trimmed_r1 * pe_factor if trimmed_r1 is not None else "NA"

            row_main = [
                sample,
                str(raw_input),
                str(trimmed_input),
                str(aligned.get("unmapped", "NA")),
                str(aligned.get("host_mapped", "NA")),
                str(clean.get("host_mapped", "NA")),
                str(dedup.get("host_mapped", "NA")),
                str(final.get("host_mapped", "NA")),
                str(aligned.get("chrM_mapped", "NA")),
                str(clean.get("chrM_mapped", "NA")),
                str(dedup.get("chrM_mapped", "NA")),
                str(final.get("chrM_mapped", "NA")),
            ]
            for virus in virus_names:
                row_main.extend(
                    [
                        str(aligned.get(f"virus_{virus}_mapped", "NA")),
                        str(clean.get(f"virus_{virus}_mapped", "NA")),
                        str(dedup.get(f"virus_{virus}_mapped", "NA")),
                        str(final.get(f"virus_{virus}_mapped", "NA")),
                    ]
                )
            out_main.write("\t".join(row_main) + "\n")
    finally:
        if out_main is not sys.stdout:
            out_main.close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
