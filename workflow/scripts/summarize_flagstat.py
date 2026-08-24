#!/usr/bin/env python3
"""Summarize flagstat files into a wide table.

Default metric is total reads ("in total").
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
import re

STEP_ORDER = ["aligned", "clean", "fixmate", "dedup", "final"]

METRIC_KEYS = {
    "total": "in total",
    "mapped": "mapped",
    "properly_paired": "properly paired",
    "duplicates": "duplicates",
}

LINE_RE = re.compile(r"^(\d+)\s+\+\s+\d+\s+(.+)$")


def parse_flagstat(path: Path) -> dict[str, int]:
    metrics: dict[str, int] = {}
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            m = LINE_RE.match(line)
            if not m:
                continue
            count = int(m.group(1))
            label = m.group(2)
            # Normalize labels that include extra context, e.g.
            # "in total (QC-passed reads + QC-failed reads)" -> "in total"
            label = label.split(" (", 1)[0]
            metrics[label] = count
    return metrics


def detect_step(name: str) -> str | None:
    if name.endswith(".aligned.flagstat.txt"):
        return "aligned"
    for step in STEP_ORDER[1:]:
        if f".aligned.{step}.flagstat.txt" in name:
            return step
    return None


def detect_step_idxstats(name: str) -> str | None:
    if name.endswith(".aligned.idxstats.txt"):
        return "aligned"
    for step in STEP_ORDER[1:]:
        if f".aligned.{step}.idxstats.txt" in name:
            return step
    return None


def parse_idxstats(path: Path, chrnames: set[str]) -> dict[str, int]:
    total_mapped = 0
    chr_mapped = 0
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
            if rname != "*":
                total_mapped += mapped
            if rname in chrnames:
                chr_mapped += mapped
    return {"total_mapped": total_mapped, "chr_mapped": chr_mapped}


def detect_sample(name: str) -> str | None:
    if ".aligned." in name:
        return name.split(".aligned.")[0]
    if name.endswith(".aligned.flagstat.txt"):
        return name[: -len(".aligned.flagstat.txt")]
    return None


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Summarize flagstat files into a table (rows=samples, columns=steps)."
    )
    parser.add_argument(
        "--results-dir",
        default="results",
        help="Results directory to scan (default: results)",
    )
    parser.add_argument(
        "--metrics",
        default="total,mapped,properly_paired,duplicates",
        help=(
            "Comma-separated metrics to report per step "
            f"(default: { 'total,mapped,properly_paired,duplicates' })"
        ),
    )
    parser.add_argument(
        "--chrnames",
        default="chrM,MT",
        help="Comma-separated contig names to count (default: chrM,MT).",
    )
    parser.add_argument(
        "--output",
        default="-",
        help="Output TSV file (default: stdout)",
    )
    parser.add_argument(
        "--fill",
        default="NA",
        help="Fill value for missing step (default: NA)",
    )
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    if not results_dir.exists():
        print(f"ERROR: results dir not found: {results_dir}", file=sys.stderr)
        return 2

    requested_metrics = [m.strip() for m in args.metrics.split(",") if m.strip()]
    for metric in requested_metrics:
        if metric not in METRIC_KEYS:
            print(f"ERROR: unknown metric: {metric}", file=sys.stderr)
            return 2

    sample_step_value: dict[str, dict[str, dict[str, int]]] = {}
    sample_step_chr: dict[str, dict[str, dict[str, int]]] = {}

    for path in results_dir.glob("*/align/*.flagstat.txt"):
        step = detect_step(path.name)
        if step is None:
            continue
        sample = detect_sample(path.name)
        if not sample:
            continue
        metrics = parse_flagstat(path)
        for metric in requested_metrics:
            key = METRIC_KEYS[metric]
            if key not in metrics:
                continue
            sample_step_value.setdefault(sample, {}).setdefault(step, {})[
                metric
            ] = metrics[key]

    chrnames = {c.strip() for c in args.chrnames.split(",") if c.strip()}
    for path in results_dir.glob("*/align/*.idxstats.txt"):
        step = detect_step_idxstats(path.name)
        if step is None:
            continue
        sample = detect_sample(path.name)
        if not sample:
            continue
        chr_metrics = parse_idxstats(path, chrnames)
        sample_step_chr.setdefault(sample, {}).setdefault(step, {}).update(chr_metrics)

    samples = sorted(sample_step_value.keys())
    header = ["sample"]
    for step in STEP_ORDER:
        for metric in requested_metrics:
            header.append(f"{step}_{metric}")
        header.append(f"{step}_mt_mapped")
        header.append(f"{step}_mt_pct_mapped")

    if args.output == "-":
        out = sys.stdout
    else:
        out = open(args.output, "w")

    try:
        out.write("\t".join(header) + "\n")
        for sample in samples:
            row = [sample]
            for step in STEP_ORDER:
                for metric in requested_metrics:
                    value = (
                        sample_step_value.get(sample, {})
                        .get(step, {})
                        .get(metric, args.fill)
                    )
                    row.append(str(value))
                chr_mapped = (
                    sample_step_chr.get(sample, {})
                    .get(step, {})
                    .get("chr_mapped", args.fill)
                )
                total_mapped = (
                    sample_step_chr.get(sample, {})
                    .get(step, {})
                    .get("total_mapped", None)
                )
                if chr_mapped == args.fill:
                    row.append(str(args.fill))
                    row.append(str(args.fill))
                else:
                    row.append(str(chr_mapped))
                    if not total_mapped:
                        row.append("0")
                    else:
                        pct = (float(chr_mapped) / float(total_mapped)) * 100.0
                        row.append(f"{pct:.4f}")
            out.write("\t".join(row) + "\n")
    finally:
        if out is not sys.stdout:
            out.close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
