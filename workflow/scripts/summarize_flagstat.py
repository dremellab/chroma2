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
            metrics[label] = count
    return metrics


def detect_step(name: str) -> str | None:
    if name.endswith(".aligned.flagstat.txt"):
        return "aligned"
    for step in STEP_ORDER[1:]:
        if f".aligned.{step}.flagstat.txt" in name:
            return step
    return None


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

    samples = sorted(sample_step_value.keys())
    header = ["sample"]
    for step in STEP_ORDER:
        for metric in requested_metrics:
            header.append(f"{step}_{metric}")

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
            out.write("\t".join(row) + "\n")
    finally:
        if out is not sys.stdout:
            out.close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
