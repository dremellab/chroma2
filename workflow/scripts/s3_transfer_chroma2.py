#!/usr/bin/env python3
"""
s3_transfer_chroma2.py

Transfer Chroma 2 pipeline outputs to S3. Called directly from the pipeline
with all config provided via command-line arguments. No external lookup tables.

Usage:
  s3_transfer_chroma2.py \\
    --workdir /path/to/pipeline/workdir \\
    --pipeline-name CHROMA \\
    --sample-set-name myrun \\
    --bucket dremel-lab-bucket \\
    --s3-prefix _HTS \\
    --storage-class GLACIER_IR \\
    --large-file-storage-class GLACIER
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path
from typing import List, Tuple, Optional


# Transfer rules for Chroma 2 outputs
CHROMA2_RULES = [
    {"kind": "path", "path": "samples.tsv", "dest": "config/samples.tsv"},
    {"kind": "path", "path": "config.yaml", "dest": "config/config.yaml"},
    {"kind": "path", "path": "contrasts.tsv", "dest": "config/contrasts.tsv"},
    {
        "kind": "path",
        "path": "results/alignmentqc/idxstats_summary.tsv",
        "dest": "qc/alignment_summary.tsv",
    },
    {
        "kind": "path",
        "path": "results/multiqc_report.html",
        "dest": "qc/multiqc_report.html",
    },
    {"kind": "dir", "dir_name": "results/multiqc_data", "dest": "qc/multiqc_data"},
    {"kind": "dir", "dir_name": "results/alignmentqc/ataqv", "dest": "qc/ataqv"},
    {
        "kind": "suffix",
        "suffixes": [".json"],
        "dest": "qc/ataqv",
        "exclude_substrings": ["fastp_report"],
    },
    {
        "kind": "suffix",
        "suffixes": [".aligned.final.bam", ".aligned.final.bam.bai"],
        "dest": "bams",
    },
    {"kind": "suffix", "suffixes": [".bw"], "dest": "bigwigs"},
    {
        "kind": "suffix",
        "suffixes": [".bb"],
        "dest": "bigbeds",
        "exclude_substrings": ["ref/"],
    },
    {
        "kind": "suffix",
        "suffixes": [".narrowPeak.gz", ".summits.bed.gz"],
        "dest": "peaks",
    },
    {"kind": "dir", "dir_name": "results/tn5_motif", "dest": "tn5_counts"},
    {"kind": "dir", "dir_name": "results/deseq2", "dest": "deseq2"},
]


def match_rule(relpath: str, rule: dict) -> Optional[str]:
    """Check if a relative path matches a rule and return the S3 destination."""
    kind = rule.get("kind")

    if kind == "path":
        if relpath == rule.get("path"):
            return rule.get("dest")

    elif kind == "suffix":
        for suffix in rule.get("suffixes", []):
            if relpath.endswith(suffix):
                if any(s in relpath for s in rule.get("exclude_substrings", [])):
                    return None
                dest_dir = rule.get("dest", "")
                base = os.path.basename(relpath)
                return f"{dest_dir}/{base}" if dest_dir else base

    elif kind == "dir":
        dir_name = rule.get("dir_name")
        if not dir_name:
            return None
        if any(s in relpath for s in rule.get("exclude_substrings", [])):
            return None
        dir_path = dir_name.strip("/")
        marker = f"/{dir_path}/"
        relpath_with_leading_slash = f"/{relpath}"
        idx = relpath_with_leading_slash.find(marker)
        if idx == -1:
            return None
        subpath = relpath_with_leading_slash[idx + len(marker) :]
        if not subpath:
            return None
        dest_dir = rule.get("dest", dir_path)
        return f"{dest_dir}/{subpath}"

    return None


def gather_files(workdir: Path) -> List[Tuple[str, str]]:
    """Walk workdir and match files to rules, returning list of (src_relpath, s3_dest) tuples."""
    entries = []
    seen = {}

    for root, dirs, files in os.walk(workdir):
        dirs.sort()
        files.sort()
        for fname in files:
            full_path = os.path.join(root, fname)
            relpath = os.path.relpath(full_path, workdir)

            for rule in CHROMA2_RULES:
                dest = match_rule(relpath, rule)
                if dest:
                    if relpath in seen and seen[relpath] != dest:
                        print(
                            f"Warning: Conflicting destinations for {relpath}: {seen[relpath]} vs {dest}",
                            file=sys.stderr,
                        )
                    if relpath not in seen:
                        seen[relpath] = dest
                        entries.append((relpath, dest))
                    break

    return entries


def get_storage_class(dest_path: str, default_class: str, large_file_class: str) -> str:
    """Assign storage class based on file type."""
    if dest_path.endswith((".bam", ".bai")):
        return large_file_class
    return default_class


def run_transfer(
    workdir: Path,
    pipeline_name: str,
    sample_set: str,
    bucket: str,
    s3_prefix: str,
    storage_class: str,
    large_file_storage_class: str,
    dry_run: bool = False,
) -> int:
    """Transfer files to S3. Returns number of successful transfers."""
    entries = gather_files(workdir)

    if not entries:
        print("No files matched Chroma 2 transfer rules", file=sys.stderr)
        return 0

    success_count = 0
    failed_count = 0

    for src_relpath, s3_dest_path in entries:
        src_file = workdir / src_relpath
        s3_path = (
            f"s3://{bucket}/{s3_prefix}/{pipeline_name}/{sample_set}/{s3_dest_path}"
        )
        storage_cls = get_storage_class(
            s3_dest_path, storage_class, large_file_storage_class
        )

        cmd = [
            "aws",
            "s3",
            "cp",
            str(src_file),
            s3_path,
            "--storage-class",
            storage_cls,
        ]

        if dry_run:
            print(f"[DRY RUN] {' '.join(cmd)}")
            success_count += 1
        else:
            print(f"Transferring {src_relpath} → {s3_path}")
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode == 0:
                success_count += 1
            else:
                print(f"Error transferring {src_file}:", file=sys.stderr)
                print(result.stderr, file=sys.stderr)
                failed_count += 1

    total = success_count + failed_count
    print(f"\nTransfer complete: {success_count}/{total} files transferred")

    return 0 if failed_count == 0 else 1


def main():
    parser = argparse.ArgumentParser(
        description="Transfer Chroma 2 pipeline outputs to S3",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--workdir",
        type=Path,
        default=Path.cwd(),
        help="Pipeline working directory (default: current directory)",
    )
    parser.add_argument(
        "--pipeline-name",
        default="CHROMA",
        help="Pipeline name for S3 path (default: CHROMA)",
    )
    parser.add_argument(
        "--sample-set-name",
        required=True,
        help="Sample set name (used in S3 path)",
    )
    parser.add_argument(
        "--bucket",
        default="dremel-lab-bucket",
        help="S3 bucket name (default: dremel-lab-bucket)",
    )
    parser.add_argument(
        "--s3-prefix",
        default="_HTS",
        help="S3 prefix within bucket (default: _HTS)",
    )
    parser.add_argument(
        "--storage-class",
        default="GLACIER_IR",
        help="Storage class for regular files (default: GLACIER_IR)",
    )
    parser.add_argument(
        "--large-file-storage-class",
        default="GLACIER",
        help="Storage class for BAM/BAI files (default: GLACIER)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show transfers without executing",
    )

    args = parser.parse_args()

    workdir = args.workdir.resolve()
    if not workdir.exists():
        print(f"Error: workdir not found: {workdir}", file=sys.stderr)
        return 1

    return run_transfer(
        workdir=workdir,
        pipeline_name=args.pipeline_name,
        sample_set=args.sample_set_name,
        bucket=args.bucket,
        s3_prefix=args.s3_prefix,
        storage_class=args.storage_class,
        large_file_storage_class=args.large_file_storage_class,
        dry_run=args.dry_run,
    )


if __name__ == "__main__":
    sys.exit(main())
