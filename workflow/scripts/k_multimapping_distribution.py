#!/usr/bin/env python3
"""
Compute the distribution of reported alignments per read from a BAM file.

This is intended for Bowtie2 -k style outputs, where each read may have
up to k alignments reported. The script counts how many alignments are
reported per read and outputs a simple distribution table.
"""

from __future__ import annotations

import argparse
import collections
import os
import sys
import tempfile
from typing import Dict, Iterable, Iterator, Optional

import pysam


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Compute distribution of number of reported alignments per read "
            "from a BAM file."
        )
    )
    p.add_argument("bam", help="Input BAM file (name-sorted)")
    p.add_argument(
        "-o",
        "--output",
        help="Output TSV path (default: stdout)",
        default="-",
    )
    p.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        help="Threads for BAM reading (default: 1)",
    )
    p.add_argument(
        "--tmpdir",
        default="/tmp",
        help="Temporary directory for coordinate sorting (default: /tmp)",
    )
    p.add_argument(
        "--include-unmapped",
        action="store_true",
        help="Include reads with zero reported alignments as n_alignments=0",
    )
    p.add_argument(
        "--proper-only",
        action="store_true",
        help="Count only alignments flagged as proper pairs (SAM flag 0x2)",
    )
    p.add_argument(
        "--bedgraph",
        default=None,
        help="Write bedGraph of average per-read k over collapsed alignments",
    )
    p.add_argument(
        "--chrom-sizes",
        default=None,
        help="Write chrom sizes from BAM header (TSV: chrom\\tlength)",
    )
    return p.parse_args()


def bam_alignments(bam: str, threads: int) -> Iterator[pysam.AlignedSegment]:
    with pysam.AlignmentFile(bam, "rb", threads=threads) as fh:
        for rec in fh.fetch(until_eof=True):
            yield rec


def _should_count(rec: pysam.AlignedSegment, proper_only: bool) -> bool:
    if rec.is_unmapped:
        return False
    if proper_only and (not rec.is_proper_pair):
        return False
    return True


def _group_k(
    paired_seen: bool,
    count_any: int,
    count_r1: int,
    count_r2: int,
) -> int:
    if paired_seen:
        return count_r1 if count_r1 > 0 else count_r2
    return count_any


def count_alignments_per_read_sorted(
    records: Iterable[pysam.AlignedSegment],
    include_unmapped: bool,
    proper_only: bool,
) -> Dict[int, int]:
    dist: Dict[int, int] = collections.Counter()

    current_qname: Optional[str] = None
    mapped_count = 0
    mapped_r1 = 0
    mapped_r2 = 0
    paired_seen = False
    saw_any = False

    def flush() -> None:
        nonlocal mapped_count, mapped_r1, mapped_r2, paired_seen, saw_any
        if not saw_any:
            return
        k = _group_k(paired_seen, mapped_count, mapped_r1, mapped_r2)
        if k > 0:
            dist[k] += 1
        elif include_unmapped:
            dist[0] += 1
        mapped_count = 0
        mapped_r1 = 0
        mapped_r2 = 0
        paired_seen = False
        saw_any = False

    for rec in records:
        qname = rec.query_name

        if current_qname is None:
            current_qname = qname
        elif qname != current_qname:
            flush()
            current_qname = qname

        saw_any = True
        if rec.is_paired:
            paired_seen = True
        if _should_count(rec, proper_only):
            mapped_count += 1
            if rec.is_read1:
                mapped_r1 += 1
            elif rec.is_read2:
                mapped_r2 += 1

    flush()
    return dist


def count_per_read(
    records: Iterable[pysam.AlignedSegment],
    proper_only: bool,
) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    paired: Dict[str, bool] = {}
    r1_counts: Dict[str, int] = {}
    r2_counts: Dict[str, int] = {}

    for rec in records:
        qname = rec.query_name
        if qname not in counts:
            counts[qname] = 0
        if rec.is_paired:
            paired[qname] = True
        if _should_count(rec, proper_only):
            counts[qname] += 1
            if rec.is_read1:
                r1_counts[qname] = r1_counts.get(qname, 0) + 1
            elif rec.is_read2:
                r2_counts[qname] = r2_counts.get(qname, 0) + 1

    if paired:
        for qname in counts.keys():
            if paired.get(qname, False):
                r1 = r1_counts.get(qname, 0)
                r2 = r2_counts.get(qname, 0)
                counts[qname] = r1 if r1 > 0 else r2
    return counts


def dist_from_counts(counts: Dict[str, int], include_unmapped: bool) -> Dict[int, int]:
    dist: Dict[int, int] = collections.Counter()
    for _, cnt in counts.items():
        if cnt > 0 or include_unmapped:
            dist[cnt] += 1
    return dist


def count_alignments_per_read_unsorted(
    records: Iterable[pysam.AlignedSegment],
    include_unmapped: bool,
    proper_only: bool,
) -> Dict[int, int]:
    counts = count_per_read(records, proper_only)
    return dist_from_counts(counts, include_unmapped)


def write_bedgraph(
    bam: str,
    counts: Dict[str, int],
    out_path: str,
    threads: int,
    proper_only: bool,
    tmpdir: str,
) -> None:
    if out_path is None:
        return
    out_fh = sys.stdout if out_path == "-" else open(out_path, "w", encoding="utf-8")
    tmp_bam: Optional[str] = None
    try:
        with pysam.AlignmentFile(bam, "rb", threads=threads) as fh:
            so = fh.header.get("HD", {}).get("SO")
        if so not in (None, "coordinate"):
            fd, tmp_bam = tempfile.mkstemp(
                prefix="k_multimapping.coord.", suffix=".bam", dir=tmpdir
            )
            os.close(fd)
            pysam.sort("-o", tmp_bam, "-@", str(threads), bam)
            bam_for_bg = tmp_bam
        else:
            bam_for_bg = bam

        with pysam.AlignmentFile(bam_for_bg, "rb", threads=threads) as fh:
            current_tid: Optional[int] = None
            events: list[tuple[int, int, int]] = []

            def flush_events(tid: int, evs: list[tuple[int, int, int]]) -> None:
                if not evs:
                    return
                evs.sort()
                ref_name = fh.get_reference_name(tid)
                sum_k = 0
                count = 0
                i = 0
                last_pos: Optional[int] = None
                last_out_start: Optional[int] = None
                last_out_end: Optional[int] = None
                last_out_val: Optional[float] = None

                def emit_segment(start: int, end: int, val: float) -> None:
                    out_fh.write(f"{ref_name}\t{start}\t{end}\t{val:.6f}\n")

                while i < len(evs):
                    pos = evs[i][0]
                    if last_pos is not None and pos > last_pos and count > 0:
                        val = sum_k / count
                        if (
                            last_out_val is not None
                            and last_out_end == last_pos
                            and abs(last_out_val - val) < 1e-9
                        ):
                            last_out_end = pos
                        else:
                            if last_out_val is not None:
                                emit_segment(last_out_start, last_out_end, last_out_val)
                            last_out_start = last_pos
                            last_out_end = pos
                            last_out_val = val
                    dsum = 0
                    dcnt = 0
                    while i < len(evs) and evs[i][0] == pos:
                        dsum += evs[i][1]
                        dcnt += evs[i][2]
                        i += 1
                    sum_k += dsum
                    count += dcnt
                    last_pos = pos
                if last_out_val is not None:
                    emit_segment(last_out_start, last_out_end, last_out_val)

            for rec in fh.fetch(until_eof=True):
                if not _should_count(rec, proper_only):
                    continue
                k = counts.get(rec.query_name, 0)
                if k <= 0:
                    continue
                if rec.reference_end is None:
                    continue
                tid = rec.reference_id
                if current_tid is None:
                    current_tid = tid
                if tid != current_tid:
                    flush_events(current_tid, events)
                    events = []
                    current_tid = tid
                start = rec.reference_start
                end = rec.reference_end
                if end <= start:
                    continue
                events.append((start, k, 1))
                events.append((end, -k, -1))

            if current_tid is not None:
                flush_events(current_tid, events)
    finally:
        if out_fh is not sys.stdout:
            out_fh.close()
        if tmp_bam is not None and os.path.exists(tmp_bam):
            os.remove(tmp_bam)


def write_tsv(dist: Dict[int, int], out_path: str) -> None:
    total_reads = sum(dist.values())
    lines = []
    lines.append("n_alignments\tread_count\tfraction")
    for n in sorted(dist.keys()):
        count = dist[n]
        frac = (count / total_reads) if total_reads else 0.0
        lines.append(f"{n}\t{count}\t{frac:.6f}")

    if out_path == "-":
        sys.stdout.write("\n".join(lines) + "\n")
    else:
        with open(out_path, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")


def write_chrom_sizes(bam: str, out_path: str) -> None:
    if out_path is None:
        return
    out_fh = sys.stdout if out_path == "-" else open(out_path, "w", encoding="utf-8")
    try:
        with pysam.AlignmentFile(bam, "rb") as fh:
            seen = set()
            for rec in fh.fetch(until_eof=True):
                if rec.is_unmapped:
                    continue
                if rec.reference_id is not None and rec.reference_id >= 0:
                    seen.add(rec.reference_id)
            for tid in sorted(seen):
                out_fh.write(f"{fh.references[tid]}\t{fh.lengths[tid]}\n")
    finally:
        if out_fh is not sys.stdout:
            out_fh.close()


def main() -> None:
    args = parse_args()
    if args.bedgraph:
        counts = count_per_read(
            bam_alignments(args.bam, args.threads), args.proper_only
        )
        dist = dist_from_counts(counts, args.include_unmapped)
        write_bedgraph(
            args.bam,
            counts,
            args.bedgraph,
            args.threads,
            args.proper_only,
            args.tmpdir,
        )
    else:
        records = bam_alignments(args.bam, args.threads)
        dist = count_alignments_per_read_sorted(
            records, args.include_unmapped, args.proper_only
        )
    write_tsv(dist, args.output)
    if args.chrom_sizes:
        write_chrom_sizes(args.bam, args.chrom_sizes)


if __name__ == "__main__":
    main()
