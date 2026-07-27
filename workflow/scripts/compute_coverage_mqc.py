#!/usr/bin/env python3
import pysam
from bisect import bisect_left
from pathlib import Path
from os.path import dirname

input_bam = snakemake.input.bam
outpfx = snakemake.params.outpfx
mapqs = sorted(snakemake.params.mapqs)

Path(dirname(outpfx)).mkdir(parents=True, exist_ok=True)

header = (
    "#rname\tstartpos\tendpos\tnumreads\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n"
)

with pysam.AlignmentFile(input_bam, "rb") as bam:
    outfiles = {mapq: open(f"{outpfx}.coverage_q{mapq}.txt", "w") for mapq in mapqs}
    for out in outfiles.values():
        out.write(header)
    try:
        for pileup_column in bam.pileup():
            if pileup_column.n == 0:
                continue
            read_mapqs = sorted(
                read.alignment.mapping_quality for read in pileup_column.pileups
            )
            counts = [len(read_mapqs) - bisect_left(read_mapqs, mapq) for mapq in mapqs]
            if not any(counts):
                continue
            chrom = bam.references[pileup_column.reference_id]
            startpos = pileup_column.reference_pos
            endpos = startpos + 1
            for mapq, high_quality_reads in zip(mapqs, counts):
                if high_quality_reads > 0:
                    outfiles[mapq].write(
                        f"{chrom}\t{startpos}\t{endpos}\t"
                        f"{high_quality_reads}\t1.0\t{high_quality_reads}\t0\t30\n"
                    )
    finally:
        for out in outfiles.values():
            out.close()
