#!/usr/bin/env python3
import pysam
from pathlib import Path
from os.path import dirname

input_bam = snakemake.input.bam
outpfx = snakemake.params.outpfx
mapqs = snakemake.params.mapqs

Path(dirname(outpfx)).mkdir(parents=True, exist_ok=True)

with pysam.AlignmentFile(input_bam, "rb") as bam:
    for mapq in mapqs:
        coverage_file = f"{outpfx}.coverage_q{mapq}.txt"
        with open(coverage_file, "w") as out:
            out.write(
                "#rname\tstartpos\tendpos\tnumreads\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n"
            )
            for pileup_column in bam.pileup():
                if pileup_column.n > 0:
                    high_quality_reads = sum(
                        1
                        for read in pileup_column.pileups
                        if read.alignment.mapping_quality >= mapq
                    )
                    if high_quality_reads > 0:
                        chrom = bam.references[pileup_column.reference_id]
                        out.write(
                            f"{chrom}\t{pileup_column.reference_pos}\t{pileup_column.reference_pos + 1}\t"
                            f"{high_quality_reads}\t1.0\t{high_quality_reads}\t0\t30\n"
                        )
