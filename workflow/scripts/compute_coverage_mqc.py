#!/usr/bin/env python3
import pysam
from bisect import bisect_left
from pathlib import Path
from os.path import dirname

input_bam = snakemake.input.bam
mapqs = sorted(snakemake.params.mapqs)
output_file = snakemake.output[0]

Path(dirname(output_file)).mkdir(parents=True, exist_ok=True)

depth_sum = {mapq: 0 for mapq in mapqs}
covered_positions = {mapq: 0 for mapq in mapqs}

with pysam.AlignmentFile(input_bam, "rb") as bam:
    for pileup_column in bam.pileup():
        if pileup_column.n == 0:
            continue
        read_mapqs = sorted(
            read.alignment.mapping_quality for read in pileup_column.pileups
        )
        for mapq in mapqs:
            high_quality_reads = len(read_mapqs) - bisect_left(read_mapqs, mapq)
            if high_quality_reads > 0:
                depth_sum[mapq] += high_quality_reads
                covered_positions[mapq] += 1

with open(output_file, "w") as out:
    out.write("mapq\tmean_depth\n")
    for mapq in mapqs:
        mean_depth = (
            depth_sum[mapq] / covered_positions[mapq] if covered_positions[mapq] else 0
        )
        out.write(f"{mapq}\t{mean_depth}\n")
