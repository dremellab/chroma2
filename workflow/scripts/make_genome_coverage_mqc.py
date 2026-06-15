#!/usr/bin/env python3
import pandas as pd
from os.path import dirname, basename
import os

output_file = snakemake.output[0]
input_files = snakemake.input

os.makedirs(dirname(output_file), exist_ok=True)

coverage_data = {}
mapq_thresholds = [0, 10, 20, 30]

for input_file in input_files:
    basename_file = basename(input_file)
    sample_name = basename_file.split(".")[0]
    mapq_match = [q for q in mapq_thresholds if f"_q{q}" in basename_file]
    if not mapq_match:
        continue
    mapq = mapq_match[0]

    df = pd.read_csv(input_file, sep="\t", skiprows=0)
    if df.shape[0] > 0:
        mean_depth = df["meandepth"].mean()
    else:
        mean_depth = 0

    if sample_name not in coverage_data:
        coverage_data[sample_name] = {}
    coverage_data[sample_name][f"MAPQ_{mapq}"] = mean_depth

sample_list = sorted(coverage_data.keys())
mapq_cols = [f"MAPQ_{q}" for q in sorted(mapq_thresholds)]

output_data = []
for sample in sample_list:
    row = {"Sample": sample}
    for col in mapq_cols:
        row[col] = coverage_data[sample].get(col, 0)
    output_data.append(row)

output_df = pd.DataFrame(output_data)

with open(output_file, "w") as f:
    f.write("# id: 'genome_coverage_mapq'\n")
    f.write("# section_name: 'Genome Coverage at MAPQ Thresholds'\n")
    f.write("# plot_type: 'table'\n")
    f.write("# pconfig:\n")
    f.write("#   id: 'genome_coverage_table'\n")
    output_df.to_csv(f, sep="\t", index=False)
