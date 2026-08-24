#!/usr/bin/env python3
import pandas as pd
from os.path import dirname
import os

from mqc_common import sample_name_from_path, write_mqc_header

output_file = snakemake.output[0]
input_files = snakemake.input

os.makedirs(dirname(output_file), exist_ok=True)

coverage_data = {}
mapq_thresholds = [0, 10, 20, 30]

for input_file in input_files:
    sample_name = sample_name_from_path(input_file)
    df = pd.read_csv(input_file, sep="\t")
    coverage_data[sample_name] = dict(zip(df["mapq"], df["mean_depth"]))

sample_list = sorted(coverage_data.keys())

output_data = []
for sample in sample_list:
    row = {"Sample": sample}
    for q in mapq_thresholds:
        row[f"MAPQ_{q}"] = coverage_data[sample].get(q, 0)
    output_data.append(row)

output_df = pd.DataFrame(output_data)

with open(output_file, "w") as f:
    write_mqc_header(
        f,
        "genome_coverage_mapq",
        "Genome Coverage at MAPQ Thresholds",
        "table",
        {"id": "genome_coverage_table"},
    )
    output_df.to_csv(f, sep="\t", index=False)
