#!/usr/bin/env python3
import pandas as pd
from os.path import dirname, basename
import os
import glob

output_file = snakemake.output[0]
input_files = snakemake.input

os.makedirs(dirname(output_file), exist_ok=True)

sample_counts = {}

for input_file in input_files:
    basename_file = basename(input_file)
    df = pd.read_csv(input_file, sep="\t", header=0)

    sample_name = basename_file.split(".")[0]
    matrix_type = "unknown"
    if "genrich" in basename_file:
        matrix_type = "genrich"
    elif "macs2" in basename_file:
        matrix_type = "macs2"
    if "virus" in basename_file:
        matrix_type += "_virus"
    elif "trna" in basename_file:
        matrix_type += "_trna"
    else:
        matrix_type += "_host"

    numeric_df = df.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
    total_count = int(numeric_df.sum().sum())

    if sample_name not in sample_counts:
        sample_counts[sample_name] = {}
    sample_counts[sample_name][matrix_type] = total_count

sample_list = sorted(sample_counts.keys())
matrix_types = set()
for counts in sample_counts.values():
    matrix_types.update(counts.keys())
matrix_types = sorted(matrix_types)

output_data = []
for sample in sample_list:
    row = {"Sample": sample}
    for mt in matrix_types:
        row[mt] = sample_counts[sample].get(mt, 0)
    output_data.append(row)

output_df = pd.DataFrame(output_data)

with open(output_file, "w") as f:
    f.write("# id: 'tn5_counts'\n")
    f.write("# section_name: 'Tn5 Cut Site Counts'\n")
    f.write("# plot_type: 'bargraph'\n")
    f.write("# pconfig:\n")
    f.write("#   id: 'tn5_counts_plot'\n")
    output_df.to_csv(f, sep="\t", index=False)
