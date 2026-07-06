#!/usr/bin/env python3
import pandas as pd
from os.path import dirname, basename
import os

output_file = snakemake.output[0]
input_files = snakemake.input

os.makedirs(dirname(output_file), exist_ok=True)

KNOWN_SAMPLES = set(snakemake.params.samples)

counts = {}  # sample -> category -> total tn5_site_count
categories = set()

for input_file in input_files:
    category = basename(dirname(input_file))
    categories.add(category)

    df = pd.read_csv(input_file, sep="\t", header=0)
    sample_columns = [c for c in df.columns if c in KNOWN_SAMPLES]

    for sample in sample_columns:
        total = int(pd.to_numeric(df[sample], errors="coerce").fillna(0).sum())
        counts.setdefault(sample, {})[category] = total

sample_list = sorted(counts.keys())
category_list = sorted(categories)

output_data = []
for sample in sample_list:
    row = {"Sample": sample}
    for cat in category_list:
        row[cat] = counts[sample].get(cat, 0)
    output_data.append(row)

output_df = pd.DataFrame(output_data, columns=["Sample"] + category_list)

with open(output_file, "w") as f:
    f.write("# id: 'tn5_counts'\n")
    f.write("# section_name: 'Tn5 Cut Site Counts'\n")
    f.write("# plot_type: 'heatmap'\n")
    f.write("# pconfig:\n")
    f.write("#   id: 'tn5_counts_plot'\n")
    output_df.to_csv(f, sep="\t", index=False)
