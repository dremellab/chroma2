#!/usr/bin/env python3
import pandas as pd
from os.path import dirname
import os
import gzip
import numpy as np

from mqc_common import sample_name_from_path, write_mqc_header

output_file = snakemake.output[0]
input_files = snakemake.input

os.makedirs(dirname(output_file), exist_ok=True)

BIN_SIZE = 50
MAX_SIZE = 2000

bins = list(range(0, MAX_SIZE + BIN_SIZE, BIN_SIZE))
bin_labels = [f"{b}-{b+BIN_SIZE}" for b in bins[:-1]]

peak_data = {}

for input_file in input_files:
    sample_name = sample_name_from_path(input_file)

    widths = []
    with gzip.open(input_file, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                start = int(parts[1])
                end = int(parts[2])
                width = end - start
                if width <= MAX_SIZE:
                    widths.append(width)

    hist, _ = np.histogram(widths, bins=bins)
    peak_data[sample_name] = hist

output_df = pd.DataFrame.from_dict(peak_data, orient="index", columns=bin_labels)
output_df.index.name = "Sample"

with open(output_file, "w") as f:
    write_mqc_header(
        f,
        "peak_size_dist",
        "Peak Size Distribution",
        "bargraph",
        {"id": "peak_size_plot"},
    )
    output_df.to_csv(f, sep="\t")
