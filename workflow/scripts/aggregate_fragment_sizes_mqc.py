#!/usr/bin/env python3
"""
Aggregate per-sample fragment size distributions into a single MultiQC-formatted TSV.
Pivots data so bins are rows and samples are columns.
"""
import pandas as pd
import os
from pathlib import Path

output_file = snakemake.output[0]
input_files = snakemake.input

# Expected bin labels (must match extract_fragment_sizes.py)
EXPECTED_BINS = [
    "0-50bp",
    "50-100bp",
    "100-150bp",
    "150-200bp",
    "200-300bp",
    "300-500bp",
    "500+bp",
]

# Read all per-sample TSVs and extract data
data_frames = {}

for input_file in input_files:
    # Extract sample name from filename: {sample}.host.fragment_sizes.tsv
    filename = os.path.basename(input_file)
    sample_name = filename.split(".")[0]

    try:
        # Read TSV, skipping comments
        df = pd.read_csv(input_file, sep="\t", comment="#")

        # Extract just the Fragment_Size_Range and Count columns
        if "Fragment_Size_Range" in df.columns and "Count" in df.columns:
            # Keep only the data rows (filter out any extra rows)
            df = df[df["Fragment_Size_Range"].isin(EXPECTED_BINS)].copy()
            df = df.set_index("Fragment_Size_Range")

            # Rename Count column to sample name
            df = df.rename(columns={"Count": sample_name})
            df[sample_name] = df[sample_name].astype(int)

            data_frames[sample_name] = df[[sample_name]]
        else:
            print(f"Warning: {input_file} missing expected columns, skipping")

    except Exception as e:
        print(f"Error reading {input_file}: {e}")
        continue

# Combine all sample dataframes
if data_frames:
    combined_df = pd.concat(data_frames.values(), axis=1)

    # Reorder columns by sample name for consistency
    combined_df = combined_df[sorted(combined_df.columns)]

    # Reorder rows to match expected bins
    combined_df = combined_df.reindex(EXPECTED_BINS)

    # Write output with MultiQC headers
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with open(output_file, "w") as f:
        # MultiQC headers for bargraph
        f.write("# id: 'fragment_size_aggregate'\n")
        f.write("# section_name: 'Fragment Size Distribution'\n")
        f.write("# plot_type: 'bargraph'\n")
        f.write("# pconfig:\n")
        f.write("#   id: 'fragment_size_aggregate_plot'\n")
        f.write("#   title: 'Fragment Size Distribution Across Samples'\n")
        f.write("#   ylab: 'Number of Fragments'\n")
        f.write("#   xlab: 'Fragment Size Range'\n")
        f.write("#   stacking: 'normal'\n")

        # Write data
        combined_df.to_csv(f, sep="\t")

    print(f"Aggregated fragment size data from {len(data_frames)} samples")
    print(f"Output written to {output_file}")
    print(f"\nFragment size distribution summary:")
    print(combined_df)

else:
    print("Error: No valid input files found")
    exit(1)
