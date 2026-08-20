#!/usr/bin/env python3
"""
Aggregate per-sample fragment size distributions into a single MultiQC-formatted TSV.
Rows are samples, columns are bin size ranges (so samples appear on the axis and
bin ranges appear as the stacked/colored series).
"""
import pandas as pd
import os

from mqc_common import sample_name_from_path, write_mqc_header

output_file = snakemake.output[0]
input_files = snakemake.input

# Expected bin labels, in size order (must match extract_fragment_sizes.py)
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
sample_rows = {}

for input_file in input_files:
    # Extract sample name from filename: {sample}.host.fragment_sizes.tsv
    sample_name = sample_name_from_path(input_file)

    try:
        # Read TSV, skipping comments
        df = pd.read_csv(input_file, sep="\t", comment="#")

        if "Fragment_Size_Range" in df.columns and "Count" in df.columns:
            counts = df.set_index("Fragment_Size_Range")["Count"]
            sample_rows[sample_name] = {
                bin_label: int(counts.get(bin_label, 0)) for bin_label in EXPECTED_BINS
            }
        else:
            print(f"Warning: {input_file} missing expected columns, skipping")

    except Exception as e:
        print(f"Error reading {input_file}: {e}")
        continue

# Combine all sample rows
if sample_rows:
    combined_df = pd.DataFrame.from_dict(
        sample_rows, orient="index", columns=EXPECTED_BINS
    )
    combined_df.index.name = "Sample"

    # Sort samples for consistent display; bin columns stay in size order
    combined_df = combined_df.sort_index()

    # Write output with MultiQC headers
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with open(output_file, "w") as f:
        write_mqc_header(
            f,
            "fragment_size_aggregate",
            "Fragment Size Distribution",
            "bargraph",
            {
                "id": "fragment_size_aggregate_plot",
                "title": "Fragment Size Distribution Across Samples",
                "ylab": "Number of Fragments",
                "xlab": "Sample",
                "stacking": "normal",
            },
        )

        # Write data
        combined_df.to_csv(f, sep="\t")

    print(f"Aggregated fragment size data from {len(sample_rows)} samples")
    print(f"Output written to {output_file}")
    print(f"\nFragment size distribution summary:")
    print(combined_df)

else:
    print("Error: No valid input files found")
    exit(1)
