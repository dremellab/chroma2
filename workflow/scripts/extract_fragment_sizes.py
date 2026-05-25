#!/usr/bin/env python3
"""
Extract paired-end fragment size statistics from a BAM file.
Outputs per-sample fragment size distribution binned by size ranges.
"""
import pysam
import sys
from collections import defaultdict

bam_file = snakemake.input.bam
output_file = snakemake.output[0]
sample_name = snakemake.wildcards.sample

# Fragment size bins: 0-50, 50-100, 100-150, 150-200, 200-300, 300-500, 500+
BIN_EDGES = [0, 50, 100, 150, 200, 300, 500, float("inf")]
BIN_LABELS = [
    "0-50bp",
    "50-100bp",
    "100-150bp",
    "150-200bp",
    "200-300bp",
    "300-500bp",
    "500+bp",
]

# Count fragments in each bin
bin_counts = defaultdict(int)
total_fragments = 0
fragment_sizes = []

# Read BAM file and extract fragment sizes
try:
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            # Only count properly paired reads (first in pair)
            if read.is_read1 and read.is_proper_pair and not read.is_unmapped:
                insert_size = abs(read.template_length)

                # Skip zero insert sizes
                if insert_size > 0:
                    fragment_sizes.append(insert_size)
                    total_fragments += 1

                    # Determine which bin this fragment belongs to
                    for i, edge in enumerate(BIN_EDGES[1:]):
                        if insert_size < edge:
                            bin_counts[BIN_LABELS[i]] += 1
                            break
except Exception as e:
    print(f"Error reading BAM file: {e}", file=sys.stderr)
    sys.exit(1)

# Calculate summary statistics
if fragment_sizes:
    fragment_sizes.sort()
    mean_size = sum(fragment_sizes) / len(fragment_sizes)
    median_size = fragment_sizes[len(fragment_sizes) // 2]
    min_size = fragment_sizes[0]
    max_size = fragment_sizes[-1]
    p25 = fragment_sizes[len(fragment_sizes) // 4]
    p75 = fragment_sizes[3 * len(fragment_sizes) // 4]
else:
    mean_size = median_size = min_size = max_size = p25 = p75 = 0

# Write output TSV with per-bin counts and summary statistics
with open(output_file, "w") as f:
    # MultiQC headers for bargraph
    f.write("# id: 'fragment_size_dist'\n")
    f.write("# section_name: 'Fragment Size Distribution'\n")
    f.write("# plot_type: 'bargraph'\n")
    f.write("# pconfig:\n")
    f.write("#   id: 'fragment_size_plot'\n")
    f.write("#   title: 'Fragment Size Distribution'\n")
    f.write("#   ylab: 'Number of Fragments'\n")
    f.write("#   xlab: 'Fragment Size Range'\n")

    # Write data rows
    f.write("Fragment_Size_Range\tCount\n")
    for bin_label in BIN_LABELS:
        count = bin_counts[bin_label]
        f.write(f"{bin_label}\t{count}\n")

    # Add summary statistics as comments for reference
    f.write(f"\n# Summary Statistics for {sample_name}\n")
    f.write(f"# Total Fragments: {total_fragments}\n")
    f.write(f"# Mean Fragment Size: {mean_size:.1f}bp\n")
    f.write(f"# Median Fragment Size: {median_size}bp\n")
    f.write(f"# Min Fragment Size: {min_size}bp\n")
    f.write(f"# Max Fragment Size: {max_size}bp\n")
    f.write(f"# 25th Percentile: {p25}bp\n")
    f.write(f"# 75th Percentile: {p75}bp\n")

print(f"Extracted {total_fragments} fragments from {bam_file}")
print(f"Output written to {output_file}")
