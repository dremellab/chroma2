#!/usr/bin/env python3
import pandas as pd
from os.path import dirname
import os

input_file = snakemake.input.summary
stats_out = snakemake.output.stats
ratio_out = snakemake.output.ratio

os.makedirs(dirname(stats_out), exist_ok=True)
os.makedirs(dirname(ratio_out), exist_ok=True)

df = pd.read_csv(input_file, sep="\t")

steps_short = {
    "raw": "aligned",
    "no_multimappers": "clean",
    "dedup": "dedup",
    "mapqfiltered": "final",
}
step_names = ["aligned", "clean", "dedup", "final"]

stats_data = []
for _, row in df.iterrows():
    sample = row["sample"]
    stats_row = {"Sample": sample}

    for col_step, step_name in steps_short.items():
        for col in df.columns:
            if col.endswith(f"_{col_step}"):
                count = row[col]
                try:
                    stats_row[f"{step_name}_reads"] = int(count)
                except (ValueError, TypeError):
                    stats_row[f"{step_name}_reads"] = 0
                break

    stats_data.append(stats_row)

stats_df = pd.DataFrame(stats_data)
stats_df = stats_df[
    ["Sample"]
    + [col for col in stats_df.columns if col != "Sample" and col.endswith("_reads")]
]

with open(stats_out, "w") as f:
    f.write("# id: 'alignment_stats'\n")
    f.write("# section_name: 'Alignment Statistics'\n")
    f.write("# plot_type: 'table'\n")
    f.write("# pconfig:\n")
    f.write("#   id: 'alignment_stats_table'\n")
    stats_df.to_csv(f, sep="\t", index=False)

ratio_data = []
for _, row in df.iterrows():
    sample = row["sample"]
    host_count = 0
    virus_count = 0

    host_final_col = "host_mapqfiltered"
    if host_final_col in df.columns:
        try:
            host_count = int(row[host_final_col])
        except (ValueError, TypeError):
            host_count = 0

    for col in df.columns:
        if col.endswith("_mapqfiltered") and not col.startswith(("host", "chrM")):
            try:
                virus_count += int(row[col])
            except (ValueError, TypeError):
                pass

    ratio_data.append({"Sample": sample, "Host": host_count, "Virus": virus_count})

ratio_df = pd.DataFrame(ratio_data)

with open(ratio_out, "w") as f:
    f.write("# id: 'host_virus_ratio'\n")
    f.write("# section_name: 'Host vs. Virus Read Distribution'\n")
    f.write("# plot_type: 'bargraph'\n")
    f.write("# pconfig:\n")
    f.write("#   id: 'host_virus_ratio_plot'\n")
    ratio_df.to_csv(f, sep="\t", index=False)
