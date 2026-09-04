#!/usr/bin/env bash
# Shared helper for bamcoverage_host/bamcoverage_virus: prints the
# --scaleFactor argument (or nothing) to add to a bamCoverage command line,
# based on the precomputed per-sample factor in norm_factors.tsv.
#
# Usage: get_bw_scale_args.sh <rule_label> <sample> <use_scaling: 0|1> <norm_factors.tsv>
set -euo pipefail

rule_label="$1"
sample="$2"
use_scaling="$3"
norm_factors="$4"

if [ "$use_scaling" != "1" ]; then
    exit 0
fi

sf=$(awk -F '\t' -v s="$sample" '$1==s{print $2}' "$norm_factors")
if [ -n "$sf" ]; then
    echo "[${rule_label}] sample=${sample} scaleFactor=${sf}" >&2
    echo "--scaleFactor $sf"
else
    echo "[${rule_label}] sample=${sample} has no precomputed norm factor; generating unscaled bigwig" >&2
fi
