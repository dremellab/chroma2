#!/usr/bin/env bash
# Shared sort -> optional score-clamp -> bedToBigBed conversion used by every
# rule in bigbed.smk. Behavior (sort key, score-clamp expressions, empty-input
# handling, bedToBigBed flag order) matches what each of those rules ran
# inline before this script existed -- this is a pure extraction, not a
# behavior change.
set -euxo pipefail

usage() {
    cat >&2 <<'EOF'
Usage: bed_to_bigbed.sh --input F --output F --chromsizes F --tmpdir D \
         --bed-type bed3|bed5|bed6+4 [--gzipped] \
         [--clamp-score --clamp-mode conditional|unconditional] \
         [--autosql F] [--tab]
EOF
    exit 2
}

input=""
output=""
chromsizes=""
tmpdir=""
bed_type=""
gzipped=0
clamp_score=0
clamp_mode="unconditional"
autosql=""
tab=0

while [ $# -gt 0 ]; do
    case "$1" in
        --input) input="$2"; shift 2 ;;
        --output) output="$2"; shift 2 ;;
        --chromsizes) chromsizes="$2"; shift 2 ;;
        --tmpdir) tmpdir="$2"; shift 2 ;;
        --bed-type) bed_type="$2"; shift 2 ;;
        --gzipped) gzipped=1; shift ;;
        --clamp-score) clamp_score=1; shift ;;
        --clamp-mode) clamp_mode="$2"; shift 2 ;;
        --autosql) autosql="$2"; shift 2 ;;
        --tab) tab=1; shift ;;
        *) echo "Unknown argument: $1" >&2; usage ;;
    esac
done

for required in input output chromsizes tmpdir bed_type; do
    if [ -z "${!required}" ]; then
        echo "Missing required argument: --${required//_/-}" >&2
        usage
    fi
done

mkdir -p "$tmpdir"
sorted_bed="$tmpdir/$(basename "$output" .bb).sorted.bed"

if [ "$gzipped" -eq 1 ]; then
    reader=(zcat "$input")
else
    reader=(cat "$input")
fi

if [ "$clamp_score" -eq 1 ]; then
    if [ "$clamp_mode" = "conditional" ]; then
        awk_prog='BEGIN{OFS="\t"} {if (NF >= 5) {score=$5; if (score == ".") score=0; if (score < 0) score=0; if (score > 1000) score=1000; $5=int(score)} print}'
    elif [ "$clamp_mode" = "unconditional" ]; then
        awk_prog='BEGIN{OFS="\t"} {score=$5; if (score == ".") score=0; if (score < 0) score=0; if (score > 1000) score=1000; $5=int(score); print}'
    else
        echo "Unknown --clamp-mode: $clamp_mode (expected conditional or unconditional)" >&2
        exit 2
    fi
    "${reader[@]}" | awk "$awk_prog" | sort -k1,1 -k2,2n > "$sorted_bed"
else
    "${reader[@]}" | sort -k1,1 -k2,2n > "$sorted_bed"
fi

if [ ! -s "$sorted_bed" ]; then
    touch "$output"
    exit 0
fi

bedtobigbed_args=()
if [ "$tab" -eq 1 ]; then
    bedtobigbed_args+=(-tab)
fi
if [ -n "$autosql" ]; then
    bedtobigbed_args+=(-as="$autosql")
fi
bedtobigbed_args+=(-type="$bed_type")

bedToBigBed "${bedtobigbed_args[@]}" "$sorted_bed" "$chromsizes" "$output"
