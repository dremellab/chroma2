###################################################################################
# Tn5 count matrix generation from BAM files
# Generates 6 types of matrices: gene, tRNA, Pol3, repeat elements, viral, rRNA
###################################################################################

# Extract count matrices configuration
COUNT_MATRICES = config.get("count_matrices", {})

# Parse enabled/disabled status for each matrix type
COUNT_MATRIX_TYPES = {
    "gene": COUNT_MATRICES.get("gene_counts", {}).get("enabled", True),
    "trna": COUNT_MATRICES.get("trna_counts", {}).get("enabled", True),
    "pol3": COUNT_MATRICES.get("pol3_counts", {}).get("enabled", True),
    "repeat_elements": COUNT_MATRICES.get(
        "repeat_element_counts", {}
    ).get("enabled", True),
    "viral": COUNT_MATRICES.get("viral_genome_counts", {}).get("enabled", True),
    "rrna": COUNT_MATRICES.get("rrna_counts", {}).get("enabled", True),
}

# Extract matrix-specific parameters
GENE_COUNTS_CONFIG = COUNT_MATRICES.get("gene_counts", {})
GENE_FLANK_SIZE = int(GENE_COUNTS_CONFIG.get("flank_size", 250))
GENE_EXCLUDE_FEATURES = GENE_COUNTS_CONFIG.get("exclude_features", [])

TRNA_COUNTS_CONFIG = COUNT_MATRICES.get("trna_counts", {})
TRNA_FLANK_SIZE = int(TRNA_COUNTS_CONFIG.get("flank_size", 100))
TRNA_MAX_SIZE = int(TRNA_COUNTS_CONFIG.get("max_size", 1000))

POL3_COUNTS_CONFIG = COUNT_MATRICES.get("pol3_counts", {})
POL3_FLANK_SIZE = int(POL3_COUNTS_CONFIG.get("flank_size", -1))
POL3_MAX_SIZE = int(POL3_COUNTS_CONFIG.get("max_size", 1000))

REPEAT_COUNTS_CONFIG = COUNT_MATRICES.get("repeat_element_counts", {})
REPEAT_FLANK_SIZE = int(REPEAT_COUNTS_CONFIG.get("flank_size", -1))
REPEAT_MAX_SIZE = int(REPEAT_COUNTS_CONFIG.get("max_size", 1000))

VIRAL_COUNTS_CONFIG = COUNT_MATRICES.get("viral_genome_counts", {})
VIRAL_BIN_SIZE = int(VIRAL_COUNTS_CONFIG.get("bin_size", 200))

RRNA_COUNTS_CONFIG = COUNT_MATRICES.get("rrna_counts", {})
RRNA_FLANK_SIZE = int(RRNA_COUNTS_CONFIG.get("flank_size", 250))

# Extract GTF file dictionaries from config
POL3_GTF_BY_TYPE = config.get("pol3_gtf", {})
REPEAT_ELEMENTS_GTF_BY_TYPE = config.get("repeat_elements_gtf", {})

# Pol3 types to process
POL3_TYPES = ["Pol3_T1", "Pol3_T2", "Pol3_T3"]
# Repeat element types to process
REPEAT_TYPES = [
    "SINE_Alu",
    "SINE_MIR",
    "LINE_L1",
    "LINE_L2",
    "LTR",
    "other_repeat_elements",
]

TN5_CALLERS = ("genrich", "macs2")

# Filter references by current host
POL3_TYPES_AVAILABLE = [
    t for t in POL3_TYPES if HOST in POL3_GTF_BY_TYPE.get(t, {})
]
REPEAT_TYPES_AVAILABLE = [
    t for t in REPEAT_TYPES if HOST in REPEAT_ELEMENTS_GTF_BY_TYPE.get(t, {})
]
VIRAL_TYPES_AVAILABLE = [v for v in VIRUS_LIST] if VIRUS_LIST else []


# Helper function to resolve GTF path
def _resolve_gtf_path(gtf_file: str) -> str:
    if os.path.isabs(gtf_file):
        return gtf_file
    return join(FASTAS_GTFS_DIR, gtf_file)


# Build output tracking lists for all matrix types

# Gene counts
GENE_BIN_OUTPUTS = []
GENE_COUNT_OUTPUTS = []
GENE_MATRIX_OUTPUTS = []
if COUNT_MATRIX_TYPES.get("gene", True) and HOST != "":
    GENE_BIN_OUTPUTS.append(
        join(REF_DIR, f"ref.tn5.host_gene_bins.{GENE_FLANK_SIZE}bp.bed")
    )
    for caller in TN5_CALLERS:
        GENE_COUNT_OUTPUTS.extend(
            expand(
                join(
                    RESULTSDIR,
                    "{sample}",
                    "tn5_motif",
                    caller,
                    "host",
                    "{sample}.host." + caller + f".tn5_gene_counts.{GENE_FLANK_SIZE}bp.tsv",
                ),
                sample=SAMPLES,
            )
        )
        GENE_MATRIX_OUTPUTS.append(
            join(RESULTSDIR, "tn5_motif", f"host.{caller}.tn5_gene_count_matrix.{GENE_FLANK_SIZE}bp.tsv")
        )

# tRNA counts
TRNA_BIN_OUTPUTS = []
TRNA_COUNT_OUTPUTS = []
TRNA_MATRIX_OUTPUTS = []
if COUNT_MATRIX_TYPES.get("trna", True) and HOST != "" and TRNAS_GTF != "":
    TRNA_BIN_OUTPUTS.append(
        join(REF_DIR, f"ref.tn5.host_trna_bins.{TRNA_FLANK_SIZE}bp.bed")
    )
    for caller in TN5_CALLERS:
        TRNA_COUNT_OUTPUTS.extend(
            expand(
                join(
                    RESULTSDIR,
                    "{sample}",
                    "tn5_motif",
                    caller,
                    "host",
                    "{sample}.host." + caller + f".tn5_trna_counts.{TRNA_FLANK_SIZE}bp.tsv",
                ),
                sample=SAMPLES,
            )
        )
        TRNA_MATRIX_OUTPUTS.append(
            join(RESULTSDIR, "tn5_motif", f"host.{caller}.tn5_trna_count_matrix.{TRNA_FLANK_SIZE}bp.tsv")
        )

# Pol3 counts (3 types: T1, T2, T3)
POL3_BIN_OUTPUTS = []
POL3_COUNT_OUTPUTS = []
POL3_MATRIX_OUTPUTS = []
if COUNT_MATRIX_TYPES.get("pol3", True) and HOST != "":
    for pol3_type in POL3_TYPES_AVAILABLE:
        POL3_BIN_OUTPUTS.append(
            join(REF_DIR, f"ref.tn5.host_pol3_{pol3_type.lower()}_bins.{POL3_FLANK_SIZE}bp.bed")
        )
        for caller in TN5_CALLERS:
            POL3_COUNT_OUTPUTS.extend(
                expand(
                    join(
                        RESULTSDIR,
                        "{sample}",
                        "tn5_motif",
                        caller,
                        "host",
                        "{sample}.host." + caller + f".tn5_pol3_{pol3_type.lower()}_counts.tsv",
                    ),
                    sample=SAMPLES,
                )
            )
            POL3_MATRIX_OUTPUTS.append(
                join(RESULTSDIR, "tn5_motif", f"host.{caller}.tn5_pol3_{pol3_type.lower()}_count_matrix.tsv")
            )

# Repeat element counts (6 types)
REPEAT_BIN_OUTPUTS = []
REPEAT_COUNT_OUTPUTS = []
REPEAT_MATRIX_OUTPUTS = []
if COUNT_MATRIX_TYPES.get("repeat_elements", True) and HOST != "":
    for repeat_type in REPEAT_TYPES_AVAILABLE:
        REPEAT_BIN_OUTPUTS.append(
            join(REF_DIR, f"ref.tn5.host_repeat_{repeat_type.lower()}_bins.{REPEAT_FLANK_SIZE}bp.bed")
        )
        for caller in TN5_CALLERS:
            REPEAT_COUNT_OUTPUTS.extend(
                expand(
                    join(
                        RESULTSDIR,
                        "{sample}",
                        "tn5_motif",
                        caller,
                        "host",
                        "{sample}.host." + caller + f".tn5_repeat_{repeat_type.lower()}_counts.tsv",
                    ),
                    sample=SAMPLES,
                )
            )
            REPEAT_MATRIX_OUTPUTS.append(
                join(RESULTSDIR, "tn5_motif", f"host.{caller}.tn5_repeat_{repeat_type.lower()}_count_matrix.tsv")
            )

# Viral genome counts (per virus)
VIRAL_BIN_OUTPUTS = []
VIRAL_COUNT_OUTPUTS = []
VIRAL_MATRIX_OUTPUTS = []
if COUNT_MATRIX_TYPES.get("viral", True) and len(VIRAL_TYPES_AVAILABLE) > 0:
    for virus in VIRAL_TYPES_AVAILABLE:
        VIRAL_BIN_OUTPUTS.append(
            join(REF_DIR, f"ref.tn5.virus_{virus}_bins.{VIRAL_BIN_SIZE}bp.bed")
        )
        for caller in TN5_CALLERS:
            VIRAL_COUNT_OUTPUTS.extend(
                expand(
                    join(
                        RESULTSDIR,
                        "{sample}",
                        "tn5_motif",
                        caller,
                        virus,
                        "{sample}." + virus + "." + caller + f".tn5_virus_counts.{VIRAL_BIN_SIZE}bp.tsv",
                    ),
                    sample=SAMPLES,
                )
            )
            VIRAL_MATRIX_OUTPUTS.append(
                join(RESULTSDIR, "tn5_motif", f"{virus}.{caller}.tn5_virus_count_matrix.{VIRAL_BIN_SIZE}bp.tsv")
            )

# rRNA counts
RRNA_BIN_OUTPUTS = []
RRNA_COUNT_OUTPUTS = []
RRNA_MATRIX_OUTPUTS = []
if COUNT_MATRIX_TYPES.get("rrna", True) and HOST != "" and CHRR_GTF != "":
    RRNA_BIN_OUTPUTS.append(
        join(REF_DIR, f"ref.tn5.host_rrna_bins.{RRNA_FLANK_SIZE}bp.bed")
    )
    for caller in TN5_CALLERS:
        RRNA_COUNT_OUTPUTS.extend(
            expand(
                join(
                    RESULTSDIR,
                    "{sample}",
                    "tn5_motif",
                    caller,
                    "host",
                    "{sample}.host." + caller + f".tn5_rrna_counts.{RRNA_FLANK_SIZE}bp.tsv",
                ),
                sample=SAMPLES,
            )
        )
        RRNA_MATRIX_OUTPUTS.append(
            join(RESULTSDIR, "tn5_motif", f"host.{caller}.tn5_rrna_count_matrix.{RRNA_FLANK_SIZE}bp.tsv")
        )

# Aggregate all count matrix outputs
COUNT_MATRIX_ALL_OUTPUTS = (
    GENE_BIN_OUTPUTS
    + GENE_COUNT_OUTPUTS
    + GENE_MATRIX_OUTPUTS
    + TRNA_BIN_OUTPUTS
    + TRNA_COUNT_OUTPUTS
    + TRNA_MATRIX_OUTPUTS
    + POL3_BIN_OUTPUTS
    + POL3_COUNT_OUTPUTS
    + POL3_MATRIX_OUTPUTS
    + REPEAT_BIN_OUTPUTS
    + REPEAT_COUNT_OUTPUTS
    + REPEAT_MATRIX_OUTPUTS
    + VIRAL_BIN_OUTPUTS
    + VIRAL_COUNT_OUTPUTS
    + VIRAL_MATRIX_OUTPUTS
    + RRNA_BIN_OUTPUTS
    + RRNA_COUNT_OUTPUTS
    + RRNA_MATRIX_OUTPUTS
)


# ====== BIN BUILDING RULES ======

if COUNT_MATRIX_TYPES.get("gene", True) and HOST != "":

    rule build_gene_bins:
        input:
            gtf=lambda wildcards: _resolve_gtf_path(
                REFERENCE_GTF_BY_HOST.get(HOST, HOST + ".gtf")
            ),
            host_regions=REF_REGIONS_HOST,
            chromsizes=join(REF_DIR, "ref.chrom.sizes.host.txt"),
        output:
            bed=join(REF_DIR, f"ref.tn5.host_gene_bins.{GENE_FLANK_SIZE}bp.bed"),
        params:
            flank_size=str(GENE_FLANK_SIZE),
            exclude_features=" ".join(GENE_EXCLUDE_FEATURES),
        threads:
            _get_threads("build_gene_bins", profile_config)
        container:
            config["containers"]["py311"]
        log:
            join(_logdir("build_gene_bins"), "gene.log")
        shell:
            r"""
            set -exo pipefail
            mkdir -p $(dirname {log})
            exec > >(tee -a {log}) 2>&1

            exclude_args=""
            if [ ! -z "{params.exclude_features}" ]; then
              exclude_args="--gene-types {params.exclude_features}"
            fi

            python {SCRIPTS_DIR}/build_tn5_tss_bins.py \
              --gtf {input.gtf} \
              --host-regions {input.host_regions} \
              --chromsizes {input.chromsizes} \
              --flank-size {params.flank_size} \
              $exclude_args \
              --output {output.bed}
            """


if COUNT_MATRIX_TYPES.get("trna", True) and HOST != "" and TRNAS_GTF != "":

    rule build_trna_bins:
        input:
            gtf=TRNAS_GTF,
            host_regions=REF_REGIONS_HOST,
            chromsizes=join(REF_DIR, "ref.chrom.sizes.host.txt"),
        output:
            bed=join(REF_DIR, f"ref.tn5.host_trna_bins.{TRNA_FLANK_SIZE}bp.bed"),
        params:
            flank_size=str(TRNA_FLANK_SIZE),
            max_size=str(TRNA_MAX_SIZE),
        threads:
            _get_threads("build_trna_bins", profile_config)
        container:
            config["containers"]["py311"]
        log:
            join(_logdir("build_trna_bins"), "trna.log")
        shell:
            r"""
            set -exo pipefail
            mkdir -p $(dirname {log})
            exec > >(tee -a {log}) 2>&1
            python {SCRIPTS_DIR}/build_tn5_midpoint_bins.py \
              --gtf {input.gtf} \
              --host-regions {input.host_regions} \
              --chromsizes {input.chromsizes} \
              --flank-size {params.flank_size} \
              --max-size {params.max_size} \
              --output {output.bed}
            """


if COUNT_MATRIX_TYPES.get("pol3", True) and HOST != "":

    rule build_pol3_bins:
        input:
            gtf=lambda wildcards: _resolve_gtf_path(
                POL3_GTF_BY_TYPE.get(wildcards.pol3_type, {}).get(HOST, "")
            ),
            host_regions=REF_REGIONS_HOST,
            chromsizes=join(REF_DIR, "ref.chrom.sizes.host.txt"),
        output:
            bed=join(
                REF_DIR,
                "ref.tn5.host_pol3_{pol3_type}_bins." + str(POL3_FLANK_SIZE) + "bp.bed",
            ),
        params:
            flank_size=str(POL3_FLANK_SIZE),
            max_size=str(POL3_MAX_SIZE),
        threads:
            _get_threads("build_pol3_bins", profile_config)
        container:
            config["containers"]["py311"]
        log:
            join(_logdir("build_pol3_bins"), "{pol3_type}.log")
        shell:
            r"""
            set -exo pipefail
            mkdir -p $(dirname {log})
            exec > >(tee -a {log}) 2>&1
            python {SCRIPTS_DIR}/build_tn5_midpoint_bins.py \
              --gtf {input.gtf} \
              --host-regions {input.host_regions} \
              --chromsizes {input.chromsizes} \
              --flank-size {params.flank_size} \
              --max-size {params.max_size} \
              --output {output.bed}
            """


if COUNT_MATRIX_TYPES.get("repeat_elements", True) and HOST != "":

    rule build_repeat_bins:
        input:
            gtf=lambda wildcards: _resolve_gtf_path(
                REPEAT_ELEMENTS_GTF_BY_TYPE.get(wildcards.repeat_type, {}).get(HOST, "")
            ),
            host_regions=REF_REGIONS_HOST,
            chromsizes=join(REF_DIR, "ref.chrom.sizes.host.txt"),
        output:
            bed=join(
                REF_DIR,
                "ref.tn5.host_repeat_{repeat_type}_bins." + str(REPEAT_FLANK_SIZE) + "bp.bed",
            ),
        params:
            flank_size=str(REPEAT_FLANK_SIZE),
            max_size=str(REPEAT_MAX_SIZE),
        threads:
            _get_threads("build_repeat_bins", profile_config)
        container:
            config["containers"]["py311"]
        log:
            join(_logdir("build_repeat_bins"), "{repeat_type}.log")
        shell:
            r"""
            set -exo pipefail
            mkdir -p $(dirname {log})
            exec > >(tee -a {log}) 2>&1
            python {SCRIPTS_DIR}/build_tn5_midpoint_bins.py \
              --gtf {input.gtf} \
              --host-regions {input.host_regions} \
              --chromsizes {input.chromsizes} \
              --flank-size {params.flank_size} \
              --max-size {params.max_size} \
              --output {output.bed}
            """


if COUNT_MATRIX_TYPES.get("viral", True) and len(VIRAL_TYPES_AVAILABLE) > 0:

    rule build_viral_bins:
        input:
            chromsizes=join(REF_DIR, "ref.chrom.sizes.{virus}.txt"),
        output:
            bed=join(
                REF_DIR, "ref.tn5.virus_{virus}_bins." + str(VIRAL_BIN_SIZE) + "bp.bed"
            ),
        params:
            bin_size=str(VIRAL_BIN_SIZE),
        threads:
            _get_threads("build_viral_bins", profile_config)
        container:
            config["containers"]["py311"]
        log:
            join(_logdir("build_viral_bins"), "{virus}.log")
        shell:
            r"""
            set -exo pipefail
            mkdir -p $(dirname {log})
            exec > >(tee -a {log}) 2>&1
            python {SCRIPTS_DIR}/build_viral_genome_bins.py \
              --chromsizes {input.chromsizes} \
              --bin-size {params.bin_size} \
              --output {output.bed}
            """


if COUNT_MATRIX_TYPES.get("rrna", True) and HOST != "" and CHRR_GTF != "":

    rule build_rrna_bins:
        input:
            gtf=CHRR_GTF,
            host_regions=REF_REGIONS_HOST,
            chromsizes=join(REF_DIR, "ref.chrom.sizes.host.txt"),
        output:
            bed=join(REF_DIR, f"ref.tn5.host_rrna_bins.{RRNA_FLANK_SIZE}bp.bed"),
        params:
            flank_size=str(RRNA_FLANK_SIZE),
        threads:
            _get_threads("build_rrna_bins", profile_config)
        container:
            config["containers"]["py311"]
        log:
            join(_logdir("build_rrna_bins"), "rrna.log")
        shell:
            r"""
            set -exo pipefail
            mkdir -p $(dirname {log})
            exec > >(tee -a {log}) 2>&1
            python {SCRIPTS_DIR}/build_tn5_tss_bins.py \
              --gtf {input.gtf} \
              --host-regions {input.host_regions} \
              --chromsizes {input.chromsizes} \
              --flank-size {params.flank_size} \
              --output {output.bed}
            """
