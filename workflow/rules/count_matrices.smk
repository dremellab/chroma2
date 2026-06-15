###################################################################################
# Tn5 count matrix generation from BAM files (Genrich only)
# Generates 6 types of matrices: gene, tRNA, Pol3, repeat elements, viral, rRNA
# Output folder: results/count_matrices
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

# Output folder for count matrices
COUNT_MATRICES_DIR = join(RESULTSDIR, "count_matrices")

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


# Get Genrich BAM for a sample
def _get_genrich_bam(sample: str, group: str = "host") -> str:
    if group == "host":
        return join(RESULTSDIR, sample, "align", f"{sample}.aligned.host.qname.bam")
    return join(RESULTSDIR, sample, "align", f"{sample}.aligned.virus.{group}.qname.bam")


# Build output tracking lists for all matrix types

# Gene counts
GENE_BIN_OUTPUTS = []
GENE_COUNT_OUTPUTS = []
GENE_MATRIX_OUTPUTS = []
if COUNT_MATRIX_TYPES.get("gene", True) and HOST != "":
    GENE_BIN_OUTPUTS.append(
        join(REF_DIR, f"ref.tn5.host_gene_bins.{GENE_FLANK_SIZE}bp.bed")
    )
    GENE_COUNT_OUTPUTS.extend(
        expand(
            join(
                COUNT_MATRICES_DIR,
                "gene",
                "{sample}.gene_counts.tsv",
            ),
            sample=SAMPLES,
        )
    )
    GENE_MATRIX_OUTPUTS.append(
        join(COUNT_MATRICES_DIR, "gene", "gene_count_matrix.tsv")
    )

# tRNA counts
TRNA_BIN_OUTPUTS = []
TRNA_COUNT_OUTPUTS = []
TRNA_MATRIX_OUTPUTS = []
if COUNT_MATRIX_TYPES.get("trna", True) and HOST != "" and TRNAS_GTF != "":
    TRNA_BIN_OUTPUTS.append(
        join(REF_DIR, f"ref.tn5.host_trna_bins.{TRNA_FLANK_SIZE}bp.bed")
    )
    TRNA_COUNT_OUTPUTS.extend(
        expand(
            join(
                COUNT_MATRICES_DIR,
                "trna",
                "{sample}.trna_counts.tsv",
            ),
            sample=SAMPLES,
        )
    )
    TRNA_MATRIX_OUTPUTS.append(
        join(COUNT_MATRICES_DIR, "trna", "trna_count_matrix.tsv")
    )

# Pol3 counts (3 types: T1, T2, T3)
POL3_BIN_OUTPUTS = []
POL3_COUNT_OUTPUTS = []
POL3_MATRIX_OUTPUTS = []
if COUNT_MATRIX_TYPES.get("pol3", True) and HOST != "":
    for pol3_type in POL3_TYPES_AVAILABLE:
        pol3_lower = pol3_type.lower()
        POL3_BIN_OUTPUTS.append(
            join(REF_DIR, f"ref.tn5.host_pol3_{pol3_lower}_bins.{POL3_FLANK_SIZE}bp.bed")
        )
        POL3_COUNT_OUTPUTS.extend(
            expand(
                join(
                    COUNT_MATRICES_DIR,
                    f"pol3_{pol3_lower}",
                    "{sample}.pol3_{pol3_lower}_counts.tsv",
                ),
                sample=SAMPLES,
            )
        )
        POL3_MATRIX_OUTPUTS.append(
            join(COUNT_MATRICES_DIR, f"pol3_{pol3_lower}", f"pol3_{pol3_lower}_count_matrix.tsv")
        )

# Repeat element counts (6 types)
REPEAT_BIN_OUTPUTS = []
REPEAT_COUNT_OUTPUTS = []
REPEAT_MATRIX_OUTPUTS = []
if COUNT_MATRIX_TYPES.get("repeat_elements", True) and HOST != "":
    for repeat_type in REPEAT_TYPES_AVAILABLE:
        repeat_lower = repeat_type.lower()
        REPEAT_BIN_OUTPUTS.append(
            join(REF_DIR, f"ref.tn5.host_repeat_{repeat_lower}_bins.{REPEAT_FLANK_SIZE}bp.bed")
        )
        REPEAT_COUNT_OUTPUTS.extend(
            expand(
                join(
                    COUNT_MATRICES_DIR,
                    f"repeat_{repeat_lower}",
                    "{sample}.repeat_{repeat_lower}_counts.tsv",
                ),
                sample=SAMPLES,
            )
        )
        REPEAT_MATRIX_OUTPUTS.append(
            join(COUNT_MATRICES_DIR, f"repeat_{repeat_lower}", f"repeat_{repeat_lower}_count_matrix.tsv")
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
        VIRAL_COUNT_OUTPUTS.extend(
            expand(
                join(
                    COUNT_MATRICES_DIR,
                    f"viral_{virus}",
                    "{sample}.viral_{virus}_counts.tsv",
                ),
                sample=SAMPLES,
            )
        )
        VIRAL_MATRIX_OUTPUTS.append(
            join(COUNT_MATRICES_DIR, f"viral_{virus}", f"viral_{virus}_count_matrix.tsv")
        )

# rRNA counts
RRNA_BIN_OUTPUTS = []
RRNA_COUNT_OUTPUTS = []
RRNA_MATRIX_OUTPUTS = []
if COUNT_MATRIX_TYPES.get("rrna", True) and HOST != "" and CHRR_GTF != "":
    RRNA_BIN_OUTPUTS.append(
        join(REF_DIR, f"ref.tn5.host_rrna_bins.{RRNA_FLANK_SIZE}bp.bed")
    )
    RRNA_COUNT_OUTPUTS.extend(
        expand(
            join(
                COUNT_MATRICES_DIR,
                "rrna",
                "{sample}.rrna_counts.tsv",
            ),
            sample=SAMPLES,
        )
    )
    RRNA_MATRIX_OUTPUTS.append(
        join(COUNT_MATRICES_DIR, "rrna", "rrna_count_matrix.tsv")
    )

# Final consolidated output folder
COUNT_MATRICES_FINAL_DIR = join(COUNT_MATRICES_DIR, "final_matrices")

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

# Add final consolidated outputs
COUNT_MATRIX_ALL_OUTPUTS.append(join(COUNT_MATRICES_FINAL_DIR, "COUNT_MATRICES_INDEX.txt"))
if COUNT_MATRIX_TYPES.get("gene", True) and HOST != "":
    COUNT_MATRIX_ALL_OUTPUTS.append(join(COUNT_MATRICES_FINAL_DIR, "gene_count_matrix.tsv"))
if COUNT_MATRIX_TYPES.get("trna", True) and HOST != "" and TRNAS_GTF != "":
    COUNT_MATRIX_ALL_OUTPUTS.append(join(COUNT_MATRICES_FINAL_DIR, "trna_count_matrix.tsv"))
if COUNT_MATRIX_TYPES.get("pol3", True) and HOST != "":
    for pol3_type in POL3_TYPES_AVAILABLE:
        COUNT_MATRIX_ALL_OUTPUTS.append(
            join(COUNT_MATRICES_FINAL_DIR, f"pol3_{pol3_type.lower()}_count_matrix.tsv")
        )
if COUNT_MATRIX_TYPES.get("repeat_elements", True) and HOST != "":
    for repeat_type in REPEAT_TYPES_AVAILABLE:
        COUNT_MATRIX_ALL_OUTPUTS.append(
            join(COUNT_MATRICES_FINAL_DIR, f"repeat_{repeat_type.lower()}_count_matrix.tsv")
        )
if COUNT_MATRIX_TYPES.get("viral", True) and len(VIRAL_TYPES_AVAILABLE) > 0:
    for virus in VIRAL_TYPES_AVAILABLE:
        COUNT_MATRIX_ALL_OUTPUTS.append(
            join(COUNT_MATRICES_FINAL_DIR, f"viral_{virus}_count_matrix.tsv")
        )
if COUNT_MATRIX_TYPES.get("rrna", True) and HOST != "" and CHRR_GTF != "":
    COUNT_MATRIX_ALL_OUTPUTS.append(join(COUNT_MATRICES_FINAL_DIR, "rrna_count_matrix.tsv"))


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


# ====== COUNTING RULES (Genrich only) ======

TN5_MAPQ_MIN = int(config.get("tn5_motif", {}).get("mapq_min", 0))
TN5_EXCLUDE_SECONDARY = bool(config.get("tn5_motif", {}).get("exclude_secondary", False))
TN5_EXCLUDE_SUPPLEMENTARY = bool(
    config.get("tn5_motif", {}).get("exclude_supplementary", False)
)
TN5_FRACTIONAL_COUNTING = bool(config.get("tn5_motif", {}).get("fractional_counting", False))


if COUNT_MATRIX_TYPES.get("gene", True) and HOST != "":

    rule count_gene_bins:
        input:
            bam=lambda wildcards: _get_genrich_bam(wildcards.sample, "host"),
            bins=join(REF_DIR, f"ref.tn5.host_gene_bins.{GENE_FLANK_SIZE}bp.bed"),
        output:
            tsv=join(COUNT_MATRICES_DIR, "gene", "{sample}.gene_counts.tsv"),
        params:
            mapq_min=str(TN5_MAPQ_MIN),
            filter_args=" ".join(
                flag
                for flag in (
                    "--exclude-secondary" if TN5_EXCLUDE_SECONDARY else "",
                    "--exclude-supplementary" if TN5_EXCLUDE_SUPPLEMENTARY else "",
                    "--fractional-counting" if TN5_FRACTIONAL_COUNTING else "",
                )
                if flag
            ),
        threads:
            _get_threads("count_gene_bins", profile_config)
        container:
            config["containers"]["pysam"]
        log:
            join(_logdir("count_gene_bins"), "{sample}.log")
        shell:
            r"""
            set -exo pipefail
            mkdir -p $(dirname {log}) $(dirname {output.tsv})
            exec > >(tee -a {log}) 2>&1
            python {SCRIPTS_DIR}/count_tn5_sites_in_bins.py \
              --bam {input.bam} \
              --bins-bed {input.bins} \
              --output {output.tsv} \
              --sample {wildcards.sample} \
              --threads {threads} \
              --mapq-min {params.mapq_min} {params.filter_args}
            """


if COUNT_MATRIX_TYPES.get("trna", True) and HOST != "" and TRNAS_GTF != "":

    rule count_trna_bins:
        input:
            bam=lambda wildcards: _get_genrich_bam(wildcards.sample, "host"),
            bins=join(REF_DIR, f"ref.tn5.host_trna_bins.{TRNA_FLANK_SIZE}bp.bed"),
        output:
            tsv=join(COUNT_MATRICES_DIR, "trna", "{sample}.trna_counts.tsv"),
        params:
            mapq_min=str(TN5_MAPQ_MIN),
            filter_args=" ".join(
                flag
                for flag in (
                    "--exclude-secondary" if TN5_EXCLUDE_SECONDARY else "",
                    "--exclude-supplementary" if TN5_EXCLUDE_SUPPLEMENTARY else "",
                    "--fractional-counting" if TN5_FRACTIONAL_COUNTING else "",
                )
                if flag
            ),
        threads:
            _get_threads("count_trna_bins", profile_config)
        container:
            config["containers"]["pysam"]
        log:
            join(_logdir("count_trna_bins"), "{sample}.log")
        shell:
            r"""
            set -exo pipefail
            mkdir -p $(dirname {log}) $(dirname {output.tsv})
            exec > >(tee -a {log}) 2>&1
            python {SCRIPTS_DIR}/count_tn5_sites_in_bins.py \
              --bam {input.bam} \
              --bins-bed {input.bins} \
              --output {output.tsv} \
              --sample {wildcards.sample} \
              --threads {threads} \
              --mapq-min {params.mapq_min} {params.filter_args}
            """


if COUNT_MATRIX_TYPES.get("pol3", True) and HOST != "":

    rule count_pol3_bins:
        input:
            bam=lambda wildcards: _get_genrich_bam(wildcards.sample, "host"),
            bins=join(
                REF_DIR,
                "ref.tn5.host_pol3_{pol3_type}_bins." + str(POL3_FLANK_SIZE) + "bp.bed",
            ),
        output:
            tsv=join(COUNT_MATRICES_DIR, "pol3_{pol3_type}", "{sample}.pol3_{pol3_type}_counts.tsv"),
        params:
            mapq_min=str(TN5_MAPQ_MIN),
            filter_args=" ".join(
                flag
                for flag in (
                    "--exclude-secondary" if TN5_EXCLUDE_SECONDARY else "",
                    "--exclude-supplementary" if TN5_EXCLUDE_SUPPLEMENTARY else "",
                    "--fractional-counting" if TN5_FRACTIONAL_COUNTING else "",
                )
                if flag
            ),
        threads:
            _get_threads("count_pol3_bins", profile_config)
        container:
            config["containers"]["pysam"]
        log:
            join(_logdir("count_pol3_bins"), "{sample}.{pol3_type}.log")
        shell:
            r"""
            set -exo pipefail
            mkdir -p $(dirname {log}) $(dirname {output.tsv})
            exec > >(tee -a {log}) 2>&1
            python {SCRIPTS_DIR}/count_tn5_sites_in_bins.py \
              --bam {input.bam} \
              --bins-bed {input.bins} \
              --output {output.tsv} \
              --sample {wildcards.sample} \
              --threads {threads} \
              --mapq-min {params.mapq_min} {params.filter_args}
            """


if COUNT_MATRIX_TYPES.get("repeat_elements", True) and HOST != "":

    rule count_repeat_bins:
        input:
            bam=lambda wildcards: _get_genrich_bam(wildcards.sample, "host"),
            bins=join(
                REF_DIR,
                "ref.tn5.host_repeat_{repeat_type}_bins." + str(REPEAT_FLANK_SIZE) + "bp.bed",
            ),
        output:
            tsv=join(COUNT_MATRICES_DIR, "repeat_{repeat_type}", "{sample}.repeat_{repeat_type}_counts.tsv"),
        params:
            mapq_min=str(TN5_MAPQ_MIN),
            filter_args=" ".join(
                flag
                for flag in (
                    "--exclude-secondary" if TN5_EXCLUDE_SECONDARY else "",
                    "--exclude-supplementary" if TN5_EXCLUDE_SUPPLEMENTARY else "",
                    "--fractional-counting" if TN5_FRACTIONAL_COUNTING else "",
                )
                if flag
            ),
        threads:
            _get_threads("count_repeat_bins", profile_config)
        container:
            config["containers"]["pysam"]
        log:
            join(_logdir("count_repeat_bins"), "{sample}.{repeat_type}.log")
        shell:
            r"""
            set -exo pipefail
            mkdir -p $(dirname {log}) $(dirname {output.tsv})
            exec > >(tee -a {log}) 2>&1
            python {SCRIPTS_DIR}/count_tn5_sites_in_bins.py \
              --bam {input.bam} \
              --bins-bed {input.bins} \
              --output {output.tsv} \
              --sample {wildcards.sample} \
              --threads {threads} \
              --mapq-min {params.mapq_min} {params.filter_args}
            """


if COUNT_MATRIX_TYPES.get("viral", True) and len(VIRAL_TYPES_AVAILABLE) > 0:

    rule count_viral_bins:
        input:
            bam=lambda wildcards: _get_genrich_bam(wildcards.sample, wildcards.virus),
            bins=join(REF_DIR, "ref.tn5.virus_{virus}_bins." + str(VIRAL_BIN_SIZE) + "bp.bed"),
        output:
            tsv=join(COUNT_MATRICES_DIR, "viral_{virus}", "{sample}.viral_{virus}_counts.tsv"),
        params:
            mapq_min=str(TN5_MAPQ_MIN),
            filter_args=" ".join(
                flag
                for flag in (
                    "--exclude-secondary" if TN5_EXCLUDE_SECONDARY else "",
                    "--exclude-supplementary" if TN5_EXCLUDE_SUPPLEMENTARY else "",
                    "--fractional-counting" if TN5_FRACTIONAL_COUNTING else "",
                )
                if flag
            ),
        threads:
            _get_threads("count_viral_bins", profile_config)
        container:
            config["containers"]["pysam"]
        log:
            join(_logdir("count_viral_bins"), "{sample}.{virus}.log")
        shell:
            r"""
            set -exo pipefail
            mkdir -p $(dirname {log}) $(dirname {output.tsv})
            exec > >(tee -a {log}) 2>&1
            python {SCRIPTS_DIR}/count_tn5_sites_in_bins.py \
              --bam {input.bam} \
              --bins-bed {input.bins} \
              --output {output.tsv} \
              --sample {wildcards.sample} \
              --threads {threads} \
              --mapq-min {params.mapq_min} {params.filter_args}
            """


if COUNT_MATRIX_TYPES.get("rrna", True) and HOST != "" and CHRR_GTF != "":

    rule count_rrna_bins:
        input:
            bam=lambda wildcards: _get_genrich_bam(wildcards.sample, "host"),
            bins=join(REF_DIR, f"ref.tn5.host_rrna_bins.{RRNA_FLANK_SIZE}bp.bed"),
        output:
            tsv=join(COUNT_MATRICES_DIR, "rrna", "{sample}.rrna_counts.tsv"),
        params:
            mapq_min=str(TN5_MAPQ_MIN),
            filter_args=" ".join(
                flag
                for flag in (
                    "--exclude-secondary" if TN5_EXCLUDE_SECONDARY else "",
                    "--exclude-supplementary" if TN5_EXCLUDE_SUPPLEMENTARY else "",
                    "--fractional-counting" if TN5_FRACTIONAL_COUNTING else "",
                )
                if flag
            ),
        threads:
            _get_threads("count_rrna_bins", profile_config)
        container:
            config["containers"]["pysam"]
        log:
            join(_logdir("count_rrna_bins"), "{sample}.log")
        shell:
            r"""
            set -exo pipefail
            mkdir -p $(dirname {log}) $(dirname {output.tsv})
            exec > >(tee -a {log}) 2>&1
            python {SCRIPTS_DIR}/count_tn5_sites_in_bins.py \
              --bam {input.bam} \
              --bins-bed {input.bins} \
              --output {output.tsv} \
              --sample {wildcards.sample} \
              --threads {threads} \
              --mapq-min {params.mapq_min} {params.filter_args}
            """


# ====== AGGREGATION RULES ======

if COUNT_MATRIX_TYPES.get("gene", True) and HOST != "":

    rule aggregate_gene_count_matrix:
        input:
            expand(join(COUNT_MATRICES_DIR, "gene", "{sample}.gene_counts.tsv"), sample=SAMPLES),
        output:
            join(COUNT_MATRICES_DIR, "gene", "gene_count_matrix.tsv"),
        threads:
            _get_threads("aggregate_gene_count_matrix", profile_config)
        container:
            config["containers"]["py311"]
        log:
            join(_logdir("aggregate_gene_count_matrix"), "aggregate.log")
        shell:
            r"""
            set -exo pipefail
            mkdir -p $(dirname {log})
            exec > >(tee -a {log}) 2>&1
            python {SCRIPTS_DIR}/build_tn5_count_matrix.py \
              --input-files {input} \
              --output {output}
            """


if COUNT_MATRIX_TYPES.get("trna", True) and HOST != "" and TRNAS_GTF != "":

    rule aggregate_trna_count_matrix:
        input:
            expand(join(COUNT_MATRICES_DIR, "trna", "{sample}.trna_counts.tsv"), sample=SAMPLES),
        output:
            join(COUNT_MATRICES_DIR, "trna", "trna_count_matrix.tsv"),
        threads:
            _get_threads("aggregate_trna_count_matrix", profile_config)
        container:
            config["containers"]["py311"]
        log:
            join(_logdir("aggregate_trna_count_matrix"), "aggregate.log")
        shell:
            r"""
            set -exo pipefail
            mkdir -p $(dirname {log})
            exec > >(tee -a {log}) 2>&1
            python {SCRIPTS_DIR}/build_tn5_count_matrix.py \
              --input-files {input} \
              --output {output}
            """


if COUNT_MATRIX_TYPES.get("pol3", True) and HOST != "":

    rule aggregate_pol3_count_matrix:
        input:
            expand(
                join(COUNT_MATRICES_DIR, "pol3_{{pol3_type}}", "{{sample}}.pol3_{{pol3_type}}_counts.tsv"),
                sample=SAMPLES,
            ),
        output:
            join(COUNT_MATRICES_DIR, "pol3_{pol3_type}", "pol3_{pol3_type}_count_matrix.tsv"),
        threads:
            _get_threads("aggregate_pol3_count_matrix", profile_config)
        container:
            config["containers"]["py311"]
        log:
            join(_logdir("aggregate_pol3_count_matrix"), "{pol3_type}.log")
        shell:
            r"""
            set -exo pipefail
            mkdir -p $(dirname {log})
            exec > >(tee -a {log}) 2>&1
            python {SCRIPTS_DIR}/build_tn5_count_matrix.py \
              --input-files {input} \
              --output {output}
            """


if COUNT_MATRIX_TYPES.get("repeat_elements", True) and HOST != "":

    rule aggregate_repeat_count_matrix:
        input:
            expand(
                join(COUNT_MATRICES_DIR, "repeat_{{repeat_type}}", "{{sample}}.repeat_{{repeat_type}}_counts.tsv"),
                sample=SAMPLES,
            ),
        output:
            join(COUNT_MATRICES_DIR, "repeat_{repeat_type}", "repeat_{repeat_type}_count_matrix.tsv"),
        threads:
            _get_threads("aggregate_repeat_count_matrix", profile_config)
        container:
            config["containers"]["py311"]
        log:
            join(_logdir("aggregate_repeat_count_matrix"), "{repeat_type}.log")
        shell:
            r"""
            set -exo pipefail
            mkdir -p $(dirname {log})
            exec > >(tee -a {log}) 2>&1
            python {SCRIPTS_DIR}/build_tn5_count_matrix.py \
              --input-files {input} \
              --output {output}
            """


if COUNT_MATRIX_TYPES.get("viral", True) and len(VIRAL_TYPES_AVAILABLE) > 0:

    rule aggregate_viral_count_matrix:
        input:
            expand(
                join(COUNT_MATRICES_DIR, "viral_{{virus}}", "{{sample}}.viral_{{virus}}_counts.tsv"),
                sample=SAMPLES,
            ),
        output:
            join(COUNT_MATRICES_DIR, "viral_{virus}", "viral_{virus}_count_matrix.tsv"),
        threads:
            _get_threads("aggregate_viral_count_matrix", profile_config)
        container:
            config["containers"]["py311"]
        log:
            join(_logdir("aggregate_viral_count_matrix"), "{virus}.log")
        shell:
            r"""
            set -exo pipefail
            mkdir -p $(dirname {log})
            exec > >(tee -a {log}) 2>&1
            python {SCRIPTS_DIR}/build_tn5_count_matrix.py \
              --input-files {input} \
              --output {output}
            """


if COUNT_MATRIX_TYPES.get("rrna", True) and HOST != "" and CHRR_GTF != "":

    rule aggregate_rrna_count_matrix:
        input:
            expand(join(COUNT_MATRICES_DIR, "rrna", "{sample}.rrna_counts.tsv"), sample=SAMPLES),
        output:
            join(COUNT_MATRICES_DIR, "rrna", "rrna_count_matrix.tsv"),
        threads:
            _get_threads("aggregate_rrna_count_matrix", profile_config)
        container:
            config["containers"]["py311"]
        log:
            join(_logdir("aggregate_rrna_count_matrix"), "aggregate.log")
        shell:
            r"""
            set -exo pipefail
            mkdir -p $(dirname {log})
            exec > >(tee -a {log}) 2>&1
            python {SCRIPTS_DIR}/build_tn5_count_matrix.py \
              --input-files {input} \
              --output {output}
            """


# ====== FINAL CONSOLIDATION & DOCUMENTATION ======


rule consolidate_count_matrices:
    """Consolidate all count matrices to final_matrices folder with index."""
    input:
        gene=join(COUNT_MATRICES_DIR, "gene", "gene_count_matrix.tsv")
        if COUNT_MATRIX_TYPES.get("gene", True) and HOST != ""
        else [],
        trna=join(COUNT_MATRICES_DIR, "trna", "trna_count_matrix.tsv")
        if COUNT_MATRIX_TYPES.get("trna", True) and HOST != "" and TRNAS_GTF != ""
        else [],
        pol3=[
            join(COUNT_MATRICES_DIR, f"pol3_{t.lower()}", f"pol3_{t.lower()}_count_matrix.tsv")
            for t in POL3_TYPES_AVAILABLE
        ]
        if COUNT_MATRIX_TYPES.get("pol3", True) and HOST != ""
        else [],
        repeat=[
            join(COUNT_MATRICES_DIR, f"repeat_{r.lower()}", f"repeat_{r.lower()}_count_matrix.tsv")
            for r in REPEAT_TYPES_AVAILABLE
        ]
        if COUNT_MATRIX_TYPES.get("repeat_elements", True) and HOST != ""
        else [],
        viral=[
            join(COUNT_MATRICES_DIR, f"viral_{v}", f"viral_{v}_count_matrix.tsv")
            for v in VIRAL_TYPES_AVAILABLE
        ]
        if COUNT_MATRIX_TYPES.get("viral", True) and len(VIRAL_TYPES_AVAILABLE) > 0
        else [],
        rrna=join(COUNT_MATRICES_DIR, "rrna", "rrna_count_matrix.tsv")
        if COUNT_MATRIX_TYPES.get("rrna", True) and HOST != "" and CHRR_GTF != ""
        else [],
    output:
        index=join(COUNT_MATRICES_FINAL_DIR, "COUNT_MATRICES_INDEX.txt"),
        gene=join(COUNT_MATRICES_FINAL_DIR, "gene_count_matrix.tsv")
        if COUNT_MATRIX_TYPES.get("gene", True) and HOST != ""
        else [],
        trna=join(COUNT_MATRICES_FINAL_DIR, "trna_count_matrix.tsv")
        if COUNT_MATRIX_TYPES.get("trna", True) and HOST != "" and TRNAS_GTF != ""
        else [],
        pol3=[
            join(COUNT_MATRICES_FINAL_DIR, f"pol3_{t.lower()}_count_matrix.tsv")
            for t in POL3_TYPES_AVAILABLE
        ]
        if COUNT_MATRIX_TYPES.get("pol3", True) and HOST != ""
        else [],
        repeat=[
            join(COUNT_MATRICES_FINAL_DIR, f"repeat_{r.lower()}_count_matrix.tsv")
            for r in REPEAT_TYPES_AVAILABLE
        ]
        if COUNT_MATRIX_TYPES.get("repeat_elements", True) and HOST != ""
        else [],
        viral=[
            join(COUNT_MATRICES_FINAL_DIR, f"viral_{v}_count_matrix.tsv")
            for v in VIRAL_TYPES_AVAILABLE
        ]
        if COUNT_MATRIX_TYPES.get("viral", True) and len(VIRAL_TYPES_AVAILABLE) > 0
        else [],
        rrna=join(COUNT_MATRICES_FINAL_DIR, "rrna_count_matrix.tsv")
        if COUNT_MATRIX_TYPES.get("rrna", True) and HOST != "" and CHRR_GTF != ""
        else [],
    threads:
        _get_threads("consolidate_count_matrices", profile_config)
    container:
        config["containers"]["py311"]
    log:
        join(_logdir("consolidate_count_matrices"), "consolidate.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {output.index}) $(dirname {log})
        exec > >(tee -a {log}) 2>&1

        outdir=$(dirname {output.index})

        # Create symlinks for all matrices
        [ -f "{input.gene}" ] && ln -sf {input.gene} $outdir/gene_count_matrix.tsv || true
        [ -f "{input.trna}" ] && ln -sf {input.trna} $outdir/trna_count_matrix.tsv || true
        [ -f "{input.rrna}" ] && ln -sf {input.rrna} $outdir/rrna_count_matrix.tsv || true

        for f in {input.pol3}; do
          [ -f "$f" ] && ln -sf "$f" "$outdir/$(basename $f)" || true
        done

        for f in {input.repeat}; do
          [ -f "$f" ] && ln -sf "$f" "$outdir/$(basename $f)" || true
        done

        for f in {input.viral}; do
          [ -f "$f" ] && ln -sf "$f" "$outdir/$(basename $f)" || true
        done

        # Generate index file
        python3 << 'PYTHON_EOF'
import os
from pathlib import Path

outdir = "{output.index}"
outdir = os.path.dirname(outdir)

# Collect all symlinked matrices
matrices = {{}}
for f in Path(outdir).glob("*.tsv"):
    matrices[f.name] = f

# Write index
with open(os.path.join(outdir, "COUNT_MATRICES_INDEX.txt"), "w") as idx:
    idx.write("=" * 80 + "\n")
    idx.write("Tn5 Count Matrices Index\n")
    idx.write("=" * 80 + "\n")
    idx.write(f"Host: {HOST}\n")
    idx.write(f"Samples: {len(SAMPLES)}\n")
    idx.write(f"Output Folder: {outdir}\n\n")
    idx.write("=" * 80 + "\n")
    idx.write("AVAILABLE MATRICES\n")
    idx.write("=" * 80 + "\n\n")

    for name in sorted(matrices.keys()):
        idx.write(f"✓ {name}\n")

    idx.write("\n" + "=" * 80 + "\n")
    idx.write("CONFIGURATION\n")
    idx.write("=" * 80 + "\n")
    idx.write("Genrich BAM files used for all counting\n")
    idx.write("Fractional counting (NH-weighted): ENABLED\n")
    idx.write(f"Mapq minimum: {TN5_MAPQ_MIN}\n\n")
    idx.write("For detailed documentation, see: COUNT_MATRICES_README.md\n")
    idx.write("=" * 80 + "\n")

PYTHON_EOF

        echo "✓ Consolidated count matrices to $outdir"
        """
