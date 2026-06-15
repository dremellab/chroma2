###################################################################################
# DESeq2 differential expression analysis on count matrices
# Supports flexible analysis of any enabled count matrix type (gene, tRNA, Pol3, etc)
###################################################################################

import shlex

# Count matrices configuration
COUNT_MATRICES_CONFIG = config.get("count_matrices", {})
COUNT_MATRICES_FINAL_DIR = join(RESULTSDIR, "count_matrices", "final_matrices")

# Build available matrices dictionary based on config
DESEQ2_AVAILABLE_MATRICES = {}

# Gene counts
if COUNT_MATRICES_CONFIG.get("gene_counts", {}).get("enabled", True) and HOST != "":
    DESEQ2_AVAILABLE_MATRICES["gene"] = join(COUNT_MATRICES_FINAL_DIR, "gene_count_matrix.tsv")

# tRNA counts
if COUNT_MATRICES_CONFIG.get("trna_counts", {}).get("enabled", True) and HOST != "" and TRNAS_GTF != "":
    DESEQ2_AVAILABLE_MATRICES["trna"] = join(COUNT_MATRICES_FINAL_DIR, "trna_count_matrix.tsv")

# Pol3 counts (Pol3_T1, Pol3_T2, Pol3_T3)
POL3_GTF_BY_TYPE = config.get("pol3_gtf", {})
POL3_TYPES = ["Pol3_T1", "Pol3_T2", "Pol3_T3"]
if COUNT_MATRICES_CONFIG.get("pol3_counts", {}).get("enabled", True) and HOST != "":
    for pol3_type in POL3_TYPES:
        if HOST in POL3_GTF_BY_TYPE.get(pol3_type, {}):
            pol3_lower = pol3_type.lower()
            DESEQ2_AVAILABLE_MATRICES[f"pol3_{pol3_lower}"] = join(
                COUNT_MATRICES_FINAL_DIR, f"pol3_{pol3_lower}_count_matrix.tsv"
            )

# Repeat element counts (6 types)
REPEAT_ELEMENTS_GTF_BY_TYPE = config.get("repeat_elements_gtf", {})
REPEAT_TYPES = [
    "SINE_Alu",
    "SINE_MIR",
    "LINE_L1",
    "LINE_L2",
    "LTR",
    "other_repeat_elements",
]
if COUNT_MATRICES_CONFIG.get("repeat_element_counts", {}).get("enabled", True) and HOST != "":
    for repeat_type in REPEAT_TYPES:
        if HOST in REPEAT_ELEMENTS_GTF_BY_TYPE.get(repeat_type, {}):
            repeat_lower = repeat_type.lower()
            DESEQ2_AVAILABLE_MATRICES[f"repeat_{repeat_lower}"] = join(
                COUNT_MATRICES_FINAL_DIR, f"repeat_{repeat_lower}_count_matrix.tsv"
            )

# Viral counts (per virus)
if COUNT_MATRICES_CONFIG.get("viral_genome_counts", {}).get("enabled", True) and len(VIRUS_LIST) > 0:
    for virus in VIRUS_LIST:
        DESEQ2_AVAILABLE_MATRICES[f"viral_{virus}"] = join(
            COUNT_MATRICES_FINAL_DIR, f"viral_{virus}_count_matrix.tsv"
        )

# rRNA counts
if COUNT_MATRICES_CONFIG.get("rrna_counts", {}).get("enabled", True) and HOST != "" and CHRR_GTF != "":
    DESEQ2_AVAILABLE_MATRICES["rrna"] = join(COUNT_MATRICES_FINAL_DIR, "rrna_count_matrix.tsv")

DESEQ2_SKIP_FEATURES = DESEQ2_CONFIG.get("skip_features", [])


# Build DESeq2 outputs if enabled
DESEQ2_OUTPUTS = []
if DESEQ2_ENABLED and DESEQ2_CONTRASTS:
    for contrast in DESEQ2_CONTRASTS:
        comparison = contrast["comparison"]
        # Generate outputs for each available matrix
        for matrix_type in DESEQ2_AVAILABLE_MATRICES.keys():
            DESEQ2_OUTPUTS.append(
                join(DESEQ2_OUTDIR, f"{matrix_type}_{comparison}", f"{matrix_type}_{comparison}.deseq2_results.tsv")
            )
            DESEQ2_OUTPUTS.append(
                join(DESEQ2_OUTDIR, f"{matrix_type}_{comparison}", f"{matrix_type}_{comparison}.volcano.png")
            )
        DESEQ2_OUTPUTS.append(join(DESEQ2_OUTDIR, comparison, f"{comparison}.deseq2_summary.html"))


if DESEQ2_ENABLED:

    rule deseq2_analysis:
        input:
            count_matrix=lambda wc: DESEQ2_AVAILABLE_MATRICES.get(wc.matrix_type, ""),
            manifest=MANIFEST_FILE,
            contrasts_file=DESEQ2_CONTRASTS_FILE,
        output:
            results_tsv=join(
                DESEQ2_OUTDIR,
                "{matrix_type}_{comparison}",
                "{matrix_type}_{comparison}.deseq2_results.tsv",
            ),
            volcano_png=join(
                DESEQ2_OUTDIR,
                "{matrix_type}_{comparison}",
                "{matrix_type}_{comparison}.volcano.png",
            ),
        params:
            group1=lambda wc: shlex.quote(_deseq2_contrast_group(wc.comparison, "group1")),
            group2=lambda wc: shlex.quote(_deseq2_contrast_group(wc.comparison, "group2")),
            matrix_type="{matrix_type}",
            comparison="{comparison}",
            alpha=str(DESEQ2_ALPHA),
            lfc_threshold=str(DESEQ2_LFC_THRESHOLD),
            min_replicates=str(DESEQ2_MIN_REPLICATES_PER_GROUP),
            size_factor_type=str(DESEQ2_SIZE_FACTOR_TYPE),
            fit_type=str(DESEQ2_FIT_TYPE),
            p_adjust_method=str(DESEQ2_P_ADJUST_METHOD),
        threads:
            _get_threads("deseq2_analysis", profile_config)
        container:
            config["containers"]["deseq2_report"]
        log:
            join(_logdir("deseq2_analysis"), "{matrix_type}_{comparison}.log")
        shell:
            r"""
            set -exo pipefail
            mkdir -p $(dirname {log}) $(dirname {output.results_tsv})
            exec > >(tee -a {log}) 2>&1

            python3 {SCRIPTS_DIR}/run_deseq2_analysis.py \
              --count-matrix {input.count_matrix} \
              --manifest {input.manifest} \
              --contrasts-file {input.contrasts_file} \
              --comparison {params.comparison} \
              --group1 {params.group1} \
              --group2 {params.group2} \
              --matrix-type {params.matrix_type} \
              --alpha {params.alpha} \
              --lfc-threshold {params.lfc_threshold} \
              --min-replicates {params.min_replicates} \
              --size-factor-type {params.size_factor_type} \
              --output-results {output.results_tsv} \
              --output-volcano {output.volcano_png}
            """


def _deseq2_contrast_group(comparison, group_key):
    """Extract group name from contrasts list"""
    for contrast in DESEQ2_CONTRASTS:
        if contrast["comparison"] == comparison:
            return contrast.get(group_key, "")
    return ""
