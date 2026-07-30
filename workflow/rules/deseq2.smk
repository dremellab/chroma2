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

# Build DESeq2 outputs if enabled
DESEQ2_OUTPUTS = []
if DESEQ2_ENABLED and DESEQ2_CONTRASTS:
    for contrast in DESEQ2_CONTRASTS:
        comparison = contrast["comparison"]
        DESEQ2_OUTPUTS.append(join(DESEQ2_OUTDIR, comparison, f"{comparison}.deseq2_report.html"))


if DESEQ2_ENABLED:

    rule deseq2_contrast_report:
        input:
            config_yaml=configfilepath,
            manifest=MANIFEST_FILE,
            contrasts=DESEQ2_CONTRASTS_FILE,
            script=join(SCRIPTS_DIR, "run_deseq2_contrast_report.R"),
            report_template=join(SCRIPTS_DIR, "deseq2_contrast_report.Rmd"),
            # Include all available matrices as inputs
            matrices=list(DESEQ2_AVAILABLE_MATRICES.values()),
        output:
            report=join(DESEQ2_OUTDIR, "{comparison}", "{comparison}.deseq2_report.html"),
        params:
            group1=lambda wc: shlex.quote(_deseq2_contrast_group(wc.comparison, "group1")),
            group2=lambda wc: shlex.quote(_deseq2_contrast_group(wc.comparison, "group2")),
            comparison=lambda wc: shlex.quote(wc.comparison),
            outdir=lambda wc: shlex.quote(join(DESEQ2_OUTDIR, wc.comparison)),
            # Map available matrices to R script parameters
            # For now, pass gene as host-matrix, others as virus-matrices (flexible approach)
            host_matrix=shlex.quote(DESEQ2_AVAILABLE_MATRICES.get("gene", "")),
            host_output=shlex.quote(
                join(DESEQ2_OUTDIR, "{comparison}", "{comparison}.gene.deseq2_results.tsv")
            ),
            host_volcano=shlex.quote(
                join(DESEQ2_OUTDIR, "{comparison}", "{comparison}.gene.volcano.png")
            ),
            virus_matrices=shlex.quote(
                ";;;".join([str(v) for k, v in DESEQ2_AVAILABLE_MATRICES.items() if k != "gene"])
            ),
            virus_labels=shlex.quote(
                ";;;".join([k for k in DESEQ2_AVAILABLE_MATRICES.keys() if k != "gene"])
            ),
            virus_outputs=shlex.quote(
                ";;;".join(
                    [
                        join(DESEQ2_OUTDIR, "{comparison}", f"{{comparison}}.{k}.deseq2_results.tsv")
                        for k in DESEQ2_AVAILABLE_MATRICES.keys()
                        if k != "gene"
                    ]
                )
            ),
            virus_volcanos=shlex.quote(
                ";;;".join(
                    [
                        join(DESEQ2_OUTDIR, "{comparison}", f"{{comparison}}.{k}.volcano.png")
                        for k in DESEQ2_AVAILABLE_MATRICES.keys()
                        if k != "gene"
                    ]
                )
            ),
        threads:
            _get_threads("deseq2_contrast_report", profile_config)
        container:
            config["containers"]["deseq2_report"]
        log:
            join(_logdir("deseq2_contrast_report"), "{comparison}.log")
        shell:
            r"""
            set -exo pipefail
            mkdir -p $(dirname {log}) {params.outdir}
            exec > >(tee -a {log}) 2>&1
            Rscript "{input.script}" \
              --config-yaml "{input.config_yaml}" \
              --manifest "{input.manifest}" \
              --contrast-file "{input.contrasts}" \
              --comparison {params.comparison} \
              --group1 {params.group1} \
              --group2 {params.group2} \
              --report-template "{input.report_template}" \
              --report-output "{output.report}" \
              --host-matrix {params.host_matrix} \
              --host-output {params.host_output} \
              --host-volcano {params.host_volcano} \
              --virus-labels {params.virus_labels} \
              --virus-matrices {params.virus_matrices} \
              --virus-outputs {params.virus_outputs} \
              --virus-volcanos {params.virus_volcanos}
            """


def _deseq2_contrast_group(comparison, group_key):
    """Extract group name from contrasts list"""
    for contrast in DESEQ2_CONTRASTS:
        if contrast["comparison"] == comparison:
            return contrast.get(group_key, "")
    return ""
