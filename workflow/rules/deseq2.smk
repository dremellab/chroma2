###################################################################################
# DESeq2 contrast reporting for host + virus aggregate genrich matrices
###################################################################################

import shlex


DESEQ2_HOST_MATRIX = join(
    RESULTSDIR,
    "tn5_motif",
    f"host.genrich.tn5_gene_count_matrix.{DESEQ2_HOST_GENE_FLANK_SIZE}bp.tsv",
)
DESEQ2_VIRUS_MATRICES = {
    virus: join(
        RESULTSDIR,
        "tn5_motif",
        f"{virus}.genrich.tn5_bin_count_matrix.{DESEQ2_VIRUS_BIN_SIZE}bp.tsv",
    )
    for virus in VIRUS_LIST
}
DESEQ2_TRNA_MATRIX = (
    join(
        RESULTSDIR,
        "tn5_motif",
        f"host.genrich.tn5_trna_gene_count_matrix.{DESEQ2_TRNA_GENE_FLANK_SIZE}bp.tsv",
    )
    if HOST != "" and TRNAS_GTF != ""
    else ""
)
DESEQ2_TRNA_LABEL = f"host.trna.genrich.{DESEQ2_TRNA_GENE_FLANK_SIZE}bp"


def _deseq2_host_tsv(comparison):
    return join(
        DESEQ2_OUTDIR,
        comparison,
        f"{comparison}.host.genrich.{DESEQ2_HOST_GENE_FLANK_SIZE}bp.deseq2.tsv",
    )


def _deseq2_virus_tsv(comparison, virus):
    return join(
        DESEQ2_OUTDIR,
        comparison,
        f"{comparison}.{virus}.genrich.{DESEQ2_VIRUS_BIN_SIZE}bp.deseq2.tsv",
    )


def _deseq2_trna_tsv(comparison):
    return join(
        DESEQ2_OUTDIR,
        comparison,
        f"{comparison}.{DESEQ2_TRNA_LABEL}.deseq2.tsv",
    )


def _deseq2_host_volcano_png(comparison):
    return join(
        DESEQ2_OUTDIR,
        comparison,
        f"{comparison}.host.genrich.{DESEQ2_HOST_GENE_FLANK_SIZE}bp.volcano.png",
    )


def _deseq2_virus_volcano_png(comparison, virus):
    return join(
        DESEQ2_OUTDIR,
        comparison,
        f"{comparison}.{virus}.genrich.{DESEQ2_VIRUS_BIN_SIZE}bp.volcano.png",
    )


def _deseq2_trna_volcano_png(comparison):
    return join(
        DESEQ2_OUTDIR,
        comparison,
        f"{comparison}.{DESEQ2_TRNA_LABEL}.volcano.png",
    )


def _deseq2_report_html(comparison):
    return join(DESEQ2_OUTDIR, comparison, f"{comparison}.report.html")


def _deseq2_output_dir(comparison):
    return join(DESEQ2_OUTDIR, comparison)


def _deseq2_contrast_record(wildcards):
    return DESEQ2_CONTRASTS_BY_COMPARISON[wildcards.comparison]


def _deseq2_host_args(wildcards):
    if HOST == "":
        return ""
    return (
        f"--host-matrix {shlex.quote(DESEQ2_HOST_MATRIX)} "
        f"--host-output {shlex.quote(_deseq2_host_tsv(wildcards.comparison))}"
    )


def _deseq2_host_volcano_args(wildcards):
    if HOST == "":
        return ""
    return f"--host-volcano {shlex.quote(_deseq2_host_volcano_png(wildcards.comparison))}"


def _deseq2_joined_volcano_paths(wildcards=None):
    if wildcards is None:
        return ""
    entries = [_deseq2_virus_volcano_png(wildcards.comparison, v) for v in VIRUS_LIST]
    if DESEQ2_TRNA_MATRIX != "":
        entries.append(_deseq2_trna_volcano_png(wildcards.comparison))
    return ";;;".join(entries) if entries else ""


def _deseq2_joined_virus_paths(wildcards=None):
    if wildcards is None:
        entries = list(DESEQ2_VIRUS_MATRICES[virus] for virus in VIRUS_LIST)
        if DESEQ2_TRNA_MATRIX != "":
            entries.append(DESEQ2_TRNA_MATRIX)
        return ";;;".join(entries) if entries else ""
    else:
        entries = list(_deseq2_virus_tsv(wildcards.comparison, virus) for virus in VIRUS_LIST)
        if DESEQ2_TRNA_MATRIX != "":
            entries.append(_deseq2_trna_tsv(wildcards.comparison))
        return ";;;".join(entries) if entries else ""


DESEQ2_OUTPUTS = []
if DESEQ2_ENABLED:
    for contrast in DESEQ2_CONTRASTS:
        comparison = contrast["comparison"]
        if HOST != "":
            DESEQ2_OUTPUTS.append(_deseq2_host_tsv(comparison))
            DESEQ2_OUTPUTS.append(_deseq2_host_volcano_png(comparison))
        if DESEQ2_TRNA_MATRIX != "":
            DESEQ2_OUTPUTS.append(_deseq2_trna_tsv(comparison))
            DESEQ2_OUTPUTS.append(_deseq2_trna_volcano_png(comparison))
        for virus in VIRUS_LIST:
            DESEQ2_OUTPUTS.append(_deseq2_virus_tsv(comparison, virus))
            DESEQ2_OUTPUTS.append(_deseq2_virus_volcano_png(comparison, virus))
        DESEQ2_OUTPUTS.append(_deseq2_report_html(comparison))


if DESEQ2_ENABLED:
    rule deseq2_contrast_report:
        input:
            config_yaml=configfilepath,
            manifest=MANIFEST_FILE,
            contrasts=DESEQ2_CONTRASTS_FILE,
            script=join(SCRIPTS_DIR, "run_deseq2_contrast_report.R"),
            report_template=join(SCRIPTS_DIR, "deseq2_contrast_report.Rmd"),
            host_matrix=([DESEQ2_HOST_MATRIX] if HOST != "" else []),
            virus_matrices=(
                [DESEQ2_VIRUS_MATRICES[virus] for virus in VIRUS_LIST]
                + ([DESEQ2_TRNA_MATRIX] if DESEQ2_TRNA_MATRIX != "" else [])
            ),
        output:
            report=_deseq2_report_html("{comparison}"),
            host_tsv=(
                _deseq2_host_tsv("{comparison}")
                if HOST != ""
                else []
            ),
            host_volcano=(
                _deseq2_host_volcano_png("{comparison}")
                if HOST != ""
                else []
            ),
            trna_tsv=(
                _deseq2_trna_tsv("{comparison}")
                if DESEQ2_TRNA_MATRIX != ""
                else []
            ),
            trna_volcano=(
                _deseq2_trna_volcano_png("{comparison}")
                if DESEQ2_TRNA_MATRIX != ""
                else []
            ),
            virus_tsv=[
                _deseq2_virus_tsv("{comparison}", virus)
                for virus in VIRUS_LIST
            ],
            virus_volcano=[
                _deseq2_virus_volcano_png("{comparison}", virus)
                for virus in VIRUS_LIST
            ],
        params:
            group1=lambda wc: shlex.quote(_deseq2_contrast_record(wc)["group1"]),
            group2=lambda wc: shlex.quote(_deseq2_contrast_record(wc)["group2"]),
            comparison=lambda wc: shlex.quote(wc.comparison),
            outdir=lambda wc: shlex.quote(_deseq2_output_dir(wc.comparison)),
            host_args=_deseq2_host_args,
            host_volcano_args=_deseq2_host_volcano_args,
            virus_labels=shlex.quote(
                ";;;".join(
                    list(VIRUS_LIST) + ([DESEQ2_TRNA_LABEL] if DESEQ2_TRNA_MATRIX != "" else [])
                )
            ),
            virus_matrices=shlex.quote(_deseq2_joined_virus_paths()),
            virus_outputs=lambda wc: shlex.quote(_deseq2_joined_virus_paths(wc)),
            virus_volcanos=lambda wc: shlex.quote(_deseq2_joined_volcano_paths(wc)),
        threads:
            _get_threads("deseq2_contrast_report", profile_config)
        container:
            config["containers"]["deseq2_report"]
        log:
            join(_logdir("deseq2_contrast_report"), "{comparison}.log")
        shell:
            r"""
            set -exo pipefail
            mkdir -p $(dirname {log})
            exec > >(tee -a {log}) 2>&1
            mkdir -p {params.outdir}
            Rscript "{input.script}" \
              --config-yaml "{input.config_yaml}" \
              --manifest "{input.manifest}" \
              --contrast-file "{input.contrasts}" \
              --comparison {params.comparison} \
              --group1 {params.group1} \
              --group2 {params.group2} \
              --report-template "{input.report_template}" \
              --report-output "{output.report}" \
              {params.host_args} \
              {params.host_volcano_args} \
              --virus-labels {params.virus_labels} \
              --virus-matrices {params.virus_matrices} \
              --virus-outputs {params.virus_outputs} \
              --virus-volcanos {params.virus_volcanos}
            """
