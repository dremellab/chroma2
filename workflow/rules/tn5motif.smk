###################################################################################
# Tn5 motif extraction from aligned BAMs (host + per-virus)
###################################################################################

def _tn5_primary_scenario_name():
    parts = [f"mapq{int(config.get('tn5_motif', {}).get('mapq_min', 0))}"]
    if config.get("tn5_motif", {}).get("exclude_secondary", False):
        parts.append("nosecondary")
    if config.get("tn5_motif", {}).get("exclude_supplementary", False):
        parts.append("nosupplementary")
    if config.get("tn5_motif", {}).get("dedup", False):
        parts.append("dedup")
    if parts == ["mapq0"]:
        return "default"
    return "_".join(parts)


TN5_WINDOW_SIZE = 2 * int(config.get("tn5_motif", {}).get("flank_size", 10)) + 1
TN5_PRIMARY_SCENARIO = _tn5_primary_scenario_name()
TN5_STRICT_SCENARIO = "mapq20_nosecondary_nosupplementary"
TN5_SCENARIOS = [TN5_PRIMARY_SCENARIO]
if TN5_PRIMARY_SCENARIO != TN5_STRICT_SCENARIO:
    TN5_SCENARIOS.append(TN5_STRICT_SCENARIO)

TN5_SCENARIO_CONFIG = {
    TN5_PRIMARY_SCENARIO: {
        "mapq_min": int(config.get("tn5_motif", {}).get("mapq_min", 0)),
        "exclude_secondary": bool(config.get("tn5_motif", {}).get("exclude_secondary", False)),
        "exclude_supplementary": bool(config.get("tn5_motif", {}).get("exclude_supplementary", False)),
    },
    TN5_STRICT_SCENARIO: {
        "mapq_min": 20,
        "exclude_secondary": True,
        "exclude_supplementary": True,
    },
}

TN5_GROUPS = (["host"] if HOST != "" else []) + VIRUS_LIST

TN5_MOTIF_OUTPUTS = []
for scenario in TN5_SCENARIOS:
    TN5_MOTIF_OUTPUTS.extend(
        expand(
            join(
                RESULTSDIR,
                "{sample}",
                "tn5_motif",
                scenario,
                "{group}",
                "{sample}.{group}." + scenario + ".tn5_sites.1bp.bed",
            ),
            sample=SAMPLES,
            group=TN5_GROUPS,
        )
    )
    TN5_MOTIF_OUTPUTS.extend(
        expand(
            join(
                RESULTSDIR,
                "{sample}",
                "tn5_motif",
                scenario,
                "{group}",
                "{sample}.{group}." + scenario + f".tn5_sites.{TN5_WINDOW_SIZE}bp.bed",
            ),
            sample=SAMPLES,
            group=TN5_GROUPS,
        )
    )
    TN5_MOTIF_OUTPUTS.extend(
        expand(
            join(
                RESULTSDIR,
                "{sample}",
                "tn5_motif",
                scenario,
                "{group}",
                "{sample}.{group}." + scenario + f".tn5_sites.{TN5_WINDOW_SIZE}bp.fa",
            ),
            sample=SAMPLES,
            group=TN5_GROUPS,
        )
    )
    TN5_MOTIF_OUTPUTS.extend(
        expand(
            join(
                RESULTSDIR,
                "{sample}",
                "tn5_motif",
                scenario,
                "{group}",
                "{sample}.{group}." + scenario + f".tn5_sites.{TN5_WINDOW_SIZE}bp.pfm.tsv",
            ),
            sample=SAMPLES,
            group=TN5_GROUPS,
        )
    )
    TN5_MOTIF_OUTPUTS.extend(
        expand(
            join(
                RESULTSDIR,
                "{sample}",
                "tn5_motif",
                scenario,
                "{group}",
                "{sample}.{group}." + scenario + f".tn5_sites.{TN5_WINDOW_SIZE}bp.logo.png",
            ),
            sample=SAMPLES,
            group=TN5_GROUPS,
        )
    )


rule extract_tn5_motifs:
    input:
        bam=join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.sorted.bam"),
        fasta=REF_FA,
        host_regions=([REF_REGIONS_HOST] if HOST != "" else []),
        virus_regions=expand(join(FASTAS_GTFS_DIR, "{virus}.fa.regions"), virus=VIRUS_LIST),
    output:
        exact_beds=expand(
            join(
                RESULTSDIR,
                "{sample}",
                "tn5_motif",
                "{scenario}",
                "{group}",
                "{sample}.{group}.{scenario}.tn5_sites.1bp.bed",
            ),
            group=TN5_GROUPS,
            allow_missing=True,
        ),
        flank_beds=expand(
            join(
                RESULTSDIR,
                "{sample}",
                "tn5_motif",
                "{scenario}",
                "{group}",
                "{sample}.{group}.{scenario}.tn5_sites." + str(TN5_WINDOW_SIZE) + "bp.bed",
            ),
            group=TN5_GROUPS,
            allow_missing=True,
        ),
        fastas=expand(
            join(
                RESULTSDIR,
                "{sample}",
                "tn5_motif",
                "{scenario}",
                "{group}",
                "{sample}.{group}.{scenario}.tn5_sites." + str(TN5_WINDOW_SIZE) + "bp.fa",
            ),
            group=TN5_GROUPS,
            allow_missing=True,
        ),
        pfms=expand(
            join(
                RESULTSDIR,
                "{sample}",
                "tn5_motif",
                "{scenario}",
                "{group}",
                "{sample}.{group}.{scenario}.tn5_sites." + str(TN5_WINDOW_SIZE) + "bp.pfm.tsv",
            ),
            group=TN5_GROUPS,
            allow_missing=True,
        ),
        logos=expand(
            join(
                RESULTSDIR,
                "{sample}",
                "tn5_motif",
                "{scenario}",
                "{group}",
                "{sample}.{group}.{scenario}.tn5_sites." + str(TN5_WINDOW_SIZE) + "bp.logo.png",
            ),
            group=TN5_GROUPS,
            allow_missing=True,
        ),
    params:
        outdir=join(RESULTSDIR, "{sample}", "tn5_motif", "{scenario}"),
        host_regions_arg=(f"--host-regions {REF_REGIONS_HOST}" if HOST != "" else ""),
        virus_regions_args=" ".join(
            [
                f"--virus-regions {virus}={join(FASTAS_GTFS_DIR, virus + '.fa.regions')}"
                for virus in VIRUS_LIST
            ]
        ),
        flank_size=str(config.get("tn5_motif", {}).get("flank_size", 10)),
        dedup_arg=("--dedup" if config.get("tn5_motif", {}).get("dedup", False) else ""),
        mapq_min=lambda wc: str(TN5_SCENARIO_CONFIG[wc.scenario]["mapq_min"]),
        exclude_secondary_arg=lambda wc: (
            "--exclude-secondary" if TN5_SCENARIO_CONFIG[wc.scenario]["exclude_secondary"] else ""
        ),
        exclude_supplementary_arg=lambda wc: (
            "--exclude-supplementary" if TN5_SCENARIO_CONFIG[wc.scenario]["exclude_supplementary"] else ""
        ),
        logo_format=str(config.get("tn5_motif", {}).get("logo_format", "png")),
    threads:
        _get_threads("extract_tn5_motifs", profile_config)
    container:
        config["containers"]["pysam"]
    log:
        join(_logdir("extract_tn5_motifs"), "{sample}.{scenario}.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        python {SCRIPTS_DIR}/extract_tn5_motifs.py \
          --bam {input.bam} \
          --fasta {input.fasta} \
          --sample {wildcards.sample} \
          --scenario-name {wildcards.scenario} \
          --outdir {params.outdir} \
          {params.host_regions_arg} \
          {params.virus_regions_args} \
          --flank-size {params.flank_size} \
          --threads {threads} \
          --mapq-min {params.mapq_min} \
          {params.dedup_arg} \
          {params.exclude_secondary_arg} \
          {params.exclude_supplementary_arg} \
          --logo-format {params.logo_format}
        """
