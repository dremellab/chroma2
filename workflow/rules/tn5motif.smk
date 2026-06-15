###################################################################################
# Tn5 motif extraction from caller-specific BAMs (host + per-virus)
###################################################################################

TN5_WINDOW_SIZE = 2 * int(config.get("tn5_motif", {}).get("flank_size", 10)) + 1
TN5_VIRUS_BIN_SIZE = int(config.get("tn5_motif", {}).get("virus_bin_size", 100))
TN5_HOST_GENE_FLANK_SIZE = int(config.get("tn5_motif", {}).get("host_gene_flank_size", 250))
TN5_TRNA_GENE_FLANK_SIZE = int(config.get("tn5_motif", {}).get("trna_gene_flank_size", 100))
TN5_TRNA_GENE_MAX_SIZE = int(config.get("tn5_motif", {}).get("trna_gene_max_size", 1000))
TN5_MAPQ_MIN = int(config.get("tn5_motif", {}).get("mapq_min", 0))

TN5_EXCLUDE_SECONDARY = bool(config.get("tn5_motif", {}).get("exclude_secondary", False))
TN5_EXCLUDE_SUPPLEMENTARY = bool(
    config.get("tn5_motif", {}).get("exclude_supplementary", False)
)
TN5_DEDUP = bool(config.get("tn5_motif", {}).get("dedup", False))
TN5_LOGO_FORMAT = str(config.get("tn5_motif", {}).get("logo_format", "png"))
TN5_GENERATE_LOGO = bool(config.get("tn5_motif", {}).get("generate_logo", False))
TN5_FRACTIONAL_COUNTING = bool(config.get("tn5_motif", {}).get("fractional_counting", False))
TN5_CALLERS = ("genrich", "macs2")


def _tn5_prefix(sample, group, caller):
    return f"{sample}.{group}.{caller}"


def _tn5_exact_bed(sample, caller, group):
    prefix = _tn5_prefix(sample, group, caller)
    return join(RESULTSDIR, sample, "tn5_motif", caller, group, f"{prefix}.tn5_sites.1bp.bed.gz")


def _tn5_flank_bed(sample, caller, group):
    prefix = _tn5_prefix(sample, group, caller)
    return join(
        RESULTSDIR,
        sample,
        "tn5_motif",
        caller,
        group,
        f"{prefix}.tn5_sites.{TN5_WINDOW_SIZE}bp.bed",
    )


def _tn5_fasta(sample, caller, group):
    prefix = _tn5_prefix(sample, group, caller)
    return join(
        RESULTSDIR,
        sample,
        "tn5_motif",
        caller,
        group,
        f"{prefix}.tn5_sites.{TN5_WINDOW_SIZE}bp.fa",
    )


def _tn5_pfm(sample, caller, group):
    prefix = _tn5_prefix(sample, group, caller)
    return join(
        RESULTSDIR,
        sample,
        "tn5_motif",
        caller,
        group,
        f"{prefix}.tn5_sites.{TN5_WINDOW_SIZE}bp.pfm.tsv.gz",
    )


def _tn5_logo(sample, caller, group):
    prefix = _tn5_prefix(sample, group, caller)
    return join(
        RESULTSDIR,
        sample,
        "tn5_motif",
        caller,
        group,
        f"{prefix}.tn5_sites.{TN5_WINDOW_SIZE}bp.logo.{TN5_LOGO_FORMAT}",
    )


def _caller_bam(sample, caller, group):
    if caller == "genrich":
        if group == "host":
            return join(RESULTSDIR, sample, "align", f"{sample}.aligned.host.qname.bam")
        return join(
            RESULTSDIR,
            sample,
            "align",
            f"{sample}.aligned.virus.{group}.qname.bam",
        )
    if caller == "macs2":
        if group == "host":
            return join(RESULTSDIR, sample, "postprocess", f"{sample}.host.bam")
        return join(
            RESULTSDIR,
            sample,
            "postprocess",
            f"{sample}.virus.{group}.bam",
        )
    raise ValueError(f"Unsupported Tn5 caller: {caller}")


TN5_MOTIF_OUTPUTS = []

if TN5_GENERATE_LOGO:
    for caller in TN5_CALLERS:
        if HOST != "":
            TN5_MOTIF_OUTPUTS.extend(
                expand(_tn5_exact_bed("{sample}", caller, "host"), sample=SAMPLES)
            )
            TN5_MOTIF_OUTPUTS.extend(
                expand(_tn5_pfm("{sample}", caller, "host"), sample=SAMPLES)
            )
            TN5_MOTIF_OUTPUTS.extend(
                expand(_tn5_logo("{sample}", caller, "host"), sample=SAMPLES)
            )
        if len(VIRUS_LIST) > 0:
            TN5_MOTIF_OUTPUTS.extend(
                expand(_tn5_exact_bed("{sample}", caller, "{virus}"), sample=SAMPLES, virus=VIRUS_LIST)
            )
            TN5_MOTIF_OUTPUTS.extend(
                expand(_tn5_pfm("{sample}", caller, "{virus}"), sample=SAMPLES, virus=VIRUS_LIST)
            )
            TN5_MOTIF_OUTPUTS.extend(
                expand(_tn5_logo("{sample}", caller, "{virus}"), sample=SAMPLES, virus=VIRUS_LIST)
            )


if TN5_GENERATE_LOGO:

    rule extract_tn5_motifs_host:
        input:
            bam=lambda wc: _caller_bam(wc.sample, wc.caller, "host"),
            fasta=REF_FA,
            host_regions=REF_REGIONS_HOST,
        output:
            exact_bed=_tn5_exact_bed("{sample}", "{caller}", "host"),
            pfm=_tn5_pfm("{sample}", "{caller}", "host"),
            logo=_tn5_logo("{sample}", "{caller}", "host"),
        params:
            outdir=join(RESULTSDIR, "{sample}", "tn5_motif", "{caller}"),
            flank_size=str(config.get("tn5_motif", {}).get("flank_size", 10)),
            dedup_arg=("--dedup" if TN5_DEDUP else ""),
            mapq_min=str(TN5_MAPQ_MIN),
            exclude_secondary_arg=("--exclude-secondary" if TN5_EXCLUDE_SECONDARY else ""),
            exclude_supplementary_arg=(
                "--exclude-supplementary" if TN5_EXCLUDE_SUPPLEMENTARY else ""
            ),
            logo_format=TN5_LOGO_FORMAT,
            fractional_arg=("--fractional-counting" if TN5_FRACTIONAL_COUNTING else ""),
        wildcard_constraints:
            caller="genrich|macs2",
        threads:
            _get_threads("extract_tn5_motifs_host", profile_config)
        container:
            config["containers"]["pysam"]
        log:
            join(_logdir("extract_tn5_motifs_host"), "{sample}.{caller}.log")
        shell:
            r"""
            set -exo pipefail
            mkdir -p $(dirname {log})
            exec > >(tee -a {log}) 2>&1
            rm -rf {params.outdir}/host
            mkdir -p {params.outdir}
            python {SCRIPTS_DIR}/extract_tn5_motifs.py \
              --bam {input.bam} \
              --fasta {input.fasta} \
              --sample {wildcards.sample} \
              --scenario-name {wildcards.caller} \
              --outdir {params.outdir} \
              --host-regions {input.host_regions} \
              --flank-size {params.flank_size} \
              --threads {threads} \
              --mapq-min {params.mapq_min} \
              {params.dedup_arg} \
              {params.exclude_secondary_arg} \
              {params.exclude_supplementary_arg} \
              --logo-format {params.logo_format} \
              --generate-logo \
              {params.fractional_arg}
            """

    rule extract_tn5_motifs_virus:
        input:
            bam=lambda wc: _caller_bam(wc.sample, wc.caller, wc.virus),
            fasta=REF_FA,
            virus_regions=join(FASTAS_GTFS_DIR, "{virus}.fa.regions"),
        output:
            exact_bed=_tn5_exact_bed("{sample}", "{caller}", "{virus}"),
            pfm=_tn5_pfm("{sample}", "{caller}", "{virus}"),
            logo=_tn5_logo("{sample}", "{caller}", "{virus}"),
        params:
            outdir=join(RESULTSDIR, "{sample}", "tn5_motif", "{caller}"),
            virus_regions_arg=lambda wc, input: f"--virus-regions {wc.virus}={input.virus_regions}",
            flank_size=str(config.get("tn5_motif", {}).get("flank_size", 10)),
            dedup_arg=("--dedup" if TN5_DEDUP else ""),
            mapq_min=str(TN5_MAPQ_MIN),
            exclude_secondary_arg=("--exclude-secondary" if TN5_EXCLUDE_SECONDARY else ""),
            exclude_supplementary_arg=(
                "--exclude-supplementary" if TN5_EXCLUDE_SUPPLEMENTARY else ""
            ),
            logo_format=TN5_LOGO_FORMAT,
            fractional_arg=("--fractional-counting" if TN5_FRACTIONAL_COUNTING else ""),
        wildcard_constraints:
            caller="genrich|macs2",
        threads:
            _get_threads("extract_tn5_motifs_virus", profile_config)
        container:
            config["containers"]["pysam"]
        log:
            join(_logdir("extract_tn5_motifs_virus"), "{sample}.{virus}.{caller}.log")
        shell:
            r"""
            set -exo pipefail
            mkdir -p $(dirname {log})
            exec > >(tee -a {log}) 2>&1
            rm -rf {params.outdir}/{wildcards.virus}
            mkdir -p {params.outdir}
            python {SCRIPTS_DIR}/extract_tn5_motifs.py \
              --bam {input.bam} \
              --fasta {input.fasta} \
              --sample {wildcards.sample} \
              --scenario-name {wildcards.caller} \
              --outdir {params.outdir} \
              {params.virus_regions_arg} \
              --flank-size {params.flank_size} \
              --threads {threads} \
              --mapq-min {params.mapq_min} \
              {params.dedup_arg} \
              {params.exclude_secondary_arg} \
              {params.exclude_supplementary_arg} \
              --logo-format {params.logo_format} \
              --generate-logo \
              {params.fractional_arg}
            """
