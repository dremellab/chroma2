###################################################################################
# Tn5 motif extraction from caller-specific BAMs (host + per-virus)
###################################################################################

TN5_WINDOW_SIZE = 2 * int(config.get("tn5_motif", {}).get("flank_size", 10)) + 1
TN5_VIRUS_BIN_SIZE = int(config.get("tn5_motif", {}).get("virus_bin_size", 100))
TN5_HOST_GENE_FLANK_SIZE = int(config.get("tn5_motif", {}).get("host_gene_flank_size", 250))
TN5_MAPQ_MIN = int(config.get("tn5_motif", {}).get("mapq_min", 0))
TN5_EXCLUDE_SECONDARY = bool(config.get("tn5_motif", {}).get("exclude_secondary", False))
TN5_EXCLUDE_SUPPLEMENTARY = bool(
    config.get("tn5_motif", {}).get("exclude_supplementary", False)
)
TN5_DEDUP = bool(config.get("tn5_motif", {}).get("dedup", False))
TN5_LOGO_FORMAT = str(config.get("tn5_motif", {}).get("logo_format", "png"))
TN5_CALLERS = ("genrich", "macs2")


def _tn5_prefix(sample, group, caller):
    return f"{sample}.{group}.{caller}"


def _tn5_exact_bed(sample, caller, group):
    prefix = _tn5_prefix(sample, group, caller)
    return join(RESULTSDIR, sample, "tn5_motif", caller, group, f"{prefix}.tn5_sites.1bp.bed")


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
        f"{prefix}.tn5_sites.{TN5_WINDOW_SIZE}bp.pfm.tsv",
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
for caller in TN5_CALLERS:
    if HOST != "":
        TN5_MOTIF_OUTPUTS.extend(
            expand(_tn5_exact_bed("{sample}", caller, "host"), sample=SAMPLES)
        )
        TN5_MOTIF_OUTPUTS.extend(
            expand(_tn5_flank_bed("{sample}", caller, "host"), sample=SAMPLES)
        )
        TN5_MOTIF_OUTPUTS.extend(
            expand(_tn5_fasta("{sample}", caller, "host"), sample=SAMPLES)
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
            expand(_tn5_flank_bed("{sample}", caller, "{virus}"), sample=SAMPLES, virus=VIRUS_LIST)
        )
        TN5_MOTIF_OUTPUTS.extend(
            expand(_tn5_fasta("{sample}", caller, "{virus}"), sample=SAMPLES, virus=VIRUS_LIST)
        )
        TN5_MOTIF_OUTPUTS.extend(
            expand(_tn5_pfm("{sample}", caller, "{virus}"), sample=SAMPLES, virus=VIRUS_LIST)
        )
        TN5_MOTIF_OUTPUTS.extend(
            expand(_tn5_logo("{sample}", caller, "{virus}"), sample=SAMPLES, virus=VIRUS_LIST)
        )


TN5_BIN_COUNT_OUTPUTS = []
for caller in TN5_CALLERS:
    if len(VIRUS_LIST) > 0:
        TN5_BIN_COUNT_OUTPUTS.extend(
            expand(
                join(
                    RESULTSDIR,
                    "{sample}",
                    "tn5_motif",
                    caller,
                    "{virus}",
                    "{sample}.{virus}."
                    + caller
                    + f".tn5_sites.bin_counts.{TN5_VIRUS_BIN_SIZE}bp.tsv",
                ),
                sample=SAMPLES,
                virus=VIRUS_LIST,
            )
        )


TN5_HOST_GENE_BIN_OUTPUTS = []
TN5_HOST_BIN_COUNT_OUTPUTS = []
TN5_AGGREGATE_COUNT_MATRIX_OUTPUTS = []
if HOST != "":
    TN5_HOST_GENE_BIN_OUTPUTS.append(
        join(REF_DIR, f"ref.tn5.host_gene_bins.{TN5_HOST_GENE_FLANK_SIZE}bp.bed")
    )
    for caller in TN5_CALLERS:
        TN5_HOST_BIN_COUNT_OUTPUTS.extend(
            expand(
                join(
                    RESULTSDIR,
                    "{sample}",
                    "tn5_motif",
                    caller,
                    "host",
                    "{sample}.host."
                    + caller
                    + ".tn5_gene_counts."
                    + str(TN5_HOST_GENE_FLANK_SIZE)
                    + "bp.tsv",
                ),
                sample=SAMPLES,
            )
        )
        TN5_AGGREGATE_COUNT_MATRIX_OUTPUTS.append(
            join(
                RESULTSDIR,
                "tn5_motif",
                f"host.{caller}.tn5_gene_count_matrix.{TN5_HOST_GENE_FLANK_SIZE}bp.tsv",
            )
        )

for caller in TN5_CALLERS:
    if len(VIRUS_LIST) > 0:
        TN5_AGGREGATE_COUNT_MATRIX_OUTPUTS.extend(
            expand(
                join(
                    RESULTSDIR,
                    "tn5_motif",
                    "{virus}."
                    + caller
                    + f".tn5_bin_count_matrix.{TN5_VIRUS_BIN_SIZE}bp.tsv",
                ),
                virus=VIRUS_LIST,
            )
        )


rule build_host_tn5_gene_bins:
    input:
        host_gtf=join(FASTAS_GTFS_DIR, HOST + ".gtf"),
        host_regions=REF_REGIONS_HOST,
        chromsizes=join(REF_DIR, "ref.chrom.sizes.host.txt"),
        trna_gtf=([TRNAS_GTF] if HOST != "" and TRNAS_GTF != "" else []),
    output:
        bed=join(
            REF_DIR,
            "ref.tn5.host_gene_bins." + str(TN5_HOST_GENE_FLANK_SIZE) + "bp.bed",
        ),
    params:
        flank_size=str(config.get("tn5_motif", {}).get("host_gene_flank_size", 250)),
        trna_gtf_arg=lambda wc, input: (
            f"--trna-gtf {input.trna_gtf[0]}" if input.trna_gtf else ""
        ),
    threads:
        _get_threads("build_host_tn5_gene_bins", profile_config)
    container:
        config["containers"]["py311"]
    log:
        join(_logdir("build_host_tn5_gene_bins"), "host.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        python {SCRIPTS_DIR}/build_host_tn5_gene_bins.py \
          --host-gtf {input.host_gtf} \
          --host-regions {input.host_regions} \
          --chromsizes {input.chromsizes} \
          {params.trna_gtf_arg} \
          --flank-size {params.flank_size} \
          --output {output.bed}
        """


rule extract_tn5_motifs_host:
    input:
        bam=lambda wc: _caller_bam(wc.sample, wc.caller, "host"),
        fasta=REF_FA,
        host_regions=REF_REGIONS_HOST,
    output:
        exact_bed=_tn5_exact_bed("{sample}", "{caller}", "host"),
        flank_bed=_tn5_flank_bed("{sample}", "{caller}", "host"),
        fasta=_tn5_fasta("{sample}", "{caller}", "host"),
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
          --logo-format {params.logo_format}
        """


rule extract_tn5_motifs_virus:
    input:
        bam=lambda wc: _caller_bam(wc.sample, wc.caller, wc.virus),
        fasta=REF_FA,
        virus_regions=join(FASTAS_GTFS_DIR, "{virus}.fa.regions"),
    output:
        exact_bed=_tn5_exact_bed("{sample}", "{caller}", "{virus}"),
        flank_bed=_tn5_flank_bed("{sample}", "{caller}", "{virus}"),
        fasta=_tn5_fasta("{sample}", "{caller}", "{virus}"),
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
          --logo-format {params.logo_format}
        """


rule count_tn5_host_genrich_bins:
    input:
        bam=join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.host.qname.bam"),
        bins=join(
            REF_DIR,
            "ref.tn5.host_gene_bins." + str(TN5_HOST_GENE_FLANK_SIZE) + "bp.bed",
        ),
    output:
        tsv=join(
            RESULTSDIR,
            "{sample}",
            "tn5_motif",
            "genrich",
            "host",
            "{sample}.host.genrich.tn5_gene_counts."
            + str(TN5_HOST_GENE_FLANK_SIZE)
            + "bp.tsv",
        ),
    params:
        mapq_min=str(TN5_MAPQ_MIN),
        exclude_secondary_arg=("--exclude-secondary" if TN5_EXCLUDE_SECONDARY else ""),
        exclude_supplementary_arg=(
            "--exclude-supplementary" if TN5_EXCLUDE_SUPPLEMENTARY else ""
        ),
    threads:
        _get_threads("count_tn5_host_genrich_bins", profile_config)
    container:
        config["containers"]["pysam"]
    log:
        join(_logdir("count_tn5_host_genrich_bins"), "{sample}.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        python {SCRIPTS_DIR}/count_tn5_sites_in_bins.py \
          --bam {input.bam} \
          --bins-bed {input.bins} \
          --output {output.tsv} \
          --sample {wildcards.sample} \
          --threads {threads} \
          --mapq-min {params.mapq_min} \
          {params.exclude_secondary_arg} \
          {params.exclude_supplementary_arg}
        """


rule count_tn5_host_macs2_bins:
    input:
        bam=join(RESULTSDIR, "{sample}", "postprocess", "{sample}.host.bam"),
        bins=join(
            REF_DIR,
            "ref.tn5.host_gene_bins." + str(TN5_HOST_GENE_FLANK_SIZE) + "bp.bed",
        ),
    output:
        tsv=join(
            RESULTSDIR,
            "{sample}",
            "tn5_motif",
            "macs2",
            "host",
            "{sample}.host.macs2.tn5_gene_counts."
            + str(TN5_HOST_GENE_FLANK_SIZE)
            + "bp.tsv",
        ),
    params:
        mapq_min=str(TN5_MAPQ_MIN),
        exclude_secondary_arg=("--exclude-secondary" if TN5_EXCLUDE_SECONDARY else ""),
        exclude_supplementary_arg=(
            "--exclude-supplementary" if TN5_EXCLUDE_SUPPLEMENTARY else ""
        ),
    threads:
        _get_threads("count_tn5_host_macs2_bins", profile_config)
    container:
        config["containers"]["pysam"]
    log:
        join(_logdir("count_tn5_host_macs2_bins"), "{sample}.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        python {SCRIPTS_DIR}/count_tn5_sites_in_bins.py \
          --bam {input.bam} \
          --bins-bed {input.bins} \
          --output {output.tsv} \
          --sample {wildcards.sample} \
          --threads {threads} \
          --mapq-min {params.mapq_min} \
          {params.exclude_secondary_arg} \
          {params.exclude_supplementary_arg}
        """


rule build_tn5_host_genrich_count_matrix:
    input:
        counts=expand(
            join(
                RESULTSDIR,
                "{sample}",
                "tn5_motif",
                "genrich",
                "host",
                "{sample}.host.genrich.tn5_gene_counts."
                + str(TN5_HOST_GENE_FLANK_SIZE)
                + "bp.tsv",
            ),
            sample=SAMPLES,
        ),
    output:
        tsv=join(
            RESULTSDIR,
            "tn5_motif",
            "host.genrich.tn5_gene_count_matrix."
            + str(TN5_HOST_GENE_FLANK_SIZE)
            + "bp.tsv",
        ),
    threads:
        _get_threads("build_tn5_host_genrich_count_matrix", profile_config)
    container:
        config["containers"]["py311"]
    log:
        join(_logdir("build_tn5_host_genrich_count_matrix"), "host.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        mkdir -p $(dirname {output.tsv})
        python {SCRIPTS_DIR}/build_tn5_count_matrix.py \
          --counts {input.counts} \
          --output {output.tsv}
        """


rule build_tn5_host_macs2_count_matrix:
    input:
        counts=expand(
            join(
                RESULTSDIR,
                "{sample}",
                "tn5_motif",
                "macs2",
                "host",
                "{sample}.host.macs2.tn5_gene_counts."
                + str(TN5_HOST_GENE_FLANK_SIZE)
                + "bp.tsv",
            ),
            sample=SAMPLES,
        ),
    output:
        tsv=join(
            RESULTSDIR,
            "tn5_motif",
            "host.macs2.tn5_gene_count_matrix."
            + str(TN5_HOST_GENE_FLANK_SIZE)
            + "bp.tsv",
        ),
    threads:
        _get_threads("build_tn5_host_macs2_count_matrix", profile_config)
    container:
        config["containers"]["py311"]
    log:
        join(_logdir("build_tn5_host_macs2_count_matrix"), "host.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        mkdir -p $(dirname {output.tsv})
        python {SCRIPTS_DIR}/build_tn5_count_matrix.py \
          --counts {input.counts} \
          --output {output.tsv}
        """


rule count_tn5_bins_genrich_virus:
    input:
        bam=join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.virus.{virus}.qname.bam"),
        fasta=REF_FA,
        virus_regions=join(FASTAS_GTFS_DIR, "{virus}.fa.regions"),
    output:
        tsv=join(
            RESULTSDIR,
            "{sample}",
            "tn5_motif",
            "genrich",
            "{virus}",
            "{sample}.{virus}.genrich.tn5_sites.bin_counts."
            + str(TN5_VIRUS_BIN_SIZE)
            + "bp.tsv",
        ),
    params:
        outdir=join(RESULTSDIR, "{sample}", "tn5_motif", "genrich"),
        virus_regions_arg=lambda wc, input: f"--virus-regions {wc.virus}={input.virus_regions}",
        mapq_min=str(TN5_MAPQ_MIN),
        exclude_secondary_arg=("--exclude-secondary" if TN5_EXCLUDE_SECONDARY else ""),
        exclude_supplementary_arg=(
            "--exclude-supplementary" if TN5_EXCLUDE_SUPPLEMENTARY else ""
        ),
    threads:
        _get_threads("count_tn5_bins_genrich_virus", profile_config)
    container:
        config["containers"]["pysam"]
    log:
        join(_logdir("count_tn5_bins_genrich_virus"), "{sample}.{virus}.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        python {SCRIPTS_DIR}/extract_tn5_motifs.py \
          --bam {input.bam} \
          --fasta {input.fasta} \
          --sample {wildcards.sample} \
          --scenario-name genrich \
          --outdir {params.outdir} \
          {params.virus_regions_arg} \
          --virus-bin-size {TN5_VIRUS_BIN_SIZE} \
          --bin-counts-only \
          --bin-counts-output {output.tsv} \
          --threads {threads} \
          --mapq-min {params.mapq_min} \
          {params.exclude_secondary_arg} \
          {params.exclude_supplementary_arg}
        """


rule count_tn5_bins_macs2_virus:
    input:
        bam=join(RESULTSDIR, "{sample}", "postprocess", "{sample}.virus.{virus}.bam"),
        fasta=REF_FA,
        virus_regions=join(FASTAS_GTFS_DIR, "{virus}.fa.regions"),
    output:
        tsv=join(
            RESULTSDIR,
            "{sample}",
            "tn5_motif",
            "macs2",
            "{virus}",
            "{sample}.{virus}.macs2.tn5_sites.bin_counts."
            + str(TN5_VIRUS_BIN_SIZE)
            + "bp.tsv",
        ),
    params:
        outdir=join(RESULTSDIR, "{sample}", "tn5_motif", "macs2"),
        virus_regions_arg=lambda wc, input: f"--virus-regions {wc.virus}={input.virus_regions}",
        mapq_min=str(TN5_MAPQ_MIN),
        exclude_secondary_arg=("--exclude-secondary" if TN5_EXCLUDE_SECONDARY else ""),
        exclude_supplementary_arg=(
            "--exclude-supplementary" if TN5_EXCLUDE_SUPPLEMENTARY else ""
        ),
    threads:
        _get_threads("count_tn5_bins_macs2_virus", profile_config)
    container:
        config["containers"]["pysam"]
    log:
        join(_logdir("count_tn5_bins_macs2_virus"), "{sample}.{virus}.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        python {SCRIPTS_DIR}/extract_tn5_motifs.py \
          --bam {input.bam} \
          --fasta {input.fasta} \
          --sample {wildcards.sample} \
          --scenario-name macs2 \
          --outdir {params.outdir} \
          {params.virus_regions_arg} \
          --virus-bin-size {TN5_VIRUS_BIN_SIZE} \
          --bin-counts-only \
          --bin-counts-output {output.tsv} \
          --threads {threads} \
          --mapq-min {params.mapq_min} \
          {params.exclude_secondary_arg} \
          {params.exclude_supplementary_arg}
        """


rule build_tn5_bins_genrich_virus_count_matrix:
    input:
        counts=expand(
            join(
                RESULTSDIR,
                "{sample}",
                "tn5_motif",
                "genrich",
                "{virus}",
                "{sample}.{virus}.genrich.tn5_sites.bin_counts."
                + str(TN5_VIRUS_BIN_SIZE)
                + "bp.tsv",
            ),
            sample=SAMPLES,
        ),
    output:
        tsv=join(
            RESULTSDIR,
            "tn5_motif",
            "{virus}.genrich.tn5_bin_count_matrix."
            + str(TN5_VIRUS_BIN_SIZE)
            + "bp.tsv",
        ),
    threads:
        _get_threads("build_tn5_bins_genrich_virus_count_matrix", profile_config)
    container:
        config["containers"]["py311"]
    log:
        join(_logdir("build_tn5_bins_genrich_virus_count_matrix"), "{virus}.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        mkdir -p $(dirname {output.tsv})
        python {SCRIPTS_DIR}/build_tn5_count_matrix.py \
          --counts {input.counts} \
          --output {output.tsv}
        """


rule build_tn5_bins_macs2_virus_count_matrix:
    input:
        counts=expand(
            join(
                RESULTSDIR,
                "{sample}",
                "tn5_motif",
                "macs2",
                "{virus}",
                "{sample}.{virus}.macs2.tn5_sites.bin_counts."
                + str(TN5_VIRUS_BIN_SIZE)
                + "bp.tsv",
            ),
            sample=SAMPLES,
        ),
    output:
        tsv=join(
            RESULTSDIR,
            "tn5_motif",
            "{virus}.macs2.tn5_bin_count_matrix."
            + str(TN5_VIRUS_BIN_SIZE)
            + "bp.tsv",
        ),
    threads:
        _get_threads("build_tn5_bins_macs2_virus_count_matrix", profile_config)
    container:
        config["containers"]["py311"]
    log:
        join(_logdir("build_tn5_bins_macs2_virus_count_matrix"), "{virus}.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        mkdir -p $(dirname {output.tsv})
        python {SCRIPTS_DIR}/build_tn5_count_matrix.py \
          --counts {input.counts} \
          --output {output.tsv}
        """
