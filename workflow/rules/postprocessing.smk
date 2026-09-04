###################################################################################
# split final BAM into host + per-virus BAMs, then generate bigWigs
###################################################################################

rule split_host_bam:
    input:
        bam=join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.final.bam"),
        bai=join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.final.bam.bai"),
        regions=REF_REGIONS_HOST,
    output:
        bam=join(RESULTSDIR, "{sample}", "postprocess", "{sample}.host.bam"),
        bai=join(RESULTSDIR, "{sample}", "postprocess", "{sample}.host.bam.bai"),
    params:
        outdir=join(RESULTSDIR, "{sample}", "postprocess"),
    threads:
        _get_threads("split_host_bam", profile_config)
    resources:
        runtime=lambda wildcards, attempt: _get_runtime_with_retries(
            "split_host_bam", profile_config, attempt
        ),
    container:
        config["containers"]["bowtie2"]
    log:
        join(_logdir("split_host_bam"), "{sample}.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        mkdir -p {params.outdir}
        host_regions=$(awk -F '\t' '{{n=split($2,a," "); for(i=1;i<=n;i++) if (a[i]!="") print a[i]}}' {input.regions} | sort -u | tr '\n' ' ')
        echo "[split_host_bam] sample={wildcards.sample} regions_file={input.regions}"
        echo "[split_host_bam] contigs: $host_regions"
        if [ ! -f {input.bam}.bai ] || [ {input.bam} -nt {input.bam}.bai ]; then
            samtools index -@ {threads} {input.bam}
        fi
        samtools view -@ {threads} -b {input.bam} $host_regions -o {output.bam}
        samtools index {output.bam} {output.bai}
        """


rule split_virus_bam:
    input:
        bam=join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.final.bam"),
        bai=join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.final.bam.bai"),
        regions=join(FASTAS_GTFS_DIR, "{virus}.fa.regions"),
    output:
        bam=join(RESULTSDIR, "{sample}", "postprocess", "{sample}.virus.{virus}.bam"),
        bai=join(RESULTSDIR, "{sample}", "postprocess", "{sample}.virus.{virus}.bam.bai"),
    params:
        outdir=join(RESULTSDIR, "{sample}", "postprocess"),
    threads:
        _get_threads("split_virus_bam", profile_config)
    resources:
        runtime=lambda wildcards, attempt: _get_runtime_with_retries(
            "split_virus_bam", profile_config, attempt
        ),
    container:
        config["containers"]["bowtie2"]
    log:
        join(_logdir("split_virus_bam"), "{sample}.{virus}.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        mkdir -p {params.outdir}
        virus_regions=$(awk -F '\t' '{{n=split($2,a," "); for(i=1;i<=n;i++) if (a[i]!="") print a[i]}}' {input.regions} | sort -u | tr '\n' ' ')
        echo "[split_virus_bam] sample={wildcards.sample} virus={wildcards.virus} regions_file={input.regions}"
        echo "[split_virus_bam] contigs: $virus_regions"
        if [ ! -f {input.bam}.bai ] || [ {input.bam} -nt {input.bam}.bai ]; then
            samtools index -@ {threads} {input.bam}
        fi
        samtools view -@ {threads} -b {input.bam} $virus_regions -o {output.bam}
        samtools index {output.bam} {output.bai}
        """


# Whether to scale bigwigs by the precomputed read-depth normalization factor
# (results/alignmentqc/norm_factors.tsv, see plans/read-depth-normalization-plan.md).
# Default on; set config normalization.use_read_depth_scaling: false to restore
# unscaled bamCoverage output.
USE_READ_DEPTH_SCALING = config.get("normalization", {}).get("use_read_depth_scaling", True)


rule bamcoverage_host:
    input:
        bam=join(RESULTSDIR, "{sample}", "postprocess", "{sample}.host.bam"),
        bai=join(RESULTSDIR, "{sample}", "postprocess", "{sample}.host.bam.bai"),
        norm_factors=join(RESULTSDIR, "alignmentqc", "norm_factors.tsv"),
    output:
        bw=join(RESULTSDIR, "{sample}", "bigwig", "{sample}.host.bw"),
    params:
        outdir=join(RESULTSDIR, "{sample}", "bigwig"),
        binsize=str(config.get("postprocessing", {}).get("host_bin_size", 50)),
        extra_args=config.get("postprocessing", {}).get("host_bw_extra_args", ""),
        use_scaling="1" if USE_READ_DEPTH_SCALING else "0",
    threads:
        _get_threads("bamcoverage_host", profile_config)
    resources:
        runtime=lambda wildcards, attempt: _get_runtime_with_retries(
            "bamcoverage_host", profile_config, attempt
        ),
    container:
        config["containers"]["deeptools"]
    log:
        join(_logdir("bamcoverage_host"), "{sample}.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        mkdir -p {params.outdir}
        scale_args=$(bash {SCRIPTS_DIR}/get_bw_scale_args.sh bamcoverage_host {wildcards.sample} {params.use_scaling} {input.norm_factors})
        bamCoverage \
          -b {input.bam} \
          -o {output.bw} \
          --binSize {params.binsize} \
          --numberOfProcessors {threads} $scale_args {params.extra_args}
        """


rule bamcoverage_virus:
    input:
        bam=join(RESULTSDIR, "{sample}", "postprocess", "{sample}.virus.{virus}.bam"),
        bai=join(RESULTSDIR, "{sample}", "postprocess", "{sample}.virus.{virus}.bam.bai"),
        norm_factors=join(RESULTSDIR, "alignmentqc", "norm_factors.tsv"),
    output:
        bw=join(RESULTSDIR, "{sample}", "bigwig", "{sample}.virus.{virus}.bw"),
    params:
        outdir=join(RESULTSDIR, "{sample}", "bigwig"),
        binsize=str(config.get("postprocessing", {}).get("virus_bin_size", 1)),
        extra_args=config.get("postprocessing", {}).get("virus_bw_extra_args", ""),
        use_scaling="1" if USE_READ_DEPTH_SCALING else "0",
    threads:
        _get_threads("bamcoverage_virus", profile_config)
    resources:
        runtime=lambda wildcards, attempt: _get_runtime_with_retries(
            "bamcoverage_virus", profile_config, attempt
        ),
    container:
        config["containers"]["deeptools"]
    log:
        join(_logdir("bamcoverage_virus"), "{sample}.{virus}.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        mkdir -p {params.outdir}
        scale_args=$(bash {SCRIPTS_DIR}/get_bw_scale_args.sh bamcoverage_virus {wildcards.sample} {params.use_scaling} {input.norm_factors})
        bamCoverage \
          -b {input.bam} \
          -o {output.bw} \
          --binSize {params.binsize} \
          --numberOfProcessors {threads} $scale_args {params.extra_args}
        """
