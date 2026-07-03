from os.path import join

MULTIQC_REPORT = join(RESULTSDIR, "multiqc", "multiqc_report.html")
MULTIQC_DATA_DIR = join(RESULTSDIR, "multiqc", "multiqc_data")
MULTIQC_CUSTOM = join(RESULTSDIR, "multiqc_extra_data", "custom")
MULTIQC_ATAQV = join(RESULTSDIR, "multiqc_extra_data", "ataqv")

MAPQ_THRESHOLDS = [0, 10, 20, 30]

rule ataqv_symlink_mqc:
    input:
        host=expand(
            join(RESULTSDIR, "{sample}", "alignmentqc", "ataqv", "{sample}.host.json"),
            sample=CASE_SAMPLES,
        ),
    output:
        expand(
            join(MULTIQC_ATAQV, "{sample}.host.ataqv.json"),
            sample=CASE_SAMPLES,
        ),
    run:
        import os
        for s in CASE_SAMPLES:
            src = join(RESULTSDIR, s, "alignmentqc", "ataqv", f"{s}.host.json")
            dst = join(MULTIQC_ATAQV, f"{s}.host.ataqv.json")
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            if os.path.exists(dst) or os.path.islink(dst):
                os.remove(dst)
            os.symlink(os.path.abspath(src), dst)

rule fragment_size:
    input:
        bam=join(RESULTSDIR, "{sample}", "postprocess", "{sample}.host.bam"),
        bai=join(RESULTSDIR, "{sample}", "postprocess", "{sample}.host.bam.bai"),
    output:
        hist=join(RESULTSDIR, "{sample}", "multiqc", "{sample}.host.bamPEFragmentSize.png"),
    threads:
        _get_threads("fragment_size", profile_config)
    container:
        config["containers"]["deeptools"]
    log:
        join(_logdir("fragment_size"), "{sample}.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log}) $(dirname {output.hist})
        exec > >(tee -a {log}) 2>&1
        bamPEFragmentSize \
          --bamfiles {input.bam} \
          --histogram {output.hist} \
          --numberOfProcessors {threads} \
          --samplesLabel {wildcards.sample}
        """

rule extract_fragment_size_data:
    input:
        bam=join(RESULTSDIR, "{sample}", "postprocess", "{sample}.host.bam"),
    output:
        data=join(RESULTSDIR, "{sample}", "multiqc", "{sample}.host.fragment_sizes.tsv"),
    container:
        config["containers"]["pysam"]
    log:
        join(_logdir("extract_fragment_size_data"), "{sample}.log")
    script:
        "../scripts/extract_fragment_sizes.py"

rule fragment_size_aggregate_mqc:
    input:
        expand(
            join(RESULTSDIR, "{sample}", "multiqc", "{sample}.host.fragment_sizes.tsv"),
            sample=CASE_SAMPLES,
        ),
    output:
        join(MULTIQC_CUSTOM, "fragment_size_mqc.tsv"),
    container:
        config["containers"]["pysam"]
    log:
        join(_logdir("fragment_size_aggregate_mqc"), "fragment_size_aggregate_mqc.log")
    script:
        "../scripts/aggregate_fragment_sizes_mqc.py"

rule frip_featurecounts:
    input:
        bam=join(RESULTSDIR, "{sample}", "postprocess", "{sample}.host.bam"),
        bai=join(RESULTSDIR, "{sample}", "postprocess", "{sample}.host.bam.bai"),
        peaks=join(RESULTSDIR, "{sample}", "peaks", "{sample}.host.macs2_peaks.narrowPeak.gz"),
    output:
        counts=join(RESULTSDIR, "{sample}", "multiqc", "{sample}.host.frip.featureCounts"),
        summary=join(RESULTSDIR, "{sample}", "multiqc", "{sample}.host.frip.featureCounts.summary"),
    threads:
        _get_threads("frip_featurecounts", profile_config)
    container:
        config["containers"]["featurecounts"]
    log:
        join(_logdir("frip_featurecounts"), "{sample}.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log}) $(dirname {output.counts})
        exec > >(tee -a {log}) 2>&1
        zcat {input.peaks} | awk 'BEGIN{{OFS="\t"}} {{print $1"_"$2"_"$3, $1, $2+1, $3, "."}}' \
          > {output.counts}.saf
        featureCounts \
          -F SAF \
          -a {output.counts}.saf \
          -o {output.counts} \
          -p --countReadPairs \
          -T {threads} \
          {input.bam}
        """

rule alignment_stats_mqc:
    input:
        summary=join(RESULTSDIR, "alignmentqc", "idxstats_summary.tsv"),
    output:
        stats=join(MULTIQC_CUSTOM, "alignment_stats_mqc.tsv"),
        ratio=join(MULTIQC_CUSTOM, "host_virus_ratio_mqc.tsv"),
    container:
        config["containers"]["pysam"]
    log:
        join(_logdir("alignment_stats_mqc"), "alignment_stats_mqc.log")
    script:
        "../scripts/make_alignment_stats_mqc.py"

rule tn5_counts_mqc:
    input:
        COUNT_MATRIX_ALL_OUTPUTS,
    output:
        join(MULTIQC_CUSTOM, "tn5_counts_mqc.tsv"),
    container:
        config["containers"]["pysam"]
    log:
        join(_logdir("tn5_counts_mqc"), "tn5_counts_mqc.log")
    script:
        "../scripts/make_tn5_counts_mqc.py"

rule peak_size_mqc:
    input:
        expand(
            join(RESULTSDIR, "{sample}", "peaks", "{sample}.host.macs2_peaks.narrowPeak.gz"),
            sample=CASE_SAMPLES,
        ),
    output:
        join(MULTIQC_CUSTOM, "peak_size_mqc.tsv"),
    container:
        config["containers"]["pysam"]
    log:
        join(_logdir("peak_size_mqc"), "peak_size_mqc.log")
    script:
        "../scripts/make_peak_size_mqc.py"

rule genome_coverage_mqc:
    input:
        bam=join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.final.bam"),
        bai=join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.final.bam.bai"),
    output:
        temp(
            expand(
                join(RESULTSDIR, "{{sample}}", "multiqc", "{{sample}}.coverage_q{mapq}.txt"),
                mapq=MAPQ_THRESHOLDS,
            )
        ),
    params:
        outpfx=lambda wc: join(RESULTSDIR, wc.sample, "multiqc", f"{wc.sample}"),
        mapqs=MAPQ_THRESHOLDS,
    container:
        config["containers"]["py311"]
    log:
        join(_logdir("genome_coverage_mqc"), "{sample}.log")
    script:
        "../scripts/compute_coverage_mqc.py"

rule genome_coverage_aggregate_mqc:
    input:
        expand(
            join(RESULTSDIR, "{sample}", "multiqc", "{sample}.coverage_q{mapq}.txt"),
            sample=CASE_SAMPLES,
            mapq=MAPQ_THRESHOLDS,
        ),
    output:
        join(MULTIQC_CUSTOM, "genome_coverage_mqc.tsv"),
    container:
        config["containers"]["pysam"]
    log:
        join(_logdir("genome_coverage_aggregate_mqc"), "genome_coverage_aggregate_mqc.log")
    script:
        "../scripts/make_genome_coverage_mqc.py"

rule multiqc:
    input:
        expand(join(RESULTSDIR, "{s}", "fastqc", "{s}.raw_R1_fastqc.zip"), s=SAMPLES),
        expand(join(RESULTSDIR, "{s}", "fastqc", "{s}.raw_R2_fastqc.zip"), s=SAMPLES),
        expand(join(RESULTSDIR, "{s}", "fastqc", "{s}.trimmed_R1_fastqc.zip"), s=SAMPLES),
        expand(join(RESULTSDIR, "{s}", "fastqc", "{s}.trimmed_R2_fastqc.zip"), s=SAMPLES),
        expand(join(RESULTSDIR, "{s}", "trim", "{s}.fastp_report.json"), s=SAMPLES),
        expand(join(MULTIQC_ATAQV, "{s}.host.ataqv.json"), s=CASE_SAMPLES),
        expand(join(RESULTSDIR, "{s}", "multiqc", "{s}.host.bamPEFragmentSize.png"), s=CASE_SAMPLES),
        expand(join(RESULTSDIR, "{s}", "multiqc", "{s}.host.frip.featureCounts.summary"), s=CASE_SAMPLES),
        join(MULTIQC_CUSTOM, "alignment_stats_mqc.tsv"),
        join(MULTIQC_CUSTOM, "host_virus_ratio_mqc.tsv"),
        join(MULTIQC_CUSTOM, "tn5_counts_mqc.tsv"),
        join(MULTIQC_CUSTOM, "peak_size_mqc.tsv"),
        join(MULTIQC_CUSTOM, "genome_coverage_mqc.tsv"),
        join(MULTIQC_CUSTOM, "fragment_size_mqc.tsv"),
    output:
        report=MULTIQC_REPORT,
        data_dir=directory(MULTIQC_DATA_DIR),
    params:
        outdir=join(RESULTSDIR, "multiqc"),
        search_dir=RESULTSDIR,
    threads:
        _get_threads("multiqc", profile_config)
    container:
        config["containers"]["multiqc"]
    log:
        join(_logdir("multiqc"), "multiqc.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log}) {params.outdir}
        exec > >(tee -a {log}) 2>&1
        multiqc {params.search_dir} \
          --outdir {params.outdir} \
          --force \
          --no-ansi \
          --ignore logs
        """
