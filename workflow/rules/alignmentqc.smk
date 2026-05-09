###################################################################################
# alignment qc summary
###################################################################################

rule alignment_flagstat_summary:
    input:
        aligned=expand(
            join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.idxstats.txt"),
            sample=SAMPLES,
        ),
        clean=expand(
            join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.clean.idxstats.txt"),
            sample=SAMPLES,
        ),
        dedup=expand(
            join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.dedup.idxstats.txt"),
            sample=SAMPLES,
        ),
        final=expand(
            join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.final.idxstats.txt"),
            sample=SAMPLES,
        ),
        raw_fastqc_r1=expand(
            join(RESULTSDIR, "{sample}", "fastqc", "{sample}.raw_R1_fastqc.zip"),
            sample=SAMPLES,
        ),
        trimmed_fastqc_r1=expand(
            join(RESULTSDIR, "{sample}", "fastqc", "{sample}.trimmed_R1_fastqc.zip"),
            sample=SAMPLES,
        ),
        regions_host=join(REF_DIR, "ref.fa.regions.host"),
        regions_viruses=join(REF_DIR, "ref.fa.regions.viruses"),
    output:
        tsv=join(RESULTSDIR, "alignmentqc", "idxstats_summary.tsv"),
    params:
        outdir=join(RESULTSDIR, "alignmentqc"),
    threads:
        _get_threads("alignment_flagstat_summary", profile_config)
    container:
        config["containers"]["py311"]
    log:
        join(_logdir("alignment_flagstat_summary"), "alignment_flagstat_summary.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        mkdir -p {params.outdir}
        python {SCRIPTS_DIR}/summarize_idxstats.py \
          --results-dir {RESULTSDIR} \
          --regions-host {input.regions_host} \
          --regions-viruses {input.regions_viruses} \
          --output {output.tsv}
        """


###################################################################################
# ATAC-seq QC with ataqv (host + per-virus, custom chromsizes)
###################################################################################

rule ataqv_host:
    input:
        bam=join(RESULTSDIR, "{sample}", "postprocess", "{sample}.host.bam"),
        bai=join(RESULTSDIR, "{sample}", "postprocess", "{sample}.host.bam.bai"),
        peaks=join(RESULTSDIR, "{sample}", "peaks", "{sample}.host.macs2_peaks.narrowPeak.gz"),
        tss=join(REF_DIR, "ref.tss.host.bed"),
        chromsizes=join(REF_DIR, "ref.chrom.sizes.host.txt"),
    output:
        json=join(RESULTSDIR, "{sample}", "alignmentqc", "ataqv", "{sample}.host.json"),
    params:
        outdir=join(RESULTSDIR, "{sample}", "alignmentqc", "ataqv"),
        tmpdir=join(TEMPDIR, "ataqv", "{sample}", "host"),
        extra_args=config.get("ataqv", {}).get("extra_args", ""),
    threads:
        _get_threads("ataqv_host", profile_config)
    container:
        config["containers"]["ataqv"]
    log:
        join(_logdir("ataqv_host"), "{sample}.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        echo "[ataqv_host] sample={wildcards.sample} bam={input.bam}"
        mkdir -p {params.outdir}
        mkdir -p {params.tmpdir}
        host_markdup_bam={params.tmpdir}/{wildcards.sample}.host.bam
        host_qname_bam={params.tmpdir}/{wildcards.sample}.host.qname.bam
        host_fixmate_bam={params.tmpdir}/{wildcards.sample}.host.fixmate.bam
        host_fixmate_sorted_bam={params.tmpdir}/{wildcards.sample}.host.fixmate.sorted.bam
        echo "[ataqv_host] fixmate + markdup on host bam"
        samtools sort -@ {threads} -n -T {params.tmpdir}/host_qname -o $host_qname_bam {input.bam}
        echo "[ataqv_host] fixmate"
        samtools fixmate -@ {threads} -m $host_qname_bam $host_fixmate_bam
        echo "[ataqv_host] coordinate sort"
        samtools sort -@ {threads} -T {params.tmpdir}/host_fixmate -o $host_fixmate_sorted_bam $host_fixmate_bam
        echo "[ataqv_host] markdup"
        samtools markdup -@ {threads} $host_fixmate_sorted_bam $host_markdup_bam
        echo "[ataqv_host] index"
        samtools index $host_markdup_bam
        echo "[ataqv_host] species=custom tss={input.tss} peaks={input.peaks} chromsizes={input.chromsizes}"
        cd {params.outdir}
        echo "[ataqv_host] running ataqv"
        ataqv {params.extra_args} --threads {threads} --tss-file {input.tss} --peak-file <(zcat {input.peaks}) \
          --autosomal-reference-file {input.chromsizes} \
          --mitochondrial-reference-name none \
          --name {wildcards.sample} \
          --metrics-file {output.json} \
          custom $host_markdup_bam 2>&1 | tee -a {log}
        echo "[ataqv_host] done"
        """


rule ataqv_virus:
    input:
        bam=join(RESULTSDIR, "{sample}", "postprocess", "{sample}.virus.{virus}.bam"),
        bai=join(RESULTSDIR, "{sample}", "postprocess", "{sample}.virus.{virus}.bam.bai"),
        tss=join(REF_DIR, "ref.tss.{virus}.bed"),
        chromsizes=join(REF_DIR, "ref.chrom.sizes.{virus}.txt"),
        autosomes=join(REF_DIR, "ref.chrom.autosomes.{virus}.txt"),
    output:
        json=join(RESULTSDIR, "{sample}", "alignmentqc", "ataqv", "{sample}.virus.{virus}.json"),
    params:
        outdir=join(RESULTSDIR, "{sample}", "alignmentqc", "ataqv"),
        tmpdir=join(TEMPDIR, "ataqv", "{sample}", "virus", "{virus}"),
        extra_args=config.get("ataqv", {}).get("extra_args", ""),
        tss_extension=str(config.get("ataqv", {}).get("virus_tss_extension", 200)),
    threads:
        _get_threads("ataqv_virus", profile_config)
    container:
        config["containers"]["ataqv"]
    log:
        join(_logdir("ataqv_virus"), "{sample}.{virus}.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        echo "[ataqv_virus] sample={wildcards.sample} virus={wildcards.virus} bam={input.bam}"
        mkdir -p {params.outdir}
        mkdir -p {params.tmpdir}
        virus_markdup_bam={params.tmpdir}/{wildcards.sample}.virus.{wildcards.virus}.bam
        virus_qname_bam={params.tmpdir}/{wildcards.sample}.virus.{wildcards.virus}.qname.bam
        virus_fixmate_bam={params.tmpdir}/{wildcards.sample}.virus.{wildcards.virus}.fixmate.bam
        virus_fixmate_sorted_bam={params.tmpdir}/{wildcards.sample}.virus.{wildcards.virus}.fixmate.sorted.bam
        echo "[ataqv_virus] fixmate + markdup on virus bam"
        samtools sort -@ {threads} -n -T {params.tmpdir}/virus_qname -o $virus_qname_bam {input.bam}
        echo "[ataqv_virus] fixmate"
        samtools fixmate -@ {threads} -m $virus_qname_bam $virus_fixmate_bam
        echo "[ataqv_virus] coordinate sort"
        samtools sort -@ {threads} -T {params.tmpdir}/virus_fixmate -o $virus_fixmate_sorted_bam $virus_fixmate_bam
        echo "[ataqv_virus] markdup"
        samtools markdup -@ {threads} $virus_fixmate_sorted_bam $virus_markdup_bam
        echo "[ataqv_virus] index"
        samtools index $virus_markdup_bam
        echo "[ataqv_virus] species=custom tss={input.tss} chromsizes={input.chromsizes} autosomes={input.autosomes}"
        cd {params.outdir}
        echo "[ataqv_virus] running ataqv"
        ataqv {params.extra_args} --threads {threads} --tss-file {input.tss} \
          --tss-extension {params.tss_extension} \
          --autosomal-reference-file {input.autosomes} \
          --mitochondrial-reference-name none \
          --name {wildcards.sample} \
          --metrics-file {output.json} \
          custom $virus_markdup_bam 2>&1 | tee -a {log}
        echo "[ataqv_virus] done"
        """


###################################################################################
# combine ataqv jsons
###################################################################################

rule ataqv_report_host:
    input:
        jsons=expand(join(RESULTSDIR, "{sample}", "alignmentqc", "ataqv", "{sample}.host.json"), sample=CASE_SAMPLES),
    output:
        report=directory(join(RESULTSDIR, "alignmentqc", "ataqv", "final_report.host")),
    params:
        genome="custom",
    threads:
        _get_threads("ataqv_report_host", profile_config)
    container:
        config["containers"]["ataqv"]
    log:
        join(_logdir("ataqv_report_host"), "ataqv_report_host.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        mkdir -p {output.report}
        echo "[ataqv_report_host] genome={params.genome}"
        echo "[ataqv_report_host] output_dir={output.report}"
        echo "[ataqv_report_host] jsons={input.jsons}"
        mkarv -f {output.report} {input.jsons}
        """


rule ataqv_report_virus:
    input:
        jsons=lambda wc: expand(
            join(RESULTSDIR, "{sample}", "alignmentqc", "ataqv", "{sample}.virus.{virus}.json"),
            sample=CASE_SAMPLES,
            virus=wc.virus,
        ),
    output:
        report=directory(join(RESULTSDIR, "alignmentqc", "ataqv", "final_report.virus.{virus}")),
    params:
        genome="custom",
    threads:
        _get_threads("ataqv_report_virus", profile_config)
    container:
        config["containers"]["ataqv"]
    log:
        join(_logdir("ataqv_report_virus"), "{virus}.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        mkdir -p {output.report}
        echo "[ataqv_report_virus] genome={params.genome}"
        echo "[ataqv_report_virus] output_dir={output.report}"
        echo "[ataqv_report_virus] jsons={input.jsons}"
        mkarv -f {output.report} {input.jsons}
        """
