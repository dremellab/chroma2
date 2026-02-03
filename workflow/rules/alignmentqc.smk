###################################################################################
# alignment qc summary
###################################################################################

rule alignment_flagstat_summary:
    input:
        aligned=expand(
            join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.flagstat.txt"),
            sample=SAMPLES,
        ),
        clean=expand(
            join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.clean.flagstat.txt"),
            sample=SAMPLES,
        ),
        fixmate=expand(
            join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.fixmate.flagstat.txt"),
            sample=SAMPLES,
        ),
        dedup=expand(
            join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.dedup.flagstat.txt"),
            sample=SAMPLES,
        ),
        final=expand(
            join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.final.flagstat.txt"),
            sample=SAMPLES,
        ),
    output:
        tsv=join(RESULTSDIR, "alignmentqc", "flagstat_summary.tsv"),
    params:
        outdir=join(RESULTSDIR, "alignmentqc"),
    threads:
        _get_threads("ataqv_sample", profile_config)
    container:
        config["containers"]["py311"]
    log:
        join(_logdir("alignment_flagstat_summary"), "alignment_flagstat_summary.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        python {SCRIPTS_DIR}/summarize_flagstat.py \
          --results-dir {RESULTSDIR} \
          --output {output.tsv}
        """


###################################################################################
# ATAC-seq QC with ataqv
###################################################################################

rule ataqv:
    input:
        bam=join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.clean.sorted.bam"),
        bai=join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.clean.sorted.bam.bai"),
        regions=REF_REGIONS_HOST,
        peaks=join(RESULTSDIR, "{sample}", "peaks", "{sample}.host.macs2_peaks.narrowPeak"),
        tss=join(REF_DIR, "ref.tss.bed"),
    output:
        json=join(RESULTSDIR, "{sample}", "alignmentqc", "ataqv", "{sample}.json"),
    params:
        outdir=join(RESULTSDIR, "{sample}", "alignmentqc", "ataqv"),
        tmpdir=join(TEMPDIR, "ataqv", "{sample}"),
        extra_args=config.get("ataqv", {}).get("extra_args", ""),
    threads:
        _get_threads("ataqv", profile_config)
    container:
        config["containers"]["ataqv"]
    log:
        join(_logdir("ataqv"), "{sample}.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        echo "[ataqv] sample={wildcards.sample} bam={input.bam}"
        mkdir -p {params.outdir}
        mkdir -p {params.tmpdir}
        host_regions=$(awk -F '\t' '{{n=split($2,a," "); for(i=1;i<=n;i++) if (a[i]!="") print a[i]}}' {input.regions} | sort -u | tr '\n' ' ')
        echo "[ataqv] host_regions=$host_regions"
        host_bam={params.tmpdir}/{wildcards.sample}.host.bam
        if [ ! -f {input.bam}.bai ] || [ {input.bam} -nt {input.bam}.bai ]; then
            echo "[ataqv] indexing input bam"
            samtools index -@ {threads} {input.bam}
        fi
        echo "[ataqv] extracting host bam"
        samtools view -@ {threads} -b {input.bam} $host_regions -o $host_bam
        samtools index $host_bam
        host_markdup_bam={params.tmpdir}/{wildcards.sample}.host.bam
        host_qname_bam={params.tmpdir}/{wildcards.sample}.host.qname.bam
        host_fixmate_bam={params.tmpdir}/{wildcards.sample}.host.fixmate.bam
        host_fixmate_sorted_bam={params.tmpdir}/{wildcards.sample}.host.fixmate.sorted.bam
        echo "[ataqv] fixmate + markdup on host bam"
        samtools sort -@ {threads} -n -T {params.tmpdir}/host_qname -o $host_qname_bam $host_bam
        samtools fixmate -@ {threads} -m $host_qname_bam $host_fixmate_bam
        samtools sort -@ {threads} -T {params.tmpdir}/host_fixmate -o $host_fixmate_sorted_bam $host_fixmate_bam
        samtools markdup -@ {threads} $host_fixmate_sorted_bam $host_markdup_bam
        samtools index $host_markdup_bam
        species=""
        if [[ "{HOST}" == hg* ]]; then
            species="human"
        elif [[ "{HOST}" == mm* ]]; then
            species="mouse"
        fi
        echo "[ataqv] species=$species tss={input.tss} peaks={input.peaks}"
        cd {params.outdir}
        ataqv {params.extra_args} --threads {threads} --tss-file {input.tss} --peak-file {input.peaks} $species $host_markdup_bam > {log} 2>&1
        mv -f {params.outdir}/{wildcards.sample}.host.bam.ataqv.json {output.json}
        """


###################################################################################
# combine ataqv jsons
###################################################################################

rule ataqv_report:
    input:
        jsons=expand(join(RESULTSDIR, "{sample}", "alignmentqc", "ataqv", "{sample}.json"), sample=SAMPLES),
    output:
        report=directory(join(RESULTSDIR, "alignmentqc", "ataqv", "final_report")),
    params:
        genome=HOST,
    threads: 1
    container:
        config["containers"]["ataqv"]
    log:
        join(_logdir("ataqv_report"), "ataqv_report.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        mkdir -p {output.report}
        echo "[ataqv_report] genome={params.genome}"
        echo "[ataqv_report] output_dir={output.report}"
        echo "[ataqv_report] jsons={input.jsons}"
        mkarv -f {output.report} {input.jsons}
        """
