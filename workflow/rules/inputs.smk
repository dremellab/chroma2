###################################################################################
# pooled input controls (host + per-virus)
###################################################################################

rule pool_host_qname_bam:
    input:
        bams=lambda wc: [
            join(RESULTSDIR, s, "align", f"{s}.aligned.host.qname.bam")
            for s in HOST_POOL_CONTROLS.get(wc.pool, [])
        ],
    output:
        bam=join(RESULTSDIR, "inputs", "host", "{pool}.qname.bam"),
    params:
        outdir=join(RESULTSDIR, "inputs", "host"),
    threads:
        _get_threads("pool_host_qname_bam", profile_config),
    container:
        config["containers"]["bowtie2"],
    log:
        join(_logdir("pool_host_qname_bam"), "{pool}.log"),
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        mkdir -p {params.outdir}
        samtools merge -f {output.bam} {input.bams}
        """


rule pool_host_filtered_bam:
    input:
        bams=lambda wc: [
            join(RESULTSDIR, s, "postprocess", f"{s}.host.bam")
            for s in HOST_POOL_CONTROLS.get(wc.pool, [])
        ],
    output:
        bam=join(RESULTSDIR, "inputs", "host", "{pool}.filtered.bam"),
        bai=join(RESULTSDIR, "inputs", "host", "{pool}.filtered.bam.bai"),
    params:
        outdir=join(RESULTSDIR, "inputs", "host"),
    threads:
        _get_threads("pool_host_filtered_bam", profile_config),
    container:
        config["containers"]["bowtie2"],
    log:
        join(_logdir("pool_host_filtered_bam"), "{pool}.log"),
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        mkdir -p {params.outdir}
        samtools merge -f {output.bam} {input.bams}
        samtools index {output.bam} {output.bai}
        """


rule pool_virus_qname_bam:
    input:
        bams=lambda wc: [
            join(RESULTSDIR, s, "align", f"{s}.aligned.virus.{wc.virus}.qname.bam")
            for s in VIRUS_POOL_CONTROLS.get(wc.pool, [])
        ],
    output:
        bam=join(RESULTSDIR, "inputs", "virus", "{pool}", "{virus}.qname.bam"),
    params:
        outdir=join(RESULTSDIR, "inputs", "virus", "{pool}"),
    threads:
        _get_threads("pool_virus_qname_bam", profile_config),
    container:
        config["containers"]["bowtie2"],
    log:
        join(_logdir("pool_virus_qname_bam"), "{pool}.{virus}.log"),
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        mkdir -p {params.outdir}
        samtools merge -f {output.bam} {input.bams}
        """


rule pool_virus_filtered_bam:
    input:
        bams=lambda wc: [
            join(RESULTSDIR, s, "postprocess", f"{s}.virus.{wc.virus}.bam")
            for s in VIRUS_POOL_CONTROLS.get(wc.pool, [])
        ],
    output:
        bam=join(RESULTSDIR, "inputs", "virus", "{pool}", "{virus}.filtered.bam"),
        bai=join(RESULTSDIR, "inputs", "virus", "{pool}", "{virus}.filtered.bam.bai"),
    params:
        outdir=join(RESULTSDIR, "inputs", "virus", "{pool}"),
    threads:
        _get_threads("pool_virus_filtered_bam", profile_config),
    container:
        config["containers"]["bowtie2"],
    log:
        join(_logdir("pool_virus_filtered_bam"), "{pool}.{virus}.log"),
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        mkdir -p {params.outdir}
        samtools merge -f {output.bam} {input.bams}
        samtools index {output.bam} {output.bai}
        """
