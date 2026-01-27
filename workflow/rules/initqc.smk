###################################################################################
# fastqc
###################################################################################

rule fastqc_raw:
    input:
        unpack(get_fastqs),
    output:
        R1_fastqc=join(RESULTSDIR, "{sample}", "fastqc", "{sample}.raw_R1_fastqc.html"),
        R1_fastqc_zip=join(RESULTSDIR, "{sample}", "fastqc", "{sample}.raw_R1_fastqc.zip"),
        R2_fastqc=join(RESULTSDIR, "{sample}", "fastqc", "{sample}.raw_R2_fastqc.html"),
        R2_fastqc_zip=join(RESULTSDIR, "{sample}", "fastqc", "{sample}.raw_R2_fastqc.zip"),
    params:
        sample="{sample}",
        outdir=join(RESULTSDIR, "{sample}", "fastqc"),
        memoryG=getmemG("fastqc"),
        peorse=get_peorse,
    threads: 1
    container:
        config["containers"]["fastqc"]
    shell:
        r"""
        set -exo pipefail
        mkdir -p {params.outdir}
        if [ "{params.peorse}" == "PE" ]; then
            fastqc {input.R1} {input.R2} \
                --threads {threads} \
                --memory {params.memoryG} \
                --outdir {params.outdir} \
                --quiet
        else
            fastqc {input.R1} \
                --threads {threads} \
                --memory {params.memoryG} \
                --outdir {params.outdir} \
                --quiet
            touch {output.R2_fastqc}
            touch {output.R2_fastqc_zip}
        fi
        """


rule fastqc_trimmed:
    input:
        R1=join(RESULTSDIR, "{sample}", "trim", "{sample}.trimmed_R1.fastq.gz"),
        R2=join(RESULTSDIR, "{sample}", "trim", "{sample}.trimmed_R2.fastq.gz"),
    output:
        R1_fastqc=join(RESULTSDIR, "{sample}", "fastqc", "{sample}.trimmed_R1_fastqc.html"),
        R1_fastqc_zip=join(RESULTSDIR, "{sample}", "fastqc", "{sample}.trimmed_R1_fastqc.zip"),
        R2_fastqc=join(RESULTSDIR, "{sample}", "fastqc", "{sample}.trimmed_R2_fastqc.html"),
        R2_fastqc_zip=join(RESULTSDIR, "{sample}", "fastqc", "{sample}.trimmed_R2_fastqc.zip"),
    params:
        sample="{sample}",
        outdir=join(RESULTSDIR, "{sample}", "fastqc"),
        memoryG=getmemG("fastqc"),
        peorse=get_peorse,
    threads: 1
    container:
        config["containers"]["fastqc"]
    shell:
        r"""
        set -exo pipefail
        mkdir -p {params.outdir}
        if [ "{params.peorse}" == "PE" ]; then
            fastqc {input.R1} {input.R2} \
                --threads {threads} \
                --memory {params.memoryG} \
                --outdir {params.outdir} \
                --quiet
        else
            fastqc {input.R1} \
                --threads {threads} \
                --memory {params.memoryG} \
                --outdir {params.outdir} \
                --quiet
            touch {output.R2_fastqc}
            touch {output.R2_fastqc_zip}
        fi
        """
