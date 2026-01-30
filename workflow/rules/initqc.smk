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
    log:
        join(_logdir("fastqc_raw"), "{sample}.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        echo "[fastqc_raw] sample={params.sample} peorse={params.peorse}"
        echo "[fastqc_raw] R1={input.R1}"
        echo "[fastqc_raw] R2={input.R2}"
        echo "[fastqc_raw] outdir={params.outdir} threads={threads} mem={params.memoryG}"
        mkdir -p {params.outdir}
        r1_base=$(basename {input.R1})
        r1_name="${{r1_base%.fastq.gz}}"
        r1_name="${{r1_name%.fq.gz}}"
        r1_name="${{r1_name%.fastq}}"
        r1_name="${{r1_name%.fq}}"
        if [ "{params.peorse}" == "PE" ]; then
            r2_base=$(basename {input.R2})
            r2_name="${{r2_base%.fastq.gz}}"
            r2_name="${{r2_name%.fq.gz}}"
            r2_name="${{r2_name%.fastq}}"
            r2_name="${{r2_name%.fq}}"
            fastqc {input.R1} {input.R2} \
                --threads {threads} \
                --outdir {params.outdir} \
                --quiet
            mv "{params.outdir}/${{r1_name}}_fastqc.html" {output.R1_fastqc}
            mv "{params.outdir}/${{r1_name}}_fastqc.zip" {output.R1_fastqc_zip}
            mv "{params.outdir}/${{r2_name}}_fastqc.html" {output.R2_fastqc}
            mv "{params.outdir}/${{r2_name}}_fastqc.zip" {output.R2_fastqc_zip}
        else
            fastqc {input.R1} \
                --threads {threads} \
                --outdir {params.outdir} \
                --quiet
            mv "{params.outdir}/${{r1_name}}_fastqc.html" {output.R1_fastqc}
            mv "{params.outdir}/${{r1_name}}_fastqc.zip" {output.R1_fastqc_zip}
            touch {output.R2_fastqc}
            touch {output.R2_fastqc_zip}
        fi
        echo "[fastqc_raw] done"
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
    log:
        join(_logdir("fastqc_trimmed"), "{sample}.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        echo "[fastqc_trimmed] sample={params.sample} peorse={params.peorse}"
        echo "[fastqc_trimmed] R1={input.R1}"
        echo "[fastqc_trimmed] R2={input.R2}"
        echo "[fastqc_trimmed] outdir={params.outdir} threads={threads} mem={params.memoryG}"
        mkdir -p {params.outdir}
        r1_base=$(basename {input.R1})
        r1_name="${{r1_base%.fastq.gz}}"
        r1_name="${{r1_name%.fq.gz}}"
        r1_name="${{r1_name%.fastq}}"
        r1_name="${{r1_name%.fq}}"
        if [ "{params.peorse}" == "PE" ]; then
            r2_base=$(basename {input.R2})
            r2_name="${{r2_base%.fastq.gz}}"
            r2_name="${{r2_name%.fq.gz}}"
            r2_name="${{r2_name%.fastq}}"
            r2_name="${{r2_name%.fq}}"
            fastqc {input.R1} {input.R2} \
                --threads {threads} \
                --outdir {params.outdir} \
                --quiet
            # mv "{params.outdir}/${{r1_name}}_fastqc.html" {output.R1_fastqc}
            # mv "{params.outdir}/${{r1_name}}_fastqc.zip" {output.R1_fastqc_zip}
            # mv "{params.outdir}/${{r2_name}}_fastqc.html" {output.R2_fastqc}
            # mv "{params.outdir}/${{r2_name}}_fastqc.zip" {output.R2_fastqc_zip}
        else
            fastqc {input.R1} \
                --threads {threads} \
                --outdir {params.outdir} \
                --quiet
            # mv "{params.outdir}/${{r1_name}}_fastqc.html" {output.R1_fastqc}
            # mv "{params.outdir}/${{r1_name}}_fastqc.zip" {output.R1_fastqc_zip}
            touch {output.R2_fastqc}
            touch {output.R2_fastqc_zip}
        fi
        echo "[fastqc_trimmed] done"
        """
