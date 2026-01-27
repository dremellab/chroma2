###################################################################################
# trimming
###################################################################################

rule fastp:
    input:
        unpack(get_fastqs),
    output:
        R1_trimmed=join(RESULTSDIR, "{sample}", "trim", "{sample}.trimmed_R1.fastq.gz"),
        R2_trimmed=join(RESULTSDIR, "{sample}", "trim", "{sample}.trimmed_R2.fastq.gz"),
        html=join(RESULTSDIR, "{sample}", "trim", "{sample}.fastp_report.html"),
        json=join(RESULTSDIR, "{sample}", "trim", "{sample}.fastp_report.json")
    params:
        sample="{sample}",
        peorse=get_peorse,
        length_required=str(config["fastp"]["length_required"]),
        qualified_quality_phred=str(config["fastp"]["qualified_quality_phred"]),
    threads:
        _get_threads("fastp", profile_config)
    container:
        config["containers"]["fastp"]
    shell:
        r"""
        set -exo pipefail
        mkdir -p "$(dirname {output.R1_trimmed})"
        if [ "{params.peorse}" == "PE" ];then
          fastp \
            -i {input.R1} \
            -I {input.R2} \
            -o {output.R1_trimmed} \
            -O {output.R2_trimmed} \
            --length_required {params.length_required} \
            --qualified_quality_phred {params.qualified_quality_phred} \
            --thread {threads} \
            --detect_adapter_for_pe \
            --html {output.html} \
            --json {output.json}
        fi
        if [ "{params.peorse}" == "SE" ];then
          fastp \
            -i {input.R1} \
            -o {output.R1_trimmed} \
            --length_required {params.length_required} \
            --qualified_quality_phred {params.qualified_quality_phred} \
            --thread {threads} \
            --html {output.html} \
            --json {output.json}
          touch {output.R2_trimmed}
        fi
        """


###################################################################################
# bowtie2 index
###################################################################################

rule bowtie2_index:
    input:
        fasta=join(REF_DIR, "ref.fa")
    output:
        expand(join(REF_DIR, "ref.fa.{ext}.bt2"), ext=["1", "2", "3", "4", "rev.1", "rev.2"])
    params:
        prefix=join(REF_DIR, "ref.fa")
    threads:
        _get_threads("bowtie2_index", profile_config)
    container:
        config["containers"]["bowtie2"]
    shell:
        r"""
        set -exo pipefail
        bowtie2-build --threads {threads} {input.fasta} {params.prefix}
        """
