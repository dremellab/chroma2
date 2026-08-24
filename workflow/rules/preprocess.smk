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
    resources:
        runtime=lambda wildcards, attempt: _get_runtime_with_retries(
            "fastp", profile_config, attempt
        ),
    container:
        config["containers"]["fastp"]
    log:
        join(_logdir("fastp"), "{sample}.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
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
    log:
        join(_logdir("bowtie2_index"), "bowtie2_index.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        bowtie2-build --threads {threads} {input.fasta} {params.prefix}
        """


###################################################################################
# build TSS bed from GTFs (host + per-virus)
###################################################################################

rule build_ref_tss_host:
    input:
        gtf=REF_GTF,
        regions=REF_REGIONS_HOST,
        fasta=REF_FA,
    output:
        tss_host=join(REF_DIR, "ref.tss.host.bed"),
        chromsizes_host=join(REF_DIR, "ref.chrom.sizes.host.txt"),
    params:
        tmpdir=join(TEMPDIR, "ref_tss", "host"),
    threads:
        _get_threads("build_ref_tss_host", profile_config)
    container:
        config["containers"]["bowtie2"]
    log:
        join(_logdir("build_ref_tss_host"), "build_ref_tss_host.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        mkdir -p {params.tmpdir}
        host_contigs_file={params.tmpdir}/host_contigs.txt
        awk -F '\t' '{{n=split($2,a," "); for(i=1;i<=n;i++) if (a[i]!="") print a[i]}}' {input.regions} | sort -u > $host_contigs_file
        if [[ "{input.gtf}" == *.gz ]]; then
            reader="zcat"
        else
            reader="cat"
        fi
        $reader {input.gtf} \
          | awk '$3 == "transcript"' \
          | awk 'BEGIN{{OFS="\t"}} {{ if($7 == "+") print $1, $4-1, $4; else if($7 == "-") print $1, $5-1, $5; }}' \
          | sort -k1,1 -k2,2n \
          | uniq \
          | awk 'NR==FNR{{keep[$1]=1;next}} ($1 in keep)' $host_contigs_file - \
          > {output.tss_host}
        if [ ! -s {input.fasta}.fai ]; then
            samtools faidx {input.fasta}
        fi
        awk 'NR==FNR{{keep[$1]=1;next}} ($1 in keep){{print $1"\t"$2}}' $host_contigs_file {input.fasta}.fai \
          > {output.chromsizes_host}
        """


rule build_ref_tss_virus:
    input:
        gtf=join(FASTAS_GTFS_DIR, "{virus}.gtf"),
        regions=join(FASTAS_GTFS_DIR, "{virus}.fa.regions"),
        fasta=join(FASTAS_GTFS_DIR, "{virus}.fa"),
    output:
        tss_virus=join(REF_DIR, "ref.tss.{virus}.bed"),
        chromsizes_virus=join(REF_DIR, "ref.chrom.sizes.{virus}.txt"),
        autosomes_virus=join(REF_DIR, "ref.chrom.autosomes.{virus}.txt"),
    params:
        tmpdir=join(TEMPDIR, "ref_tss", "virus", "{virus}"),
    threads:
        _get_threads("build_ref_tss_virus", profile_config)
    container:
        config["containers"]["bowtie2"]
    log:
        join(_logdir("build_ref_tss_virus"), "{virus}.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        mkdir -p {params.tmpdir}
        virus_contigs_file={params.tmpdir}/virus_contigs.txt
        awk -F '\t' '{{n=split($2,a," "); for(i=1;i<=n;i++) if (a[i]!="") print a[i]}}' {input.regions} | sort -u > $virus_contigs_file
        if [[ "{input.gtf}" == *.gz ]]; then
            reader="zcat"
        else
            reader="cat"
        fi
        $reader {input.gtf} \
          | awk '$3 == "transcript"' \
          | awk 'BEGIN{{OFS="\t"}} {{ if($7 == "+") print $1, $4-1, $4; else if($7 == "-") print $1, $5-1, $5; }}' \
          | sort -k1,1 -k2,2n \
          | uniq \
          | awk 'NR==FNR{{keep[$1]=1;next}} ($1 in keep)' $virus_contigs_file - \
          > {output.tss_virus}
        if [ ! -s {input.fasta}.fai ]; then
            samtools faidx {input.fasta}
        fi
        awk '{{print $1"\t"$2}}' {input.fasta}.fai \
          > {output.chromsizes_virus}
        awk '{{print $1}}' {input.fasta}.fai \
          > {output.autosomes_virus}
        """
