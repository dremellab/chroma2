###################################################################################
# MACS2 peak calling (ATAC-seq)
###################################################################################

rule macs2_atac_callpeak:
    input:
        bam=join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.final.bam"),
        bai=join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.final.bam.bai"),
    output:
        narrowpeak=join(RESULTSDIR, "{sample}", "peaks", "{sample}.macs2_peaks.narrowPeak"),
        summits=join(RESULTSDIR, "{sample}", "peaks", "{sample}.macs2_summits.bed"),
        xls=join(RESULTSDIR, "{sample}", "peaks", "{sample}.macs2_peaks.xls"),
    params:
        outdir=join(RESULTSDIR, "{sample}", "peaks"),
        name="{sample}.macs2",
        peorse=get_peorse,
        genome_size=(
            "hs" if HOST.lower() == "hg38" else
            "mm" if HOST.lower() == "mm39" else
            str(config.get("macs2", {}).get("genome_size", "hs"))
        ),
        qvalue=str(config.get("macs2", {}).get("qvalue", 0.01)),
        shift=str(config.get("macs2", {}).get("shift", -100)),
        extsize=str(config.get("macs2", {}).get("extsize", 200)),
        keep_dup=str(config.get("macs2", {}).get("keep_dup", "all")),
        extra_args=config.get("macs2", {}).get("extra_args", ""),
    threads: 1
    container:
        config["containers"]["macs2"]
    log:
        join(_logdir("macs2_atac_callpeak"), "{sample}.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        mkdir -p {params.outdir}
        if [ "{params.peorse}" == "PE" ]; then
            fmt="BAMPE"
        else
            fmt="BAM"
        fi
        macs2 callpeak \
          -t {input.bam} \
          -f $fmt \
          -g {params.genome_size} \
          -n {params.name} \
          --outdir {params.outdir} \
          --keep-dup {params.keep_dup} \
          --call-summits \
          --nomodel \
          --shift {params.shift} \
          --extsize {params.extsize} \
          -q {params.qvalue} \
          {params.extra_args} \
          > {log} 2>&1
        """


###################################################################################
# Genrich peak calling (ATAC-seq)
###################################################################################

rule genrich_atac_callpeak:
    input:
        bam=join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.final.qname.bam"),
    output:
        narrowpeak=join(RESULTSDIR, "{sample}", "peaks", "{sample}.genrich.narrowPeak"),
    params:
        outdir=join(RESULTSDIR, "{sample}", "peaks"),
        peorse=get_peorse,
        qvalue=str(config.get("genrich", {}).get("qvalue", 0.05)),
        remove_dups=str(config.get("genrich", {}).get("remove_dups", True)),
        junctions=str(config.get("genrich", {}).get("junctions", True)),
        exclude_chr=config.get("genrich", {}).get("exclude_chr", "chrM"),
        blacklist=config.get("genrich", {}).get("blacklist", ""),
        fraglen=str(config.get("genrich", {}).get("fraglen", 500)),
        minlen=str(config.get("genrich", {}).get("minlen", 150)),
        maxlen=str(config.get("genrich", {}).get("maxlen", 1000)),
        extra_args=config.get("genrich", {}).get("extra_args", ""),
    threads: 1
    container:
        config["containers"]["genrich"]
    log:
        join(_logdir("genrich_atac_callpeak"), "{sample}.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        mkdir -p {params.outdir}
        rm_flag=""
        if [ "{params.remove_dups}" == "True" ]; then
            rm_flag="-r"
        fi
        junc_flag=""
        if [ "{params.junctions}" == "True" ]; then
            junc_flag="-j"
        fi
        bl_arg=""
        if [ -n "{params.blacklist}" ]; then
            bl_arg="-E {params.blacklist}"
        fi
        Genrich \
          -t {input.bam} \
          -o {output.narrowpeak} \
          $junc_flag \
          $rm_flag \
          -a {params.fraglen} \
          -l {params.minlen} \
          -g {params.maxlen} \
          -e {params.exclude_chr} \
          $bl_arg \
          -q {params.qvalue} \
          {params.extra_args} \
          > {log} 2>&1
        """
