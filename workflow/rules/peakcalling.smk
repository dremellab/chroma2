###################################################################################
# MACS2 peak calling (ATAC-seq) - host + per-virus
###################################################################################

rule macs2_atac_callpeak_host:
    input:
        bam=join(RESULTSDIR, "{sample}", "postprocess", "{sample}.host.bam"),
        bai=join(RESULTSDIR, "{sample}", "postprocess", "{sample}.host.bam.bai"),
        control=lambda wc: _opt_input(get_host_control_filtered_bam(wc)),
    output:
        narrowpeak=join(RESULTSDIR, "{sample}", "peaks", "{sample}.host.macs2_peaks.narrowPeak.gz"),
        summits=join(RESULTSDIR, "{sample}", "peaks", "{sample}.host.macs2_summits.bed.gz"),
        xls=temp(join(RESULTSDIR, "{sample}", "peaks", "{sample}.host.macs2_peaks.xls")),
    params:
        outdir=join(RESULTSDIR, "{sample}", "peaks"),
        name="{sample}.host.macs2",
        peorse=get_peorse,
        control_arg=lambda wc, input: (
            f"-c {input.control}"
            if config.get("peakcalling", {}).get("use_host_input", False) and input.control
            else ""
        ),
        genome_size=_macs2_genome_size(
            HOST, str(config.get("macs2", {}).get("genome_size", "hs"))
        ),
        qvalue=str(config.get("macs2", {}).get("qvalue", 0.01)),
        shift=str(config.get("macs2", {}).get("shift", -100)),
        extsize=str(config.get("macs2", {}).get("extsize", 200)),
        keep_dup=str(config.get("macs2", {}).get("keep_dup", "all")),
        extra_args=config.get("macs2", {}).get("extra_args", ""),
    threads:
        _get_threads("macs2_atac_callpeak_host", profile_config)
    container:
        config["containers"]["macs2"]
    log:
        join(_logdir("macs2_atac_callpeak_host"), "{sample}.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        mkdir -p {params.outdir}
        if [ "{params.peorse}" == "PE" ]; then
            fmt="BAMPE"
        else
            fmt="BAM"
        fi
        macs2 callpeak \
          -t {input.bam} {params.control_arg} \
          -f $fmt \
          -g {params.genome_size} \
          -n {params.name} \
          --outdir {params.outdir} \
          --keep-dup {params.keep_dup} \
          --call-summits \
          --nomodel \
          --shift {params.shift} \
          --extsize {params.extsize} \
          -q {params.qvalue} {params.extra_args}
        gzip -f {params.outdir}/{params.name}_peaks.narrowPeak
        gzip -f {params.outdir}/{params.name}_summits.bed
        """


rule macs2_atac_callpeak_virus:
    input:
        bam=join(RESULTSDIR, "{sample}", "postprocess", "{sample}.virus.{virus}.bam"),
        bai=join(RESULTSDIR, "{sample}", "postprocess", "{sample}.virus.{virus}.bam.bai"),
        fasta=join(FASTAS_GTFS_DIR, "{virus}.fa"),
        control=lambda wc: _opt_input(get_virus_control_filtered_bam(wc)),
    output:
        narrowpeak=join(RESULTSDIR, "{sample}", "peaks", "{sample}.virus.{virus}.macs2_peaks.narrowPeak.gz"),
        summits=join(RESULTSDIR, "{sample}", "peaks", "{sample}.virus.{virus}.macs2_summits.bed.gz"),
        xls=temp(join(RESULTSDIR, "{sample}", "peaks", "{sample}.virus.{virus}.macs2_peaks.xls")),
    params:
        outdir=join(RESULTSDIR, "{sample}", "peaks"),
        name="{sample}.virus.{virus}.macs2",
        peorse=get_peorse,
        control_arg=lambda wc, input: (
            f"-c {input.control}"
            if config.get("peakcalling", {}).get("use_virus_input", True) and input.control
            else ""
        ),
        qvalue=str(config.get("macs2", {}).get("qvalue", 0.01)),
        shift=str(config.get("macs2", {}).get("shift", -100)),
        extsize=str(config.get("macs2", {}).get("extsize", 200)),
        keep_dup=str(config.get("macs2", {}).get("keep_dup", "all")),
        extra_args=config.get("macs2", {}).get("extra_args", ""),
    threads:
        _get_threads("macs2_atac_callpeak_virus", profile_config)
    container:
        config["containers"]["macs2"]
    log:
        join(_logdir("macs2_atac_callpeak_virus"), "{sample}.{virus}.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        mkdir -p {params.outdir}
        if [ "{params.peorse}" == "PE" ]; then
            fmt="BAMPE"
        else
            fmt="BAM"
        fi
        virus_genome_size=$(awk 'BEGIN{{n=0}} /^>/{{next}} {{n+=length($0)}} END{{print n}}' {input.fasta})
        macs2 callpeak \
          -t {input.bam} {params.control_arg} \
          -f $fmt \
          -g $virus_genome_size \
          -n {params.name} \
          --outdir {params.outdir} \
          --keep-dup {params.keep_dup} \
          --call-summits \
          --nomodel \
          --shift {params.shift} \
          --extsize {params.extsize} \
          -q {params.qvalue} {params.extra_args}
        gzip -f {params.outdir}/{params.name}_peaks.narrowPeak
        gzip -f {params.outdir}/{params.name}_summits.bed
        """


###################################################################################
# Genrich peak calling (ATAC-seq) - host + per-virus
###################################################################################

rule genrich_atac_callpeak_host:
    input:
        bam=join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.host.qname.bam"),
        control=lambda wc: _opt_input(get_host_control_qname_bam(wc)),
    output:
        narrowpeak=join(RESULTSDIR, "{sample}", "peaks", "{sample}.host.genrich.narrowPeak.gz"),
    params:
        outdir=join(RESULTSDIR, "{sample}", "peaks"),
        peorse=get_peorse,
        narrowpeak_uncompressed=lambda wc, output: output.narrowpeak[:-3],
        control_arg=lambda wc, input: (
            f"-c {input.control}"
            if config.get("peakcalling", {}).get("use_host_input", False) and input.control
            else ""
        ),
        qvalue=str(config.get("genrich", {}).get("qvalue", 0.05)),
        remove_dups=str(config.get("genrich", {}).get("remove_dups", True)),
        junctions=str(config.get("genrich", {}).get("junctions", True)),
        exclude_chr=config.get("genrich", {}).get("exclude_chr", "chrM"),
        exclude_arg=(
            f"-e {config.get('genrich', {}).get('exclude_chr', 'chrM')}"
            if config.get("genrich", {}).get("exclude_chr", "chrM")
            else ""
        ),
        blacklist=config.get("genrich", {}).get("blacklist", ""),
        fraglen=str(
            config.get("genrich", {}).get(
                "host_fraglen", config.get("genrich", {}).get("fraglen", 500)
            )
        ),
        mval=str(config.get("genrich", {}).get("host_mval", 30)),
        minlen=str(config.get("genrich", {}).get("minlen", 150)),
        maxlen=str(config.get("genrich", {}).get("maxlen", 1000)),
        extra_args=config.get("genrich", {}).get("extra_args", ""),
    threads:
        _get_threads("genrich_atac_callpeak_host", profile_config)
    container:
        config["containers"]["genrich"]
    log:
        join(_logdir("genrich_atac_callpeak_host"), "{sample}.log")
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
          -t {input.bam} {params.control_arg} \
          -o {params.narrowpeak_uncompressed} \
          $junc_flag \
          $rm_flag \
          -m {params.mval} \
          -a {params.fraglen} \
          -l {params.minlen} \
          -g {params.maxlen} \
          {params.exclude_arg} \
          $bl_arg \
          -q {params.qvalue} {params.extra_args} \
          2>&1 | tee -a {log}
        gzip -f {params.narrowpeak_uncompressed}
        """


rule genrich_atac_callpeak_virus:
    input:
        bam=join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.virus.{virus}.qname.bam"),
        control=lambda wc: _opt_input(get_virus_control_qname_bam(wc)),
    output:
        narrowpeak=join(RESULTSDIR, "{sample}", "peaks", "{sample}.virus.{virus}.genrich.narrowPeak.gz"),
    params:
        outdir=join(RESULTSDIR, "{sample}", "peaks"),
        peorse=get_peorse,
        narrowpeak_uncompressed=lambda wc, output: output.narrowpeak[:-3],
        control_arg=lambda wc, input: (
            f"-c {input.control}"
            if config.get("peakcalling", {}).get("use_virus_input", True) and input.control
            else ""
        ),
        qvalue=str(config.get("genrich", {}).get("qvalue", 0.05)),
        remove_dups=str(config.get("genrich", {}).get("remove_dups", True)),
        junctions=str(config.get("genrich", {}).get("junctions", True)),
        exclude_chr=config.get("genrich", {}).get("virus_exclude_chr", ""),
        exclude_arg=(
            f"-e {config.get('genrich', {}).get('virus_exclude_chr', '')}"
            if config.get("genrich", {}).get("virus_exclude_chr", "")
            else ""
        ),
        blacklist=config.get("genrich", {}).get("virus_blacklist", ""),
        fraglen=str(
            config.get("genrich", {}).get(
                "virus_fraglen", config.get("genrich", {}).get("fraglen", 500)
            )
        ),
        mval=str(config.get("genrich", {}).get("virus_mval", 5)),
        minlen=str(config.get("genrich", {}).get("minlen", 150)),
        maxlen=str(config.get("genrich", {}).get("maxlen", 1000)),
        extra_args=config.get("genrich", {}).get("extra_args", ""),
    threads:
        _get_threads("genrich_atac_callpeak_virus", profile_config)
    container:
        config["containers"]["genrich"]
    log:
        join(_logdir("genrich_atac_callpeak_virus"), "{sample}.{virus}.log")
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
          -t {input.bam} {params.control_arg} \
          -o {params.narrowpeak_uncompressed} \
          $junc_flag \
          $rm_flag \
          -m {params.mval} \
          -a {params.fraglen} \
          -l {params.minlen} \
          -g {params.maxlen} \
          {params.exclude_arg} \
          $bl_arg \
          -q {params.qvalue} {params.extra_args} \
          2>&1 | tee -a {log}
        gzip -f {params.narrowpeak_uncompressed}
        """
