###################################################################################
# bowtie2 alignment
###################################################################################

rule bowtie2_align:
    input:
        index=expand(join(REF_DIR, "ref.fa.{ext}.bt2"), ext=["1", "2", "3", "4", "rev.1", "rev.2"]),
        R1=join(RESULTSDIR, "{sample}", "trim", "{sample}.trimmed_R1.fastq.gz"),
        R2=join(RESULTSDIR, "{sample}", "trim", "{sample}.trimmed_R2.fastq.gz"),
    output:
        bam=temp(join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.sorted.bam")),
        idxstats=join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.idxstats.txt"),
        flagstat=join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.flagstat.txt"),
    params:
        prefix=join(REF_DIR, "ref.fa"),
        outdir=join(RESULTSDIR, "{sample}", "align"),
        tmpdir=join(TEMPDIR, "align", "{sample}"),
        peorse=get_peorse,
        bowtie2_preset=config.get("bowtie2_align", {}).get("preset", "--very-sensitive"),
        bowtie2_k=str(config.get("bowtie2_align", {}).get("k", 10)),
        bowtie2_maxins=str(config.get("bowtie2_align", {}).get("max_insert", 2000)),
        bowtie2_no_mixed="--no-mixed" if config.get("bowtie2_align", {}).get("no_mixed", True) else "",
        bowtie2_no_discordant="--no-discordant" if config.get("bowtie2_align", {}).get("no_discordant", True) else "",
        bowtie2_extra=config.get("bowtie2_align", {}).get("extra_args", ""),
    threads:
        _get_threads("bowtie2_align", profile_config)
    resources:
        runtime=lambda wildcards, attempt: 960 * attempt,
    container:
        config["containers"]["bowtie2"]
    log:
        join(_logdir("bowtie2_align"), "{sample}.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        echo "[bowtie2_align] sample={wildcards.sample} peorse={params.peorse}"
        echo "[bowtie2_align] R1={input.R1}"
        echo "[bowtie2_align] R2={input.R2}"
        echo "[bowtie2_align] index_prefix={params.prefix} threads={threads}"
        echo "[bowtie2_align] outdir={params.outdir} tmpdir={params.tmpdir}"
        mkdir -p {params.outdir}
        mkdir -p {params.tmpdir}
        unsorted_bam={params.tmpdir}/{wildcards.sample}.aligned.bam
        if [ "{params.peorse}" == "PE" ]; then
            echo "[bowtie2_align] running bowtie2 (PE)"
            bowtie2 -x {params.prefix} \
              -1 {input.R1} -2 {input.R2} \
              -p {threads} \
              {params.bowtie2_preset} \
              -k {params.bowtie2_k} \
              -X {params.bowtie2_maxins} \
              {params.bowtie2_no_mixed} {params.bowtie2_no_discordant} \
              {params.bowtie2_extra} \
              2> {log} \
            | samtools view -@ {threads} -bS - \
            > $unsorted_bam
        else
            echo "[bowtie2_align] running bowtie2 (SE)"
            bowtie2 -x {params.prefix} \
              -U {input.R1} \
              -p {threads} \
              {params.bowtie2_preset} \
              -k {params.bowtie2_k} \
              -X {params.bowtie2_maxins} \
              {params.bowtie2_no_mixed} {params.bowtie2_no_discordant} \
              {params.bowtie2_extra} \
              2> {log} \
            | samtools view -@ {threads} -bS - \
            > $unsorted_bam
        fi
        echo "[bowtie2_align] sorting bam"
        samtools sort -@ {threads} -T {params.tmpdir}/sort -o {output.bam} $unsorted_bam
        echo "[bowtie2_align] indexing bam"
        samtools index {output.bam}
        echo "[bowtie2_align] idxstats"
        samtools idxstats {output.bam} > {output.idxstats}
        echo "[bowtie2_align] flagstat"
        samtools flagstat {output.bam} > {output.flagstat}
        echo "[bowtie2_align] done"
        """


###################################################################################
# bowtie2 filtering
###################################################################################

rule bowtie2_filter:
    input:
        bam=join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.sorted.bam"),
    output:
        clean_bam=temp(join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.clean.sorted.bam")),
        clean_idxstats=join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.clean.idxstats.txt"),
        clean_flagstat=join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.clean.flagstat.txt"),
        fixmate_bam=temp(join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.fixmate.sorted.bam")),
        fixmate_idxstats=join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.fixmate.idxstats.txt"),
        fixmate_flagstat=join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.fixmate.flagstat.txt"),
        dedup_bam=temp(join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.dedup.bam")),
        dedup_idxstats=join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.dedup.idxstats.txt"),
        dedup_flagstat=join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.dedup.flagstat.txt"),
        final_bam=join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.final.bam"),
        final_idxstats=join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.final.idxstats.txt"),
        final_flagstat=join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.final.flagstat.txt"),
        final_bai=join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.final.bam.bai"),
    params:
        outdir=join(RESULTSDIR, "{sample}", "align"),
        tmpdir=join(TEMPDIR, "align", "{sample}", "filter"),
        peorse=get_peorse,
        flag_include=str(config.get("bowtie2_align", {}).get("filter", {}).get("include_flag", 2)),
        flag_exclude=str(config.get("bowtie2_align", {}).get("filter", {}).get("exclude_flag", "0x904")),
        mapq_arg=(
            "" if config.get("bowtie2_align", {}).get("filter", {}).get("mapq", 30) in ["", None]
            else f"-q {config.get('bowtie2_align', {}).get('filter', {}).get('mapq', 30)}"
        ),
        exclude_rname_expr=(
            " && ".join(
                [f'rname!=\"{rname}\"' for rname in config.get("bowtie2_align", {}).get("filter", {}).get("exclude_rnames", ["chrM", "MT"])]
            )
            if config.get("bowtie2_align", {}).get("filter", {}).get("exclude_rnames", ["chrM", "MT"])
            else "1"
        ),
        markdup_enable=str(config.get("bowtie2_align", {}).get("filter", {}).get("markdup_enable", True)),
        markdup_remove_arg="-r" if config.get("bowtie2_align", {}).get("filter", {}).get("markdup_remove", True) else "",
    threads:
        _get_threads("bowtie2_filter", profile_config)
    container:
        config["containers"]["bowtie2"]
    log:
        join(_logdir("bowtie2_filter"), "{sample}.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        echo "[bowtie2_filter] sample={wildcards.sample} peorse={params.peorse}"
        echo "[bowtie2_filter] flags: include={params.flag_include} exclude={params.flag_exclude} mapq={params.mapq_arg}"
        echo "[bowtie2_filter] exclude_rnames_expr={params.exclude_rname_expr}"
        mkdir -p {params.outdir}
        mkdir -p {params.tmpdir}

        if [ "{params.peorse}" == "PE" ]; then
            echo "[bowtie2_filter] filtering proper pairs"
            rm -f {params.tmpdir}/clean.*.bam
            samtools view -@ {threads} -b -f {params.flag_include} -F {params.flag_exclude} {input.bam} \
              | samtools sort -@ {threads} -T {params.tmpdir}/clean -o {output.clean_bam} -
        else
            echo "[bowtie2_filter] filtering (SE, no include flag)"
            rm -f {params.tmpdir}/clean.*.bam
            samtools view -@ {threads} -b -F {params.flag_exclude} {input.bam} \
              | samtools sort -@ {threads} -T {params.tmpdir}/clean -o {output.clean_bam} -
        fi
        if [ ! -s {output.clean_bam}.bai ]; then
            samtools index {output.clean_bam}
        fi
        samtools idxstats {output.clean_bam} > {output.clean_idxstats}
        samtools flagstat {output.clean_bam} > {output.clean_flagstat}

        if [ "{params.peorse}" == "PE" ]; then
            echo "[bowtie2_filter] fixmate"
            rm -f {params.tmpdir}/clean_qname.*.bam
            samtools sort -@ {threads} -n -T {params.tmpdir}/clean_qname -o {params.tmpdir}/clean.qname.bam {output.clean_bam}
            samtools fixmate -@ {threads} -m {params.tmpdir}/clean.qname.bam {params.tmpdir}/fixmate.bam
            rm -f {params.tmpdir}/fixmate_sort.*.bam
            samtools sort -@ {threads} -T {params.tmpdir}/fixmate_sort -o {output.fixmate_bam} {params.tmpdir}/fixmate.bam
        else
            cp -f {output.clean_bam} {output.fixmate_bam}
        fi
        if [ ! -s {output.fixmate_bam}.bai ]; then
            samtools index {output.fixmate_bam}
        fi
        samtools idxstats {output.fixmate_bam} > {output.fixmate_idxstats}
        samtools flagstat {output.fixmate_bam} > {output.fixmate_flagstat}

        if [ "{params.markdup_enable}" == "True" ]; then
            echo "[bowtie2_filter] markdup"
            samtools markdup -@ {threads} {params.markdup_remove_arg} {output.fixmate_bam} {output.dedup_bam}
        else
            cp -f {output.fixmate_bam} {output.dedup_bam}
        fi
        if [ ! -s {output.dedup_bam}.bai ]; then
            samtools index {output.dedup_bam}
        fi
        samtools idxstats {output.dedup_bam} > {output.dedup_idxstats}
        samtools flagstat {output.dedup_bam} > {output.dedup_flagstat}

        echo "[bowtie2_filter] mapq + rname filtering"
        samtools view -@ {threads} -b {params.mapq_arg} {output.dedup_bam} \
          | samtools view -@ {threads} -b -e '{params.exclude_rname_expr}' \
          > {output.final_bam}
        if [ ! -s {output.final_bam}.bai ]; then
            samtools index {output.final_bam}
        fi
        samtools idxstats {output.final_bam} > {output.final_idxstats}
        samtools flagstat {output.final_bam} > {output.final_flagstat}
        echo "[bowtie2_filter] done"
        """


###################################################################################
# split aligned.sorted.bam into host/virus and qname-sort for Genrich
###################################################################################

rule split_host_qname_bam:
    input:
        bam=join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.sorted.bam"),
        regions=REF_REGIONS_HOST,
    output:
        qname_bam=join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.host.qname.bam"),
    params:
        tmpdir=join(TEMPDIR, "align", "{sample}", "split_host_qname"),
    threads:
        _get_threads("split_host_qname_bam", profile_config)
    container:
        config["containers"]["bowtie2"]
    log:
        join(_logdir("split_host_qname_bam"), "{sample}.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        mkdir -p {params.tmpdir}
        host_regions=$(awk -F '\t' '{{n=split($2,a," "); for(i=1;i<=n;i++) if (a[i]!="") print a[i]}}' {input.regions} | sort -u | tr '\n' ' ')
        samtools view -@ {threads} -b {input.bam} $host_regions \
          | samtools sort -@ {threads} -n -T {params.tmpdir}/qname -o {output.qname_bam} -
        """


rule split_virus_qname_bam:
    input:
        bam=join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.sorted.bam"),
        regions=join(FASTAS_GTFS_DIR, "{virus}.fa.regions"),
    output:
        qname_bam=join(RESULTSDIR, "{sample}", "align", "{sample}.aligned.virus.{virus}.qname.bam"),
    params:
        tmpdir=join(TEMPDIR, "align", "{sample}", "split_virus_qname", "{virus}"),
    threads:
        _get_threads("split_virus_qname_bam", profile_config)
    container:
        config["containers"]["bowtie2"]
    log:
        join(_logdir("split_virus_qname_bam"), "{sample}.{virus}.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        mkdir -p {params.tmpdir}
        virus_regions=$(awk -F '\t' '{{n=split($2,a," "); for(i=1;i<=n;i++) if (a[i]!="") print a[i]}}' {input.regions} | sort -u | tr '\n' ' ')
        samtools view -@ {threads} -b {input.bam} $virus_regions \
          | samtools sort -@ {threads} -n -T {params.tmpdir}/qname -o {output.qname_bam} -
        """
