###################################################################################
# Convert browser-facing BED-like outputs to BigBed for faster UCSC loading
###################################################################################

BIG_NARROWPEAK_AS = join(RESOURCES_DIR, "bigNarrowPeak.as")

BIGBED_OUTPUTS = []
if HOST != "":
    BIGBED_OUTPUTS.append(join(REF_DIR, "ref.tss.host.bb"))
    BIGBED_OUTPUTS.extend(
        expand(join(RESULTSDIR, "{sample}", "peaks", "{sample}.host.macs2_summits.bb"), sample=CASE_SAMPLES)
    )
    BIGBED_OUTPUTS.extend(
        expand(join(RESULTSDIR, "{sample}", "peaks", "{sample}.host.macs2_peaks.bb"), sample=CASE_SAMPLES)
    )
    BIGBED_OUTPUTS.extend(
        expand(join(RESULTSDIR, "{sample}", "peaks", "{sample}.host.genrich.bb"), sample=CASE_SAMPLES)
    )

if len(VIRUS_LIST) > 0:
    BIGBED_OUTPUTS.extend(expand(join(REF_DIR, "ref.tss.{virus}.bb"), virus=VIRUS_LIST))
    BIGBED_OUTPUTS.extend(
        expand(
            join(RESULTSDIR, "{sample}", "peaks", "{sample}.virus.{virus}.macs2_summits.bb"),
            sample=CASE_SAMPLES,
            virus=VIRUS_LIST,
        )
    )
    BIGBED_OUTPUTS.extend(
        expand(
            join(RESULTSDIR, "{sample}", "peaks", "{sample}.virus.{virus}.macs2_peaks.bb"),
            sample=CASE_SAMPLES,
            virus=VIRUS_LIST,
        )
    )
    BIGBED_OUTPUTS.extend(
        expand(
            join(RESULTSDIR, "{sample}", "peaks", "{sample}.virus.{virus}.genrich.bb"),
            sample=CASE_SAMPLES,
            virus=VIRUS_LIST,
        )
    )


rule tss_host_to_bigbed:
    input:
        bed=join(REF_DIR, "ref.tss.host.bed"),
        chromsizes=join(REF_DIR, "ref.chrom.sizes.host.txt"),
    output:
        bb=join(REF_DIR, "ref.tss.host.bb"),
    params:
        tmpdir=join(TEMPDIR, "bigbed", "ref", "host"),
    threads:
        _get_threads("tss_host_to_bigbed", profile_config)
    container:
        config["containers"]["bedToBigBed"]
    log:
        join(_logdir("tss_host_to_bigbed"), "host.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        bash {SCRIPTS_DIR}/bed_to_bigbed.sh \
          --input {input.bed} --output {output.bb} --chromsizes {input.chromsizes} \
          --tmpdir {params.tmpdir} --bed-type bed3
        """


rule tss_virus_to_bigbed:
    input:
        bed=join(REF_DIR, "ref.tss.{virus}.bed"),
        chromsizes=join(REF_DIR, "ref.chrom.sizes.{virus}.txt"),
    output:
        bb=join(REF_DIR, "ref.tss.{virus}.bb"),
    params:
        tmpdir=join(TEMPDIR, "bigbed", "ref", "{virus}"),
    threads:
        _get_threads("tss_virus_to_bigbed", profile_config)
    container:
        config["containers"]["bedToBigBed"]
    log:
        join(_logdir("tss_virus_to_bigbed"), "{virus}.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        bash {SCRIPTS_DIR}/bed_to_bigbed.sh \
          --input {input.bed} --output {output.bb} --chromsizes {input.chromsizes} \
          --tmpdir {params.tmpdir} --bed-type bed3
        """


rule macs2_summits_host_to_bigbed:
    input:
        bed=join(RESULTSDIR, "{sample}", "peaks", "{sample}.host.macs2_summits.bed.gz"),
        chromsizes=join(REF_DIR, "ref.chrom.sizes.host.txt"),
    output:
        bb=join(RESULTSDIR, "{sample}", "peaks", "{sample}.host.macs2_summits.bb"),
    params:
        tmpdir=join(TEMPDIR, "bigbed", "{sample}", "host", "macs2_summits"),
    threads:
        _get_threads("macs2_summits_host_to_bigbed", profile_config)
    container:
        config["containers"]["bedToBigBed"]
    log:
        join(_logdir("macs2_summits_host_to_bigbed"), "{sample}.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        bash {SCRIPTS_DIR}/bed_to_bigbed.sh \
          --input {input.bed} --output {output.bb} --chromsizes {input.chromsizes} \
          --tmpdir {params.tmpdir} --bed-type bed5 --gzipped --clamp-score --clamp-mode conditional
        """


rule macs2_summits_virus_to_bigbed:
    input:
        bed=join(RESULTSDIR, "{sample}", "peaks", "{sample}.virus.{virus}.macs2_summits.bed.gz"),
        chromsizes=join(REF_DIR, "ref.chrom.sizes.{virus}.txt"),
    output:
        bb=join(RESULTSDIR, "{sample}", "peaks", "{sample}.virus.{virus}.macs2_summits.bb"),
    params:
        tmpdir=join(TEMPDIR, "bigbed", "{sample}", "{virus}", "macs2_summits"),
    threads:
        _get_threads("macs2_summits_virus_to_bigbed", profile_config)
    container:
        config["containers"]["bedToBigBed"]
    log:
        join(_logdir("macs2_summits_virus_to_bigbed"), "{sample}.{virus}.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        bash {SCRIPTS_DIR}/bed_to_bigbed.sh \
          --input {input.bed} --output {output.bb} --chromsizes {input.chromsizes} \
          --tmpdir {params.tmpdir} --bed-type bed5 --gzipped --clamp-score --clamp-mode conditional
        """


rule narrowpeak_host_to_bigbed:
    input:
        bed=join(RESULTSDIR, "{sample}", "peaks", "{sample}.host.{caller}.narrowPeak.gz"),
        chromsizes=join(REF_DIR, "ref.chrom.sizes.host.txt"),
        autosql=BIG_NARROWPEAK_AS,
    output:
        bb=join(RESULTSDIR, "{sample}", "peaks", "{sample}.host.{caller}.bb"),
    params:
        tmpdir=join(TEMPDIR, "bigbed", "{sample}", "host", "{caller}"),
    threads:
        _get_threads("narrowpeak_host_to_bigbed", profile_config)
    wildcard_constraints:
        caller="macs2_peaks|genrich"
    container:
        config["containers"]["bedToBigBed"]
    log:
        join(_logdir("narrowpeak_host_to_bigbed"), "{sample}.{caller}.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        bash {SCRIPTS_DIR}/bed_to_bigbed.sh \
          --input {input.bed} --output {output.bb} --chromsizes {input.chromsizes} \
          --tmpdir {params.tmpdir} --bed-type bed6+4 --gzipped --clamp-score --clamp-mode unconditional \
          --tab --autosql {input.autosql}
        """


rule narrowpeak_virus_to_bigbed:
    input:
        bed=join(RESULTSDIR, "{sample}", "peaks", "{sample}.virus.{virus}.{caller}.narrowPeak.gz"),
        chromsizes=join(REF_DIR, "ref.chrom.sizes.{virus}.txt"),
        autosql=BIG_NARROWPEAK_AS,
    output:
        bb=join(RESULTSDIR, "{sample}", "peaks", "{sample}.virus.{virus}.{caller}.bb"),
    params:
        tmpdir=join(TEMPDIR, "bigbed", "{sample}", "{virus}", "{caller}"),
    threads:
        _get_threads("narrowpeak_virus_to_bigbed", profile_config)
    wildcard_constraints:
        caller="macs2_peaks|genrich"
    container:
        config["containers"]["bedToBigBed"]
    log:
        join(_logdir("narrowpeak_virus_to_bigbed"), "{sample}.{virus}.{caller}.log")
    shell:
        r"""
        set -exo pipefail
        mkdir -p $(dirname {log})
        exec > >(tee -a {log}) 2>&1
        bash {SCRIPTS_DIR}/bed_to_bigbed.sh \
          --input {input.bed} --output {output.bb} --chromsizes {input.chromsizes} \
          --tmpdir {params.tmpdir} --bed-type bed6+4 --gzipped --clamp-score --clamp-mode unconditional \
          --tab --autosql {input.autosql}
        """
