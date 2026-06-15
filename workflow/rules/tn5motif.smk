###################################################################################
# Tn5 motif extraction from caller-specific BAMs (host + per-virus)
###################################################################################

TN5_WINDOW_SIZE = 2 * int(config.get("tn5_motif", {}).get("flank_size", 10)) + 1
TN5_VIRUS_BIN_SIZE = int(config.get("tn5_motif", {}).get("virus_bin_size", 100))
TN5_HOST_GENE_FLANK_SIZE = int(config.get("tn5_motif", {}).get("host_gene_flank_size", 250))
TN5_TRNA_GENE_FLANK_SIZE = int(config.get("tn5_motif", {}).get("trna_gene_flank_size", 100))
TN5_TRNA_GENE_MAX_SIZE = int(config.get("tn5_motif", {}).get("trna_gene_max_size", 1000))
TN5_MAPQ_MIN = int(config.get("tn5_motif", {}).get("mapq_min", 0))

TN5_EXCLUDE_SECONDARY = bool(config.get("tn5_motif", {}).get("exclude_secondary", False))
TN5_EXCLUDE_SUPPLEMENTARY = bool(
    config.get("tn5_motif", {}).get("exclude_supplementary", False)
)
TN5_DEDUP = bool(config.get("tn5_motif", {}).get("dedup", False))
TN5_LOGO_FORMAT = str(config.get("tn5_motif", {}).get("logo_format", "png"))
TN5_GENERATE_LOGO = bool(config.get("tn5_motif", {}).get("generate_logo", False))
TN5_FRACTIONAL_COUNTING = bool(config.get("tn5_motif", {}).get("fractional_counting", False))
TN5_CALLERS = ("genrich", "macs2")


def _tn5_prefix(sample, group, caller):
    return f"{sample}.{group}.{caller}"


def _tn5_exact_bed(sample, caller, group):
    prefix = _tn5_prefix(sample, group, caller)
    return join(RESULTSDIR, sample, "tn5_motif", caller, group, f"{prefix}.tn5_sites.1bp.bed.gz")


def _tn5_flank_bed(sample, caller, group):
    prefix = _tn5_prefix(sample, group, caller)
    return join(
        RESULTSDIR,
        sample,
        "tn5_motif",
        caller,
        group,
        f"{prefix}.tn5_sites.{TN5_WINDOW_SIZE}bp.bed",
    )


def _tn5_fasta(sample, caller, group):
    prefix = _tn5_prefix(sample, group, caller)
    return join(
        RESULTSDIR,
        sample,
        "tn5_motif",
        caller,
        group,
        f"{prefix}.tn5_sites.{TN5_WINDOW_SIZE}bp.fa",
    )


def _tn5_pfm(sample, caller, group):
    prefix = _tn5_prefix(sample, group, caller)
    return join(
        RESULTSDIR,
        sample,
        "tn5_motif",
        caller,
        group,
        f"{prefix}.tn5_sites.{TN5_WINDOW_SIZE}bp.pfm.tsv.gz",
    )


def _tn5_logo(sample, caller, group):
    prefix = _tn5_prefix(sample, group, caller)
    return join(
        RESULTSDIR,
        sample,
        "tn5_motif",
        caller,
        group,
        f"{prefix}.tn5_sites.{TN5_WINDOW_SIZE}bp.logo.{TN5_LOGO_FORMAT}",
    )


def _caller_bam(sample, caller, group):
    if caller == "genrich":
        if group == "host":
            return join(RESULTSDIR, sample, "align", f"{sample}.aligned.host.qname.bam")
        return join(
            RESULTSDIR,
            sample,
            "align",
            f"{sample}.aligned.virus.{group}.qname.bam",
        )
    if caller == "macs2":
        if group == "host":
            return join(RESULTSDIR, sample, "postprocess", f"{sample}.host.bam")
        return join(
            RESULTSDIR,
            sample,
            "postprocess",
            f"{sample}.virus.{group}.bam",
        )
    raise ValueError(f"Unsupported Tn5 caller: {caller}")


TN5_MOTIF_OUTPUTS = []
for caller in TN5_CALLERS:
    if HOST != "":
        TN5_MOTIF_OUTPUTS.extend(
            expand(_tn5_exact_bed("{sample}", caller, "host"), sample=SAMPLES)
        )
        TN5_MOTIF_OUTPUTS.extend(
            expand(_tn5_pfm("{sample}", caller, "host"), sample=SAMPLES)
        )
        if TN5_GENERATE_LOGO:
            TN5_MOTIF_OUTPUTS.extend(
                expand(_tn5_logo("{sample}", caller, "host"), sample=SAMPLES)
            )
    if len(VIRUS_LIST) > 0:
        TN5_MOTIF_OUTPUTS.extend(
            expand(_tn5_exact_bed("{sample}", caller, "{virus}"), sample=SAMPLES, virus=VIRUS_LIST)
        )
        TN5_MOTIF_OUTPUTS.extend(
            expand(_tn5_pfm("{sample}", caller, "{virus}"), sample=SAMPLES, virus=VIRUS_LIST)
        )
        if TN5_GENERATE_LOGO:
            TN5_MOTIF_OUTPUTS.extend(
                expand(_tn5_logo("{sample}", caller, "{virus}"), sample=SAMPLES, virus=VIRUS_LIST)
            )
