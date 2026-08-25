# 📥 CHROMA2 Input Requirements

## 1. Sample Manifest (Required for Initialization) {#sample-manifest}

`samples.tsv` is a tab-delimited file, copied into the working directory at `init` time via `-s/--manifest` (or taken from the pipeline's checked-in `config/samples.tsv` template if omitted).

| Column             | Required                                     | Description                                                                                                                         |
| ------------------ | -------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| `sampleName`       | **Required**, unique                         | Sample identifier, used throughout output file/directory naming                                                                     |
| `groupName`        | **Required**, non-empty                      | Biological group/condition — used for DESeq2 contrasts (`contrasts.tsv`) and count-matrix grouping                                  |
| `batch`            | Optional                                     | Batch label; blank values are normalized to `batch1`                                                                                |
| `path_to_R1_fastq` | **Required**, must exist/be readable         | Absolute path to R1 FASTQ (`.fastq.gz`)                                                                                             |
| `path_to_R2_fastq` | Optional                                     | Absolute path to R2 FASTQ; leave blank for single-end samples — presence/absence determines PE vs. SE per sample                    |
| `role`             | **Required**, one of `case`, `control`       | Whether this sample is an experimental case or an input/IgG-style control used as a peak-calling background                         |
| `target`           | **Required**, one of `host`, `virus`, `both` | For control samples, which alignment space (host and/or virus) this control applies to                                              |
| `host_input_pool`  | Optional                                     | Pool name grouping this control sample's **host** reads with other controls into a shared input BAM, used as `-c` for MACS2/Genrich |
| `virus_input_pool` | Optional                                     | Same, for **virus**-space pooling                                                                                                   |

### Example Sample Manifest

```text
sampleName	groupName	batch	path_to_R1_fastq	path_to_R2_fastq	role	target	host_input_pool	virus_input_pool
KOS_1h_rep1	KOS_1h		/data/KOS_1h_rep1_R1.fastq.gz	/data/KOS_1h_rep1_R2.fastq.gz	case	both
KOS_1h_rep2	KOS_1h		/data/KOS_1h_rep2_R1.fastq.gz	/data/KOS_1h_rep2_R2.fastq.gz	case	both
KOS_4h_rep1	KOS_4h		/data/KOS_4h_rep1_R1.fastq.gz	/data/KOS_4h_rep1_R2.fastq.gz	case	both
Input_ctrl_rep1	Input	batch1	/data/Input_ctrl_rep1_R1.fastq.gz	/data/Input_ctrl_rep1_R2.fastq.gz	control	virus		VPOOL_ALL
Input_ctrl_rep2	Input	batch1	/data/Input_ctrl_rep2_R1.fastq.gz	/data/Input_ctrl_rep2_R2.fastq.gz	control	virus		VPOOL_ALL
```

### Supported Library Types

Both paired-end and single-end libraries are supported, determined per sample by whether `path_to_R2_fastq` is populated — there is no separate config flag.

### Validation Rules

Applied when Snakemake first parses the manifest (during `dryrun`/`run`/`runlocal`, not at `init`):

- `sampleName` values must be unique across the manifest.
- `groupName` must be non-empty for every sample.
- `path_to_R1_fastq` must point to an existing, readable file for every sample.
- `role` must be exactly `case` or `control`.
- `target` must be exactly `host`, `virus`, or `both`.
- Any missing required column raises an error before any rules run.

### Input-Pooling Mechanism

Control samples (`role: control`) are not aligned/filtered any differently from case samples — the pooling happens downstream, when their host and/or virus alignments are merged by pool name:

1. All `control` samples sharing the same `host_input_pool` value have their host-space BAMs (`aligned.host.qname.bam` for Genrich, `postprocess/*.host.bam` for MACS2) merged into one pooled control BAM per pool.
2. The same happens per virus for `virus_input_pool`.
3. A `case` sample's peak-calling rule looks up the pool assigned to whichever control samples are configured, and passes that pooled BAM as `-c` — but **only if** `peakcalling.use_host_input`/`use_virus_input` is `true` in `config.yaml` (see [§3](#reference-configuration) below; defaults are host input **off**, virus input **on**).
4. If no pool resolves for a given case sample/organism, peak calling proceeds without a control (`-c` omitted).

## 2. Working Directory (`--workdir`) {#working-directory}

Set via `-w/--workdir` on every `chroma2` invocation. Must not already exist at `init` time; all other runmodes require it to already exist with a valid `config.yaml` inside.

## 3. Reference Configuration (Host, Additives, and Viruses) {#reference-configuration}

### Parameters

- `--host|-g` — a single host genome name. Must have an entry in `reference_gtf` (currently `mm39`, `hg38`, `hg38_basic`, `hs1`). Default: `hg38_basic`.
- `--additives|-a` — comma-separated list of additive reference contigs (e.g. `ERCC`, `BAC16Insert`). Default: none.
- `--viruses|-v` — comma-separated list of viral genome accessions to co-align against. Default: `KT899744.1`.

At least one of host or viruses must be non-empty, or the pipeline exits with an error.

### Supported Host Genomes

The four `--host` values are distinct genome builds, not just naming variants — different FASTA, different annotation source, and (except `hg38_basic`) a custom `chrR` contig added specifically so reads from ribosomal DNA repeats, which collapse/misassemble in a standard reference, map correctly:

| Host         | Assembly                       | Chromosomes                    | Annotation source                                                     | Genes  | `chrR`? |
| ------------ | ------------------------------ | ------------------------------ | --------------------------------------------------------------------- | ------ | ------- |
| `hg38_basic` | GRCh38 primary assembly (NCBI) | chr1–22, X, Y, chrM (25)       | GENCODE v38 (comprehensive: gene/transcript/exon/CDS/UTR/codons)      | 60,649 | No      |
| `hg38`       | GRCh38, rDNA-enriched          | chr1–22, X, Y, chrM, chrR (26) | NCBI RefSeq, rDNA-mapping genome (minimal: gene/transcript/exon only) | 28,091 | Yes     |
| `mm39`       | GRCm39, rDNA-enriched          | chr1–19, X, Y, chrM, chrR (23) | rDNA-mapping genome (mouse)                                           | 78,277 | Yes     |
| `hs1`        | T2T-CHM13v2.0, rDNA-enriched   | chr1–22, X, Y, chrM, chrR (26) | rDNA-mapping genome (T2T)                                             | 28,363 | Yes     |

**Choosing between `hg38` and `hg38_basic`** (the only build offered in both forms): `hg38_basic` is the general-purpose choice — comprehensive GENCODE gene models, full CDS/UTR detail, best for standard RNA-seq/splice-junction/protein-coding analysis (this is why it's the default). `hg38` trades that annotation depth for the `chrR` contig, which matters when rRNA quantification is a goal — without a dedicated rDNA-mapping chromosome, reads from the highly repetitive rDNA locus map ambiguously or get lost against the standard assembly. `mm39` and `hs1` only ship as their rDNA-enriched build, with `chrR` included by default.

### Reference Data Paths

| Type                             | Default Location                                                                                      |
| -------------------------------- | ----------------------------------------------------------------------------------------------------- |
| Host/additive/virus FASTA + GTFs | `fastas_gtfs_dir` in `config.yaml`, default `/project/dremel_lab/workflows/reference_data/fasta_gtf/` |
| Kraken2 database                 | `kraken2_db` in `config.yaml`                                                                         |

### Per-Host Annotation GTFs

CHROMA2 resolves several optional/required annotation GTFs per selected host genome from `config.yaml`:

```yaml
reference_gtf:
  mm39: "mm39.no_tRNA.gtf"
  hg38: "hg38-rDNA_v1.0_enriched.canonical_chromosomes_plus_chrR.gtf"
  hg38_basic: "hg38_basic.canonical_chromosomes.gtf"
  hs1: "hs1-rDNA_v1.0_enriched.canonical_plus_chrR.gtf"
trnas_gtf: { mm39: "...", hg38: "...", hg38_basic: "...", hs1: "..." }
chrr_gtf: { mm39: "...", hg38: "...", hg38_basic: "...", hs1: "..." }
pol3_gtf: { Pol3_T1: { ... }, Pol3_T2: { ... }, Pol3_T3: { ... } }
repeat_elements_gtf:
  {
    SINE_Alu: { ... },
    SINE_MIR: { ... },
    LINE_L1: { ... },
    LINE_L2: { ... },
    LTR: { ... },
    other_repeat_elements: { ... },
  }
```

- **`reference_gtf`** — main annotation, required for every selected host; used for the composite `ref.gtf`, TSS beds, and the gene count matrix.
- **`trnas_gtf`** — tRNA gene annotations; used for the tRNA count matrix. If unresolved for the selected host, the tRNA matrix is silently unavailable.
- **`chrr_gtf`** — rRNA ("chrR") gene annotations; used for the rRNA count matrix.
- **`pol3_gtf`** — three separate Pol III gene-class GTFs (T1/T2/T3); each unresolved class is dropped from the available Pol III matrices.
- **`repeat_elements_gtf`** — six separate repeat-class GTFs (RepeatMasker-derived); each unresolved class is dropped from the available repeat-element matrices.

### Validation Rules

Before any pipeline rules run, CHROMA2 prints a validation report confirming every GTF required for the selected host/viruses exists and is readable, and exits non-zero (with the specific missing paths listed) if any is not.

## 4. Contrasts Manifest (Optional, for DESeq2) {#contrasts-manifest}

`contrasts.tsv` is a 2-column tab-delimited file, auto-copied into the workdir at `init` from the pipeline's `config/contrasts.tsv`. Only consumed if `deseq2.enabled: true`.

```text
group1	group2
KOS_1h	KOS_4h
KOS_1h	KOS_8h
KOS_4h	KOS_8h
```

- Header must be exactly `group1` / `group2`.
- Each row references two `groupName` values from `samples.tsv`; both must exist and must differ from each other.
- Each group referenced in a contrast must have at least `deseq2.min_replicates_per_group` samples (default `2`).
- Each row generates a comparison slug `{group1}_vs_{group2}`; duplicate slugs are rejected.

## 5. `config.yaml` Reference

Rendered from `config/config.yaml` at `init`/`reconfig` time — `WORKDIR`, `HOST`, `ADDITIVES`, `VIRUSES`, `TEMP_DIR`, `REFS_DIR`, `KRAKEN2_DB` placeholders are substituted automatically; everything else is a plain default you can hand-edit after `init`.

::: {.scrollable-table}
| Section | Key(s) | Default | Purpose |
|---|---|---|---|
| top-level | `workdir`, `tempdir`, `samples`, `host`, `additives`, `viruses`, `scriptsdir`, `resourcesdir`, `fastas_gtfs_dir` | templated | Core paths and reference-genome selection |
| top-level | `kraken2_db` | templated | Kraken2 database path (present in config; not consumed by any current rule) |
| `count_matrices.gene_counts` | `enabled`, `flank_size`, `exclude_features` | `true`, `250`, `[rRNA, tRNA, Mt_tRNA]` | Gene-TSS Tn5 count matrix — bin width and excluded biotypes |
| `count_matrices.trna_counts` | `enabled`, `flank_size`, `max_size` | `true`, `100`, `1000` | tRNA gene-midpoint count matrix |
| `count_matrices.pol3_counts` | `enabled`, `flank_size`, `max_size` | `true`, `-1` (full feature span), `1000` | Pol III count matrices (3 gene classes) |
| `count_matrices.repeat_element_counts` | `enabled`, `flank_size`, `max_size` | `true`, `-1`, `1000` | Repeat-element count matrices (6 classes) |
| `count_matrices.viral_genome_counts` | `enabled`, `bin_size` | `true`, `200` | Fixed-width viral genome bin matrix |
| `count_matrices.rrna_counts` | `enabled`, `flank_size` | `true`, `250` | rRNA-TSS count matrix |
| `fastp` | `length_required`, `qualified_quality_phred` | `15`, `15` | Read trimming thresholds |
| `bowtie2_align` | `preset`, `k`, `max_insert`, `no_mixed`, `no_discordant`, `extra_args` | `--very-sensitive`, `20`, `2000`, `true`, `true`, `""` | Alignment parameters |
| `bowtie2_align.filter` | `include_flag`, `exclude_flag`, `mapq`, `exclude_rnames`, `markdup_enable`, `markdup_remove` | `2`, `0x904`, `30`, `[chrM, MT]`, `true`, `true` | Post-alignment BAM filtering |
| `macs2` | `qvalue`, `shift`, `extsize`, `keep_dup`, `extra_args` | `0.01`, `-100`, `200`, `all`, `""` | MACS2 ATAC-seq peak-calling parameters |
| `peakcalling` | `use_host_input`, `use_virus_input` | `false`, `true` | Whether pooled input-control BAMs are passed to peak callers |
| `genrich` | `qvalue`, `remove_dups`, `junctions`, `exclude_chr`, `virus_exclude_chr`, `blacklist`, `virus_blacklist`, `host_fraglen`, `virus_fraglen`, `host_mval`, `virus_mval`, `minlen`, `maxlen`, `extra_args` | `0.05`, `true`, `true`, `chrM,MT`, `""`, `""`, `""`, `200`, `100`, `30`, `5`, `150`, `1000`, `""` | Genrich peak-calling parameters, separately tunable for host vs. virus |
| `postprocessing` | `host_bin_size`, `virus_bin_size`, `host_bw_extra_args`, `virus_bw_extra_args` | `50`, `1`, `""`, `""` | bigWig generation bin sizes |
| `tn5_motif` | `flank_size`, `virus_bin_size`, `host_gene_flank_size`, `trna_gene_flank_size`, `trna_gene_max_size`, `mapq_min`, `dedup`, `exclude_secondary`, `exclude_supplementary`, `fractional_counting`, `logo_format`, `generate_logo` | `10`, `100`, `250`, `100`, `1000`, `0`, `false`, `false`, `true`, `true`, `png`, `false` | Tn5 cut-site extraction / motif-logo config; also read by count-matrix rules for shared filters |
| `ataqv` | `extra_args`, `virus_tss_extension` | `""`, `200` | ataqv QC parameters |
| `deseq2` | `enabled`, `contrasts_tsv`, `outdir`, `alpha`, `lfc_threshold`, `min_replicates_per_group`, `min_total_count`, `min_features_per_matrix`, `size_factor_type`, `fit_type`, `shrink_type`, `p_adjust_method`, `cooks_cutoff`, `independent_filtering`, `vst_blind`, `report_top_n`, `report_label_top_n`, `report_max_table_rows`, `html_self_contained`, `design_factors`, `skip_features` | `false`, `""`, `results/deseq2`, `0.05`, `1.0`, `2`, `10`, `50`, `poscounts`, `parametric`, `apeglm`, `BH`, `true`, `true`, `true`, `50`, `15`, `50000`, `true`, `[group]`, `[tRNA]` | Full DESeq2 differential-accessibility configuration (`design_factors` is currently locked to `["group"]` — additional covariates are not yet implemented) |
| `containers` | one Docker URI per tool | `seqinfomics/*` / `staphb/*` images | Container images for every rule (fastp, fastqc, multiqc, featurecounts, bowtie2, bedToBigBed, py311, pysam, macs2, genrich, deeptools, ataqv, deseq2_report, aws) |
| S3 (top-level) | `push_to_s3`, `s3_pipeline_name`, `s3_sample_set_name`, `s3_aws_credentials_file`, `s3_bucket`, `s3_output_prefix`, `s3_default_storage_class`, `s3_large_file_storage_class` | `false`, `CHROMA`, `""`, ..., `dremel-lab-bucket`, `_HTS`, `GLACIER_IR`, `GLACIER` | S3 archival toggle and destination settings — see [S3 Configuration](s3_configuration.md) |
:::

---

## Summary

### Required Inputs

1. A `samples.tsv` manifest with valid `sampleName`/`groupName`/`role`/`target` values and readable FASTQ paths.
2. A host genome and/or one or more viruses selected via `-g`/`-v` at `init`.

### Optional but Recommended Inputs

1. `contrasts.tsv` + `deseq2.enabled: true` for automated differential-accessibility reporting — see [Outputs — DESeq2 Reports](outputs.md#deseq2-reports).
2. `host_input_pool`/`virus_input_pool` assignments on control samples for background-corrected peak calling.
3. S3 configuration if you want final outputs archived automatically — see [S3 Configuration](s3_configuration.md).
