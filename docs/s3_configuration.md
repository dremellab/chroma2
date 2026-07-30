# ☁️ S3 Output Deposition

CHROMA2 can optionally archive a curated set of final outputs to Amazon S3 as the last step of a run, using the lab's shared Globus-managed S3 credentials.

## Prerequisites

### S3 Bucket Access

Your S3 push destination is the lab-managed bucket (`dremel-lab-bucket` by default). Confirm with your PI/lab admin that you have write access via the shared credentials before enabling `push_to_s3`.

### S3 Credentials (Lab-Managed)

CHROMA2 does not manage AWS credentials itself — it reads them from a shared credentials file (`s3_aws_credentials_file` in `config.yaml`, default `/project/dremel_lab/scripts/aws/credentials`) and authenticates using the `s3-globus-user` AWS profile. If your transfer fails with an authentication error, check with lab IT whether this credentials file needs rotation.

## Configuration in `config.yaml`

### S3 Settings Block

```yaml
push_to_s3: false # master on/off switch
s3_pipeline_name: "CHROMA" # pipeline name segment in the S3 path
s3_sample_set_name: "" # REQUIRED (non-empty) if push_to_s3 is true
s3_aws_credentials_file: "/project/dremel_lab/scripts/aws/credentials"
s3_bucket: "dremel-lab-bucket"
s3_output_prefix: "_HTS"
s3_default_storage_class: "GLACIER_IR"
s3_large_file_storage_class: "GLACIER"
```

### Configuration Details

::: {.scrollable-table}
| Key | Type | Default | Required? | Description |
|---|---|---|---|---|
| `push_to_s3` | bool | `false` | — | Master toggle; when `false` the transfer rule is a no-op (it always touches `.s3_transfer.done` so it never blocks the pipeline's final target) |
| `s3_pipeline_name` | string | `"CHROMA"` | when enabled | First path segment after the S3 prefix — note this defaults to `CHROMA`, not `CHROMA2` |
| `s3_sample_set_name` | string | `""` | **yes, when enabled** | Pipeline errors out if this is empty and `push_to_s3: true` |
| `s3_aws_credentials_file` | path | lab-managed path | when enabled | Passed as `AWS_SHARED_CREDENTIALS_FILE` |
| `s3_bucket` | string | `dremel-lab-bucket` | when enabled | Destination S3 bucket |
| `s3_output_prefix` | string | `_HTS` | when enabled | Top-level S3 key prefix |
| `s3_default_storage_class` | string | `GLACIER_IR` | when enabled | Storage class for all files except BAM/BAI |
| `s3_large_file_storage_class` | string | `GLACIER` | when enabled | Storage class specifically for `.bam`/`.bai` files |
:::

## S3 Path Hierarchy

```text
s3://{bucket}/{s3_output_prefix}/{s3_pipeline_name}/{s3_sample_set_name}/{category}/{file}
```

Worked example (`bucket=dremel-lab-bucket`, `s3_output_prefix=_HTS`, `s3_pipeline_name=CHROMA`, `s3_sample_set_name=my_atac_run`):

```text
s3://dremel-lab-bucket/_HTS/CHROMA/my_atac_run/
├── config/
│   ├── samples.tsv
│   ├── config.yaml
│   └── contrasts.tsv
├── qc/
│   ├── alignment_summary.tsv
│   ├── multiqc_report.html
│   ├── multiqc_data/...
│   ├── custom/...
│   └── ataqv/...             # aggregate ataqv reports + per-sample .json files
├── bams/
│   └── {sample}.aligned.final.bam(.bai)
├── bigwigs/
│   └── {sample}.{host,virus}.bw
├── bigbeds/
│   └── *.bb                  # (ref/ TSS bigbeds excluded)
├── peaks/
│   └── *.narrowPeak.gz, *.summits.bed.gz
├── tn5_counts/
│   └── ...                   # contents of results/count_matrices/ — per-sample counts,
│                              # per-category *_count_matrix.tsv, and final_matrices/
└── deseq2/
    └── {group1}_vs_{group2}/...
```

## Files Transferred to S3

::: {.scrollable-table}
| Category | Source | Destination prefix | Storage Class |
|---|---|---|---|
| Config | `samples.tsv`, `config.yaml`, `contrasts.tsv` | `config/` | default |
| Alignment QC | `results/alignmentqc/idxstats_summary.tsv` | `qc/alignment_summary.tsv` | default |
| MultiQC | `results/multiqc/multiqc_report.html`, `results/multiqc/multiqc_data/` | `qc/` | default |
| Custom QC data | `results/multiqc_extra_data/custom/` | `qc/custom/` | default |
| ataqv reports | `results/alignmentqc/ataqv/`, plus any other `*.json` (excluding fastp report JSONs) | `qc/ataqv/` | default |
| Final BAMs | `*.aligned.final.bam(.bai)` | `bams/` | **large-file class** |
| Coverage tracks | `*.bw` | `bigwigs/` | default |
| Browser tracks | `*.bb` (excludes files under `ref/`) | `bigbeds/` | default |
| Peak calls | `*.narrowPeak.gz`, `*.summits.bed.gz` | `peaks/` | default |
| Tn5 count matrices | `results/count_matrices/` | `tn5_counts/` | default |
| DESeq2 reports | `results/deseq2/` | `deseq2/` | default |
:::

## Enabling and Running S3 Transfer

### Step 1: Configure

Edit `$WORKDIR/config.yaml`:

```yaml
push_to_s3: true
s3_sample_set_name: "my_atac_run"
```

### Step 2: Dry Run

```bash
chroma2 -w $WORKDIR -m dryrun
```

Confirm `s3_transfer_if_enabled` appears in the job list and depends on the full set of expected final outputs (it will not run until everything upstream — alignment QC, all enabled count matrices, DESeq2 reports if enabled, MultiQC — is complete).

### Step 3: Run

```bash
chroma2 -w $WORKDIR -m run
```

Verify after completion:

```bash
aws s3 ls s3://dremel-lab-bucket/_HTS/CHROMA/my_atac_run/ --recursive
```

## Troubleshooting S3 Transfer

### Transfer Failed or Hung

- Check `$WORKDIR/pipeline.failed` and the SLURM log for the `s3_transfer_if_enabled` rule.
- Common issues: expired/rotated shared credentials, missing bucket write permission, or `s3_sample_set_name` left empty (the rule errors out immediately in that case).

### Partial Transfer (Some Files Missing)

- Re-check the [Files Transferred](#files-transferred-to-s3) table above — only the listed categories are transferred; anything else is expected to be absent from S3.
- You can re-run the transfer manually against an already-completed workdir:

```bash
python workflow/scripts/s3_transfer_chroma2.py \
  --workdir $WORKDIR \
  --pipeline-name CHROMA \
  --sample-set-name my_atac_run \
  --bucket dremel-lab-bucket \
  --s3-prefix _HTS \
  --storage-class GLACIER_IR \
  --large-file-storage-class GLACIER
```

Add `--dry-run` to preview the exact `aws s3 cp` commands without executing them.

### Cost Overruns

- BAM/BAI files use the more deeply archived `GLACIER` class by default (vs. `GLACIER_IR` for everything else) specifically to control cost for the largest files.
- Review `s3_default_storage_class`/`s3_large_file_storage_class` in `config.yaml` if your access pattern needs faster retrieval (e.g. `STANDARD_IA`) at higher storage cost.

## Best Practices

1. Use a descriptive, unique `s3_sample_set_name` per run — it's the only thing disambiguating separate runs under the same pipeline/bucket prefix.
2. Always `dryrun` after changing S3 settings to confirm the transfer rule is picked up before submitting a real run.
3. Treat S3-archived outputs as cold storage — `GLACIER`/`GLACIER_IR` retrieval is not instantaneous; don't rely on this as your primary/working copy of results.

## See Also

- [Usage — Step 5: S3 Deposition](usage.md)
- [Outputs — S3 Transfer Sentinel](outputs.md#s3-transfer-sentinel-optional)
