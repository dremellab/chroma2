# 🧭 Using CHROMA2 on Rivanna

Running CHROMA2 is a three-step process: **initialize** a working directory, **dry-run** it to validate the plan, then **execute** it. This page walks through each step, plus the optional DESeq2 and S3 stages.

Before you start, make sure you have:

- A **sample manifest** (`samples.tsv`) — see [Inputs](inputs.md#sample-manifest) for the required 9-column format.
- Decided on your **host genome**, **additives**, and **viruses** — see [Inputs](inputs.md#reference-configuration) for the supported values.

---

## Step 1: Initialization (`--runmode=init`)

`init` creates the working directory, renders `config.yaml` from the pipeline's template, and copies in your sample manifest and (if present) `contrasts.tsv`.

**Required arguments:**

- `-w/--workdir` — path to the new working directory (must not already exist)
- `-m/--runmode init`

**Common optional arguments** (only meaningful with `init`/`reconfig`):

- `-g/--host` — host genome, default `hg38_basic` (also available: `hg38`, `mm39`, `hs1`)
- `-a/--additives` — comma-separated additive contigs, e.g. `ERCC,BAC16Insert` (default: none)
- `-v/--viruses` — comma-separated viral accessions, default `KT899744.1`
- `-s/--manifest` — absolute path to your `samples.tsv`, copied into the workdir

```bash
chroma2 -w /scratch/$USER/analysis/my_atac_run \
  -m init \
  -g hg38_basic \
  -v KT899744.1 \
  -s /scratch/$USER/rawdata/my_atac_run/samples.tsv
```

After `init`, review and, if needed, hand-edit `$WORKDIR/config.yaml` and `$WORKDIR/samples.tsv` before proceeding — `init` only templates paths and genome selection; it does not validate the manifest contents (that happens when Snakemake first parses it during `dryrun`/`run`).

## Step 2: Dry Run (`--runmode=dryrun`, alias `dry`)

```bash
chroma2 -w /scratch/$USER/analysis/my_atac_run -m dryrun
```

Runs `snakemake --dry-run` against the workdir's cluster profile and logs to `$WORKDIR/logs/dryrun.<timestamp>.log`. This is the point where manifest-column validation, GTF-availability checks, and DESeq2 contrast validation (if enabled) all run — fix any errors here before submitting a real job. See [Prerequisites §8](prereq.md#8-understanding-snakemake-dry-run-output) for how to read the output.

## Step 3: Execution

CHROMA2 supports two execution modes:

### `--runmode=run` — submit to SLURM

```bash
chroma2 -w /scratch/$USER/analysis/my_atac_run -m run
```

This generates `$WORKDIR/run_head_job.sbatch` (a driver job that itself invokes `snakemake --profile config/rivanna`, requesting 2 CPUs / 40G / 3 days by default, account/partition read from `config.yaml`'s `slurm_account`/`slurm_partition` keys — default `dremel_lab`/`standard`) and submits it with `sbatch --parsable`. Snakemake then submits one SLURM job per pipeline rule under that umbrella, per the resource settings in `config/rivanna/config.yaml`.

`run` refuses to start if `$WORKDIR/pipeline.running` already exists (i.e. a run is already active) — use `unlock` first if a previous run was interrupted uncleanly.

### `--runmode=runlocal` (alias `local`) — run on an interactive node

```bash
salloc --account=dremel_lab --partition=standard --cpus-per-task=4 --mem=32G --time=8:00:00
chroma2 -w /scratch/$USER/analysis/my_atac_run -m runlocal
```

Runs Snakemake directly in the foreground (still submitting one SLURM job per rule via the same profile) — requires an active `$SLURM_JOB_ID` (i.e. you must be inside an `salloc`/interactive session first). Streams output to `$WORKDIR/logs/runlocal.<timestamp>.log`. `Ctrl-C`/`SIGTERM` is caught and marks the run `canceled` cleanly rather than leaving stale state.

### Tracking Run State

Every `run`/`runlocal` invocation writes and updates:

- `$WORKDIR/pipeline.running` — present while the pipeline is actively executing; contains a live "N/M steps complete" progress snapshot, refreshed roughly every 60 seconds.
- `$WORKDIR/pipeline.completed` / `pipeline.failed` / `pipeline.canceled` — written on exit, whichever applies. The `failed` marker includes a digest of which rule(s) failed.
- `$WORKDIR/pipeline.status.json` — machine-readable sidecar mirroring the same state, for scripts/dashboards that want to poll run status without parsing log files.
- `$WORKDIR/logs/events.log` and `logs/events.jsonl` — structured event log (human-readable and JSONL) of every state transition.

You can check whether a run is active, completed, or failed at any time just by looking at which `pipeline.*` file exists in the workdir — no need to `squeue`/parse Snakemake logs.

## Step 4 (Optional): DESeq2 Differential Accessibility

If `$WORKDIR/contrasts.tsv` exists and `deseq2.enabled: true` in `config.yaml`, CHROMA2 runs a DESeq2 differential-accessibility report for every contrast row, against every enabled Tn5 count-matrix type. There is no separate CLI flag — `contrasts.tsv` is copied into the workdir automatically at `init` (from `config/contrasts.tsv` in the pipeline install) and just needs to match your experiment's group names. See [Inputs](inputs.md#contrasts-manifest) for the format and [Outputs](outputs.md#deseq2-reports) for the report contents.

## Step 5 (Optional): S3 Deposition

If `push_to_s3: true` and `s3_sample_set_name` is set (non-empty) in `config.yaml`, CHROMA2 transfers final outputs (configs, QC reports, count matrices, DESeq2 reports) to S3 as the last step of the run. See [S3 Configuration](s3_configuration.md) for full setup instructions.

---

## Monitoring SLURM Jobs

```bash
squeue -u $USER
```

Each pipeline rule's SLURM stdout/stderr lands under `$WORKDIR/logs/`, named by rule and job ID (per `config/rivanna/config.yaml`'s `slurm-logdir: logs` setting). To find the log for a specific failed rule:

```bash
find $WORKDIR/logs -name "*<rule_name>*" -newer $WORKDIR/pipeline.running
```

## Other Runmodes

`unlock`, `reconfig`, `reset`, `recluster`, `touch`, and `printbinds` are all documented in [Prerequisites §9](prereq.md#other-runmodes) — use them only when you understand the tradeoffs, since several are destructive or overwrite manual edits.

---

## Summary

A typical CHROMA2 session is: `init` → edit `config.yaml`/`samples.tsv` as needed → `dryrun` (fix any validation errors) → `run` (or `runlocal` on an interactive node) → monitor via `pipeline.*` markers and `squeue` → once complete, explore [Outputs](outputs.md) and, if enabled, the DESeq2 and MultiQC reports.
