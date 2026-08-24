# ⚙️ Prerequisites for Running CHROMA2 on Rivanna

CHROMA2 is designed to run on the University of Virginia **Rivanna** (Afton/`shen`) HPC cluster via SLURM, using Apptainer/Singularity containers for every tool. This page covers everything you need before your first run.

---

## 1. Obtain a Rivanna Account

You need an active Rivanna account under the `dremel_lab` allocation. Contact your PI or UVA Research Computing if you don't have one yet.

## 2. Connect via UVA VPN

If you're off-Grounds, connect to the UVA VPN before SSHing into Rivanna.

## 3. Accessing Rivanna from the Command Line (ssh)

```bash
ssh <your_computing_id>@login.hpc.virginia.edu
```

## 4. Load the Mamba Environment

CHROMA2 runs from the lab's shared `pipelines` conda/mamba environment, which provides Snakemake, the `chroma2` wrapper's Python dependencies, and the `SINGULARITY_*`/`DEFAULT_SHARED_SIF_DIR` environment variables the wrapper requires (see [§6](#apptainer-images-and-cache-layout)).

```bash
source ~/.sh_common   # or the lab's equivalent environment bootstrap script
mamba activate pipelines
```

## 5. Verify CHROMA2 Installation

```bash
chroma2 --help
```

Expected output:

```text
##########################################################################################

CHROMA2 Snakemake wrapper (Rivanna)

USAGE (both forms supported):
  chroma2 -w/--workdir=<WORKDIR> -m/--runmode=<RUNMODE>
  chroma2 -w/--workdir <WORKDIR> -m/--runmode <RUNMODE>

Required Arguments:
1.  WORKDIR     : Absolute or relative path to output folder with write permissions.
2.  RUNMODE     : Valid options:
    * init      : initialize workdir
    * dryrun    : dry run snakemake to generate DAG
    * run       : run with slurm (submit head job via sbatch)
    * runlocal  : run on interactive node without sbatch
    ADVANCED RUNMODES (use with caution!!)
    * unlock    : unlock WORKDIR if locked by snakemake
    * reconfig  : recreate config file in WORKDIR (debugging option)
    * reset     : DELETE workdir and re-init (debugging option)
    * recluster : re-copy cluster profile files into WORKDIR/config
    * touch     : touch outputs
    * printbinds: print apptainer bind args (from profile)
    * local     : same as runlocal
    * dry       : same as dryrun

Optional Arguments:
--host|-g       : host genome (e.g. hg38_basic or mm39_basic)                 (--runmode=init only)
--additives|-a  : comma-separated list of additives (e.g. ERCC,BAC16Insert)   (--runmode=init only)
--viruses|-v    : comma-separated list of viruses                       (--runmode=init only)
--manifest|-s   : absolute path to samples.tsv (copied to WORKDIR)      (--runmode=init only)
--sifdir|-i     : path for pre-downloaded SIF images (default: DEFAULT_SHARED_SIF_DIR)
--help|-h       : print this help
--version|-V    : print version and exit

Example:
  chroma2 -w /my/output/folder -m init
  chroma2 --workdir /my/output/folder --runmode dryrun
  chroma2 --workdir=/my/output/folder --runmode=run

##########################################################################################
```

Check the installed version at any time with `chroma2 --version`.

## 6. Apptainer Images and Cache Layout {#apptainer-images-and-cache-layout}

Every rule in CHROMA2 runs inside an Apptainer (Singularity) container (`seqinfomics/*` and `staphb/*` images pulled from Docker Hub — see the `containers:` block in `config.yaml`). Two things need to be true before you can run the pipeline:

1. **A shared, pre-pulled SIF image cache must exist.** By default the wrapper reads its location from the `DEFAULT_SHARED_SIF_DIR` environment variable, which is set automatically by the `pipelines` mamba environment. You can override it per-invocation with `-i/--sifdir <path>`.
2. **Four Singularity environment variables must be set**: `SINGULARITY_CACHEDIR`, `SINGULARITY_TMPDIR`, `SINGULARITY_PULLFOLDER`, and `DEFAULT_SHARED_SIF_DIR`. These are also provided by the `pipelines` environment. If any of them is missing, `chroma2 run`/`runlocal` will exit immediately with an error telling you which variable is unset — if that happens after activating `pipelines`, email `cud2td@virginia.edu`.

Every per-rule SLURM job additionally sets up its own scratch-local Apptainer cache/tmp/pull directory tree under `$SCRATCH/singularity/{cache,tmp,sif}` at job start (see `config/rivanna/jobscript.sh`), falling back to the shared SIF directory if a compute node can't reach it.

To inspect the exact bind-mount arguments a workdir's cluster profile will pass to Apptainer, run:

```bash
chroma2 -w <WORKDIR> -m printbinds
```

## 7. Reference Bundle Created Inside Each Work Directory

At `init` time, CHROMA2 does **not** build reference indices — that happens on the first real `dryrun`/`run`/`runlocal` invocation, when Snakemake's `init.smk` logic concatenates the selected host/additive/virus FASTAs and GTFs and builds the Bowtie2 index. The resulting reference bundle lives under `$WORKDIR/ref/`:

```text
$WORKDIR/ref/
├── ref.fa                  # concatenated host + additive + virus FASTA
├── ref.fa.regions          # combined contig list
├── ref.fa.host.regions     # host-only contig list
├── ref.fa.virus.regions    # virus-only contig list (per configured virus)
├── ref.gtf                 # concatenated annotation GTF
├── ref.fa.{1,2,3,4}.bt2    # Bowtie2 index
├── ref.fa.rev.{1,2}.bt2
├── ref.tss.host.bed        # host TSS positions (from ref.gtf `transcript` features)
├── ref.tss.<virus>.bed     # per-virus TSS positions
└── ...                     # chrom-sizes / autosomes files used by ataqv
```

Source FASTA/GTF paths per host/virus/additive come from `reference_gtf`, `trnas_gtf`, `chrr_gtf`, `pol3_gtf`, and `repeat_elements_gtf` in `config.yaml`, resolved against the shared reference data directory (`refs_dir`, default `/project/dremel_lab/workflows/reference_data/fasta_gtf/`). CHROMA2 validates that every required GTF for your selected host/viruses exists and is readable _before_ running any rules, and will exit with a clear report of any missing files.

## 8. Understanding Snakemake Dry-Run Output {#8-understanding-snakemake-dry-run-output}

Before submitting a real run, always dry-run first:

```bash
chroma2 -w <WORKDIR> -m dryrun
```

This calls `snakemake --dry-run` under the hood and logs to `$WORKDIR/logs/dryrun.<timestamp>.log`. Read the log for:

- **The job list and counts** — one line per rule showing how many jobs of that rule will run (e.g. `macs2_atac_callpeak_host: 4`). Sanity-check this against your sample count and enabled features (peak calling ×2 callers ×2 organisms, count matrices ×6 types, etc.).
- **Missing input errors** — if a `MissingInputException` appears, check `samples.tsv` paths and `config.yaml` GTF resolution first.
- **A locked working directory warning** — if `.snakemake/lock` exists from an interrupted previous run, `dryrun` will warn; use `chroma2 -w <WORKDIR> -m unlock` to clear it before proceeding.

## 9. Other Runmodes (Advanced) {#other-runmodes}

These runmodes exist for debugging and recovery — use with caution:

- **`unlock`** — releases a stale Snakemake lock (`.snakemake/lock`) left behind by an interrupted run. Refuses to run if the pipeline is currently marked as actively running.
- **`reconfig`** — regenerates `$WORKDIR/config.yaml` from the pipeline's template, **overwriting any manual edits you made to it**. Only the `host`/`additives`/`viruses`/paths get re-templated; use it to pick up a new host/virus selection without a full `reset`.
- **`recluster`** — re-copies the `config/` cluster-profile tree (SLURM resource settings, job script) from the pipeline install into `$WORKDIR/config`, without touching `config.yaml` or `samples.tsv`. Use this after a pipeline upgrade that changed per-rule SLURM resource defaults.
- **`reset`** — **deletes the entire working directory** and re-runs `init`. Irreversible; only use on workdirs you're certain have no results worth keeping.
- **`touch`** — runs `snakemake --touch` to mark existing outputs as up to date without re-running anything (useful after restoring outputs from backup or fixing a false timestamp mismatch).
- **`printbinds`** — prints the Apptainer bind-mount arguments configured in the workdir's cluster profile, with no side effects.

## 10. Cluster Profile

CHROMA2 auto-detects the SLURM cluster name via `scontrol show config`. On Rivanna (`ClusterName: shen`), it selects the `rivanna` profile (`config/rivanna/`, containing `config.yaml` with per-rule SLURM resource overrides and `jobscript.sh`, the template used for every per-rule job submission) and the `standard` partition, and prepends the lab's shared script/env `bin` directories to `PATH`. On any other cluster, the wrapper prints a warning and you'll need to hand-edit `config/` and `config.yaml` for compatibility — CHROMA2 currently ships only the `rivanna` profile.

---

**Next:** continue to [Usage](usage.md) to initialize and run your first CHROMA2 analysis.
