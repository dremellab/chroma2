#!/bin/bash
#SBATCH --cpus-per-task={threads}
#SBATCH --mem={resources.mem_mb}
#SBATCH --time={resources.runtime}
#SBATCH --account=dremel_lab
#SBATCH --output=logs/{rule}.{jobid}.slurm-%j.out
#SBATCH --error=logs/{rule}.{jobid}.slurm-%j.err

source ~/.bashrc
export PROFILE={profile}

SCRATCH_ROOT="${SCRATCH:-/scratch/$USER}"
if [[ ! -d "$SCRATCH_ROOT" ]]; then
  mkdir -p "$SCRATCH_ROOT" 2>/dev/null || SCRATCH_ROOT="$PWD/.harold_runtime"
fi

SINGULARITY_RUNTIME_ROOT="${SCRATCH_ROOT}/singularity"
APPTAINER_CACHE_DIR="${SINGULARITY_RUNTIME_ROOT}/cache"
APPTAINER_TMP_DIR="${SINGULARITY_RUNTIME_ROOT}/tmp"
FALLBACK_SIF_DIR="${SINGULARITY_RUNTIME_ROOT}/sif"

mkdir -p "$APPTAINER_CACHE_DIR" "$APPTAINER_TMP_DIR"

if [[ -n "${APPTAINER_SIF_DIR:-}" ]]; then
  if [[ ! -d "$APPTAINER_SIF_DIR" ]]; then
    echo "WARNING: APPTAINER_SIF_DIR='${APPTAINER_SIF_DIR}' is unavailable on this node; caching to ${FALLBACK_SIF_DIR} instead." >&2
    APPTAINER_SIF_DIR="$FALLBACK_SIF_DIR"
    mkdir -p "$APPTAINER_SIF_DIR"
  fi
else
  APPTAINER_SIF_DIR="$FALLBACK_SIF_DIR"
  mkdir -p "$APPTAINER_SIF_DIR"
fi

export APPTAINER_CACHEDIR="$APPTAINER_CACHE_DIR"
export SINGULARITY_CACHEDIR="$APPTAINER_CACHE_DIR"
export APPTAINER_TMPDIR="$APPTAINER_TMP_DIR"
export SINGULARITY_TMPDIR="$APPTAINER_TMP_DIR"
export APPTAINER_PULLDIR="$APPTAINER_SIF_DIR"
export TMPDIR="$APPTAINER_TMP_DIR"

{exec_job}
