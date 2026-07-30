# ❓ Help and Support

## 🧩 Common Issues

1. **`ValueError` about a missing manifest column.** `samples.tsv` must have exactly these 9 tab-delimited columns: `sampleName, groupName, batch, path_to_R1_fastq, path_to_R2_fastq, role, target, host_input_pool, virus_input_pool`. See [Inputs — Sample Manifest](inputs.md#sample-manifest) for the full spec, and check for stray spaces or a copy-pasted header with the wrong delimiter.

2. **"Unsupported host/virus" or a missing-GTF error at dry-run time.** CHROMA2 validates that every required GTF for your selected `-g/--host` and `-v/--viruses` exists before running any rules. Confirm your host name matches a key in `reference_gtf` (`mm39`, `hg38`, `hg38_basic`, `hs1`) and that the reference data directory (`fastas_gtfs_dir`) is reachable. See [Inputs — Reference Configuration](inputs.md#reference-configuration).

3. **`role`/`target` validation error.** `role` must be exactly `case` or `control`; `target` must be exactly `host`, `virus`, or `both` — check for typos or mixed case.

4. **SLURM submission failures (`sbatch` errors, account/association errors).** Confirm your `dremel_lab` SLURM account association is active (`sacctmgr show associations user=$USER`), and that `slurm_account`/`slurm_partition` in `config.yaml` match a partition you have access to. If your account association itself is broken, that's a Research Computing ticket, not a CHROMA2 issue.

5. **`chroma2` exits immediately complaining about `SINGULARITY_*`/`DEFAULT_SHARED_SIF_DIR`.** These must be set by the `pipelines` mamba environment — re-activate it (`mamba activate pipelines`) and retry. If the error persists, email `cud2td@virginia.edu`.

6. **DESeq2 contrast validation errors** ("group not found" / "insufficient replicates"). Every `group1`/`group2` value in `contrasts.tsv` must exactly match a `groupName` in `samples.tsv`, and each group needs at least `deseq2.min_replicates_per_group` (default `2`) samples. See [Inputs — Contrasts Manifest](inputs.md#contrasts-manifest).

7. **A previous run left `.snakemake/lock` behind and `dryrun`/`run` refuses to proceed.** Run `chroma2 -w <WORKDIR> -m unlock` first (only if no run is currently active).

## 💬 Getting Help

Contact `cud2td@virginia.edu` for pipeline support. When reporting an issue, please include:

- The exact `chroma2` command you ran (workdir, runmode, and any flags).
- The relevant section of `$WORKDIR/pipeline.failed` (if present) and/or the SLURM log for the failing rule under `$WORKDIR/logs/`.
- Your `config.yaml` and `samples.tsv` (or the relevant excerpt) if the issue looks input-related.

## 🧠 Additional Resources

- [Snakemake documentation](https://snakemake.readthedocs.io/)
- [MultiQC documentation](https://multiqc.info/docs/)
- [MACS2 documentation](https://github.com/macs3-project/MACS)
- [Genrich documentation](https://github.com/jsh58/Genrich)
- [ataqv documentation](https://github.com/ParkerLab/ataqv)
- [UVA Research Computing / Rivanna docs](https://www.rc.virginia.edu/userinfo/hpc/)
