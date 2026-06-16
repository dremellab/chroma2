## Problem

`sbatch` job submission is failing for user `cud2td` with:

```text
sbatch: error: slurm_job_submit: user specified partition: standard
sbatch: error: slurm_job_submit: resulting partitions standard-rivanna,standard-afton,standard-rivanna-largemem,standard-afton-largemem
sbatch: error: Batch job submission failed: Invalid account or account/partition combination specified
```

## Batch script used

```bash
#!/bin/bash
#SBATCH --job-name=20260225_KSclinicalRNASeq
#SBATCH --output=/scratch/cud2td/20260225_KSclinicalRNASeq/20260225_KSclinicalRNASeq.out
#SBATCH --error=/scratch/cud2td/20260225_KSclinicalRNASeq/20260225_KSclinicalRNASeq.err
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=40G
#SBATCH --partition=standard
#SBATCH --account=dremel_lab

source /project/dremel_lab/scripts/.sh_common
mamba activate pipelines

harold --workdir=/scratch/cud2td/20260225_KSclinicalRNASeq --runmode=run
```

## What I checked

### `sinfo`

```text
PARTITION                 AVAIL  TIMELIMIT  NODES  STATE NODELIST
standard                     up 7-00:00:00      1  down* udc-aw36-22c1
standard                     up 7-00:00:00     10  drain udc-aw36-20c1,udc-aw36-21c[0-1],udc-aw36-22c0,udc-aw36-23c[0-1],udc-ba02-3c[0-1],udc-ba02-4c[0-1]
standard                     up 7-00:00:00      1   idle udc-aw36-20c0
parallel                     up 3-00:00:00      4  drain udc-ba02-3c[0-1],udc-ba02-4c[0-1]
gpu                          up 3-00:00:00      1  inval udc-ba27-26
gpu                          up 3-00:00:00      1  drain udc-an25-16
standard-rivanna             up 7-00:00:00      2  drain udc-aw36-21c1,udc-aw36-22c0
standard-rivanna-largemem    up 7-00:00:00      1  down* udc-aw36-22c1
standard-afton               up 7-00:00:00      5  drain udc-aw36-21c0,udc-ba02-3c[0-1],udc-ba02-4c[0-1]
standard-afton               up 7-00:00:00      1   idle udc-aw36-20c0
standard-afton-largemem      up 7-00:00:00      1  drain udc-aw36-20c1
gpu-a6000                    up 3-00:00:00      1  drain udc-an25-16
gpu-v100                     up 3-00:00:00      1  inval udc-ba27-26
interactive-rivanna          up   12:00:00      0    n/a
interactive-afton            up   12:00:00      0    n/a
gpu-a100-40                  up 3-00:00:00      1  inval udc-an26-7
interactive-rtx2080          up 3-00:00:00      0    n/a
interactive-rtx3090          up 3-00:00:00      0    n/a
```

### `scontrol show partition`

Relevant sections show that the partitions appear to allow all accounts:

```text
PartitionName=standard
   AllowGroups=ALL AllowAccounts=ALL AllowQos=ALL
   QoS=standard
```

```text
PartitionName=standard-rivanna
   AllowGroups=ALL AllowAccounts=ALL AllowQos=ALL
   QoS=standard
```

```text
PartitionName=standard-afton
   AllowGroups=ALL AllowAccounts=ALL AllowQos=ALL
   QoS=standard
```

### `sshare -A dremel_lab`

```text
Account                    User  RawShares  NormShares    RawUsage  EffectvUsage  FairShare
-------------------- ---------- ---------- ----------- ----------- ------------- ----------
```

### `sshare -A dremel_lab -u cud2td`

```text
Account                    User  RawShares  NormShares    RawUsage  EffectvUsage  FairShare
-------------------- ---------- ---------- ----------- ----------- ------------- ----------
```

### `sshare -u cud2td`

```text
Account                    User  RawShares  NormShares    RawUsage  EffectvUsage  FairShare
-------------------- ---------- ---------- ----------- ----------- ------------- ----------
root                                          0.000000     5086517      1.000000
 hpc_build                               1    0.500000     5086517      1.000000
```

### `sacctmgr`

```bash
sacctmgr -nP show user cud2td withassoc
```

No output.

```bash
sacctmgr -nP show assoc where user=cud2td format=cluster,account,partition,qos,defaultqos
```

No output.

### `id`

```text
uid=15786556(cud2td) gid=100(users) groups=100(users),508(parallel),526(standard),567(gpu),628(interactive),660(largemem),647008(uvaRaveAppArmor2),647167(UV_Contractor),647538(dremel_lab)
```

## Interpretation

This looks like a SLURM accounting association problem rather than a partition-definition problem:

- user `cud2td` is in Unix group `dremel_lab`
- the `dremel_lab` account appears to exist
- but `sshare -A dremel_lab -u cud2td` is empty
- and `sacctmgr` shows no associations for `cud2td`

## Request

Please check whether user `cud2td` is associated with the correct SLURM account for batch submission on the `standard` partitions, and either:

- add `cud2td` to the correct SLURM account (possibly `dremel_lab`), or
- tell me which SLURM account I should use for `sbatch`

## Date

Observed on 2026-03-18.
