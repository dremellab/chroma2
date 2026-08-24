# 🧬 CHROMA2 — ATAC-seq Host/Virus Co-analysis Pipeline

**CHROMA2** is a modular and reproducible **ATAC-seq analysis pipeline**, designed to characterize chromatin accessibility across a host genome and one or more co-infecting or integrated viral genomes in a single, coherent analysis.

CHROMA2 automates the process from raw FASTQ files to peak calls, Tn5 insertion count matrices, and differential-accessibility reports, aggregating extensive quality-control metrics into an interactive MultiQC dashboard. It is designed for **virus-host chromatin accessibility studies**, timecourse infection experiments, and any ATAC-seq experiment that needs a composite host+additive+virus reference in a single run.

---

## 🚀 Key Features

- **CHROMA2 supports composite host/additive/virus references.** It builds a single Bowtie2 index that concatenates a host genome (e.g. `hg38`, `mm39`, `hs1`), optional additives (`ERCC`, `BAC16Insert`), and one or more viral genomes, so accessibility can be quantified consistently across all of them in one alignment pass.

- **CHROMA2 runs dual, independently tunable peak callers.** MACS2 (ATAC-centric `--nomodel --shift -100 --extsize 200` cut-site model) and Genrich (native ATAC-seq mode on qname-sorted, unfiltered BAMs) are both run for every sample, each with separately configurable parameters for host vs. viral genomes.

- **CHROMA2 supports pooled input/IgG-style controls.** Control samples in the sample manifest can be grouped into named host and virus "input pools" via `role`/`target`/`host_input_pool`/`virus_input_pool` columns; the pooled control BAM is automatically passed as `-c` to MACS2 and Genrich for the case samples that reference it.

- **CHROMA2 is powered by Snakemake for workflow orchestration.** Each stage — trimming, alignment, filtering, peak calling, quantification, and reporting — is managed through clearly defined, container-based rules that ensure reproducibility and modular execution.

- **CHROMA2 generates six types of Tn5 insertion count matrices.** Gene-TSS, tRNA, Pol III (3 gene classes), repeat elements (6 classes: SINE/Alu, SINE/MIR, LINE/L1, LINE/L2, LTR, other), viral genome bins, and rRNA — all built from precise Tn5 insertion sites with fractional (NH-weighted) counting for multi-mapping reads, so repeat-rich and multi-copy loci are counted fairly.

- **CHROMA2 can run DESeq2 differential accessibility directly.** Given a `contrasts.tsv` manifest, CHROMA2 runs a DESeq2 analysis per contrast against every enabled count-matrix type, producing a self-contained HTML report per comparison with volcano plots — no separate manual invocation required.

- **CHROMA2 produces a comprehensive MultiQC report.** This report aggregates FastQC, fastp, alignment/filtering statistics, host-vs-virus read distribution, fragment-size distributions, FRiP scores, ataqv accessibility metrics, peak-size distributions, genome coverage breadth, and Tn5 count-matrix summaries into one interactive dashboard.

- **CHROMA2 includes seamless HPC integration.** The pipeline runs on clusters like Rivanna through a single wrapper script (`chroma2`) that supports initialization, dry-run validation, and SLURM-based execution, with live progress tracking during a run.

- **CHROMA2 tracks and reports its own run state.** Every run writes `pipeline.{running,completed,failed,canceled}` markers and a `pipeline.status.json` sidecar to the working directory, with structured event logs (`logs/events.log`/`.jsonl`), so status can be checked without parsing raw Snakemake output.

- **CHROMA2 supports optional S3 cloud deposition.** Final outputs (configs, QC reports, count matrices, DESeq2 reports) can be automatically transferred to Amazon S3 for cloud storage and collaboration, with configurable storage classes to control archival cost.

---

## 🧬 Workflow Overview

An ATAC-seq experiment probes which regions of the genome are accessible to the Tn5 transposase — a proxy for open chromatin and, indirectly, for regulatory activity. In a host/virus co-infection or co-culture context, the same question applies to the pathogen genome: which regions of a viral episome or integrated provirus are accessible at a given timepoint? CHROMA2 automates the full path from raw reads to interpretable accessibility signal for both compartments simultaneously.

Reads are trimmed with fastp before alignment against a composite reference assembled from the selected host, additive, and viral genomes. Bowtie2 alignments are then carried through a multi-stage filtering pipeline — proper-pair filtering, fixmate, duplicate marking, MAPQ and blacklist-contig filtering — to produce a clean alignment set, which is then split into host-only and virus-only BAMs for all downstream, organism-specific analysis.

From the filtered alignments, CHROMA2 calls peaks with both MACS2 and Genrich, computes genome-wide accessibility tracks (bigWig and BigBed), and extracts precise Tn5 insertion sites to build count matrices across genes, tRNAs, Pol III loci, repeat element classes, viral genome bins, and rRNA loci. These count matrices feed directly into per-contrast DESeq2 differential-accessibility reports, while ataqv and a battery of custom QC metrics feed into a single MultiQC dashboard — giving a complete picture of both signal quality and biological result from one pipeline run.

---
