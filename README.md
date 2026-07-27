# CHROMA2

**CHROMA2** is a modular, reproducible **ATAC-seq analysis pipeline** built with Snakemake, purpose-built for **host + virus co-analysis**: it aligns reads against a composite reference (a host genome such as `hg38`/`mm39`/`hs1`, optional additives like `ERCC`/`BAC16Insert`, and one or more viral genomes) in a single coherent run, so chromatin accessibility can be examined across host and pathogen genomes together.

The pipeline runs dual peak callers (MACS2 and Genrich, each with optional pooled input-control support), generates six types of Tn5 insertion count matrices (gene, tRNA, Pol III, repeat elements, viral, and rRNA), performs DESeq2 differential-accessibility reporting per contrast, and produces a consolidated MultiQC report. It is driven end-to-end by the `chroma2` wrapper script for SLURM/HPC execution on Rivanna, with optional S3 archival of final outputs.

📖 **[Read the full documentation](https://dremellab.github.io/chroma2/dev/docs/index.html)** for prerequisites, usage, pipeline architecture, inputs/outputs, and more.
