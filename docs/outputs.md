# 📊 CHROMA2 Output Files

## Output Directory Structure

```text
$WORKDIR/
├── config.yaml                    # rendered pipeline config
├── samples.tsv                    # sample manifest
├── contrasts.tsv                  # DESeq2 contrasts manifest (if present)
├── pipeline.{running,completed,failed,canceled}  # run-state marker (whichever applies)
├── pipeline.status.json           # machine-readable run-state sidecar
├── logs/                          # per-rule SLURM logs, events.log/.jsonl, dryrun/run logs
├── ref/                           # composite reference bundle (see prereq.md §7)
└── results/
    ├── preprocess/                # fastp trimmed reads + HTML/JSON reports
    ├── initqc/                    # FastQC on raw + trimmed reads
    ├── align/                     # sorted/filtered/qname BAMs + idxstats/flagstat
    ├── alignmentqc/
    │   ├── idxstats_summary.tsv   # aggregate alignment stats across all filter stages
    │   └── ataqv/                 # per-sample ataqv JSON + mkarv aggregate HTML reports
    ├── inputs/                    # pooled input-control BAMs (host + per-virus)
    ├── peaks/                     # MACS2 + Genrich narrowPeak/summits (gzipped)
    ├── postprocess/                # host/virus-split filtered BAMs
    ├── bigwig/                    # per-organism coverage tracks
    ├── bigbed/                    # UCSC browser tracks (TSS, peaks, summits)
    ├── tn5_motif/                 # optional Tn5 insertion-site PFM/logo outputs
    ├── count_matrices/
    │   ├── gene/                  # gene-TSS matrix (per-sample + aggregated)
    │   ├── trna/                  # tRNA matrix
    │   ├── pol3_pol3_t{1,2,3}/    # Pol III matrices (3 types)
    │   ├── repeat_{sine_alu,sine_mir,line_l1,line_l2,ltr,other_repeat_elements}/  # repeat matrices (6 types)
    │   ├── viral_{virus}/         # per-virus genome-bin matrix
    │   ├── rrna/                  # rRNA matrix
    │   └── final_matrices/        # consolidated symlinks + COUNT_MATRICES_INDEX.txt
    ├── deseq2/
    │   └── {group1}_vs_{group2}/  # one self-contained HTML report per contrast
    ├── multiqc/                   # multiqc_report.html + multiqc_data/
    └── multiqc_extra_data/
        └── custom/                # custom-content TSVs feeding MultiQC (kept out of multiqc's own scan dir)
```

## Reference Files

| File                              | Description                                                              |
| --------------------------------- | ------------------------------------------------------------------------ |
| `ref/ref.fa`                      | Concatenated host + additive + virus FASTA                               |
| `ref/ref.gtf`                     | Concatenated annotation GTF                                              |
| `ref/ref.fa.*.bt2`                | Bowtie2 index                                                            |
| `ref/ref.fa.{host,virus}.regions` | Per-organism contig lists (used to split BAMs downstream)                |
| `ref/ref.tss.{host,virus}.bed`    | TSS positions per organism, derived from `ref.gtf` `transcript` features |

## Per-Sample Output Files

### Trimmed Reads (`results/preprocess/`)

`{sample}.trimmed_R{1,2}.fastq.gz` plus fastp `{sample}.fastp.{html,json}` reports.

### Alignment Files (`results/align/`)

| File                                       | Stage                                                              | Notes                                                                       |
| ------------------------------------------ | ------------------------------------------------------------------ | --------------------------------------------------------------------------- |
| `{sample}.aligned.sorted.bam`              | raw sorted alignment                                               | temporary, deleted after filtering                                          |
| `{sample}.aligned.clean.bam`               | proper-pair, unmapped/secondary/supplementary filtered             | temporary                                                                   |
| `{sample}.aligned.fixmate.bam`             | mate info fixed, re-sorted                                         | temporary                                                                   |
| `{sample}.aligned.dedup.bam`               | duplicates marked/removed (`samtools markdup`)                     | temporary                                                                   |
| `{sample}.aligned.final.bam`               | MAPQ + excluded-contig filtered                                    | **canonical filtered alignment**, used everywhere downstream except Genrich |
| `{sample}.aligned.host.qname.bam`          | host contigs only, from the **unfiltered** sorted BAM, name-sorted | used by Genrich and the Genrich branch of Tn5 motif/count-matrix rules      |
| `{sample}.aligned.virus.{virus}.qname.bam` | same, per virus                                                    |                                                                             |

Each stage above also has a matching `.idxstats.txt`/`.flagstat.txt` side-output, aggregated by `alignmentqc.smk` into `results/alignmentqc/idxstats_summary.tsv`.

### Coverage Tracks

`results/bigwig/{sample}.{host,virus}.bw` — `deeptools bamCoverage` on the postprocess-split filtered BAM (50bp bins host, 1bp bins virus).

### Quality Control

**ataqv** (`results/alignmentqc/ataqv/`): per-case-sample `{sample}.{host,virus}.ataqv.json`, plus aggregate `final_report.host/` and `final_report.virus.{virus}/` interactive HTML directories (via `mkarv`), scoped to case samples only.

**FastQC** (`results/initqc/`): `{sample}.{raw,trimmed}_R{1,2}_fastqc.{html,zip}`.

## Aggregate Output Files

### Count Matrices {#count-matrices}

CHROMA2 generates up to **6 types of Tn5 insertion count matrices**, each capturing insertion events at a different genomic scale, all built exclusively from the Genrich-caller qname BAMs with **fractional (NH-weighted) counting** for multi-mapping reads (a single-mapped read contributes `+1.0`; an `N`-way multi-mapped read contributes `1/N` at each of its `N` locations — this avoids over-counting insertions in repetitive/multi-copy loci).

| #   | Matrix             | Path (per-type dir under `results/count_matrices/`)                                       | Binning                                                          | Use case                                           |
| --- | ------------------ | ----------------------------------------------------------------------------------------- | ---------------------------------------------------------------- | -------------------------------------------------- |
| 1   | Gene               | `gene/gene_count_matrix.tsv`                                                              | protein-coding TSS ± 250bp (excludes rRNA/tRNA/Mt_tRNA biotypes) | Promoter accessibility                             |
| 2   | tRNA               | `trna/trna_count_matrix.tsv`                                                              | gene-body midpoint ± 100bp, capped at 1000bp                     | tRNA locus accessibility                           |
| 3   | Pol III ×3         | `pol3_pol3_t{1,2,3}/pol3_pol3_t{1,2,3}_count_matrix.tsv`                                  | full feature span, capped at 1000bp                              | Pol III gene accessibility (5S/7SL RNA loci, etc.) |
| 4   | Repeat elements ×6 | `repeat_{sine_alu,sine_mir,line_l1,line_l2,ltr,other_repeat_elements}/*_count_matrix.tsv` | full feature span, capped at 1000bp                              | Heterochromatin/transposon accessibility           |
| 5   | Viral genome       | `viral_{virus}/viral_{virus}_count_matrix.tsv` (per virus)                                | fixed non-overlapping 200bp bins across the whole genome         | Viral genome accessibility/infection dynamics      |
| 6   | rRNA               | `rrna/rrna_count_matrix.tsv`                                                              | TSS ± 250bp                                                      | Ribosomal gene accessibility                       |

Every enabled matrix's aggregated TSV is also symlinked into `results/count_matrices/final_matrices/`, alongside a human-readable `COUNT_MATRICES_INDEX.txt` summarizing which matrices are available for the current host/virus selection and the counting parameters used.

**Shared TSV format:**

```text
feature_id    feature_name    [feature_type]    [strand]    [position]    sample1    sample2    ...
```

- `feature_id` — unique identifier (gene_id, repeat_id, bin_id, etc.)
- `feature_name` / optional type/strand/position metadata columns
- One column per sample, containing Tn5 insertion counts (integer if all contributing reads were uniquely mapped, fractional otherwise)

**Per-sample vs. aggregated:** `{type}/{sample}.{type}_counts.tsv` holds one sample's counts (useful for per-sample QC); `{type}/{type}_count_matrix.tsv` is the column-wise concatenation across all samples (no normalization applied — normalize downstream as needed, e.g. via DESeq2's size factors or a simple library-size division).

**Example: loading in R**

```r
library(data.table)
gene_counts <- fread("results/count_matrices/final_matrices/gene_count_matrix.tsv")
counts_matrix <- as.matrix(gene_counts[, -(1:5)])  # drop metadata columns
normalized <- sweep(counts_matrix, 2, colSums(counts_matrix), "/")
```

### Alignment Quality Summary

`results/alignmentqc/idxstats_summary.tsv` — per-sample read counts at every filter stage (raw sorted → clean → dedup → final) plus host/virus split, feeding the MultiQC `alignment_stats` and `host_virus_ratio` custom-content panels.

### Peak Calls

`results/peaks/{sample}.{host,virus}.macs2.narrowPeak.gz` (+ `.summits.bed.gz`) and `results/peaks/{sample}.{host,virus}.genrich.narrowPeak.gz`. UCSC-browser-ready `.bb` equivalents live under `results/bigbed/`.

### Tn5 Motif Outputs (Optional)

Only generated if `tn5_motif.generate_logo: true`. Per sample × caller × organism: `{sample}.{group}.{caller}.tn5_sites.1bp.bed.gz`, `...pfm.tsv.gz` (position-frequency matrix), `...logo.{png|...}` (sequence logo image) under `results/tn5_motif/`.

## DESeq2 Reports {#deseq2-reports}

Only generated if `deseq2.enabled: true` and `contrasts.tsv` is present. One self-contained HTML report per contrast at `results/deseq2/{group1}_vs_{group2}/{group1}_vs_{group2}.deseq2_report.html`, analyzing every available count-matrix type (gene matrix as the primary analysis with a volcano plot; other enabled matrix types analyzed alongside, with any `deseq2.skip_features` — tRNA by default — flagged and skipped with a stated reason directly in the report).

## MultiQC Report

`results/multiqc/multiqc_report.html` + `multiqc_data/` — aggregates FastQC, fastp, alignment/filtering stats, host-vs-virus read distribution, ataqv metrics, fragment-size distributions, FRiP scores, peak-size distributions, genome coverage breadth, and Tn5 count-matrix summaries. The custom-content TSVs feeding the non-standard panels live in `results/multiqc_extra_data/custom/` (deliberately outside MultiQC's own output-scan directory).

## S3 Transfer Sentinel (Optional)

`$WORKDIR/.s3_transfer.done` — only meaningfully populated if `push_to_s3: true` and `s3_sample_set_name` is set; always touched (as a no-op) otherwise so it doesn't block the pipeline's final target. See [S3 Configuration](s3_configuration.md).

## Output Naming Conventions

### Sample Names

All per-sample outputs are prefixed with the exact `sampleName` from `samples.tsv`.

### Organism/Genomic Region Names

Host-space outputs use `.host.`; virus-space outputs use `.virus.{accession}.` (or `.{accession}.` for count-matrix directories), where `{accession}` is the exact viral accession supplied via `-v/--viruses` at `init`.

## How to Use These Outputs

### For Differential Accessibility Analysis

1. Enable `deseq2` and provide `contrasts.tsv` (see [Inputs](inputs.md#contrasts-manifest)) for automated per-contrast reports, or
2. Load `results/count_matrices/final_matrices/*.tsv` directly into your own DESeq2/edgeR/limma workflow.

### For Quality Assessment

Start with `results/multiqc/multiqc_report.html`, then drill into `results/alignmentqc/ataqv/final_report.*/` for per-sample accessibility QC detail.

### For Visualization

Load `results/bigwig/*.bw` and `results/bigbed/*.bb` into IGV or the UCSC Genome Browser.

### For Downstream Integration

`results/count_matrices/final_matrices/` is the single directory to point any downstream R/Python analysis at — it contains every enabled matrix type with a stable naming scheme.

## Frequently Asked Questions

**Q: Why are peak calls generated by two different tools (MACS2 and Genrich)?**
A: They use different underlying models and BAM inputs (MACS2 on the fully filtered BAM, Genrich on the unfiltered qname-sorted BAM with its own duplicate handling) — comparing both gives a sanity check on peak-calling robustness, and downstream FRiP/count-matrix rules can be pointed at whichever caller's output fits your analysis.

**Q: Why do count matrices only use the Genrich BAM, not the MACS2 BAM?**
A: Consistency — using one canonical BAM source for all six matrix types avoids introducing a peak-caller-dependent bias into the count-matrix comparisons across feature classes.

**Q: Some Pol III/repeat-element/tRNA matrices are missing for my host genome — why?**
A: Those GTFs are configured per host in `config.yaml`; if a GTF isn't defined for your selected host, that matrix type is unavailable and won't appear in `final_matrices/`. Check `COUNT_MATRICES_INDEX.txt` for exactly which matrices were built for your run.

## Troubleshooting

**Missing outputs:** check `$WORKDIR/pipeline.failed` (if present) for which rule failed, and the corresponding SLURM log under `$WORKDIR/logs/`.

**Missing count matrices:** check `results/count_matrices/final_matrices/COUNT_MATRICES_INDEX.txt` to confirm the matrix type was enabled and a GTF was resolved for your host/virus selection.

**Zero counts across all samples in a matrix:** likely a coordinate mismatch between the bin file and BAM contigs, or genuinely low Tn5 signal at those features (e.g. rRNA/repeat loci can be low-accessibility by biology, not by bug) — cross-check against `results/alignmentqc/idxstats_summary.tsv` for overall alignment health first.
