# S3 Output Hierarchy: HAROLD vs Chroma 2

## HAROLD S3 Hierarchy

HAROLD pipelines use the following S3 directory structure for output archival:

```
s3://bucket/s3_prefix/sample_set_name/
├── config/
│   ├── samples.tsv                     # Sample metadata
│   ├── config.yaml                     # Pipeline config
│   └── rivanna/config.yaml             # HPC-specific config
├── qc/
│   ├── alignment_summary.tsv           # Aggregated alignment stats
│   ├── multiqc_report.html             # MultiQC HTML report
│   └── multiqc_data/                   # MultiQC data directory (preserves structure)
├── bigwigs/                            # BigWig files (flattened)
│   └── *.bw, *.bb                      # Flattened file names
├── bams/                               # BAM files (flattened, excludes intermediates)
│   └── *.bam, *.bai
├── SJ/                                 # Splice junction files (flattened)
│   └── SJ.out.tab
└── counts/                             # Count matrices (structure preserved)
    └── normalized_counts/
        └── (all files including quarto output)
```

### Key HAROLD Rules:

1. **Config files** (path rule): Direct mapping to config/
2. **Alignment summary** (path rule): Direct mapping to qc/
3. **BigWigs/BigBeds** (suffix rule): Flattened into single directory, preserves filenames
4. **BAM/BAI files** (suffix rule): Flattened into single directory, excludes "Aligned.out.bam" intermediates
5. **Splice junctions** (suffix rule): Flattened into single directory, excludes "\_STARpass1" intermediates
6. **MultiQC files** (dir rule): Directory structure preserved (qc/multiqc_data/)
7. **Counts matrices** (dir rule): Full directory structure preserved

### Storage Classes:

- **Default**: GLACIER_IR (infrequent access, retrievable in hours)
- **Large files** (BAM/BAI): GLACIER (deep archive, retrievable in hours to days)

---

## Proposed Chroma 2 S3 Hierarchy

Adapting HAROLD's pattern for Chroma 2's ATAC-seq outputs:

```
s3://bucket/s3_prefix/sample_set_name/
├── config/
│   ├── samples.tsv                     # Sample metadata
│   ├── config.yaml                     # Pipeline config
│   └── contrasts.tsv                   # Contrast definitions (if present)
├── qc/
│   ├── alignment_summary.tsv           # Aggregated alignment stats (flagstats + fastqc)
│   ├── multiqc_report.html             # MultiQC HTML (if enabled)
│   ├── multiqc_data/                   # MultiQC data directory (if enabled)
│   └── ataqv/                          # ATAC-seq QC reports (flattened)
│       ├── final_report.host
│       ├── final_report.virus.{virus}  # Per-virus if applicable
│       └── {sample}.*.json
├── bams/                               # BAM files - sample-level (flattened)
│   └── {sample}.aligned.final.bam*
├── bigwigs/                            # BigWig files (flattened)
│   └── {sample}.*.bw
├── bigbeds/                            # BigBed files (flattened)
│   └── {sample}.*.bb
├── peaks/                              # Peak files (flattened)
│   └── {sample}.*.narrowPeak.gz
│   └── {sample}.*.summits.bed.gz
├── tn5_counts/                         # Tn5 motif counting outputs (structure preserved)
│   ├── *_tn5_gene_count_matrix.*.tsv
│   ├── *_tn5_bin_count_matrix.*.tsv
│   ├── *_tn5_trna_gene_count_matrix.*.tsv
│   └── (all supporting files)
└── deseq2/                             # DESeq2 contrast reports (structure preserved)
    ├── {comparison}/
    │   ├── {comparison}.*.deseq2.tsv
    │   ├── {comparison}.*.volcano.png
    │   └── {comparison}.report.html
    └── (more contrasts as subdirs)
```

### Chroma 2 Transfer Rules:

| Rule Type | Pattern                          | Source                                     | Destination                | Storage Class | Notes                           |
| --------- | -------------------------------- | ------------------------------------------ | -------------------------- | ------------- | ------------------------------- |
| path      | samples.tsv                      | `samples.tsv`                              | `config/samples.tsv`       | GLACIER_IR    | Sample metadata                 |
| path      | config.yaml                      | `config.yaml`                              | `config/config.yaml`       | GLACIER_IR    | Pipeline config                 |
| path      | contrasts.tsv                    | `contrasts.tsv`                            | `config/contrasts.tsv`     | GLACIER_IR    | Optional, if exists             |
| path      | alignmentqc/idxstats_summary.tsv | `results/alignmentqc/idxstats_summary.tsv` | `qc/alignment_summary.tsv` | GLACIER_IR    | Aggregated stats                |
| path      | multiqc_report.html              | `results/multiqc_report.html`              | `qc/multiqc_report.html`   | GLACIER_IR    | Optional, if exists             |
| dir       | multiqc_data                     | `results/multiqc_data/`                    | `qc/multiqc_data/`         | GLACIER_IR    | Optional, if exists             |
| suffix    | .aligned.final.bam, .bai         | `results/*/align/`                         | `bams/`                    | GLACIER       | Flattened, final alignments     |
| suffix    | .bw                              | `results/*/bigwig/`                        | `bigwigs/`                 | GLACIER_IR    | Flattened, all samples          |
| suffix    | .bb                              | `results/*/bigwig/`                        | `bigbeds/`                 | GLACIER_IR    | Flattened, all samples          |
| suffix    | .narrowPeak.gz, .summits.bed.gz  | `results/*/peaks/`                         | `peaks/`                   | GLACIER_IR    | Flattened, all peak types       |
| suffix    | .json (from ataqv)               | `results/*/alignmentqc/ataqv/`             | `qc/ataqv/`                | GLACIER_IR    | Flattened JSON reports          |
| path      | ataqv final_report.\*            | `results/alignmentqc/ataqv/final_report.*` | `qc/ataqv/final_report.*`  | GLACIER_IR    | Aggregate reports               |
| dir       | tn5_motif                        | `results/tn5_motif/`                       | `tn5_counts/`              | GLACIER_IR    | Preserves matrix file structure |
| dir       | deseq2                           | `results/deseq2/`                          | `deseq2/`                  | GLACIER_IR    | Preserves comparison structure  |

### Key Design Decisions:

1. **Sample-level outputs are flattened** (BAM, BigWig, peaks) for easy discovery
2. **Aggregate outputs preserve structure** (DESeq2 contrasts, Tn5 matrices) to maintain organization
3. **QC reports consolidated** in `qc/` directory (including ataqv reports) similar to HAROLD
4. **Config files consolidated** in `config/` directory for easy access
5. **Storage classes**:
   - **GLACIER**: All BAM/BAI files (large files, infrequent access, deep archive)
   - **GLACIER_IR**: All other files (config, QC reports, BigWigs, peaks, count matrices, DESeq2 reports)
6. **Excluded outputs**: FastQC reports, pooled input controls, and reference files are not uploaded to S3

---

## Implementation Notes

The script (`s3_transfer_chroma2.py`) will:

- Walk the local results directory recursively
- Match files against these rules in order
- Construct S3 paths following the hierarchy above
- Assign storage classes based on file type:
  - **GLACIER**: BAM and BAI files (large, infrequent access)
  - **GLACIER_IR**: All other files (config, reports, matrices, counts)
- Handle errors gracefully (retain local files if transfer fails)
- Verify checksums before deleting local copies
