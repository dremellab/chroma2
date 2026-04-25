## dev version

- feat(deseq2): add DESeq2 contrast reporting module
- add `deseq2_contrast_report` rule driven by a user-supplied contrasts TSV; generates per-comparison HTML reports and differential accessibility TSVs for host gene bins and per-virus bins
- add full config schema (`deseq2` block) with validated parameters and init-time contrast validation
- add `deseq2_report` container and Rivanna resource settings
- feat(trna): refactor tRNA analysis into a dedicated count matrix pipeline (closes #11)
- separate tRNA genes from protein-coding genes into their own bin reference (`ref.tn5.host_trna_gene_bins.{flank}bp.bed`) using gene body center ± configurable flank (default 100 bp) instead of TSS
- add per-sample tRNA count rules for genrich and macs2 callers and aggregate tRNA count matrices
- remove tRNA genes from the protein-coding host gene bin pipeline (no regression)
- add `trna_gene_flank_size` config key (default 100 bp) under `tn5_motif`
- GTF filtering for tRNA gene_type is case-insensitive (tRNA/trna/TRNA all match)
- add Rivanna SLURM resource settings for all 5 new tRNA rules
- chore(ci): remove `style-files` pre-commit hook (Rscript unavailable in this environment)

- feat(bigbed): add UCSC BigBed outputs for browser-facing BED tracks
- convert host/virus TSS BEDs, MACS2 summit BEDs, and MACS2/Genrich narrowPeak outputs to `.bb`
- add `bedToBigBed` container support and Rivanna resource settings for the new conversion rules
- feat(tn5-motif): add Tn5 motif extraction outputs to the workflow
- generate caller-specific host and viral Tn5 BED/FASTA/PFM/logo outputs from the Genrich and MACS2 input BAMs
- add viral 100 bp bin counts, host TSS-centered gene bins, per-sample Tn5 count TSVs, and aggregate count matrices across samples
- add Tn5 motif config defaults, helper scripts, `pysam`/`py311` container support, and Rivanna resources for motif extraction and counting jobs
