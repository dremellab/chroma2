## dev version

- feat(bigbed): add UCSC BigBed outputs for browser-facing BED tracks
- convert host/virus TSS BEDs, MACS2 summit BEDs, and MACS2/Genrich narrowPeak outputs to `.bb`
- add `bedToBigBed` container support and Rivanna resource settings for the new conversion rules
- feat(tn5-motif): add Tn5 motif extraction outputs to the workflow
- generate caller-specific host and viral Tn5 BED/FASTA/PFM/logo outputs from the Genrich and MACS2 input BAMs
- add viral 100 bp bin counts, host TSS-centered gene bins, per-sample Tn5 count TSVs, and aggregate count matrices across samples
- add Tn5 motif config defaults, helper scripts, `pysam`/`py311` container support, and Rivanna resources for motif extraction and counting jobs
