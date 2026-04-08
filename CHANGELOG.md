## dev version

- feat(bigbed): add UCSC BigBed outputs for browser-facing BED tracks
- convert host/virus TSS BEDs, MACS2 summit BEDs, and MACS2/Genrich narrowPeak outputs to `.bb`
- add `bedToBigBed` container support and Rivanna resource settings for the new conversion rules
- feat(tn5-motif): add Tn5 motif extraction outputs to the workflow
- add scenario-aware Tn5 BED/FASTA/PFM/logo generation and a dedicated `extract_tn5_motifs.py` helper
- add Tn5 motif config defaults, `pysam` container support, and Rivanna resources for motif extraction jobs
