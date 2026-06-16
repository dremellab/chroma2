# GTF-Based Count Matrix Expansion — Feature Plan

**Status**: v1.4 | **Updated**: 2026-06-04  
**Focus**: hg38 and hs1 only (mm39 deferred for future)  
**No backward compatibility** — fresh implementation  
**Ready for iterative updates**

---

## Feature Implementation Status by Genome

| Feature Type | hg38 | hs1 (T2T) | Notes |
|---|---|---|---|
| **tRNA** | ✅ Ready | ✅ Ready (enriched) | hg38: `hg38.tRNA.hg38chroms.gtf` (basic); hs1: enriched GTF with tRNAscan scores, pseudo-marking, nomenclature (230 novel discoveries included) |
| | | | **Status**: Both genomes have usable tRNA GTFs; hs1 has additional quality metadata |
| | | | **Gaps**: (1) hs1 enriched file needs copying to main folder; (2) hg38 lacks scoring—keep as-is OR run tRNAscan (user choice) |
| | | | **Low-confidence tRNA**: hs1 includes ~230 score < 50 genes (flagged `.low_confidence`); kept by default; configurable filter in Phase 4 |
| | | | **Pseudo-tRNA**: hs1 marks non-functional copies with `pseudo_` prefix; kept by default; configurable filter in Phase 4 |
| **POL3 Class 1** | ✅ Ready | ✅ Ready | T1_Genes (RNA5S): hg38 and hs1 both available; hg38 file lifted to hs1 via liftOver |
| | | | **Status**: Both genomes ready; needs renaming from T1_Genes.* to POL3_class1.* naming |
| | | | **Content**: RNA5S ribosomal RNA genes |
| **POL3 Class 2** | ✅ Ready | ✅ Ready | T2_noAlu (POL3-dependent, non-Alu): hg38 and hs1 both available; hg38 file lifted to hs1 |
| | | | **Status**: Both genomes ready; needs renaming from T2_noAlu.* to POL3_class2.* naming |
| | | | **Content**: POL3-dependent genes excluding Alu elements |
| **POL3 Class 3** | ✅ Ready | ✅ Ready | T3_Genes (RMRP, other POL3): hg38 and hs1 both available; hg38 file lifted to hs1 |
| | | | **Status**: Both genomes ready; needs renaming from T3_Genes.* to POL3_class3.* naming |
| | | | **Content**: Other POL3-dependent genes (RMRP, etc.) |
| **ALU Elements** | ❌ Not ready | ⚠️ Partial | hg38.repeats.hg38chroms.gtf mixed (needs ALU extraction script); hs1 has chm13v2.0_Alu.gtf (836M, needs chromosome filtering) |
| | | | **Status**: hs1 type-specific GTF exists; hg38 needs extraction from mixed repeats file |
| | | | **Issues**: hg38 needs filtering script to extract ALU; hs1 needs chromosome filtering |
| **LINE Elements** | ❌ Not ready | ⚠️ Partial | hg38.repeats.hg38chroms.gtf mixed (needs L1 extraction); hs1 has chm13v2.0_L1.gtf (664M, needs chromosome filtering) |
| | | | **Status**: hs1 type-specific GTF exists; hg38 needs extraction |
| | | | **Issues**: hg38 needs filtering script; hs1 needs chromosome filtering |
| **SINE Elements** | ❌ Not ready | ⚠️ Partial | hg38.repeats.hg38chroms.gtf mixed (needs SINE extraction); hs1 has chm13v2.0_SINE.gtf (1.2G, needs chromosome filtering) |
| | | | **Status**: hs1 type-specific GTF exists; hg38 needs extraction |
| | | | **Issues**: hg38 needs filtering script; SINE/ALU overlap (ALU is SINE subtype) |
| **LTR Elements** | ❌ Not ready | ⚠️ Partial | hg38.repeats.hg38chroms.gtf mixed (needs LTR extraction); hs1 has chm13v2.0_HERV.gtf (7.7M, needs chromosome filtering) |
| | | | **Status**: hs1 HERV GTF exists; hg38 needs extraction |
| | | | **Issues**: hg38 needs filtering script; HERVs vs all LTRs distinction |

---

## Data Files Quick Reference

### hg38

| Feature | Current File | Status | Location | Action Required |
|---|---|---|---|---|
| tRNA | hg38.tRNA.hg38chroms.gtf | ✅ Ready | `/project/dremel_lab/workflows/reference_data/fasta_gtf/` | None — already prepared |
| repeats (mixed) | hg38.repeats.hg38chroms.gtf | ❌ Needs extraction | `/project/dremel_lab/workflows/reference_data/fasta_gtf/` | Extract ALU, LINE, SINE, LTR type-specific GTFs (create filtering script or obtain files) |
| POL3 Class 1 | T1_Genes.hg38.gtf | ✅ Ready | `/project/dremel_lab/workflows/reference_data/fasta_gtf/` | Rename to `hg38.POL3_class1.hg38chroms.gtf` |
| POL3 Class 2 | T2_noAlu.hg38.gtf | ✅ Ready | `/project/dremel_lab/workflows/reference_data/fasta_gtf/` | Rename to `hg38.POL3_class2.hg38chroms.gtf` |
| POL3 Class 3 | T3_Genes.hg38.gtf | ✅ Ready | `/project/dremel_lab/workflows/reference_data/fasta_gtf/` | Rename to `hg38.POL3_class3.hg38chroms.gtf` |
| chrR (rRNA genes) | hg38.chrR.gtf | ✅ Ready | `/project/dremel_lab/workflows/reference_data/fasta_gtf/` | Not used for tRNA bins; for reference only |

**Overall hg38 Status**: POL3 classes ✅ ready (need renaming); tRNA ✅ ready; repeats ❌ need type-specific filtering/extraction (unlike hs1 which has pre-filtered type-specific GTFs)

---

### hs1 (T2T-CHM13v2.0)

| Feature | Current File | Status | Location | Action Required |
|---|---|---|---|---|
| tRNA | all_tRNA_enriched.gtf | ✅ Ready | `/project/dremel_lab/workflows/reference_data/fasta_gtf/temp/hs1/Human_hs1-rDNA_genome_v1.0/trnascan_out/` | Copy to main folder as `hs1.tRNA.hs1chroms.gtf` |
| POL3 Class 1 | T1_Genes.hs1.gtf | ✅ Ready | `/project/dremel_lab/workflows/reference_data/fasta_gtf/` | Rename to `hs1.POL3_class1.hs1chroms.gtf` |
| POL3 Class 2 | T2_noAlu.hs1.gtf | ✅ Ready | `/project/dremel_lab/workflows/reference_data/fasta_gtf/` | Rename to `hs1.POL3_class2.hs1chroms.gtf` |
| POL3 Class 3 | T3_Genes.hs1.gtf | ✅ Ready | `/project/dremel_lab/workflows/reference_data/fasta_gtf/` | Rename to `hs1.POL3_class3.hs1chroms.gtf` |
| ALU | chm13v2.0_Alu.gtf (836M) | ⚠️ Needs chrom filter | `/project/dremel_lab/workflows/reference_data/fasta_gtf/temp/hs1/` | Filter to hs1 standard chromosomes (chr1-22, X, Y, M) |
| LINE | chm13v2.0_L1.gtf (664M) | ⚠️ Needs chrom filter | `/project/dremel_lab/workflows/reference_data/fasta_gtf/temp/hs1/` | Filter to hs1 standard chromosomes |
| SINE | chm13v2.0_SINE.gtf (1.2G) | ⚠️ Needs chrom filter | `/project/dremel_lab/workflows/reference_data/fasta_gtf/temp/hs1/` | Filter to hs1 standard chromosomes |
| LTR/HERV | chm13v2.0_HERV.gtf (7.7M) | ⚠️ Needs chrom filter | `/project/dremel_lab/workflows/reference_data/fasta_gtf/temp/hs1/` | Filter to hs1 standard chromosomes |

**Overall hs1 Status**: POL3 classes ✅ ready (need renaming); tRNA ✅ ready (need copying); repeats ⚠️ type-specific GTFs exist (created by Python script), need chromosome filtering

---

## Final Count Matrices Generated by Pipeline

The following **count matrices** will be generated from Genrich BAM input after the feature expansion is complete. Each matrix contains Tn5 insertion site counts per sample, aggregated across all samples for downstream analysis (DESeq2, etc.).

### Host-Based Count Matrices (9 total)

1. **Gene count matrix** — Tn5 insertion counts within full gene bodies (±250bp flanking region)
2. **tRNA count matrix** — Tn5 insertion counts within tRNA gene bodies (±100bp flanking region)
3. **POL3 Class 1 count matrix** — Tn5 insertion counts within RNA5S ribosomal RNA genes
4. **POL3 Class 2 count matrix** — Tn5 insertion counts within POL3-dependent genes (non-Alu)
5. **POL3 Class 3 count matrix** — Tn5 insertion counts within other POL3-dependent genes
6. **ALU element count matrix** — Tn5 insertion counts within ALU repeat elements
7. **LINE element count matrix** — Tn5 insertion counts within LINE/L1 repeat elements
8. **SINE element count matrix** — Tn5 insertion counts within SINE repeat elements
9. **LTR/HERV count matrix** — Tn5 insertion counts within LTR/HERV retrotransposon elements

### Viral-Based Count Matrices

10. **Viral genome bin count matrix** (per virus) — Tn5 insertion counts in 200bp non-overlapping intervals across each viral genome

---

**Total output**: 9 host feature matrices + 1 per virus = **(9 + number of viruses)** count matrices
