# GTF-Based Count Matrix Feature Expansion

## Overview
Expand the Tn5 count matrix pipeline from single-type (tRNA genes) to multi-type support, including POL3 elements and repeat-related annotations. Each type will follow the established tRNA pattern: bin creation → per-sample counting → matrix aggregation.

---

## Currently Implemented

### tRNA Gene Bins
- **Status**: ✅ Implemented (commits 88859f4, 7c1f635)
- **Pattern**: Gene body center ± configurable flank (default 100 bp)
- **Files**:
  - Input: `{host}.tRNAs.{host}chroms.gtf`
  - Output bin: `ref.tn5.host_trna_gene_bins.{flank}bp.bed`
  - Per-sample counts: `{sample}.host.{caller}.tn5_trna_gene_counts.{flank}bp.tsv`
  - Matrix: `host.{caller}.tn5_trna_gene_count_matrix.{flank}bp.tsv`
- **Callers**: Genrich, MACS2
- **Config**:
  - `trnas_gtf`: Host-specific GTF mapping
  - `tn5_motif.trna_gene_flank_size`: Bin width (default: 100 bp)
  - `include_trnas_gtf_in_ref`: Include in reference (default: true)

---

## Proposed Features

### 1. POL3 BED-Based Bins (3 variants)

#### 1.1 POL3 Core Elements
- **Type**: Direct BED input (no GTF parsing needed)
- **Description**: RNA Polymerase III-transcribed genes without flanking regions
- **Files**:
  - Input: `{host}.POL3.core.bed` 
  - Output bin: `ref.tn5.host_pol3_core_bins.bed` (direct copy/normalize)
  - Per-sample counts: `{sample}.host.{caller}.tn5_pol3_core_counts.tsv`
  - Matrix: `host.{caller}.tn5_pol3_core_count_matrix.tsv`
- **Config**:
  - `pol3_core_bed`: Host-specific BED mapping
  - `include_pol3_core_in_analysis`: Toggle (default: true if file exists)

#### 1.2 POL3 with 500bp Flank
- **Type**: BED input with flank expansion
- **Description**: POL3 elements ± 500 bp flanking regions
- **Files**:
  - Input: `{host}.POL3.core.bed`
  - Output bin: `ref.tn5.host_pol3_500bp_bins.bed`
  - Per-sample counts: `{sample}.host.{caller}.tn5_pol3_500bp_counts.tsv`
  - Matrix: `host.{caller}.tn5_pol3_500bp_count_matrix.tsv`
- **Config**:
  - `pol3_core_bed`: Host-specific BED mapping
  - `tn5_motif.pol3_flank_size`: Flank width (default: 500 bp)
  - `include_pol3_500bp_in_analysis`: Toggle

#### 1.3 POL3 with 1000bp Flank
- **Type**: BED input with extended flank
- **Description**: POL3 elements ± 1000 bp flanking regions
- **Files**:
  - Input: `{host}.POL3.core.bed`
  - Output bin: `ref.tn5.host_pol3_1000bp_bins.bed`
  - Per-sample counts: `{sample}.host.{caller}.tn5_pol3_1000bp_counts.tsv`
  - Matrix: `host.{caller}.tn5_pol3_1000bp_count_matrix.tsv`
- **Config**:
  - `include_pol3_1000bp_in_analysis`: Toggle

---

### 2. Repeat-Related GTF Bins (4 variants)

#### 2.1 ALU Elements
- **Type**: GTF with gene_type filtering
- **Description**: Alu repeat elements (young, active)
- **Pattern**: Repeat body center ± configurable flank (default: 100 bp)
- **Files**:
  - Input: `{host}.repeats.gtf` (filtered for ALU)
  - Output bin: `ref.tn5.host_alu_gene_bins.{flank}bp.bed`
  - Per-sample counts: `{sample}.host.{caller}.tn5_alu_gene_counts.{flank}bp.tsv`
  - Matrix: `host.{caller}.tn5_alu_gene_count_matrix.{flank}bp.tsv`
- **Config**:
  - `repeats_gtf`: Host-specific GTF mapping
  - `tn5_motif.alu_gene_flank_size`: Bin width (default: 100 bp)
  - `include_alu_in_analysis`: Toggle

#### 2.2 LINE Elements (L1)
- **Type**: GTF with gene_type filtering
- **Description**: LINE-1 long interspersed elements
- **Pattern**: Element center ± configurable flank (default: 100 bp)
- **Files**:
  - Input: `{host}.repeats.gtf` (filtered for LINE)
  - Output bin: `ref.tn5.host_line_gene_bins.{flank}bp.bed`
  - Per-sample counts: `{sample}.host.{caller}.tn5_line_gene_counts.{flank}bp.tsv`
  - Matrix: `host.{caller}.tn5_line_gene_count_matrix.{flank}bp.tsv`
- **Config**:
  - `tn5_motif.line_gene_flank_size`: Bin width (default: 100 bp)
  - `include_line_in_analysis`: Toggle

#### 2.3 SINE Elements
- **Type**: GTF with gene_type filtering
- **Description**: Short interspersed nuclear elements (e.g., Alu in primates, B1 in mice)
- **Pattern**: Element center ± configurable flank (default: 50 bp)
- **Files**:
  - Input: `{host}.repeats.gtf` (filtered for SINE)
  - Output bin: `ref.tn5.host_sine_gene_bins.{flank}bp.bed`
  - Per-sample counts: `{sample}.host.{caller}.tn5_sine_gene_counts.{flank}bp.tsv`
  - Matrix: `host.{caller}.tn5_sine_gene_count_matrix.{flank}bp.tsv`
- **Config**:
  - `tn5_motif.sine_gene_flank_size`: Bin width (default: 50 bp)
  - `include_sine_in_analysis`: Toggle

#### 2.4 LTR Elements (Retrotransposons)
- **Type**: GTF with gene_type filtering
- **Description**: Long terminal repeat retrotransposons
- **Pattern**: LTR element center ± configurable flank (default: 100 bp)
- **Files**:
  - Input: `{host}.repeats.gtf` (filtered for LTR)
  - Output bin: `ref.tn5.host_ltr_gene_bins.{flank}bp.bed`
  - Per-sample counts: `{sample}.host.{caller}.tn5_ltr_gene_counts.{flank}bp.tsv`
  - Matrix: `host.{caller}.tn5_ltr_gene_count_matrix.{flank}bp.tsv`
- **Config**:
  - `tn5_motif.ltr_gene_flank_size`: Bin width (default: 100 bp)
  - `include_ltr_in_analysis`: Toggle

---

## Implementation Roadmap

### Phase 1: Consolidate Bin Building (Low Priority, High Value)
Create a generic bin builder that handles both GTF-based and BED-based inputs:
- **New script**: `build_generic_bins.py`
- **Input modes**:
  - GTF with gene_type filter → center-based bins
  - BED direct copy → identity mapping
  - BED with flank expansion → widened bins
- **Supports**: metadata preservation (gene_id, gene_name, gene_type, strand, center)

### Phase 2: Consolidate Counting (Medium Priority, High Reusability)
Generalize counting to all bin types (already flexible via CLI args):
- Current `count_tn5_sites_in_bins.py` works generically
- Create Snakemake wrapper rules for each type
- Handle both Genrich and MACS2 callers

### Phase 3: Consolidate Matrix Building (Medium Priority)
Generalize count matrix aggregation:
- Current `build_tn5_count_matrix.py` works generically
- Create Snakemake wrapper rules for each type

### Phase 4: Configuration & Integration (Medium Priority)
- Add config keys for each type (include toggles, flank sizes)
- Integrate into DESeq2 analysis (skip_features list)
- Update CHANGELOG and documentation

### Phase 5: Refactor tn5motif.smk (Optional, High Quality)
- Replace repetitive rule blocks with a rule factory pattern
- Parameterize by bin type, caller, flank size
- Reduces ~600+ lines to ~100-150 lines of template code

---

## Configuration Schema Addition

```yaml
# Add to config.yaml under tn5_motif block:
tn5_motif:
  # existing keys...
  
  # POL3 Configuration
  pol3_core_bed:
    hg38: "hg38.POL3.core.bed"
    mm39: "mm39.POL3.core.bed"
  pol3_flank_size: 500  # for 500bp variant
  include_pol3_core_in_analysis: true
  include_pol3_500bp_in_analysis: true
  include_pol3_1000bp_in_analysis: false  # optional
  
  # Repeat Configuration
  repeats_gtf:
    hg38: "hg38.repeats.gtf"
    mm39: "mm39.repeats.gtf"
  alu_gene_flank_size: 100
  line_gene_flank_size: 100
  sine_gene_flank_size: 50
  ltr_gene_flank_size: 100
  
  # Repeat Analysis Toggles
  include_alu_in_analysis: true
  include_line_in_analysis: false  # optional
  include_sine_in_analysis: false  # optional
  include_ltr_in_analysis: false   # optional
```

---

## Snakemake Integration Points

### Bin Creation Rules
- `build_host_tn5_trna_gene_bins` → generalize
- `build_host_tn5_pol3_core_bins` (new)
- `build_host_tn5_pol3_500bp_bins` (new)
- `build_host_tn5_pol3_1000bp_bins` (new)
- `build_host_tn5_alu_gene_bins` (new)
- `build_host_tn5_line_gene_bins` (new)
- `build_host_tn5_sine_gene_bins` (new)
- `build_host_tn5_ltr_gene_bins` (new)

### Counting Rules (per caller per type)
- `count_tn5_host_{caller}_{type}_bins` (8 rule definitions × 2 callers = 16 rules)
  - Types: trna_gene, pol3_core, pol3_500bp, pol3_1000bp, alu_gene, line_gene, sine_gene, ltr_gene

### Matrix Building Rules (per caller per type)
- `build_tn5_host_{caller}_{type}_count_matrix` (8 rule definitions × 2 callers = 16 rules)

### DESeq2 Integration
- Update `deseq2.smk` to conditionally include tRNA, ALU, LINE, SINE, LTR matrices
- Skip types based on `skip_features` list
- Generate volcano plots for enabled types

### Outputs (rule all)
- Add conditional outputs for enabled feature types
- Default: tRNA only (backward compatibility)
- User can enable additional types via config

---

## Testing Strategy

1. **Unit tests**: Expand existing test suite
   - `test_build_host_trna_gene_bins.py` → `test_build_generic_bins.py`
   - Add test cases for BED input, flank expansion, repeat filtering

2. **Integration tests**: 
   - Small test BAM with known Tn5 sites
   - Verify counts match manual BED overlap
   - Validate matrix schema across all types

3. **Data validation**:
   - GTF parsing error handling
   - BED coordinate edge cases (negative starts, beyond chromosome bounds)
   - Empty feature sets (graceful skip)

---

## Backward Compatibility

✅ **Fully backward compatible**
- tRNA remains the only enabled type by default
- Existing configs work unchanged
- New types opt-in via `include_{type}_in_analysis: true`
- DESeq2 analysis unchanged if repeat types not enabled
- CHANGELOG documents breaking changes (none)

---

## Open Questions / Clarifications Needed

1. **POL3 inputs**: Confirm 3 BED files (core, 500bp-flanked, 1000bp-flanked) or prefer single BED file with flank parameter?
2. **Repeat GTF**: Single combined `repeats.gtf` with gene_type filtering, or separate files per repeat type?
3. **Flank sizes**: Are the defaults (100 bp for tRNA/ALU/LINE/LTR, 50 bp for SINE) reasonable?
4. **Priority order**: Which types should be implemented first? (suggest: POL3 → ALU → LINE → SINE → LTR)
5. **Data availability**: Do reference genomes (hg38, mm39) already have these GTF/BED files in `REFS_DIR`?
