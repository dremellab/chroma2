# Count Matrices Generation Plan

## Overview

Generate 6 types of count matrices capturing Tn5 insertions at different genomic scales:

---

## 1. Gene Count Matrix

**Description:** Tn5 insertions around TSS of protein-coding genes (±250 bp)

**Config:**

```yaml
count_matrices:
  gene_counts:
    enabled: true
    flank_size: 250 # ±250 bp around TSS
    exclude_features: [rRNA, tRNA, Mt_tRNA] # Exclude specific biotypes
```

**Workflow:**

```
GTF (protein_coding)
  → build_tn5_tss_bins.py (host_gene_flank_size)
  → per-sample count_tn5_sites_in_bins.py
  → aggregate build_tn5_count_matrix.py
  → host.gene.tn5_count_matrix.tsv
```

**Rules to add:**

- `build_host_gene_tss_bins` - calls `build_tn5_tss_bins.py`
- `count_tn5_host_gene_bins` (per sample)
- `aggregate_tn5_host_gene_count_matrix`

---

## 2. tRNA Count Matrix

**Description:** Tn5 insertions within tRNA gene bodies (±100 bp around midpoint)

**Config:**

```yaml
count_matrices:
  trna_counts:
    enabled: true
    flank_size: 100 # ±100 bp around midpoint
    max_size: 1000 # Feature size filter
```

**Workflow:**

```
tRNA GTF
  → build_tn5_midpoint_bins.py (trna_gene_flank_size, trna_gene_max_size)
  → per-sample count_tn5_sites_in_bins.py
  → aggregate build_tn5_count_matrix.py
  → host.trna.tn5_count_matrix.tsv
```

**Status:** Rules already partially in place (build_host_tn5_trna_gene_bins). Just needs config toggle.

---

## 3. POL3 Count Matrix (RNA5S rRNA)

**Description:** Tn5 insertions within Pol III RNA5S rRNA genes (full feature, max 1000 bp)

**Config:**

```yaml
count_matrices:
  pol3_counts:
    enabled: true
    flank_size: -1 # Full feature span (no flanking)
    max_size: 1000 # Feature size filter
    # Generates 3 separate matrices:
    # - pol3_T1 (T1_Genes)
    # - pol3_T2 (T2_noAlu)
    # - pol3_T3 (T3_Genes)
```

**Workflow (per Pol3 type - T1, T2, T3):**

```
T1_Genes.GTF (from config.pol3_gtf)
  → build_tn5_midpoint_bins.py (flank_size=-1, max_size=1000)
  → per-sample count_tn5_sites_in_bins.py
  → aggregate build_tn5_count_matrix.py
  → host.pol3_T1.tn5_count_matrix.tsv
  → host.pol3_T2.tn5_count_matrix.tsv
  → host.pol3_T3.tn5_count_matrix.tsv
```

**Rules to add:**

- `build_host_pol3_XX_bins` (X3: T1, T2, T3)
- `count_tn5_host_pol3_XX_bins` (per sample, X3)
- `aggregate_tn5_host_pol3_XX_count_matrix` (X3)

---

## 4. Repeat Element Count Matrices

**Description:** Tn5 insertions within repeat elements (full feature, max 1000 bp)

**Config:**

```yaml
count_matrices:
  repeat_element_counts:
    enabled: true
    flank_size: -1 # Full feature span
    max_size: 1000 # Feature size filter
    # Generates 6 separate matrices:
    # - SINE_Alu
    # - SINE_MIR
    # - LINE_L1
    # - LINE_L2
    # - LTR
    # - other_repeat_elements
```

**Workflow (per repeat type - 6 total):**

```
SINE_Alu.GTF (from config.repeat_elements_gtf)
  → build_tn5_midpoint_bins.py (flank_size=-1, max_size=1000)
  → per-sample count_tn5_sites_in_bins.py
  → aggregate build_tn5_count_matrix.py
  → host.repeat_SINE_Alu.tn5_count_matrix.tsv
  → [5 more for other types]
```

**Rules to add:**

- `build_host_repeat_element_XX_bins` (X6: one per repeat type)
- `count_tn5_host_repeat_element_XX_bins` (per sample, X6)
- `aggregate_tn5_host_repeat_element_XX_count_matrix` (X6)

---

## 5. Viral Genome Bin Count Matrix

**Description:** Tn5 insertions in 200 bp fixed non-overlapping bins across each viral genome

**Config:**

```yaml
count_matrices:
  viral_genome_counts:
    enabled: true
    bin_size: 200 # Fixed bin size
```

**Workflow (per virus):**

```
Virus FASTA + chromsizes
  → build_viral_genome_bins.py (bin_size=200)
  → per-sample count_tn5_sites_in_bins.py
  → aggregate build_tn5_count_matrix.py
  → virus.tn5_bin_count_matrix.200bp.tsv (per virus)
```

**Scripts to add:**

- `build_viral_genome_bins.py` - NEW: creates fixed-size non-overlapping bins across genome

**Rules to add:**

- `build_virus_XX_bins` (per virus in VIRUS_LIST)
- `count_tn5_virus_XX_bins` (per sample, per virus)
- `aggregate_tn5_virus_XX_count_matrix` (per virus)

---

## 6. rRNA Gene Count Matrix

**Description:** Tn5 insertions around TSS of rRNA genes (±250 bp)

**Config:**

```yaml
count_matrices:
  rrna_counts:
    enabled: true
    flank_size: 250 # ±250 bp around TSS
```

**Workflow:**

```
rRNA GTF (from future config.rrna_gtf)
  → build_tn5_tss_bins.py (flank_size=250)
  → per-sample count_tn5_sites_in_bins.py
  → aggregate build_tn5_count_matrix.py
  → host.rrna.tn5_count_matrix.tsv
```

**Rules to add:**

- `build_host_rrna_tss_bins`
- `count_tn5_host_rrna_bins` (per sample)
- `aggregate_tn5_host_rrna_count_matrix`

**Prerequisites:**

- Need to configure rRNA GTF in config.yaml (analogous to tRNA GTF)

---

## Implementation Phases

### Phase 1: Config Structure

1. Add `count_matrices` section to config.yaml with toggles and parameters
2. Extract toggles and parameters in init.smk or tn5motif.smk

### Phase 2: New Scripts

1. `build_viral_genome_bins.py` - creates fixed-size bins
2. Update config to include rRNA GTF paths

### Phase 3: Snakemake Rules

1. Add conditional rules for each matrix type (if enabled)
2. Create bin-building rules
3. Create per-sample counting rules
4. Create aggregation rules

### Phase 4: Output Organization

```
results/tn5_motif/
├── gene.tn5_gene_count_matrix.250bp.tsv
├── trna.tn5_trna_count_matrix.100bp.tsv
├── pol3_t1.tn5_pol3_t1_count_matrix.tsv
├── pol3_t2.tn5_pol3_t2_count_matrix.tsv
├── pol3_t3.tn5_pol3_t3_count_matrix.tsv
├── repeat_sine_alu.tn5_repeat_sine_alu_count_matrix.tsv
├── repeat_sine_mir.tn5_repeat_sine_mir_count_matrix.tsv
├── repeat_line_l1.tn5_repeat_line_l1_count_matrix.tsv
├── repeat_line_l2.tn5_repeat_line_l2_count_matrix.tsv
├── repeat_ltr.tn5_repeat_ltr_count_matrix.tsv
├── repeat_other.tn5_repeat_other_count_matrix.tsv
├── virus_rinderpest.tn5_bin_count_matrix.200bp.tsv
├── virus_ebola.tn5_bin_count_matrix.200bp.tsv
└── rrna.tn5_rrna_count_matrix.250bp.tsv
```

---

## Key Reusable Components

**Scripts (already exist or will be minimal new):**

- `build_tn5_tss_bins.py` - TSS-based (gene, rRNA)
- `build_tn5_midpoint_bins.py` - midpoint-based (tRNA, Pol3, repeat elements)
- `build_viral_genome_bins.py` - NEW: fixed-size bins
- `count_tn5_sites_in_bins.py` - reused for all (already handles BAM counting)
- `build_tn5_count_matrix.py` - reused for all (already handles aggregation)

**Config parameters (all per-matrix):**

- `enabled` - toggle on/off
- `flank_size` - for TSS/midpoint modes (-1 = full feature)
- `max_size` - feature size filter
- `bin_size` - fixed bin mode (viral genomes)

---

## Conditional Rule Logic

```python
# In init.smk or tn5motif.smk
COUNT_MATRIX_TYPES = {
    'gene': config.get('count_matrices', {}).get('gene_counts', {}).get('enabled', True),
    'trna': config.get('count_matrices', {}).get('trna_counts', {}).get('enabled', True),
    'pol3': config.get('count_matrices', {}).get('pol3_counts', {}).get('enabled', False),
    'repeat_elements': config.get('count_matrices', {}).get('repeat_element_counts', {}).get('enabled', False),
    'viral': config.get('count_matrices', {}).get('viral_genome_counts', {}).get('enabled', False),
    'rrna': config.get('count_matrices', {}).get('rrna_counts', {}).get('enabled', False),
}

# Then in rules:
if COUNT_MATRIX_TYPES['gene']:
    rule build_host_gene_tss_bins:
        ...
```

---

## Next Steps

1. **Design config.yaml structure** for count_matrices section
2. **Write build_viral_genome_bins.py** script
3. **Update tn5motif.smk** with conditional rules and parameters
4. **Test** with subset of samples
5. **Document** final matrix outputs and interpretation

---

## CLARIFICATIONS

### chrR GTFs are rRNA GTFs

The `chrR_gtf` (Chromatin Reference) in config.yaml already contains rRNA genes for all
three genomes (hg38, mm39, hs1). These are the same GTFs used for rRNA count matrix.
No additional rRNA GTF configuration needed.

### Actual Prerequisites (Simplified)

1. Add `count_matrices` config section with toggles and parameters
2. Write `build_viral_genome_bins.py` script
3. Add ~20 Snakemake rules (3 per matrix type: build → count → aggregate)

That's it. All GTF files already exist and are configured.
