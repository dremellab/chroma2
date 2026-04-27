"""Unit tests for build_host_trna_gene_bins.py"""

import sys
import tempfile
import textwrap
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from build_host_trna_gene_bins import (
    collect_trna_genes,
    gene_center_from_fields,
    load_chromsizes,
    load_host_contigs,
    write_bins,
)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _fields(chrom, start, end, strand, gene_type, gene_id):
    """Build a minimal GTF-like field list (0-indexed, 9 fields)."""
    attrs = f'gene_id "{gene_id}"; gene_type "{gene_type}";'
    return [chrom, ".", "gene", str(start), str(end), ".", strand, ".", attrs]


def _write_regions(tmp_path, contigs):
    p = tmp_path / "ref.fa.regions"
    with open(p, "w") as f:
        for name, chrs in contigs.items():
            f.write(f"{name}\t{' '.join(chrs)}\n")
    return str(p)


def _write_chromsizes(tmp_path, sizes):
    p = tmp_path / "chrom.sizes"
    with open(p, "w") as f:
        for chrom, size in sizes.items():
            f.write(f"{chrom}\t{size}\n")
    return str(p)


def _write_gtf(tmp_path, records):
    """records: list of (chrom, start_1based, end_1based, strand, gene_type, gene_id)"""
    p = tmp_path / "test.gtf"
    with open(p, "w") as f:
        for chrom, start, end, strand, gtype, gid in records:
            attrs = f'gene_id "{gid}"; gene_type "{gtype}";'
            f.write(f"{chrom}\t.\tgene\t{start}\t{end}\t.\t{strand}\t.\t{attrs}\n")
    return str(p)


# ---------------------------------------------------------------------------
# gene_center_from_fields
# ---------------------------------------------------------------------------


def test_gene_center_plus_strand():
    # GTF: start=1001 end=1100 (1-based) → 0-based: 1000..1099 → center = 1049
    fields = _fields("chr1", 1001, 1100, "+", "tRNA", "g1")
    assert gene_center_from_fields(fields) == 1049


def test_gene_center_minus_strand():
    # Same coordinates, minus strand → center must be the same (strand-independent)
    fields_plus = _fields("chr1", 1001, 1100, "+", "tRNA", "g1")
    fields_minus = _fields("chr1", 1001, 1100, "-", "tRNA", "g1")
    assert gene_center_from_fields(fields_plus) == gene_center_from_fields(fields_minus)


def test_gene_center_single_base():
    # start == end → center == start (0-based)
    fields = _fields("chr1", 500, 500, "+", "tRNA", "g1")
    assert gene_center_from_fields(fields) == 499


def test_gene_center_even_span():
    # start=1 end=4 → 0-based 0..3 → center = (0+3)//2 = 1
    fields = _fields("chr1", 1, 4, "+", "tRNA", "g1")
    assert gene_center_from_fields(fields) == 1


# ---------------------------------------------------------------------------
# collect_trna_genes: filtering
# ---------------------------------------------------------------------------


def test_trna_gene_type_included(tmp_path):
    gtf = _write_gtf(tmp_path, [("chr1", 1001, 1100, "+", "tRNA", "g1")])
    contigs = {"ref": ["chr1"]}
    regions = _write_regions(tmp_path, contigs)
    host_contigs = load_host_contigs(regions)
    genes = collect_trna_genes(gtf, host_contigs, {"trna"})
    assert "g1" in genes


def test_protein_coding_excluded(tmp_path):
    gtf = _write_gtf(tmp_path, [("chr1", 1001, 1100, "+", "protein_coding", "g1")])
    contigs = {"ref": ["chr1"]}
    regions = _write_regions(tmp_path, contigs)
    host_contigs = load_host_contigs(regions)
    genes = collect_trna_genes(gtf, host_contigs, {"trna"})
    assert "g1" not in genes


def test_case_insensitive_filter_upper(tmp_path):
    gtf = _write_gtf(tmp_path, [("chr1", 1001, 1100, "+", "TRNA", "g1")])
    contigs = {"ref": ["chr1"]}
    regions = _write_regions(tmp_path, contigs)
    host_contigs = load_host_contigs(regions)
    genes = collect_trna_genes(gtf, host_contigs, {"trna"})
    assert "g1" in genes


def test_case_insensitive_filter_mixed(tmp_path):
    gtf = _write_gtf(tmp_path, [("chr1", 1001, 1100, "+", "tRNA_gene", "g1")])
    contigs = {"ref": ["chr1"]}
    regions = _write_regions(tmp_path, contigs)
    host_contigs = load_host_contigs(regions)
    genes = collect_trna_genes(gtf, host_contigs, {"trna_gene"})
    assert "g1" in genes


def test_off_host_contig_excluded(tmp_path):
    gtf = _write_gtf(tmp_path, [("chrUn_random", 1001, 1100, "+", "tRNA", "g1")])
    contigs = {"ref": ["chr1"]}
    regions = _write_regions(tmp_path, contigs)
    host_contigs = load_host_contigs(regions)
    genes = collect_trna_genes(gtf, host_contigs, {"trna"})
    assert "g1" not in genes


# ---------------------------------------------------------------------------
# write_bins: bin coordinates
# ---------------------------------------------------------------------------


def test_bin_centered_correctly(tmp_path):
    # gene: chr1:1001-1100 (1-based) → center=1049 (0-based), flank=100
    # expected: start=949, end=1150
    gtf = _write_gtf(tmp_path, [("chr1", 1001, 1100, "+", "tRNA", "g1")])
    contigs = {"ref": ["chr1"]}
    regions = _write_regions(tmp_path, contigs)
    chromsizes_file = _write_chromsizes(tmp_path, {"chr1": 200000})
    host_contigs = load_host_contigs(regions)
    chromsizes = load_chromsizes(chromsizes_file)
    genes = collect_trna_genes(gtf, host_contigs, {"trna"})
    out = tmp_path / "out.bed"
    write_bins(str(out), chromsizes, 100, genes)
    row = out.read_text().strip().split("\t")
    assert row[0] == "chr1"
    assert int(row[1]) == 949  # 1049 - 100
    assert int(row[2]) == 1150  # 1049 + 100 + 1


def test_bin_clamped_at_zero(tmp_path):
    # gene near start: center=50, flank=100 → start would be -50 → clamped to 0
    gtf = _write_gtf(tmp_path, [("chr1", 1, 101, "+", "tRNA", "g1")])
    contigs = {"ref": ["chr1"]}
    regions = _write_regions(tmp_path, contigs)
    chromsizes_file = _write_chromsizes(tmp_path, {"chr1": 200000})
    host_contigs = load_host_contigs(regions)
    chromsizes = load_chromsizes(chromsizes_file)
    genes = collect_trna_genes(gtf, host_contigs, {"trna"})
    out = tmp_path / "out.bed"
    write_bins(str(out), chromsizes, 100, genes)
    row = out.read_text().strip().split("\t")
    assert int(row[1]) == 0


def test_bin_clamped_at_chrom_end(tmp_path):
    # gene near end of chr: center=9950, flank=100, chromsize=10000
    # end would be 10051 → clamped to 10000
    gtf = _write_gtf(tmp_path, [("chr1", 9901, 10000, "+", "tRNA", "g1")])
    contigs = {"ref": ["chr1"]}
    regions = _write_regions(tmp_path, contigs)
    chromsizes_file = _write_chromsizes(tmp_path, {"chr1": 10000})
    host_contigs = load_host_contigs(regions)
    chromsizes = load_chromsizes(chromsizes_file)
    genes = collect_trna_genes(gtf, host_contigs, {"trna"})
    out = tmp_path / "out.bed"
    write_bins(str(out), chromsizes, 100, genes)
    row = out.read_text().strip().split("\t")
    assert int(row[2]) == 10000


def test_bin_id_uses_center_suffix(tmp_path):
    gtf = _write_gtf(tmp_path, [("chr1", 1001, 1100, "+", "tRNA", "g1")])
    contigs = {"ref": ["chr1"]}
    regions = _write_regions(tmp_path, contigs)
    chromsizes_file = _write_chromsizes(tmp_path, {"chr1": 200000})
    host_contigs = load_host_contigs(regions)
    chromsizes = load_chromsizes(chromsizes_file)
    genes = collect_trna_genes(gtf, host_contigs, {"trna"})
    out = tmp_path / "out.bed"
    write_bins(str(out), chromsizes, 100, genes)
    row = out.read_text().strip().split("\t")
    assert row[3].endswith("|center")


# ---------------------------------------------------------------------------
# integration smoke test
# ---------------------------------------------------------------------------


def test_integration_two_genes(tmp_path):
    records = [
        ("chr1", 1001, 1100, "+", "tRNA", "g1"),
        ("chr2", 5001, 5100, "-", "tRNA_gene", "g2"),
    ]
    gtf = _write_gtf(tmp_path, records)
    contigs = {"ref": ["chr1", "chr2"]}
    regions = _write_regions(tmp_path, contigs)
    chromsizes_file = _write_chromsizes(tmp_path, {"chr1": 200000, "chr2": 200000})
    host_contigs = load_host_contigs(regions)
    chromsizes = load_chromsizes(chromsizes_file)
    genes = collect_trna_genes(gtf, host_contigs, {"trna", "trna_gene"})
    assert len(genes) == 2
    out = tmp_path / "out.bed"
    write_bins(str(out), chromsizes, 100, genes)
    lines = out.read_text().strip().splitlines()
    assert len(lines) == 2
    chroms = {line.split("\t")[0] for line in lines}
    assert chroms == {"chr1", "chr2"}
