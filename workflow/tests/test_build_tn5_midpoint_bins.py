"""Unit tests for build_tn5_midpoint_bins.py's write_bins() size filter.

Pure stdlib -- no pysam/pandas needed, runs under any Python.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from build_tn5_midpoint_bins import write_bins  # noqa: E402


def _meta(chrom, start, end, strand="+"):
    return {
        "chrom": chrom,
        "strand": strand,
        "center": (start + end) // 2,
        "start": start,
        "end": end,
        "gene_name": "gene_name",
        "gene_type": "protein_coding",
    }


def _perline_meta(
    chrom, start, end, feature_type, gene_type, repeat_name="", strand="+"
):
    return {
        "chrom": chrom,
        "strand": strand,
        "center": (start + end) // 2,
        "start": start,
        "end": end,
        "gene_name": "gene_name",
        "gene_type": gene_type,
        "feature_type": feature_type,
        "repeat_name": repeat_name,
        "repeat_family": "",
        "sw_score": "",
    }


def _bin_ids_written(path, chromsizes, flank_size, max_size, genes):
    write_bins(str(path), chromsizes, flank_size, max_size, genes)
    return {line.split("\t")[3] for line in path.read_text().splitlines() if line}


def test_max_size_filter_keeps_feature_exactly_at_the_limit(tmp_path):
    # meta start/end are 0-based *inclusive*, so a feature spanning
    # start=1000..end=1999 is exactly 1000bp long (1999 - 1000 + 1).
    genes = {"g_at_limit": _meta("chr1", 1000, 1999)}
    chromsizes = {"chr1": 1_000_000}

    bin_ids = _bin_ids_written(
        tmp_path / "out.bed", chromsizes, flank_size=-1, max_size=1000, genes=genes
    )

    assert bin_ids == {"g_at_limit_gene_name"}


def test_max_size_filter_drops_feature_one_bp_over_the_limit(tmp_path):
    # Regression test: feature_size used to be computed as end - start
    # (missing the +1 for inclusive coordinates), so a feature that's
    # actually max_size+1 bp long was incorrectly computed as exactly
    # max_size and retained instead of filtered out.
    genes = {
        "g_at_limit": _meta("chr1", 1000, 1999),  # 1000bp -- must be kept
        "g_over_limit": _meta("chr1", 1000, 2000),  # 1001bp -- must be dropped
    }
    chromsizes = {"chr1": 1_000_000}

    bin_ids = _bin_ids_written(
        tmp_path / "out.bed", chromsizes, flank_size=-1, max_size=1000, genes=genes
    )

    assert bin_ids == {"g_at_limit_gene_name"}


def test_max_size_filter_unaffected_well_above_and_below_the_limit(tmp_path):
    genes = {
        "g_small": _meta("chr1", 1000, 1099),  # 100bp, well under
        "g_huge": _meta("chr1", 1000, 4999),  # 4000bp, well over
    }
    chromsizes = {"chr1": 1_000_000}

    bin_ids = _bin_ids_written(
        tmp_path / "out.bed", chromsizes, flank_size=-1, max_size=1000, genes=genes
    )

    assert bin_ids == {"g_small_gene_name"}


def test_perline_bin_id_uses_feature_type_not_gene_type(tmp_path):
    # Regression test: bin_id used to be built from meta["gene_type"] instead
    # of meta["feature_type"], contradicting the documented
    # feature_type_chrom_start_end_name format. The two happen to be equal
    # whenever a GTF line has no gene_type/gene_biotype/transcript_type
    # attribute (gene_type falls back to the raw feature column), which is
    # why this went unnoticed -- so this test gives the two fields different
    # values, as a per-line GTF row carrying an explicit biotype attribute
    # would. Two rows share coordinates/repeat_name but have different
    # feature_type values and the SAME gene_type; with the bug they'd
    # collide on one bin_id, with the fix they're distinct.
    genes = {
        "g1": _perline_meta(
            "chr1",
            1000,
            1099,
            feature_type="exon",
            gene_type="protein_coding",
            repeat_name="AluY",
        ),
        "g2": _perline_meta(
            "chr1",
            1000,
            1099,
            feature_type="intron",
            gene_type="protein_coding",
            repeat_name="AluY",
        ),
    }
    chromsizes = {"chr1": 1_000_000}

    bin_ids = _bin_ids_written(
        tmp_path / "out.bed", chromsizes, flank_size=-1, max_size=-1, genes=genes
    )

    assert bin_ids == {"exon_chr1_1001_1100_AluY", "intron_chr1_1001_1100_AluY"}


if __name__ == "__main__":
    import pytest

    sys.exit(pytest.main([__file__, "-v"]))
