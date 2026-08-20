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


if __name__ == "__main__":
    import pytest

    sys.exit(pytest.main([__file__, "-v"]))
