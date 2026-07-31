"""Unit tests for extract_tn5_motifs.py's Tn5 cut-site extraction.

Requires pysam; run inside the pysam apptainer/singularity image, e.g.:
  apptainer exec .../pysam_0.22.1.sif python3 -m pytest workflow/tests/test_extract_tn5_motifs.py
"""

import sys
from pathlib import Path

import pysam
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from extract_tn5_motifs import (  # noqa: E402
    alignment_passes_filters,
    cut_sites_from_record,
    cut_sites_from_records,
    single_cut_site,
)

HEADER = pysam.AlignmentHeader.from_dict(
    {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 10000}]}
)


def _make_record(
    query_name,
    flag,
    reference_start,
    cigarstring,
    mapping_quality=60,
    next_reference_start=0,
    template_length=0,
):
    rec = pysam.AlignedSegment(HEADER)
    rec.query_name = query_name
    rec.flag = flag
    rec.reference_id = 0
    rec.reference_start = reference_start
    rec.mapping_quality = mapping_quality
    rec.cigarstring = cigarstring
    rec.next_reference_id = 0
    rec.next_reference_start = next_reference_start
    rec.template_length = template_length
    return rec


def test_single_cut_site_forward_strand():
    rec = _make_record("r1", flag=0, reference_start=1000, cigarstring="50M")

    sites = single_cut_site(rec, 1)

    assert len(sites) == 1
    site = sites[0]
    assert site.chrom == "chr1"
    assert site.start == 1004
    assert site.end == 1005
    assert site.strand == "+"
    assert site.name == "r1|read1"


def test_single_cut_site_reverse_strand():
    rec = _make_record("r1", flag=16, reference_start=1000, cigarstring="50M")

    sites = single_cut_site(rec, 1)

    assert len(sites) == 1
    site = sites[0]
    assert site.start == rec.reference_end - 5 == 1045
    assert site.strand == "-"


def test_single_cut_site_returns_empty_when_shifted_cut_is_negative():
    rec = _make_record("r1", flag=16, reference_start=0, cigarstring="3M")

    assert single_cut_site(rec, 1) == []


def test_alignment_passes_filters_respects_mapq_secondary_supplementary():
    rec = _make_record("r1", flag=0, reference_start=1000, cigarstring="50M")
    rec.mapping_quality = 5

    assert not alignment_passes_filters(
        rec, mapq_min=10, exclude_secondary=False, exclude_supplementary=False
    )
    assert alignment_passes_filters(
        rec, mapq_min=5, exclude_secondary=False, exclude_supplementary=False
    )

    rec.mapping_quality = 60
    rec.flag = 0x100  # secondary
    assert not alignment_passes_filters(
        rec, mapq_min=0, exclude_secondary=True, exclude_supplementary=False
    )
    assert alignment_passes_filters(
        rec, mapq_min=0, exclude_secondary=False, exclude_supplementary=False
    )


def test_cut_sites_from_record_returns_empty_for_unmapped():
    rec = _make_record("r1", flag=4, reference_start=1000, cigarstring="50M")

    assert (
        cut_sites_from_record(
            rec, mapq_min=0, exclude_secondary=False, exclude_supplementary=False
        )
        == []
    )


def test_cut_sites_from_records_derives_each_mates_site_from_its_own_alignment():
    # Proper pair, FR orientation, with a gap between the mates (the ATAC-seq
    # norm for fragments longer than one read): R1 forward at 1000-1050, R2
    # reverse at 1100-1150. This is the exact shape of the previously-fixed
    # R1-only bug: deriving both cut sites from R1 alone would put the "far"
    # site at R1.reference_end - 5 == 1045, not R2's real 1145.
    r1 = _make_record(
        "pair1",
        flag=99,
        reference_start=1000,
        cigarstring="50M",
        next_reference_start=1100,
        template_length=150,
    )
    r2 = _make_record(
        "pair1",
        flag=147,
        reference_start=1100,
        cigarstring="50M",
        next_reference_start=1000,
        template_length=-150,
    )

    sites = cut_sites_from_records(
        [r1, r2], mapq_min=0, exclude_secondary=False, exclude_supplementary=False
    )

    assert len(sites) == 2
    by_name = {site.name: site for site in sites}
    assert by_name["pair1|read1"].start == 1004
    assert by_name["pair1|read1"].strand == "+"
    assert by_name["pair1|read2"].start == 1145
    assert by_name["pair1|read2"].strand == "-"
    # Guard against regressing to R1-only derivation for the mate's site.
    assert by_name["pair1|read2"].start != r1.reference_end - 5


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
