"""Unit tests for extract_tn5_motifs.py's Tn5 cut-site extraction.

Requires pysam; run inside the pysam apptainer/singularity image, e.g.:
  apptainer exec .../pysam_0.22.1.sif python3 -m pytest workflow/tests/test_extract_tn5_motifs.py
"""

import io
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from extract_tn5_motifs import (  # noqa: E402
    CutSite,
    alignment_passes_filters,
    cut_sites_from_record,
    fetch_window_sequence,
    parse_args,
    single_cut_site,
    write_flank_bed,
)

SCRIPT_PATH = Path(__file__).parent.parent / "scripts" / "extract_tn5_motifs.py"

HEADER = pysam.AlignmentHeader.from_dict(
    {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 10000}]}
)
# A short contig, sized like a real viral genome, for boundary-condition tests
# where a cut site's offset math can plausibly push it past the contig end.
SHORT_HEADER = pysam.AlignmentHeader.from_dict(
    {"HD": {"VN": "1.6"}, "SQ": [{"SN": "hsv1", "LN": 1010}]}
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


def test_single_cut_site_returns_empty_when_shifted_cut_exceeds_chrom_length():
    # Regression test: single_cut_site() only checked cut < 0, never the
    # upper boundary. A forward-strand read whose alignment ends exactly at
    # the contig length still pushes the +4 Tn5-offset cut site past it --
    # this used to crash downstream (--bin-counts-only indexes bin_counts[chrom]
    # with no bounds check) instead of being dropped like the negative case.
    rec = pysam.AlignedSegment(SHORT_HEADER)
    rec.query_name = "r1"
    rec.flag = 0  # forward strand
    rec.reference_id = 0
    rec.reference_start = 1007  # SHORT_HEADER's hsv1 is 1010bp
    rec.cigarstring = "3M"  # reference_end == 1010, a valid alignment
    assert rec.reference_end == 1010

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


def test_cut_sites_from_record_labels_mates_correctly_when_processed_independently():
    # Regression test: mate labeling used to rely on grouping adjacent
    # same-query_name BAM records, which only works on name-sorted input. The
    # real --bam inputs are coordinate-sorted, so mates are essentially never
    # adjacent -- this simulates that reality by calling cut_sites_from_record
    # on each mate as a fully independent, non-adjacent call (no shared state,
    # no batching) and asserting the SAM is_read1/is_read2 flags alone are
    # enough to label each site correctly.
    #
    # Proper pair, FR orientation, with a gap between the mates (the ATAC-seq
    # norm for fragments longer than one read): R1 forward at 1000-1050, R2
    # reverse at 1100-1150. This is also the exact shape of the previously-fixed
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

    r1_sites = cut_sites_from_record(
        r1, mapq_min=0, exclude_secondary=False, exclude_supplementary=False
    )
    r2_sites = cut_sites_from_record(
        r2, mapq_min=0, exclude_secondary=False, exclude_supplementary=False
    )

    assert len(r1_sites) == 1 and len(r2_sites) == 1
    assert r1_sites[0].name == "pair1|read1"
    assert r1_sites[0].start == 1004
    assert r1_sites[0].strand == "+"
    assert r2_sites[0].name == "pair1|read2"
    assert r2_sites[0].start == 1145
    assert r2_sites[0].strand == "-"
    # Guard against regressing to R1-only derivation for the mate's site.
    assert r2_sites[0].start != r1.reference_end - 5


def test_cut_sites_from_record_labels_unpaired_read_as_read1():
    rec = _make_record("single_end", flag=0, reference_start=1000, cigarstring="50M")

    sites = cut_sites_from_record(
        rec, mapq_min=0, exclude_secondary=False, exclude_supplementary=False
    )

    assert len(sites) == 1
    assert sites[0].name == "single_end|read1"


def test_write_flank_bed_clamps_negative_start_near_chromosome_beginning():
    # Regression test: write_flank_bed() used to write start - flank_size
    # unclamped, producing a negative BED coordinate for a cut site near the
    # start of a chromosome.
    site = CutSite("chr1", start=2, end=3, name="r1|read1", score=60, strand="+")
    fh = io.StringIO()

    write_flank_bed(fh, site, flank_size=10, chrom_len=10000)

    fields = fh.getvalue().rstrip("\n").split("\t")
    assert int(fields[1]) == 0  # clamped, not 2 - 10 == -8
    assert int(fields[2]) == 13  # end + flank_size, well within chrom_len


def test_write_flank_bed_clamps_end_past_chromosome_length():
    site = CutSite("chr1", start=995, end=996, name="r1|read1", score=60, strand="+")
    fh = io.StringIO()

    write_flank_bed(fh, site, flank_size=10, chrom_len=1000)

    fields = fh.getvalue().rstrip("\n").split("\t")
    assert int(fields[1]) == 985  # start - flank_size, well within chrom_len
    assert int(fields[2]) == 1000  # clamped, not 996 + 10 == 1006


def test_fetch_window_sequence_respects_passed_in_chrom_len(tmp_path):
    fasta_path = tmp_path / "ref.fa"
    fasta_path.write_text(">chr1\n" + "ACGT" * 5 + "\n")  # chr1, 20bp
    pysam.faidx(str(fasta_path))

    with pysam.FastaFile(str(fasta_path)) as fasta:
        chrom_len = fasta.get_reference_length("chr1")
        assert chrom_len == 20

        # cut_start=2, flank_size=5 -> left=-3, out of bounds -> None, matching
        # the site write_flank_bed's sibling test above clamps instead of drops.
        assert (
            fetch_window_sequence(
                fasta,
                "chr1",
                cut_start=2,
                strand="+",
                flank_size=5,
                chrom_len=chrom_len,
            )
            is None
        )

        seq = fetch_window_sequence(
            fasta,
            "chr1",
            cut_start=10,
            strand="+",
            flank_size=5,
            chrom_len=chrom_len,
        )
        assert seq is not None
        assert len(seq) == 11


def test_bin_counts_only_skips_out_of_bounds_cut_site_without_crashing(tmp_path):
    # End-to-end regression test for the --bin-counts-only IndexError: a read
    # whose computed cut site lands past the (short, virus-genome-like)
    # contig's length must be silently dropped, not crash the whole run.
    chrom = "hsv1"
    # chrom_len is an exact multiple of bin_size so the boundary read's
    # out-of-range cut site maps to bin_idx == len(bin_counts[chrom]) -- an
    # actual IndexError on the unpatched code, not just an incorrect count
    # absorbed into a trailing partial bin (which a non-multiple chrom_len
    # would silently do instead, masking the crash this test exists to catch).
    chrom_len = 1000
    bin_size = 100

    fasta_path = tmp_path / "hsv1.fa"
    fasta_path.write_text(f">{chrom}\n{'A' * chrom_len}\n")
    pysam.faidx(str(fasta_path))

    regions_path = tmp_path / "hsv1.fa.regions"
    regions_path.write_text(f"{chrom}\t{chrom}\n")

    header_dict = {"HD": {"VN": "1.6"}, "SQ": [{"SN": chrom, "LN": chrom_len}]}
    header = pysam.AlignmentHeader.from_dict(header_dict)

    def make_record(query_name, reference_start, cigarstring):
        rec = pysam.AlignedSegment(header)
        rec.query_name = query_name
        rec.flag = 0  # forward strand, not paired
        rec.reference_id = 0
        rec.reference_start = reference_start
        rec.mapping_quality = 60
        rec.cigarstring = cigarstring
        rec.next_reference_id = 0
        return rec

    # Valid read: cut site at 504, inside bin 5 (500-600).
    valid = make_record("valid_read", reference_start=500, cigarstring="50M")
    # Boundary read: alignment ends exactly at chrom_len (a valid alignment),
    # but the +4 Tn5 offset pushes the cut site to chrom_len+1 -- out of range.
    boundary = make_record(
        "boundary_read", reference_start=chrom_len - 3, cigarstring="3M"
    )
    assert boundary.reference_end == chrom_len

    bam_path = tmp_path / "reads.bam"
    with pysam.AlignmentFile(str(bam_path), "wb", header=header_dict) as bam:
        bam.write(valid)
        bam.write(boundary)

    bin_counts_output = tmp_path / "bin_counts.tsv"
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--bam",
            str(bam_path),
            "--fasta",
            str(fasta_path),
            "--sample",
            "sample_a",
            "--scenario-name",
            "macs2",
            "--outdir",
            str(tmp_path / "out"),
            "--virus-regions",
            f"{chrom}={regions_path}",
            "--virus-bin-size",
            str(bin_size),
            "--bin-counts-only",
            "--bin-counts-output",
            str(bin_counts_output),
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0, result.stderr

    rows = {}
    lines = bin_counts_output.read_text().splitlines()
    for line in lines[1:]:
        fields = line.split("\t")
        rows[int(fields[5])] = int(fields[-1])  # bin_index -> tn5_site_count

    total_bins = -(-chrom_len // bin_size)  # ceil division, matches the script
    assert len(rows) == total_bins
    assert rows[5] == 1  # valid_read's cut site (504) landed in bin 5 (500-600)
    assert sum(rows.values()) == 1  # boundary_read contributed nothing


def test_parse_args_rejects_negative_mapq_min(capsys):
    # Every other required arg is supplied so the only possible source of a
    # SystemExit is the --mapq-min type check itself -- otherwise argparse's
    # unrelated "missing required arguments" error also exits with code 2,
    # and this test would pass even if the negative-value check were reverted.
    with pytest.raises(SystemExit) as exc_info:
        parse_args(
            [
                "--bam",
                "in.bam",
                "--fasta",
                "ref.fa",
                "--sample",
                "sample_a",
                "--scenario-name",
                "macs2",
                "--outdir",
                "/tmp/out",
                "--mapq-min",
                "-1",
            ]
        )
    assert exc_info.value.code == 2
    assert "--mapq-min" in capsys.readouterr().err


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
