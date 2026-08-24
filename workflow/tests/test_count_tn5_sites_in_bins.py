"""Unit tests for count_tn5_sites_in_bins.py.

The pure bin/counting logic (load_bins, overlapping_bin_ids, write_counts,
validate_output) needs no BAM library and runs under any Python. The
end-to-end test drives the real script against a synthetic BAM and needs
pysam; run inside the pysam apptainer/singularity image, e.g.:
  apptainer exec .../pysam_0.22.1.sif python3 -m pytest workflow/tests/test_count_tn5_sites_in_bins.py
"""

import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from count_tn5_sites_in_bins import (  # noqa: E402
    get_nh_value,
    load_bins,
    overlapping_bin_ids,
    parse_args,
    validate_output,
    write_counts,
)
from gtf_common import (  # noqa: E402
    BIN_CORE_COLUMNS,
    GROUPED_BIN_METADATA_COLUMNS,
    PERLINE_BIN_METADATA_COLUMNS,
)

SCRIPT_PATH = Path(__file__).parent.parent / "scripts" / "count_tn5_sites_in_bins.py"


def _write_bins_bed(path, rows):
    with open(path, "w", encoding="utf-8") as handle:
        for chrom, start, end, bin_id in rows:
            handle.write(
                "\t".join(
                    [
                        chrom,
                        str(start),
                        str(end),
                        bin_id,
                        f"{bin_id}_gene",
                        f"{bin_id}_name",
                        "protein_coding",
                        "+",
                        str(start),
                    ]
                )
                + "\n"
            )


def test_load_bins_indexes_by_chrom_and_tracks_max_width(tmp_path):
    bed = tmp_path / "bins.bed"
    _write_bins_bed(
        bed,
        [
            ("chr1", 1000, 1010, "bin1"),
            ("chr1", 1140, 1150, "bin2"),
            ("chr2", 500, 600, "bin3"),
        ],
    )

    bins_by_chrom, starts_by_chrom, counts, max_width = load_bins(str(bed))

    assert set(bins_by_chrom.keys()) == {"chr1", "chr2"}
    assert starts_by_chrom["chr1"] == [1000, 1140]
    assert counts == {"bin1": 0, "bin2": 0, "bin3": 0}
    assert max_width == 100


def test_overlapping_bin_ids_distinguishes_near_and_far_bins(tmp_path):
    bed = tmp_path / "bins.bed"
    _write_bins_bed(
        bed,
        [
            ("chr1", 1000, 1010, "near_bin"),
            ("chr1", 1140, 1150, "far_bin"),
        ],
    )
    bins_by_chrom, starts_by_chrom, _counts, max_width = load_bins(str(bed))

    # The R1-side cut site (1004) must land in near_bin only; the R2-side cut
    # site (1145) must land in far_bin only -- this is the bin-assignment half
    # of the previously-fixed R1-only cut-site bug (extract_tn5_motifs.py):
    # a regression there would put the far site around 1045, which falls in
    # neither bin here.
    assert overlapping_bin_ids(
        "chr1", 1004, bins_by_chrom, starts_by_chrom, max_width
    ) == ["near_bin"]
    assert overlapping_bin_ids(
        "chr1", 1145, bins_by_chrom, starts_by_chrom, max_width
    ) == ["far_bin"]
    assert (
        overlapping_bin_ids("chr1", 1045, bins_by_chrom, starts_by_chrom, max_width)
        == []
    )
    assert (
        overlapping_bin_ids("chr9", 5, bins_by_chrom, starts_by_chrom, max_width) == []
    )


def test_get_nh_value_defaults_to_one_when_tag_missing():
    class FakeRecord:
        def get_tag(self, name):
            raise KeyError(name)

    assert get_nh_value(FakeRecord()) == 1


def test_get_nh_value_returns_tag_value():
    class FakeRecord:
        def get_tag(self, name):
            assert name == "NH"
            return 3

    assert get_nh_value(FakeRecord()) == 3


def test_write_counts_and_validate_output_round_trip(tmp_path):
    bed = tmp_path / "bins.bed"
    output = tmp_path / "counts.tsv"
    _write_bins_bed(bed, [("chr1", 1000, 1010, "bin1"), ("chr1", 1140, 1150, "bin2")])
    counts = {"bin1": 1.0, "bin2": 0.5}

    rows_written = write_counts(str(output), "sample_a", str(bed), counts)

    assert rows_written == 2
    lines = output.read_text().splitlines()
    assert lines[0].split("\t") == [
        "sample",
        "chrom",
        "start",
        "end",
        "bin_id",
        "gene_id",
        "gene_name",
        "gene_type",
        "strand",
        "tss",
        "tn5_site_count",
    ]
    assert lines[1].split("\t")[-1] == "1"  # round(1.0)
    assert lines[2].split("\t")[-1] == "0"  # round(0.5) -> banker's rounding to 0

    validate_output(str(output), rows_written, expected_rows=2)


def _write_perline_bins_bed(path, rows):
    """rows: list of (chrom, start, end, bin_id, gene_id, gene_name, gene_type,
    strand, reference_pos, feature_type, repeat_name, repeat_family, sw_score)
    -- the 13-column per-line/repeatMasker-mode schema build_tn5_midpoint_bins.py
    emits when called without --gene-types."""
    with open(path, "w", encoding="utf-8") as handle:
        for row in rows:
            handle.write("\t".join(str(field) for field in row) + "\n")


def test_write_counts_detects_perline_repeatmasker_schema(tmp_path):
    # Regression test for the schema-mismatch bug: write_counts() used to
    # assume every bins BED was 9-column grouped-mode, silently mislabeling
    # (and truncating) the 13-column per-line/repeatMasker-mode schema that
    # build_tn5_midpoint_bins.py emits by default (no --gene-types). This
    # exercises the real repeat-element path end-to-end through write_counts.
    bed = tmp_path / "bins.bed"
    output = tmp_path / "counts.tsv"
    _write_perline_bins_bed(
        bed,
        [
            (
                "chr1",
                3000,
                3100,
                "SINE_Alu_chr1_3000_3100_AluY",
                "SINE_Alu_3000_3100_AluY",
                "SINE_Alu_3000_3100_AluY",
                "SINE/Alu",
                "+",
                3049,
                "SINE/Alu",
                "AluY",
                "SINE/Alu",
                500,
            ),
        ],
    )
    counts = {"SINE_Alu_chr1_3000_3100_AluY": 2.0}

    rows_written = write_counts(str(output), "sample_a", str(bed), counts)

    assert rows_written == 1
    lines = output.read_text().splitlines()
    assert lines[0].split("\t") == (
        ["sample"]
        + BIN_CORE_COLUMNS
        + PERLINE_BIN_METADATA_COLUMNS
        + ["tn5_site_count"]
    )
    row = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
    # None of these should have shifted, and none should have been dropped --
    # the two bugs the original 9-column-only write_counts() had.
    assert row["gene_id"] == "SINE_Alu_3000_3100_AluY"
    assert row["gene_name"] == "SINE_Alu_3000_3100_AluY"
    assert row["gene_type"] == "SINE/Alu"
    assert row["strand"] == "+"
    assert row["reference_pos"] == "3049"
    assert row["feature_type"] == "SINE/Alu"
    assert row["repeat_name"] == "AluY"
    assert row["repeat_family"] == "SINE/Alu"
    assert row["sw_score"] == "500"
    assert row["tn5_site_count"] == "2"


def test_write_counts_empty_bins_bed_falls_back_to_grouped_schema(tmp_path):
    bed = tmp_path / "bins.bed"
    bed.write_text("")
    output = tmp_path / "counts.tsv"

    rows_written = write_counts(str(output), "sample_a", str(bed), {})

    assert rows_written == 0
    header = output.read_text().splitlines()[0].split("\t")
    assert header == (
        ["sample"]
        + BIN_CORE_COLUMNS
        + GROUPED_BIN_METADATA_COLUMNS
        + ["tn5_site_count"]
    )


def test_write_counts_raises_on_unrecognized_schema(tmp_path):
    bed = tmp_path / "bins.bed"
    bed.write_text("chr1\t100\t200\tbin1\tf1\tf2\tf3\tf4\tf5\tf6\n")  # 10 fields
    output = tmp_path / "counts.tsv"

    with pytest.raises(ValueError, match="Unrecognized bins-BED schema"):
        write_counts(str(output), "sample_a", str(bed), {"bin1": 1.0})


def test_validate_output_raises_on_row_count_mismatch(tmp_path):
    output = tmp_path / "counts.tsv"
    output.write_text("header\nrow1\nrow2\n")

    with pytest.raises(RuntimeError, match="row count mismatch"):
        validate_output(str(output), rows_written=2, expected_rows=3)


def _make_synthetic_bam(path, header_dict, records):
    import pysam

    with pysam.AlignmentFile(str(path), "wb", header=header_dict) as bam:
        for rec in records:
            bam.write(rec)
    pysam.index(str(path))


def test_end_to_end_counts_each_mates_cut_site_into_its_own_bin(tmp_path):
    pysam = pytest.importorskip("pysam")

    header_dict = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 10000}]}
    header = pysam.AlignmentHeader.from_dict(header_dict)

    def make_record(query_name, flag, reference_start, cigarstring):
        rec = pysam.AlignedSegment(header)
        rec.query_name = query_name
        rec.flag = flag
        rec.reference_id = 0
        rec.reference_start = reference_start
        rec.mapping_quality = 60
        rec.cigarstring = cigarstring
        rec.next_reference_id = 0
        return rec

    r1 = make_record("pair1", flag=99, reference_start=1000, cigarstring="50M")
    r2 = make_record("pair1", flag=147, reference_start=1100, cigarstring="50M")

    bam_path = tmp_path / "reads.bam"
    _make_synthetic_bam(bam_path, header_dict, [r1, r2])

    bed = tmp_path / "bins.bed"
    _write_bins_bed(
        bed, [("chr1", 1000, 1010, "near_bin"), ("chr1", 1140, 1150, "far_bin")]
    )
    output = tmp_path / "counts.tsv"

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--bam",
            str(bam_path),
            "--bins-bed",
            str(bed),
            "--output",
            str(output),
            "--sample",
            "sample_a",
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0, result.stderr
    rows = {
        line.split("\t")[4]: line.split("\t")[-1]
        for line in output.read_text().splitlines()[1:]
    }
    assert rows["near_bin"] == "1"
    assert rows["far_bin"] == "1"


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
                "--bins-bed",
                "bins.bed",
                "--output",
                "out.tsv",
                "--sample",
                "sample_a",
                "--mapq-min",
                "-1",
            ]
        )
    assert exc_info.value.code == 2
    assert "--mapq-min" in capsys.readouterr().err


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
