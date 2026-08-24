"""Unit tests for build_tn5_count_matrix.py"""

import subprocess
import sys
from pathlib import Path

SCRIPT_PATH = Path(__file__).parent.parent / "scripts" / "build_tn5_count_matrix.py"


HEADER = (
    "sample\tchrom\tstart\tend\tbin_id\tgene_id\tgene_name\tgene_type\tstrand\t"
    "tss\ttn5_site_count\n"
)


def _write_count_tsv(path, sample, counts):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(HEADER)
        for idx, count in enumerate(counts, start=1):
            handle.write(
                "\t".join(
                    [
                        sample,
                        "chr1",
                        str(idx * 100),
                        str(idx * 100 + 201),
                        f"bin{idx}",
                        f"gene{idx}",
                        f"gene{idx}",
                        "tRNA",
                        "+",
                        str(idx * 100 + 100),
                        str(count),
                    ]
                )
                + "\n"
            )


def _run_matrix(counts, output):
    return subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--counts",
            *[str(path) for path in counts],
            "--output",
            str(output),
        ],
        text=True,
        capture_output=True,
    )


def test_matrix_build_succeeds_with_valid_inputs(tmp_path):
    count_a = tmp_path / "sample_a.tsv"
    count_b = tmp_path / "sample_b.tsv"
    output = tmp_path / "matrix.tsv"
    _write_count_tsv(count_a, "sample_a", [3, 5])
    _write_count_tsv(count_b, "sample_b", [7, 11])

    result = _run_matrix([count_a, count_b], output)

    assert result.returncode == 0, result.stderr
    assert "Created Tn5 count matrix" in result.stderr
    assert output.stat().st_size > 0
    lines = output.read_text().splitlines()
    assert lines[0].endswith("\tsample_a\tsample_b")
    assert lines[1].endswith("\t3\t7")
    assert lines[2].endswith("\t5\t11")


def test_matrix_build_reports_missing_inputs(tmp_path):
    count_a = tmp_path / "sample_a.tsv"
    missing = tmp_path / "missing.tsv"
    output = tmp_path / "matrix.tsv"
    _write_count_tsv(count_a, "sample_a", [3])

    result = _run_matrix([count_a, missing], output)

    assert result.returncode != 0
    assert "Invalid Tn5 count matrix input TSVs" in result.stderr
    assert "missing:" in result.stderr
    assert str(missing) in result.stderr
    assert not output.exists()


def test_matrix_build_reports_empty_inputs(tmp_path):
    count_a = tmp_path / "sample_a.tsv"
    empty = tmp_path / "empty.tsv"
    output = tmp_path / "matrix.tsv"
    _write_count_tsv(count_a, "sample_a", [3])
    empty.touch()

    result = _run_matrix([count_a, empty], output)

    assert result.returncode != 0
    assert "Invalid Tn5 count matrix input TSVs" in result.stderr
    assert "empty:" in result.stderr
    assert str(empty) in result.stderr
    assert not output.exists()


def test_matrix_build_skips_header_only_inputs(tmp_path):
    count_a = tmp_path / "sample_a.tsv"
    header_only = tmp_path / "header_only.tsv"
    output = tmp_path / "matrix.tsv"
    _write_count_tsv(count_a, "sample_a", [3])
    header_only.write_text(HEADER)

    result = _run_matrix([count_a, header_only], output)

    assert result.returncode == 0, result.stderr
    assert "Created Tn5 count matrix" in result.stderr
    lines = output.read_text().splitlines()
    assert lines[0].endswith("\tsample_a")
    assert lines[1].endswith("\t3")


def test_matrix_build_writes_header_only_output_when_all_inputs_are_header_only(
    tmp_path,
):
    header_only_a = tmp_path / "header_only_a.tsv"
    header_only_b = tmp_path / "header_only_b.tsv"
    output = tmp_path / "matrix.tsv"
    header_only_a.write_text(HEADER)
    header_only_b.write_text(HEADER)

    result = _run_matrix([header_only_a, header_only_b], output)

    assert result.returncode == 0, result.stderr
    lines = output.read_text().splitlines()
    assert len(lines) == 1
    assert lines[0].split("\t") == [
        "chrom",
        "start",
        "end",
        "bin_id",
        "gene_id",
        "gene_name",
        "gene_type",
        "strand",
        "tss",
    ]
