"""Unit tests for the chroma2 wrapper's logging + pipeline.* state-marker helpers."""

import importlib.machinery
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

WRAPPER_PATH = Path(__file__).parent.parent.parent / "chroma2"

_loader = importlib.machinery.SourceFileLoader("chroma2_wrapper", str(WRAPPER_PATH))
_spec = importlib.util.spec_from_loader("chroma2_wrapper", _loader)
chroma2_wrapper = importlib.util.module_from_spec(_spec)
_loader.exec_module(chroma2_wrapper)


# ---------------------------------------------------------------------------
# logging helpers
# ---------------------------------------------------------------------------


def test_log_helpers_use_fixed_width_tags(capsys):
    chroma2_wrapper.log_info("hello")
    chroma2_wrapper.log_step("hello")
    chroma2_wrapper.log_ok("hello")
    chroma2_wrapper.log_warn("hello")
    chroma2_wrapper.log_error("hello")
    chroma2_wrapper.log_next("hello")
    chroma2_wrapper.log_divider()

    lines = capsys.readouterr().out.splitlines()
    assert lines == [
        "INFO  hello",
        "STEP  hello",
        "OK    hello",
        "WARN  hello",
        "ERROR hello",
        "NEXT  hello",
        "-" * 66,
    ]


# ---------------------------------------------------------------------------
# pipeline.* state markers
# ---------------------------------------------------------------------------


def test_write_pipeline_state_marker_writes_marker_and_sidecar(tmp_path):
    (tmp_path / "pipeline.running").write_text("stale progress")

    chroma2_wrapper.write_pipeline_state_marker(
        tmp_path,
        "completed",
        reason="test_reason",
        job_id="123",
        content="all done",
        runmode="run",
    )

    assert (tmp_path / "pipeline.completed").read_text() == "all done"
    assert not (tmp_path / "pipeline.running").exists()
    assert not (tmp_path / "pipeline.failed").exists()
    assert not (tmp_path / "pipeline.canceled").exists()

    sidecar = json.loads((tmp_path / "pipeline.status.json").read_text())
    assert sidecar["pipeline"] == "chroma2"
    assert sidecar["state"] == "completed"
    assert sidecar["reason"] == "test_reason"
    assert sidecar["slurm_job_id"] == "123"
    assert sidecar["runmode"] == "run"
    assert sidecar["num_failed_rules"] == 0
    assert sidecar["failed_rules"] == []
    assert "timestamp_utc" in sidecar


def test_write_pipeline_state_marker_removes_siblings_but_creates_new_first(tmp_path):
    chroma2_wrapper.write_pipeline_state_marker(tmp_path, "running", reason="start")
    assert (tmp_path / "pipeline.running").exists()

    chroma2_wrapper.write_pipeline_state_marker(
        tmp_path, "failed", reason="boom", num_failed_rules=2
    )
    assert not (tmp_path / "pipeline.running").exists()
    assert (tmp_path / "pipeline.failed").exists()


def test_check_not_already_running_raises_when_marker_present(tmp_path):
    (tmp_path / "pipeline.running").touch()
    try:
        chroma2_wrapper.check_not_already_running(tmp_path)
        assert False, "expected SystemExit"
    except SystemExit as e:
        assert e.code == 1


def test_check_not_already_running_noop_when_no_marker(tmp_path):
    chroma2_wrapper.check_not_already_running(tmp_path)  # should not raise


def test_warn_if_already_running_prints_warn_without_raising(tmp_path, capsys):
    (tmp_path / "pipeline.running").touch()
    chroma2_wrapper.warn_if_already_running(tmp_path)
    out = capsys.readouterr().out
    assert out.startswith("WARN  ")

    capsys.readouterr()
    chroma2_wrapper.warn_if_already_running(tmp_path.parent / "does-not-exist")
    assert capsys.readouterr().out == ""


# ---------------------------------------------------------------------------
# failed-rule digest
# ---------------------------------------------------------------------------

FIXTURE_LOG = """\
Some unrelated Snakemake output
Error in rule align:
    jobid: 5
    message: SLURM-job '999' failed, SLURM status is: 'FAILED'.
    log: logs/align/sampleA/999.log (check log file(s) for error details)

Error in rule peakcall:
    jobid: 7
    message: SLURM-job '1000' failed, SLURM status is: 'FAILED'.
    log: logs/peakcall/sampleB/1000.log (check log file(s) for error details)

Finished jobid: 7
More trailing output
"""


def test_extract_failed_rules_digest_drops_jobs_that_succeeded_on_retry(tmp_path):
    logfile = tmp_path / "snakemake.log"
    logfile.write_text(FIXTURE_LOG)

    digest = chroma2_wrapper.extract_failed_rules_digest(logfile)

    assert len(digest) == 1
    entry = digest[0]
    assert entry["rule"] == "align"
    assert entry["jobid"] == "5"
    assert entry["slurm_job"] == "999"
    assert entry["state"] == "FAILED"
    assert entry["sample"] == "sampleA"
    assert entry["log"] == "logs/align/sampleA/999.log"


def test_extract_failed_rules_digest_missing_file_returns_empty(tmp_path):
    assert chroma2_wrapper.extract_failed_rules_digest(tmp_path / "nope.log") == []


def test_format_digest_line():
    line = chroma2_wrapper.format_digest_line(
        {
            "rule": "align",
            "jobid": "5",
            "slurm_job": "999",
            "state": "FAILED",
            "sample": "sampleA",
            "log": "x.log",
        }
    )
    assert (
        line
        == "FAILED: rule=align jobid=5 slurm_job=999 state=FAILED sample=sampleA log=x.log"
    )


# ---------------------------------------------------------------------------
# progress snapshots
# ---------------------------------------------------------------------------


def test_progress_snapshot_and_derive_final_progress(tmp_path):
    logfile = tmp_path / "run.log"
    logfile.write_text("5 of 10 steps (50%) done\n8 of 10 steps (80%) done\n")

    snapshot = chroma2_wrapper.progress_snapshot(logfile)
    assert "Progress : 8 / 10 steps complete (80%)" in snapshot
    assert "Remaining: 2 steps" in snapshot

    final = chroma2_wrapper.derive_final_progress(logfile)
    assert "Progress : 10 / 10 steps complete (100%)" in final
    assert "Remaining: 0 steps" in final


def test_progress_snapshot_no_match_returns_none(tmp_path):
    logfile = tmp_path / "run.log"
    logfile.write_text("nothing useful here\n")
    assert chroma2_wrapper.progress_snapshot(logfile) is None


# ---------------------------------------------------------------------------
# --internal-finalize CLI mode
# ---------------------------------------------------------------------------


def test_internal_finalize_completed(tmp_path):
    logfile = tmp_path / "slurm-1.out"
    logfile.write_text("10 of 10 steps (100%) done\n")

    result = subprocess.run(
        [
            sys.executable,
            str(WRAPPER_PATH),
            "--internal-finalize",
            "--workdir",
            str(tmp_path),
            "--state",
            "completed",
            "--reason",
            "snakemake_succeeded",
            "--job-id",
            "42",
            "--log-file",
            str(logfile),
        ],
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "pipeline.completed").exists()
    sidecar = json.loads((tmp_path / "pipeline.status.json").read_text())
    assert sidecar["state"] == "completed"
    assert sidecar["slurm_job_id"] == "42"


def test_internal_finalize_failed_includes_digest(tmp_path):
    logfile = tmp_path / "slurm-2.out"
    logfile.write_text(FIXTURE_LOG)
    (tmp_path / "pipeline.running").write_text(
        "Progress : 3 / 10 steps complete (30%)\n"
    )

    result = subprocess.run(
        [
            sys.executable,
            str(WRAPPER_PATH),
            "--internal-finalize",
            "--workdir",
            str(tmp_path),
            "--state",
            "failed",
            "--reason",
            "snakemake_failed",
            "--job-id",
            "43",
            "--log-file",
            str(logfile),
        ],
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, result.stderr
    sidecar = json.loads((tmp_path / "pipeline.status.json").read_text())
    assert sidecar["state"] == "failed"
    assert sidecar["num_failed_rules"] == 1
    assert "align" in sidecar["failed_rules"][0]
    content = (tmp_path / "pipeline.failed").read_text()
    assert "3 / 10 steps complete" in content
    assert "1 job(s) failed" in content


def test_internal_finalize_canceled_no_logfile(tmp_path):
    result = subprocess.run(
        [
            sys.executable,
            str(WRAPPER_PATH),
            "--internal-finalize",
            "--workdir",
            str(tmp_path),
            "--state",
            "canceled",
            "--reason",
            "signal_SIGTERM",
            "--job-id",
            "44",
        ],
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "pipeline.canceled").exists()
