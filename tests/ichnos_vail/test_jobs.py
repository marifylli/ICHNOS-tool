"""Background jobs: commands run in order, detached, with status and log on disk."""
import json
import os
import sys
import time

import pytest

from ichnos_vail import jobs


def wait(job, states=jobs.FINAL_STATES, timeout=30):
    end = time.time() + timeout
    while time.time() < end:
        status = jobs.read_status(job)
        if status.get("state") in states:
            return status
        time.sleep(0.1)
    raise AssertionError(f"job still {jobs.read_status(job)} after {timeout} s")


def py(code, *args):
    return [sys.executable, "-c", code, *args]


def test_commands_run_in_order_and_the_log_keeps_their_output(tmp_path):
    job = jobs.create_job("two steps", [py("print('first')"), py("print('second')")], jobs_dir=tmp_path)
    assert jobs.read_status(job)["state"] == "created"
    jobs.start_job(job)
    status = wait(job)
    assert status["state"] == "finished" and status["exit_code"] == 0 and status["step"] == 2
    log = jobs.tail(job)
    assert log.index("first") < log.index("second")


def test_the_first_failing_command_stops_the_job(tmp_path):
    marker = tmp_path / "third-ran"
    job = jobs.create_job("fails", [py("pass"), py("import sys; print('boom'); sys.exit(3)"),
                                    py("import pathlib, sys; pathlib.Path(sys.argv[1]).touch()", str(marker))],
                          jobs_dir=tmp_path / "jobs")
    jobs.start_job(job)
    status = wait(job)
    assert (status["state"], status["exit_code"], status["step"]) == ("failed", 3, 2)
    assert "boom" in jobs.tail(job) and not marker.exists()


def test_arguments_are_passed_literally_never_through_a_shell(tmp_path):
    out = tmp_path / "argv.json"
    tricky = "a b; rm -rf / $(echo x) 'q\""
    job = jobs.create_job("argv", [py("import json, sys; open(sys.argv[1], 'w').write(json.dumps(sys.argv[2:]))",
                                      str(out), tricky)], jobs_dir=tmp_path / "jobs")
    jobs.start_job(job)
    wait(job)
    assert json.loads(out.read_text()) == [tricky]


def test_a_job_can_be_stopped_with_everything_it_started(tmp_path):
    job = jobs.create_job("long", [py("import time; print('started', flush=True); time.sleep(60)")],
                          jobs_dir=tmp_path)
    jobs.start_job(job)
    end = time.time() + 15
    while "started" not in jobs.tail(job) and time.time() < end:
        time.sleep(0.1)
    jobs.stop_job(job)
    status = wait(job, timeout=15)
    assert status["state"] == "stopped"


def test_a_running_job_whose_process_vanished_is_reported_interrupted(tmp_path):
    job = jobs.create_job("ghost", [py("pass")], jobs_dir=tmp_path)
    (job / "status.json").write_text(json.dumps({"state": "running", "pid": 2**22 + 12345}))
    assert jobs.read_status(job)["state"] == "interrupted"


def test_progress_counts_finished_batches_from_the_log(tmp_path):
    job = jobs.create_job("cells", [py("pass")], jobs_dir=tmp_path, progress={"kind": "batches", "total": 4})
    (job / "log.txt").write_text("batch 0: manifest rows 0-19\nbatch 1: manifest rows 20-39\n")
    assert jobs.progress(job) == pytest.approx(0.25)
    assert jobs.progress(jobs.create_job("other", [py("pass")], jobs_dir=tmp_path)) is None


def test_invalid_commands_are_refused_and_jobs_are_listed_newest_first(tmp_path):
    with pytest.raises(ValueError):
        jobs.create_job("bad", ["python -c pass"], jobs_dir=tmp_path)
    first = jobs.create_job("one", [py("pass")], jobs_dir=tmp_path)
    second = jobs.create_job("two", [py("pass")], jobs_dir=tmp_path)
    listed = [row["path"] for row in jobs.list_jobs(tmp_path)]
    assert listed[:2] == sorted([first, second], reverse=True)


def test_the_jobs_folder_follows_the_environment(tmp_path, monkeypatch):
    monkeypatch.setenv("ICHNOS_VAIL_JOBS_DIR", str(tmp_path / "here"))
    job = jobs.create_job("env", [py("pass")])
    assert job.parent == tmp_path / "here" and os.path.isdir(job)
