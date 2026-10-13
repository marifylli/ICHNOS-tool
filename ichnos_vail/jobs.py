"""Background jobs: a list of commands run one after another, detached from ichnos_vail.

A job is a folder under the jobs directory:

  job.json     name, creation time, the commands (argument lists, never a shell string),
               the files the job is expected to write, and how to measure progress
  status.json  state (running / finished / failed / stopped / interrupted), current
               step, exit code, process id, start and end times
  log.txt      everything the commands printed

The commands run in a separate process (``python -m ichnos_vail.chain <job>``)
in its own session, so closing the browser or restarting ichnos_vail does not
stop an 11-hour segmentation run; ichnos_vail only reads these files.
"""
from __future__ import annotations

import json
import os
import re
import signal
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_JOBS_DIR = ROOT / "outputs" / "ichnos_vail_jobs"


def default_jobs_dir() -> Path:
    """Where jobs are kept: $ICHNOS_VAIL_JOBS_DIR if set, else outputs/ichnos_vail_jobs in the repository."""
    return Path(os.environ.get("ICHNOS_VAIL_JOBS_DIR", DEFAULT_JOBS_DIR))
FINAL_STATES = {"finished", "failed", "stopped", "interrupted"}


def _slug(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "-", text).strip("-").lower()[:40] or "job"


def _write_json(path: Path, data: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, ensure_ascii=False))
    os.replace(tmp, path)  # atomic: ichnos_vail never reads a half-written file


def _read_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return {}


def create_job(name: str, commands: list[list[str]], jobs_dir: Path | None = None,
               outputs: list[str] | None = None, progress: dict | None = None) -> Path:
    """Write a job folder; nothing runs until ``start_job``."""
    jobs_dir = Path(jobs_dir) if jobs_dir is not None else default_jobs_dir()
    if not commands or not all(isinstance(c, list) and c and all(isinstance(a, str) for a in c) for c in commands):
        raise ValueError("commands must be a non-empty list of argument lists of strings")
    jobs_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    folder = jobs_dir / f"{stamp}_{_slug(name)}"
    n = 1
    while folder.exists():
        n += 1
        folder = jobs_dir / f"{stamp}_{_slug(name)}-{n}"
    folder.mkdir()
    _write_json(folder / "job.json", {"name": name, "created": datetime.now().isoformat(timespec="seconds"),
                                      "commands": commands, "outputs": outputs or [], "progress": progress or {}})
    _write_json(folder / "status.json", {"state": "created"})
    return folder


def start_job(job: Path) -> int:
    """Launch the job's runner in its own session; return its process id."""
    process = subprocess.Popen([sys.executable, "-m", "ichnos_vail.chain", str(job)], cwd=ROOT,
                               stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                               start_new_session=True)
    status = _read_json(job / "status.json")
    if status.get("state") in (None, "created"):
        _write_json(job / "status.json", {"state": "running", "pid": process.pid, "step": 0,
                                          "started": time.time()})
    return process.pid


def _alive(pid: int | None) -> bool:
    if not pid:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    # a finished child of this process lingers as a zombie until reaped
    try:
        done, _ = os.waitpid(pid, os.WNOHANG)
        return done == 0
    except ChildProcessError:
        return True


def read_status(job: Path) -> dict:
    """The job's status; a 'running' job whose process has vanished is reported 'interrupted'."""
    status = _read_json(job / "status.json")
    if status.get("state") in ("running", "stopping") and not _alive(status.get("pid")):
        time.sleep(0.2)  # the runner may be writing its final status right now
        status = _read_json(job / "status.json")
        if status.get("state") == "stopping":  # stopped from ichnos_vail; the runner died with its children
            status = {**status, "state": "stopped"}
        elif status.get("state") == "running":  # killed from outside, machine slept or rebooted
            status = {**status, "state": "interrupted"}
    return status


def read_job(job: Path) -> dict:
    return _read_json(job / "job.json")


def tail(job: Path, lines: int = 60) -> str:
    log = job / "log.txt"
    if not log.exists():
        return ""
    with log.open("rb") as handle:
        handle.seek(0, os.SEEK_END)
        size = handle.tell()
        handle.seek(max(0, size - 64_000))
        text = handle.read().decode("utf-8", errors="replace")
    return "\n".join(text.splitlines()[-lines:])


def progress(job: Path) -> float | None:
    """Fraction done for jobs that declare how to count it (batches of the segmentation run)."""
    spec = read_job(job).get("progress", {})
    if spec.get("kind") == "batches" and spec.get("total"):
        log = job / "log.txt"
        text = log.read_text(errors="replace") if log.exists() else ""
        started = len(re.findall(r"^batch \d+: manifest rows", text, flags=re.M))
        if read_status(job).get("state") == "finished":
            return 1.0
        return min(max(started - 1, 0) / spec["total"], 1.0)
    return None


def stop_job(job: Path) -> None:
    """Stop the runner and everything it started (they share its process group)."""
    status = read_status(job)
    pid = status.get("pid")
    if status.get("state") != "running" or not pid:
        return
    _write_json(job / "status.json", {**status, "state": "stopping"})
    try:
        os.killpg(os.getpgid(pid), signal.SIGTERM)
    except (ProcessLookupError, PermissionError):
        pass


def list_jobs(jobs_dir: Path | None = None, limit: int = 30) -> list[dict]:
    jobs_dir = Path(jobs_dir) if jobs_dir is not None else default_jobs_dir()
    if not jobs_dir.exists():
        return []
    rows = []
    for folder in sorted((p for p in jobs_dir.iterdir() if (p / "job.json").exists()), reverse=True)[:limit]:
        job, status = read_job(folder), read_status(folder)
        rows.append({"path": folder, "name": job.get("name", folder.name), "created": job.get("created", ""),
                     "state": status.get("state", "?"), "step": status.get("step"),
                     "steps": len(job.get("commands", []))})
    return rows
