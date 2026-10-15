"""Run a job's commands in order: ``python -m ichnos_vail.chain <job folder>``.

Started by ``jobs.start_job`` in its own session. Each command's output goes
to log.txt; the first command that fails stops the job. Never uses a shell.
"""
from __future__ import annotations

import os
import subprocess
import sys
import time
from pathlib import Path

from .jobs import ROOT, _read_json, _write_json


def _wait_for_launcher(job: Path, timeout: float = 5.0) -> None:
    """Let start_job record 'running' first, so it can never overwrite this runner's final state."""
    end = time.time() + timeout
    while _read_json(job / "status.json").get("state") in (None, "created") and time.time() < end:
        time.sleep(0.02)


def run(job: Path) -> int:
    _wait_for_launcher(job)
    spec = _read_json(job / "job.json")
    commands = spec.get("commands", [])
    status = {"state": "running", "pid": os.getpid(), "step": 0, "steps": len(commands), "started": time.time()}
    _write_json(job / "status.json", status)
    env = {**os.environ, "PYTHONUNBUFFERED": "1"}
    with (job / "log.txt").open("a", encoding="utf-8") as log:
        for index, command in enumerate(commands):
            status.update(step=index + 1)
            _write_json(job / "status.json", status)
            log.write(f"\n$ {' '.join(command)}\n")
            log.flush()
            code = subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT,
                                  stdin=subprocess.DEVNULL, env=env).returncode
            if code != 0:
                stopping = _read_json(job / "status.json").get("state") == "stopping"
                status.update(state="stopped" if stopping else "failed", exit_code=code, ended=time.time())
                log.write(f"\n[{'stopped' if stopping else f'failed with exit code {code}'} at step {index + 1}]\n")
                _write_json(job / "status.json", status)
                return code or 1
    status.update(state="finished", exit_code=0, ended=time.time())
    _write_json(job / "status.json", status)
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: python -m ichnos_vail.chain <job folder>")
    raise SystemExit(run(Path(sys.argv[1])))
