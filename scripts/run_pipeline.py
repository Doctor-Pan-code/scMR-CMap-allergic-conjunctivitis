"""Sequential runner for Days 3-7 plus the two report builders.

Waits for the Day 2 process to exit, then runs each downstream day in order.
A failing day is recorded and the chain continues so that the log shows every
independent failure in one pass instead of one at a time.
"""
from __future__ import annotations

import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = ROOT / "scripts"
LOGS = ROOT / "logs"
PY = sys.executable

STEPS = [
    ("day3", "day3_coloc_smr.py"),
    ("day4", "day4_scrna_signature.py"),
    ("day5", "day5_cmap_reversal.py"),
    ("day6", "day6_phewas_safety.py"),
    ("day7", "day7_integration.py"),
    ("report_plan", "build_report_plan.py"),
    ("report_master", "build_report_master.py"),
]


def stamp() -> str:
    return datetime.now().strftime("%H:%M:%S")


def say(msg: str) -> None:
    line = f"[{stamp()}] {msg}"
    print(line, flush=True)
    with open(LOGS / "pipeline.log", "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def wait_for(pid: int, timeout_s: int = 7200) -> bool:
    """Block until the given PID is gone. Returns True if it exited in time."""
    t0 = time.time()
    while time.time() - t0 < timeout_s:
        r = subprocess.run(["tasklist", "/FI", f"PID eq {pid}", "/NH"],
                           capture_output=True, text=True)
        if str(pid) not in (r.stdout or ""):
            return True
        time.sleep(10)
    return False


def run(name: str, script: str) -> int:
    say(f"=== {name}: {script} ===")
    t0 = time.time()
    log_path = LOGS / f"pipeline_{name}.log"
    with open(log_path, "w", encoding="utf-8") as fh:
        proc = subprocess.run([PY, "-X", "utf8", str(SCRIPTS / script)],
                              cwd=str(ROOT), stdout=fh, stderr=subprocess.STDOUT)
    dt = time.time() - t0
    say(f"=== {name}: exit={proc.returncode} in {dt/60:.1f} min (log: {log_path.name}) ===")
    return proc.returncode


def main() -> None:
    pid = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    if pid:
        say(f"waiting for Day 2 process {pid} to exit")
        ok = wait_for(pid)
        say(f"Day 2 process gone: {ok}")
    codes = {}
    for name, script in STEPS:
        codes[name] = run(name, script)
    say("pipeline summary: " + ", ".join(f"{k}={v}" for k, v in codes.items()))
    bad = [k for k, v in codes.items() if v != 0]
    say("all steps succeeded" if not bad else f"failed steps: {bad}")


if __name__ == "__main__":
    main()
