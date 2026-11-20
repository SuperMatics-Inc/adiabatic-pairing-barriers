from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
LANES = (
    ("LINEAR6", ROOT / "linear6_replay.py", r"(\d+)/(\d+) mutants rejected"),
    ("LINEAR7", ROOT / "linear7_replay.py", r"all (\d+) mutants rejected with exit 1"),
)


def run(script: Path, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-B", str(script), *args],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )


def main() -> int:
    total_rejected = 0
    total_mutants = 0
    for lane, script, count_pattern in LANES:
        baseline = run(script)
        print(f"{lane} baseline: exit {baseline.returncode}")
        if baseline.stdout.strip():
            print(baseline.stdout.strip())
        if baseline.stderr.strip():
            print(baseline.stderr.strip())
        if baseline.returncode != 0:
            return baseline.returncode or 1

        mutant_run = run(script, "--all-mutants")
        if mutant_run.stdout.strip():
            print(mutant_run.stdout.strip())
        if mutant_run.stderr.strip():
            print(mutant_run.stderr.strip())
        if mutant_run.returncode != 0:
            print(f"{lane} all-mutants: exit {mutant_run.returncode}")
            return mutant_run.returncode
        match = re.search(count_pattern, mutant_run.stdout)
        if not match:
            print(f"{lane}: could not find the mutant rejection count")
            return 1
        if lane == "LINEAR6":
            rejected, count = map(int, match.groups())
        else:
            rejected = count = int(match.group(1))
        if rejected != count:
            print(f"{lane}: only {rejected}/{count} mutants rejected")
            return 1
        total_rejected += rejected
        total_mutants += count
        print(f"{lane} all-mutants: exit 0; {rejected}/{count} rejected")

    print(f"COMBINED: {total_rejected}/{total_mutants} mutants rejected; all replay modes passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
