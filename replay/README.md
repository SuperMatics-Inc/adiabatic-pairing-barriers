# Replay bundle

`linear6_replay.py` and `linear7_replay.py` are byte-for-byte copies of the replay scripts named in `SOURCE-PINS.json`. The pin file records each source path and SHA-256. `run_all.py` runs each lane's baseline checks and its complete mutant suite, then reports the combined rejection count. The lane scripts are not modified by the driver.

From this directory, run:

```powershell
$env:UV_OFFLINE = "1"
$env:UV_CACHE_DIR = "D:/dev/.cache/uv"
uv run --no-project --with numpy --with scipy --with mpmath --with python-flint --with pymupdf python -B run_all.py
```

A successful run exits zero only if both baseline replays and both all-mutant runs exit zero and every listed mutant is rejected. Individual mutant probes intentionally exit one when rejected; each lane's all-mutants mode interprets that result and exits zero on a complete rejection set.
