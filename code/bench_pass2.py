"""SCF Phase 3 pass-2 driver: run evaluation arms against the FROZEN
thresholds/g* in results/benchmark_freeze.json. Resume-safe and chunked for
contended machines; run repeatedly until it reports all cells complete.

Usage: python3 code/bench_pass2.py [--workers 1] [--time-budget 280]
                                   [--reps-chunk 60]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import benchmarks as B  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def _cell_done(job) -> bool:
    raw_path = B.RAW_DIR / job["config_name"] / f"{job['arm']}.parquet"
    if not raw_path.exists():
        return False
    try:
        have = pd.read_parquet(raw_path, columns=["rep"])
        return int(have["rep"].nunique()) >= job["spec"]["reps"]
    except Exception:
        return False


def main():
    args = sys.argv[1:]
    workers = int(args[args.index("--workers") + 1]) if "--workers" in args else 1
    budget = float(args[args.index("--time-budget") + 1]) \
        if "--time-budget" in args else None
    chunk = int(args[args.index("--reps-chunk") + 1]) \
        if "--reps-chunk" in args else None

    freeze = json.loads((ROOT / "results" / "benchmark_freeze.json").read_text())
    jobs = B.build_jobs("pass2", freeze)
    import os as _os
    flt = _os.environ.get("SCF_CONFIGS")
    if flt:
        keep = set(flt.split(","))
        jobs = [j for j in jobs if j["config_name"] in keep]
    flt_a = _os.environ.get("SCF_ARMS")
    if flt_a:
        keepa = set(flt_a.split(","))
        jobs = [j for j in jobs if j["arm"] in keepa]
    import time as _time
    if chunk is not None:
        for j in jobs:
            rp = B.RAW_DIR / j["config_name"] / f"{j['arm']}.parquet"
            have = 0
            if rp.exists():
                try:
                    have = int(pd.read_parquet(rp,
                                               columns=["rep"])["rep"].nunique())
                except Exception:
                    have = 0
            j["spec"] = {**j["spec"],
                         "reps": min(j["spec"]["reps"], have + chunk)}
    print(f"pass2: {len(jobs)} cells", flush=True)
    t0 = _time.perf_counter()
    if budget is not None or workers <= 1:
        for j in jobs:
            if budget is not None and _time.perf_counter() - t0 > budget:
                break
            B.run_cell(j)
    else:
        B.run_jobs(jobs, workers=workers)
    remaining = [f"{j['config_name']}/{j['arm']}" for j in jobs
                 if not _cell_done(j)]
    if remaining:
        print(f"[incomplete] {len(remaining)}: " + ", ".join(remaining),
              flush=True)
    else:
        print("[complete] all pass-2 cells done", flush=True)


if __name__ == "__main__":
    main()
