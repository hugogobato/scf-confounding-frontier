"""SCF Phase 3 pass-1 driver: run calibration arms (null, perm_null), derive
frozen thresholds and predicted frontiers g*, write
results/benchmark_freeze.json BEFORE any pass-2 arm is generated.

Thresholds frozen per benchmark config:
  S2_mc95        : q95 of t_maxz over the matched-null reps (gate statistic)
  ucm_q95 / js_q95 : q95 over permuted-Y reps (permutation calibration of the
                     two APPROXIMATE baselines)
  g_star         : smallest g with F12-law-predicted alarm power >= 0.8 at
                   size calibrated to S2_mc95 (predicted_g_star)

Usage: python3 code/bench_freeze.py [--reps-scale 1.0] [--workers 3]
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


def main():
    args = sys.argv[1:]
    reps_scale = 1.0
    workers = 3
    budget = None
    if "--reps-scale" in args:
        reps_scale = float(args[args.index("--reps-scale") + 1])
    if "--workers" in args:
        workers = int(args[args.index("--workers") + 1])
    if "--time-budget" in args:
        budget = float(args[args.index("--time-budget") + 1])

    import time as _time
    jobs = B.build_jobs("pass1", None)
    chunk = None
    if "--reps-chunk" in args:
        chunk = int(args[args.index("--reps-chunk") + 1])
    for j in jobs:
        j["spec"] = {**j["spec"],
                     "reps": max(50, int(j["spec"]["reps"] * reps_scale))}
        if chunk is not None:
            import pandas as _pd
            rp = B.RAW_DIR / j["config_name"] / f"{j['arm']}.parquet"
            have = 0
            if rp.exists():
                try:
                    have = int(_pd.read_parquet(rp,
                                                columns=["rep"])["rep"].nunique())
                except Exception:
                    have = 0
            j["spec"]["reps"] = min(j["spec"]["reps"], have + chunk)
    print(f"pass1: {len(jobs)} cells", flush=True)
    t_start = _time.perf_counter()
    incomplete = []
    if budget is not None or workers <= 1:
        for j in jobs:
            if budget is not None and \
                    _time.perf_counter() - t_start > budget:
                incomplete.append(j)
                continue
            B.run_cell(j)
    else:
        B.run_jobs(jobs, workers=workers)
    remaining = [j for j in jobs if not _cell_done(j)]
    if remaining:
        print(f"[budget] {len(remaining)} pass-1 cells still incomplete: "
              + ", ".join(f"{j['config_name']}/{j['arm']}" for j in remaining),
              flush=True)
        return
    freeze_configs(cfg := B.arm_specs(), freeze := {"version": cfg["version"],
                                                    "ledger_hash": cfg["ledger_hash"],
                                                    "configs": {}})
    (ROOT / "results").mkdir(exist_ok=True)
    out = ROOT / "results" / "benchmark_freeze.json"
    out.write_text(json.dumps(freeze, indent=1))
    print("wrote", out)


def _cell_done(job) -> bool:
    import pandas as pd

    raw_path = B.RAW_DIR / job["config_name"] / f"{job['arm']}.parquet"
    if not raw_path.exists():
        return False
    try:
        have = pd.read_parquet(raw_path, columns=["rep"])
        return int(have["rep"].nunique()) >= job["spec"]["reps"]
    except Exception:
        return False


def freeze_configs(cfg, freeze):
    for name in cfg["configs"]:
        null_df = pd.read_parquet(B.RAW_DIR / name / "null.parquet")
        zcols = [c_ for c_ in null_df.columns if c_.startswith("zeta")]
        zetas = null_df[zcols].to_numpy()
        coord_scales = np.sqrt((zetas ** 2).mean(axis=0))
        T_null = (np.max(np.abs(zetas) / coord_scales[None, :], axis=1)
                  if len(zcols) else null_df["T_bench"].to_numpy())
        mc95 = float(np.quantile(T_null, 0.95))
        s2_size_bonf = float(null_df["rej_s2_bonf"].mean())
        s0_size = float(null_df["rej_s0_tw99"].mean())
        b1_size = float(null_df["rej_b1_f95"].mean())
        perm_df = pd.read_parquet(B.RAW_DIR / name / "perm_null.parquet")
        bench = B.Bench(name)
        entry = {
            "mc95_S2": round(mc95, 4),
            "coord_scales": [round(float(v), 4) for v in coord_scales],
            "size_s2_analytic_null": round(s2_size_bonf, 4),
            "size_s0_null": round(s0_size, 4),
            "size_b1_null": round(b1_size, 4),
            "ucm_q95_perm": round(float(np.quantile(
                perm_df["ucm_rho"].to_numpy(), 0.95)), 6),
            "js_q95_perm": round(float(np.quantile(
                perm_df["js_asym"].to_numpy(), 0.95)), 6),
            "se2_bench": round(bench.se2, 6),
            "ktop_alarm": int(bench.ktop),
            "r_inj": int(bench.r_inj),
            "l_hat_used": [round(float(v), 3)
                           for v in bench.l_hat[:bench.r_inj]],
            "g_star": None,
        }
        entry["g_star"] = B.predicted_g_star(bench, mc95,
                                             coord_scales=coord_scales)
        freeze["configs"][name] = entry
        print(name, json.dumps(entry))


if __name__ == "__main__":
    main()
