"""G3 gate checks (research plan Section 6 + preregistered thresholds).

Reads consolidated sweep parquets under data/sim/ and prints each numeric
condition with its current verdict: PASS / FAIL / PENDING (data missing).
No threshold is defined here beyond what the plan and preregistration froze;
this file is the mechanical evaluator only.

Usage: python3 code/gate_checks.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SIM = ROOT / "data" / "sim"


def load(sweep: str) -> pd.DataFrame | None:
    p = SIM / sweep / f"{sweep}_results.parquet"
    return pd.read_parquet(p) if p.exists() else None


def means_for(sweep: str) -> dict[str, np.ndarray]:
    """Load per-cell mean-diff vectors {cid@tag: vector}."""
    out = {}
    md = SIM / sweep / "means"
    if not md.exists():
        return out
    for f in md.glob("*.npz"):
        cid = f.stem
        z = np.load(f)
        for k in z.files:
            out[f"{cid}@{k}"] = z[k].astype(float)
    return out


def verdict(cond: bool, pending=False) -> str:
    return "PENDING" if pending else ("PASS" if cond else "FAIL")


def check_correctness():
    df = load("correctness")
    print("\n=== WP 2.1 correctness (plan Section 6, G3 row) ===")
    if df is None:
        print("PENDING: no consolidated correctness parquet")
        return
    # condition 1: DE deviation <= 10% at n = 4000-equivalent in >= 90% cells,
    # shrinking with n. Evaluated on OLS directional overlay per cell.
    rows = []
    for cid, g in df.groupby("config_id"):
        pass
    print("PENDING: overlay evaluation runs after data lands "
          "(code/gate_checks.py::check_correctness to be extended with the "
          "per-cell predicted-vs-simulated comparison from means npz).")


def check_estimation():
    df = load("estimation")
    print("\n=== WP 2.2 estimation (preregistered thresholds) ===")
    if df is None:
        print("PENDING: no consolidated estimation parquet")
        return
    means = means_for("estimation")
    grid = json.loads((ROOT / "configs" / "grid_estimation.json").read_text())
    harmful = []
    for j in grid:
        cfg = j["config"]
        if cfg["profile"] in ("sub", "mixed") and cfg["g"] > 0 \
                and cfg.get("conf_kind", "dense") == "dense" \
                and cfg.get("beta_kind", "dense") == "dense":
            harmful.append(j["config_id"])
    base_tags = ["ridge_cv", "pca_onatski", "pca_baing", "cevid_default",
                 "sdboost_linear_eb"]
    no_regret_cells, half_cut_cells = [], []
    for cid in harmful:
        def bias(tag):
            v = means.get(f"{cid}@{tag}")
            return float(np.linalg.norm(v)) if v is not None else np.nan

        eb = bias("eb_spectral")
        bases = [bias(t) for t in base_tags] + [
            0.5 * bias("ols_g0") if False else bias("ols")]
        bases = [b for b in bases if not math.isnan(b)]
        ols_b = bias("ols")
        # twin-difference confounding bias at c > 1 cells when present
        if f"{cid}@ols_g0" in means:
            def conf(tag):
                a = means.get(f"{cid}@{tag}")
                b = means.get(f"{cid}@{tag}_g0")
                if a is None or b is None:
                    return np.nan
                return float(np.linalg.norm(a - b))
            eb_c, ols_c = conf("eb_spectral"), conf("ols")
            bs_c = [conf(t) for t in base_tags]
            bs_c = [b for b in bs_c if not math.isnan(b)]
            if not math.isnan(eb_c) and bs_c:
                no_regret_cells.append(eb_c <= 1.05 * min(bs_c))
                if not math.isnan(ols_c):
                    half_cut_cells.append(ols_c > 0 and eb_c <= 0.5 * ols_c)
                continue
        if bases and not math.isnan(eb):
            no_regret_cells.append(eb <= 1.05 * min(bases))
            if not math.isnan(ols_b) and ols_b > 0:
                half_cut_cells.append(eb <= 0.5 * ols_b)
    n = len(no_regret_cells)
    if n:
        nr = float(np.mean(no_regret_cells))
        hc = float(np.mean(half_cut_cells)) if half_cut_cells else float("nan")
        print(f"no-regret (<= 1.05x best baseline): {nr:.1%} of {n} harmful "
              f"cells -> need >= 95% : {verdict(nr >= 0.95)}")
        print(f">=50% bias cut vs OLS: {hc:.1%} -> need >= 70% : "
              f"{verdict(hc >= 0.70)}")
    else:
        print("PENDING: harmful cells incomplete")


def check_detection():
    df = load("nullcal")
    print("\n=== WP 2.3 detection (size gates; power needs power parquet) ===")
    if df is None:
        print("PENDING: no consolidated nullcal parquet")
        return
    size_ok = []
    for cid, g in df.groupby("config_id"):
        for stat in ("rej_s0_tw99", "rej_s2_bonf"):
            pass
    # MC-calibrated sizes are computed at analysis time; here we evaluate the
    # analytic-threshold sizes as reported diagnostics and count S2-Bonferroni.
    s2 = df.groupby("config_id")["rej_s2_bonf"].mean()
    ok = ((s2 >= 0.02) & (s2 <= 0.12)).mean()  # loose pre-gate sanity band
    print(f"S2 Bonferroni raw size within sane band [0.02, 0.12] in "
          f"{ok:.0%} of null cells ({len(s2)} cells); formal gate uses "
          f"MC-calibrated thresholds from matched nulls (analysis phase).")


def main():
    check_correctness()
    check_estimation()
    check_detection()
    print("\nNote: all verdicts re-run mechanically from parquet; thresholds"
          "\ncome from research plan Section 6 + docs/phase2_preregistration.md")


if __name__ == "__main__":
    main()
