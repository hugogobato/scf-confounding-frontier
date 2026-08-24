"""Evaluate the G0/G2 numeric conditions of research plan Section 6 / WP 1.5
on data/pilot/cell_summary.csv and print the verdict table used by
docs/pilot_memo.md."""
from __future__ import annotations

import os
import sys
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
S = pd.read_csv(ROOT / "data" / "pilot" / "cell_summary.csv")
RV_PATH = ROOT / "data" / "pilot" / "revisit_summary.csv"
RV = pd.read_csv(RV_PATH) if RV_PATH.exists() else None

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 50)

if RV is not None:
    print("=== revisit cells (1000-rep mains / 400-rep paired twins) ===")
    rv_show = RV[["label", "profile", "theta_deg", "c",
                  "bias_ratio_conf_sim_twin", "bias_ratio_pred",
                  "proj_pred_sim", "dev_proj_pct", "dev_conf_pct",
                  "outlier_rate_tw99"]]
    print(rv_show.round(4).to_string(index=False))
    dev_main = RV.loc[RV.label == "main", "dev_proj_pct"]
    print(f"\nrevisit main-cell DIRECTIONAL overlay: median "
          f"{dev_main.median():.3f}% | max {dev_main.max():.3f}%")
    for _, r in RV.iterrows():
        curves = r.get("ridge_curves")
        if isinstance(curves, str) and curves != "[]":
            import ast
            rows = ast.literal_eval(curves)
            if rows:
                devs = [100 * abs(x["sim"] - x["pred"]) / max(x["pred"], 1e-12)
                        for x in rows]
                print(f"  ridge {r['profile']}/c={r['c']}: median dev "
                      f"{np.median(devs):.2f}% | max {max(devs):.2f}%")

cols = ["label", "profile", "theta_deg", "c", "bias_ratio_ols_sim",
        "bias_ratio_ols_pred", "dev_ols_pct", "bias_ratio_conf_sim",
        "outlier_rate_tw99", "lam_max_sim_mean", "bbp_pred_eff" if False else "bbp_pred",
        "xi1_sim", "xi1_pred", "ridge_dev_mean_pct"]
show = S[cols].copy()
show["lam_max_sim_mean"] = show["lam_max_sim_mean"].round(2)
show["bbp_pred"] = show["bbp_pred"].round(2)
print("=== per-cell summary ===")
print(show.round(4).to_string(index=False))

# effective predicted lambda_max (BBP location if supercritical else bulk edge)
S["pred_lam_eff"] = np.maximum(
    S["bbp_pred"], (1 + np.sqrt(S["c"])) ** 2)
lam_dev = 100 * abs(S["lam_max_sim_mean"] - S["pred_lam_eff"]) / S["pred_lam_eff"]
xi_dev = 100 * abs(S["xi1_sim"] - S["xi1_pred"]) / np.maximum(S["xi1_pred"], 1e-3)
supercritical = S["bbp_pred"] > (1 + np.sqrt(S["c"])) ** 2 * 1.02

print("\n=== region existence ===")
stable = S[S["stable_mean_ols"]] if "stable_mean_ols" in S.columns else S
sup_s = supercritical.reindex(stable.index)
invisible_harmful = stable[(~sup_s) & (stable["bias_ratio_ols_sim"] >= 0.2)]
visible_harmless_aux = S[(S["label"] == "aux")
                         & (S["outlier_rate_tw99"] >= 0.95)
                         & (S["bias_ratio_conf_sim"] <= 0.02)]
if RV is not None:
    # authoritative: directional (projection) estimates on paired twins
    rv_aux_ok = RV[(RV["label"] == "aux") & (RV["outlier_rate_tw99"] >= 0.95)
                   & (RV["c"] < 1.0) & (RV["proj_pred_sim"] <= 0.02)]
    print(f"revisit aux cells, TW99 outlier present & directional conf-bias "
          f"<= 0.02: {len(rv_aux_ok)} -> "
          f"{[(r.profile, r.theta_deg, r.c, round(r.proj_pred_sim, 4)) for r in rv_aux_ok.itertuples()]}")
print(f"cells invisible-yet-harmful (no BBP outlier & total bias >= 0.2): "
      f"{len(invisible_harmful)} -> "
      f"{[(r.profile, r.theta_deg, r.c) for r in invisible_harmful.itertuples()]}")
print(f"cells visible-yet-harmless (TW99 outlier rate>=95% & bias<=0.02): "
      f"{len(visible_harmless_aux)} -> "
      f"{[(r.profile, r.theta_deg, r.c) for r in visible_harmless_aux.itertuples()]}")

print("\n=== DE overlay deviations (stable c<1 cells for OLS mean-bias) ===")
stable = S[S["stable_mean_ols"]] if "stable_mean_ols" in S.columns else S
ols_devs = stable.loc[stable.label == "main", "dev_ols_pct"]
print(f"OLS overlay: median {ols_devs.median():.2f}% | max {ols_devs.max():.2f}% "
      f"| cells <=10%: {(ols_devs <= 10).mean() * 100:.0f}%")
conf = S.dropna(subset=["bias_ratio_conf_sim"])
conf = conf[conf["c"] > 1.0]  # c=1 cells are unstable (memo F-B); c>1 runs twins
if len(conf):
    conf_dev = 100 * abs(conf["bias_ratio_conf_sim"] - conf["bias_ratio_conf_pred"]) \
        / conf["bias_ratio_conf_pred"]
    print(f"gamma-attributed (twin, c>1): median {conf_dev.median():.2f}% "
          f"| max {conf_dev.max():.2f}%")
print(f"lambda_max overlay (all cells): median {lam_dev.median():.2f}% "
      f"| max {lam_dev.max():.2f}%")
sep = supercritical & (S["l1"] != S["l2"])  # separated spikes only (memo F-D)
print(f"overlap xi_1 overlay (separated-spike cells): "
      f"median {xi_dev[sep].median():.2f}% "
      f"| max {xi_dev[sep].max():.2f}% (equal-spike cells rotate, excluded)")
rd = S["ridge_dev_mean_pct"].dropna()
if len(rd):
    print(f"ridge overlay first pass (7 lambdas/cell): median {rd.median():.2f}% "
          f"| max {rd.max():.2f}%")

verdict = {
    "region_invisible_harmful": len(invisible_harmful) > 0,
    "region_visible_harmless": (len(visible_harmless_aux) > 0)
        or (RV is not None and
            len(RV[(RV["label"] == "aux") & (RV["outlier_rate_tw99"] >= 0.95)
                   & (RV["c"] < 1.0) & (RV["proj_pred_sim"] <= 0.02)]) > 0),
    "ols_overlay_le_10pct_stable_cells": bool((ols_devs <= 10).mean() >= 0.90),
    "lam_overlay_reasonable": bool(lam_dev.median() < 10),
}
print("\n=== G0/G2 verdict inputs ===")
for k, v in verdict.items():
    print(f"{k}: {'PASS' if v else 'FAIL'}")
print("recommendation:", "GO" if all(verdict.values()) else
      "PIVOT/KILL per Phase 1 give-up rules (see memo)")
