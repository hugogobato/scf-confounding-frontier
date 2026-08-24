"""Generate pilot figures from data/pilot (parquet + cell_summary.csv only).

Figures regenerate from parquet via one command:
    python3 code/make_pilot_figures.py
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PILOT = ROOT / "data" / "pilot"
FIGS = ROOT / "figures"
FIGS.mkdir(exist_ok=True)

THETA_ORDER = [0.0, 30.0, 90.0]
PROFILES = ["sub", "mixed", "super", "onespike"]


def load():
    s = pd.read_csv(PILOT / "cell_summary.csv")
    raw = pd.read_parquet(PILOT / "pilot_results.parquet",
                          columns=["config_id", "estimator", "rep", "lam_max_cov",
                                   "overlap1", "tw_stat", "outlier99"])
    rv_path = PILOT / "revisit_summary.csv"
    rv = pd.read_csv(rv_path) if rv_path.exists() else None
    return s, raw, rv


def fig_phase_regions(s: pd.DataFrame):
    main = s[s["label"] == "main"]
    aux = s[s["label"] == "aux"]
    profiles = ["sub", "mixed", "super"]
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5), sharex=True)
    for col, prof in enumerate(profiles):
        sub = main[main["profile"] == prof].sort_values("theta_deg")
        stable = sub[sub["stable_mean_ols"]] if "stable_mean_ols" in sub.columns else sub
        unstable = sub[~sub["stable_mean_ols"]] if "stable_mean_ols" in sub.columns else sub.iloc[0:0]
        for row, (metric, title) in enumerate([
                ("bias_ratio_ols_sim", "OLS total bias ratio"),
                ("bias_ratio_conf_sim", "confounding-attributed bias ratio")]):
            ax = axes[row, col]
            vals = stable[metric].to_numpy(float)
            sc = ax.scatter(stable["c"], stable["theta_deg"], c=vals,
                            cmap="viridis", vmin=0.0,
                            vmax=max(0.6, np.nanmax(vals)),
                            s=260, marker="o", edgecolor="k")
            # unstable c=1.0 cells: mean-bias undefined (E[(X'X)^-1] diverges)
            if len(unstable):
                ax.scatter(unstable["c"], unstable["theta_deg"], marker="x",
                           c="lightgray", s=200, edgecolor="k", linewidths=1.5,
                           label="c=1: mean-bias undefined")
            # outlier visibility marker: TW99 outlier rate >= 95%
            out = sub["outlier_rate_tw99"].to_numpy(float) >= 0.95
            for x, y, o in zip(sub["c"], sub["theta_deg"], out):
                ax.annotate("^" if o else "o", (x, y), ha="center", va="center",
                            color="red" if o else "white", fontsize=11, weight="bold")
            if len(unstable) and row == 0 and col == 0:
                ax.legend(fontsize=7, loc="upper left")
            ax.set_xscale("log")
            ax.set_xticks([0.2, 1.0, 5.0])
            ax.set_xticklabels(["0.2", "1", "5"])
            ax.set_yticks(THETA_ORDER)
            ax.set_ylim(-15, 105)
            if row == 0:
                ax.set_title(f"profile = {prof}", fontsize=12)
            if col == 0:
                ax.set_ylabel(f"{title}\ntheta (deg)")
            else:
                ax.set_ylabel("theta (deg)")
            ax.set_xlabel("c = p/n")
            plt.colorbar(sc, ax=ax, fraction=0.046)
    fig.suptitle(
        "Pilot phase regions (n=2000, r=3, g=1). Red ^ = BBP outlier detected "
        "(TW99 rate>=95%), white o = none.\nTop row: total OLS bias; bottom row: "
        "gamma-attributed bias (twin arms at c>1).", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(FIGS / "pilot_phase_regions.png", dpi=160)
    plt.close(fig)


def fig_de_overlay(s: pd.DataFrame, raw: pd.DataFrame, rv=None):
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9))

    # (a) OLS bias: predicted vs simulated
    ax = axes[0, 0]
    ok = s.dropna(subset=["bias_ratio_ols_sim", "bias_ratio_ols_pred"])
    for lab, mk in [("main", "o"), ("aux", "^")]:
        subm = ok[ok["label"] == lab]
        ax.scatter(subm["bias_ratio_ols_pred"], subm["bias_ratio_ols_sim"],
                   marker=mk, c=np.log10(subm["c"]), cmap="coolwarm",
                   s=70, edgecolor="k", zorder=3)
    if rv is not None:
        ax.scatter(rv["bias_ratio_pred"], rv["bias_ratio_conf_sim_twin"],
                   marker="*", c="limegreen", s=140, edgecolor="k",
                   zorder=4, label="revisit (1000-rep twins)")
        ax.legend(fontsize=8)
    lim = [1e-3, 2]
    ax.plot(lim, lim, "k--", lw=1, label="y = x")
    ax.set_xlim(lim); ax.set_ylim(lim); ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("predicted ||E[beta_hat] - beta|| (DE)")
    ax.set_ylabel("simulated")
    ax.set_title("(a) OLS bias overlay, all cells")
    ax.legend()

    # (b) ridge curves for two representative cells (revisit means carry the
    # per-lambda keys; first-pass npz predate the key fix)
    ax = axes[0, 1]
    lams = np.array([0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0])
    picks = [("mixed", 0.8, 0.0), ("super", 5.0, 0.0)]
    import sys
    sys.path.insert(0, str(ROOT / "code"))
    from de_formulas import ridge_bias_vector, ridge_capture
    rv_means_dir = PILOT / "revisit_means"
    for k, (prof, c, th) in enumerate(picks):
        row = s[(s["profile"] == prof) & (s["c"] == c)
                & (s["theta_deg"] == th) & (s["label"] == "main")]
        if row.empty:
            continue
        cid = row["config_id"].iloc[0]
        mp = None
        for cand in (rv_means_dir / f"{cid}.npz", PILOT / "means" / f"{cid}.npz"):
            if Path(cand).exists():
                m = np.load(cand)
                if f"ridge_fixed@{lams[0]}" in m:
                    mp = m
                    break
        if mp is None:
            continue
        sim = [np.linalg.norm(mp[f"ridge_fixed@{lm}"].astype(float))
               for lm in lams]
        ls = row[[f"l{j}" for j in (1, 2, 3)]].iloc[0].to_numpy(float)
        gam = np.array([np.cos(np.deg2rad(th)), np.sin(np.deg2rad(th)), 0.0])
        if c <= 1.0:
            pred = [float(np.linalg.norm(ridge_bias_vector(ls, gam[:len(ls)], lm)))
                    for lm in lams]
        else:
            pred = []
            for lm in lams:
                capr = ridge_capture(ls, lm, c)
                coef = capr * np.sqrt(ls) / (1.0 + ls) * gam[:len(ls)]
                pred.append(float(np.linalg.norm(coef)))
        ax.plot(lams, sim, "o-", color=f"C{k}", label=f"sim {prof}, c={c}")
        ax.plot(lams, pred, "--", color=f"C{k}", label=f"pred {prof}, c={c}")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("ridge lambda (Sigma-scale)")
    ax.set_ylabel("||bias||")
    ax.set_title("(b) Ridge bias vs lambda")
    ax.legend(fontsize=8)

    # (c) BBP outlier location
    ax = axes[1, 0]
    ols = raw[raw["estimator"] == "ols"]
    agg = ols.groupby("config_id").agg(lam=("lam_max_cov", "mean")).reset_index()
    m = s.merge(agg, on="config_id")
    m["bbp_pred_eff"] = m.apply(
        lambda r_: max(r_["bbp_pred"], (1 + np.sqrt(r_["c"])) ** 2), axis=1)
    ax.scatter(m["bbp_pred_eff"], m["lam"], c=np.log10(m["c"]),
               cmap="coolwarm", s=60, edgecolor="k", zorder=3)
    lim = [(1 + np.sqrt(0.2)) ** 2 * 0.7, 400]
    ax.plot(lim, lim, "k--", lw=1)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("predicted lambda_max (BBP mu(l_max), bulk edge if subcritical)")
    ax.set_ylabel("simulated mean lambda_max")
    ax.set_title("(c) Outlier location overlay")

    # (d) BGN overlap of top spike
    ax = axes[1, 1]
    mm = m.dropna(subset=["xi1_sim"])
    sup = mm["bbp_pred_eff"] > (1 + np.sqrt(mm["c"])) ** 2 * 1.02
    ax.scatter(mm.loc[sup, "xi1_pred"], mm.loc[sup, "xi1_sim"],
               c=np.log10(mm.loc[sup, "c"]), cmap="coolwarm",
               s=60, edgecolor="k", zorder=3, label="supercritical cells")
    ax.scatter(mm.loc[~sup, "xi1_pred"], mm.loc[~sup, "xi1_sim"],
               facecolors="none", edgecolor="gray", label="subcritical cells")
    ax.plot([0, 1], [0, 1], "k--", lw=1)
    ax.set_xlim(-0.03, 1.03); ax.set_ylim(-0.03, 1.03)
    ax.set_xlabel("predicted xi_1"); ax.set_ylabel("simulated mean overlap_1")
    ax.set_title("(d) BGN eigenvector overlap")
    ax.legend(fontsize=8)

    fig.suptitle("DE overlays vs simulation (WP 1.5)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIGS / "pilot_de_overlay.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    s, raw, rv = load()
    fig_phase_regions(s)
    fig_de_overlay(s, raw, rv)
    print("figures written to", FIGS)
