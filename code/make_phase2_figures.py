"""Preregistered Phase 2 figures from results/*.csv (derived from parquet
by code/phase2_analysis.py; regeneration order documented in the execution
memo). Outputs to figures/ per docs/phase2_preregistration.md names."""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
FIG = ROOT / "figures"


def de_overlay_grid():
    ov = pd.read_csv(RES / "correctness_overlays.csv")
    fig, ax = plt.subplots(figsize=(7, 4))
    for nval, g in ov.groupby("n"):
        m = g.groupby("c").ols_reldev.median()
        ax.plot(m.index, m.values * 100, marker="o", label=f"n={nval}")
    ax.axhline(10, ls="--", c="k", lw=0.8)
    ax.text(0.02, 10.6, "gate: 10%", fontsize=8)
    ax.set_xlabel("c = p/n")
    ax.set_ylabel("median |sim - pred| / pred (%)")
    ax.set_title("OLS mean-bias DE overlay (per-cell, median over profiles)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG / "de_overlay_grid.pdf")
    plt.close(fig)


def power_surface_vs_frontier():
    ps = pd.read_csv(RES / "power_surface.csv")
    strata = ps.groupby(["c", "profile", "theta"])
    fig, axes = plt.subplots(1, len(ps[["c", "profile"]].drop_duplicates()),
                             figsize=(11, 3.2), squeeze=False,
                             sharey=True)
    col = 0
    for (c, prof), g in ps.groupby(["c", "profile"]):
        for th, gg in g.groupby("theta"):
            gg = gg.sort_values("g")
            ax = axes[0][min(col, axes.shape[1] - 1)]
            ax.plot(gg.g, gg.emp_power_S2, marker="o", ms=3,
                    label=f"theta={th:.2f}")
            gp = gg.g_pred_S2.iloc[0]
            if str(gp) != "inf" and gp == gp:
                ax.axvline(gp, ls=":", lw=0.9)
        axes[0][min(col, axes.shape[1] - 1)].set_title(f"c={c}, {prof}",
                                                       fontsize=9)
        col += 1
    axes[0][0].set_ylabel("S2_cal empirical power")
    axes[0][0].set_ylim(-0.03, 1.05)
    for a in axes[0]:
        a.set_xlabel("g")
    fig.suptitle("Power surfaces vs F12-law predicted frontier "
                 "(dotted vertical)", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "power_surface_vs_frontier.pdf")
    plt.close(fig)


def alignment_stress():
    al = pd.read_csv(RES / "alignment_stress.csv")
    fig, ax = plt.subplots(figsize=(6, 3.5))
    for colname, lbl in [("pow_S2_cal", "S2 cal"), ("pow_S1_analytic",
                                                   "S1 analytic"),
                         ("pow_B1", "B1 F-PCS"), ("pow_S0", "S0 scree")]:
        if colname in al:
            ax.plot(al.theta, al[colname], marker="o", ms=3, label=lbl)
    ax.set_xlabel("theta")
    ax.set_ylabel("power (n=2000, c=0.8, mixed, g=1)")
    ax.set_title("Alignment stress: S2 tracks supercritical-aligned mass")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "alignment_stress.pdf")
    plt.close(fig)


def crossover_curves():
    cv = pd.read_csv(RES / "crossover_curves.csv")
    fig, axes = plt.subplots(1, cv.c.nunique(), figsize=(10, 3.2),
                             sharey=True, squeeze=False)
    tags = ["ols", "pca_onatski", "eb_spectral", "eb_cv_tau"]
    for i, (c, g) in enumerate(cv.groupby("c")):
        ax = axes[0][i]
        for t in tags:
            if t in g:
                ax.plot(g.g, g[t], marker="o", ms=3, label=t)
        ax.set_xscale("log", base=2)
        ax.set_yscale("log")
        ax.set_title(f"c={c}", fontsize=9)
        ax.set_xlabel("g")
    axes[0][0].set_ylabel("mean-bias norm")
    axes[0][0].legend(fontsize=7)
    fig.suptitle("Crossover strip: bias vs confounding strength", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "crossover_curves.pdf")
    plt.close(fig)


def bias_phase_diagram_empirical():
    det = pd.read_csv(RES / "estimation_cell_detail.csv")
    det = det.dropna(subset=["ols_conf"])
    det["half_cut"] = det.eb_conf <= 0.5 * det.ols_conf
    piv = det.pivot_table(index="c", columns="theta",
                          values="half_cut", aggfunc="mean")
    fig, ax = plt.subplots(figsize=(5.5, 3.4))
    im = ax.imshow(piv.values, cmap="RdYlGn", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(piv.columns)))
    ax.set_xticklabels([f"{t:.2f}" for t in piv.columns])
    ax.set_yticks(range(len(piv.index)))
    ax.set_yticklabels(piv.index)
    ax.set_xlabel("theta")
    ax.set_ylabel("c")
    ax.set_title("Share of harmful cells with >=50% OLS bias cut (eb_spectral)",
                 fontsize=9)
    fig.colorbar(im, ax=ax, shrink=0.85)
    fig.tight_layout()
    fig.savefig(FIG / "bias_phase_diagram_empirical.pdf")
    plt.close(fig)


def lecam_probe_auc():
    lc = pd.read_csv(RES / "lecam_probe_auc.csv")
    fig, ax = plt.subplots(figsize=(6, 3.5))
    for (c, prof), g in lc.groupby(["c", "profile"]):
        blind = g.g_pred.astype(str) == "inf"
        lbl = f"c={c} {prof}" + (" (blind)" if blind.all() else "")
        ax.plot(g.g, g.auc_gbm, marker="o", ms=3, label=lbl,
                ls="--" if blind.all() else "-")
    ax.axhline(0.55, ls="--", c="k", lw=0.8)
    ax.text(0.02, 0.57, "declaration threshold 0.55", fontsize=7)
    ax.set_xlabel("g")
    ax.set_ylabel("GBM probe AUC vs matched nulls")
    ax.set_title("Le Cam probe across the s-grid")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG / "lecam_probe_auc.pdf")
    plt.close(fig)


if __name__ == "__main__":
    de_overlay_grid()
    power_surface_vs_frontier()
    alignment_stress()
    crossover_curves()
    bias_phase_diagram_empirical()
    lecam_probe_auc()
    print("figures written:", sorted(p.name for p in FIG.glob('*.pdf')))
