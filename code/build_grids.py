"""SCF Phase 2 frozen grids, AMENDMENT v2 (docs/phase2_preregistration.md,
deviation D6). Emits job dicts consumed identically by local execution
(code/run_sweep.py) and Colab shards (code/make_shards.py).

Grid edits AFTER data generation are forbidden; this file is hashed into
shard manifests. Memory ceiling (binding): no cell with n*p > 2.4e8 or
min(n, p) > 8000.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

THETA_MAIN = math.pi / 6
# Amendment v2 (D6a): c = 10 dropped everywhere.
C_LIST = [0.1, 0.2, 0.5, 0.8, 2.0, 5.0]
C_LIST_DET = C_LIST + [1.0]  # detection-only cells may use c = 1.0
PROFILES6 = [(1, "sub"), (1, "super"), (5, "sub"), (5, "mixed"), (5, "super"),
             (25, "mixed")]


def profile_l(name: str, c: float, r: int) -> tuple[float, ...]:
    sc = math.sqrt(c)
    if name == "sub":
        return tuple([0.5 * sc] * r)
    if name == "mixed":
        return tuple([3.0 * sc] + [0.5 * sc] * (r - 1))
    if name == "super":
        return tuple([3.0 * sc] * r)
    raise ValueError(name)


def ceiling_ok(n: int, p: int) -> bool:
    return n * p <= 2.4e8 and min(n, p) <= 8000


def _job(cfg_kwargs: dict, mode: str, sweep: str, out_root: str,
         reps: int) -> dict:
    from simulator import Config

    cfg = Config(**cfg_kwargs)
    cid = cfg.cid
    return {
        "config": cfg_kwargs,
        "config_id": cid,
        "mode": mode,
        "sweep": sweep,
        "reps": reps,
        "raw_path": f"{out_root}/{sweep}/raw/{cid}.parquet",
        "means_path": f"{out_root}/{sweep}/means/{cid}.npz",
    }


def _base(n, p, r, prof, theta, label="main", **kw):
    kw.setdefault("q_fixed", True)
    kw.setdefault("twin_gamma0", p / n > 1.0)
    return dict(n=n, p=p, r=r, l=profile_l(prof, p / n, r), theta=theta,
                profile=prof, label=label, **kw)


def grid_correctness(out_root="data/sim") -> list[dict]:
    jobs = []
    # mainA: n = 500, full profiles, reps 1000
    for c in C_LIST:
        p = int(round(c * 500))
        for r, prof in PROFILES6:
            jobs.append(_job(_base(500, p, r, prof, THETA_MAIN),
                             "correctness", "correctness", out_root, 600))
    # mainB: n = 2000, full profiles, reps 350 (D6b)
    for c in C_LIST:
        p = int(round(c * 2000))
        for r, prof in PROFILES6:
            jobs.append(_job(_base(2000, p, r, prof, THETA_MAIN),
                             "correctness", "correctness", out_root, 200))
    # alignS: alignment sensitivity at n = 2000, r = 5 profiles only (D6b)
    for th in (0.0, math.pi / 2):
        for c in C_LIST:
            p = int(round(c * 2000))
            for r, prof in ((5, "sub"), (5, "mixed")):
                jobs.append(_job(_base(2000, p, r, prof, th),
                                 "correctness", "correctness", out_root, 150))
    # deepN: skinny-c branch only at n = 8000 (D6b)
    for c in (0.1, 0.2):
        p = int(round(c * 8000))
        for r, prof in ((1, "sub"), (5, "mixed"), (1, "super")):
            jobs.append(_job(_base(8000, p, r, prof, THETA_MAIN),
                             "correctness", "correctness", out_root, 150))
    return jobs


def grid_nullcal(out_root="data/sim") -> list[dict]:
    def null_job(n, p, r, prof, th, reps):
        return _job(_base(n, p, r, prof, th, g=0.0, twin_gamma0=False),
                    "nullcal", "nullcal", out_root, reps)

    jobs = []
    # core: includes detection-only c = 1.0
    for c in C_LIST_DET:
        p = int(round(c * 2000))
        for prof in ("sub", "mixed"):
            jobs.append(null_job(2000, p, 3, prof, THETA_MAIN, 1200))
    # power-grid matched theta nulls
    for c in (0.2, 0.8, 2.0):
        p = int(round(c * 2000))
        for prof in ("sub", "mixed", "super"):
            for th in (THETA_MAIN, math.pi / 2):
                jobs.append(null_job(2000, p, 3, prof, th, 500))
    # n-ladder: size at small n, trend-only at n = 4000 (D6c)
    for c in (0.2, 2.0):
        p = int(round(c * 500))
        for prof in ("sub", "mixed"):
            jobs.append(null_job(500, p, 3, prof, THETA_MAIN, 4000))
    jobs.append(null_job(4000, 800, 3, "mixed", THETA_MAIN, 400))
    seen, uniq = set(), []
    for j in jobs:
        if j["config_id"] not in seen:
            seen.add(j["config_id"])
            uniq.append(j)
    return uniq


def grid_power(out_root="data/sim") -> list[dict]:
    g_grid = [0.15, 0.4, 0.8, 1.6, 3.2]  # D6d/D7
    jobs = []
    for c in (0.2, 0.8, 2.0):
        p = int(round(c * 2000))
        for prof in ("sub", "mixed"):  # D7: super anchored by nullcal/scree
            for th in (THETA_MAIN, math.pi / 2):  # D6d
                for g in g_grid:
                    jobs.append(_job(_base(2000, p, 3, prof, th, g=g),
                                     "power", "power", out_root, 300))
    return jobs


def grid_alignment(out_root="data/sim") -> list[dict]:
    p = int(round(0.8 * 2000))
    thetas = [math.pi * k / 12 for k in range(12)]
    return [
        _job(_base(2000, p, 3, "mixed", th, g=1.0),
             "alignment", "alignment", out_root, 800)  # D6d/D7
        for th in thetas
    ]


def grid_estimation(out_root="data/sim") -> list[dict]:
    jobs = []
    # D7: c-set {0.2,0.8,2,5}; theta {pi/6, pi/2} at both n;
    #     reps 200 (n=500) / 140 (n=2000), c=5&n=2000 at 100
    for n, base_reps in ((500, 200), (2000, 140)):
        for c in (0.2, 0.8, 2.0, 5.0):
            p = int(round(c * n))
            if not ceiling_ok(n, p):
                continue
            reps = base_reps if not (c == 5.0 and n == 2000) else 100
            for r in (1, 5, 25):
                for prof in ("sub", "mixed"):
                    for th in (THETA_MAIN, math.pi / 2):
                        jobs.append(_job(_base(n, p, r, prof, th),
                                         "estimation", "estimation", out_root,
                                         reps))
    # deep evidence point on the cheap branch (D6e)
    jobs.append(_job(_base(8000, 1600, 5, "mixed", THETA_MAIN),
                     "estimation", "estimation", out_root, 80))
    # rung 4 baseline-favorable: sparse confounding + aligned beta
    for c in (0.2, 0.8):
        p = int(round(c * 2000))
        jobs.append(_job(_base(2000, p, 3, "mixed", THETA_MAIN,
                               conf_kind="sparse", beta_kind="aligned"),
                         "estimation", "estimation", out_root, 300))
    return jobs


def grid_crossover(out_root="data/sim") -> list[dict]:
    jobs = []
    for c in (0.2, 0.8, 2.0):
        p = int(round(c * 2000))
        for g in (0.25, 0.5, 1.0, 2.0, 4.0):
            l = tuple([3.0 * math.sqrt(c)] + [0.5 * math.sqrt(c)] * 2)
            jobs.append(_job(dict(n=2000, p=p, r=3, l=l, theta=THETA_MAIN,
                                  g=g, profile="mixed", label="crossover",
                                  q_fixed=True,
                                  twin_gamma0=(p > 2000)),
                             "crossover", "crossover", out_root, 150))  # D6e/D7
    return jobs


ROBUST_VARIANTS = {
    "V0_gauss": {},
    "V1_t5": {"error_law": "t5"},
    "V2_rademacher_half": {"loading_kind": "rademacher_half"},
    "V3_hetero_u": {"hetero_u": True},
    "V4_corr_f": {"corr_factors": True},
    "V6_sparse_conf": {"conf_kind": "sparse"},
}


def grid_robustness(out_root="data/sim") -> list[dict]:
    jobs = []
    for cname, extra in ROBUST_VARIANTS.items():
        for c in (0.2, 2.0):  # D6e: two aspect anchors
            p = int(round(c * 2000))
            for prof in ("sub", "mixed"):
                kw = dict(extra)
                jobs.append(_job(_base(2000, p, 5, prof, THETA_MAIN,
                                       label=cname, **kw),
                                 "robustness", "robustness", out_root, 250))
    # V5 r-misspecification
    for dr in (-1, 1):
        for c in (0.2, 0.8):
            p = int(round(c * 2000))
            jobs.append(_job(_base(2000, p, 5, "mixed", THETA_MAIN,
                                   label=f"V5_r{dr:+d}", r_misspec=dr),
                             "robustness", "robustness", out_root, 125))
    return jobs


def grid_m2(out_root="data/sim") -> list[dict]:
    jobs = []
    for c in (0.2, 0.8, 2.0):
        p = int(round(c * 2000))
        jobs.append(_job(_base(2000, p, 5, "mixed", THETA_MAIN,
                               m2_treatment=True, m2_tau=1.0, delta_g=0.3,
                               label="m2_weak"),
                         "m2", "m2", out_root, 150))  # D6f/D7
    return jobs


def grid_scaling(out_root="data/sim") -> list[dict]:
    """Timing/memory envelope study: 2 reps per size, correctness mode."""
    sizes = [(1000, 1000), (2000, 8000), (4000, 8000), (8000, 1600)]
    jobs = []
    for n, p in sizes:
        if not ceiling_ok(n, p):
            continue
        c = p / n
        l = tuple([3.0 * math.sqrt(c)] + [0.5 * math.sqrt(c)] * 4)
        jobs.append(_job(dict(n=n, p=p, r=5, l=l, theta=THETA_MAIN, g=1.0,
                              profile="mixed", label="scaling", q_fixed=True,
                              twin_gamma0=(c > 1.0)),
                         "correctness", "scaling", out_root, 2))
    return jobs


def all_grids(out_root="data/sim") -> dict[str, list[dict]]:
    return {
        "correctness": grid_correctness(out_root),
        "nullcal": grid_nullcal(out_root),
        "power": grid_power(out_root),
        "alignment": grid_alignment(out_root),
        "estimation": grid_estimation(out_root),
        "crossover": grid_crossover(out_root),
        "robustness": grid_robustness(out_root),
        "m2": grid_m2(out_root),
        "scaling": grid_scaling(out_root),
    }


def main():
    grids = all_grids()
    out = Path("configs")
    out.mkdir(exist_ok=True)
    summary = {}
    for name, jobs in grids.items():
        path = out / f"grid_{name}.json"
        path.write_text(json.dumps(jobs, indent=1))
        total_reps = sum(j["reps"] for j in jobs)
        summary[name] = {"cells": len(jobs), "reps": total_reps,
                         "file": str(path)}
        print(f"{name}: {len(jobs)} cells, {total_reps} reps -> {path}")
    blob = json.dumps(summary, sort_keys=True)
    print("grids hash:", hashlib.sha256(blob.encode()).hexdigest()[:12])


if __name__ == "__main__":
    main()
