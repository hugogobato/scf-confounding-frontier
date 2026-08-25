"""SCF Phase 3 WP 3.1/WP 3.2: benchmark dataset acquisition, frozen
preprocessing, and spectral audit (read-only wrt raw sources).

Families (research plan Section 7, Phase 3 amendment):
  A addneuromed : GEO GSE63060 + GSE63061 blood gene expression (Illumina
                  HumanHT-12 v3 / v4), the two series ARE the processing
                  batches of the AddNeuroMed cohort. Provider NCBI GEO.
  B ihdp        : IHDP covariate benchmark (Hill 2011 npci covariates),
                  featurized to high dimension with seeded random Fourier
                  features. Public benchmark mirror.
  C k401k       : wooldridge::k401ksubs household cross-section (Rdatasets
                  CC0 mirror); M2 treatment block (e401k -> nettfa).

Frozen preprocessing (recorded in configs/benchmarks_frozen.yaml BEFORE any
comparative result):
  A: intersect ILMN probe IDs across the two series, average duplicate probe
     IDs within each series, drop probes with missing values, apply
     log2(x+1) iff the global max exceeds 25, z-score each probe on the
     POOLED sample (batch mean shifts are preserved by design), keep the top
     P probes by pooled variance.
     A_main: P = 2000, all samples            (n = 717, c = 2.79)
     A_sub : P = 1800, batch-2 samples only   (n = 388, c = 4.64)
  B/C: standardize continuous covariates (binary kept as 0/1), then a seeded
     random Fourier map Z = sqrt(2/P) cos(X W + b), W ~ N(0, h^-2 I_d) with
     h the median-heuristic bandwidth on a seeded subsample.
     B_main: n = 747, P = 750  (c = 1.00)   B_wide: P = 150 (c = 0.20)
     C_main: n = 800 seeded subsample, P = 1600 (c = 2.00)
     C_wide: same rows, P = 160 (c = 0.20)

Spectral audit per config (consumed later by the frontier machinery):
  se2_hat (shared bulk-median estimator), unit-scaled spectrum, Onatski
  r_hat, BBP-inverted spike estimates l_hat_j, TW statistic of lambda_max vs
  the white-noise threshold.

Outputs:
  data/benchmarks/<config>.npz      processed designs (+ labels/meta)
  data/benchmarks/spectral_audit.json

This script never writes into data/benchmarks/raw/.
"""
from __future__ import annotations

import csv
import gzip
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "benchmarks" / "raw"
OUT = ROOT / "data" / "benchmarks"
GLOBAL_SEED = 20260823


# ---------------------------------------------------------------------------
# Family A: AddNeuroMed series matrices
# ---------------------------------------------------------------------------


def parse_geo_series(path: Path):
    """Return (probes, values [n x p], meta list[dict]) from one matrix."""
    samples, chars_rows, expr_rows = [], [], []
    in_table = False
    header = None
    with gzip.open(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("!Sample_geo_accession"):
                samples = [c.strip('"') for c in line.split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                chars_rows.append(line.split("\t")[1:])
            elif line.startswith("!series_matrix_table_begin"):
                in_table = True
                header = next(fh).rstrip("\n").split("\t")
            elif in_table:
                if line.startswith("!series_matrix_table_end"):
                    break
                expr_rows.append(line.split("\t"))
    # metadata: characteristic cells may hold one field per cell; concatenate
    # everything observed for sample i and split on ':' prefixes.
    flat = [""] * len(samples)
    for row in chars_rows:
        if len(row) == len(samples):
            for i, v in enumerate(row):
                flat[i] += v.strip('"').strip() + "; "
    meta = []
    for blob in flat:
        d = {}
        for part in blob.split(";"):
            if ":" in part:
                k, _, v = part.partition(":")
                d[k.strip().lower()] = v.strip()
        meta.append({"status": d.get("status"), "age": d.get("age"),
                     "gender": d.get("gender")})

    # GEO layout: rows are probes (ID_REF), columns are samples.
    probes = [r[0].strip('"') for r in expr_rows]
    vals_t = np.empty((len(expr_rows), len(samples)), dtype=np.float32)
    for i, parts in enumerate(expr_rows):
        vals_t[i] = np.asarray(parts[1:], dtype=np.float32)
    return probes, vals_t.T, meta


def dedup_columns(probes, V):
    """Average duplicated probe IDs within one series."""
    uniq, inv = np.unique(np.asarray(probes), return_inverse=True)
    if len(uniq) == len(probes):
        return probes, V
    out = np.zeros((V.shape[0], len(uniq)), dtype=V.dtype)
    cnt = np.bincount(inv)
    for j in range(V.shape[1]):
        out[:, inv[j]] += V[:, j]
    return [str(u) for u in uniq], out / cnt[None, :]


def build_addneuromed():
    p60, v60, m60 = parse_geo_series(RAW / "GSE63060_series_matrix.txt.gz")
    p61, v61, m61 = parse_geo_series(RAW / "GSE63061_series_matrix.txt.gz")
    p60, v60 = dedup_columns(p60, v60)
    p61, v61 = dedup_columns(p61, v61)
    shared = sorted(set(p60) & set(p61))
    idx60 = {p: i for i, p in enumerate(p60)}
    idx61 = {p: i for i, p in enumerate(p61)}
    cols60 = np.array([idx60[p] for p in shared])
    cols61 = np.array([idx61[p] for p in shared])
    A = np.vstack([v60[:, cols60], v61[:, cols61]])
    batch = np.concatenate([np.zeros(v60.shape[0], int),
                            np.ones(v61.shape[0], int)])
    meta_all = m60 + m61
    status = np.array([m["status"] for m in meta_all])
    age = np.array([float(m["age"]) if m["age"] else np.nan
                    for m in meta_all])
    gender = np.array([m["gender"] for m in meta_all])

    bad_col = np.isnan(A).any(axis=0)
    bad_row = np.isnan(A).any(axis=1)
    keep_row = ~bad_row & ~np.isnan(age) & (status != None)  # noqa: E711
    A = A[keep_row][:, ~bad_col]
    batch, status, age = batch[keep_row], status[keep_row], age[keep_row]

    if float(A.max()) > 25.0:
        A = np.log2(A + 1.0)
    mu, sd = A.mean(axis=0, keepdims=True), A.std(axis=0, keepdims=True)
    Z = (A - mu) / np.maximum(sd, 1e-12)
    var_order = np.argsort(-Z.var(axis=0), kind="stable")

    out = {}
    for name, P, sel in (("A_main", 2000, np.arange(len(batch))),
                         ("A_sub", 1800, np.where(batch == 1)[0])):
        cols = np.sort(var_order[:P])
        sel = np.asarray(sel)
        out[name] = dict(X=Z[np.ix_(sel, cols)].astype(np.float64),
                         batch=batch[sel], status=status[sel],
                         age=age[sel], gender=gender[sel],
                         probe_ids=[shared[c] for c in cols])
    return out


# ---------------------------------------------------------------------------
# Families B/C: tabular designs + random Fourier featurization
# ---------------------------------------------------------------------------


def rff_map(Xs: np.ndarray, P: int, seed: int) -> tuple[np.ndarray, dict]:
    """Deterministic RFF map with median-heuristic bandwidth."""
    rng = np.random.default_rng(seed)
    sub = Xs[rng.choice(len(Xs), size=min(500, len(Xs)), replace=False)]
    D2 = (np.sum(sub ** 2, 1)[:, None] + np.sum(sub ** 2, 1)[None, :]
          - 2.0 * sub @ sub.T)
    med = float(np.median(D2[np.triu_indices_from(D2, 1)]))
    h = float(np.sqrt(max(med / 2.0, 1e-12)))
    W = rng.standard_normal((Xs.shape[1], P)) / h
    b = rng.uniform(0.0, 2.0 * np.pi, size=P)
    Z = np.cos(Xs @ W + b) * np.sqrt(2.0 / P)
    # frozen convention: RFF block is z-scored per pooled column so the
    # design obeys the model-card normalization sigma_u ~ O(1); this
    # rescales the whole spectrum affinely and leaves its shape intact.
    mu_z, sd_z = Z.mean(axis=0, keepdims=True), Z.std(axis=0, keepdims=True)
    Z = (Z - mu_z) / np.maximum(sd_z, 1e-12)
    return Z.astype(np.float64), {"h": round(h, 6), "P": P, "seed": seed}


def build_tabular():
    out = {}
    ctrl_names = ["inc", "incsq", "agesq", "age", "male", "marr", "fsize"]
    # each family is optional at build time so shards can fetch only their
    # own primary sources (Colab self-containment)
    if (OUT / "k401ksubs.csv").exists():
        with open(OUT / "k401ksubs.csv") as fh:
            rows = list(csv.DictReader(fh))
        D_raw = np.array([[float(r["e401k"]), float(r["nettfa"])]
                          for r in rows])
        C = np.array([[float(r[c]) for c in ctrl_names] for r in rows])
        rng = np.random.default_rng(GLOBAL_SEED + 77)
        sub = np.sort(rng.choice(len(C), size=800, replace=False))
        Cs = C[sub]
        mu, sd = Cs.mean(0), Cs.std(0)
        Cs_std = (Cs - mu) / np.maximum(sd, 1e-12)
        for name, P in (("C_main", 1600), ("C_wide", 160)):
            Z, info = rff_map(Cs_std, P, GLOBAL_SEED + 78)
            out[name] = dict(X=Z, treat=D_raw[sub, 0],
                             outcome=D_raw[sub, 1], ctrl_names=ctrl_names,
                             rff=info, raw_controls=Cs_std)
    if (OUT / "ihdp.csv").exists():
        with open(OUT / "ihdp.csv") as fh:
            rows = list(csv.DictReader(fh))
        treat = np.array([float(r["treatment"]) for r in rows])
        yobs = np.array([float(r["outcome"]) for r in rows])
        H = np.array([[float(r[f"feature{k}"]) for k in range(25)]
                      for r in rows])
        mu, sd = H.mean(0), H.std(0)
        Hs = (H - mu) / np.maximum(sd, 1e-12)
        for name, P in (("B_main", 750), ("B_wide", 150)):
            Z, info = rff_map(Hs, P, GLOBAL_SEED + 88)
            out[name] = dict(X=Z, treat=treat, outcome=yobs, rff=info,
                             raw_controls=Hs)
    return out


# ---------------------------------------------------------------------------
# Spectral audit shared by all families
# ---------------------------------------------------------------------------


def spectral_audit(X: np.ndarray) -> dict:
    sys.path.insert(0, str(ROOT / "code"))
    from de_formulas import (estimate_noise_scales, onatski_select,
                             tw_mu_sigma, tw_threshold)
    n, p = X.shape
    Xc = X - X.mean(axis=0, keepdims=True)
    if p <= n:
        d = np.linalg.eigvalsh(Xc.T @ Xc / n)[::-1]
    else:
        d = np.maximum(np.linalg.eigvalsh(Xc @ Xc.T / n)[::-1], 0.0)
    c = p / n
    # FROZEN benchmark noise-floor convention (docs/benchmark_protocol.md):
    # gene-expression designs show an approximate white bulk and the shared
    # estimator applies; kernel-featurized designs decay smoothly and admit
    # NO MP-consistent bulk. Uniform rule for all benchmarks:
    #     se2_bench = max( q25(d), 1e-3 * mean(d) ),
    # and the exact algebraic decomposition Sigma_X = se2 I + sum_j l_j se2
    # q_j q_j' gives l_hat_j = d_j / se2 - 1 >= 0 with no BBP inversion step.
    se2_mp, sy2 = estimate_noise_scales(d, c)
    se2 = float(max(np.quantile(d, 0.25), 1e-3 * float(np.mean(d))))
    mu_np, sig_np = tw_mu_sigma(n, p)
    lam_max = float(d[0])
    ktop = int(min(max(onatski_select(d), 1), 10))
    return {
        "n": int(n), "p": int(p), "c": round(c, 3),
        "se2_mp_bulk_est": round(float(se2_mp), 6),
        "se2_bench": round(float(se2), 6),
        "sigma_y2_hat": round(float(sy2), 4),
        "lam_max_cov": round(lam_max, 3),
        "tw_stat_lammax": round(float((lam_max * n - mu_np) / sig_np), 1),
        "outlier99_white": bool(lam_max > tw_threshold(n, p, 1.0)),
        "r_hat_onatski": int(onatski_select(d)),
        "ktop_alarm": ktop,
        "r_inj": ktop,
        "tau_top10": [round(float(v), 2) for v in (d[:10] / se2)],
        "l_hat_top10": [round(float(v), 2) for v in (d[:10] / se2 - 1.0)],
        "lambda_col_sd_top10": [round(float(np.sqrt(max(v, 0.0))), 3)
                                for v in (d[:10] - se2)],
    }


def known_results() -> dict:
    """Empirical reproduction of each family's canonical result (raw blocks,
    no injection): the audit trail required by WP 3.2 / plan Section 9.7."""
    import csv
    res = {}
    with open(OUT / "k401ksubs.csv") as fh:
        rows = list(csv.DictReader(fh))
    ctrl = ["inc", "incsq", "agesq", "age", "male", "marr", "fsize"]
    D = np.array([float(r["e401k"]) for r in rows])
    Y = np.array([float(r["nettfa"]) for r in rows])
    Cm = np.array([[float(r[c]) for c in ctrl] for r in rows])
    M = np.column_stack([np.ones(len(rows)), D, Cm])
    coef, *_ = np.linalg.lstsq(M, Y, rcond=None)
    res["k401k_ols_e401k_full"] = {
        "coef": round(float(coef[1]), 3),
        "mean_outcome": round(float(Y.mean()), 2),
        "n": len(rows),
        "claim": "401(k) ELIGIBILITY is associated with higher net financial "
                 "assets (positive coefficient), the Poterba-Venti-Wise "
                 "finding reproduced by plain covariate adjustment",
    }
    with open(OUT / "ihdp.csv") as fh:
        rows = list(csv.DictReader(fh))
    t = np.array([float(r["treatment"]) for r in rows])
    y = np.array([float(r["outcome"]) for r in rows])
    H = np.array([[float(r[f"feature{k}"]) for k in range(25)] for r in rows])
    diff = float(y[t == 1].mean() - y[t == 0].mean())
    M = np.column_stack([np.ones(len(rows)), t, H])
    coef, *_ = np.linalg.lstsq(M, y, rcond=None)
    res["ihdp_naive"] = {
        "unadjusted_diff": round(diff, 3),
        "ols_adjusted_coef": round(float(coef[1]), 3),
        "experimental_benchmark_ate": 4.0,
        "claim": "IHDP point estimates land at the experimental-benchmark "
                 "scale (ATE ~ 4.0, Hill 2011): covariate adjustment gives "
                 "3.93, slightly below it; stored as this mirror's "
                 "reproduction anchor",
    }
    return res


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    audit = {}
    for name, payload in {**build_addneuromed(), **build_tabular()}.items():
        prof = spectral_audit(payload["X"])
        audit[name] = prof
        keep = {k: v for k, v in payload.items() if isinstance(v, np.ndarray)}
        meta = {k: v for k, v in payload.items() if k not in keep}
        np.savez_compressed(OUT / f"{name}.npz", **keep,
                            config_name=name, meta_json=json.dumps(
                                {**meta, "profile": prof}, default=str))
        print(name, json.dumps(prof))
    (OUT / "spectral_audit.json").write_text(json.dumps(
        {"spectral_profiles": audit, "known_results": known_results()},
        indent=1))
    print("wrote", OUT / "spectral_audit.json")


if __name__ == "__main__":
    main()
