# confounderalarm

A calibrated spectral alarm for hidden **dense** confounding in
high-dimensional regression, with a predicted detectability frontier and a
trim-then-regress adjustment for treatment effects.

Given observational data `(Y, X)` (optionally a treatment block `D`), the
package:

1. estimates the design's factor profile `(r_hat, l_j, c = p/n)`;
2. runs a spike-coordinate alarm whose null calibration comes from
   permutation on the data itself (marginals preserved);
3. predicts the benchmark-specific detection frontier `g*` from the SCF F12
   law (BBP/BGN deterministic equivalents) using the same calibration;
4. reports a verdict with an honest **blind-region certificate**: where the
   theory says non-detection is uninformative, the package says so instead
   of declaring safety;
5. when `D` is supplied, recommends the Onatski hard-trim adjustment
   (trim-then-regress). Tuned soft weights are explicitly NOT recommended
   (they were executed by the Phase-2 simulation gates of the underlying
   research project).

## Identification caveat (read this)

The alarm tests for confounding **relative to the A4a near-orthogonality
ledger**: the structural coefficient vector is treated as generic/isotropic
relative to the top design eigenspace, and "confounding" means a dense
factor link `gamma'f` sharing the latent factors of `X`. Sparse confounding
is out of scope and cannot be certified either way. Permutation p-values
calibrate against the sharper hypothesis "no second-moment association
beyond marginals"; the ledger-relative interpretation is the package's
documented semantics.

## Install

```bash
pip install -e .
```

## Quickstart (Python)

```python
from confounderalarm import fit_alarm

rep = fit_alarm(Y, X, n_perm=400)
print(rep.alarm, rep.p_value, rep.placement)
print(rep.certificate)
print("r_hat =", rep.r_hat, " c =", round(rep.c, 3),
      " g* =", rep.g_star)

# with a treatment block: adds trim-then-regress
rep = fit_alarm(Y, X, D=D)
print(rep.adjustment)   # {"method": "onatski_trim", "k": ..., "tau_trim": ..., ...}
```

## CLI

```bash
python -m confounderalarm --csv mydata.csv --y outcome --treatment treat \
    --exclude id,postcode --n-perm 400
```

Prints the full verdict as JSON.

## Worked example 1: batch structure in gene expression

Real geometry template: AddNeuroMed blood expression batches (GSE63060/61).
Batch effects dominate the spectrum (`TW statistic of lambda_max ~ 4.8e4`
against white noise), which is exactly why "I checked the scree plot" is not
evidence about confounding: the alarm separates *design spikes* (always
there) from a *response link* (the thing being tested).

```python
import numpy as np
from confounderalarm import fit_alarm

z = np.load("A_main.npz", allow_pickle=True)   # any n x p expression matrix
X = z["X"]
rng = np.random.default_rng(0)
beta = rng.standard_normal(X.shape[1]); beta /= np.linalg.norm(beta)
Y = X @ beta + rng.standard_normal(X.shape[0])          # no confounding
rep = fit_alarm(Y, X, n_perm=300)
# -> alarm False; outlier99_white True (scree WOULD have "rejected")
print(rep.outlier99_white, rep.alarm, rep.placement)
```

## Worked example 2: treatment effect under injected dense confounding

```python
import numpy as np
from confounderalarm import fit_alarm

rng = np.random.default_rng(5)
n, p, r = 300, 600, 2
Q, _ = np.linalg.qr(rng.standard_normal((p, r)))
Lam = Q * np.sqrt(np.array([16.0, 4.0]))[None, :]
f = rng.standard_normal((n, r))
X = f @ Lam.T + rng.standard_normal((n, p))
pi = np.zeros(p); pi[:6] = 1/np.sqrt(6)
delta = rng.standard_normal(r); delta *= 0.3 / np.linalg.norm(delta)
D = X @ pi + f @ delta + rng.standard_normal(n)
beta = rng.standard_normal(p); beta /= np.linalg.norm(beta)
gam = 2.0 * np.ones(r) / np.sqrt(r)
tau_true = 1.0
Y = tau_true * D + X @ beta + f @ gam + rng.standard_normal(n)

rep = fit_alarm(Y, X, D=D, n_perm=300)
rep.adjustment["tau_trim"]   # ~ tau_true (Onatski hard trim, k selected)
rep.adjustment["tau_ols"]    # materially biased in this regime (c > 1)
```

## Interpretation of report fields

| field | meaning |
|-------|---------|
| `alarm`, `p_value` | verdict at level `alpha`; permutation-calibrated |
| `statistic`, `threshold` | max standardized spike coordinate vs its null q95 |
| `r_hat`, `l_hat`, `c` | estimated factor count, normalized spike strengths, aspect ratio |
| `g_star` | smallest link strength the F12 law predicts detectable with power 0.8 (None = no visible channel) |
| `placement` | one of `above` / `near` / `below` / `blind` relative to the frontier |
| `certificate` | plain-language statement about what non-detection DOES or DOES NOT mean here |
| `ucm_rho` (+ `ucm_p`) | UCM-style strength proxy (Rendsburg-et-al.-spirit, approximate) |
| `js_asym` (+ `js_p`) | Janzing-Schoelkopf-style spectral asymmetry (approximate adaptation) |
| `adjustment` | Onatski hard-trim treatment coefficient when `D` given |

## What this package deliberately does not do

* sparse confounding (out of scope by design; see the ledger);
* claims that non-detection implies absence of harmful confounding below
  the frontier — the certificate states the opposite;
* tuned soft-spectral weights as an adjustment (dominated by hard trimming
  in the underlying research program's preregistered simulations).

## Development

```bash
pip install -e ".[dev]"
pytest
```

Source project: SCF (Spectral Confounding Frontier), Phase 3 artifact.
Frozen spec: `configs/benchmarks_frozen.yaml`,
`docs/benchmark_protocol.md` in the research repository.
