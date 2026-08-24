# Robustness and Scaling Memo (WP 2.4) - SKELETON, data pending

Status: pre-registered expectations frozen 2026-08-24 (before sweep launch).

## Variant grid (D6e/D7)

V0 gaussian reference; V1 t5 errors; V2 rademacher_half loadings (half-
support, columns rescaled to exact Lambda'Lambda = diag(sigma_u^2 l_j));
V3 row-heteroskedastic u (variance factor (1+chi1^2)/2, mean 1); V4
correlated factors (Omega random SPD on [0.5, 1.5]); V5 r-misspecification
(+/-1 consumed by the SEB spike-profile prior); V6 sparse-confounding.
Cells: c in {0.2, 2.0} x {sub, mixed} x r = 5, theta = pi/6, n = 2000,
125 reps per variant. M2 block: weak-treatment alignment delta_g = 0.3,
tau = 1, c in {0.2, 0.8, 2}, 150 reps.

## Preregistered pass rule

Qualitative phase-diagram structure survives all perturbations; diagnostic
size < 0.15 under V3-V4 or a documented robust variant fixes it.
Fail: advantage/calibration is a Gaussian-homoskedastic artifact only ->
restrict claims to that regime explicitly (PIVOT within G3).

## Expected sensitivities (pre-data)

- V4 correlated factors break Var(f) = I: bias functionals shift by the
  Omega-scale; the capture-law overlay is expected to degrade gracefully;
  S2 size may inflate (whitening assumes isotropic factors).
- V2 breaks Haar-isotropy of rows: overlap geometry changes; SEB tuner's
  l_hat extraction degrades first at equal-spike profiles.
- M2: tau_trim estimators should dominate tau_ols under dense confounding;
  absolute tau errors are expected O(0.05-0.4) depending on c.

## Scaling study

Timing/memory envelope at (n, p) in {(1000,1000), (2000,8000), (4000,8000),
(8000,1600)} x 2 reps -> package compute-envelope table.

## Results

PENDING.

## Verdict

PENDING.
