# Estimation Gate Memo (WP 2.2) - SKELETON, data pending

Status: pre-registered expectations frozen 2026-08-24 (before sweep launch).
Results fill from data/sim/estimation/ + crossover/.

## Frozen estimator roster (12; lava_default dropped by D6e)

ols, ridge_cv (exact-LOO over LAM grid), pca_onatski, pca_baing,
pca_oracle_r, cevid_default (Trim transform, tau = median singular value,
OLS after transform - linear special case of Cevid et al. 2020 eq. 3.3),
sdboost_linear_eb (Nava et al.: EB variance components on the XX' spectrum,
LAVA weights, boosting path, BLUP-corrected CV stopping; shared-subspace
fold approximation D4), eb_spectral (OURS: SEB tuner on soft-trim family
with F12-Erratum-1-consistent gamma^2 inversion), eb_cv_tau (ablation:
prediction-LOO tuning), eb_oracle_tau (ablation: true-(l,gamma^2) objective),
oracle_gamma (upper bound).

## Preregistered gate conditions (plan Section 6 G3 + preregistration)

1. No-regret: bias(eb_spectral) <= 1.05 x best baseline in >= 95% of harmful
   cells (baseline set above). Primary bias functional at c > 1: twin
   difference (gamma-attributed); at c <= 1: mean-diff norm.
2. >= 50% bias reduction vs OLS in >= 70% of harmful cells.
3. Practical effect: >= 15% relative MSE reduction vs best baseline in the
   mixed/subcritical region (secondary metric, reported).
4. Attribution expectation: EB-vs-CV ablation gap largest near frontier cells;
   perturbation-free SEB relies on dense-beta geometry (A4a) - rung 4 cells
   must show PCA-k/OLS winning with eb_spectral not claiming otherwise.

## Give-up rules armed (binding)

- KILL if adaptive worse than strongest baseline in majority of harmful cells.
- INCREMENTAL-ONLY if gains < 10% bias reduction everywhere.

## Results

PENDING - populate from consolidated parquet via code/gate_checks.py
(check_estimation implements conditions 1-2 mechanically).

## Verdict

PENDING.
