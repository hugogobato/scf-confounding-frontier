# Phase 2 Pre-Registration (SCF, gate G3)

Status: FROZEN 2026-08-24 before any Phase 2 sweep execution. Thresholds are
copied from research plan Section 6 (G3 row) and WP 2.2/2.3; grids inherit the
pilot memo amendments. Any deviation after this freeze must appear in the
relevant gate memo's deviation register with reason and timing.

Ledger basis: model card v1.0 + assumption ledger v1.0 (hash recorded per
row). Additive clarifications frozen alongside this document:
docs/detection_statistics.md (statistic family) and model-card-compatible
DGP extensions defined in code/simulator.py Config (conf_kind sparse,
error_law variants, heteroskedastic flag, correlated-factors flag, M2
treatment block). These instantiate ledger-documented robust variants; no
A1-A6 row is weakened.

## Deviation register (pre-data, updated as deviations arise)

D1 (2026-08-24, before any sweep data): F12 null-variance law corrected from
(1 + c d) to (sigma_eps2 + d/c)/sigma_y2; see docs/detection_statistics.md
Erratum 1. Affects S2 calibration and the SEB tuner's g^2 inversion only;
no grid or threshold changes.

D2 (2026-08-24, before any sweep data): gate rejections for detection use
MC-calibrated thresholds from matched g = 0 configs; the analytic S1
threshold is co-recorded but known miscalibrated at finite n
(detection_statistics.md Erratum 2). This was predeclared in the original
freeze ("MC-calibrated variant must pass"); nothing new is unlocked.

D3 (2026-08-24, before any sweep data): estimation metrics consume the RAW
centered response (pilot convention); response standardization happens
inside detection statistics only. Rationale: standardizing Y rescales fitted
coefficients by 1/sigma_y and would contaminate bias ratios; caught by
tests/test_phase2.py::test_eb_spectral_reduces_harmful_bias_vs_ols.

D4: SDBoost linear special case uses a shared-left-singular-subspace
approximation inside its BLUP-corrected CV folds (no fold-wise spectra);
documented in code, cost-motivated, discrimination-equivalent by design.

D5 (2026-08-24, before any sweep data): S1 aug_bbp demoted from primary to
co-recorded diagnostic after pre-data diagnostics showed no discrimination
beyond S0+S2 (detection_statistics.md Erratum 3); gate statistics for WP 2.3
become S2 maxz_cal (MC-calibrated) plus probe/baselines; undetectability of
subcritical dense confounding promoted to a measured headline object with a
T2/T3 theory question attached. No grid changes; falsification rules updated
in the spec file.

D6 (2026-08-24, before any sweep data): COMPUTE-DRIVEN GRID AMENDMENT v2,
invoking the plan's own risk-register recovery action ("reduce grid
adaptively") after cost calibration against WP 1.5 timings showed the
originally frozen grids at 1.5-2x the 40-notebook budget. Changes, all
pre-data:
  (a) c = 10 dropped from all grids (c = 5 retains the p >> n regime;
      c = 10 cells were the most expensive per bit of information);
  (b) correctness: reps 1000 -> 350 at n = 2000 (main), alignment slice
      restricted to r = 5 profiles with reps 200, n = 8000 tier restricted
      to the skinny-c branch c in {0.1, 0.2} (spectrum rides the p-side
      Gram there; at c >= 0.8 an 8000-row cell costs minutes PER REP);
  (c) null calibration: core reps 10000 -> 2500 (size SE 0.0044, ample for
      the [0.035, 0.065] gate), theta-matched nulls 800, n-ladder keeps
      n = 500 at 4000 reps and adds ONE n = 4000 trend cell at 400 reps
      (trend-only, explicitly NOT size-gating); n = 8000 dropped from
      nullcal (minutes-per-rep spectra x 10k reps is infeasible);
  (d) power surface: theta set {pi/6, pi/2}, g-grid 7 points
      {0.125, 0.3, 0.6, 1.0, 1.6, 2.5, 4.0}, reps 450 (power SE 0.023);
      alignment sweep reps 5000 -> 1500;
  (e) estimation: reps 300 (n = 500) / 250 (n = 2000) with c = 5 column at
      150; crossover reps 300; robustness reps 250 on c in {0.2, 2.0};
      lava_default estimator dropped from the roster (redundant with
      cevid_default for the transform+OLS class);
  (f) M2 block reps 250.
Projected total with the calibrated cost model (code/make_shards.py prints
it): ~480-520 core-hours across <= 36 shards of <= 5.5 h target. Residual
overshoot risk is absorbed by the resume-safe design (any interrupted cell
reruns from zero loss; shards stop cleanly before the Colab wall limit).

D8 (2026-08-24, at consolidation, before any analysis): grid-file config_id
columns were computed by a pre-freeze Config.cid implementation; all joins
recompute authoritative ids from config dicts under the pinned tag (the ids
that seed every rep). configs/cid_remap.json is the audit mapping; no config
dict changed and no data was touched.

D9 (2026-08-24, at analysis): plan wording "n = 4000-equivalent" interpreted
as the largest full-coverage tier (n = 2000), n = 8000 skinny-c tier reported
alongside.

D10 (2026-08-24, at analysis): null-size gate evaluated noise-aware (pooled
chi-square on standardized split-half deviations) alongside the frozen raw
arithmetic, because D6c rep cuts made the raw band comparable to MC error.

D11 (2026-08-24, at analysis): subspace-overlap functionals (WP 2.1 overlay
iv) not stored by the runner schema; recorded scope cut.

D12 (2026-08-24, at analysis): Le Cam declaration "below the claimed
frontier" operationalized as S-blind strata (infinite predicted frontier);
per-cell labels co-recorded.

D13 (2026-08-24, at analysis): probe feature map truncates log-eigenvalues
at 10 plus 4 scalars; r-dependent zero padding omitted (no gate depends on
truncated dimensions).

## Inherited amendments (from docs/pilot_memo.md)

1. c = 1.0 is EXCLUDED from all mean-bias claims (E[(X'X)^{-1}] diverges);
   c = 1.0 remains eligible for eigenstructure-only and detection-only cells.
   Mean-bias grids use c in {0.1, 0.2, 0.5, 0.8} on the p < n branch and
   c in {2, 5, 10} on the p > n branch.
2. Loading directions Q fixed per config (q_fixed=True everywhere); bias
   functionals are loading-conditional (model card Section 4).
3. Twin gamma=0 arms on common seeds run for every c > 1 estimation cell and
   every tiny-signal cell; rel_bias_conf is the primary bias metric there.
4. Equal-spike profiles (sub, super) report SUBSPACE overlaps
   (sum of top-j squared overlaps), not per-vector xi_j.
5. Capture law cap_j = (1 + l_j)/(c + l_j) used for c > 1 OLS overlays
   (pilot-validated conjecture; its derivation stays open as T1 input).

## Memory ceiling (binding, extends plan Section 10.1 rule 3)

No cell may have n*p > 2.4e8 entries (X at float64 <= 1.92 GB) or
min(n, p)^3 eigendecomposition beyond 8000^3. Consequence: n = 8000 rows are
restricted to c in {0.1, 0.2, 0.5, 0.8, 2} (p <= 16000). Cells violating the
ceiling are dropped from grids and listed as dropped-by-ceiling in the
correctness memo.

## WP 2.1 Correctness sweep (frozen grid)

Estimation correctness cells (DE overlays):
- c in {0.1, 0.2, 0.5, 0.8, 2.0, 5.0, 10.0}
- (r, profile) in {(1, sub), (1, super), (5, sub), (5, mixed), (5, super),
  (25, mixed)} [six spike profiles per plan; r=25 all-sub/all-super dropped:
  r=25 equal spikes at l = 3 sqrt(c) make Sigma nearly rank-deficient
  numerically at c < 1 and duplicate the mixed cell's information; documented
  economy]
- theta = pi/6 for all cells; theta in {0, pi/2} additionally at n = 2000
  (alignment sensitivity of the overlay)
- n in {500, 2000, 8000 subject to ceiling}: full 7c x 6 profiles at n = 500
  and n = 2000; at n = 8000 only c in {0.1, 0.2, 0.5, 0.8, 2.0}
- g = 1, Gaussian, q_fixed, twins for c > 1, 1000 reps
- Cell counts: 7*6*2 + 7*6 (theta extra at n=2000) + 5*6 = 84 + 42 + 30 = 156

Metrics per cell: relative deviation |simulated - predicted| / predicted for
(i) OLS mean-bias norm (capture law at c > 1, exact identity at c < 1),
(ii) ridge curves over the frozen lambda grid LAM = (0.01, 0.03, 0.1, 0.3, 1,
3, 10) (provisional lambda-interpolated capture at c > 1, flagged open),
(iii) lambda_max vs BBP/bulk edge, (iv) subspace overlaps vs BGN sums.

Null calibration subset (detection statistics, stats-only reps):
- 20 null configs: c in {0.1, 0.2, 0.5, 0.8, 1.0, 2.0, 5.0} x profile in
  {sub, mixed} x n in {2000} plus c in {0.2, 2.0} x n in {500, 8000-at-c<=2};
  gamma = 0; statistics S0-S3 + B1 per docs/detection_statistics.md;
  10,000 reps each (size SE ~ 0.002).

Pass/fail (plan G3 correctness layer, verbatim thresholds):
- DE deviation <= 10% at n = 4000-equivalent in >= 90% of estimation cells,
  shrinking with n;
- null size in [0.035, 0.065] in >= 95% of null cells (MC-calibrated variant
  must pass; analytic-threshold size reported honestly alongside);
- fail: systematic deviation > 25% in >= 30% of cells after two fix
  iterations -> Phase 2 give-up rule 1.

## WP 2.2 Estimation contribution sweep (frozen)

Estimator roster (all tuned with logged budgets; common seeds per rep):
- ols (min-norm at c > 1)
- ridge_cv: 5-fold CV over LAM grid above
- pca_onatski, pca_baing: hard trim, k by Onatski ratio / Bai-Ng selector
- pca_oracle_r: hard trim at true r (diagnostic upper bound of trims)
- cevid_default: Trim transform d~_i = min(d_i, median(d)) then OLS on
  transformed data (linear special case; Cevid et al. 2020 eq. 3.3 with
  tau = median singular value). For A4b sparse-beta cells only: cevid_lasso =
  same transform + Lasso (lambda by 40-point path, EBIC-free 5-fold CV).
- sdboost_linear_eb: Nava et al. linear special case - LAVA-type spectral
  loss w_i = sigma_e^2/(sigma_r^2 d_i^2 + sigma_e^2) with variance components
  tuned by Gaussian marginal likelihood (their eq. for ell(theta)), estimator
  = spectral-loss least squares (m -> infinity limit of their boosting
  recursion with linear base learners; equivalence documented in code)
- eb_spectral (OURS): soft-trim weights w_tau(d) = min(1, tau/d); tau chosen
  by the SEB tuner = minimization of DE-predicted causal MSE
  sum_j [w_tau(d_j) capture_j sqrt(l_j)/(1+l_j) ghat_j]^2 + variance terms
  (var_j = sigma_eps_hat^2 w_tau(d_j)^2 / (n d_j) + rowspace artifact at
  c > 1), with (l_j, ghat_j^2) estimated from the spectrum by BBP inversion
  and cross-moment mixture shrinkage; details in estimators.eb_spectral
- eb_oracle_tau: same family, tau minimizing true known-gamma predicted MSE
  (ablation: tuner-quality upper bound)
- eb_cv_tau: same family, tau by 5-fold prediction-CV (ablation: no-EB)
- oracle_gamma: beta_hat_ols minus true-direction bias removal (upper bound)

Grids:
- Harmful region: profiles {sub, mixed}, c in {0.2, 0.5, 0.8, 2.0, 5.0},
  r in {1, 5, 25}, theta in {0, pi/6, pi/2}, n in {500, 2000}, 1000 reps,
  g = 1. (n = 8000 slice: c in {0.2, 0.8, 2.0}, theta = pi/6, r = 5 only,
  400 reps - cost control, predeclared.)
- Crossover strip (figures): c in {0.2, 0.8, 2.0} x g in {0.25, 0.5, 1, 2, 4}
  x theta in {pi/6}, r = 3 profile mixed, n = 2000, 1000 reps.
- Rung 4 baseline-favorable: sparse-confounder DGP (conf_kind="sparse",
  ||b||_0 = ceil(sqrt(p)) nonzeros in b, same ||b|| as dense twin), aligned
  beta (beta_kind="aligned": unit vector with cosine 0.9 to top factor
  direction), c in {0.2, 0.8}, n = 2000, 500 reps. Expectation: PCA-k and
  OLS win; eb_spectral must not claim otherwise (its honest-regime signal is
  checked: alarm statistic reported but method output not penalized).
- Rung 3 null: gamma = 0 arms of the harmful-region cells (twins) double as
  the no-advantage check.

Metrics (predeclared priority): rel_bias_conf (primary), rel_bias total,
MSE, runtime_s. Coverage deferred to Phase 3 package work (variance formula
not yet derived; noted as scope cut, not a silent omission).

Gate thresholds (plan Section 6, binding):
- No-regret: bias(eb_spectral) <= 1.05 x best-baseline bias in >= 95% of
  harmful cells. Baseline set: {ridge_cv, pca_onatski, pca_baing,
  cevid_default, sdboost_linear_eb}.
- >= 50% bias reduction vs OLS in >= 70% of harmful cells.
- Practical-effect threshold: >= 15% relative MSE reduction vs best baseline
  in the mixed/subcritical region.
- Attribution: EB-vs-CV ablation gap largest near the frontier cells;
  documented expectation, checked descriptively.

Fail rules: plan Phase 2 give-up rules 2-3 verbatim (loses to strongest
baseline in majority of harmful cells -> KILL; gains < 10% bias reduction
everywhere -> INCREMENTAL-ONLY).

## WP 2.3 Detection contribution sweep (frozen)

Power surface: s-grid g in {0, 0.125, 0.25, 0.5, 0.75, 1, 1.5, 2, 3}
x (c, profile) in ({0.2, 0.8, 2.0} x {sub, mixed, super}) x theta in
{0, pi/6, pi/2}, r = 3, n = 2000; 1000 reps per alternative cell;
27 null configs at 10,000 reps (shared with WP 2.1 null subset where
overlapping; union deduplicated by config id).
Statistics: S1 aug_bbp (analytic + MC threshold), S2 maxz_cal, S0 scree,
B1 f_test_pcs; B2 ucm_strength_boot on the headline subset only (c = 0.8
mixed/sub, all g, theta = pi/6; 200 bootstrap replicas).

Frontier comparison: empirical power-1/2 contour per (c, profile, theta)
vs F13/F12 predictions; pass if within factor 1.5 for S1-analytic and
S2 (plan: power >= 0.8 at s <= 1.5 x predicted frontier; size in
[0.035, 0.065] in >= 95% of null cells).

Alignment stress: 12-point theta sweep on (c = 0.8, mixed, r = 3, n = 2000,
g = 1), 5000 reps per point, S1/S2/B1 power + H0 size-invariance check.

Le Cam probe: pooled discrimination on (c in {0.2, 0.8}) x profile in
{sub, mixed} x g-grid; 2000 + 2000 datasets per cell; GBM + MMD AUC;
declaration threshold AUC <= 0.55 both probes below the claimed frontier.

Fail rules: plan give-up rule 4 (uncalibratable size) and the falsification
list in docs/detection_statistics.md.

## WP 2.4 Robustness and scaling (frozen)

Core reference block: c in {0.2, 0.8, 2.0} x profile in {sub, mixed} x r = 5
x theta = pi/6 x n = 2000 x 500 reps under each variant:
V0 gaussian (reference), V1 t5 errors (u_ij, eps scaled t_5),
V2 bernoulli-loadings (entries +/- sqrt(p)/sqrt(nnz) with nnz = p/2 per
column, preserving Lambda Lambda' scale approximately),
V3 heteroskedastic-u (row variances (1 + chi1^2)/2, i.e., E = 1),
V4 correlated factors (f ~ N(0, Omega), Omega random SPD with eigenvalues
uniform[0.5, 1.5]; breaks Var(f) = I normalization deliberately),
V5 r-misspecification (methods receive r_hat +/- 1 where they consume r),
V6 sparse-confounder (as rung 4 above).
Pass: qualitative phase-diagram structure survives; diagnostic size < 0.15
under V3-V4 or a documented robust variant fixes it.
Scaling study: timing/memory of spectrum + estimator suite at
(n, p) in {(1000, 1000), (2000, 8000), (4000, 8000), (8000, 8000),
(8000, 16000)}, single rep, Colab-only for the last two; deliverable =
compute-envelope table for the package.

## Compute and sharding

Cost calibration: micro-pilot on local workstation measures per-rep cost
model t(n, p, mode) for mode in {estimation, stats-only}; shards assembled to
target <= 5 h wall each on 2 Colab vCPU (weaker than local cores; safety
factor 2.5x applied to local timings). Budget: <= 36 notebooks total,
leaving headroom within the plan's 40-notebook ceiling for reruns. Local
execution limited to: null-calibration stats-only cells, cost pilots, and
anything measured <= 2 h wall at <= 6 workers (plan Section 10.2 rule).

Consolidation: shard outputs land in data/sim/<sweep>/shards/ with per-shard
manifests (sha256 of each parquet); code/consolidate_shards.py merges,
verifies checksums, and emits the sweep-level parquet + completeness report
(all config ids x reps present).

## Analysis and figure regeneration

Every figure script reads parquet only: figures/de_overlay_grid.pdf,
power_surface_vs_frontier.pdf, lecam_probe_auc.pdf, alignment_stress.pdf,
bias_phase_diagram_empirical.pdf, crossover_curves.pdf. Gate verdicts come
from code/gate_checks.py reading the consolidated parquets and printing the
Section 6 numeric conditions verbatim.
