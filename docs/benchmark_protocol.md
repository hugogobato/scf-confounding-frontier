# Benchmark Protocol (WP 3.2, Phase 3)

Status: FROZEN 2026-08-24 before any comparative result. Config source of
truth: `configs/benchmarks_frozen.yaml` (v1); this doc is its readable
companion. Any post-data change requires a numbered deviation entry here
(D-B1, ...).

## What is being tested

The calibrated confounding alarm and its benchmark-specific predicted
frontier g*, on real covariate geometry where ground truth is created by
injection (semi-synthetic design in the template of Ulmer et al. 2025):
real X, known beta/gamma/f mechanism, known tau. Predeclared findings:

* **PF-1 (primary):** the gate alarm straddles each main family's own
  predicted frontier: size-calibrated under negative controls (rejection in
  [0.02, 0.10]) with power >= 0.8 at 2x g* and <= 0.25 at 0.5x g*.
* **PF-2:** naive spectrum-gazing false-alarms: TW99 scree inspection
  rejects "white noise" on every unmodified real design before any
  confounding is injected (audit table already quantifies this).
* **PF-3:** trim-then-regress (Onatski k, the post-G3 recipe) recovers
  tau_true = 1 sign/magnitude under 1x-frontier injected dense confounding
  on C_main where raw-OLS tau error inflates materially.
* **PF-4 (secondary):** the response-aware UCM strength proxy tracks
  injected g monotonically while its permutation threshold is invalid as a
  calibrated test on real geometry (informative ordering, no calibration).

## Designs and injection

Six processed designs (three families x main/sensitivity cuts; full audit in
`docs/data_audit.md`): A_main/A_sub (AddNeuroMed blood expression, geo batch
structure), B_main/B_wide (IHDP, random-Fourier featurized), C_main/C_wide
(k401ksubs, RFF featurized). Injection follows simulator.gen_data semantics:

    X_obs = Xc_base + f Lam',   Y = X_obs beta + f gam + eps,
    beta Haar ||beta|| = 1, f ~ N(0, I_r), eps ~ N(0, 1),
    Lambda columns = first r_inj sample eigenvectors, column sd sqrt(d_j - se2),

so Cov(X_obs, Y) = Sigma_obs beta + Lambda gamma exactly as in M1 and all
Phase-1/2 functionals apply with the observed second moment. Matched-null
twins share (beta, f, eps) seeds and drop only the f@gam link from Y (the
design keeps its realized factor structure, Phase-2 twin convention). M2
block (C_main): D = X pi + f delta + nu with sparse pi, ||delta|| = 0.3,
tau_true = 1.

Noise-floor convention F-BENCH: se2 = max(q25(d), 1e-3 mean(d));
l_j = d_j/se2 - 1 (exact decomposition, no BBP inversion); r_inj =
ktop_alarm = clip(Onatski, 1, 10). Default gamma direction: equal mass over
the r_inj coordinates.

## Gate alarm (deviation D-B0, amended BEFORE any calibration data)

Smoke testing exposed two failures of the frozen Phase-2 S2 statistic on
real geometry, both pre-data: (1) gamma-blind response-scale suppression;
(2) MP-white-bulk misestimation on family A's smooth spectrum. Final frozen
form ("S2-bench"):

    zeta_j = sqrt(n) v_j' (Xc'Yc/n)_j / sqrt(d_j)     (raw centered Y)
    T      = max_{j < ktop} |zeta_j| / s_j,
    s_j    = twin-estimated per-coordinate null scales (pass 1, then frozen),
    reject iff T > mc95(T | matched-null pool)              (frozen).

The H1 mean shift v_j'Lambda gamma enters zeta in absolute units and is
immune to response-scale dilution; the twin scales absorb every bulk-shape
factor empirically. The F12 law survives as the PREDICTION layer: slopes

    slope_j = dir_j omega_j sqrt(n se2 l_j / d_pred_j) / s_j,
    d_pred_j = bbp_location(l_j, c, se2), omega_j = capture weight (c > 1),

feed the same Gaussian-max MC construction as
`phase2_analysis.predicted_frontier_g`, giving g* = smallest g with
predicted power >= 0.8. PF-1 tests exactly this prediction layer against
the empirical curve.

## Two-pass order (binding)

1. **Pass 1 (calibration only):** arms `null` (600 reps/config) and
   `perm_null` (400 perms/config). Produces results/benchmark_freeze.json:
   per-config coord scales s_j, mc95, permutation q95 thresholds for the two
   APPROXIMATE baselines (ucm_rho, js_asym), analytic sizes, and g*. No
   pass-2 arm may be generated before this file exists.
2. **Pass 2 (evaluation):** positive arms at {0.5x, 1x, 2x} g*, negative
   controls (permuted-Y evaluation rows are part of pass 1's threshold pool;
   split-half placebo; batch-free subset = A_sub), sensitivity (alignment
   top/weak, r_inj +/- 1, heteroskedastic eps), M2 block (C_main), evaluated
   against the frozen values.

## Head-to-head methods (WP 3.3 deliverable: calibration-vs-informality)

On identical data and controls: our S2-bench alarm (twin-calibrated);
UCM-strength proxy rho_hat (response-aware, Rendsburg-et-al.-spirit,
APPROXIMATE flag; permutation q95 threshold); Janzing-Schoelkopf-style
spectral asymmetry score js_asym (rank-one response-removal eigenvalue
drops, APPROXIMATE transcription of janzing2018detecting; permutation q95);
S0/TW99 scree (analytic white threshold); partial-F on PC scores (analytic
F95). Reported per method: realized rejection under every negative control
and power at matched alternatives.

## Fair-comparison and honesty rules (inherited from plan Section 8.3)

Common seeds across arms within a config; tuning budgets fixed (no tuning
for any detection method); failed/degenerate reps preserved; primary metric
predeclared (rejection indicators and T statistics; tau absolute errors for
the M2 block); MC uncertainty as 95% intervals across reps; anomalies
investigated before being called findings; contradictions with trusted
benchmark results treated as bugs until proven otherwise.

## Deviation register

* **D-B0** (2026-08-24, pre-pass-1): alarm statistic replaced by the
  S2-bench form above after smoke testing; full rationale in the YAML
  deviations block. Two intermediate candidates were implemented and
  rejected during the same smoke loop (SEB-inversion rescaling: inherits the
  blindness it tries to fix; realized-sd var_cal: fixes dilution but not
  the family-A bulk-shape mismatch). Phase 2 code paths untouched; Phase 2
  conclusions unaffected (its regimes had moderate spikes and weak links,
  where the frozen statistic's frontier ratios were measured, not assumed).
