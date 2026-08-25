# Detection Statistics Specification (SCF Phase 2)

Status: FROZEN 2026-08-24, before any WP 2.1-2.3 sweep data was generated.
ERRATUM 4 (2026-08-25, session 4; full forensics in
docs/theory_T3a_eigenvalue_contiguity.md Section 5): the S1 statistic t_aug
in code/detection.py silently mis-brackets its secular bisection whenever
disc = (1+lo)^2 - 4(lo + s_total) < 0, which is COMMON at c <= 1 designs:
the sqrt(max(disc,1e-300)) fallback collapses hi below d[0], the bisection
runs on an inverted interval, and the returned value is a smooth surrogate
correlated with lambda_max(Sigmahat) rather than lambda_max(M_aug). All
stored t_aug / pow_S1_cal artifacts read as "surrogate-S1". Gate verdicts
are unaffected (S1 demoted pre-data by D5/E3; MC thresholds self-consistent
under H0/H1). The TRUE augmented alarm behaves very differently - it
consistently detects subcritical confounding (direction depends on
geometry); see the T3(a) document. Any future S1 usage must adopt the
validated robust bracketing (tests/test_theory_T3a.py::_aug_stats).
ERRATUM 1 (2026-08-24, same day, during unit testing, BEFORE any sweep data):
the original draft of this document stated Var(z_j) = sigma_eps2 (1 + c d_j);
the correct law is Var(z_j) = (sigma_eps2 + d_j/c)/sigma_y2^0... precisely

    Var(z_j) = [sigma_eps2 + d_j / c] / sigma_y2,

with sigma_y2 = tr(Sigma)/p + sigma_eps2 under A4a (Y is standardized). The
draft inverted the rowspace factor: n d_j/p = d_j/c, not c d_j. Caught by
tests/test_phase2.py::test_maxz_null_variance; measured mean z_1^2 of 5.45
against the wrong prediction 1.94 and the correct prediction ~5.4 at
(n, p) = (900, 300). Consequence recorded for the frontier interpretation:
dense-beta leakage through spike coordinates GROWS like d_j/c, so calibrated
spike-coordinate tests pay a heavy tax at small c (large p); the operational
detection frontier is statistic-dependent, which strengthens the decoupling
story rather than weakening it.

ERRATUM 2 (2026-08-24, BEFORE any sweep data): the v1 analytic threshold for
S1 (plug-in secular root + white-Johnstone TW95 width) is miscalibrated at
finite n: measured sd(t_aug - lam0_plugin)/width ~ 8.5 with a negative mean
offset at n = 600, p = 240 (200 null reps), because plug-in root noise
dominates the genuine TW fluctuation. Per the frozen falsification rules the
analytic variant is REPORTED as-is; gate decisions use the MC-calibrated
thresholds from matched g = 0 configs (mc_thresholds in detection.py).
Deriving the correct finite-n width for the deformed augmented root is a T2
work item.

ERRATUM 3 (2026-08-24, BEFORE any sweep data): pre-data diagnostics (60-rep
cells at n=600, p=240) show (i) S1 aug_bbp carries NO discrimination beyond
S0+S2 in any probed cell - its realized root tracks the design spike
(correlation with lambda_max of Sigmahat near 1) and its H1-vs-H0 shifts are
within null noise; plausibly because the statistic is self-normalizing:
sigma-growing confounding inflates sigma_y^2, which deflates the standardized
couplings and the plug-in root adapts; (ii) in SUBCRITICAL cells at g = 2,
BOTH S1 and S2 are blind (t_maxz 1.196 -> 1.243, t_aug unchanged) while
||b||^2 can even SHRINK (0.913 -> 0.787) because sigma_y grows faster than
the cross-moment mean - direct finite-n evidence for the invisible-yet-
harmful region being genuinely undetectable from (Y, X) second moments by
this whole statistic class, consistent with the OMH contiguity reading and
with the identification subtlety (plan Section 3.2). CONSEQUENCES, frozen:
(a) S1 aug_bbp is DEMOTED from primary to co-recorded diagnostic; the gate
statistics are S2 maxz_cal (supercritical-aligned detector, MC-calibrated)
plus the probe/baselines; (b) the WP 2.3 pass rule "frontier within factor
1.5" applies to the S2 MC-calibrated power-1/2 contour against the F12-law
prediction WITH the d/c leakage tax; (c) the undetectability claim now
carries numerical evidence at the pilot scale and is a headline object, to
be quantified by the Le Cam probe across the s-grid; (d) deriving whether
ANY second-moment statistic can detect subcritical dense confounding (or
whether sigma-normalization makes it impossible) is promoted to a T2/T3
theory question with real content.

This document is the single source of truth for every detection statistic used in
Phase 2; it extends ledger Section "Detection problem under A4a" and de
formula sheet items F9/F12/F13. Any change after data generation requires a
deviation entry in the relevant gate memo.

Conventions follow docs/model_card.md: sigma_u = 1 default; population
eigenvalues tau_j = 1 + l_j of Sigma_X; c = p/n; sample covariance Sigmahat =
X_c' X_c / n with column-centered X_c; Y standardized to variance 1 per
dataset: Ytilde = (Y - mean(Y)) / sd(Y). Write b := X_c' Ytilde / n (the
sample cross-moment vector) and w_j := v_j' b for sample eigenvectors v_j,
d_j eigenvalues of Sigmahat. z_j := sqrt(n) w_j / sqrt(d_j).

Under H0 (gamma = 0), A4a (Haar beta independent of everything), A2/A3:

    E[w_j] = d_j (v_j' beta) ~ N(0, d_j^2 / p),
    Var contribution of eps/f noise: sigma_eps^2 d_j / n,

so

    z_j -> N(0, 1 + c d_j)   approximately at finite n.   (F12)

The c d_j term is the beta-mass leakage amplified by the spike: it vanishes
only when c d_j << 1 or n/p fixed and d_j bounded. This is a REAL null
inflation for naive max-z tests on spiked coordinates; all max-z statistics
below divide it out using estimated d_j and c = p/n. Unit test
test_maxz_null_variance checks this to Monte-Carlo tolerance.

## Statistic family (frozen)

S0 scree_tw99 (practitioner baseline): lambda_max(X'X/n)/sigma_u^2 vs TW99
white-noise threshold tw_threshold(n, p, sigma_u). This is the pilot's
"visibility" statistic; it ignores Y entirely.

S1 aug_bbp (primary; ours): largest root of the augmented second-moment
matrix of (Ytilde, X_c):

    M_aug = [[1, b'], [b, Sigmahat]],  T_aug = lambda_max(M_aug).

Computed without forming M_aug: eigenvalues of Sigmahat plus one secular root
solving lam = 1 + sum_j (v_j'b)^2 / (lam - d_j) above the bulk edge (bisection
on [edge+eps, edge + ||b||^2 + 1]). Null location: plug-in root lam0 solving
the same equation with b replaced by its H0-equivalent magnitude
sqrt(beta' Sigma^2 beta / sigma_y^2) estimated by sqrt(tr(Sigmahat^2)/p /
sigma_y_hat^2), sigma_y_hat^2 = tr(Sigmahat)/p + sigma_eps_hat^2, where
sigma_eps_hat^2 = min(1, median(d_j)) heuristic documented in estimators.py
(median eigenvalue of a white-design covariance is ~1; conservative choice
min with 1 protects against strong-spike designs). Analytic threshold:
T_aug > lam0 + kappa * width_np with Johnstone width at dims (n, p+1)
rescaled by 1/n, kappa = TW1_Q95 = 0.9793 (size-0.05 intent); an empirically
calibrated variant T_aug > quantile_{0.95}(null reps) is co-recorded. Both
rejections stored per rep (`rej_s1_analytic`, `rej_s1_mc`).

DE prediction for power (F13): under H1 the cross-moment mean shifts by
E[b] = (Sigma beta + Lambda gamma)/sigma_y, so the secular-root location moves
by Delta(lam) approx gamma' Lambda'(lam I - Sigma)^{-1} Lambda gamma /
sigma_y^2 evaluated at lam near lam0 (plus the O(1/sqrt(p)) cross term which
vanishes under A4a). Unlike S2, this responds to ALL confounding mass,
subcritical-aligned included, because (lam I - Sigma)^{-1} does not project
onto supercritical directions only. Frontier prediction s_detect computed by
inverting Delta(g) = kappa * width.

S2 maxz_cal (supercritical-aligned detector; matches F9's s_eff):
    T_z = max_{j <= ktop} |z_j| / sqrt(1 + c d_j),
ktop = max(r_hat_onatski(Sigmahat), 1) capped at 10. Threshold: two-sided
Gaussian with Bonferroni over ktop effective coordinates at alpha = 0.05:
|T_z| > z_{1 - 0.05/(2 ktop)}. Correlations between z_j (shared beta draw)
make Bonferroni conservative; the MC-calibrated empirical 95th percentile of
T_z under H0 is co-recorded (`rej_s2_mc`). Power follows the supercritical
projection of gamma only: mu_j = sqrt(n l_j)/(sqrt(1+l_j) sigma_y) g dir_j,
so subcritical-aligned confounding is invisible to S2 by design.

S3 onatski_ratio_whitened: Onatski ED ratio statistic applied to the spectrum
of the whitened cross-moment Gram G := Sigmahat^{+1/2} (b b') Sigmahat^{+1/2}
* n restricted to the top-(ktop+5) coordinates; rhat = onatski_select on its
descending eigenvalues; reject if rhat >= 1. Implemented as a rank-one
special case; retained because the ledger named it, expected to be dominated
by S1/S2 (documented expectation, not a gate statistic).

Baselines:
B1 f_test_pcs: partial F-statistic of Ytilde on the top-ktop PC scores
(ktop as in S2); classical residual F-test with (ktop, n - ktop) degrees of
freedom. Under dense alternatives it pools all directions weakly.
B2 ucm_strength_boot: Rendsburg-et-al.-style confounding-strength estimate
(hatted via the PE-decomposition implemented in estimators.ucm_strength,
documented approximation of their estimator) with a bootstrap threshold:
refit on B = 200 parametric-bootstrap replicas (gamma = 0 resampled residuals);
reject if point estimate exceeds the 95th percentile. Costly; run only on the
WP 2.3 headline cells, flagged optional per shard.
B3 scree inspection = S0 (same object, reported as baseline arm).

## Numerical Le Cam probe (WP 2.3 action 3)

Two-sample discrimination between pooled H0 and H1 datasets at matched
configuration. Feature map (fixed, frozen here; dimension 12 + 4 r):

    feat(Y, X) = [log eigenvalues d_1..d_10 of Sigmahat, tr Sigmahat,
                  lambda_max(M_aug), ||b||^2,
                  top-min(r,4) values of |z_j|/sqrt(1+c d_j)].

Classifier: gradient-boosted trees (sklearn HistGradientBoostingClassifier,
default params, 50% train/calibration split) AND a Gaussian-kernel MMD with
median heuristic on standardized features. Report AUC of each on held-out
data. Operational declaration: computational undetectability at level
AUC <= 0.55 (both probes). LABEL: computational probe only, not an
information-theoretic statement; OMH (2013) contiguity below the BBP
threshold implies no max-eigenvalue alarm works there, while the probe bounds
what flexible spectral-summary discriminators can extract.

## Alignment stress test (WP 2.3 action 4)

theta swept from 0 (gamma aligned with strongest factor) to pi/2 (aligned
with weakest reported factor direction e_2) at fixed (l, c, g, n). Under H0
theta is undefined (gamma = 0), so SIZE must be theta-invariant by
construction; any size drift flags a bug, not a phenomenon. POWER is recorded
for S1, S2, B1 separately: theory predicts S2 power decays like cos(theta)
(supercritical-aligned mass shrinks) while S1 power depends on total mass and
decays much more slowly; the empirical contrast is the operational content of
the decoupling claim C1/C2 boundary mapping.

## What would falsify the detection ansatz (predeclared)

1. Size of S1 (analytic threshold) outside [0.02, 0.15] in >= 30% of null
   cells after the plug-in refinement (one fix iteration allowed: recalibrate
   kappa analytically, not per-cell).
2. Empirical power-1/2 contour of S1 off the F13 frontier by more than
   factor 3 in >= 30% of probed cells.
3. S2 calibrated size outside [0.035, 0.065] systematically (would mean the
   1 + c d variance law is wrong beyond finite-n corrections).
