# Assumption Ledger (SCF Phase 1)

Ledger version: v1.0 (2026-08-23). The simulation config hash is sha256(model_card.md || assumption_ledger.md), truncated to 12 hex chars; it appears in every parquet row. Simulations must use the WEAK variant wherever weak and strong variants are distinguished.

Purpose: every claim C1-C7 of the research plan traces to a subset of these rows; every row names the claims it carries and its failure mode. This document exists because of the identification subtlety (research plan Section 3.2): from second moments of (Y, X) alone, only beta* = beta + Sigma_X^{-1} Lambda gamma is identified, not (beta, gamma).

## A1. Dense factor structure

- Statement: X = Lambda f + u with f in R^r, r fixed (or slowly growing, r = o(p^{1/2})), loadings Lambda = Q D with Q p x r Haar-isotropically distributed (dense rows), D = sigma_u diag(sqrt(l_1), ..., sqrt(l_r)), l_1 >= ... >= l_r >= 0 fixed as n, p -> inf with c = p/n in (0, inf).
- Weak variant: r fixed at 3 for all Phase 1-2 simulations. Strong variant: r growing; deferred to Phase 4 theory and never used in simulations before then.
- Role: produces the spectral profile (spikes or hidden bulk mass) that both the bias functional and the detection statistic act on.
- Failure mode if violated: sparse or discrete loadings change the overlap geometry (no isotropic delocalization); results become loading-dependent.
- Carries: C1, C2, C3, C4.

## A2. Error law

- Statement: u_ij, eps_i i.i.d. N(0, sigma^2) at finite n, all mutually independent and independent of f.
- Weak variant: Gaussian everywhere in Phases 1-3 gates. Strong/universal variant: sub-Gaussian or finite-moment errors with universality argued empirically only (plan C5 cut by default).
- Role: licenses exact finite-n distributional facts (TW fluctuation scale, Gaussian resolvent identities) and keeps unit tests deterministic-tolerance based.
- Failure mode: heavy tails inflate lambda_max beyond BBP/TW predictions; heteroskedasticity shifts the bulk edge and breaks naive TW calibration (robustness WP 2.4 owns this).
- Carries: C1 (DE overlay accuracy), C2 (size calibration). Claims C1/C3 qualitative regions should survive non-Gaussianity; that is tested, not assumed.

## A3. Independence of (f, u, eps)

- Statement: f_i, u_i, eps_i mutually independent across i and across sources.
- Role: makes Cov(X, Y) = Sigma_X beta + Lambda gamma an exact identity and makes the OLS/ridge bias functionals well defined without distributional assumptions on f beyond Var(f) = I.
- Failure mode: dependence between f and u redefines Sigma_X and moves the spikes; out of scope.
- Carries: everything.

## A4. Identification / alignment ledger (the material one)

Primary variant A4a (near-orthogonality, Cevid-style perturbed design):
- Statement: the structural direction beta is drawn independently of (Q, gamma, l) with E[beta] = 0, ||beta|| = 1, and delocalized coordinates (e.g., uniform on the sphere S^{p-1}); equivalently, beta has asymptotically negligible inner product with any fixed-r subspace spanned by population eigenvector directions u_j. Formally checkable by a referee: |<beta, u_j>| = O_p(p^{-1/2}) for each j <= r.
- Consequence: in whitened cross-moments the beta-signal spreads over the bulk while any confounding component Lambda gamma concentrates on spike directions; testing H0: gamma = 0 against H1: gamma != 0 is well posed with Lambda treated as a nuisance estimated from X.
- Failure mode: if beta correlates with top eigenvectors (aligned beta), part of the signal mimics a spike and size inflates; the alignment stress sweep (WP 2.3 action 4) quantifies exactly how much violation the test tolerates. Estimation claims C1/C3 do not need A4a; only C2 does.

Sensitivity variants:
- A4b (perturbed sparsity): beta = beta_s + Lambda gamma / sigma_u^2 with beta_s k-sparse, k = o(n/log p); the testing problem becomes detecting the dense perturbation component; same statistic family applies to residuals after Lasso partialling.
- A4c (multiple responses): m >= r responses Y^(t) = X beta_t + gamma_t' f + eps_t sharing the factor loadings; identification of the confounding subspace from the m x p cross-moment matrix when m >= r; primary use: strengthens power, no new estimator family.
- A4d (independent causal mechanisms, Janzing-Scholkopf UCM logic): beta independent of the spectrum of Sigma_X in the mechanism-independence sense; different estimator family (spectral concentration scores); Phase 1 uses it only as a comparison baseline concept, not as an assumption we maintain.

Standing rule: any C2 statement must name its variant. Default is A4a. If a result needs A4b/c/d instead, it says so inline.

## A5. Homoskedastic noise

- Statement: Var(u_ij) = sigma_u^2 constant in i, j; Var(eps_i) = sigma_eps^2 constant in i.
- Weak variant: enforced in all gate simulations. Robust variants (heteroskedastic u) are Phase 2 WP 2.4 stress tests with predeclared degradation tolerance (size < 0.15 or a documented robust fix).
- Failure mode: heteroskedasticity invalidates raw TW thresholds (bulk edge moves); biwhitening-type normalization is the documented recovery path.
- Carries: C2 size calibration, C4 alarm calibration.

## A6. Number of factors

- Statement: either r known (theory statements) or r estimated by the Onatski (2010) ratio rule applied to the sample spectrum of X (simulations using PCA-k baselines).
- Role: PCA-k trim needs k; misspecification by +-1 is a robustness cell, not a gate blocker.
- Failure mode: over-estimated k keeps confounding directions inside the retained subspace (bias survives); under-estimated k removes signal. Both effects are visible in the phase diagram and reported.
- Carries: C3 (baseline fairness), ablation attribution.

## Detection problem under A4a (C2's formal object)

Hypotheses about the observable joint law of (Y, X):
- H0: gamma = 0 (Y = X beta + eps).
- H1: gamma != 0 with ||Lambda gamma|| = Theta(1) and effective detection spike s(l, c, theta, g) >= 0 (defined below).

Statistic family (WP 1.1 action 3): let Sigmahat = X'X/n and W = XY/n (cross moment). Family members:
- S1: lambda_max( Sigmahat^{-1} W W' Sigmahat^{-1} / n ) with Tracy-Widom-calibrated threshold; the whitened cross-moment largest-eigenvalue test.
- S2: linear spectral statistics of Sigmahat^{-1} W W' with data-chosen test functions (EB-selected; Phase 2).
- S3: Onatski ratio statistic for the number of spiked directions in the whitened cross-moment matrix.
Baselines for comparison (not part of the family): residual F-test, UCM bootstrap (Rendsburg et al.), scree inspection.

Effective detection spike (working definition, DE level, refined in WP 1.3/F7): s_eff := || P_spike Sigma_X^{-1/2} Lambda gamma ||^2 where P_spike projects onto supercritical population eigenvector directions; intuition: only the supercritical-aligned part of b = Lambda gamma produces a coherent rank-r mean shift in the whitened cross-moment, the rest hides in the bulk. Under subcritical profiles s_eff = 0 at leading order regardless of g, which is the formal content of "invisible".

Frontier mapping: s_detect(c, theta, g) is the smallest g such that power(S1) >= 1/2 at size 0.05 given (l, c, theta); s_remove(c, delta) the smallest g such that worst-case rel_bias of the tuned spectral estimator <= delta. Claim C2 asserts these boundaries are computable and distinct; claim C1 asserts the corresponding bias boundary differs from detectability.

## Mechanical verification hooks

1. Config hash: sha256(card+ledger)[:12] recorded per row (see code/configs.py).
2. Grep rule: simulation code may reference only assumption symbols A1-A6, A4a-d; any new modeling knob requires a ledger edit first (checked by review at each gate memo).
3. Weak-variant rule: grid configs instantiate weak variants; strong variants appear only in docs/theory sections.
