# T3: Scoped Impossibility and the Visibility Boundary

Status: statements frozen 2026-08-25; proof skeletons grounded in the local
copy lit/omh2013.pdf (Onatski-Moreira-Hallin 2013, arXiv:1306.4867;
verified title/ids). Numerical anchors already on disk: results/power_surface.csv
(blind strata 9/9 confirmed), results/lecam_probe_auc.csv (c-split).

## T3(a): Eigenvalue-alarm impossibility below the BBP threshold (ADAPT)

### What "eigenvalue alarm" must mean (scope fix from the plan wording)

Pure-X eigenvalue statistics T(X) = spec(X_c'X_c/n) are INDEPENDENT of
gamma under every model in the ledger class - their impossibility is
trivial (data-processing identity, no asymptotics needed) and is recorded
as Remark 0 rather than a theorem. The substantive claim concerns
AUGMENTED-eigenvalue alarms: statistics measurable wrt the spectrum of

    M_aug = [[1, b'], [b, Sigma_hat]],   b = X_c'Ytilde/n,

(S1 t_aug and relatives), i.e., eigenvalues of the sample second-moment
matrix of the standardized pair (Ytilde, X_c).

### Statement target (to prove)

Fix c in a compact subset of (0, inf) \ {1}, r fixed, spike profile with
l_j <= sqrt(c) for all j (subcritical), sigma_u = sigma_eps = 1, A4a.
Let H0: gamma = 0 and H1: gamma = g * dir with ||dir|| = 1. Then:

(i) the joint laws of spec(M_aug) under H0 and H1 are mutually contiguous
    whenever g^2 * K(l, c) -> 0, where K(l,c) := sum_j l_j dir_j^2 /
    ((sqrt(c) - l_j)^2-ish) is the OMH-type distance-to-threshold
    functional (exact form derived in step A2 below);
(ii) consequently no test based on spec(M_aug) has asymptotic power exceeding
    size uniformly over such alternatives (Le Cam two-point + contiguity);
(iii) at the boundary g^2 K -> const > 0 the log likelihood ratio converges
    to a Gaussian process in the direction index - the OMH Theorem-7
    analogue - giving a strictly positive power envelope INSIDE the
    subcritical region near threshold (this is what makes the scoped claim
    sharp rather than pessimistic; it matches the measured S2/S1 behavior
    at supercritical-aligned cells only).

### Proof skeleton (grounded)

Step A1 (augmented-model reduction). M_aug is the top-left (p+1)-corner of
the Gram of the standardized pair; its population version is
Sigma_aug = [[sigma_y^{-2} sigma_eps^2 + m'm, m'],[m, Sigma]] with
m = (Sigma beta + Lambda gamma)/sigma_y. Under A4a the beta part of m is
isotropic O(1/sqrt(p)); the gamma part is the rank-one coherent piece
Lambda gamma / sigma_y. Thus H0-vs-H1 is exactly an OMH rank-one
alternative problem on the FIRST COORDINATE of a spiked population matrix,
with effective spike strength

    h_eff^2 = ||Lambda gamma||^2 / (sigma_y^2 * scale(M_aug-null edge)),

and the response STANDARDIZATION enters through sigma_y in the denominator
(the D5/E3 subtlety: sigma_y grows with gamma, deflating h_eff).

Step A2 (contiguity region). Apply the OMH machinery (their Theorem 7:
log L(h; .) = -1/2[Delta_p(z0(h)) - ln(1 - h/cp)] + o_p(1), uniform in
h <= h-bar < threshold) with their spherical alternative replaced by the
first-coordinate rank-one deformation of Sigma_aug. The needed
generalization: their null is sphericity; ours has bulk-plus-spikes. The
adaptation argument: inside the bulk, spike directions behave like
regular MP mass for z away from detached eigenvalues (KY anisotropic law,
arXiv:1410.3516, Theorem 3.7); the LR statistic depends on the spectrum
only through MP-small-deviation functionals that are continuous wrt the
bulk measure, so the same Delta_p functional arises with c replaced by the
effective aspect and thresholds by the deformed edges. Deliverable: exact
K(l, c) and the boundary curve g*(c, l-profile).

Step A3 (transfer to all eigenvalue alarms). Le Cam's lemma converts LR
contiguity into power bounds for ANY measurable function of spec(M_aug);
monotone data-processing covers coarser alarms (S0, S3-whitened variants
after checking measurability).

Guardrails (specific traps):
* DO NOT drop sigma_y from h_eff: the whole point of the D5/E3 finding is
  that self-normalization pushes h_eff down as g grows; an analysis in raw
  coordinates "proves" detection where none exists.
* The beta-isotropy term contributes O(1/p) shifts that are NOT uniformly
  negligible across j near the bulk edge - track them (they set the
  finite-p floor seen in the probe at low c).
* Contiguity gives impossibility for FIXED alternatives; the uniform-over-
  g statement needs the standard two-point packing argument - include it.

Numerical falsifiers (frozen): results/power_surface.csv blind strata -
max emp power <= 0.25 at max g in every stratum predicted blind (already
9/9); if any eigenvalue alarm beats 0.25 inside the claimed contiguity
region at n >= 2000, this theorem program is wrong.

## T3(b): Visibility boundary for quadratic functionals in (c, g) (NEW + adapt)

Empirical driver (Phase 2, gate_verdicts.WP23_lecam): GBM-probe AUC stays
chance-level up to g ~ 1.6 at c = 0.8 but reaches 1.0 at c = 0.2 - the
||b||^2 sigma-inflation signature is readable exactly when the beta-mass
floor is small RELATIVE to the link-induced variance shift.

Statement target: delta-method phase curve for the probe feature
q(b) = ||b||^2 = b'b under Y-standardization:

    H0:  q ~ q0 + noise with E[q0] = tr(Sigma^2 + sigma_eps^2 Sigma)/p /
         sigma_y^2-scale floor growing like 1 + c d_max-ish,
    H1:  E[q] shifts by delta(g, c) = 2 g^2 ||Lambda-dir weighting||^2 /
         sigma_y^2(g)^2 + higher order,

so the DETECTABLE-BY-PROBE region is {delta(g,c) >= kappa * sd_0(q)}:
a curve g_vis(c) increasing in c, matching the measured ordering
g_vis(0.2) << g_vis(0.8). Deliverable: closed-form delta(g,c) and sd_0(q)
under A4a + Gaussianity; validate against lecam_probe_auc.csv within MC
tolerance INCLUDING the c-ordering; then state the two-class map:
{eigenvalue alarms: impossible below BBP}, {quadratic probes: possible iff
above the visibility curve}. Guardrail: sd_0(q) is governed by FOURTH-
moment geometry of beta (Haar) times spectral tails - do not substitute
Gaussian-beta simplifications silently (A4a fixes Haar); keep both variants
side by side and simulate the Haar one.

Stop rules: T3 total 4 weeks (per research plan). Fallbacks documented
there (r = 1 subfamily; Gaussian-design-first for T3(b)).
