# T3(a): Augmented-Eigenvalue Alarms - Detachment Boundary, Falsified
# Impossibility Claim, and the S1 Pipeline Erratum

Status: CLOSED-AS-REVISED 2026-08-25 (session 4). This document REPLACES the
first draft of the T3(a) theory note written earlier the same session. The
draft claimed OMH-type mutual contiguity of spec(M_aug) below a population
detachment boundary. Its own pre-registered falsifier FALSIFIED that claim:
the true lambda_max(M_aug) separates CONSISTENTLY (AUC -> 1 with n) deep
inside the no-detachment region, in a direction that FLIPS with the
confounding geometry. What survives is proven and reported here: (i) the
population detachment boundary (Theorem A, elementary, correct);
(ii) the retraction with full forensic trail (Section 4); (iii) a pipeline
erratum of independent consequence: the frozen S1 statistic t_aug has been
computing a degenerate surrogate, NOT lambda_max(M_aug), since Phase 1
(Section 5). Grounding: lit/omh2013.pdf, lit/bgn2011.pdf,
lit/knowles_yin_anisotropic.pdf.

## 0. Scope and statistic-class map

Remark 0 (trivial layer, unchanged). Pure-X statistics spec(X_c'X_c/n) are
independent of gamma; power equals size identically. S0 belongs here (its
"power" under gamma != 0 reflects design spikes only).

Class under study: tests measurable wrt spec(M_aug),
M_aug = [[1, b'],[b, Sigmahat]], b = X_c'Ytilde/n. The intended flagship
member is lambda_max(M_aug); S1 was supposed to realize it.

## 1. Setup

Model card M1 conventions; sigma_u = sigma_eps = 1 unless carried;
gamma = g dir; Haar beta (A4a). CRITICAL NORMALIZATION (erratum vs draft):
plim sd^2(Y) = beta'Sigma_X beta + ||gamma'||^2 + sigma_eps^2
             = A + g^2,          A := 1 + sigma_eps^2,
because the response carries the RAW factor link gamma'f (variance
||gamma||^2 = g^2), NOT its loading-filtered image (that image ||Lambda
gamma||^2 = g^2 omega enters Cov(X,Y) instead). Verified directly:
E[sd_y^2] at (c=0.2, g=3.2) measured 12.2065 vs A + g^2 = 12.24.

Population of the standardized pair z_i = (Ytilde_i, x_i'):
Omega(g) = [[1, u'], [u, Sigma_X]],  u = (Sigma_X beta + Lambda gamma)/sigma_y(g).
H0 vs H1 differ by a rank-2 population perturbation; the bulk and the
design spikes tau_j = 1 + l_j are common. Population secular equation for
detached levels (tau outside spec(Sigma_X)), deterministic limit:

    tau = 1 + sum_j kappa_j/(tau - tau_j),
    kappa_j := l_j g_j-units: l_j g^2 dir_j^2 /(A + g^2).

## 2. Theorem A: population detachment boundary (PROVED, elementary)

Let b_pm = (1 -/+ sqrt(c))^2. Define omega_e := sum_j l_j dir_j^2/(c +
2 sqrt(c) - l_j) (spike denominators positive for subcritical profiles).
(i) TOP: f(tau) = tau - 1 - sum_j kappa_j/(tau - tau_j) is strictly
increasing on (b_+, inf); a detached level exists iff

    sup_g Theta_top(g) = omega_e > B := b_+ - 1 = c + 2 sqrt(c).     (D_top)

For scaling sub-profiles (l_j = 0.5 sqrt(c)): omega_e = w0/(c + 1.5 sqrt(c))
with w0 = dir-weighted profile constant; at the Phase-2 cells:
omega_e/B = 0.235 (c=0.2), 0.081 (c=0.8), 0.036 (c=2.0): NEVER detached,
at any g. (Mass concentrated AT the edge can satisfy (D_top) only for
c <~ 0.26.)
(ii) BOTTOM: on (-inf, b_-) exactly one root exists for every nonzero
coupling; writing s(g) = g^2 omega/((A + g^2)-normalized coherent share),
the root sticks to b_- as s -> 0 and moves down continuously. Wake-up is
possible iff (1 - b_-) < W(b_-), W(t) := sum_j l_j dir_j^2/(tau_j - t);
equivalently for equal spikes (1 - b_-)((1 + l) - b_-) < omega. This holds
AUTOMATICALLY for c > 4 and fails at all three Phase-2 sub cells
(c = 0.2: 0.64 > 0.22; c = 0.8: 1.31 > 0.45; c = 2.0: 1.27 > 0.71).
(iii) Inside the no-detachment region the limiting ESDs coincide
(finite-rank invariance of the absolutely continuous part).
All three parts elementary; proofs unchanged from the draft (they were
never the problem).

## 3. What the falsifiers actually found (the draft's death certificate)

Pre-registered falsifier design: inside the no-detachment region, an
OMH-adapted contiguity claim implies lambda_max(M_aug) cannot separate
consistently; test at (c = 0.2, sub profile l = 0.5 sqrt(c), r = 3,
n in {800, 1600, 3200}, theta = pi/6, g = 1.6 and 3.2). Outcome
(Mann-Whitney AUC of lambda_max(M_aug), H1 vs H0, fresh draws):

    n=800:  AUC 0.286 (g=3.2), 0.351 (g=1.6)
    n=1600: AUC 0.180 (g=3.2), 0.314 (g=1.6)
    n=3200: AUC 0.051 (g=3.2), 0.133 (g=1.6)

Separation is CONSISTENT (tending to perfect, DOWNWARD: H1 roots smaller)
- the opposite of contiguity. Geometry control at the T3(b) visibility
fixed point (theta = pi/2, l = (0.5, 0.5), so omega = omega* = 1/A): the
cross-moment mass sum_j w_j^2 is FLAT (0.7014 vs 0.7066) and every bulk
functional (tr, tr d^2, d_top) is flat, yet lambda_max(M_aug) now separates
UPWARD, AUC 0.998 (n=1600) / 1.000 (n=3200) at g = 3.2. Additional facts:
the direction and existence of separation depend on the coupling geometry,
not merely on omega; the effect persists at fixed g with growing n.

CONCLUSION (forced): the laws of spec(M_aug) under H0/H1 are NOT mutually
contiguous in the no-detachment region; no OMH-style impossibility holds
for this class. Mechanistically, the response-standardization rescales the
entire first-coordinate coupling (u and the sampling variability of b
alike) through sigma_y(g), so the augmented spectrum carries a
self-normalization signature that has no analogue in OMH's sphericity-plus-
rank-one alternative; the maximal-invariant structure their proof relies on
does not transfer. The precise analytic law for E[lambda_max(M_aug)] shift
(composition of the coherent-mean channel, the beta-leakage rank-one spike
nu_beta(g) = beta'Sigma^2 beta/sigma_y^2 ~ 1/(A+g^2), and the unchanged
eps-channel) is registered as an OPEN derivation target; its fixed-point
structure is NOT the naive omega* one (falsified above).

Consequence for the project narrative (recorded verbatim): "subcritical
dense confounding is invisible to eigenvalue alarms" is FALSE for the true
augmented alarm; it was an artifact of testing a broken surrogate
(Section 5) plus the genuinely ceiling-bound S2. The invisible-yet-harmful
phenomenon remains real for the S2 family (docs/theory_T2_frontier.md,
Proposition C) and for any statistic class actually shown contiguous; the
paper's decoupling claim must cite those scopes, not spec(M_aug).

## 4. Retracted claims (kept for the audit trail)

From the draft of this note, retracted in place:
* "(i) joint laws of spec(M_aug) mutually contiguous whenever g^2 K -> 0"
  - FALSIFIED by the Section-3 measurements.
* "(ii) no spec(M_aug)-measurable test has power exceeding size" - falls
  with (i).
* "(iii) OMH Theorem-7 Gaussian-process analogue with positive envelope
  inside the layer" - unproven and moot for the class; the Laplace-layer
  transfer was never established (the draft's Section on Theorem B had
  tagged it [adapt] without proof).
* The draft's Corollary A1 saturation arithmetic used sigma_y^2 = A +
  g^2 omega (wrong normalization, Section 1); its qualitative never/edge
  conclusions for (D_top) happen to survive under the corrected
  normalization because omega_e <= B holds a fortiori, but the printed
  constants g_top = 13.4 etc. are void.
Kept: Theorem A (all three parts, with corrected normalization), Remark 0,
the low-edge wake criterion, and the pipeline erratum below.

## 5. Pipeline erratum: frozen S1/t_aug is a degenerate surrogate

Discovery path: the falsifier's lambda_max(M_aug) (validated against the
explicitly formed (p+1)x(p+1) matrix to 1e-8, permanent unit test) disagreed
with the pipeline's t_aug on identical datasets (e.g. H0 mean 2.185 vs
1.808 at (c=0.2, sub, n=2000)).

Root cause, code/detection.py compute_stats: the upper bracket for the
secular bisection is computed as the larger root of the conservative
quadratic lam^2 - (1+lo) lam + (lo + s_total) = 0 via disc = (1+lo)^2 -
4(lo + s_total). When s_total > (1 - lo)^2/4 this DISC IS NEGATIVE, the
sqrt(max(disc, 1e-300)) silently returns ~0, and hi collapses BELOW d[0]
(e.g. hi ~= 1.53 < lo = 2.07 at the anchor cell). The subsequent bisection
on the inverted interval converges below the bulk edge; the returned
"t_aug" is a smooth but meaningless surrogate strongly correlated with
lambda_max(Sigmahat) - exactly the behavior Erratum 3 documented
pre-data ("tracks the design spike") without identifying the bug.

Blast radius (checked):
* All stored t_aug values (nullcal/power/alignment/bench sweeps) are
  affected at cells where disc < 0; this includes most c <= 1 designs
  (threshold s_total > (1-d_1)^2/4 is easily exceeded there) and rarely
  triggers at c >= 2 (large d_1).
* Gate verdicts do NOT change: S1 was demoted to diagnostic BEFORE data
  (D5/E3) on independent grounds; its MC thresholds were computed from the
  same surrogate under H0 and H1 (self-consistent calibration); no gate
  rule used raw S1 decisions.
* Frozen figures/tables involving pow_S1_cal / t_aug must be read as
  "surrogate-S1", annotated accordingly in any manuscript use.
Fix (one line, robust bracketing): hi = d[0] + s_total/max(d_min,eps) + 1
(or the exact bound sum w_j^2/(lam-d_j) <= s_total/(lam-d_max)); recommend
adopting the validated root solver from tests/test_theory_T3a.py::_aug_stats
in any future S1 usage. NOT applied retroactively to frozen artifacts.

## 6. Validation ledger

| Check | Status |
|-------|--------|
| Secular-root identity vs explicit M_aug eigh (permanent unit test) | PASS |
| Theorem A(i) monotonicity/(D_top) constants, corrected normalization | PASS |
| Theorem A(ii) bottom wake criterion (corrected algebra) | PASS |
| F-B': consistent downward separation, c=0.2 sub (AUC table above) | PASS (falsifies draft claim) |
| F-B": fixed-point geometry control: flat Q & bulk, upward root shift | PASS (falsifies probe-disguise guess) |
| F-C: low-edge channel at c=5 sub opens (down-AUC 0.81 at g=1) | PASS (prediction P-1 confirmed) |
| Frozen blind strata csv assertion (surrogate-S1 + S2) | PASS (as recorded artifacts) |
| Bracket-bug reproduction (disc<0 -> hi < d[0]) | PASS (documented) |

Honest status: T3(a) closes with a NEGATIVE result on its planned theorem,
one proved population-boundary theorem, one confirmed structural prediction
(low-edge channel), one superseded-class-map correction, and one pipeline
erratum with bounded blast radius. The OMH machinery remains unused for
this problem; any future impossibility claim must first exhibit a statistic
class whose laws it can actually prove contiguous.
