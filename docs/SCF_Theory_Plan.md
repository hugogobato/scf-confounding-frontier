# SCF Theory Plan (Phase 4)

Status: COMPLETE (2026-08-25, session 4). All work packages closed; one
package (T3a) closed AS REVISED after its own falsifier overturned the
planned theorem - details in docs/theory_T3a_eigenvalue_contiguity.md.
Theory never blocked submission (give-up rule 1); nothing outstanding.

This file operationalizes the research plan's Phase-4 target table into
work packages. It follows the theory-development workflow: validity audit
first, one frozen assumption ledger, honesty tags, one verifiable unit per
session, a numerical falsifier per number-predicting claim, and a specific
guardrail per package.

## Preflight certificate (validity audit, compressed)

Audited specification (from `docs/model_card.md`, `docs/assumption_ledger.md`,
both hash-pinned into every simulation config):

- Observable space: (Y, X) in R^n x R^{n x p}; latent f in R^{n x r}.
- Model class M1: X = Lambda f + u, Y = beta'X + gamma'f + eps, with
  Lambda = Q diag(sigma_u sqrt(l_j)), Q Haar-conditioned (q_fixed), dense
  beta (A4a: Haar, ||beta|| = 1), Gaussian errors (A2/A3/A5).
- Estimand: the loading-conditional bias functional Bias(T) =
  E[T(Y,X) | Lambda, Q] - beta (NOT the Haar-unconditional mean, which
  vanishes by symmetry; model card Section 4 records why the conditional
  functional is the meaningful one).
- Identification: (beta, gamma) not separately identified from second
  moments (plan Section 3.2); all theory targets below are functionals of
  (Sigma_X, Lambda gamma) or of the testing problem relative to ledger A4,
  hence constant on observational equivalence classes. No target claims to
  recover gamma itself.
- Consistency/non-vacuity: witness configurations exist for both qualitative
  regions (invisible-yet-harmful, visible-yet-harmless) and were CONFIRMED
  numerically in WP 1.5 and Phase 2 (gate verdicts on disk).

Checks performed: typing of every object entering T1-T7 (dimensions,
normalizations per model card Section 2); quantifier audit of each target
statement below; degenerate-limit tests (r = 1; c -> 1 boundary; g = 0);
consistency with known impossibility results (OMH-style power limits for
max-eigenvalue alarms below the BBP threshold, respected rather than
contradicted by T3(a)'s scope).

Verdict: PASS (no material issue; caveats carried as explicit scope notes:
T3(b)'s visibility boundary is derived first for Gaussian designs;
T4's optimality class must exclude sign-flipping weight tricks).

## Frozen assumption ledger reference

Single source of truth: `docs/assumption_ledger.md` (hash 11b162ac814d with
the model card). Every package cites ledger rows; none restates them.

## Theory-asset inventory

| Source | Headline content | Contributes to | Limitation we must carry |
|--------|------------------|----------------|--------------------------|
| Baik-Ben Arous-Peche (2005) | BBP transition, outlier locations | T1 (locations), T2 | real case only as stated; no eigenvector joint law |
| Benaych-Georges-Nadakuditi (2011, 2012) | spike overlaps xi(l,c), singular values | T1 overlaps layer, T2 | rectangular case separate; fluctuations not needed |
| Johnstone (2001) | TW1 null law for lambda_max | T2 finite-n null | white-noise null; our H0 is spiked-design + beta leakage |
| Onatski-Moreira-Hallin (2013) | asymptotic power of sphericity-type tests; contiguity machinery | T3(a) impossibility | their alternative is sphericity, not dense-factor links |
| Dobriban-Wager (2018); Hastie et al. (2022) | exact ridge/min-norm risk via resolvent calculus, c>1 | T1 (resolvent skeleton), T4 (min-norm geometry) | Frobenius/prediction risk objects, NOT loading-conditional directional means with gamma weights - the glue |
| Knowles-Yin (2017); Bloemendal et al. local laws | anisotropic local laws for spiked sample covariance | T1 rigorous backbone if pushed beyond DE level | heavy; we target DE-level theorem with sharp constants instead |
| Cevid et al. (2020) | spectral transform class, perturbed sparse model | T4 (estimator class definition) | no optimality statement to inherit |
| Nava et al. (2026) | SDBoost EB tuning, boosting path | T6 collapse lemma | linear special case only |

Single most valuable deliverable: **T1** (the capture-law theorem) - it is
the paper's spine, fully de-risked by the Phase-2 overlay (91.7%/100% tiers
within 10%, ridge median deviation 0.0%), and its r = 1 piece is provable in
isolation for an early complete win.

## Work packages

### TP-1 (T1) Capture-law theorem  [tag: adapt]

EXECUTION STATUS (2026-08-25, session 2): T1.b CLOSED at r = 1. The
Wishart-route joint-moment assembly that failed at (6.708, 5) was diagnosed
as a bookkeeping bug (A_e limit is sqrt(l) t/sqrt(n), not sqrt(l t/n); the
t = 1 cell had masked it), after which the collapse sqrt(l)t/(1+t(1+l)) =
sqrt(l)/(c+l) is an EXACT algebraic identity in t - no joint-moment
calculation needed; all cross-moment means vanish exactly by Haar-probe odd
symmetry. Full proof + companion artifact theorem R2 (E[q'Pi q] ->
(c-1)/(c+l)) written into docs/theory_T1_capture_law.md Sections 4.1-4.9;
adversarially audited by an independent agent (checklist A-N all PASS,
85+ fresh falsification cells, verdict PROOF STANDS; one Section-4.9
exposition erratum and one variance-bound typo found and fixed in place).
Fixed-r extension fully reduced (Woodbury + vanishing lemmas L1-L3,
Section 5), numerically confirmed at r = 2 incl. a subcritical spike.
T1.c ridge interpolation ALSO CLOSED same session: corrected shifted-
resolvent capture cap(lam) = (1+l)m_bar/(1+(1+l)m_bar) with m_bar the root
of lam m^2 + (lam+c-1)m - 1 = 0 (Section 6); decisive reconciliation
falsified the PROVISIONAL xi-split form (max err 0.081 + lambda drift vs
0.003 for the proved form); code updated (ridge_capture now proved form,
old kept as ridge_capture_superseded); lam->0 collapse to the min-norm law
is exact, replacing the old test that enshrined the superseded constant.
Remaining: componentwise r>=2 write-up (clerical).

Target statement (r = 1 first, then fixed r):

(T1.a, provable) c <= 1: E[beta_hat_OLS - beta | Lambda] =
Sigma_X^{-1} Lambda gamma EXACTLY at every n (joint-Gaussian independence
argument; already written in de_formulas.ols_bias_vector docstring; needs
formal write-up).

(T1.b, adapt) c > 1, sigma_u = 1 normalization, min-norm OLS beta_hat =
X'(XX')^{-1}Y: the loading-conditional mean bias decomposes as

    E[beta_hat - beta | Q] = sum_j (cap_j - 1) <beta, q_j> q_j
                             - (1/c - 1) beta_perp
                             + sum_j cap_j sqrt(l_j)/(1+l_j) gamma_j q_j,

with cap_j = (1 + l_j)/(c + l_j) and beta_perp = beta -
sum_j <beta,q_j> q_j. Consequently the RMS-over-A4a-beta directional norm
matches minnorm_total_bias_norm (already implemented and overlay-validated).

Sub-tasks (one session each):
1. Reduce to r = 1, write the three terms as expectations of explicit
   resolvent functionals of the spiked sample covariance; prove the
   rowspace artifact term -(1/c - 1)<beta_perp, .> from
   E[P_null] = (1 - 1/c) I + O(spike corrections) (DW/Hastie skeleton).
2. Prove the spike-direction coefficient: show
   E[<beta_hat - beta, q_j>] = (cap_j - 1)<beta,q_j> +
   cap_j sqrt(l_j)/(1+l_j) gamma_j using the leave-one-out/spiked-resolvent
   identity for q_j' Sigma_hat^+ q_j and the BGN overlap; THIS is the new
   glue. Guardrail: the naive guess xi/nu + (1-xi)*E[1/T_bulk] gives the
   SUPERSEDED formula (xi + (1-xi)/c), which simulation REJECTED (0.607 vs
   measured 0.657 at (l, c) = (3 sqrt 5, 5)); any derivation that passes
   through independent leaked-mass averaging is WRONG - the leaked mass
   correlates with the spike through the shared null space. Collapse check:
   l -> 0 must return cap = 1/c; l -> infinity must return cap = 1.
3. Extend r = 1 -> fixed r (orthogonality of distinct spike coordinates via
   asymptotic eigenvector orthogonality; cross terms o(1)).
4. Ridge-lambda interpolation: same three terms with
   cap_j(lam) = xi nu/(nu+lam) + (1-xi)(1/c)(1 - lam m_inv(lam, 1/c))
   (already coded as ridge_capture, flagged PROVISIONAL); prove it or scope
   it to lam >= c-limited regime. Numerical hook: ridge overlay median
   deviation currently 0.0% - any proof contradicting it is wrong.
Numerical falsifier: results/correctness_overlays.csv cells at n = 8000,
r in {1, 5}: predicted-vs-simulated directional means within MC tolerance
(frozen tolerance 10%; measured max ~4%).
Stop rule: 4 weeks; fallback Gaussian-design-only version, general case
conjectural.

### TP-2 (T2) Frontier achievability + finite-n null law for S2  [adapt]

CLOSED 2026-08-25 (session 4; docs/theory_T2_frontier.md). The null
inflation mechanism is a marginal-scale effect, not correlation: calibrated
coordinates are CONDITIONALLY INDEPENDENT Gaussians (Lemma A - the plan's
"correlated chi-square" guess superseded), pairwise covariances vanish
exactly, and kappa_n,k factorizes into an explicit scale-calibration ratio
(times a documented positive mixture-spread correction). At c > 1 the
frozen se2/sigma_y^2 estimator conventions mis-scale systematically on the
n-side spectrum (closed forms via the restricted-MP median: kappa_top =
3.50 at (5, sub); measured sd(z_cal) 4.21 incl. spread) - recorded as an
erratum-quality finding with NO gate impact (MC calibration was primary).
Achievability: limiting experiment = independent Gaussian shifts; NP
envelope + product-form power attained by S2-MC; slope law carries the BGN
overlap weight xi^{1/2} and is quadratic-in-g in second moments (sign of
v_j'q_j is rep-random - means carry no signal), satiating via
sigma_y^2 = A + g^2. Subcritical ceiling proposition makes "infinite
frontier" precise. Falsifier: tests/test_theory_T2.py.
Stop rule: met.

### TP-3 (T3) Scoped impossibility + visibility boundary  [adapt + NEW]

(a) CLOSED-AS-REVISED 2026-08-25 (session 4;
docs/theory_T3a_eigenvalue_contiguity.md). The planned OMH-contiguity
impossibility for spec(M_aug) was FALSIFIED BY ITS OWN PRE-REGISTERED
FALSIFIER: the true lambda_max(M_aug) separates CONSISTENTLY inside the
no-detachment region at c = 0.2 (AUC 0.29 -> 0.05 as n goes 800 -> 3200 at
g = 3.2), with a direction that FLIPS with geometry (upward, AUC -> 1,
even at the T3(b) visibility fixed point omega = omega* while every bulk
functional stays flat). The response-standardization channel breaks the
maximal-invariant structure OMH rely on; the retraction, the corrected
normalization sigma_y^2 = A + g^2 (the response carries the RAW factor
link), the PROVED population detachment boundary (Theorem A: top
detachment needs omega_e > c + 2 sqrt(c), never satisfied by scaling sub
profiles; bottom wake only for c > 4 - confirmed by fresh simulation),
and a PIPELINE ERRATUM of independent consequence are all documented:
the frozen t_aug silently mis-brackets its secular bisection (disc < 0
collapses hi below d[0]) and has been a degenerate surrogate correlated
with lambda_max(Sigmahat) all along - gate verdicts unaffected (S1 demoted
pre-data; MC thresholds self-consistent), but "invisible to eigenvalue
alarms" must be re-scoped to the S2 family and probe-blind bands in any
manuscript use.
(b) CLOSED 2026-08-25 (session 3; docs/theory_T3_visibility_boundary.md):
the probe feature is q(b) = ||b||^2 with b the self-standardized cross-
moment. Closed forms (all verified at n=2000 to <=2%): v0 = 1+sigma_eps^2,
M0 = 1+c(1+sigma_eps^2), floor m0 = M0/v0, mean-shift curve
E[q|g] = (M0 + g^2(omega + c))/(v0 + g^2) with omega = sum l_j dir_j^2.
THE VISIBILITY LAW: delta(g) satiates at |delta|max = |omega_star - omega|,
omega_star = 1/(1+sigma_eps^2); quadratic probes are BLIND at every g iff
|omega - omega_star| <= A sd_0(q0). For scaling profiles l ~ sqrt(c) the
sub-profile hits the fixed point exactly AT c = 1 (kappa_sub = 0.5(1-sqrt(c))):
the Marchenko-Pastur boundary is where quadratic probes are born blind.
Sign map (sub deflate / mixed inflate / c=1 null) matches forensics;
frozen-csv ordering and plateaus reproduced as the class map, with the
honest caveat that small-g frozen AUCs exceed the Gaussian mean-shift
envelope (GBM exploits shape beyond first moment; recorded, not tuned).
Falsifier: tests/test_theory_T3b.py.

### TP-4 (T4) Hard-trim dominance  [adapt]

CLOSED 2026-08-25 (session 4; docs/theory_T4_hard_trim_dominance.md).
Class corrected during falsification: spectral REGRESSION maps
beta_hat_m = V diag(m(d_hat)) V' b_raw, m >= 0 (projection-type maps are
not estimators of beta - caught by the first run, documented). The whole
dominance question reduces to per-direction transmission coefficients
T_j(m) inherited from closed T1 results: T(ols) = cap_j (c>1), T(ridge) =
cap_j(lam), T(trim|retained supercritical) = xi -> 1, T(trim|subcritical or
dropped) = 0. Theorem E: on all-subcritical cells the minimal
confounding-attributed floor is ZERO, attained exactly by pi = 0 trims;
any soft family transmits strictly and loses at every g bounded from zero;
Onatski selection consistency gives data-driven attainment. Formalizes the
G3 kill (oracle tau no-regret 6.2%, pca_onatski wins 89/97) as weight-SHAPE
necessity, not tuning failure. Falsifier: tests/test_theory_T4.py (fresh
twins; csv anchors matched to +0.2%/+2.5%).
Stop rule: met.

### TP-5 (T6) SDBoost collapse lemma  [NEW small, elementary]

Under A4a the EB variance-components ratio sr2/se2 -> 0 and all LAVA
weights -> 1, so the Nava et al. linear special case reduces EXACTLY to
OLS along its boosting path at any stopping time. Evidence: byte-identical
mean-bias vectors in 42/97 harmful cells. One-week write-up; falsifier =
equality check over returned means npz (already on disk).

### TP-6 (T7) Trimmed-tau DE for M2  [NEW small]

CLOSED 2026-08-25 (session 3; docs/theory_T7_trimmed_tau.md): exact FWL
identity for the trim arm and exact Sherman-Morrison identity
tau_hat_ols = d'G^-1Y/(1+d'G^-1D) for the min-norm arm (verified to 15
digits); plim DEs proved at elementary level. The six-fold OLS inflation is
the shrinkage term -tau/(1+Lambda_D) with Lambda_D = (1-1/c)||pi||^2 +
sum_j delta_j^2 t(1+t)/(1+t(1+l_j)) + sigma_nu^2 t (predicted 1.558 vs
measured 1.5548 +/- 0.072 at c=2); the trim is immune because [D,S] has
full column rank. Trim flatness across c follows from the delta-gamma
survival channel rho_j in [1/(1+l_j), 1] and a c-flat denominator v_D.
CLT stated with separated-spikes scope guardrail. Falsifier:
tests/test_theory_T7.py.

## Claims register

| ID | Claim | Tag | Numerical hook |
|----|-------|-----|----------------|
| C-T1a | exact c<=1 identity | provable (PROVED) | identity unit test (1e-8) |
| C-T1b | capture decomposition c>1 + artifact R2; fixed-r componentwise theorem PROVED (Section 5 write-up complete, session 4) | PROVED | correctness_overlays.csv <= 10% at n=8000; tests/test_theory_T1.py |
| C-T1c | ridge interpolation | PROVED at r=1 (closed 2026-08-25): shifted-resolvent cap(lam) = (1+l)m_bar/(1+(1+l)m_bar), m_bar root of lam m^2+(lam+c-1)m-1 (Stieltjes derivation closed session 4); superseded xi-split FALSIFIED by decisive reconciliation | ridge overlays median deviation 0.0%; tests/test_theory_T1.py ridge cells |
| C-T2 | kappa_n,k scale law + frontier achievability + subcritical ceiling | CLOSED (independence lemma supersedes correlated-chi-square guess; c>1 estimator mis-scale recorded) | tests/test_theory_T2.py; frozen ratios <= 1.5 |
| C-T3a | eigenvalue-alarm impossibility below BBP | REVISED: planned impossibility FALSIFIED by own falsifier; detachment boundary proved; true augmented alarm detects consistently with geometry-dependent direction; S1 pipeline erratum issued | tests/test_theory_T3a.py |
| C-T3b | visibility boundary phase curve | CLOSED (mean-shift class map + saturation ceiling; small-g shape envelope caveat) | lecam_probe_auc.csv ordering; tests/test_theory_T3b.py |
| C-T4 | hard-trim dominance on subcritical cells within nonnegative spectral regression maps | CLOSED (transmission law anchored in T1.b/T1.c/BGN; strict-loss corollary) | tests/test_theory_T4.py; estimation_cell_detail.csv anchors |
| C-T6 | sdboost collapse | provable | npz byte-equality 42/97 |
| C-T7 | trimmed-tau plim + CLT | CLOSED (both arms; OLS shrinkage mechanism verified) | m2_treatment.csv flatness; tests/test_theory_T7.py |

## Dependency graph

TP-1.1 -> TP-1.2 -> {TP-1.3, TP-6(T7)} ; TP-2 and TP-3a independent of TP-1
(start parallel); TP-3b needs TP-1.2 vocabulary only; TP-4 needs TP-1.2.
TP-5 standalone.

## Risks

R1: the T1.2 glue may resist at DE-rigor level; fallback = state as
deterministic-equivalent theorem with the local-law step cited as
conjectural (plan rule: theory never blocks submission).
R2: T3(b) delta-method variance asymptotics under A4a may need
Gaussian-X-first restriction.
R3: T4's weight-class definition may admit pathological exclusions; iterate
class definition against the ablation grid before writing the proof.

## Immediate next actions

1. Manuscript assembly: the theory layer is complete. Class map for the
   detection section must cite the REVISED scopes: pure-X trivial,
   S2-family ceiling (T2 Prop C), probe visibility law (T3b), augmented-
   alarm consistency with geometry-dependent direction (T3a revised) -
   do NOT quote "invisible to eigenvalue alarms" without the S1-erratum
   caveat.
2. Optional refinements (non-blocking, queued): analytic law for
   E[lambda_max(M_aug)] shift under H1 (T3a open mechanism target);
   subcritical eigenvector-CLT constant for the S2 ceiling; finite-n
   score-capture interpolation for T7's rho_j.

STATUS SNAPSHOT (2026-08-25, end of session 4): ALL Phase-4 packages
CLOSED. Session-4 closures - T2 (kappa law + achievability), T4
(dominance), clerical T1 items (fixed-r componentwise theorem write-up;
Stieltjes micro-gap). T3(a) closed AS REVISED: the planned OMH
impossibility was falsified by its own falsifier; surviving content =
proved detachment boundary + low-edge channel prediction + retraction
trail + S1 pipeline erratum (t_aug surrogate bug, gate verdicts
unaffected). Suite state: tests/test_theory_T{1,2,3a,3b,4,7}.py all green
alongside the Phase 1/2 suites.
