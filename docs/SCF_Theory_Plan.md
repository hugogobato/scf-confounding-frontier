# SCF Theory Plan (Phase 4)

Status: ACTIVE for T1, T2, T3(a) (G3 cleared 2026-08-24); T3(b), T6, T7 run
alongside Phase 3 (numerical falsifiers on disk); T4 after T1; T5 cut by
default; theory never blocks submission (give-up rule 1).

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

Target: (a) under H0, t_maxz = O_p(max over k of standardized coordinate)
has a finite-n inflation factor kappa_n,k over the Bonferroni Gaussian scale
(measured ~1.7x at Phase-2 scales); derive kappa from the max of correlated
chi-square-type coordinates with the D10 noise-aware variance law.
(b) power >= 1/2 at s_eff = s_detect(l,c) as predicted by the F12 slope
construction; achievability = upper bound matching the OMH power template.
Guardrail: the composite alternative reduces to aligned mass only through
the capture-weighted slopes; dropping omega_j breaks the frontier ratio
calibration (measured 0.93/1.02/1.25 across strata).
Numerical falsifier: WP 2.3 power_surface.csv g80/g_pred ratios within 1.5.
Stop rule: 3 weeks.

### TP-3 (T3) Scoped impossibility + visibility boundary  [adapt + NEW]

(a) No eigenvalue-alarm statistic detects sub-BBP dense confounding under
response standardization: two-point argument in the contiguity regime
(OMH technique), scoped to statistics measurable wrt the eigenvalues of the
augmented moment matrix. Status: subcritical blindness confirmed 9/9 strata.
(b) NEW object from the Le Cam probe split: delta-method phase curve in
(c, g) for concentrated quadratic functionals of b = X'Ytilde/n explaining
why ||b||^2-class probes are blind at c = 0.8 up to g ~ 1.6 but informative
at c = 0.2 (relative shift vs beta-mass floor growing with p). Deliverable:
predicted AUC-vs-g curve matched to results/lecam_probe_auc.csv within MC
tolerance including the c-ordering. Guardrail: the sigma_y-standardization
makes the H0/H1 pair NONstandard (sigma-normalization subtlety flagged since
E3/D5); a proof that ignores the random denominator is void.
Stop rule: 4 weeks total; fallbacks documented in the research plan.

### TP-4 (T4) Hard-trim dominance  [adapt]

Under separated spikes and A4a, Onatski-truncated OLS attains the minimal
directional mean-bias floor within estimators measurable wrt spectral
coordinates with NONNEGATIVE weights. Formalizes the G3 kill (oracle tau
no-regret 6.2%; pca_onatski wins 89/97). Class definition guardrail: must
exclude sign-flipping tricks; Onatski selection consistency enters.
Collapse check: k = r recovers pca_oracle_r behavior; single-spike case
must reduce to cap-law arithmetic. Numerical falsifier: ablation grid -
every fixed soft-weight family loses somewhere the theory says it must.
Stop rule: 3 weeks; fallback restricted class (diagonal monotone weights).

### TP-5 (T6) SDBoost collapse lemma  [NEW small, elementary]

Under A4a the EB variance-components ratio sr2/se2 -> 0 and all LAVA
weights -> 1, so the Nava et al. linear special case reduces EXACTLY to
OLS along its boosting path at any stopping time. Evidence: byte-identical
mean-bias vectors in 42/97 harmful cells. One-week write-up; falsifier =
equality check over returned means npz (already on disk).

### TP-6 (T7) Trimmed-tau DE for M2  [NEW small]

plim and CLT for the Onatski-trimmed treatment coefficient under dense
confounding via Frisch-Waugh: tau_trim inherits the capture-law bias floor;
tau errors flat in c while OLS inflates six-fold (results/m2_treatment.csv).
Feeds the package's adjustment documentation. 1-2 weeks.

## Claims register

| ID | Claim | Tag | Numerical hook |
|----|-------|-----|----------------|
| C-T1a | exact c<=1 identity | provable (PROVED) | identity unit test (1e-8) |
| C-T1b | capture decomposition c>1 + artifact R2 | PROVED at r=1 (elementary Wishart route, audited); fixed-r write-up clerical | correctness_overlays.csv <= 10% at n=8000; tests/test_theory_T1.py |
| C-T1c | ridge interpolation | PROVED at r=1 (closed 2026-08-25): shifted-resolvent cap(lam) = (1+l)m_bar/(1+(1+l)m_bar), m_bar root of lam m^2+(lam+c-1)m-1; superseded xi-split FALSIFIED by decisive reconciliation (max err 0.081 vs 0.003) and preserved as ridge_capture_superseded | ridge overlays median deviation 0.0%; tests/test_theory_T1.py ridge cells |
| C-T2 | frontier achievability + kappa_n,k | adapt | power_surface ratios <= 1.5 |
| C-T3a | eigen-alarm impossibility below BBP | adapt | blind strata power <= 0.25 at max g |
| C-T3b | visibility boundary phase curve | NEW/adapt | lecam_probe_auc.csv match incl. c-ordering |
| C-T4 | hard-trim dominance | adapt | estimation_cell_detail.csv dominance pattern |
| C-T6 | sdboost collapse | provable | npz byte-equality 42/97 |
| C-T7 | trimmed-tau plim + CLT | NEW | m2_treatment.csv flatness |

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

1. TP-6 (T7 trimmed-tau DE): unblocked end-to-end - consume the proved
   capture law via Frisch-Waugh; feeds the package docs.
2. TP-3b visibility-boundary curve (highest novelty per week).
3. Fixed-r componentwise write-up for T1 (clerical; lemmas L1-L3 proved).
4. Optional: rerun ridge overlay figures with the proved cap(lam) and
   record the (small) prediction deltas in a memo line.
