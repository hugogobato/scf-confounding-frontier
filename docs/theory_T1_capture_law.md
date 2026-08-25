# T1: The Capture-Law Theorem (directional mean-bias DE for min-norm OLS)

Status: derivation document for TP-1 (research plan Phase 4 / SCF_Theory_Plan).
Conventions follow docs/model_card.md (sigma_u = 1; population eigenvalues
tau_j = 1 + l_j of Sigma_X; c = p/n; A4a dense Haar beta with ||beta|| = 1;
loading-conditional functional E[ . | Q]).

## 0. What is proved, what is adapted, what is checked

* T1.a (c <= 1): EXACT finite-n identity - PROVED here (Section 3), already
  unit-tested at 1e-8 (tests/test_identities.py).
* T1.b (c > 1): the three-term decomposition with cap_j = (1 + l_j)/(c + l_j)
  - stated as a deterministic-equivalent theorem; the proof reduces to two
  resolvent lemmas (Section 4) whose verification is the 'adapt' glue
  (anisotropic local laws for spiked sample covariance, Knowles-Yin
  backbone). The DE content is FULLY VALIDATED by the frozen Phase-2
  overlay (results/correctness_overlays.csv: <= 10% in 91.7% of n = 2000
  cells, 100% at n = 8000) and by the pilot at 0.5-1% accuracy.
* Guardrail history (why naive routes fail): TWO candidate derivations were
  tried and REJECTED before this formulation -
  (i) the independent-bulk split q' Sig+q -> xi/nu + (1 - xi) * E_bulk[T^-1]
      produces the SUPERSEDED formula xi + (1 - xi)/c (bgn_capture_superseded),
      which simulation rejected (0.607 predicted vs 0.657 measured at
      (l, c) = (3 sqrt 5, 5));
  (ii) a Sherman-Morrison conditional-expectation chain that treats the
      noise direction u_q as uncorrelated with the Gram resolvent - it
      drops the Uq-cross-terms inside G and gives denominator (l + c - 1)
      instead of (l + c).
  Both failures come from the same source: the leaked spike mass and the
  null space are correlated through the shared noise. Any correct proof
  must keep that correlation (this is exactly what the local-law/QVE
  machinery does).

## 1. Setup and notation (r = 1 first)

X = sqrt(l) f q' + U, f ~ N(0, I_n), U iid N(0,1), q unit, Y =
X beta + gamma f + eps. Min-norm OLS beta_hat = X'(XX')^{-1}Y (p > n).
Gram G := XX' (n x n, invertible a.s.), Sigma_hat^+ := X'G^{-1}X/n.

Decompose the bias along the spike direction:

    <beta_hat - beta, q> = A + B + C,
    A := <(P - I) beta, q>          (fit artifact; P = Sigma_hat^+ n-scale)
    B := gamma * q' X'G^{-1} f      (confounding channel)
    C := q' X'G^{-1} eps            (mean-zero noise)

## 2. What the elementary algebra can and cannot do (full autopsy)

The confounding channel is

    B/gamma = (sqrt(l) f + Uq)' G^{-1} f,
    G = XX' = l ff' + sqrt(l)(f q'U' + U q f') + H,   H = UU'.

CRITICAL CORRECTION (recorded because two derivation attempts died here):
G is NOT H + l ff'. The CROSS TERMS sqrt(l) f q'U' + sqrt(l) U q f' are
O_p(sqrt(l n p)) - the same order as H - and NEVER vanish. Consequently:

* Attempt (ii-a) "Sherman-Morrison with G = H + aa'" is VOID from the
  start; any constant derived through it (e.g. denominator l + c - 1) is
  meaningless.
* Attempt (iii) raw-resolvent reduction q'Sigma_hat^+q -> 1/(c+l):
  FALSIFIED numerically (sim 0.0052 vs 0.1667 at (l,c)=(4,2)); the channel
  involves cancellations the raw resolvent does not see.
* A Schur-complement computation in the (e = f/||f||)-basis (kappa := U'e,
  eta := U - e kappa', D-tilde := ||kappa||^2 - ||P_r kappa||^2 with
  P_r the rowspace projector of eta) correctly reproduces the auxiliary
  laws T = f'H^{-1}f -> (p-n+1)^{-1} ~ (n(c-1))^{-1} and
  V = (Uq)'H^{-1}f = o_P(1), but it CANNOT reach B, because applying SM to
  G requires the cross terms again. There is no shortcut around them.

CORRECT TARGET LEMMA (R'), stated directly on the channel:

    LEMMA R': E[<beta_hat - beta, q> | Q] = sqrt(l)/(c + l) gamma + o(1),
    with companion artifact statement R2: E[q' Pi q] = (c - 1)/(c + l).

STATUS OF R': VALIDATED NUMERICALLY, NOT YET PROVED. Micro-grid (r = 1,
gamma = 1.3, n = 160, 80 reps):

    (l, c)     sim bias coef   theory sqrt(l)/(c+l)*gamma
    (4, 2)     +0.4388         +0.4333
    (9, 2)     +0.3609         +0.3545
    (6.708, 5) +0.2896         +0.2876
    (0.5, 2)   +0.3660         +0.3677
    (2, 1.5)   +0.5243         +0.5253
(max deviation 1.8%, consistent with MC noise; extends the Phase-2
full-grid overlay validation of minnorm_total_bias_norm.)

Guardrail history now FOUR documented failed routes:
  (i) independent-bulk xi-split -> superseded formula xi + (1-xi)/c;
  (ii) Sherman-Morrison dropping Gram cross-correlations;
  (iii) raw-resolvent reduction -> Lemma R, falsified;
  (iv) SM/Schur-completum elementary chain -> invalid because the Gram
       cross terms are first-order; any "exact identity" that starts by
       deleting them is wrong.
Common root cause: X = a q' + U makes the spike interact with the noise
through FIRST-ORDER cross terms; only resolvent machinery that treats the
full deformed Gram (local laws / QVE for spiked sample covariance)
qualifies as a proof vehicle.

Guardrail history now THREE documented failed naive routes:
  (i) independent-bulk xi-split -> superseded formula xi + (1-xi)/c;
  (ii) Sherman-Morrison chain dropping Uq-Gram cross-correlations ->
       denominator (l + c - 1);
  (iii) raw-resolvent reduction -> Lemma R, falsified above.
Common root cause: leaked spike mass, noise direction, and null space are
mutually correlated through the shared U; only identities that keep the
exact cancellation (like R') survive.

Collapse checks (must hold in any proof):
  l -> 0:  cap -> 1/c (bulk capture n/p);
  l -> inf: cap -> 1 (spike resolved);
  c -> 1+: 1/(c + l) must match the c <= 1 exact identity's implied
  resolvent limit (continuous transition; the c <= 1 branch has NO
  null-space artifact term).

## 3. T1.a: the c <= 1 exact identity (PROOF)

For p <= n, Sigma_hat invertible a.s. and Y = X beta + w with
w = gamma f + eps INDEPENDENT of X (joint Gaussianity of (f, u, eps) makes
w _|_ X even though X depends on f: condition on f; then X|f and eps are
independent, and gamma f is a constant shift). Then

    beta_hat - beta = Sigma_hat^{-1} X' w,
    E[beta_hat - beta | X] = Sigma_hat^{-1} X' E[w] = 0 ... wait no:

E[w] = 0 but the bias comes from Corr(w, X): write w = X a + zeta with
a = Sigma_X^{-1} Lambda gamma and zeta = w - Xa INDEPENDENT of X
(Gaussian projection). Then

    E[beta_hat] - beta = E[Sigma_hat^{-1} X'] X a = a   (since Sigma_hat^{-1} X'X = I),

EXACTLY at every n. Hence plim/E-identity: E[beta_hat] - beta =
Sigma_X^{-1} Lambda gamma = sum_j sqrt(l_j)/(1 + l_j) gamma_j q_j.
QED. (This is F1; unit-tested.)

Note the c <= 1 case has NO fit-artifact term: P = I identically.

## 4. Proof route for R'/R2 (adapt level; the only surviving vehicle)

Route A (deformed-Gram resolvent calculus). Grounded citations (working
copies in lit/, verified 2026-08-25):

  * Knowles-Yin, "Anisotropic local laws for random matrices"
    (arXiv:1410.3516; PTRF 169): anisotropic local law
    Sigma^{-1}G(z) - Pi(z)Sigma^{-1} = O_prec(Psi(z)) (their Theorem 3.6);
    outside-spectrum version (their Theorem 3.7); Pi(z) defined via the
    population QVE solution m(z) (their Section 3). Their Remark 3.8
    explicitly flags outlier regions as requiring the companion analysis -
    TP-1's evaluation point z = 0-below-the-bulk ALONG THE SPIKE DIRECTION
    is exactly such a region, so the missing glue is: extend their
    outlier-side resolvent to spike-ALIGNed test vectors and take z -> 0.
  * BGN (arXiv:0910.2120): overlap constants xi(l, c) entering the
    statement's interpretation; not needed inside Route A's proof chain.
  * Dobriban-Wager (arXiv:1507.03003) / Hastie et al. (arXiv:1903.08560):
    ridge-resolvent fixed points; used for the T1.c interpolation layer.

Proof program:

  A1. Express <beta_hat - beta, q> through an augmented resolvent identity
      that keeps ALL Gram cross terms (exact finite-n algebra).
  A2. Apply KY Theorem 3.7 to replace deformed resolvents by DEs on
      regular directions; handle the spike-aligned direction with the
      outlier extension of Remark 3.8 (the new content).
  A3. Evaluate at z = 0, simplify to sqrt(l)/(c + l); control errors at
      order n^{-1/2}.

Numerical falsifier (pre-proof screening): the micro-grid table in
Section 2 IS the falsifier for R'; tolerance 5% at n = 400 across the grid,
shrinking with n. Currently PASSING at n = 160 within 1.8% (results/
theory_T1_checks.csv). Any future "proof" whose constants disagree with
this table is wrong - investigate, do not explain away.

## 5. r > 1 extension sketch

Distinct spike directions decouple: q_i' Sigma_hat^+ q_j = o(1) for i != j
(eigenvector delocalization/orthogonality of distinct spikes, BGN 2011
Corollary-level fact), so the coefficient vector enters componentwise with
l_j-specific caps; the beta_perp rowspace artifact is common. Cross terms
between the artifact and confounding channels vanish at A4a (E<beta,q_j> =
0 exactly by independence). This gives the vector statement implemented in
minnorm_bias_vector / minnorm_total_bias_norm.

## 6. Ridge interpolation (T1.c)

Replace Sigma_hat^+ by (Sigma_hat + lam I/n)^{-1}; the same QVE evaluation
at z = -lam/n gives cap_j(lam); consistency anchor lim_{lam->0} matches
Lemma R. Empirically the ridge overlays are EXACT to MC noise (median dev
0.0%), which any formula must reproduce. ridge_capture currently carries
the xi-split form flagged PROVISIONAL - if its predictions match sim while
Lemma-R-based algebra differs, reconcile BEFORE proving either (both cannot
be right unless they coincide; check numerically first).

## 7. Validation ledger for this document

| Check | Status |
|-------|--------|
| T1.a identity unit test 1e-8 | green since Phase 1 |
| overlay <= 10% share (n=2000) | 91.7% PASS (Phase 2 gate) |
| overlay <= 10% share (n=500/8000) | 100% PASS |
| ridge overlay median deviation | 0.0% PASS |
| R' micro-grid r=1 (5 points, n=160) | PASS, max dev 1.8% |
| raw-resolvent Lemma R | FALSIFIED (recorded as guardrail) |
| SM/Schur elementary routes | INVALID - Gram cross terms first-order |

Honest status: T1.b is a VALIDATED deterministic-equivalent theorem whose
rigorous proof (Route A) remains open work; per the plan's binding rule,
theory never blocks submission - the paper states the capture law at DE
level with the overlay evidence and marks the local-law step as the
formalization gap.
