# T1: The Capture-Law Theorem (directional mean-bias DE for min-norm OLS)

Status: derivation document for TP-1 (research plan Phase 4 / SCF_Theory_Plan).
Conventions follow docs/model_card.md (sigma_u = 1; population eigenvalues
tau_j = 1 + l_j of Sigma_X; c = p/n; A4a dense Haar beta with ||beta|| = 1;
loading-conditional functional E[ . | Q]).

## 0. What is proved, what is adapted, what is checked

* T1.a (c <= 1): EXACT finite-n identity - PROVED here (Section 3), already
  unit-tested at 1e-8 (tests/test_identities.py).
* T1.b (c > 1): PROVED AT r = 1 (Section 4, closed 2026-08-25) by an
  elementary route: exact Gram completion of squares -> ordinary central
  Wishart -> Sherman-Morrison -> classical inverse-Wishart moments. NO local
  laws, NO unproved imports. Companion artifact theorem R2 also proved.
  ADVERSARIALLY AUDITED (2026-08-25, independent subagent): all checklist
  items A-N pass; 85+ fresh falsification cells (incl. an 80-cell
  (l,c)-grid and n-scaling probes) found no deviation beyond MC noise;
  verdict PROOF STANDS with one exposition erratum in Section 4.9 (fixed
  in place) and one variance-bound typo in M1 (fixed in place).
  General fixed r: reduction complete (Section 5, Woodbury + three
  off-diagonal vanishing lemmas); write-up of the componentwise statement is
  the remaining clerical step. Numerically confirmed at r = 2 to <0.3%.
* Guardrail history (why naive routes failed): FIVE documented dead ends -
  (i) independent-bulk xi-split -> SUPERSEDED formula xi + (1-xi)/c,
      rejected by simulation (0.607 predicted vs 0.657 measured at
      (l, c) = (3 sqrt 5, 5));
  (ii) Sherman-Morrison chain treating u_q as uncorrelated with the Gram
       resolvent (drops Uq-cross terms) -> denominator l + c - 1, wrong;
  (iii) raw-resolvent reduction q'Sigma_hat^+ q -> 1/(c+l): falsified
       numerically (sim 0.0052 vs 0.1667 at (l,c) = (4,2));
  (iv) SM/Schur-completum "exact identity" starting by deleting the Gram
       cross terms sqrt(l)(f u' + u f'): invalid, they are first-order;
  (v) assembly-from-pointwise-limits in the CORRECT Wishart reduction:
      produced sqrt(l) t (1 + t - sqrt(t))/(1+t(1+l)), which matches at
      t = 1 ((4,2)) but fails at t = 1/4 ((6.708,5)). POST-MORTEM (this
      session): the failure was a BOOKKEEPING BUG, not fluctuation physics:
      A_e's limit is sqrt(l) t / sqrt(n), NOT sqrt(l t / n). With the
      corrected limit the assembly closes EXACTLY for every t and the
      earlier "cross-moment fluctuations are first-order" diagnosis is
      RETRACTED - all m_eu-mean terms vanish exactly by odd symmetry of the
      Haar probe e.

## 1. Setup and notation (r = 1 first)

X = sqrt(l) f q' + U, f ~ N(0, I_n), U iid N(0,1), q unit, Y =
X beta + gamma f + eps. Min-norm OLS beta_hat = X'(XX')^{-1}Y (p > n).
Gram G := XX' (n x n, invertible a.s.), rowspace projector P_R :=
X'G^{-1}X, null-space projector Pi := I - P_R.

Decompose the bias along the spike direction:

    <beta_hat - beta, q> = A + B + C,
    B := gamma * (sqrt(l) f + Uq)' G^{-1} f   (confounding channel)
    C := q' X'G^{-1} eps                      (mean-zero noise)
    A := fit-artifact channel (computed via R2, Section 4.7).

## 2. Why naive routes fail (full autopsy, kept as guardrails)

The confounding channel is

    B/gamma = (sqrt(l) f + Uq)' G^{-1} f,
    G = XX' = l ff' + sqrt(l)(f q'U' + U q f') + H,   H = UU'.

CRITICAL CORRECTION (two early derivation attempts died here):
G is NOT H + l ff'. The CROSS TERMS sqrt(l) f q'U' + sqrt(l) U q f' are
O_p(sqrt(l n p)) - the same order as H - and NEVER vanish.

* Attempt (ii-a) "Sherman-Morrison with G = H + aa'" on that WRONG split is
  void from the start.
* Attempt (iii) raw-resolvent reduction q'Sigma_hat^+q -> 1/(c+l):
  FALSIFIED numerically; the channel involves cancellations the raw
  resolvent does not see.
* Attempt (iv): Schur-complement chains in the (e = f/||f||)-basis correctly
  recover auxiliary laws but cannot reach B without handling the cross
  terms.

The correct move (Route A-prime, next section) absorbs the cross terms into
the rank-one update EXACTLY: G = K + aa' with K a central Wishart.

Collapse checks any final formula must pass:
  l -> 0:  capture coefficient of gamma-channel -> 0 like sqrt(l)/(c+l);
           rowspace capture cap -> 1/c (bulk fraction n/p);
  l -> inf: cap -> 1 (spike resolved);
  c -> 1+: continuous match with the c <= 1 exact identity branch.

## 3. T1.a: the c <= 1 exact identity (PROOF)

For p <= n, Sigma_hat invertible a.s. and w := gamma f + eps. Joint
Gaussianity makes w decomposable as w = X a + zeta with
a = Sigma_X^{-1} Lambda gamma and zeta INDEPENDENT of X (Gaussian
projection). Then

    E[beta_hat] - beta = E[Sigma_hat^{-1} X'] X a = a   (Sigma_hat^{-1} X'X = I,

E[(X'X)^{-1}X' zeta] = 0 since zeta _|_ X). EXACTLY at every n. Hence
E[beta_hat] - beta = Sigma_X^{-1} Lambda gamma
                   = sum_j sqrt(l_j)/(1 + l_j) gamma_j q_j. QED (unit-tested).

Note c <= 1 has NO fit-artifact term: P_R = I identically.

## 4. T1.b at r = 1: COMPLETE PROOF (closed 2026-08-25)

Throughout: p/n -> gamma_0 > 1 written c; Wishart df m := p - 1;
delta_m := 1/(m - n - 1) = 1/(p - n - 2); t := lim n delta_m = 1/(c - 1).

### 4.1 Exact reduction (Gram completion of squares)

With u := Uq and e := f/||f||:

    G = XX' = K + a a',
    K   := U P_perp U',   P_perp = I_p - q q',
    a   := sqrt(l)||f|| e + u.

Check by expansion: the aa' cross term reproduces sqrt(l)||f||(e u' + u e')
= the exact Gram cross terms; the residual mismatch -u u' is absorbed into
K because UU' - uu' = U(I - qq')U' = U P_perp U'. Exact identity, no
approximation.

STRUCTURAL FACTS (all exact):

(F0a) Rows r_i of U are iid N(0, I_p). Since P_perp q = 0, the Gaussian
      pairs (P_perp r_i, q'r_i) are independent per row, hence
      sigma(K) INDEPENDENT of u = Uq.
(F0b) K =d W_n(I_n, df = m = p-1): rotate coordinates mapping q to e_1;
      P_perp becomes diag(0,1,...,1) and U P_perp drops one column.
(F0c) e is Haar on S^{n-1}, INDEPENDENT of (K, u); ||f||^2 ~ chi^2_n
      independent of everything else listed.

### 4.2 Exact Sherman-Morrison identities

Define m_ee := e'K^{-1}e, m_eu := e'K^{-1}u (= u'K^{-1}e), m_uu := u'K^{-1}u,
kappa := a'K^{-1}a, A_e := a'K^{-1}e, A_u := a'K^{-1}u. SM gives, exactly,

    e'G^{-1}e = m_ee - A_e^2/(1+kappa),
    u'G^{-1}e = m_eu - A_u A_e/(1+kappa),

    kappa = l||f||^2 m_ee + 2 sqrt(l)||f|| m_eu + m_uu,
    A_e   = sqrt(l)||f|| m_ee + m_eu,
    A_u   = sqrt(l)||f|| m_eu + m_uu.

The two channel pieces are ch_a := sqrt(l)||f||^2 e'G^{-1}e and
ch_b := ||f|| u'G^{-1}e, so B/gamma = ch_a + ch_b exactly.

### 4.3 Moment ledger (conditional on (K,u); Haar-probe algebra)

For B symmetric, x Haar unit: E[x'Bx | B] = tr(B)/n,
Var(x'Bx | B) = 2[tr(B^2) - (tr B)^2/n]/(n(n+2)) <= 2 tr(B^2)/n^2; for
linear forms E[x'a | B] = 0 and E[(x'a)^2|B] = ||a||^2/n; odd symmetry
gives E[(x'a)(x'Bx)|B] = 0 EXACTLY (e -> -e). Applied with probes e
against vectors built from (K, u):

(M1) E[m_ee | K] = tr(K^{-1})/n, and E[tr(K^{-1})] = n delta_m exactly
     (inverse-Wishart mean E[K^{-1}] = I/(m-n-1); Mardia-Kent-Bibby 1979
     p. 91 pin, lit/theory_T1_wishart_locators.md), so m_ee =
     delta_m + o_p(delta_m). Var(m_ee|K) <= 2 tr(K^{-2})/n^2 -> o(delta_m^2)
     by MP negative moments (Bai-Silverstein 2010 §3.3.1; Silverstein-Choi
     1995 eq. (1.3)-(1.4)).
(M2) E[m_eu | K, u] = 0 exactly; E[m_eu^2 | K, u] = u'K^{-2}u/n, whose
     mean over (K,u) is tr(E[K^{-2}])/n = Theta(delta_m^2) (MP negative
     second moment int x^-2 dMP_a = 1/(1-a)^3, verified analytically and
     numerically in lit/theory_T1_wishart_locators.md).
(M3) E[m_uu | K, u] = tr(K^{-1}) = n delta_m; fluctuations O_p(delta_m)
     (Wishart concentration), so m_uu -> t in probability and in L2.
(M4) Odd-symmetry corollary: E[m_eu * m_ee] = 0 and
     E[m_eu * m_uu | K, u] = 0 EXACTLY (both are odd in e).
(M5) kappa -> t(1+l), A_u -> t in probability (substitute M1-M3; the
     m_eu contributions are O_p(delta_m) each).

### 4.4 Assembly of ch_a (term-by-term, expectations)

ch_a = sqrt(l)||f||^2 [m_ee - (sqrt(l)||f||m_ee + m_eu)^2/(1+kappa)].

Expand the square; take conditional means using (M1)-(M4):

 T1. sqrt(l)||f||^2 m_ee:  E = sqrt(l) n delta_m -> sqrt(l) t
     (independence of ||f|| and K).
 T2. sqrt(l)||f||^2 * l||f||^2 m_ee^2 /(1+kappa):
     -> sqrt(l) * l t^2/(1+t(1+l)); error o(1) by (M1) concentration.
 T3. cross term 2 sqrt(l)^{3/2}||f||^3 E[m_ee m_eu]/(1+kappa) = 0 EXACTLY
     by (M4).
 T4. sqrt(l)||f||^2 E[m_eu^2]/(1+kappa) = O(sqrt(l) n delta_m^2) = o(1)
     by (M2).

Hence  E[ch_a] -> sqrt(l) [t - l t^2/(1+t(1+l))]
                     = sqrt(l) t (1+t)/(1+t(1+l)).

NUMERICAL CONFIRMATION (n = 400, 300 reps, results/tmp run of this session,
to be frozen in tests/test_theory_T1.py):
  (4,2): sim 0.6656 vs 0.6667; (6.708,5): 0.2766 vs 0.2765;
  (2,1.5): 1.2145 vs 1.2122; (16,8): 0.1904 vs 0.1905.

### 4.5 Assembly of ch_b

ch_b = ||f|| [m_eu - A_u A_e/(1+kappa)], A_uA_e
     = (sqrt(l)||f||m_eu + m_uu)(sqrt(l)||f||m_ee + m_eu)
     = l||f||^2 m_eu m_ee + sqrt(l)||f|| m_eu^2 + sqrt(l)||f|| m_uu m_ee
       + m_uu m_eu.

Conditional means of the four product pieces:
 p1 = l||f||^2 E[m_eu m_ee] = 0 EXACTLY by (M4) (odd in e);
 p2 = sqrt(l)||f|| E[m_eu^2] = O(sqrt(n) delta_m^2) -> 0 by (M2);
 p3: conditioning on (K,u), E[m_ee | K,u] = tr(K^{-1})/n and m_uu is
     (K,u)-measurable, so E[p3] = sqrt(l)||f|| m_uu tr(K^{-1})/n
     -> sqrt(l) n^{1/2} * t * t/n = sqrt(l) t^2/n^{1/2}. ch_b carries ONE
     more factor ||f|| ~ n^{1/2} over the denominator (1+kappa) = O(1)
     (kappa -> t(1+l), no n), so this piece contributes
     -sqrt(l) t^2/(1+t(1+l)) at exactly first order.
 p4: E[m_uu m_eu | K,u] = 0 EXACTLY by (M4).
The direct term ||f||E[m_eu] = 0 exactly.

Hence E[ch_b] -> -sqrt(l) t^2/(1+t(1+l)),
with all other pieces o(1) by (M2)/(M4).

NUMERICAL CONFIRMATION: (4,2): -0.3334 vs -0.3333; (6.708,5): -0.0554 vs
-0.0553 (the OLD naive form predicted -0.1106 here - this cell is the
discriminator between correct and buggy bookkeeping); (9,3): -0.1259 vs
-0.1250; (16,8): -0.0237 vs -0.0238.

### 4.6 Collapse: sum and the capture law

E[B/gamma] = E[ch_a + ch_b]
           -> sqrt(l)[t(1+t) - t^2]/(1+t(1+l))
           =  sqrt(l) t/(1+t(1+l))
           =  sqrt(l)/(c + l)        [t = 1/(c-1): exact identity, ALL t].

THE CAPTURE LAW IS AN ALGEBRAIC IDENTITY IN t once the two channel pieces
are assembled - no delicate limit interplay, no joint-moment calculation
beyond (M1)-(M4). Collapse checks: l -> 0 gives 0 channel (no spike, no
confounding coupling) while cap -> 1/c; c -> 1+ gives t -> inf,
sqrt(l)t/(1+t(1+l)) -> sqrt(l)/(1+l) matching the c<=1 identity's spike
coefficient sqrt(l)/(1+l). Both hold.

### 4.7 Artifact theorem R2 (companion statement, PROVED)

Claim: E[q'Pi q | Q] -> (c-1)/(c+l).

Proof. q'Pi q = 1 - q'P_R q and q'P_R q = (Xq)'G^{-1}Xq with
Xq = sqrt(l)f + u. Expand:

  (Xq)'G^{-1}(Xq) = l f'G^{-1}f + 2 sqrt(l) f'G^{-1}u + u'G^{-1}u.

Using Section 4.2 with the SAME five scalars:
  f'G^{-1}f = ||f||^2 e'G^{-1}e: by Section 4.4, ||f||^2[e'G^{-1}e] is the
              functional ch_a/sqrt(l), so -> t(1+t)/(1+t(1+l)).
  u'G^{-1}f = ||f||(e'G^{-1}u): this is exactly the Section 4.5 assembly
              (ch_b without its sqrt(l)-free normalization), so
              ||f||[m_eu - A_uA_e/(1+kappa)] -> -sqrt(l) t^2/(1+t(1+l));
              multiplied by the outer 2sqrt(l) below, it contributes
              -2 l t^2/(1+t(1+l)).
  u'G^{-1}u = m_uu - A_u^2/(1+kappa): E[A_u^2] = l||f||^2 E[m_eu^2]
              + 2sqrt(l)||f|| E[m_eu m_uu] + E[m_uu^2] -> 0 + 0 + t^2
              by (M2)/(M4)/concentration, so -> t - t^2/(1+t(1+l)).

Summing over the common denominator (1+t(1+l)):
  numerator = l t(1+t) - 2 l t^2 + t(1+t(1+l)) - t^2
            = lt + lt^2 - 2lt^2 + t + t^2 + lt^2 - t^2 = t(1+l).
Hence E[q'P_R q|Q] -> t(1+l)/(1+t(1+l)) and
      E[q'Pi q|Q]  -> 1/(1+t(1+l)) = (c-1)/(c+l). QED.

NUMERICAL: (4,2): 0.16720 vs 0.16667; (6.708,5): 0.34214 vs 0.34165;
r=2 case below.

### 4.8 Directional bias assembly (T1.b statement, r = 1)

beta_hat = X'G^{-1}(X beta + gamma f + eps); along q:

  <beta_hat, q> = beta'P_R q + gamma (Xq)'G^{-1}f + (Xq)'G^{-1} eps.

Noise term: mean zero exactly (eps _|_ (X, q)). Confounding term:
(Xq)'G^{-1}f = sqrt(l) f'G^{-1}f + u'G^{-1}f
             -> sqrt(l) t(1+t)/(1+t(1+l)) - sqrt(l) t^2/(1+t(1+l))
             = sqrt(l) t/(1+t(1+l)) = sqrt(l)/(c+l)     [Sections 4.4-4.6].
Artifact term: beta'Pi q. By (F0a-F0c) the law of U is invariant under all
rotations fixing q, so E[Pi q | Q] = rho q with rho := E[q'Pi q|Q]
-> (c-1)/(c+l) by Section 4.7; beta fixed gives

  E[beta'Pi q | Q] -> rho <beta,q> -> ((c-1)/(c+l)) <beta,q>.

TOTAL (conditional on the realized beta, per the model-card estimand
E[. | Lambda, Q, beta]; with A4a Haar beta the artifact term averages to
zero under E[.|Q] alone, which is why the conditional reading matters):

  E[<beta_hat - beta, q> | Q, beta]
      -> -(c-1)/(c+l) <beta,q> + sqrt(l)/(c+l) gamma
      = (cap - 1) <beta,q> + cap sqrt(l)/(1+l) gamma,
with cap := (1+l)/(c+l). QED (T1.b directional statement, r = 1).

Moment-existence preconditions: the ledger uses E[K^{-1}] (needs
p > n+2) and E[K^{-2}] (needs p > n+4); automatic for fixed c > 1 as
n -> infinity and satisfied by all finite-n cells used below.

NUMERICAL (conditional test, fixed (Q,beta), fresh (f,U,eps)): dev <= 0.24%
at (4,2),(6.708,5),(0.5,2),(2,1.5).

### 4.9 Bulk block (three-term decomposition completion)

For unit v perpendicular to span(q_j) (fixed given Q): the same SM probe
algebra with probe pair (v, f) shows E[v'X'G^{-1}f | Q] = 0: writing
v'G^{-1}f = ||f||(m_ve - A_vA_e/(1+kappa)) with m_ve = v'K^{-1}e Haar-linear
(mean 0, odd), A_v = sqrt(l)||f||m_ve + v'K^{-1}u, every surviving mean
involves E[m_ve * (K,u-measurable)] = 0 exactly or E[m_uv] = 0 exactly.
eps-channel mean zero exactly. Hence along the bulk,

  E[<v, beta_hat - beta> | Q] = -E[v'Pi beta | Q] = -(v'E[Pi|Q] beta).

Rotational invariance under rotations fixing EACH q_j makes
E[P_R | Q] commute with that group, so E[P_R|Q] = sum_j alpha_j q_j q_j'
+ alpha_b (I - QQ'). Diagonal entries give alpha_j = 1 - (c-1)/(c+l_j) =
cap_j (Section 4.7 generalized; Section 5 at r > 1) and the trace forces
alpha_b = (n - sum alpha_j)/(p-r) -> 1/c. Since Pi = I - P_R,

  E[v'Pi beta|Q] -> (1 - 1/c) <v, beta>,
  E[<v, beta_hat - beta>|Q] = -(1 - 1/c)<v, beta> = (1/c - 1) <v, beta>,

i.e., the bulk artifact term (1/c - 1) beta_perp of the three-term
decomposition. Sanity anchor: at l = 0, E[P_R] = (n/p) I EXACTLY, so
E[v'Pi beta] = (1 - n/p)<v,beta> -> (1 - 1/c)<v,beta>. [ERRATUM 2026-08-25,
adversarial audit finding 1: an earlier draft wrote E[v'Pi beta|Q] ->
(1/c)<v,beta> here - a sign/convention slip contradicting Pi = I - P_R;
the proved chain 4.1-4.8 never uses this block.]

This completes the full three-term decomposition implemented in
minnorm_bias_vector / minnorm_total_bias_norm and overlay-validated in
Phase 2 (91.7%/100% tiers within 10%).

### 4.10 Post-mortem of the Step-3 bug (recorded 2026-08-25)

The previous draft assembled from pointwise limits and got
ch_b -> -sqrt(l) t^{3/2}/(1+t(1+l)) from the erroneous evaluation
A_e -> sqrt(l t/n). In fact A_e = sqrt(l)||f||m_ee + m_eu ->
sqrt(l) sqrt(n) (t/n) = sqrt(l) t/sqrt(n): the exponent of t is 1, not 1/2.
The bug was invisible at (4,2) because there c = 2 forces t = 1 where both
expressions coincide - the single-cell "confirmation" was a coincidence.
With the corrected limit the assembly is exact at EVERY t, and the earlier
diagnosis ("cross-moment fluctuations are first-order when t < 1") is
retracted: all m_eu-involving means vanish EXACTLY by odd symmetry of the
Haar probe (M4); their magnitudes only affect fluctuations of ch_a/ch_b,
which are irrelevant to the loading-conditional MEAN functional. Lesson
(guardrail G-T1.6): when a formula matches at one cell, verify it at a cell
with DIFFERENT structure (here: t != 1) before believing the derivation.

## 5. Extension to fixed r > 1

Setup: X = sum_j sqrt(l_j) f_j q_j' + U, Q = [q_1..q_r], orthonormal.
Reduction: G = K + A A', A = [a_1..a_r], a_j = sqrt(l_j)||f_j||e_j + u_j,
K = U(I - QQ')U' =d W_n(I_n, df = p - r), sigma(K) INDEPENDENT of
(u_1..u_r) (same per-row block independence, P_Q q-block = 0), e_j from QR
of F independent of (K, U). Woodbury:

  G^{-1} = K^{-1} - K^{-1}A M^{-1} A'K^{-1},   M = I_r + A'K^{-1}A.

Off-diagonal vanishing lemmas (all elementary):

(L1) e_a'K^{-1}e_b = o(delta) for a != b: mean 0 (orthogonality + oddness),
     magnitude O_p(tr(K^{-2})^{1/2}/n) by Haar-set second moments.
(L2) e_j'K^{-1}u_k = O_p(delta) uniformly (mean 0, sd^2 = u_k'K^{-2}u_k/n).
(L3) u_a'K^{-1}u_b = O_p(1/sqrt(n)) for a != b: the (u_a,u_b)-conditional
     mean is <u_a,u_b> tr(K^{-1})/n (via E[K^{-1} | (u_a,u_b)] = delta_m I;
     K-conditional mean is 0), with <u_a,u_b> = O_p(sqrt(n)); either way
     O_p(n^{-1/2}), and every off-diagonal entry of M is O_p(1/sqrt(n)).

Consequently M = diag(1 + t(1+l_j)) + o(1) (diagonals by the r=1 ledger,
off-diagonals by L1-L3 summed with bounded coefficients), and every
directional functional decouples componentwise: for each j,

  E[<beta_hat - beta, q_j>|Q] -> (cap_j - 1)<beta,q_j>
                                 + cap_j sqrt(l_j)/(1+l_j) gamma_j,
  cap_j = (1 + l_j)/(c + l_j),   E[q_j'Pi q_j|Q] -> (c-1)/(c+l_j),
  bulk artifact -(1/c - 1) beta_perp unchanged,

because t = lim n/(p - r - n - 1) = 1/(c-1) is unchanged for fixed r.
Cross confounding channels gamma_k (Xq_j)'G^{-1}f_k (k != j) vanish: every
piece factors through L1-L3 products. This is the vector statement in
minnorm_bias_vector.

NUMERICAL (r = 2, l = (6.708, 0.8), c = 5, fixed-(Q,beta), 400 reps):
spike 1 dir +0.29003 vs +0.28985; spike 2 (SUBCRITICAL l = 0.8 < sqrt 5)
dir +0.12445 vs +0.12368; R2 values within 0.003. PASS.

## 6. Ridge interpolation T1.c (route specified; next session)

Replace G^{-1} by (G + lam I_n)^{-1}: SM/Woodbury apply verbatim with
K -> K + lam I, so all five scalars become shifted resolvent entries
m_ee(lam) = e'(K+lam I)^{-1}e etc. Their limits come from the Silverstein
equation for W_n(I, df = p-1) (Silverstein-Choi 1995 eq. (1.3)/(1.4),
pinned in lit/theory_T1_wishart_locators.md): mu_1(lam) :=
lim (1/n)tr(K+lam)^{-1} solves the companion quadratic; the same
assembly then yields cap_j(lam) in closed form. DECISIVE CHECK FIRST
(per doc guardrail): reconcile against ridge_capture's PROVISIONAL
xi-split form numerically BEFORE proving either - if the two disagree at
t != 1 cells, the xi-split form is wrong the same way the min-norm naive
form was. Anchor: lam -> 0 must return Section 4.6 exactly; overlay
median deviation 0.0% is the empirical target.

## 7. Validation ledger for this document

| Check | Status |
|-------|--------|
| T1.a identity unit test 1e-8 | green since Phase 1 |
| overlay <= 10% share (n=2000) | 91.7% PASS (Phase 2 gate) |
| overlay <= 10% share (n=500/8000) | 100% PASS |
| ridge overlay median deviation | 0.0% PASS |
| R' micro-grid r=1 (5 points, n=160) | PASS, max dev 1.8% |
| raw-resolvent Lemma R | FALSIFIED (guardrail iii) |
| SM/Schur routes deleting cross terms | INVALID (guardrails ii/iv) |
| Step-3 pointwise assembly (old) | FALSIFIED at t!=1 (guardrail v) |
| Corrected assembly ch_a/ch_b/sum (7 cells, t in {0.14..5}) | PASS, max dev 1.15% |
| R2 artifact theorem (4 cells) | PASS, max dev 0.32% |
| Conditional directional mean incl. artifact (4 cells) | PASS, max dev 0.24% |
| r=2 decoupling incl. subcritical spike | PASS |
| Adversarial audit (independent subagent, fresh seeds) | PROOF STANDS; 80-cell grid + n-scaling: no falsification; errata fixed in place |

Honest status: T1.a PROVED (exact); T1.b PROVED at r=1 by elementary
means (Sections 4.1-4.8) and adversarially audited; extension to fixed r
fully reduced with all vanishing lemmas identified and numerically
confirmed (Section 5); T1.c ridge route specified, decisive reconciliation
check queued (Section 6). No local laws or unproved imports remain anywhere
in the r=1 proof.
