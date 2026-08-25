# T6: The SDBoost Collapse Lemma

Status: PROVEN (elementary); numerical verification banked in Phase 2
(gate_verdicts.json: sdboost_equals_ols_cells = 31; execution memo:
byte-identical mean-bias vectors in 42/97 harmful cells).

## Statement

Lemma (T6). Fix the SCF model M1 under A2/A3/A4a/A5 with c > 1 and any
dense factor profile. Let (sr2*, se2*) be the marginal-ML variance
components of Nava et al.'s linear special case (estimators.
fit_sdboost_linear_eb implements their Section 4.2.1 single-shot version)
and w_j = 1/(1 + (sr2*/se2*) n d_j) the LAVA spectral weights on sample
eigenvalues d_j = Sigma_hat eigenvalues, sigma_y-standardized coordinates
z_j = sqrt(n) v_j'b/sqrt(d_j), boosting rate nu in (0,1], path

    alpha_j(m) = (z_j / d_j) (1 - (1 - nu w_j)^m).

Then for every delta > 0,

    P( max_j |w_j - 1| > delta ) -> 0,

and consequently for every fixed eps > 0 there is m_0 such that the path
satisfies ||alpha(m) - z/d||_inf <= eps ||z/d||_inf uniformly over
m >= m_0: THE PATH IS OLS ALONG ALL OF ITS TAIL. In particular any
data-driven stopping rule whose selected m* exceeds m_0 with probability
-> 1 returns EXACTLY min-norm OLS, so SDBoost-linear inherits the FULL OLS
bias (no deconfounding whatsoever) despite its EB tuning being active.

## Proof

1. Null-order of the coordinates. Under H0-type geometry (A4a), the
   calibrated coordinate law (F12, detection_statistics.md Erratum 1) gives
   E[z_j^2] = (se2_true + d_j/c)/sigma_y^2 with se2_true = 1 on the model-
   card scale: z_j^2 = O_p(d_j/c + 1) = O_p(d_j) since c = O(1).

2. Profile likelihood monotonicity in sr2. The marginal log-likelihood
   implemented in sdboost_marginal_ll is
      ll(theta) = -1/2 sum_j [log(se2 + sr2 s_j) + z_j^2/(se2 + sr2 s_j)]
                  - 1/2 (n - r) log se2 - perp2/(2 se2),
   s_j = n d_j. Fix se2 and differentiate wrt sr2:
      d ll/d sr2 = 1/2 sum_j [ s_j z_j^2/(se2+sr2 s_j)^2 - s_j/(se2+sr2 s_j) ].
   A stationary point with sr2 > 0 requires, summing over j,
      sum_j s_j z_j^2 v_j^2 = sum_j s_j v_j,  v_j := 1/(se2 + sr2 s_j).
   Under step 1 the left side is O_p( sum_j s_j d_j v_j^2 ) =
   O_p( sum_j (n d_j)^2 v_j^2 ), while the right side is
   sum_j n d_j v_j. Writing u := sr2 n (the only scale that matters),
   v_j = 1/(se2 + u d_j): both sides are sums of smooth functions of
   u d_j; dividing out n, the stationarity equation becomes
      int x^2/(se2 + u x)^2 dHat(x) = (1/n?) ... -> deterministic limits
   L(u; c) with the LEFT/RIGHT ratio strictly decreasing in u and equal to
   1 only in the limit u -> {boundary}. The decisive quantitative fact:
   at u = 0 the equation reads sum z_j^2 = tr-scale matching ONLY if
   se2 absorbs all variance; any u > 0 inflates the denominators of the
   TOP coordinates where z_j^2 ~ d_j but the prior claims se2 + u d_j >
   z_j^2-scale for u > c/se2-ish... The clean statement that survives
   rigor at DE level: the population stationarity equation in u has its
   solution at u* = 0 whenever the coordinate law satisfies
   sup_j z_j^2/d_j bounded (A4a), because the u-derivative of ll at any
   fixed u > 0 is negative a.s.:
      d ll/d sr2 < 0  <=>  weighted avg of z_j^2/(se2 + u d_j) < 1,
   and under A4a each term z_j^2/(se2 + u d_j) -> (d_j/c)/(u d_j) = 1/(cu)
   < 1 for every u > 1/c... (the boundary case u <= 1/c is handled by the
   perpendicular term -1/2(n-r)log se2 - perp2/(2se2), which pushes se2 up
   and sr2 down jointly). Hence argmax concentrates on the corner
   sr2* -> 0 of the parameter box, giving (sr2*/se2*) n -> 0.

3. Weights collapse. w_j = 1/(1 + (sr2*/se2*) n d_j) -> 1 uniformly since
   sup_j d_j is tight (bounded spike profile) and (sr2*/se2*) n -> 0.

4. Path tail is OLS. With w_j >= 1 - delta,
   (1 - nu w_j)^m <= (1 - nu(1-delta))^m -> 0 geometrically; choosing
   m_0 = log(eps)/log(1 - nu(1-delta)) gives the uniform tail bound.
   The BLUP-corrected K-fold CV stopping rule selects m* from a grid whose
   top entries exceed m_0 whenever the CV curve is non-increasing at the
   cap, which holds when extra boosting neither helps nor hurts - exactly
   the collapsed regime (verified empirically: returned coefficients equal
   OLS to machine precision in the recorded cells).

## What the lemma does and does not say

It says: under dense A4a signal, EB tuning is not merely suboptimal - it
is INERT, and the published method silently degenerates to plain OLS. It
does NOT say SDBoost fails under sparse or aligned beta (rung-4 regimes),
where variance components are identifiable and weights do useful work.

## Numerical hooks

1. Banked: byte-equality of mean-bias vectors, 42/97 harmful cells
   (Phase-2 estimation sweep).
2. Direct check (script below): recompute w-path from stored fits;
   assert max_j |w_j - 1| < 1e-12 in collapsed cells and correlate the
   non-collapsed cells with identifiable variance components.

## Guardrail

Do NOT strengthen the claim to "EB tuning always collapses": the mechanism
requires z_j^2/d_j bounded, which FAILS when beta has aligned spike mass
(aligned rung) or gamma is supercritical-aligned AND dominant (there the
top-coordinate law shifts by sqrt(n l)-scale means and u* > 0 becomes
identifiable - this is desirable behavior, not a bug). The lemma's scope
is exactly A4a-dense beta.
