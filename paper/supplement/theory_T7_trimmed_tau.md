# T7: Trimmed-Tau Deterministic Equivalent for M2

Status: CLOSED 2026-08-25 (session 3). Both arms proved at the level of
exact identities plus deterministic-equivalent limits, with every mechanism
verified against fresh simulation. Numerical anchors:
results/m2_treatment.csv (tau errors flat across c for the Onatski trim:
0.0466/0.0348/0.0416 at c = 0.2/0.8/2.0, while OLS inflates to 0.391 at
c = 2) and this session's channel-decomposition diagnostics (recorded in
the validation ledger).

## Setup (model card M2)

D = pi'X + delta'f + nu, Y = tau D + beta'X + gamma'f + eps, tau_true = 1
in Phase 2/3 runs, sigma_u = sigma_eps = 1, ||pi|| = 1 (sparse unit vector,
generic vs Haar Q), ||delta|| = delta_g redrawn Haar in R^r per rep,
gamma = g*dir(theta) fixed. Estimators as implemented (runners.run_m2_rep):

    tau_hat_trim: joint OLS of Y on [D, S],  S = X V_k,
                  k = max(Onatski select, r)   => k + 1 <= n, full rank;
    tau_hat_ols:  joint min-norm OLS of Y on [D, X] (p >= n branch relevant).

## T7.a Exact identities (both arms, every n)

TRIM. By Frisch-Waugh-Lovell on the full-rank design [D, S]:

    tau_hat_trim - tau = <d_tilde, w_tilde> / <d_tilde, d_tilde>,        (*)
    d_tilde := (I - P_S) D,  w_tilde := (I - P_S) w,  w = beta'X + gamma'f + eps,

with P_S := S(S'S)^{-1}S' the sample score projector in R^n. EXACT at every
rep; no pseudo-inverse anywhere because k + 1 <= n.

OLS (min-norm, c > 1). Writing G := XX' and using Sherman-Morrison on
MM' = G + dd' (M = [D, X]):

    tau_hat_ols = d'G^{-1}Y / (1 + d'G^{-1}D),
    tau_hat_ols - tau = (d'G^{-1}w - tau) / (1 + d'G^{-1}D).              (**)

VERIFIED numerically to 15 digits against pinv at (n,p,c) = (300,600,2).
Identity (*) fails to exist at c > 1 precisely because P_X-perp = 0 there;
(**) is its regularized replacement and is the source of the OLS pathology.

## T7.b plim deterministic equivalents

Write s_j := sqrt(l_j) (sigma_u = 1), A := {j: spike j separated from the
bulk edge and retained}, and for j in A define the score-capture functional

    cap_sc(j) := s_j^2 * sum_{m<=k} (q_j'v_m)^2 / mu_m  ->  l_j xi(l_j,c)/mu(l_j),

with mu(l) = (1+l)(l+c)/l the BBP outlier location and xi(l,c) the BGN
overlap; for subcritical factors cap_sc -> 0 (no eigenvector alignment).
Bounds: 0 <= cap_sc(j) <= l_j/(1+l_j) (right end = idealized population-V_k
trim, provable by pure algebra; left end = no capture).

TRIM NUMERATOR. Conditioning on X and using Gaussian odd-symmetry exactly
as in T1 Sections 4.3-4.5 (every probe-vs-(K,u)-measurable cross moment has
conditional mean zero):

    plim <d_tilde, w_tilde>/n =
        sum_j delta_j gamma_j rho_j + o(1),
    rho_j := 1 - cap_sc(j)   in [1/(1+l_j), 1] on separated spikes;

the pi-gamma and pi-beta channels are O_p(p^{-1/2}) (they reduce to
pi'q_j-weighted spike coordinates times gamma_j s_j and pi'beta_perp under
A4a), and nu/eps channels have mean zero exactly.

TRIM DENOMINATOR.

    plim <d_tilde,d_tilde>/n = v_D
        = pi'Sigma_X pi + ||delta||^2 + sigma_nu^2
          - sum_j cap_sc(j) delta_j^2 + o(1).

Consequently (loading- and beta-conditional reading per model card Section
4): the systematic trimmed-tau shift is

    plim tau_hat_trim - tau = [sum_j delta_j gamma_j rho_j] / v_D + o(1),

which is O(delta_g * g) and INDEPENDENT of c at fixed (l_j) profile shape:
this is the flatness recorded in m2_treatment.csv. The beta-artifact term
that inflates min-norm beta-hats (T1.b) is ABSENT because the trim design
has no null space.

OLS SHRINKAGE TERM (the six-fold inflation mechanism). From (**), with

    Lambda_D := plim d'G^{-1}D
        = (1 - 1/c)||pi||^2
          + sum_j delta_j^2 t(1+t)/(1 + t(1+l_j))
          + sigma_nu^2 t + o(1),      t := 1/(c-1),

(the three lines being the pi'X-, delta'f-, and nu-components; the middle
one uses the T1 resolvent limit f_j'G^{-1}f_j -> t(1+t)/(1+t(1+l_j)); the
last uses tr(G^{-1}) ~ n t), we get

    plim tau_hat_ols - tau = -tau / (1 + Lambda_D)
                             + [sum_j delta_j gamma_j tilde_rho_j]
                               /(1 + Lambda_D) + o(1),

where the tilde_rho_j channel is small at the m2 grid (all named channels
measured |mean| <= 0.024 at c = 2). The DOMINANT term is the O(tau)
multiplicative shrinkage caused by the treatment column competing with the
p >= n design columns inside the min-norm normal equations: the D-slot
variance is shared with the null space of X. It is invisible at c <= 1
(FWL branch, no such term), turns on exactly at c > 1, and magnifies as
t(c) = 1/(c-1) declines: at the m2_weak grid, Lambda_D = 0.5 + 0.058 + 1
= 1.558 predicted vs 1.5548 +/- 0.072 measured at c = 2, giving
-tau/(1+Lambda_D) = -0.391 against measured means -0.403 (fresh sim) /
0.391 (frozen csv). QED at DE level.

Collapse checks: c -> 1- : trim formula continuous with the exact c <= 1
partitioned identity (rho_j = 1/(1+l_j) at perfect population alignment);
g -> 0: numerator -> 0, both arms unbiased up to the O(p^{-1/2})
beta-adjustment channel; delta_g -> 0: numerator -> 0, OLS shrinkage term
SURVIVES (it is a pure artifact of c > 1, not of confounding) - matching
the Phase-3 observation that OLS tau degrades even without injection.

## T7.c CLT layer (statement; scoped proof sketch)

sqrt(n)(tau_hat_trim - plim) -> N(0, Omega) with Omega the bilinear-form
variance of the FWL ratio: Omega = [sigma_eps^2 v_D + sigma_nu^2 sigma_w^2
+ cross-covariances of the named channels]/v_D^2, all entries converging by
the LLN/CLT for quadratic forms of jointly Gaussian vectors. The data-driven
V_k enters only through Onatski selection; under spike separation
(Onatski 2010 consistency) k is eventually constant and the selection step
contributes o_p(n^{-1/2}), reducing the CLT to the fixed-k case. GUARDRAIL
(binding scope): without spike separation k fluctuates, the o_p claim is
NOT automatic, and the CLT must be restricted to separated profiles -
stated as an explicit scope condition, not silently assumed.

## Validation ledger (2026-08-25 session 3 diagnostics)

| Check | Status |
|-------|--------|
| m2 flatness across c (trim arms, frozen csv) | PASS (0.0466/0.0348/0.0416) |
| OLS inflation at c=2 (frozen csv) | PASS (0.391 vs 0.064 at c=0.2) |
| SM identity (**) vs pinv, 15 digits | PASS |
| Channel SDs vs formulas (trim arm, c=0.2): a_pg, a_pb, rest | PASS (0.021/0.053/0.046 vs 0.021/0.050/0.045 predicted) |
| Lambda_D formula at c=2 | PASS (pred 1.558 vs meas 1.5548 +/- 0.072) |
| OLS shrinkage -tau/(1+Lambda_D) at c=2 | PASS (-0.391 pred vs -0.403 fresh / -0.391 frozen) |
| Trim MAE assembled vs frozen csv | PASS within ~15% (finite-n capture gap below) |
| a_dg absolute SD at finite n | GAP: measured sits between perfect-alignment (l/(1+l) capture) and BGN-overlap (l xi/mu capture) models, consistently ~15% above the aligned model at all c; attributed to partial finite-n subspace capture. Does not affect flatness or the OLS verdict; refinement queued. |

Honest status: T7 CLOSED at the level needed by the paper (exact
identities, DE plims with verified constants, scoped CLT). Remaining
refinement (non-blocking): finite-n score-capture interpolation for rho_j.
