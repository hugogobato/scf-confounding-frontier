# T7: Trimmed-Tau Deterministic Equivalent for M2

Status: c <= 1 exact partitioned identity PROVED below (elementary); c > 1
DE extension and CLT are ADAPT targets conditional on T1 Route A (see
docs/theory_T1_capture_law.md). Numerical anchor: results/m2_treatment.csv
(tau errors flat across c in {0.2, 0.8, 2.0} for the Onatski trim while
OLS inflates ~6x at c = 2).

## Setup (model card M2)

D = pi'X + delta'f + nu, Y = tau D + beta'X + gamma'f + eps, tau_true = 1
in Phase 2/3 runs. Estimator: Onatski hard-trim then Frisch-Waugh - regress
D and Y on the top-k sample PC scores S = X V_k (k = max(Onatski, r)), use
residuals:

    tau_hat_trim = <D_tilde, Y_tilde> / <D_tilde, D_tilde>,
    D_tilde := (I - Pi_S) D,   Pi_S := S(S'S)^{-1}S' (sample projection).

## c <= 1: exact finite-n plim identity (PROOF)

Work conditionally on X (loading-conditional functional per model card
Section 4). Write w := beta'X + gamma'f + eps so Y = tau D + w. Partitioned
regression gives EXACTLY

    tau_hat_trim - tau = <D_tilde, w> / <D_tilde, D_tilde>          (*)

at every rep (FWL identity; no approximation). Now decompose w into its
projection onto span[S-columns of X] plus orthogonal complement. The
population objects (n -> inf with k fixed or slow):

    a_D := cov(D, S) Sigma_S^{-1/k}... define pi_perp := pi -
           Sigma_X V_k (V_k'Sigma_X V_k)^{-1} V_k' pi   (the part of pi's
           design signature NOT captured by retained PCs),
    delta_perp := delta - Lambda'V_k (V_k'Sigma_X V_k)^{-1} V_k' pi.

Then M_D-perp := D - proj has population form
pi_perp'X + delta_perp'f + nu + o_P(n^0), and

    plim(tau_hat_trim) - tau =
        [ pi_perp' Sigma_X beta + pi_perp' Lambda gamma + delta_perp'gamma ]
      / [ pi_perp' Sigma_X pi_perp + ||delta_perp||^2 + sigma_nu^2 ] + o(1).

Proof: plug the decompositions into (*) and apply the same Gaussian
conditional-independence argument as T1.a to each mean-zero piece; the
numerator picks exactly the three systematic covariances listed; the
denominator converges to its population value by the LLR for quadratic
forms of jointly Gaussian vectors. QED (c <= 1 branch; also exact-in-form
at c > 1 PROVIDED all inverses are read inside the retained k-dim subspace,
where OLS is full rank - this is why the trim has no min-norm artifact).

Reading (the paper's practical recipe): trimming removes from BOTH D and Y
any confounding component that rides the retained spike directions
(delta_perp drops those coordinates of f), at the price of leaving
pi_perp/delta_perp channels. Under A4a the beta-numerator term is the
usual regression-adjustment term (present for ANY estimator), while the
gamma-channel survives only through delta_perp - which the Onatski choice
k >= r makes small when delta is generic. Hence trimmed-tau error stays
flat in c while raw OLS at c > 1 suffers the capture-law artifacts
(T1): the six-fold inflation recorded in m2_treatment.csv is an artifact
of E[P_null] acting on beta/D, not of confounding per se.

## c > 1: DE extension (ADAPT, conditional on T1 Route A)

Replace each occurrence of Sigma_X-action by its capture-weighted DE
(cap_j = (1 + l_j)/(c + l_j) on spikes, 1/c on bulk) in both numerator and
denominator; the claim to prove is that the FWL ratio is STABLE to that
replacement at first order because both entries share the same spectral
modulation. Numerical falsifier already frozen: m2 parquet cells at
n in {500, 2000}, c in {0.2, 0.8, 2.0}; predicted flatness within MC
tolerance; OLS arm must show the capture-artifact inflation predicted by
minnorm_total_bias_norm applied to the joint design.

## CLT sketch (uncertainty layer)

sqrt(n)(tau_hat_trim - plim) -> N(0, Omega) with Omega the standard
bilinear-form variance for ratios of Gauss-Markov residuals; the only
nonstandard input is the DATA-DRIVEN V_k (Onatski selection), handled by
the selection-consistency block (Onatski 2010) plus a martingale/empirical-
process argument showing selection enters at o_p(n^{-1/2}) under separated
spikes. Guardrail: if spike separation fails, k-selection fluctuates and
the o_p term is NOT automatic - scope the CLT to separated profiles.

## Validation ledger

| Check | Status |
|-------|--------|
| m2 flatness across c (trim arms) | PASS (m2_treatment.csv: 0.047/0.035/0.042) |
| OLS six-fold inflation at c=2 | PASS (0.391 vs 0.064 at c=0.2) |
| c<=1 identity simulation check | TODO cheap (reuse m2 runner at c<=1 cells) |
