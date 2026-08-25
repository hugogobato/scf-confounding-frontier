# T3(b): Visibility Boundary for Quadratic Probes in (c, g)

Status: CLOSED 2026-08-25 (session 3). Closed-form phase curve derived,
validated against fresh simulation of the exact Phase-2 pipeline and
against results/lecam_probe_auc.csv (systematic structure); permanent
falsifier in tests/test_theory_T3b.py. Companion to
docs/theory_T3_scoped_impossibility.md (T3(a), eigenvalue-alarm side).

## Empirical driver

Phase-2 Le Cam probe (GBM over 14 frozen features): chance-level at
c = 0.8 up to g = 3.2 but AUC -> 1 at c = 0.2. Forensics (this session):
the ONLY informative feature is f12 = q(b) := ||b||^2 with
b = X_c' Ytilde / n, Ytilde = Y_centered / sd_y (self-standardized).
Under H1 the AUC of f12 moves BELOW 1/2 (0.073 at c=0.2,g=0.8): ||b||^2
DEFLATES - the sigma-inflation signature. Everything below explains the
magnitude, the sign, the c-ordering, and a genuine saturation ceiling.

## Setup and objects

M1 model, sigma_u = sigma_eps = 1 (general sigma_eps carried through),
gamma = g*dir(theta), omega := ||Lambda dir||^2 = sum_j l_j dir_j^2
(the per-unit-g confounding mass). Probe quantities:

    q_raw := ||X_c' Y_c / n||^2,      q := q_raw / sd_y^2,
    sd_y^2 := ||Y_c||^2/n,            v(g) := plim sd_y^2.

## Three exact deterministic limits (each verified to <= 2% at n = 2000)

(L1) H0 response variance: v0 := v(0) = beta'Sigma_X beta + sigma_eps^2
     -> 1 + sigma_eps^2 (beta Haar: beta'Sigma_X beta -> 1).
     Measured 1.9948 / 2.0114 at c = 0.2 / 0.8.

(L2) H0 raw floor: M0 := plim E[q_raw | H0] = 1 + c(1 + sigma_eps^2) + o(1).
     Derivation: E[q_raw] = E[beta' Sigmahat^2 beta] + E||X'eps/n||^2;
     Haar-beta averaging makes the first term a trace identity,
     E[beta'Sigmahat^2 beta] = tr(E[Sigmahat^2])/p = 1 + c + o(1)
     (from E[tr Sigmahat^2] = [(tr Sigma)^2 + (n+1) tr Sigma^2]/n), and
     E||X'eps/n||^2 = sigma_eps^2 tr(Sigmahat)/n = sigma_eps^2 c + o(1).
     Measured 1.3986 / 2.6168 vs 1.4 / 2.6.

(L3) H1 numerator gain: plim E[q_raw | gamma = g dir]
     = M0 + g^2 (omega + c) + o(g^2):
     the coherent piece ||Lambda gamma||^2/n-scale = omega g^2 PLUS the
     leakage ||U' (gamma f)||^2 / n^2 -> c g^2 (the factor vector projects
     on the noise part of the design with per-coordinate variance p/n).
     Hence the standardized mean

     E[q | g]  ->  [M0 + g^2 (omega + c)] / (v0 + g^2)
                =  m0 - g^2 kappa / (v0 + g^2),
     m0 := M0/v0 = 1/(1+sigma_eps^2) + c,   (measured 0.7005 / 1.3006
                                             vs 0.70 / 1.30 - exact)

## THE VISIBILITY CURVE

Define the headroom

    kappa(profile, theta) := m0 - (omega + c)
                           = omega_star - omega,
    omega_star := 1/(1 + sigma_eps^2)        [= 1/2 at sigma_eps = 1].

All c-dependence cancels: the boundary is governed ONLY by whether the
per-unit-g confounding mass omega (= dir'L L' dir) sits below or above the
fixed point omega_star. Mean-shift phase curve:

    delta(g) := E[q|g] - m0 = -g^2 (omega_star - omega)/(v0 + g^2),

with a SATURATION CEILING |delta|max = |omega_star - omega| as g -> inf:
quadratic probes cannot push the standardized floor past omega + c no
matter how strong the confounding.

VISIBILITY CONDITION (two-class map completing T3(a)):
Let A be the separation (in sd_0 units) required for the target AUC
(A ~ 2.5 for AUC ~ 0.96 under Gaussian q; Mann-Whitney AUC =
Phi(-|delta|/(sqrt(2) sd_0)) for a deflation-side signal).

  BLIND CLASS:   |omega - omega_star| <= A sd_0(q0)  =>
                 NO quadratic functional of (b, sd_y) detects at any g
                 (supremum of the standardized shift is bounded).
  VISIBLE CLASS: otherwise g_vis solves g^2(kappa - A sd_0) = A sd_0 v0:
                 g_vis = sqrt( A sd_0 v0 / (kappa - A sd_0) ).

Measured sd_0(q0) = 0.0331 (c=0.2) / 0.0422 (c=0.8); kappa = 0.276 / 0.053.
Ceiling ratios |delta|max/sd_0 = 8.3 (visible class, AUC -> 1 reached) and
1.26 (blind class, AUC stuck near 0.18-0.25 even at g = 3.2): matches the
frozen csv ordering g_vis(0.2) << g_vis(0.8) and the measured plateaus.

WHY THE c-ORDERING (the headline): for the scaling profiles used
everywhere in Phase 2-3, l_j proportional to sqrt(c), so
omega = w0 sqrt(c) with w0 = 0.5 (sub) or 2.375 (mixed) at theta = pi/6.
The sub-profile hits the fixed point omega = omega_star EXACTLY AT c = 1:
kappa_sub(c) = 0.5(1 - sqrt(c)), vanishing at the Marchenko-Pastur
boundary. The probe is born blind as c -> 1 from either side, and for
c > 1 the sign flips (omega > omega_star: ||b||^2 INFLATES instead).
The mixed profile crosses at c = 0.044, so it is visible (positive-signed)
across the whole tested grid - exactly the frozen csv pattern
(mixed AUCs > 1/2 and rising; sub AUCs < 1/2-side and stuck at c = 0.8).

Sign bookkeeping (matches forensics): sub at c < 1: deflation (f12-AUC
below 1/2); mixed: inflation (f12-AUC above 1/2); sub at c = 1: neither.

## Scope notes and guardrails honored

* sd_0(q0) DE: sd_0^2 = 2 tr(tilde_Sigma^2)/p^2
  + 4 sigma_eps^4 beta'Sigmatilde^3 beta/n + o(...), with tilde_Sigma the
  E[Sigmahat^2] operator (trace form 1 + c + o()); leading constant
  validated against measured 0.033/0.042 within ~25%; exact finite-n
  constant left as a refinement (does not affect the class map).
* Guardrail from the plan (do NOT substitute Gaussian-beta for Haar):
  honored - the derivation USES Haar averaging (beta draws fresh per rep
  in both the frozen pipeline and the falsifier); a Gaussian-beta variant
  would change only O(r/p) terms here, recorded, not assumed.
* The frozen csv small-g AUCs at c = 0.2 (e.g. 0.845 at g = 0.15) EXCEED
  the Gaussian mean-shift envelope (predicted ~0.53): the GBM also
  exploits distributional SHAPE beyond the first moment at low signal.
  The mean-shift curve is therefore the LOWER ENVELOPE / class-map
  theory; per-point AUC prediction at small g is explicitly out of scope.
  Recorded honestly rather than tuned away.

## Validation ledger

| Check | Status |
|-------|--------|
| v0 = 1 + sigma_eps^2 | PASS (1.9948 / 2.0114 vs 2) |
| M0 = 1 + c(1+sigma_eps^2) | PASS (1.3986 / 2.6168 vs 1.4 / 2.6) |
| m0 floor = 1/(1+sigma_eps^2) + c | PASS (0.7005 / 1.3006 vs 0.70 / 1.30) |
| kappa = omega_star - omega | PASS (impl. 0.258-0.277 / 0.050 vs pred 0.276 / 0.053) |
| Saturation ceiling E[q|g->inf] = omega + c | PASS (trend: 0.465 -> 0.424; 1.2546 -> 1.247) |
| Sign map (sub deflate / mixed inflate / c=1 null) | PASS (forensic AUC sides) |
| Class map vs frozen lecam_probe_auc.csv | PASS (ordering + plateaus; small-g envelope caveat above) |
| Permanent falsifier | tests/test_theory_T3b.py |
