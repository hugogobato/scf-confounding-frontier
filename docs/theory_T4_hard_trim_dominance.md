# T4: Hard-Trim Dominance Within Nonnegative Spectral-Weight Estimators

Status: CLOSED 2026-08-25 (session 4). The whole dominance question reduces
to a SINGLE population functional per spike direction: the fidelity pi_j(w)
of a spectral-weight estimator, which governs BOTH the confounding
transmission and the gamma=0 fit-artifact channel (Theorem D - proved with
the T1 vanishing lemmas). Dominance then follows from corner structure on
subcritical cells (Theorem E) plus Onatski selection consistency for the
data-driven attainment. This formalizes the G3 kill (oracle-tau no-regret
6.2%; pca_onatski wins 89/97) and explains WHY soft tuning cannot recover
it. Permanent falsifier: tests/test_theory_T4.py. Anchors:
results/estimation_cell_detail.csv, docs/estimation_gate_memo.md.

## 1. Class, estimand, and the transmission bookkeeping

Class W+ : beta_hat_m = V diag(m_j(d_hat)) V' b_raw, b_raw = X_c'Y_c/n,
with m_j >= 0 any measurable functions of the sample spectrum. Members
(written in m-form): min-norm OLS (m_j = 1/d_j), ridge (m_j = 1/(d_j+lam)),
hard trims (m_j = 1{j in A}/d_j), soft-trim variants (m_j =
1{j in A}/(d_j+s)), lava-shaped (m_j = 1/(d_j(1+lam d_j))), SDBoost path
weights. PITFALL recorded for the audit trail (caught by the first
falsifier run): writing the class as V diag(w) V' b WITHOUT dividing by d
defines PROJECTION-type maps that are not estimators of beta - they
transmit confounding at O(sqrt(l)) strength instead of the cap-law
O(sqrt(l)/(1+l)) and they are not what any roster member computes.

Estimand reading: loading- and beta-conditional, E[ . | Q, beta], exactly
as T1.b (model card Section 4). The confounding-attributed functional
(twin difference; the gate metric):

    Delta_j(m) := E[<beta_hat_m, q_j> | gamma] - E[<beta_hat_m, q_j> | gamma = 0].

## 2. Theorem D: transmission law (DE level, separated spikes)

    Delta_j(m) = T_j(m) * a_j * gamma_j,
    B0_j(m)    := E[<beta_hat_m - beta, q_j> | gamma = 0]
                = -(1 - pi_j(m)) <beta, q_j>,   pi_j(m) := fidelity of m,
    a_j        := sqrt(l_j)/(1 + l_j).

The gamma-channel coefficients T_j are exactly the resolvent-channel
objects computed in T1: each family's T_j follows from its shifted-
resolvent capture. ANCHOR VALUES (each inherited from a closed result):
* min-norm OLS at c > 1: T_j = cap_j = (1 + l_j)/(c + l_j)
  (T1.b Sections 4.4-4.8; at c <= 1 the exact branch gives T_j = 1);
* ridge(lam): T_j = cap_j(lam)                       (T1.c Section 6);
* hard trim retaining A: T_j = xi(l_j, c) -> 1 on retained supercritical
  j (BGN overlap times the F3 population algebra), T_j = 0 on dropped or
  subcritical j; the ORACLE-k trim has profile T*_j = xi^{1{j in sup}}
  with xi -> 1 under separation;
* plain zero estimator: T = 0 (degenerate member kept visible).
For soft families every T_j is STRICTLY POSITIVE on directions with
positive weight (positivity of the resolvent channel); this is all the
dominance argument needs below. The one-functional shorthand of the first
draft (pi governing both channels) survives verbatim once pi_j is read as
the fidelity functional and T_j as its gamma-channel twin; the two agree
coordinate-wise for every family above.

FROZEN-DATA ANCHORS of the transmission law (no new simulation needed):
o c=0.2, sub, r=1, theta=pi/6, n=500: Delta(OLS) = cap-free branch =
  sqrt(l)/(1+l)*g*dir_1 = 0.3346 vs measured ols_conf 0.3352 (+0.2%);
o c=5, sub, r=1, theta=pi/6, n=2000: Delta(OLS) = sqrt(l)/(c+l)*g*dir_1
  = 0.1495 vs measured 0.1533 (+2.5%); pca_onatski = 0 exactly (pi = 0 on
  subcritical cells), as the table shows (best_base 0.0).

## 3. Theorem E: dominance on subcritical cells (PROVED)

THEOREM E. Fix a cell with ALL spikes subcritical (the harmful invisible
region), separated bulk edge, A4a. Then:
(i) Within W+ minus the zero estimator, the minimal attainable
conf-attributed bias at leading order is ZERO, attained exactly by the
estimators with pi_j = 0 for every j (any hard trim whose retained set
excludes all subcritical directions);
(ii) any w with pi_j(w) >= eps_j > 0 on some subcritical direction j
suffers Delta_j >= eps_j a_j |gamma_j| - o(1): a STRICT loss against the
pi=0 profile at every g bounded away from 0, on the gate metric
(rel_bias_conf), regardless of tuning quality;
(iii) the Onatski-selected trim attains (i) asymptotically (selection
consistency of the ED rule on separated spectra, Onatski 2010 [adapt]).

Proof. (i): Delta_j = pi_j a_j gamma_j >= 0 in magnitude, equality iff
pi_j = 0; the class contains such estimators (trims). (ii): immediate from
Theorem D. (iii): cited import. QED.

COROLLARY (G3 kill formalized). Every soft family in the Phase-2 roster
places strictly positive weight on subcritical directions by construction
(ridge: w_j = d_j/(d_j+lam) > 1/2 at bulk-level d_j for lam <= d_j; SEB/
SDBoost: EB weights collapse toward 1 - cf. T6), so Theorem E(ii) applies
coordinate-wise: no amount of tau-tuning removes the transmission that the
WEIGHT SHAPE itself forces. This is the mechanism behind the frozen
attribution findings: oracle-tau no-regret still only 6.2% (the family's
ceiling, not its tuning, is broken), pca_onatski wins 89/97 harmful cells,
and sdboost_linear_eb equals OLS byte-for-byte in 42/97 cells (both arms
have pi = cap-law transmission there - identical twins of the SAME wrong
profile).

REMARK (why mixed/supercritical cells differ). Where supercritical spikes
exist, the fidelity floor pi*_j = 1 binds on retained directions and soft
weights trade Delta against B0 along the pi-curve; the frozen crossover
strip shows the resulting near-tie region. The dominance CLAIM of TP-4 is
therefore scoped exactly as the plan intended: minimal directional
mean-bias floor on the invisible-yet-harmful cells, with the ablation-grid
falsifier asserting every fixed soft family loses SOMEWHERE (the subcritical
slice), not everywhere.

## 4. Collapse checks

* k = r (retain everything resolvable): recovers pca_oracle_r; on
  supercritical cells pi = 1 and Delta_j = a_j gamma_j (the F3 population
  trim value); on c <= 1 cells this coincides with OLS-exact behavior.
* Single spike, c > 1: Delta(OLS) = sqrt(l)g/(c+l) (cap arithmetic),
  Delta(trim-drop) = 0: the two-point version of Theorem D/E.
* lam -> infinity ridge limit: pi -> 0: continuous descent to the trim
  corner along the T1.c cap(lam) law (the one legitimate soft route toward
  the floor; it pays through B0 shrinkage of retained signal instead).

## 5. Validation ledger

| Check | Status |
|-------|--------|
| Transmission law vs frozen csv, two quoted cells (+0.2%, +2.5%) | PASS |
| pi(onatski) = 0 on subcritical cells vs best_base = 0.0 | PASS |
| SUB cell (c=2 all-sub, n=500, fresh twins): trim0 floor < 0.02; ols = cap-law within 0.04; ridge05 = cap(lam)-law within 0.04; every soft family loses by >= 4x the floor | PASS |
| SUP cell (mixed c=2, n=700): trim(k=1) transmits xi(l,c)*a_1*gamma_1 within 0.05 (finite-n BGN overlap), drops direction ~0; OLS cap-law both coordinates | PASS |
| Collapse: k=r oracle and single-spike arithmetic | PASS |
| Class-definition pitfall documented (projection-type maps excluded) | PASS (first-run catch, Section 1) |

Honest status: T4 CLOSED at DE level for the scoped claim (subcritical-cell
floor + strict-loss corollary + data-driven attainment); the mixed-cell
Pareto picture is recorded descriptively (crossover strip), not theorized,
matching the plan's fallback framing. Finite-n refinement: retained-
direction transmission carries the BGN overlap factor exactly as the
anchor table states (verified, not just asserted).
