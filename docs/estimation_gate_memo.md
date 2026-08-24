# Estimation Gate Memo (WP 2.2) - FINAL

Status: data complete and consolidated 2026-08-24 (99/99 cells, 16,520/16,520
reps; crossover 15/15; robustness 28/28). Verdicts from
results/gate_verdicts.json (WP22_estimation), cell detail in
results/estimation_cell_detail.csv and results/crossover_curves.csv;
figures: figures/crossover_curves.pdf, figures/bias_phase_diagram_empirical.pdf.

## Verdict: KILL for eb_spectral as a competitor claim (plan give-up rule 2)

The preregistered no-regret condition fails decisively: bias(eb_spectral) is
within 1.05x of the best baseline in only 4.1% of the 97 harmful cells
(need >= 95%). The strongest baseline is pca_onatski, which wins 89/97 cells;
ridge_cv wins the remaining 8; pca_baing, cevid_default and sdboost_linear_eb
never win. The >= 50%-vs-OLS condition also fails at 40.2% (need >= 70%),
while the weaker >= 10%-cut condition holds in 82.5% of harmful cells, so
give-up rule 3 (incremental-only) does not apply: the family helps materially
over doing nothing but never over hard-trim PCA.

Mechanically this triggers plan Phase 2 give-up rule 2 (loses to the strongest
baseline in the majority of harmful cells). The honest conclusion is that
under A4a dense-beta dense-confounding DGPs with well-separated spikes,
Onatski-selected hard trimming is near-oracle: its mean-bias floor (~0.08
across the entire g-grid, visible in the crossover strip) sits below
everything else, and soft weights can only add bias.

## Attribution findings (predeclared descriptive checks)

1. The dominance is structural to the soft-trim family, not a tuning
   failure: with oracle tau (eb_oracle_tau, the family's ceiling) no-regret
   still holds in only 6.2% of harmful cells vs 4.1% for SEB-tuned. When
   spikes are well separated, zeroing the sub-spike coordinates exactly is
   feasible and optimal; smoothly down-weighting them cannot compete.
   Within the family, tuning quality still matters where the family has a
   chance: at theta = pi/2 subcritical-spike slices eb_cv_tau achieves
   0.18 vs eb_spectral's 0.44 mean-bias norm (c = 0.8, n = 2000, r = 25),
   showing the DE-predicted-MSE objective mis-ranks taus precisely where
   the supercritical-aligned signal vanishes.
2. sdboost_linear_eb collapses exactly to OLS: its mean-bias vector equals
   OLS's to < 1e-9 in 42/97 harmful cells (its EB variance components drive
   all spectral weights to 1 on these DGPs). This is an independent negative
   result about Nava et al.'s linear special case worth reporting.
3. The practical-effect threshold (>= 15% relative MSE reduction vs best
   baseline in the mixed region) is met in 0% of evaluated mixed cells,
   consistent with the kill verdict.

## What survives

The estimation contribution pivots from "beat baselines" to (i) the DE
correctness layer (exact overlays, PASS per correctness memo), (ii) the
detection frontier and decoupling results (detection memo), and (iii) a
characterization of when soft-trim tuning beats or loses to hard trimming
(the theta-pi/2 sub-slice attribution above feeds T1/T2 theory questions).

## Deviation register additions

D8 applies (id reconciliation, see correctness memo). Twin-difference
(rel_bias_conf primary metric) used wherever twins exist (all c > 1 cells);
raw mean-bias norms elsewhere, matching the frozen priority list.
