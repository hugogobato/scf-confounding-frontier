# Data Feasibility Audit (WP 3.1, Phase 3)

Status: COMPLETE (2026-08-24). Read-only audit of candidate benchmark
families per the research plan; preprocessing is specified here and frozen in
`configs/benchmarks_frozen.yaml` before any comparative result (WP 3.2).
Machine-readable source: `data/benchmarks/spectral_audit.json`.

## Families selected

| ID | Family | Provider / license | Access | n | p (processed) | Role |
|----|--------|--------------------|--------|---|---------------|------|
| A | AddNeuroMed blood gene expression, batches 1-2 | NCBI GEO (public submission), GSE63060 + GSE63061, Illumina HumanHT-12 v3/v4 | series matrix files, sha256 pinned under `data/benchmarks/raw/SHA256SUMS` | 717 pooled (329 + 388) | 2000 (A_main) / 1800 (A_sub) | M1-style design with real batch factors; PF-1/PF-2 primary |
| B | IHDP covariate benchmark (Hill 2011 npci covariates) | public benchmark mirror (github.com/gpeng9/ihdp-causality), research use | CSV snapshot `data/benchmarks/ihdp.csv` | 747 | 750 (B_main) / 150 (B_wide) | high-dimensional featurized causal design; injected dense confounder |
| C | wooldridge::k401ksubs household cross-section | Rdatasets mirror (CC0), Vincent Arel-Bundock | CSV snapshot `data/benchmarks/k401ksubs.csv` | 800 seeded subsample of 9275 | 1600 (C_main) / 160 (C_wide) | M2 treatment block (e401k -> nettfa); PF-3 |

Rejected candidates: GSE46861 (n = 12, too small for the validated DE
regime); GSE54275 (multi-platform pooling confounds the batch semantics we
want to exploit); ACIC mirrors (no stable unauthenticated copy found at audit
time; IHDP covers the same methodological slot).

## Frozen preprocessing (hash-locked in WP 3.2; summary)

Family A: intersect ILMN probe IDs across the two series; average duplicate
probe IDs within a series; drop probes with missing values; log2(x+1) iff the
global max exceeds 25 (values are LUMI-processed and already on log scale;
rule never fires on these files); z-score each probe on the POOLED sample so
batch mean shifts are preserved by design; keep the top-P probes by pooled
variance. A_main: P = 2000, all samples. A_sub: single-batch sensitivity
design (GSE63061 only), P = 1800. Metadata kept: status (AD / MCI /
control), age, gender.

Families B/C: standardize continuous covariates (binary kept as 0/1), then a
seeded random Fourier feature map Z = sqrt(2/P) cos(X W + b), bandwidth by
median heuristic on a seeded subsample, W fixed by seed; RFF block z-scored
per pooled column (restores the model-card scale sigma_u = O(1); affine
rescaling only, spectral shape unchanged). B: 25 raw covariates -> P in
{750, 150}. C: controls {inc, incsq, agesq, age, male, marr, fsize} on an
800-row seeded subsample -> P in {1600, 160}; e401k is the treatment,
nettfa the outcome; p401k and pira are EXCLUDED as post-treatment variables.

## Empirical spectral profiles (the frontier is benchmark-specific)

Noise-floor convention F-BENCH (frozen): kernel-featurized designs decay
smoothly and admit no MP-consistent bulk, so for ALL benchmarks
sigma_u^2_bench = max(q25(d), 1e-3 mean(d)) and spike strengths follow the
exact algebraic decomposition Sigma_X = se2 I + sum_j l_j se2 q_j q_j',
l_hat_j = d_j / se2 - 1. The shared MP bulk estimator is co-reported as
se2_mp_bulk_est for comparison (it agrees to within 2x on family A and
collapses on families B/C, which is exactly why the uniform convention is
needed).

| Config | n | p | c | r_hat_onatski | ktop_alarm = r_inj | top tau_j (d/se2) | TW stat of lam_max vs white noise |
|--------|---|---|---|---------------|--------------------|-------------------|-----------------------------------|
| A_main | 717 | 2000 | 2.79 | 2 | 2 | 141937, 1757, 781, ... | 48430.8 |
| A_sub  | 388 | 1800 | 4.64 | 1 | 1 | 846, 316, 236, ...   | 155.4 |
| B_main | 747 | 750  | 1.00 | 0 | 1 | 11050, 8977, 8147, ... | 1444.8 |
| B_wide | 747 | 150  | 0.20 | 0 | 1 | 81, 70, 58, ...       | 308.8 |
| C_main | 800 | 1600 | 2.00 | 0 | 1 | 142846, 138585, ...   | 8360.0 |
| C_wide | 800 | 160  | 0.20 | 0 | 1 | 31666, 29443, ...     | 1191.0 |

Reading: every UNMODIFIED design is massively outlier-positive against the
white-noise TW99 threshold (PF-2's premise quantified: scree inspection
rejects "white" on all six designs before any confounding is injected). The
dominant spikes are processing/batch structure (family A), income-age kernel
structure (families B/C). Injection geometry: Lambda uses the first r_inj =
ktop_alarm sample eigenvector directions with per-column sd
sqrt(d_j - se2); default gamma direction places equal mass on those r_inj
coordinates.

## Known results reproduced on unmodified data (plan Section 9.7)

1. k401k: OLS coefficient of e401k on nettfa given the control set is
   +9.41 (outcome mean 19.07, n = 9275): 401(k) ELIGIBILITY associates with
   HIGHER net financial assets, the Poterba-Venti-Wise sign reproduced by
   plain covariate adjustment.
2. IHDP (this mirror's realization): unadjusted treated-control difference
   4.02; OLS-adjusted coefficient 3.93; experimental-benchmark ATE ~ 4.0
   (Hill 2011). Point estimates land at the benchmark scale with adjustment
   slightly below it.
3. AddNeuroMed: batch structure is the dominant spectral feature
   (lam_max TW statistic 48431 pooled; 155 within batch 2 alone); diagnosis
   (status) metadata retained for label-based sanity checks but no
   published effect size is claimed or tested.

## Pass rule check

Plan requires >= 2 families with n_eff/p compatible with the theory
(p/n in [0.1, 10], clear factor structure). PASS: three families, six
configs, c in [0.20, 4.64], all with strong realized spike structure and
explicit batch/treatment semantics. Compute footprint: eigendecompositions
are computed ONCE per config (designs are fixed across injections), so every
later rep costs O(np); no Colab sharding needed (all cells < 2 min local).
