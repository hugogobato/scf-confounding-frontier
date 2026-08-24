# Phase 2 Execution Memo (gate G3 infrastructure)

Date: 2026-08-24. Scope: Phase 2 build-out per research plan Section 7
(Phase 2), launched on the Phase 1 GO verdict. Everything in this memo
happened BEFORE any Phase 2 sweep data was generated; all thresholds and
grids are frozen with a deviation register (D1-D7) in
docs/phase2_preregistration.md.

## What exists now

1. Frozen pre-registration: docs/phase2_preregistration.md (grids v2 after
   compute-driven amendments D6/D7, estimator roster, thresholds, give-up
   rules, memory ceiling).
2. Detection specification + errata: docs/detection_statistics.md. The three
   errata are pre-data findings from the verification loop:
   - E1: F12 null-variance law corrected to (sigma_eps2 + d_j/c)/sigma_y2;
     the draft had the rowspace factor inverted.
   - E2: v1 analytic S1 threshold miscalibrated at finite n (measured
     sd/width ~ 8.5); MC-calibrated thresholds are the gate statistics,
     as predeclared in the original freeze.
   - E3: S1-aug demoted to diagnostic (tracks design spectrum, no
     discrimination beyond S0+S2 in probed cells); subcritical dense
     confounding empirically invisible to the whole statistic class at
     pilot scale (t_maxz/t_aug unchanged at g = 2 while bias is O(1)) -
     this CERTIFIES the invisible-yet-harmful region operationally and
     promotes "is any second-moment test possible below the BBP threshold?"
     to a real theory question (T2/T3 input).
3. Code (all unit-tested; pytest 35 green):
   - code/de_formulas.py: F12 (corrected), F13 aug-secular machinery,
     bbp_invert, MP quantiles, noise-scale estimation shared by detection
     and tuning, weight families (trim/lava), SDBoost path, SEB objective,
     UCM proxy.
   - code/estimators.py: full frozen roster incl. faithful Cevid Trim
     default (median singular value) and SDBoost linear special case
     (EB variance components on the XX' spectrum + boosting path +
     BLUP-corrected CV stopping, shared-subspace fold approximation D4),
     our SEB soft-trim estimator, ablations, oracle.
   - code/detection.py: S0/S1/S2/S3/B1 + probe features + Le Cam probe
     (GBM + MMD) + MC threshold helper.
   - code/simulator.py: extended DGP surface (sparse/rademacher loadings,
     t5 errors, heteroskedastic u, correlated factors, r-misspecification,
     M2 treatment block), backward-compatible with Phase 1.
   - code/runners.py + run_jobs: 9 sweep modes, checkpoint/resume per cell,
     worker caps, single-thread BLAS.
   - code/build_grids.py -> configs/grid_*.json (hash a920e554e692).
4. Compute plan: data/sim/shard_manifest.json + notebooks/colab_shard_01..38
   (172 projected wall-hours at 2 workers/notebook, all <= 5.5 h target,
   resume-safe, download-fallback included). Budget respected: 38 <= 40.

## How to run Phase 2 (operator instructions)

0. Code distribution: the public repository
   github.com/hugogobato/scf-confounding-frontier carries the harness; each
   notebook clones it at the pinned tag `phase2-freeze` and verifies the
   five harness files against generation-time sha256 prefixes before
   running (mismatch aborts). If the harness ever changes: commit, move/
   re-create the tag, re-run code/make_shards.py, and redistribute
   notebooks - never mix notebooks with an unpinned checkout.
1. Upload notebooks colab_shard_XX.ipynb to Colab accounts (any grouping of
   the 38; each is self-contained given internet access for the clone).
   Run them; download the produced scf_shard_XX.zip archives.
2. Drop all archives into one folder, then:
       python3 code/consolidate_shards.py <folder-with-zips>
   It verifies checksums against embedded manifests, ingests payloads into
   data/sim/<sweep>/raw|means/, writes consolidated parquets +
   completeness.csv per sweep and data/sim/consolidation_report.json.
3. Re-run any incomplete shards (resume-safe: completed cells are skipped).
4. Gate evaluation:
       python3 code/gate_checks.py
   followed by filling the four memos' Results sections
   (docs/{correctness,estimation_gate,detection_gate,robustness}_memo.md).
5. Figures regenerate exclusively from parquet by the figure scripts
   (to be added with the analysis pass; preregistered list in the
   pre-registration document).

## Design decisions recorded this session (beyond the errata)

- Estimation metrics consume the RAW centered response (D3): standardizing Y
  rescales fitted coefficients by 1/sigma_y and would contaminate bias
  ratios. Detection statistics standardize internally (scale-invariant).
- ridge_cv and eb_cv_tau use exact LOO identities over their grids (any
  diagonal spectral fit has hat-diagonal h_ii = sum_j w_j u_ij^2), replacing
  fold-wise eigendecompositions; verified against brute-force refits.
- SDBoost stopping uses BLUP-corrected K-fold CV with fold random-effect
  operators built from the full-data left-singular subspace (Woodbury);
  documented approximation, cost-motivated.
- n = 8000 cells are affordable only on the skinny-c branch (c <= 0.2)
  because the spectrum rides the p-side Gram there; at c >= 0.8 a single
  rep costs minutes (eigh(8000^3)). This drove D6b/D6c.
- c = 10 dropped project-wide (D6a); c = 5 retains the p >> n regime where
  the capture-law story lives.

## Open items carried forward

- T1: derivation of cap_j = (1+l_j)/(c+l_j) and its ridge-lambda
  interpolation (blocks nothing; provisional curve flagged in overlays).
- T2/T3 (promoted by E3): is ANY second-moment test detectable below the
  BBP threshold given sigma-normalization? The nullcal/power/probe data will
  quantify the empirical side.
- Analysis-phase artifacts: overlay/figure scripts, MC-threshold analysis,
  probe training runs (all operate on returned parquet only).

## FINALIZATION RECORD (2026-08-24, post-data)

### Consolidation and integrity

All 38 shard archives returned and ingested: sha256 verification passed for
every payload file against embedded manifests; completeness 8 sweeps at
100% (504 cells total, ~131k reps); the single incomplete cell was m2's wide
twin arm (F-1 below). Zone.Identifier artifacts from the Windows download
hop were removed before ingestion. Deviations D8-D13 were registered at
consolidation/analysis time (see preregistration register).

Post-freeze harness fixes applied on main (tag `phase2-freeze` remains the
exact code the Colab shards ran):
- F-1: run_m2_rep added an identity matrix sized min(n,p) to a p x p Gram,
  crashing every p > n M2 cell (this is why shard 08 produced no wide-m2
  output; not a wall-limit issue). Fixed, cell rerun locally under the
  plan's local-execution allowance (~57 core-minutes), consolidated; all
  nine sweeps now report 100% completeness (504/504 cells, ~131k reps).
- F-2: run_cell now creates output directories itself (Colab shards call it
  without run_jobs' pre-creation). Caught by the end-to-end clone test
  before distribution.
- F-3: consolidate_shards.py made layout-tolerant (archives store members
  relative to data/, manifest keys /content-relative).

### Analysis pipeline

code/phase2_analysis.py computes all gate verdicts mechanically from
consolidated parquet into results/gate_verdicts.json plus per-topic CSVs;
code/make_phase2_figures.py renders the six preregistered figures from those
CSVs (regeneration order: consolidate -> analyze -> figures).

### G3 verdict summary (full detail in the four gate memos)

1. WP 2.1 correctness: PASS. OLS mean-bias DE overlay within 10% in 91.7%
   of n = 2000 cells, 100% at n in {500, 8000}; ridge overlays exact.
   The capture law (c > 1) and exact identity (c < 1) both hold at MC noise.
2. WP 2.2 estimation: KILL for eb_spectral as competitor (give-up rule 2):
   no-regret in 4.1% of harmful cells vs >= 95% required; pca_onatski wins
   89/97. Attribution: the SEB tuner loses the worst slices while eb_cv_tau
   does not; sdboost_linear_eb collapses exactly to OLS in 42/97 cells.
   Contribution pivots to correctness + detection layers.
3. WP 2.3 detection: size PASS noise-aware (chi2 p = 0.41, max|z| = 2.91;
   raw arithmetic FAIL co-reported per D10); power frontier PASS (3/3
   supercritical strata within factor 1.5, median ratio 1.02); blind region
   confirmed 9/9 strata with the sharpest alignment signature (S2 power
   1.000 everywhere except 0.022 exactly at theta = pi/2); Le Cam probe
   scopes the invisibility claim by c (chance-level probes at c = 0.8 up to
   g ~ 1.6, but AUC up to 1.0 at c = 0.2 via the ||b||^2 sigma-inflation
   signature). S1 analytic falsification rule fired as pre-declared; D5
   consequence stands. T2/T3 question upgraded: characterize the
   detectability phase diagram of second-moment functionals in (c, g).
4. WP 2.4 robustness/scaling/M2: estimation verdict variant-invariant
   (PIVOT recorded); scaling envelope table delivered; diagnostic-size
   sub-rule scope cut (D-note); M2 finalized after F-1 rerun.

Phase 2 closes with G3 = GO on the correctness + detection frontier program
(the paper's spine), NO-GO on soft-trim tuning as an estimator contribution,
and a materially sharpened theory agenda (T1 capture-law derivation now has
exact full-scale validation; T2/T3 has empirical phase structure to explain).
