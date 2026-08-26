# Claim Map (WP 5.2)

Every empirical claim in `paper.tex` maps to the figure/artifact/config that regenerates it. Paths are relative to the repository root (`scf/`). Grids are hash-locked (`configs/grid_*.json`, hash a920e554e692); every result row carries the sha256 of model_card.md + assumption_ledger.md.

| Paper claim (section) | Figure / Table | Artifact(s) | Config source |
|---|---|---|---|
| Capture-law overlay accuracy: 91.67% of 60 cells <=10% at n=2000; 100% at n=500/8000; tier medians 3.73/1.03/0.68%; maxima 13.6/4.5/0.9%; failing-cell diagnosis (Sec 3.5) | Fig. fig:overlay | `results/correctness_overlays.csv`, `data/sim/correctness/correctness_results.parquet` | `configs/grid_correctness.json` |
| Ridge interpolation + superseded-form falsification (Sec 3, Remark) | - | `code/de_formulas.py::ridge_capture` / `ridge_capture_superseded`; reconciliation cells in `docs/theory_T1_capture_law.md` Sec 6 | - |
| Trimmed-tau attenuation: Lambda_D = 1.558 vs 1.5548+/-0.072; fresh-run mean -0.403; six-fold OLS inflation 0.064 -> ~0.40 (Sec 3.4) | - | `results/m2_treatment.csv`, `data/sim/m2/m2_results.parquet`, `docs/theory_T7_trimmed_tau.md` | `configs/grid_m2.json` |
| Alarm calibration: pooled chi2 p=0.4114, max|z|=2.91, 96.8% in +/-2 sigma; raw band share 80.65%; Bonferroni median size 0.1625; S0 analytic size 1.0 (Sec 4.2) | - | `results/null_sizes.csv`, `data/sim/nullcal/nullcal_results.parquet` | `configs/grid_nullcal.json` |
| Power frontier ratios 0.93/1.25/1.02 (median 1.02); worst cell c=0.8 ratio 1.246; 9/9 blind strata <=0.25 power at g=3.2 (Sec 4.3) | Fig. fig:power | `results/power_surface.csv`, `results/frontier_check.csv`, `data/sim/power/power_results.parquet` | `configs/grid_power.json` |
| Alignment stress: S2 power 1.000 across theta except 0.022 at pi/2; F baseline 0.262 (Sec 5.1) | Fig. fig:alignment | `results/alignment_stress.csv`, `data/sim/alignment/alignment_results.parquet` | alignment sweep configs |
| Detachment boundary ratios omega_e/B = 0.235/0.081/0.036; wake margins 0.64 vs 0.22, 1.42 vs 0.45, 1.27 vs 0.71 (Prop detach + App C) | - | `docs/theory_T3a_eigenvalue_contiguity.md`; falsifiers in `tests/test_theory_T3a.py` | - |
| Le Cam probe split: c=0.8 chance-level (AUC 0.50-0.58) through g=0.8, informative ~g>=1.6; c=0.2 AUC 0.84 at g=0.15 rising to 1.0; MMD deferred (Sec 5.2) | Fig. fig:lecam | `results/lecam_probe_auc.csv` | frozen probe spec in `docs/detection_statistics.md` |
| Visibility anchors: v0/M0/m0 measured vs theory; sd0 0.0331/0.0422; headrooms 0.276/0.053; ceiling ratios 8.3 vs 1.26 (Sec 5.2) | - | `docs/theory_T3_visibility_boundary.md` validation table | - |
| Soft-trim kill: no-regret 4.12% vs required 95%; oracle ceiling 6.2%; half-cut 40.21% vs 70%; Onatski best in 89/97; ridge 8; SDBoost==OLS 31/97 (Sec 6) | Figs fig:crossover, fig:phasediagram | `results/estimation_cell_detail.csv`, `results/crossover_curves.csv`, `data/sim/estimation/estimation_results.parquet`, `data/sim/crossover/crossover_results.parquet` | `configs/grid_estimation.json`, `configs/grid_crossover.json` |
| Robustness: 28 variants, ordering intact; r+/-1 median bias ratio 1.48; injection powers 0.28-0.44 (Sec 6, Sec 8) | - | `results/robustness_variants.csv`, `data/sim/robustness/robustness_results.parquet` | `configs/grid_robustness.json` |
| Benchmark premise PF-2: scree rejects all six unmodified designs and matched nulls; TW statistics up to 48430.8 (Sec 7) | Tab tab:pf1 context | `results/benchmark_arms.csv`, `data/benchmarks/spectral_audit.json` | `configs/benchmarks_frozen.yaml` |
| PF-1 straddle: power(0.5x g*) in [0.1425, 0.2050], power(2x g*) in [0.975, 0.9975]; control sizes 0.000-0.0033 (Sec 7) | Tab tab:pf1, Fig fig:bench | `results/bench_pf1.csv` | benchmark freeze `results/benchmark_freeze.json` |
| PF-3: tau trim 0.173 vs OLS 0.421 at 1x frontier; ratio 2.43; sign 100%; null cost 0.128 vs 0.110; ridge 0.238 (Sec 7) | - | `results/bench_m2.csv` | benchmark freeze |
| PF-4: strength proxy rejects H0 out-of-sample at 0.07-0.92; asymmetry score 0.61-1.00; partial-F 0.25-0.93; scree 1.000; alarm 0.000-0.003 (Sec 7) | Tab tab:calib | `results/calibration_informality.csv`, `results/bench_replication.csv` | benchmark freeze |
| G4 adjudication transparency: mechanical REVIEW (band fails from below only), adjudicated GO; zero anti-conservative rejections in 3000+ controls (Sec 7, Sec 8) | - | `results/bench_g4_verdict.json`, `results/bench_g4_adjudication.json` | - |
| Sensitivity: align_top power 0.05 vs spread 0.48 at same g*; hetero_eps 0.47-0.58; r_inj+1 powers 0.28-0.44 (Sec 7, Sec 8) | - | `results/bench_sensitivity.csv` | sensitivity battery |
| Replication: largest statistic shift 7.7e-13 on independent hardware (Sec 7, App D) | - | `results/bench_replication.csv` | replication check |

Theory claims map to permanent falsifiers: capture law -> `tests/test_theory_T1.py` (39 checks); trimmed law -> `tests/test_theory_T7.py`; visibility law -> `tests/test_theory_T3b.py` (T7+T3b suites total 48 checks); detachment boundary -> `tests/test_theory_T3a.py`; identity layer -> `tests/test_identities.py`.
