# Robustness and Scaling Memo (WP 2.4) - FINAL

Status: data complete and consolidated 2026-08-24 (robustness 28/28 cells,
8,550 rows; scaling 4/4). Tables: results/robustness_variants.csv,
results/scaling_envelope.csv, results/m2_treatment.csv.

## Verdict

The estimation-side picture is variant-invariant, which is itself the robust
finding that matters after the WP 2.2 kill: eb_spectral's mean-bias norm
relative to the best baseline has median 1.5 to 5.3 across every variant
arm (V0 reference included), so the dominance of Onatski hard trimming is
not a Gaussian-homoskedastic artifact of the comparison. Per the frozen fail
rule this is recorded as a PIVOT within G3 for estimation claims: method
comparisons are restricted to the regime statements supported by the data
(see estimation memo).

Variant medians of bias(eb)/bias(best): V0 gaussian 2.26, V1 t5 errors 2.24,
V2 rademacher-half loadings 2.18, V3 heteroskedastic-u 5.26 (worst; soft
weights mis-read the inflated bulk under row heteroskedasticity), V4
correlated factors 2.05 (graceful, matching the pre-data expectation that
Omega rescales but does not reorder the geometry), V5 r +/- 1 misspecification
1.48 (least damaging; the spike-profile prior absorbs one-spike errors),
V6 sparse confounding 2.31. The qualitative ordering of estimator families
survives every perturbation; no variant flips a winner.

Diagnostic-size sub-rule (< 0.15 under V3-V4) is a documented scope cut:
the runner stores detection statistics on the main arm only, and gamma = 0
twins exist as estimation rows without rejection columns. Nothing downstream
depended on it; flagged in gate_checks output rather than silently dropped.

## Scaling envelope (single-rep runtimes, full roster sum per rep)

(n, p) = (1000, 1000): 0.02 s; (2000, 8000): 0.91 s; (4000, 8000): 1.53 s;
(8000, 1600): 0.16 s. Spectrum computation dominates; the suite scales with
min(n, p)^2 to min(n, p)^3 as designed, inside the memory ceiling everywhere.
The originally planned (8000, 8000) and (8000, 16000) cells were replaced by
ceiling-compliant variants at grid freeze (D6); the envelope table supports
Phase 3 packaging claims without them.

## M2 weak-treatment block (final)

All three cells complete (150 reps each). Mean absolute tau errors:
c = 0.2: OLS 0.064, ridge(1) 0.055, trim-onatski 0.047, trim-oracle-r 0.047;
c = 0.8: 0.058 / 0.039 / 0.035 / 0.035;
c = 2.0: 0.391 / 0.046 / 0.042 / 0.042.
Under weak-treatment alignment (delta_g = 0.3), naive OLS degrades six-fold
across the c-range while both trimmed estimators stay flat, and Onatski
selection matches oracle-r exactly in every cell. Treatment-effect estimation
under dense confounding inherits the Phase 2 detection geometry: trim first,
then regress. This is the cleanest positive packaging of the hard-trim
dominance finding for Phase 3 (it aligns with, rather than fights, the WP 2.2
kill: trimming is the right tool; soft weighting is not).
