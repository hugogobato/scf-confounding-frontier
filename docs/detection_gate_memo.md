# Detection Gate Memo (WP 2.3) - SKELETON, data pending

Status: pre-registered expectations frozen 2026-08-24 (before sweep launch),
amended by Errata 1-3 of docs/detection_statistics.md (all pre-data).

## Frozen statistics and decision rules

Gate statistics: S2 maxz_cal with MC-calibrated thresholds from matched
g = 0 nulls (mc_thresholds); S0 scree TW99 as the practitioner-visibility
anchor; B1 partial F; UCM proxy on headline cells. S1 aug_bbp co-recorded as
diagnostic only (Erratum 3). Le Cam probe: GBM + median-heuristic MMD on the
frozen 14-feature map; computational undetectability declared at AUC <= 0.55
for BOTH probes.

## Pre-data empirical expectations (from unit-test-scale diagnostics)

1. Supercritical-aligned cells: S2 fires strongly once calibrated
   (measured t_maxz 1.09 -> 5.85 at g=2, n=600 c=0.4 cell).
2. Subcritical cells at g = 2: S2 AND S1-aug blind (t_maxz 1.20 -> 1.24);
   ||b||^2 can SHRINK under strong subcritical confounding because sigma_y
   grows faster than the cross-moment mean - the invisible-yet-harmful
   region is expected to certify as undetectable-by-this-class.
3. The F12 law predicts the calibrated power frontier for S2 including the
   d_j/c leakage tax: power should DECREASE in c for fixed supercritical l
   (leakage grows), a counterintuitive, falsifiable signature.

## Preregistered pass/fail

- Size gate: MC-calibrated size in [0.035, 0.065] in >= 95% of null cells.
- Power gate: S2 power >= 0.8 at s <= 1.5x predicted frontier; empirical
  power-1/2 contour within factor 1.5 of the F12-based prediction.
- Probe gate: AUC <= 0.55 below the claimed frontier (undetectability side).
- Fail (plan give-up rule 4): uncalibratable size (> 2 alpha everywhere after
  robust variants) -> KILL detection component, keep C1/C3.
- Fail: empirical frontier off by > factor 3 -> ansatz wrong for the test.

## Alignment stress expectation

Size theta-invariant by construction (gamma = 0); S2 power decays with
theta as the supercritical projection shrinks; B1 decays slower; the S2-vs-B1
contrast is the operational decoupling evidence.

## Results

PENDING - populate from data/sim/nullcal/, power/, alignment/ + probe runs
(probe training happens at analysis time from stored features f0-f13).

## Verdict

PENDING.
