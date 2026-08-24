# Correctness Gate Memo (WP 2.1) - SKELETON, data pending

Status: pre-registered expectations frozen 2026-08-24 (before sweep launch).
Result tables fill from data/sim/correctness/ once shard archives return.
Verdicts are computed mechanically by code/gate_checks.py.

## Preregistered grid (amendment v2, D6/D7)

102 cells: mainA n=500 x 6c x 6 profiles @1000 reps; mainB n=2000 @200 reps;
alignS theta in {0, pi/2} r=5 profiles @150; deepN n=8000 skinny-c branch
(c in {0.1, 0.2}) @150. c = 10 dropped (D6a). Twins on c > 1 cells.

## Preregistered pass rule (plan Section 6 G3 row)

DE-vs-sim deviation <= 10% at n = 4000-equivalent in >= 90% of estimation
cells, shrinking with n. Null size (MC-calibrated) in [0.035, 0.065] in
>= 95% of null cells. Fail: systematic deviation > 25% in >= 30% of cells
after two fix iterations -> Phase 2 give-up rule 1.

## Predictions to overlay

- OLS mean-bias: exact identity at c < 1 (F1); capture law cap_j =
  (1+l_j)/(c+l_j) at c > 1 (pilot-validated conjecture F8). Deviations here
  measure MC noise and code correctness, not asymptotics, at c < 1.
- Ridge curves over LAM grid: population identity at c <= 1; provisional
  lambda-interpolated capture at c > 1 (flagged open item F8/T1; deviations
  there are expected and informative, NOT gate-blocking beyond the plan's
  two-fix-iteration allowance for missing finite-size terms).
- lambda_max vs BBP/bulk edge; subspace overlaps vs BGN sums (equal-spike
  profiles use subspace masses per inherited amendment 4).

## Known pre-data caveats carried into evaluation

- Erratum 1 (F12): spike-coordinate leakage scales d_j/c; affects S2
  calibration and SEB g^2 inversion only.
- Erratum 2: analytic S1 threshold miscalibrated at finite n (reported,
  not gating).
- Erratum 3: S1-aug demoted to diagnostic; subcritical blindness expected -
  nullcal cells should CONFIRM size-calibration of S2/S0 despite blindness
  to alternatives (size is a null property).

## Results

PENDING - populate after consolidation:
| cell class | median dev | max dev | n-shrink | verdict |
|---|---|---|---|---|

## Verdict

PENDING against the numeric conditions above.
