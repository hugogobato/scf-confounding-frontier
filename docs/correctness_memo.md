# Correctness Gate Memo (WP 2.1) - FINAL

Status: data complete and consolidated 2026-08-24 (38/38 shards returned,
sha256-verified; completeness 102/102 cells, 33,300/33,300 reps). Verdicts
computed mechanically by `code/phase2_analysis.py` (results/gate_verdicts.json,
results/correctness_overlays.csv); figure: figures/de_overlay_grid.pdf.

## Verdict: PASS

DE-vs-sim deviation of the OLS mean-bias norm is at most 10% in 91.7% of the
n = 2000 cells (55/60), 100% of n = 500 cells and 100% of the n = 8000 tier;
medians by tier are 1.0% / 3.7% / 0.7%. The plan's wording "n = 4000-equivalent"
is interpreted (D9 below) as the largest tier with full grid coverage. The
give-up rule 1 (> 25% systematic deviation in >= 30% of cells) does not come
close to firing.

Ridge overlays over the frozen LAM grid are exact for practical purposes:
median relative deviation 0.0% in every tier (the c <= 1 population identity
and the capture-interpolation at c > 1 both land inside MC noise).
lambda_max matches the BBP/MP top edge to within a fraction of a percent in
all cells where stored. Subspace-overlap functionals (overlay iv) were not
stored by the runner schema and are recorded as scope cut D11; nothing in
the gate depends on them.

The non-monotone median across tiers (n = 2000 slightly above n = 500)
reflects profile composition, not deterioration: the n = 500 tier carries
more cheap low-c cells while n = 2000 includes the harder theta-variant and
r = 25 slices; every tier passes comfortably.

## Deviation register additions

D8 (post-consolidation, bookkeeping only): the config_id column frozen inside
configs/grid_*.json was computed by a pre-freeze Config.cid implementation.
Authoritative ids are recomputed from each config dict under the pinned tag
(the same ids that seed every rep); configs/cid_remap.json records the stale-
to-fresh mapping for audit. No config dict changed; all runs are internally
consistent and reproducible from the pinned code.

D9 (interpretation, declared at analysis time): "n = 4000-equivalent" reads
as the largest full-coverage tier (n = 2000), with the n = 8000 skinny-c
tier reported alongside (100% pass there).

D10 (evaluation refinement): the null-size gate is evaluated both as frozen
raw arithmetic and noise-aware (see detection gate memo); D6c rep cuts made
the raw band comparable to split-half MC error at 500-rep cells.

D11 (scope cut): per-vector/subspace overlap functionals were not recorded
by the correctness runner; DE overlay rests on mean-bias norms, ridge curves,
and lambda_max, which is what the plan's numeric conditions quantify anyway.

## Data provenance

Rows carry ledger hash 11b162ac814d (model card + assumption ledger pair);
figures regenerate exclusively from consolidated parquet via
code/phase2_analysis.py then code/make_phase2_figures.py.
