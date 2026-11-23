# Drawert et al. 2017: sampled parameter points and the authors' stochastic summaries

- **File**: `baseline_points.csv`, the 2132 parameter points of the Baseline model, with the
  authors' 30-run stochastic summaries for each point (no intervention).
- **Source**: figshare, "Simulation Code from Using stochastic epidemiological models to
  evaluate conservation strategies for endangered amphibians",
  https://doi.org/10.6084/m9.figshare.5328127.v1 (The Royal Society, published 2017-08-21),
  file `rsif20170480_si_002.gz` (MD5 48829907e9c66c9cae0587be3d1c5a24), member
  `SI_simulation_code_and_data/baseline_params.csv` (SHA-256
  4cb20874cf8ba90a91eab50a88a7fe0159370015c35a4dfceac0d28e54e88c1b). Retrieved 2026-11-13.
  A copy of the archive, extracted, is in `lab/papers/drawert2017_si/`.
- **Paper**: Drawert B, Griesemer M, Petzold LR, Briggs CJ (2017) Using stochastic
  epidemiological models to evaluate conservation strategies for endangered amphibians.
  *J R Soc Interface* 14:20170480. https://doi.org/10.1098/rsif.2017.0480
- **Licence**: CC BY 4.0 (the figshare record). Attribution: Drawert et al. 2017.

## What was changed

Only columns were dropped; every kept value is the original string. The dropped columns are
constant across all rows (`param_nu1`, `param_sigma1`, `param_phi`, `param_Ar`, `param_Ac`,
`param_rho`, `param_initZoospores` = 0; `param_initSporangia` = 100; `run_tag` = baseline;
`param_ensemble_size` = 30; `status` = COMPLETE; all intervention columns = 0), always
`NULL` (`data_s_star_value`, `simulation_data`) or bookkeeping (`last_updated`, all
2013-11-26). The extraction keeps the columns below in this order.

## Fields

| Field | Meaning |
|:------|:--------|
| `ID` | point id, 1–2132 (ESM Table S9 cites ids 120, 1094, 182, 982, 192) |
| `param_Nfrogs` | initial number of frogs `N` |
| `param_gamma` | zoospore–frog encounter rate `γ` (frog⁻¹ day⁻¹; pool volume 1) |
| `param_eta` | zoospores released per sporangium per day `η` |
| `param_nu` | encystment probability `ν0` |
| `param_f` | self-reinfection fraction `f` |
| `param_sigma` | sporangium loss rate `σ0` (day⁻¹) |
| `param_mu` | zoospore death rate `μ` (day⁻¹) |
| `data_meanDaysToSmax`, `data_stdDaysToSmax` | first output day on which fewer than `N` frogs are alive (1000 if never) |
| `data_meanFracSurviving`, `data_stdFracSurviving` | fraction of the `N` frogs alive at the end of a run |
| `data_meanDaysToZ0`, `data_stdDaysToZ0` | day of fungal extinction (1000 if none) |
| `data_meanDaysToF0`, `data_stdDaysToF0` | day on which both frogs and zoospores are gone (1000 if never) |
| `data_fracPersistence`, `data_fracFrogsExtinct`, `data_fracFungusExtinct` | fraction of the 30 runs with each outcome |
| `data_lambda` | linearised Bd growth rate λ (day⁻¹), the dominant eigenvalue of the Baseline model |
| `simulation_time` | CPU seconds the authors' C code took for the 30 runs |

Standard deviations are population sds (divisor 30), as in `print_stats()` of
`WithinHostModel.c`. Outputs were taken once a day, so the day columns are whole days.

## Checks

- The table reproduces the paper's Baseline numbers (§4): 65.0% of points end in frog
  extinction in all 30 runs (paper "65%"), 24.4% in fungal extinction ("24%"), 0.14% in
  persistence ("about 0.01%") and 10.5% are variable ("approximately 10%"). The groups of
  Fig. 3 have 1386 (A, 65.0%), 109 (B), 221 (C), 396 (D) and 20 (E) points.
- `data_lambda` equals `load_growth_rate(DrawertParameters(:drawert2017; point = id))` to
  the 6 significant figures printed, and lies in 0.0351–0.3499 (the ESM's 0.035–0.35).
- Under the distance of ESM eq. (10), with the range of log c2 and log c3 as scale, no two
  points are closer than 0.0901 (ε = 0.09).

## Reproduction output (ours, not the authors')

Written by `scripts/drawert2017_*.jl`; see the script headers for designs and seeds.

- `baseline_sample.csv`: our 30-run summaries at a stratified sample of 116 points.
- `variants_sample.csv`: Nu and Sigma models, 10 runs per point.
- `interventions_sample.csv`: one-time interventions, 30 runs per arm.
- `reproduction.toml`: the summary read by the tests and the validation targets
  (`scripts/drawert2017_summarise.jl`).
- `deterministic.csv` (checks on all 2132 points) is not committed. It is 380 kB, and
  `scripts/drawert2017_deterministic.jl` regenerates it in a few minutes.

## Outcome table for decision analysis (ours, not the authors')

Written by `scripts/drawert2017_decision_table.jl` (design and options in its header); read with
`load_drawert_decisions()`, whose docstring defines the columns.

- `decision_runs.csv`: one row per run, 100 draws × 10 actions × 30 replicates of the Baseline
  model. `Inf` in `first_loss` means that no frog was lost. `uninfected` and `infected` count live
  frogs without and with sporangia (load threshold 0, `--threshold` to change it) at the end of the
  run and on days 120 and 365; culled frogs are in neither.
- `decision_runs_seconds.csv`: the wall time of each run, in the same row order. It is apart
  because it is not reproducible; the other files are, with the same Julia version.
- `decision_draws.csv`: the parameters of the draws, and `weight`, the number of the 2132 points
  in the draw's group per sampled draw. The draws oversample the smaller groups, so use the
  weights to describe all 2132 points. The points are a space-filling design, not a posterior.
- `decision_design.toml`: Julia version, actions, variants, replicates and rows.

**Seeds.** `drawert_run_seed(draw_id, replicate)`, the same for every action: common random
numbers up to the first intervention, one stream per run, no per-frog streams.

**Pulses that are not applied.** As in the authors' code, a pulse is skipped once the frogs or Bd
have gone (`applied` is false). By group, the share skipped is, on days 14, 28 and 56: A 0%, 5%,
61%; B 0%, 0%, 1%; C 0%, 11%, 26%; D 0%, 43%, 72%; E 0%. In Groups C and D (Bd has usually died
out) this makes a late cull look much less harmful than it would be for a manager who did not
know that Bd had gone: in those runs no frog is culled, and `handled` is 0.

**Cost.** 0.78 s per run over 50 random points, 10 actions and one replicate (7.8 s per point;
median 3.2, maximum 46.5; 6 threads on an Apple M-series laptop). It varies by group: A 0.87, B
0.18, C 0.02, D 0.003 and E 4.7 s per run. The Nu variant takes 8.7 s per run on the
same design (87 s per point; median 18, maximum 1073), because surviving runs go on to day 1000
with large pools: about 73 CPU-hours for 100 draws. The Sigma variant takes 11.5 s per run
(115 s per point; median 20, maximum 1041), or about 95 CPU-hours; other jobs were running, so that
figure is an upper estimate. The Nu and Sigma tables are not stored; make them with `--variants=nu,sigma`.
