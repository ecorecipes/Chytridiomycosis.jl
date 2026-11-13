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
