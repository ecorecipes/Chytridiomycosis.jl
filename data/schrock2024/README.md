# Schrock et al. 2025: Bd aggregation group summaries (NOT BUNDLED YET)

- **Paper**: Schrock SAR, Walsman JC, DeMarchi J, et al. (2025) Do fungi look like
  macroparasites? Quantifying the patterns and mechanisms of aggregation for host–fungal
  parasite relationships. *Proc R Soc B* 292:20242013. https://doi.org/10.1098/rspb.2024.2013
  (CC BY 4.0). Preprint: bioRxiv https://doi.org/10.1101/2024.08.29.609018 (2024), the
  `schrock2024` key of `lab/papers/references.tsv`.
- **Data**: Dryad https://doi.org/10.5061/dryad.mgqnk998g, version of 2025-01-09, licence
  **CC0 1.0** (from the Dryad API record, 2026-11-14). Files: `analysis_aggregation_dataset.csv`
  (858,877 bytes, SHA-256
  `2f6468d581d936755ad18675d9e8ecbb794aff3d6e4ac92e06a7244f16951a3e`), `abund_and_prev.csv`
  and `README.md`. The raw swab loads are not deposited (the authors' grouping script says
  "Cannot provide raw data"); the CSV holds one row per group with the mean and variance of
  the load with and without zeros, the distribution fits and the aggregation metrics.
- **Code**: Zenodo https://doi.org/10.5281/zenodo.14291806 (MIT), read (not run) on
  2026-11-14: `01_Aggregation_Groups.R` builds the group summaries, `03_Fig.1_TPL_Plot.R`
  fits the per-dataset Taylor's power law slopes.

## Why the file is not here

On 2026-11-14 the Dryad file endpoints refused the package's automated client (HTTP 401 on
the API download, 403 on the web download), so the file could not be retrieved and checked.
It is CC0, so it may be bundled once someone downloads it by hand.

## To add it

1. Download `analysis_aggregation_dataset.csv` from the Dryad page above.
2. Check its SHA-256 against the value above, and save it here unmodified as
   `data/schrock2024/analysis_aggregation_dataset.csv`.
3. `load_schrock2024()` then reads it, and
   `schrock_taylor_slopes(load_schrock2024())` gives the per-dataset slopes. The pending
   validation targets `schrock2024_tpl_sierra`, `schrock2024_tpl_brazil` and
   `schrock2024_tpl_adult` can then be given calls, for example
   `schrock_taylor_slopes(load_schrock2024())[:sierra_nevada].slope` (published 1.90).

## Columns used

`dataset` (`brazil`, `east_bay`, `sierra_nevada`, `serdp` = eastern US), `species`,
`life_stage`, `site`, `N`, `num_infected`, `mean_natural_with_zeros`,
`var_natural_with_zeros` and `phase` (Sierra only). Missing values are `NA`,
`Not Applicable` or `Not Available`.
