# Schrock et al. 2025: Bd aggregation group summaries

- **Paper**: Schrock SAR, Walsman JC, DeMarchi J, et al. (2025) Do fungi look like
  macroparasites? Quantifying the patterns and mechanisms of aggregation for host–fungal
  parasite relationships. *Proc R Soc B* 292:20242013. https://doi.org/10.1098/rspb.2024.2013
  (CC BY 4.0). Preprint: bioRxiv https://doi.org/10.1101/2024.08.29.609018 (2024), the
  `schrock2024` key of `lab/papers/references.tsv`.
- **Data**: Dryad https://doi.org/10.5061/dryad.mgqnk998g, version of 2025-01-09, licence
  **CC0 1.0** (from the Dryad API record, 2026-11-14).
- **Code**: Zenodo https://doi.org/10.5281/zenodo.14291806 (MIT), read (not run) on
  2026-11-14: `01_Aggregation_Groups.R` builds the group summaries, `03_Fig.1_TPL_Plot.R`
  fits the per-dataset Taylor's power law slopes, `04_TPL_Life_Stage.R` the life-stage
  slopes and `05_Fig.2_Distribution_Fit_Boxplot.R` the distribution comparison. No script
  for the epizoological-phase slopes is in the deposit.

## File

`analysis_aggregation_dataset.csv` (858,877 bytes, SHA-256
`2f6468d581d936755ad18675d9e8ecbb794aff3d6e4ac92e06a7244f16951a3e`), bundled unmodified.
The Dryad file endpoints refused the package's automated client on 2026-11-14 (HTTP 401/403);
the file was downloaded by hand from the Dryad page on 2026-11-15 and its SHA-256 checked
against the value recorded on 2026-11-14 from the Dryad API. The other files of the deposit
(`abund_and_prev.csv`, `README.md`) are not bundled; the hand-downloaded copies are in
`lab/papers/schrock2024_si/`.

The raw swab loads are not deposited (the authors' grouping script says "Cannot provide raw
data"); the CSV holds one row per group (1,861 groups) with the mean and variance of the load
with and without zeros, the authors' distribution fits and AICs, and aggregation metrics.

## Columns used

- `load_schrock2024`: `dataset` (`brazil`, `east_bay`, `sierra_nevada`, `serdp` = eastern US),
  `species`, `life_stage`, `site`, `N`, `num_infected`, `mean_natural_with_zeros`,
  `var_natural_with_zeros` and `phase` (Sierra only).
- `schrock_family_fractions`: `num_infected`, `gamma_aic`, `lnorm_aic`, `weibull_aic`,
  `exp_aic`.

Missing values are `NA`, `Not Applicable` or `Not Available`. Quoted fields can contain line
breaks (647 rows; the Sierra phase labels, e.g. `"enzootic\n(Northern Sierra)"`), which the
loader handles.
