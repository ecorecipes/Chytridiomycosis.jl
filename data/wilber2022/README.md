# Wilber et al. 2022: host density and Bd outcomes across Sierra Nevada lakes

- **Paper**: Wilber MQ, Knapp RA, Smith TC, Briggs CJ (2022). Host density has limited effects
  on pathogen invasion, disease-induced declines and within-host infection dynamics across a
  landscape of disease. *J Anim Ecol* 91:2451–2464. https://doi.org/10.1111/1365-2656.13823
- **Paper copy used**: the accepted manuscript from the NSF Public Access Repository,
  https://par.nsf.gov/servlets/purl/10401712 (retrieved 2026-11-14; 30 pages, without the
  supporting information). Kept as `lab/papers/wilber2022_si/nsfpar/accepted.pdf`.
  Appendices S1–S4 were not obtained: the Wiley supplement links returned HTTP 403, and the
  PMC author manuscript (PMC12594071, which lists `NIHMS2113987-supplement-Supplement.pdf`) is
  behind a browser check.
- **Source of the data and code**: Zenodo, https://doi.org/10.5281/zenodo.7114634
  (`mqwilber/host_density_and_bd` v1.0.0, GitHub repository id 496634407, release of
  2022-09-26), file `mqwilber/host_density_and_bd-v1.0.0.zip`. Retrieved 2026-11-14 through the
  Zenodo REST API. The zip and its extracted files are kept in
  `lab/papers/wilber2022_si/zenodo7114634/`.
- **Licence**: the Zenodo record gives the licence id `other-open` with no text, and the
  repository has no licence file. No CC0 or CC BY licence is attached, so **the data are not
  bundled here**.

## Using the deposit

[`wilber2022_path()`](../../src/models/wilber2022.jl) finds the extracted deposit:

1. `ENV["CHYTRIDIOMYCOSIS_WILBER2022"]`, if set, should point at the extracted root
   (the directory holding `code/`, `data/` and `results/`); otherwise
2. `papers/wilber2022_si/zenodo7114634/extracted/mqwilber-host_density_and_bd-4903d73` under
   any ancestor of the package directory (the lab workspace layout).

The tests and validation targets that need the data are skipped when it is absent.

Files read (all from the deposit's `data/`):

| file | used for |
|---|---|
| `formatted/combined_Bd_data_for_analysis_allVES.csv` | swabs (Bd load per swab) and visual-encounter counts per lake and year; Hypothesis 1 |
| `candidate_datasets/analysisI_failedinvasions_excluded.csv` | the lakes of Hypothesis 1 (≥ 3 non-missing years: 79 lakes) |
| `environmental_covariates/summer_lake_temperatures.csv` | mean June–September maximum air temperature, °C (gridMET), per lake and year |
| `environmental_covariates/cdec_snowpack_april1.csv` | 1 April snow water equivalent, mm (nearest CDEC station), per lake and year |
| `candidate_datasets/analysisII_magnitudeofdecline_with{,out}subadults.csv` | Hypothesis 2 (decline) |
| `candidate_datasets/analysisIII_loadanddensity.csv` | Hypothesis 3 (enzootic path analysis) |

Not in the deposit, so not reproducible: `data/archival/sites_nov2020.csv` (lake areas, needed
for the density covariate of Hypothesis 1), `data/archival/humphreys_translocation_2010.csv`
(the translocation of Hypothesis 3), and the fitted Stan objects (`results/*.rds`). The deposit
does hold `results/looic_model_comparison_results.csv` (the LOOIC values behind Table 1) and
`results/effect_sizes_across_rhos_rarifyFALSE.jpg` (Fig. 3 at higher resolution, digitised for
`data/compendium/wilber2022.toml`).

No code from the deposit was run; it was read to reimplement the analyses.
