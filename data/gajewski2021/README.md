# Gajewski et al. 2021: Bd optical density at constant and varying temperatures

- **Paper**: Gajewski Z, Stevenson LA, Pike DA, Roznik EA, Alford RA, Johnson LR (2021).
  Predicting the growth of the amphibian chytrid fungus in varying temperature environments.
  *Ecol Evol* 11:17920–17931. https://doi.org/10.1002/ece3.8379 (CC BY 4.0). Full text from
  PMC8717292 and Appendix S1 (`ECE3-11-17920-s001.pdf`) from the Europe PMC supplementary
  files, retrieved 2026-11-19 and kept in `lab/papers/gajewski2021/`.
- **Source of the data and code**: Zenodo, https://doi.org/10.5281/zenodo.4525954
  (`gzach93/Predict_Bd_Growth_VaryEnvir` V1.0, GitHub commit c4ade0e), file
  `gzach93/Predict_Bd_Growth_VaryEnvir-V1.0.zip`, retrieved 2026-11-19 through the Zenodo REST
  API and extracted to
  `lab/papers/gajewski2021/zenodo/gzach93-Predict_Bd_Growth_VaryEnvir-c4ade0e/`.
- **Licence**: the Zenodo record gives the licence id `other-open` with no text, and the
  repository has no licence file. No CC0 or CC BY licence is attached, so **the data are not
  bundled here**. `gajewski2021_path()` finds the extracted folder in the lab workspace, or
  the folder named by `ENV["CHYTRIDIOMYCOSIS_GAJEWSKI2021"]`; tests that need it are skipped
  otherwise.
- **What the package reads**: `optical_QLD.csv` (Stevenson et al. 2013, doi:10.1371/
  journal.pone.0073830: one Queensland isolate at 13–28 °C, days 0–14; column `adjusted` is the
  OD the authors fitted) and `Vary_OD.csv` (Stevenson et al. 2020: the Dry Low Air, Wet High
  Frog and Wet Low Frog regimes, days 0–12). These are the same experiments as Stevenson et
  al.; do not count them again as independent data.
- **Not read**: `TemperatureData.csv` (the incubator loggers). The package uses the authors'
  piecewise fits to it from `functions.R`, which are in the compendium
  (`data/compendium/gajewski2021.toml`, sets `regime_*`).
