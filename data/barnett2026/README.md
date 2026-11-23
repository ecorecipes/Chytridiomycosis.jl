# Barnett et al. 2026 (bioRxiv preprint): imperfect Bd prophylaxis

- **Paper**: Barnett KM et al. (2026) From protection to amplification: imperfect
  chytridiomycosis prophylaxis increases infections in wild amphibians. bioRxiv **preprint**
  (not peer reviewed), https://doi.org/10.64898/2026.05.15.725113 (v2, CC BY-NC-ND 4.0).
- **Source of the code and data**: https://github.com/kmbarnett100/protection_to_amplification,
  commit 29f4af6 (2026-08-12), retrieved 2026-11-23 and kept in
  `lab/papers/barnett2026/github/`: the NetLogo 6.3.0 model, the BehaviorSpace output tables,
  the field and challenge data, and the R code.
- **Licence**: the repository has no licence file, so **its files are not bundled here**.
  `barnett2026_path()` finds the folder in the lab workspace, or the folder named by
  `ENV["CHYTRIDIOMYCOSIS_BARNETT2026"]`; the field-trial tests and validation targets are
  skipped otherwise.
- `reproduction.toml` is written by `scripts/barnett2026_abm.jl`. For each of the authors' 12
  output files it holds summary statistics (means and medians) of their runs and of the
  package's runs with the same parameters, the z score of the difference in means, and the
  Kolmogorov–Smirnov distance. It contains no individual run.
- **Older and newer files**: the two-way and adverse files come from an older version of the
  NetLogo reporters (columns `metas_prev`, `avg_spn_inten_metas`), in which intensity and
  prevalence are 0 when fewer than two metamorphs carry zoosporangia; the newer ones use −100.
  The figure code in the repository reads two files under other names
  (`anti_transmission_vs_coverage.csv`) and has column names of the older version for it.
