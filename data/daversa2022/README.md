# Daversa et al. 2022: three-host Bd mesocosms

- **Files** (unmodified): `design.csv` (initial numbers per tank), `endnums.csv` (numbers alive
  at the end) and `Multihost_mesocosms_2017_infection.csv` (one row per sampled animal, with its
  Bd load in GE).
- **Source**: Dryad, https://doi.org/10.5061/dryad.2v6wwpznk, fetched from Dryad's Zenodo mirror
  (record 6456515) on 2026-11-12, because the Dryad file API needs a login token. The full
  deposit and the authors' R code (Zenodo, https://doi.org/10.5281/zenodo.5068208, MIT) are in
  the lab workspace: `lab/data/daversa2022/` and `lab/papers/daversa2022_si1.Rmd`.
- **Paper**: Daversa DR, Bosch J, Manica A, Garner TWJ, Fenton A (2022) *Am Nat* 200(4).
  https://doi.org/10.1086/720638
- **Licence**: CC0 1.0.

`load_daversa2022()` repeats the data preparation in Part 4 of the R code. Day 0 is 2017-07-04.
It drops the `accidental` samples, calls an animal infected when `bd.load > 0`, and computes
loss rates μ = ln(H(0)/H(T))/T.

`fig4cd_R0tot_digitised.csv` is not part of the deposit. It holds the per-tank posterior
medians of R0_TOT (black points), digitised on 2026-11-12 from the raster of Fig. 4d (page 11,
386 ppi; one pixel is 0.0086 in R0).

- Points were matched to tanks by their frog:toad ratio on the x-axis.
- Tanks 10 and 19 both have a ratio of 0.6 and medians of 1.046 and 1.033, which cannot be
  told apart, so both are given 1.040.
- The five 12:6:6 tanks (16, 20, Ri12, Ri18, SN) plot at the same height, 1.580.
