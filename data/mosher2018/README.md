# Mosher et al. 2018: boreal toad and Bd detection histories

- **Files**: the two deposited CSV files, unmodified.
  - `boreal_toad_bd_histories.csv`: encounter histories of boreal toads (*Anaxyrus boreas
    boreas*) and Bd at 83 sites in the southern Rocky Mountains, 2001–2010, with 10 survey
    slots per year.
  - `boreal_toad_bd_covariates.csv`: site covariates and the day of year of each survey.
- **Source**: Dryad, https://doi.org/10.5061/dryad.1c990 (version 2, published 2018-01-16).
  The files were downloaded by the user on 2026-11-20 (the Dryad file API needs a login),
  and the originals are kept in `lab/papers/mosher2018/dryad/`.
- **Paper**: Mosher BA, Bailey LL, Muths E, Huyvaert KP (2018). Host–pathogen metapopulation
  dynamics suggest high elevation refugia for boreal toads. *Ecol Appl* 28:926–937.
  https://doi.org/10.1002/eap.1699
- **Licence**: CC0 1.0. This is the licence on the Dryad landing page, as
  `"license": "https://spdx.org/licenses/CC0-1.0.html"` in its metadata.
- **Reader**: `load_mosher2018()` (`src/fitting/mosher2018.jl`).

## Layout

**`boreal_toad_bd_histories.csv`**

- **Row 1** gives the year, once every 20 columns (2001 … 2010).
- **Row 2** is the header: `Site`, then `Toad 1`, `Bd 1`, …, `Toad 10`, `Bd 10` for each
  year.
- **Then one row per site.**
  - Codes: `1` detected, `0` not detected, `.` no survey (for Bd, no swab).
  - Bd records exist only for 2003–2010, and only where a toad record exists.
  - One survey has `Toad 0`, `Bd 0`; with `pB = rBa = 0` this is equivalent to a toad
    non-detection.

**`boreal_toad_bd_covariates.csv`**

- **Site covariates**: `Elevation (m)`; `Permanence` (`permanent`/`semi-permanent`);
  `Other Amphibians` (`yes`/`no`, whether tiger salamanders or boreal chorus frogs were
  ever detected); `Region` (`NE`, `NW`, `SE`, `SW`, separated by I-70 and the Continental
  Divide).
- **`jdate1` … `jdate100`**: the day of year of survey `j` of year `y` is in column
  `10(y − 1) + j`, and `0.5` means no survey. Site 5 has a date at survey 2 of 2008 with no
  record; the reader ignores it.

## Site 33 and the paper's 82 sites

The paper analysed 82 sites, and the deposit has 83. Site 33 is the only site with no
permanence or other-amphibian value, and toads were never detected there. Dropping it gives
494 site-years with at least one survey. That is the sample size implied by the AICc and
deviance values of the paper's Table 3. `load_mosher2018()` therefore drops it by default;
`drop_incomplete = false` keeps it.
