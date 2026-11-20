# Hardy et al. 2023: boreal toad capture histories, Colorado, 1998-2020

- **Files**: the deposited files, unmodified: `Site_k_CMR.csv` (capture histories) and
  `Site_k_ENV.csv` (site covariates) for k = 1…5, and the deposit's `README_Hardy_etal.txt`.
- **Source**: Dryad, https://doi.org/10.5061/dryad.fttdz08wv (version 2, published
  2022-10-07). Retrieved 2026-11-19.
- **Paper**: Hardy BM, Muths E, Lambert BA, Schneider SC, Funk WC, Bailey LL (2023).
  Compensatory recruitment unlikely in high-elevation amphibian populations challenged with
  disease. *J Appl Ecol* 60:400–411. https://doi.org/10.1111/1365-2664.14310
- **Licence**: CC0 1.0 (public domain dedication).

## Sites

The files are not named by site. They are matched to Table 1 by elevation and the year of
Bd arrival:

| file | site | elevation (m) | Bd from | toads |
|---|---|---|---|---|
| 1 | Delilah Pond | 2755 | 2015 | 86 |
| 2 | Collegiate Camp | 2951 | 2015 | 266 |
| 3 | Denny Creek | 3238 | 2016 | 271 |
| 4 | South Cottonwood | 3183 | 2014 | 796 |
| 5 | Morgan's Gulch | 3366 | 2014 | 175 |

## Columns

**`_CMR`**

- Column 1 is a row number (or the PIT tag at site 1).
- Then there are 69 columns `YY_v` or `YYYY_v`: years 1998–2020, three visits each. `1`
  means caught, `0` not caught and `.` not sampled.
- Next come five group indicators, `End` (`;`) and `Concatenated`. The last holds the PIT tag
  in a comment and the Program MARK history string. `load_hardy2023` checks that string
  against the 69 columns.

All rows of a site share the same sampling pattern.

**`_ENV`**

- `Year` runs from 1996 to 2020.
- Drought is the number of days without precipitation in the active season. Winter onset is
  the day of year of the first sustained killing frost. Both are `N/A` after 2017.
- `Bd` is 0 or 1, and `Elevation` is in m.

The headers of `Site_1_ENV.csv` and `Site_2_ENV.csv` say "Drought, Winter Start", but those
columns hold the winter-onset day (about 270–325) and then the drought days (about 60–140).
That is the reverse of the headers and of the other three files. `load_hardy2023` assigns
the columns by value. The appendix −2 log L of f(DROUGHT) and f(WINSTART) are reproduced
only this way.
