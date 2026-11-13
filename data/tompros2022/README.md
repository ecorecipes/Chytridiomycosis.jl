# Tompros et al. 2022: Bsal transmission in eastern newt mesocosms

- **Paper**: Tompros A, Dean AD, Fenton A, Wilber MQ, Carter ED, Gray MJ (2022).
  Frequency-dependent transmission of *Batrachochytrium salamandrivorans* in eastern newts.
  *Transbound Emerg Dis* 69:731–741. https://doi.org/10.1111/tbed.14043
- **Source**: the paper's data availability statement points to the University of Tennessee
  TRACE repository, https://dx.doi.org/10.7290/xw8yy2RBeu
  (https://trace.tennessee.edu/handle/20.500.14382/16899). The item's single file is
  `utk.ir.fg_2341.zip`, holding the data, the R code and a README. Retrieved 2026-11-12
  through the repository's REST API. The zip, the TRACE `license.txt` and the extracted files
  are kept in `lab/papers/tompros2022_trace/`.
- **Licence**: the paper says the data "are openly available". The deposit README gives
  "Licenses/restrictions placed on the data: NA". The TRACE deposit agreement grants the
  University of Tennessee non-exclusive rights to distribute the item. No open licence (such
  as CC BY) is attached to the item. The article itself is CC BY 4.0. These files are
  redistributed here for reproduction, with attribution: Tompros et al. 2021, TRACE,
  doi:10.7290/xw8yy2RBeu.

## Files

| file | origin |
|---|---|
| `counts_days0to5.csv` | `SIDat.Trans` in the deposit's `DataFile2_Days0to5.RData`, written as CSV |
| `counts_day8.csv` | `SIDat.Day8.NoDeaths` in `DataFile3_Day8_NoDeaths.RData`, written as CSV |
| `swabs.csv` | the deposit's `DataFile1.csv`, keeping 10 of its 13 columns (the Ct values are dropped) |

The two `.RData` files were read with `load()` in R into a new environment, checked to hold
only data frames of atomic columns, and written with `write.csv`. Quotes were then removed.
No code from the deposit was run.

**`counts_days0to5.csv`**: `Tanks` (1–20), `Days` (0, 2, 5), `S` and `I` (numbers of newts
qPCR-negative and positive after the decision rules). These are the data of the paper's
Table 2 fits.

**`counts_day8.csv`**: the same on days 5 and 8 for the 10 tanks with no deaths or recoveries
by day 8 (section 2.6).

**`swabs.csv`**: one row per newt and swab. `Sample` (`HDT n`), `Treatment` (`2.5x10^6` for
the newts exposed to zoospores before the experiment, otherwise `Susceptible`), `Tank`,
`InfectionPrev` (initial prevalence), `Day` (2–62, or `Necropsy`), `Pos` and `Neg` (raw
qPCR result), `DDE` (1 if the newt died during the experiment), `DaysSurvival` and
`Quantity` (qPCR load). 2911 rows, 290 newts. `tompros_counts()` applies the authors'
rules (`CodeFile1_PreProcess.R`) and rebuilds both count files exactly.

Read by `load_tompros2022()` and `load_tompros_swabs()`.
