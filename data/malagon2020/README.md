# Malagon et al. 2020: newt contacts and Bsal cohabitation experiments

- **Paper**: Malagon DA, Melara LA, Prosper OF, Lenhart S, Carter ED, Fordyce JA, Peterson AC,
  Miller DL, Gray MJ (2020). Host density and habitat structure influence host contact rates
  and *Batrachochytrium salamandrivorans* transmission. *Sci Rep* 10:5584.
  https://doi.org/10.1038/s41598-020-62351-x
- **Source**: Supplementary File 2 of the article (a zip of `Repository_Files/`: data CSVs, an R
  script and MATLAB code), from the publisher's site. The zip and the extracted files are
  kept in `lab/papers/malagon2020_si2.zip` and `lab/papers/malagon2020_si2/`. Retrieved
  2026-11-12.
- **Licence**: the article and its supplementary information are CC BY 4.0. The deposit's
  README (`MALAGON_DATASET_ReadmeTemplate_2020-02-13.txt`) lists no licence or restrictions
  ("None"). Attribution: Malagon et al. 2020.

## Files

| file | origin |
|---|---|
| `contactrates.csv` | the deposit's `contactrates.csv`, unmodified (it starts with a UTF-8 byte-order mark) |
| `survival.csv` | the deposit's `DansSurvival.csv`, unmodified |
| `cohabitation_swabs.csv` | transcribed from the matrix `YMAT` in the deposit's `MasterRun.m` |
| `contact_fit_inputs.csv` | transcribed from `total_contacts` in the deposit's `MasterRun.m` (plants in) and `objfun.m` (plants out) |

**`contactrates.csv`** (Experiment 1; 20–22 °C; 1-m² pools): `Density` (Two, Four or Eight
newts), `Plants` (No, Yes: nine artificial plants), `Contacts/Hr` (contacts counted in one
hour of video), `Day`, `Diel` (1–4: 0–6, 6–12, 12–18 and 18–24 h) and `Pool` (1–6). 672 rows.

**`survival.csv`** (Experiment 2; 14 °C): one row per newt. `Treatment` is the disease state
of the infected partner (Day 12, 18 or 24) or `5x10^6` for the newts exposed to zoospores.
`Days Survival ` is censored at 90 days, `Dead ` is 1 for a death, `Group` is Susceptible
or Infected, and `CT` is the contact treatment (`1second`, `10minutes`, `30minutes`). Some
column names end in a space.

**`cohabitation_swabs.csv`**: the cumulative proportion of the five susceptible newts per
treatment that were qPCR-positive on two swabs, on days 6, 12, 18 and 24 after cohabitation.
The authors' grid search (`SolveSEImodel6.m`) fits these. The mapping of the six rows of
`YMAT` to disease state and cohabitation time (12/10, 12/30, 18/10, 18/30, 24/10, 24/30)
follows the row order of Supplementary Table S1. The refit reproduces every γ and SSR of the
table with this mapping.

**`contact_fit_inputs.csv`**: total contacts per hour at 2, 4 and 8 newts m⁻², and at
133 m⁻² (the Experiment 2 containers, 2 newts in 0.015 m²), as fitted by the authors' code.
The first three values are the back-transformed means `exp(mean(log(x + 1))) − 1` of the
pool-day means of `contactrates.csv`. `contact_means()` recomputes them; the plants-out
values are rounded to one decimal in the code. The 133 value is 7466 without plants
(56.1 contacts newt⁻¹ h⁻¹; the SI says 112 per newt per hour). With plants it is 2331.384,
which is close to 7466 times the plants/no-plants ratio at 8 newts; the code does not
explain it.

Read by `load_malagon2020()`.
