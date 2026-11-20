# Watering timing and extinction risk (Turner 2025)

```@meta
CurrentModule = Chytridiomycosis
```

A reimplementation of the stochastic population viability model of Turner et al. (2025),
*Conserv Sci Pract* 7:e70022
([doi:10.1111/csp2.70022](https://doi.org/10.1111/csp2.70022)), for the southern bell frog
*Litoria raniformis* in a regulated floodplain wetland that receives environmental water in
August (~10 °C), September (~15 °C) or October (~20 °C).

!!! note "Sources"
    The model follows the authors' R code (Zenodo
    [doi:10.5281/zenodo.10921748](https://doi.org/10.5281/zenodo.10921748), CC BY 4.0) and
    was checked against the published paper, which describes it in outline.

!!! note "Bd's role"
    The authors converted watering month into adult and metamorph survival through Bd
    load–survival relationships and stored the result as precomputed arrays. The code that
    produced them is not in the archive, and the paper describes the step without printing
    its regression coefficients or temperature series, so here **Bd enters only through
    those survival series**, which are bundled (`data/turner2025`).

## Model

Six female age classes (eggs, then ages 1–5), annual step, 52 years. Each year the
survivals, hatching, metamorphosis and fecundities are redrawn; density dependence acts on
first-year, juvenile and adult survival and on egg laying; in a dry (unwatered) year there is
no breeding. A replicate is extinct when no female aged 1 or more is left.

The reimplementation reproduces the authors' archived outputs (e.g. extinction 0.16, 0.06
and 0.04 for August, September and October with independent dry years, and the median
minimum population sizes of all nine scenarios). The paper's values (0.13 and 0.05 for
August and September; 0.025 and 0.0004 with at most two dry years; the quasi-extinction
probabilities of its Fig. 3) are not reproduced, and the archived outputs do not reproduce
them either: the paper's figures match a table hard-coded in the authors' script and
appear to come from an earlier version of the code or arrays. October (0.04) agrees.

!!! note "As coded"
    The authors' `pop.func` reads the survival arrays with the replicate count (10,000) as
    the row index, so every replicate uses the last replicate's survival series. That is the
    default here; pass the full array ([`load_turner2025_survival_array`](@ref)) to use one
    row per replicate. Wetness sequences, not archived as code, are drawn by
    [`turner_watering_sequence`](@ref).

```@docs
TurnerParameters
turner_egg_multiplier
turner_juvenile_multiplier
turner_tadpole_multiplier
turner_egg_ceiling
turner_watering_sequence
turner2025_path
load_turner2025
load_turner2025_survival_array
turner2025_archive_path
TurnerPVAProblem
TurnerPVASolution
turner_extinction_probability
turner_minimum_population
turner_quasi_extinction_probability
```
