# Age truncation and metapopulation viability (Heard 2024)

```@meta
CurrentModule = Chytridiomycosis
```

A reimplementation of the models of Heard, Scroggie, Hollanders and Scheele (2024), *J Anim
Ecol* 93:1670–1683
([doi:10.1111/1365-2656.14177](https://doi.org/10.1111/1365-2656.14177)), for the green
and golden bell frog (*Litoria aurea*) and the growling grass frog (*L. raniformis*).

!!! warning "Checked against the authors' code, not the paper"
    The full text could not be retrieved. Everything here follows the authors' R code and
    data (Zenodo [doi:10.5281/zenodo.10457229](https://doi.org/10.5281/zenodo.10457229),
    CC BY 4.0) and the source of the R package `steps` 1.3.0 that they used. The compendium
    records are therefore unverified. Only the abstract was read.

!!! note "Bd's role"
    Bd acts **only through an imposed change in adult (and sub-adult) survival**, from the
    pre-Bd estimate (museum specimens) to the post-Bd estimate (remnant populations). The
    models have no infection, load or transmission process.

## Survival from age structure

Ages (lines of arrested growth) of 403 frogs, bundled in `data/heard2024`. Age − 1 is
geometric with annual survival `S` (a catch curve, Scroggie 2012):

```
P(a) = (1 − S) S^a,   logit S = α + β₁ [L. raniformis] + β₂ [post-Bd] + u_site
```

with Normal site effects and the priors of the authors' NIMBLE code. The posterior gives
0.32 (post-Bd) against 0.65 (pre-Bd) for *L. raniformis* and 0.36 against 0.69 for *L.
aurea*; the authors simulated with 0.319, 0.641, 0.352 and 0.673. Survival roughly halves.

```@docs
HeardAgeRecord
load_heard2024
HeardCatchCurve
heard_loglikelihood
heard_survival
HeardCatchCurveFit
heard_survival_summary
```

[`parameter_names`](@ref) also has a method for a [`HeardCatchCurve`](@ref).

## Spatial metapopulation

A female-only two-stage (sub-adult, adult) matrix model in every habitat cell of a 25 × 25
landscape of 65 m cells, with global environmental stochasticity, demographic
stochasticity, exponential-kernel dispersal and a ceiling carrying capacity, run for 50
years over a grid of recruitment (`mets`), fecundity CV and mortality CV scenarios and 20
landscapes. The reimplementation reproduces the authors' archived persistence fractions
(0.18 and 0.82 for *L. raniformis*, 0.27 and 0.87 for *L. aurea*, post- and pre-Bd).

!!! note "Quirks of `steps` kept by default"
    With `as_coded = true` the simulation keeps two behaviours of `steps` 1.3.0: a cell in
    the last column of the grid disperses from a neighbourhood shifted to column 0 of the
    next row, and after the ceiling the whole landscape is rounded by one multinomial draw,
    which re-allocates individuals among patches and stages every year. `as_coded = false`
    removes both.

```@docs
HeardParameters
heard_scenario_grid
heard_projection_matrix
heard_survival_sd
heard_lambda
HeardMetapopulationProblem
HeardMetapopulationSolution
heard_persistence
heard_persistence_grid
heard2024_path
heard2024_landscapes
```
