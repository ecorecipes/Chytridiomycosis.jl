# Hurdle model of Bd across species (Grear 2025)

```@meta
CurrentModule = Chytridiomycosis
```

Grear et al. (2025, *Ecological Monographs* 95:e70001,
[doi:10.1002/ecm.70001](https://doi.org/10.1002/ecm.70001)) swabbed 11,183 amphibians of 54
species across the United States. They fitted a Bayesian **hurdle-lognormal** model with brms
to the *Bd* quantity per swab. It has two parts:

- a logistic part for whether a swab is positive (prevalence);
- a lognormal part for the load of positive swabs (intensity).

Both parts use the same covariates:

- age and sex;
- the mean daily maximum temperature over the 30 days before sampling, and its square;
- precipitation over the previous 10 days;
- random intercepts for site, NOAA region, and species (phylogenetically correlated).

!!! warning "Source"
    The full text could not be fetched by an automated client, so only the abstract was read.
    The model, priors and data here come from the authors' USGS code release
    ([doi:10.5066/P1MXUJZM](https://doi.org/10.5066/P1MXUJZM), CC0). Its data file is bundled
    in `data/grear2025`. No printed coefficient is reproduced; those validation targets are
    pending.

What is checked:

- The abstract's counts of swabs (11,183) and species (54) are reproduced. The site count is
  not: the data have 594 sites, the abstract says 609.
- The abstract says infection peaks at intermediate temperatures. An approximate posterior
  mode of a reduced model agrees: temperature effects are concave, peaking at 14.9 °C for
  prevalence and 10.9 °C for intensity.

## Hurdle regression layer

[`HurdleRegressionProblem`](@ref) adds covariates and random intercepts to the two-part model
of [`LoadHurdleModel`](@ref), in brms' `hurdle_lognormal` parametrisation. There, `hu` is the
probability of a *zero*, so a prevalence coefficient is minus the corresponding `hu`
coefficient.

The problem is a LogDensityProblems target with an analytic gradient, so any gradient-based
sampler can draw from the full posterior. [`posterior_mode`](@ref) finds its mode with L-BFGS.

```@docs
HurdleRegressionProblem
hurdle_coefficients
posterior_mode(::Any)
phylogenetic_correlation
```

## Grear et al. 2025

```@docs
grear2025_path
Grear2025Data
load_grear2025
grear2025_summary
grear2025_phylogeny_path
grear2025_phylogeny
Grear2025Problem
grear2025_fit
grear2025_temperature_optimum
```
