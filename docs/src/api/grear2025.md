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

The model is checked against the published full text, and the authors' USGS code release
([doi:10.5066/P1MXUJZM](https://doi.org/10.5066/P1MXUJZM), CC0) supplies the data file,
which is bundled in `data/grear2025`. Their formula and priors agree with eqs 1–10 of the
paper. Two differences between paper and code are switches:

- **The σ prior** (`sigma_prior`). The paper puts the half-Student-t(3, 0, 2.5) prior on the
  variance σ² (`:as_printed`). The code leaves brms' default, which is on σ (`:as_coded`, the
  default).
- **The phylogeny** (`tree` in [`grear2025_phylogeny`](@ref)). The paper used 30 trees, one
  per imputed data set. The code's file has 100 trees, and its loop over them is broken.

Missing sex, life stage and genus-only *Taricha* were multiply imputed in the paper (30 data
sets). [`Grear2025Problem`](@ref) draws one such data set by the authors' rule
([`grear2025_impute`](@ref)); `missing_covariates = :drop` drops those swabs instead.

What is checked:

- **Counts.** Swabs (11,183), species (54), sites (593, in the Results) and genus-only records
  (42) are reproduced.
  - The abstract's 609 sites contradicts the Results and is not found anywhere in the data
    or code.
  - The modelled sample (11,172 against 11,175) and the missing sex and life-stage counts
    (191 and 123 against 159 and 91) are discrepant.
- **Temperature.** Both parts are concave in temperature, and prevalence peaks at 14.6 °C,
  against "approximately 15 °C" in the paper.
- **Pending.** The printed 95% credible intervals of the 14 fixed effects, the intensity
  optimum (13 °C) and the PSIS-LOO table need the full posterior
  (`scripts/grear2025_nuts.jl`). Those targets carry the printed values.

!!! warning "The posterior mode does not converge"
    Many sites, and some species, have only negative swabs. Their hurdle random effects and
    the random-effect standard deviations keep growing, so [`posterior_mode`](@ref) does not
    reach a finite mode. After 3000 iterations the site standard deviation of the prevalence
    part exceeds 70 on the logit scale.

    The prevalence fixed effects settle within a few hundred iterations. The intensity ones
    drift: the intensity optimum is 11.9, 14.0 and 15.9 °C after 500, 1000 and 3000
    iterations. For that reason only prevalence quantities are compared numerically, and the
    printed credible intervals are compared only with a full posterior.

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
grear2025_impute
Grear2025Problem
grear2025_fit
grear2025_temperature_optimum
grear2025_printed_intervals
grear2025_coefficient_check
```
