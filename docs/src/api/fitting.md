# Fitting to field data

```@meta
CurrentModule = Chytridiomycosis
```

!!! note "New analysis"
    The load equation and these fits are new work in this package, not a published analysis.
    The data are from Knapp et al. (2022), CC BY 4.0.

The frog load equation `dL/dt = ι + r0 L − s(t) L/(h + L)` fitted by maximum likelihood to the
Dusy Basin 2012 skin swabs of Knapp et al. (2022)
([doi:10.7717/peerj.12712](https://doi.org/10.7717/peerj.12712)). Four mechanisms for the extra
clearance `s(t)` are compared. A [`KnappLoadModel`](@ref) implements the LogDensityProblems
interface (`logdensity`, `dimension`, `capabilities`), so it can also be passed to
LogDensityProblems-based samplers. See the tutorial
[Fitting the frog load model to a field trial](../tutorials/03_knapp_fit.md), and
[Bayesian fit to field data](knapp_bayes.md) for priors, NUTS and PSIS-LOO.

## Data

```@docs
KnappSwab
knapp2022_path
load_knapp2022
```

## Model and fitting

```@docs
KnappLoadModel
parameter_names
knapp_start
predicted_loads
maximum_likelihood
profile_likelihood
```
