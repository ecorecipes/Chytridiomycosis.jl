# Bd growth in varying temperatures (Gajewski 2021)

```@meta
CurrentModule = Chytridiomycosis
```

This is a reference implementation of the model of Gajewski et al. (2021), *Ecol Evol*
11:17920–17931 ([doi:10.1002/ece3.8379](https://doi.org/10.1002/ece3.8379)): the optical
density (OD) of a Bd culture follows a logistic curve with a delay, whose growth rate is a
thermal performance curve (TPC). The authors fitted five TPCs to Stevenson et al.'s (2013)
constant-temperature data, and used each fit to predict growth under three fluctuating
regimes. Equations and values were checked against the full text, Appendix S1, and the
authors' code (Zenodo [doi:10.5281/zenodo.4525954](https://doi.org/10.5281/zenodo.4525954)).

## Model

Time in days, temperature `T` in °C:

```
Y(t) = Y0                                          t < d     (Eq. 1)
     = K Y0 / (Y0 + (K − Y0) exp(−r(T)(t − d)))    t ≥ d
dx/dt = r(T(t)) x (1 − x/K) 1(t > d),  x(0) = Y0             (Eq. 3)
```

The five curves are [`Briere`](@ref) (Brière-2), [`Ratkowsky`](@ref), [`SharpeSchoolfieldIkemoto`](@ref),
[`Logan10`](@ref) and [`Stinner`](@ref), in the forms of the authors' code.

!!! warning "Differences between the paper and the code"
    - **Likelihood.** The paper prints a lognormal likelihood for every fit. The code uses it
      for Brière-2 and Ratkowsky only, and a normal truncated at zero for Ikemoto, Logan-10 and
      Stinner. [`GajewskiModel`](@ref) follows the code by default (`likelihood = :as_coded`);
      with the data, −2 log L at the posterior mode plus 2k matches the penalized deviances of
      Table 1 to within 8 units for all five fits. The ΔDIC of Table 1, which favours Stinner,
      therefore compares deviances of different likelihoods.
    - **Ikemoto delay.** The Ikemoto JAGS model uses `exp(d)` as the delay; Table 2 prints `d`
      (0.007), and the predictions use `d`. See `ikemoto_delay` in [`GajewskiParameters`](@ref).
    - **Kelvin.** The fits use `T + 273.15`, the predictions `T + 273.16`, and the last segment
      of the Wet Low Frog regime omits the offset, so the kelvin-based curves see 25 K there
      (`as_coded` in [`gajewski_prediction`](@ref)).
    - **Logan-10 α.** Appendix S1 prints its median as 0.134, outside its own HDI
      (0.885, 0.957); with the data the posterior mode is 0.914. Targets built on the printed
      α are discrepant.
    - **Table B.1.** The logistic fits without a TPC to each regime are not reproduced from the
      deposited data: the mode of the authors' model is far from the printed medians.

The summaries of the curve at the printed posterior medians ([`gajewski_summary`](@ref))
reproduce Table 1 (adjusted limits, maximum r), Table 2 (optima) and Table B.2 (median rate in
each regime) to 0.2 °C and 1 %; they are not the posterior medians of the summaries that the
paper reports, so the agreement is approximate except for Stinner, whose optimum is a
parameter.

```@docs
Stinner
Logan10
SharpeSchoolfieldIkemoto
delayed_logistic
gajewski_curve
GajewskiParameters
optical_density
GajewskiRegime
GajewskiGrowthProblem
gajewski_prediction
regime_growth_rate
adjusted_limits
gajewski_summary
```

## Data and posterior

!!! note "Data not bundled"
    The Zenodo release has no licence, so its data are not bundled.
    [`gajewski2021_path`](@ref) finds them in the lab workspace or in
    `ENV["CHYTRIDIOMYCOSIS_GAJEWSKI2021"]`. The constant-temperature data are Stevenson et
    al.'s (2013) experiment, not an independent one. Validation targets that need the data are
    recorded as pending; `test/gajewski2021.jl` runs them when the data are present.

`scripts/gajewski2021_mcmc.jl` samples the posterior of each fit by random-walk Metropolis
from its mode.

```@docs
gajewski2021_path
load_gajewski2021
GajewskiModel
GajewskiLogisticModel
gajewski_start
posterior_mode(::Union{GajewskiModel,GajewskiLogisticModel})
```

[`parameter_names`](@ref) also has methods for [`GajewskiModel`](@ref) and
[`GajewskiLogisticModel`](@ref).
