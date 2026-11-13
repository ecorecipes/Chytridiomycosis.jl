# Bsal transmission functions (Malagon 2020, Tompros 2022)

```@meta
CurrentModule = Chytridiomycosis
```

Is transmission of *Batrachochytrium salamandrivorans* (Bsal) density-dependent or
frequency-dependent? The answer decides whether thinning a population can stop an outbreak.
Two experiments on adult eastern newts (*Notophthalmus viridescens*) bear on it:

- **Malagon et al. (2020)**, *Sci Rep* 10:5584
  ([doi:10.1038/s41598-020-62351-x](https://doi.org/10.1038/s41598-020-62351-x)), counted
  contacts between uninfected newts at three densities. They fitted a saturating (Holling
  type II) contact function and put it in an SEI model of a cohabitation experiment.
- **Tompros et al. (2022)**, *Transbound Emerg Dis* 69:731–741
  ([doi:10.1111/tbed.14043](https://doi.org/10.1111/tbed.14043)), followed epidemics in
  twenty 1-m² tanks with 2 to 32 newts. They fitted eight transmission functions and found
  that frequency dependence fits best.

This page gives the transmission functions as reusable types, the two experiments' models,
and refits of both analyses to the authors' data. Equations and values were checked against
the published articles, Malagon's Supplementary Information and Tompros's appendix. The
data are the authors' deposits: Malagon's Supplementary File 2 (`data/malagon2020`) and the
Tompros TRACE deposit
([doi:10.7290/xw8yy2RBeu](https://dx.doi.org/10.7290/xw8yy2RBeu); `data/tompros2022`).

!!! note "Reproduction"
    These are reproductions of published analyses, not new models. What reproduces and what
    does not is recorded in the [validation targets](validation.md) (sources `malagon2020`
    and `tompros2022`).

## Transmission functions

A transmission function gives `φ(S, I, N)`, the total rate of new infections when `S`
susceptible and `I` infectious hosts are in a population of `N`. Tompros et al. write it as
`φ = α S I / N`, where `α` is the per capita contact rate times the probability that a
contact infects. Frequency dependence keeps `α` constant. Density dependence makes it grow
with `N`. The saturating forms fall between the two.

```@docs
TransmissionFunction
DensityDependent
FrequencyDependent
PowerLaw
QuadraticSaturation
StructuralSaturation
PrevalenceDependent
PrevalencePowerLaw
PrevalenceSaturation
HollingContact
infection_rate
force_of_infection
transmission_coefficient
```

[`basic_reproduction_number`](@ref) has a method for a transmission function: the
coefficient divided by the rate at which infectious hosts are removed.

### In a Schmidt-style S–L–I model

Any of these can replace the mass-action term of the fire salamander model of Schmidt et
al. (2017). With `DensityDependent(p.β)` the result is the published model. The newt
estimates are per newt per 1-m² tank per day, so they do not transfer to that model's units
(adult females per hectare per year) without a scaling argument, which is not made here.

```@docs
TransmissionSLIProblem
```

## Tompros et al. (2022): eight candidate functions

In each tank, death and recovery are ignored over the first five days. The number of
susceptibles then falls by one at rate `φ(S, N − S, N)`, and its probability distribution
follows the master equation (appendix eq. 2):

```
dp_S/dt = φ(S+1, N−S−1, N) p_{S+1} − φ(S, N−S, N) p_S
```

The likelihood multiplies, for each tank and each interval between swab days (0, 2 and 5),
the probability of the observed number of susceptibles given the numbers at the start of the
interval.

```@docs
SusceptibleChainProblem
susceptible_distribution
expected_susceptibles
```

### Refit

| Function | NLL (paper) | NLL (refit) | AIC (paper) | AIC (refit) |
|---|---|---|---|---|
| frequency-dependent | 63.2 | 63.154 | 128.3 | 128.31 |
| quadratic saturation | 62.2 | 62.162 | 130.3 | 130.32 |
| power law | 63.2 | 63.215 | 132.4 | 132.43 |
| structural saturation | 62.5 | 62.433 | 132.9 | 132.87 |
| simple prevalence | 68.7 | 68.672 | 139.3 | 139.34 |
| power law (prevalence) | 68.4 | 68.383 | 142.8 | 142.77 |
| quadratic saturation (prevalence) | 68.7 | 68.672 | 143.3 | 143.34 |
| density-dependent | failed to converge | 77.414 | — | 156.83 |

The refit reproduces the ranking and the estimates of Table 2, including `β = 0.230`
day⁻¹ and its interval [0.181, 0.293]. It also reproduces the hold-out "AIC" of section 3.3
(65.9 and 70.6) and the paper's R₀ = 4.9. It also rebuilds the authors' processed counts
from the raw swab results, using their false-negative and false-positive rules.

!!! warning "Differences from the paper"
    - **Density dependence.** The paper says its fit failed to converge. Here it has a
      finite maximum (β = 0.0117 newt⁻¹ day⁻¹) that is 28.5 AIC units worse than frequency
      dependence. That supports the paper's conclusion.
    - **κ₀ of quadratic saturation.** The refit gives 0.0227, against 0.026 printed. κ₀ is
      barely identified (Fig. 4 puts 95% of its profile in [0, 34.2]); κ₁, κ₂ and the NLL
      agree.
    - **Mortality.** The paper reports 45% total mortality. The deposit gives 42.9% of the
      newts exposed by contact and 53.8% of all newts.
    - **The interval for β** is the central 95% of the normalised likelihood, a credible
      interval under a flat prior on the grid, not a likelihood-ratio interval.

```@docs
TomprosTank
tompros2022_path
load_tompros2022
load_tompros_swabs
clean_infection_series
tompros_status
tompros_counts
TomprosModel
transmission_function
tompros_model_selection
likelihood_interval
tompros_negloglik
TomprosParameters
```

[`maximum_likelihood`](@ref) and [`parameter_names`](@ref) have methods for a
[`TomprosModel`](@ref), and [`basic_reproduction_number`](@ref) for
[`TomprosParameters`](@ref).

## Malagon et al. (2020): contacts and cohabitation

Experiment 1 filmed newts at 2, 4 and 8 per m², with and without plants. The contact rate
per newt was fitted to

```
f(N) = c N / (K + N)
```

using the three densities plus a fourth point from the cohabitation containers, which held 2
newts in 0.015 m² (about 133 m⁻²). In Experiment 2 an infected newt was housed with a
susceptible one for 1 second, 10 minutes or 30 minutes. The SEI model of that experiment is

```
dS/dt = −β c S I_c / (K + N)
dE/dt =  β c S I_c / (K + N) − γ E
dI/dt =  γ E − d I
```

where the infected newts `I_c` are present only during the cohabitation. `β` (the probability
that a contact transmits) and `γ` (the rate of becoming qPCR-positive) were chosen on a grid
to match the cumulative proportion of newts positive at days 6, 12, 18 and 24.

The refit reproduces:

- the means behind Fig. 1;
- the four ANOVA F statistics;
- `c` and `K` of Supplementary Table S1 (2153.29 and 79.55 without plants, 587.31 and 52.78
  with plants);
- all twelve `γ` and sums of squares of Table S1, and eleven of its twelve `β`;
- the statements about Fig. 6;
- R₀ = 1.9 for complex habitat.

!!! warning "What the R₀ of Malagon et al. measures"
    - **Fig. 6 is not a population.** The caption describes a population after one infected
      newt arrives. The authors' code instead simulates the cohabitation experiment: 66.5
      infected and 66.5 susceptible newts per m², together for 10 minutes. There is no
      transmission afterwards.
    - **R₀ is not an epidemic R₀.** R₀ = −ln(1 − z)/z is applied to the cumulative
      proportion infected. For plants, the curve's asymptote (0.761) gives 1.9. Without
      plants, the asymptote is 0.988, which gives 4.5. The paper's 0.95, which gives 3.2, is
      the value at about one month.
    - **Neither number is the R₀ of an epidemic in a newt population.** With the fitted
      contact rate and `d` as the removal rate, `β c N/(K + N)/d` is about 5.7 × 10⁴ at
      133 newts m⁻² (derived here, not published).
    - **Units and other slips.** Table S1's `c` is per day, although the text says per hour.
      The code's bounds on `c` and `K` are the reverse of the printed ones. In the table's
      row "Plants In, day 18, 10 min", β is printed as 0.93, but the grid minimum is at 0.98,
      which matches the printed sum of squares.

```@docs
MalagonParameters
contact_rate
CohabitationProblem
cohabitation_curve
final_size_reproduction_number
malagon2020_path
load_malagon2020
contact_means
contact_anova
fit_contact_function
fit_cohabitation
malagon_table_s1
```

## A nonparametric transmission term

The wiki plans to estimate the transmission function nonparametrically, as a term in the
compartmental skeleton (PartiallySpecifiedModels.jl). That is not implemented here.
