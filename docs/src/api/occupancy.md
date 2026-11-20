# Occupancy models (Mosher 2018, Adams 2010)

```@meta
CurrentModule = Chytridiomycosis
```

Occupancy models estimate where a species, or a pathogen, occurs when surveys can miss it.
This page has a reusable occupancy layer and two published Bd applications:

- **Mosher et al. (2018)**, *Ecol Appl* 28:926–937
  ([doi:10.1002/eap.1699](https://doi.org/10.1002/eap.1699)): a dynamic two-species
  occupancy model of boreal toads and Bd at historic breeding sites in the southern Rocky
  Mountains, 2001–2010, with colonisation and extinction of each species conditional on the
  other and on elevation.
- **Adams et al. (2010)**, *Ecol Appl* 20:289–302
  ([doi:10.1890/08-2319.1](https://doi.org/10.1890/08-2319.1)): single-season occupancy of
  Bd in native-amphibian and bullfrog populations, where each animal tested is a survey.

!!! warning "What could and could not be checked"
    Neither article could be read (both are paywalled).
    - **Mosher et al. 2018.** The equations come from the first author's PhD dissertation
      (Mosher 2017, Colorado State University). Its Chapter 2 is the pre-publication version
      of the paper, and its Chapter 1 gives the transition and detection matrices. The
      dissertation analysed 82 sites; the Dryad deposit
      ([doi:10.5061/dryad.1c990](https://doi.org/10.5061/dryad.1c990), CC0) has 83. The
      deposit could not be downloaded, so the paper's estimates are **not yet reproduced**:
      its targets are pending, and there is no reader for its CSV files yet.
    - **Adams et al. 2010.** Only the abstract and the two Ecological Archives appendices
      (figshare, CC0) were read. No detection histories are deposited. The covariate models
      of the paper are not implemented; its two occupancy estimates are pending targets.

## Occupancy layer

Single-season occupancy (MacKenzie et al. 2002): a site is occupied with probability `ψ`,
and each survey detects the species with probability `p` if it is there. Dynamic multistate
occupancy is a hidden Markov model over seasons, with the surveys of a season as repeated
observations of a fixed state (robust design). It uses [`hmm_loglikelihood`](@ref) and
[`robust_design_likelihoods`](@ref).

```@docs
occupancy_loglikelihood
occupancy_mle
dynamic_occupancy_loglikelihood
simulate_dynamic_occupancy
```

### Two species: host and pathogen

The conditional parameterisation of Richmond et al. (2010) has four states: `A` (species A
only), `B` (B only), `U` (neither) and `AB` (both). In Mosher et al. species A is the toad
and B is Bd. With rows the state in year `t` and columns the state in year `t + 1`:

```
      A                B                U                AB
A     (1−εA)(1−γBA)    εA γBA           εA(1−γBA)        (1−εA)γBA
B     γAB εB           (1−γAB)(1−εB)    (1−γAB)εB        γAB(1−εB)
U     γA(1−γB)         (1−γA)γB         (1−γA)(1−γB)     γA γB
AB    (1−εAB)εBA       εAB(1−εBA)       εAB εBA          (1−εAB)(1−εBA)
```

In one survey, a site in state `AB` is recorded as `AB` with probability `rA rBA`, as `A`
with `rA(1 − rBA)`, as `B` with `(1 − rA) rBa` and as `U` otherwise. Bd is found only on
swabs of toads that were found, so `pB = rBa = 0`. The Bd-only state is then never
observed, and some parameters are poorly identified (Mosher 2017, Ch. 1).

```@docs
two_species_initial
two_species_transition
two_species_detection
two_species_code
two_species_records
```

## Mosher et al. 2018: boreal toads and Bd

[`Mosher2018Problem`](@ref) projects the expected proportions of sites in each state. With
the dissertation's "SRM boreal toad–Bd" scenario (compendium base set), toad occupancy falls
from 0.90 to 0.71 over four years; with its rare, susceptible host it falls from 0.20 to
0.07. Both match dissertation Table 1.3, which checks the transition matrix.

```@docs
Mosher2018Parameters
initial_state
transition_matrix
detection_matrix
Mosher2018Problem
simulate_mosher2018
```

### Fitting

[`Mosher2018Model`](@ref) is the likelihood with covariates on the logit scale: elevation
(linear or quadratic), other amphibians, habitat permanence and region for colonisation and
extinction; day of year for toad detection; elevation for Bd detection. Initial occupancy is
constant, as in the source. Fit it with [`maximum_likelihood`](@ref).

!!! note "Choices where the source is silent"
    - When both members of a pair (`γA`/`γAB`, `γB`/`γBA`, `εA`/`εAB`, `εB`/`εBA`,
      `pA`/`rA`) have the same covariate structure, they share slopes and differ only in
      intercept. With this rule, the parameter counts of all eleven models listed in
      dissertation Table 2.3 are reproduced. The exception is the null model: Table 2.3
      prints `K = 14`, and this rule gives 13.
    - Covariates are standardised. This does not change the maximised likelihood.
    - AICc uses the number of sites as the sample size.
    - A survey with no recorded day of year gets the mean day of year.

```@docs
Mosher2018Data
mosher2018_path
Mosher2018Model
mosher2018_loglikelihood
mosher2018_site_parameters
```

[`maximum_likelihood`](@ref) and [`parameter_names`](@ref) also have methods for a
[`Mosher2018Model`](@ref).

## Adams et al. 2010: Bd in amphibian populations

Here a "site" is a host population and a "survey" is one animal tested. The likelihood is
[`occupancy_loglikelihood`](@ref). The abstract reports Bd in an estimated 53.4% of 78
native-amphibian populations and 66.4% of 40 bullfrog populations (compendium `adams2010`).
These are pending targets until the detection histories are available.

Appendix A notes that the number of animals that could be sampled is finite, so the
probability that Bd occurs and could be detected is `ψ[1 − (1 − p′)^N]`, not `ψ`.

```@docs
finite_population_occurrence
```
