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

!!! note "What was checked"
    - **Mosher et al. 2018.** The equations were checked against the accepted manuscript of
      the article, and its models were refitted to the Dryad data
      ([doi:10.5061/dryad.1c990](https://doi.org/10.5061/dryad.1c990), CC0, bundled in
      `data/mosher2018`).
      - The fits reproduce the AICc of 10 of the 12 models in Table 3, to 0.01, and the
        estimates quoted in the text.
      - The deposit has 83 sites; dropping Site 33, which has no covariates, gives the paper's
        82.
      - The transition and detection matrices are written out only in the first author's
        dissertation (Mosher 2017, Ch. 1). The article describes them in words.
    - **Adams et al. 2010.** Checked against the published article. No detection histories
      are deposited, so the models cannot be refitted. The targets use quantities derived
      from the printed coefficients (Tables 4 and 6), and three of them disagree with numbers
      printed in the text.

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

### Data and refit

[`Mosher2018Model`](@ref) is the likelihood with covariates on the logit scale:

- elevation (linear or quadratic), other amphibians, habitat permanence or region for
  colonisation and extinction;
- day of year (quadratic) for toad detection;
- elevation for Bd detection.

Initial occupancy is constant. When both members of a pair (`γA`/`γAB`, `γB`/`γBA`,
`εA`/`εAB`, `εB`/`εBA`, `pA`/`rA`) have the same covariate structure, they share slopes, as
the caption of Table 3 says. Fit the model with [`maximum_likelihood`](@ref).
[`mosher2018_models`](@ref) lists the 12 models printed in Table 3, and
`scripts/mosher2018_refit.jl` fits them all in about two minutes.

The refit reproduces:

- the AICc of 10 of the 12 printed models, and the deviance (2175.3) and `K` (25) of the
  best;
- `εA = εB = 0.00` in the all-constant model; in the best model, `εA = 0.00`,
  `εB = 0.02` and `γA = 1.00`;
- the elevation patterns of Fig. 5: `εAB` and `εBA` peak at about 0.3 at intermediate
  elevations, toads recolonise high sites (`γAB` increases with elevation), and Bd colonises
  low and mid-elevation sites.

!!! warning "Not reproduced"
    - **Model 3** (`γB`, `γBA`, `εBA` quadratic). Our deviance is 0.22 lower than the
      printed one, so the printed fit seems to fall short of the maximum.
    - **The null model.** Table 3 prints `K = 14` and deviance 2232.15. Ours has one
      constant per rate pair: `K = 13` and a lower deviance, 2230.93. The paper does not
      describe its null model further.

!!! note "Choices where the paper is silent"
    - **AICc sample size.** It is the number of site-years with at least one survey (494 at
      the 82 sites). The paper does not state it. This value reproduces every AICc–deviance
      pair in Table 3, which number of sites (82) does not.
    - **Standardised covariates.** This does not change the maximised likelihood.
    - **Site 33 dropped.** This gives the paper's 82 sites; see `data/mosher2018/README.md`.

```@docs
Mosher2018Data
mosher2018_path
load_mosher2018
Mosher2018Model
mosher2018_loglikelihood
mosher2018_site_parameters
mosher2018_parameters_at
mosher2018_models
mosher2018_model
mosher2018_model_selection
```

[`maximum_likelihood`](@ref) and [`parameter_names`](@ref) also have methods for a
[`Mosher2018Model`](@ref).

## Adams et al. 2010: Bd in amphibian populations

Here a "site" is a host population and a "survey" is one animal swabbed once. The
probability `p` that an animal from a population with Bd tests positive combines prevalence
and test sensitivity. Native anurans (72 populations analysed) and bullfrogs (40) were
fitted separately in Program MARK, with logit-linear covariates:

- **Native anurans, best model.** `p` depends on larvae, date, elevation, date × elevation
  and the human footprint index; `ψ` depends on *Rana luteiventris*.
- **Bullfrogs, best model.** `p` depends on larvae, date and stagnant water; `ψ` depends on
  latitude.

[`SingleSeasonOccupancyModel`](@ref) is this likelihood with any design matrices. No
detection histories are deposited, so the paper's models cannot be refitted. Its printed
coefficients (Tables 4 and 6) are in the compendium.

From those coefficients:

- **Bullfrogs.** `ψ = 0.664` (CI 0.429–0.838) is reproduced, and so are all four odds
  ratios quoted in the text.
- **Native anurans.** The *R. luteiventris* odds ratio (4.2) is reproduced. Three printed
  numbers do not follow:

  | Quantity | Printed | From Table 4 |
  |---|---|---|
  | `ψ` | 53.4% | 0.741 (non-*R. luteiventris* populations) |
  | larval odds ratio | 0.10 | 0.088 |
  | odds ratio per human-footprint point | 2.37 | 1.54; 2.37 is close to 2.39, the ratio per standard deviation |

```@docs
SingleSeasonOccupancyModel
adams2010_coefficients
adams2010_odds_ratio
adams2010_occupancy
```

Appendix A: because a population has a finite number `N` of animals, the probability that
Bd occurs and could be detected is `ψ[1 − (1 − p′)^N]`, not `ψ`. The difference is small
for `p′ > 0.05` and `N > 50`.

```@docs
finite_population_occurrence
```
