# Stochastic individual-based model (Drawert 2017)

```@meta
CurrentModule = Chytridiomycosis
```

This is a reference implementation of the model of Drawert, Griesemer, Petzold and Briggs
(2017, *J R Soc Interface* 14:20170480,
[10.1098/rsif.2017.0480](https://doi.org/10.1098/rsif.2017.0480)). It takes the
[Briggs et al. (2010)](briggs2010.md) load model and makes it individual-based and
demographically stochastic: each frog carries a whole number of sporangia, and the pond holds a
whole number of zoospores. The paper asks which one-time intervention best helps a frog
population survive the first Bd outbreak: treating frogs with antifungals, culling frogs, or
disinfecting the pond.

It was checked against the paper, its supplementary text and the authors' C code and
parameter table (figshare [10.6084/m9.figshare.5328127](https://doi.org/10.6084/m9.figshare.5328127),
CC BY 4.0). Time is in days.

## The model

For each live frog `i` with load `S_i`, and the zoospores `Z` in a pool of volume 1, the
deterministic version is (eqs 3.1–3.2)

```math
\frac{dS_i}{dt} = \nu(S_i)\,\gamma Z + \nu(S_i)\,\eta f S_i - \sigma(S_i) S_i, \qquad
\frac{dZ}{dt} = \sum_i \bigl[\eta(1-f) S_i - \gamma Z\bigr] - \mu Z .
```

A frog dies when its load exceeds `Smax = 10 000`, and its sporangia die with it.

- Zoospores meet each frog at rate `γ` and encyst with probability `ν(S)`.
- Sporangia release `η` zoospores a day. A fraction `f` of them meet their own frog again,
  and the rest go to the pool.
- Sporangia are lost at rate `σ(S)`, and zoospores die at rate `μ`.

There are three variants:

- **Baseline**: `ν = ν0` and `σ = σ0`. These are the Briggs et al. (2010) equations.
- **Nu**: encystment falls with load, `ν(S) = ν0 exp(−ν1 S)`.
- **Sigma**: sporangium loss rises with load, `σ(S) = σ0 + σ1 S`.

In the Nu and Sigma variants, `ν1` or `σ1` is chosen so that loads would level off at
`S* × Smax`, with `S* = 1.2` or `2` (just above the lethal load).

The stochastic version (ESM Table S3) turns each term into an event: a zoospore encysting, a
sporangium reinfecting its own frog, a zoospore being released, a sporangium being lost, a
zoospore dying. [`DrawertIBMProblem`](@ref) simulates it exactly, event by event, as the
authors' C code does. It starts with one frog carrying 100 sporangia, every other frog
uninfected and an empty pool. A run ends in one of three ways:

- frog extinction: all frogs have died;
- fungal extinction: no sporangia or zoospores are left;
- persistence: frogs and fungus are both present after 1000 days.

The paper samples 2132 parameter points over wide ranges (Table 1). It kept only points
where Bd grows at 0.035–0.35 per day, so that every frog dies in the deterministic model,
and ran each point 30 times. It then sorted the points into five groups by the mean and
spread of the fraction of frogs that survive (Fig. 3):

- A: all frogs always die;
- B: almost all die;
- C: highly variable outcomes;
- D: most survive;
- E: partial die-offs.

The authors' table of points and results is bundled ([`load_drawert2017`](@ref)).

```@docs
DrawertParameters
model_variant
encystment_probability
sporangium_loss_rate
drawert_variant
equilibrium_load_ratio
load_growth_rate(::DrawertParameters)
nondimensional_parameters
nondimensional_ranges
nondimensional_distance
DrawertIBMProblem
DrawertODEProblem
```

## Interventions

Each intervention is a single pulse on day 14, 28 or 56, treating 25%, 50% or 75% (§3.3):

- cleaning frogs clears all sporangia from a random set of frogs;
- culling removes a random set of frogs;
- cleaning the environment removes zoospores from the pool.

The number of frogs (or zoospores) treated is binomial. The paper's Fig. 4 counts the
fraction of points at which the mean surviving fraction rises (or falls) by at least 0.1,
0.3 or 0.5. The paper's main finding is that cleaning frogs helps at 20–30% of points and
rarely harms, culling is riskier, and cleaning the environment does little.

The authors' code removes *every* zoospore when it cleans the environment, whatever the
fraction (`clean_env()` sets `Z = 0` before drawing the binomial). This fits the paper's
remark that this treatment is insensitive to the fraction removed. `as_coded = true`
reproduces the code; the default follows the paper's description.

```@docs
DrawertIntervention
drawert_interventions
drawert_summary
outcome_class
survival_group
intervention_effect
```

## Data

```@docs
drawert2017_path
load_drawert2017
load_drawert2017_reproduction
drawert2017_draws
```

## Outcome table for decision analysis

!!! note "Package additions"
    The treatment mortality, several interventions in one run, the table below and the survey
    functions are additions of this package, not published results. With the defaults the
    model is the published one.

`scripts/drawert2017_decision_table.jl` runs the stochastic model over draws × actions ×
replicates, and writes one row per run to `data/drawert2017/decision_runs.csv`. Utilities,
costs and decision rules are left to the user. The seed depends only on the draw and the
replicate, so all actions share random numbers up to the first intervention. `saveat` days
before that give the state on which a survey would be taken.

```@docs
load_drawert_decisions
load_drawert_decision_draws
drawert_run_seed
```

## Adapter for DecisionAnalysis.jl

[DecisionAnalysis.jl](https://github.com/ecorecipes/DecisionAnalysis.jl) is a weak dependency:
this package does not need it, and when both are loaded an extension gives its adapter verbs
methods for [`DrawertIBMProblem`](@ref) and [`DrawertODEProblem`](@ref). Then

```julia
using Chytridiomycosis, DecisionAnalysis
draws = load_drawert_decision_draws()
belief = SampleBelief(draws; weights = [d.weight for d in draws])
prob = DrawertIBMProblem(DrawertParameters(first(draws)); tfinal = 120.0)
actions = (; none = nothing, clean = DrawertIntervention(:clean_frogs, 14, 0.75))
table = outcomes(prob, actions, belief; replicates = 6, seed = 2017)
```

runs the model over draws × actions × replicates with common random numbers and returns an
`OutcomeTable`, to which a utility is applied afterwards. The four verbs are:

| Verb | For the Drawert problems |
|---|---|
| `with_parameters(prob, θ)` | the fields of [`DrawertParameters`](@ref) named in the draw `θ` replace those of `prob`; other columns (`draw_id`, `group`, `weight`) are ignored; `ν1` and `σ1` are not recomputed |
| `apply_action(prob, action)` | `nothing` (no intervention), a [`DrawertIntervention`](@ref), a named tuple `(; kind, day, fraction)`, or a list of those (several pulses; individual-based model only) |
| `simulate(prob; rng)` | `solve(prob; rng)`; the deterministic model ignores `rng` |
| `outcome(prob, run)` | counts: `initial_frogs`, `survivors`, the outcome class as three 0/1 columns, and for the individual-based model `uninfected`, `infected`, `handled`, `treatment_deaths` and the counts on each day of `evaluate` |

Outcomes are counts of frogs. No cost or utility is built in, and no intervention is always
among the actions that the user lists. The stored table is available in the same form:

```@docs
drawert_outcome_table
drawert_survey_study
```

The tutorial [Which response to a Bd outbreak, with new actions?](../tutorials/11_decision_adapter.md)
uses the adapter.

## Survey data from a run

```@docs
simulate_survey(::Chytridiomycosis.Random.AbstractRNG, ::DrawertParameters, ::NamedTuple, ::SurveyObservation, ::SurveyPond)
survey_loglikelihood
```

## What reproduces

The heavy runs are `scripts/drawert2017_*.jl`. Their output is stored in
`data/drawert2017/reproduction.toml`, which the tests and validation targets read.

- **Deterministic model, all 2132 points.** λ equals the authors' column (to 1e-6).
- **Deterministic extinction.** Every point ends in frog extinction, in the Baseline and in the
  four Nu/Sigma variants, "by construction" as the paper says.
- **Stochastic Baseline, 116 points.** We drew a stratified random sample of points and ran
  each 30 times, then compared it with the authors' 30 runs:
  - 113/116 points have the same Fig. 2a outcome class and 105/116 the same Fig. 3 group;
  - no Fisher test on frog-extinction counts rejects, and 5/116 Welch tests on the mean
    surviving fraction reject at 5%;
  - the weighted fractions are 65.4% frog extinction (paper 65%), 24.2% fungal extinction
    (24%), 10.4% variable (~10%) and 65.4% in Group A (65%).
- **Nu and Sigma.** Load regulation lets frogs and fungus coexist (18–32% of runs at
  S* = 1.2 on a random 12-point sample, against none in the Baseline).
- **Interventions** (Fig. 4, Fig. S14), on 66 sample points and two full 27-strategy grids:
  - cleaning frogs never harms, and helps at most Group B/C points;
  - culling helps and harms in Groups B/C, and harms in Group D, more when larger and earlier;
  - cleaning the environment changes outcomes no more often than a placebo (a second
    no-intervention ensemble).
- **Not reproduced.** The overall Fig. 4 bars and Fig. 2a's Nu/Sigma bars need the full design.

Three printed statements disagree with the authors' own data or with the rest of the paper:

- **ESM eq. (11) for ν1** lacks the logarithm and the `ηf` term of the S* formula it should
  invert. With the printed form, frogs would survive deterministically at 270 points.
  `drawert_variant(p, :nu; as_printed = true)` uses it.
- **"Death of all frogs within 14–275 days" (§3.2)** fits the day of the *first* death:
  rounded up to whole days, it spans exactly 14–275. The last death comes as late as day 567.
- **"Coexisting in about 0.01% of the cases"**: the authors' table has 3 of 2132 points
  (0.14%).

ESM eq. (10) also prints `log(c_range)` where the authors' points are spaced by the range of
`log c` (no two closer than 0.0901, with ε = 0.09).
