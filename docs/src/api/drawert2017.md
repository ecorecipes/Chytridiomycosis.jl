# Stochastic individual-based model (Drawert 2017)

```@meta
CurrentModule = Chytridiomycosis
```

This is a reference implementation of the model of Drawert, Griesemer, Petzold and Briggs
(2017, *J R Soc Interface* 14:20170480,
[10.1098/rsif.2017.0480](https://doi.org/10.1098/rsif.2017.0480)). It takes the
[Briggs et al. (2010)](published.md) load model and makes it individual-based and
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
```

## What reproduces

The heavy runs are in `scripts/drawert2017_*.jl`, and their results are stored in
`data/drawert2017/reproduction.toml`. What was found:

- **Growth rate.** The linearised growth rate λ matches the authors' column at all 2132
  points, to the six figures they printed.
- **Deterministic model.** Every point ends in frog extinction, in the Baseline and in all
  four Nu and Sigma variants, as the paper says it should "by construction".
- **Stochastic Baseline against the authors' runs.** On a stratified sample of 116 points
  with 30 runs each, our runs match the authors' 30 runs point by point: the same outcome
  class and group at most points, and two-sample tests reject at about their nominal rate.
  The weighted fractions of points ending in frog extinction, fungal extinction or variable
  outcomes are close to the paper's 65%, 24% and 10%.

Three printed statements disagree with the authors' own data or with the rest of the paper:

- **ESM eq. (11) for ν1** lacks the logarithm and the `ηf` term of the S* formula it should
  invert. `drawert_variant(p, :nu; as_printed = true)` uses the printed form.
- **"Death of all frogs within 14–275 days" (§3.2)** fits the day of the *first* death:
  rounded up to whole days, it spans exactly 14–275. The last death comes as late as day 567.
- **"Coexisting in about 0.01% of the cases"**: the authors' table has 3 of 2132 points
  (0.14%).

ESM eq. (10) also prints `log(c_range)` where the authors' points are spaced by the range of
`log c` (no two closer than 0.0901, with ε = 0.09).
