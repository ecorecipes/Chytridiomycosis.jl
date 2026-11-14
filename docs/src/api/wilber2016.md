# Integral projection models (Wilber 2016)

```@meta
CurrentModule = Chytridiomycosis
```

An integral projection model (IPM) follows the distribution of a continuous trait through
discrete time. In a host–parasite IPM the trait is pathogen load, so the model sits between an
individual-based load model such as [Briggs et al. (2010)](published.md) or
[Drawert et al. (2017)](drawert2017.md), which follows every host, and a compartmental model,
which follows only how many hosts are infected. This page has two parts: generic IPM
machinery, and the reference implementation of Wilber, Langwig, Kilpatrick, McCallum and
Briggs (2016, *Methods Ecol Evol* 7:1182–1194,
[10.1111/2041-210X.12561](https://doi.org/10.1111/2041-210X.12561)), checked against the
paper, the authors' R code and results
([github.com/mqwilber/ipm_for_parasites](https://github.com/mqwilber/ipm_for_parasites)) and
their data (Dryad [10.5061/dryad.06sd9](https://doi.org/10.5061/dryad.06sd9), CC0, bundled).

## Generic machinery

The state is one uninfected class `S` and the infected density `I(x)` over load `x`
(eqns 1–2 without reproduction):

```math
S(t+1) = s_0(1-\phi) S(t) + \int s(x)\,l(x)\,I(x,t)\,dx, \qquad
I(x',t+1) = \int s(x)\,(1-l(x))\,G(x',x)\,I(x,t)\,dx + s_0\,\phi\,S(t)\,G_0(x').
```

An [`IPMKernel`](@ref) holds the vital-rate functions: survival `s(x)`, growth `G(x′, x)`
(a distribution for each `x`), loss of infection `l(x)`, initial load `G0(x′)` and
uninfected survival `s0`. [`discretize`](@ref) evaluates it on an [`IPMMesh`](@ref) by the
midpoint rule (or exact bin integrals, `integration = :cdf`). Hosts whose next load falls
outside the mesh are lost (`eviction = :none`, as in the paper) or kept in the edge bins
(`eviction = :boundary`). [`eviction`](@ref) measures the loss and
[`mesh_convergence`](@ref) checks the mesh size. The infection probability `φ` is passed to
[`projection_matrix`](@ref) separately, because it may depend on the state.

The matrix analyses work on any projection matrix: [`population_growth_rate`](@ref)
(the dominant eigenvalue `λ`), [`stable_distribution`](@ref),
[`sensitivity_matrix`](@ref) and [`elasticity_matrix`](@ref), and
[`absorption_times`](@ref) (mean and variance of the time to death, from the fundamental
matrix). [`basic_reproduction_number`](@ref) has a method for a [`DiscretizedIPM`](@ref):
`R0` by the next-generation method of Klepac and Caswell (2011).

```@docs
IPMMesh
meshpoints
meshwidth
meshedges
IPMKernel
DiscretizedIPM
discretize
projection_matrix
population_growth_rate(::AbstractMatrix)
stable_distribution
sensitivity_matrix
elasticity_matrix
absorption_times
eviction
mesh_convergence
load_summary
```

The machinery is written here rather than taken from IntegralProjectionModels.jl
(ecorecipes), which sources a dependency by relative path (so CI cannot install it), takes
Plots as a hard dependency and requires SciMLBase 2.

## The Bd–*Rana muscosa* model

Twenty adult frogs were exposed to about 10⁶ Bd zoospores, housed singly at 4, 12, 20 and
26 °C (five per temperature; 26 °C was not analysed) and swabbed every 3 days. One time step
is 3 days. Load `x` is the natural log of zoospore equivalents (ZE) per swab; uninfected is a
separate class, not `x = 0`. The vital rates (Table 1) are regressions on `x` and
temperature `T`:

- survival: `logit s(x) = b00 + b10 x`, fitted at 20 °C (no frog died at 4 or 12 °C);
- growth: normal with mean `b01 + b11 x + b21 T` and variance `v01 exp(2 c01 x)`;
- loss of infection: `logit l(x, T) = b02 + b12 x + b22 T`;
- initial load: normal with mean `b03 + b13 T` and variance `v03 exp(2 c03 T)`;
- infection: `logit φ(T) = b04 + b14 T`; frogs were housed singly, so this is infection from
  the environment.

The kernel is discretised by the midpoint rule with 100 bins on `[−5, 18]`.

```@docs
WilberParameters
load_growth_rate(::WilberParameters)
BdIPMProblem
BdIPMSolution
to_discrete_problem(::BdIPMProblem)
LoadTransmission
```

[`infection_probability`](@ref) has methods for `WilberParameters` (`φ(T)`) and for a
[`LoadTransmission`](@ref) and a state.

### Analyses

```@docs
population_growth_rate(::WilberParameters)
stable_load_distribution
mean_time_to_death
vital_rate_elasticities
```

[`basic_reproduction_number`](@ref) for `WilberParameters` gives the authors' Appendix S1
`R0`, with transmission between frogs at `β = 9.82 × 10⁻⁴` per frog per step among 100
susceptible frogs, either from the full IPM or from its collapse to an SIS model.

### Data and refits

```@docs
wilber2016_path
load_wilber2016
fit_wilber2016
WilberFit
WilberRegression
growth_rate_bootstrap
```

## What reproduces

Paper and deposit are compared in `test/wilber2016.jl` and the validation targets
`wilber2016_*`. "Deposit" means the authors' own result files, decoded from the vector
graphics of their `results/*.pdf` (to about 10⁻⁵) or read from their CSVs.

- **Vital-rate fits.** Every coefficient and covariance of the deposited fit is reproduced
  by [`fit_wilber2016`](@ref) to at least 8 significant figures, as are the growth, loss,
  initial-load and infection tests of the Results (χ² = 196.36, 13.56, 9.8, 14.555, 23.701,
  6.00, 6.0361; t₄₁ = 2.53).
- **Survival.** Table 1's survival coefficients (11.824, −0.8605) and the survival test
  (χ² = 12.197) are those of a fit to **all three temperatures**, not the 20 °C fit the text
  and Fig. 2a describe. The 20 °C fit gives 11.597, −0.8873 and χ² = 8.58, and is what the
  authors' analyses used. `WilberParameters(:wilber2016)` is Table 1 as printed;
  `set = "deposited"` is what the figures were made with.
- **Density-independent IPM** (deposited coefficients, 100 bins): λ for 4–20 °C (0.99896 to
  0.93690), the mean and variance-to-mean ratio of the stable load distribution (2.61 at
  4 °C to 0.85 at 20 °C, Fig. 5 inset), the peaks of Fig. 5 and the mean times to death of
  Fig. 4a all agree with the deposit to plotting precision.
- **R0** (Appendix S1 code, 50 bins): full and collapsed `R0` agree with the deposited plot
  to 4 × 10⁻⁵; `R0 = 1` near 12.2 °C and peaks at 1.68 near 17 °C.
- **Density-dependent transmission** (eqn 13 as coded, 40 × 40 grid, 12–20 °C): the
  population after 120 days and the fraction uninfected agree with the deposited CSVs to
  10⁻¹¹.
- **Bootstrap** (Fig. 4b): medians agree within sampling error; the very wide interquartile
  ranges come from a slip in the authors' code, which never resamples the survival slope
  (`as_coded = true`). With both survival coefficients resampled the ranges are much narrower.

Printed statements that do not hold, or hold only for the code:

- Eqn 13 integrates `x I(x, t)` (log loads); the code uses `ln(Σ eˣ n + 1)`, the log of the
  total zoospore load. The conclusion that density-dependent transmission "had little effect"
  holds as coded; with eqn 13 as printed, `β` matters a lot when environmental infection is
  rare (at 20 °C and `ω = 0.01`, the summer decline rises from 16% to 79% across the `β`
  range, against 16% to 32% as coded).
- "Populations at 12 °C will experience a maximum of a 20% population decline … with 70%
  prevalence": the deposited simulation gives at most 25% and 77%.
- The text describes 20 frogs at three temperatures, five per temperature; the data have five
  frogs at each of four temperatures, and 26 °C is silently excluded.
- The paper scans `ω` over 0.01–0.6; the deposited simulation stops at 0.3.

## How it relates to the Briggs and Drawert models

All three describe Bd load on *R. muscosa*/*R. sierrae*. Briggs et al. (2010) and Drawert et
al. (2017) grow sporangia exponentially (rate `ηνf − σ` per day,
[`load_growth_rate`](@ref)) until a hard lethal load `Smax = 10⁴` sporangia; transmission goes
through a zoospore pool. The IPM is fitted to qPCR loads (ZE per swab; the two load units
have no published conversion) and replaces both mechanisms with regressions:

- **Growth.** The fitted mean growth (`b11 ≈ 0.8 < 1`) is not exponential: loads are pulled
  towards `x* = (b01 + b21 T)/(1 − b11)`, about 9.2 ln ZE (10⁴ ZE) at 20 °C and 1.9 ln ZE at
  4 °C. At 20 °C, ln load rises by 0.62 per day from `x = 0` (1 ZE), above the 0.035–0.35
  per day of sporangia in the Drawert design, and the rate falls to zero at `x*`.
- **Death.** Instead of a threshold, survival declines smoothly with load (the deposited
  fit): a frog at 10⁴ ZE survives a 3-day step with probability 0.97, and survival is ½ at
  13.1 ln ZE (4.7 × 10⁵ ZE), whatever the temperature.
- **Loss of infection** has no counterpart in the Briggs model, and transmission is an
  environmental probability rather than a pool.

`test/wilber2016.jl` checks that with Briggs-style vital rates (exponential growth at the
linearised rate of a Drawert Group A point, death above `Smax`) the IPM machinery reproduces
exponential load growth and the time to reach `Smax`.
