# Load distributions and extinction risk (Louca 2014)

```@meta
CurrentModule = Chytridiomycosis
```

Louca, Lampo and Doebeli (2014, *Proc R Soc B* 281:20132783,
[doi:10.1098/rspb.2013.2783](https://doi.org/10.1098/rspb.2013.2783)) built a generic,
continuous-time model of Bd in a frog population with tadpoles and adults. It tracks the whole
distribution of the degree of infection (DOI, zoospores on a host) in each stage, and a pool
of free zoospores. They derived a basic reproduction quotient `Q_o`, sampled a 25-dimensional
parameter space and ranked the parameters by their "risk effect" on host extinction.

The implementation was checked against the paper and its electronic supplementary material
(ESM), which the first author's site appends to the article (the journal's ESM link returned
HTTP 403). No code or data were deposited. Time is in days and a year is `Y = 365` days.

## The model

``P_s(t, n)`` is the number of hosts of stage `s` (A adults, T tadpoles) with DOI `n` at time
`t`, and `Z(t)` the free zoospores. A zoospore that meets a host raises its DOI by one; the
within-host dynamics move the DOI up or down one at a time (main eqs 2.1–2.2, ESM eqs 1–3):

```math
\begin{aligned}
\frac{dP_A(n)}{dt} &= [\beta_A Z + I_{A,u}(n-1)]P_A(n-1) + I_{A,d}(n+1)P_A(n+1)
  - [\beta_A Z + I_{A,u}(n) + I_{A,d}(n)]P_A(n) \\
  &\quad - [\mu_A(N_A) + g_A(n)]P_A(n) + \frac{\delta_{n0}}{\tau_m}\sum_{n'}\sigma_m(n')P_T(n'),\\
\frac{dP_T(n)}{dt} &= (\text{the same transitions}) - \Big[\mu_T(N_T) + g_T(n) + \frac{1}{\tau_m}\Big]P_T(n)
  + \delta_{n0}\, r_{br}(t)\sum_{n'}\sigma_{br}(n')P_A(n'),\\
\frac{dZ}{dt} &= \sum_s \eta_s N_s S_s - \mu_Z (Z - Z_o).
\end{aligned}
```

`N_s` is the number of hosts and `S_s` their mean DOI. The functional responses are those of
ESM §1.4:

- **Within-host growth** (eq. 7): `I_u(n) = 2 max(0, μ_S + λ)(√(n+1) − 1)` and
  `I_d(n) = μ_S n − 2 min(0, μ_S + λ)(√(n+1) − 1)`, so `I′(0) = λ` and, for `λ > 0`,
  infections settle at the within-host capacity `K_S = 4λ(λ + μ_S)/μ_S²`
  ([`louca_load_rates`](@ref), [`louca_within_host_capacity`](@ref)).
- **Disease mortality** (eq. 6): `g_s(n) = (F_s − T_s) ln 2/(Y T_s) · n/(F_s − n)`, which halves
  yearly survival at the tolerance `T_s` and kills at the lethal DOI `F_s = 10⁴`.
- **Fecundity and metamorphosis**: `σ_br(n) = 1/(1 + n/T_br)`, `σ_m(n) = 1/(1 + n/T_m)`.
- **Demography** (eqs 5, 9): `μ_s(N) = −ln σ_s/Y + c_s N_s`, with `c_s` chosen so that the
  disease-free equilibrium under constant breeding is `(K_A, K_T)` ([`louca_competition`](@ref)).
- **Breeding** (eq. 4): a wrapped Cauchy pulse centred on `t = 0` that integrates to `r` a year,
  with full width at half maximum `τ_br` ([`louca_breeding_rate`](@ref)).
- **Transmission**: `β_s(t) = β̄_s [1 + b cos(2π(t − φ)/Y)]`.

**Tolerant tadpoles** (ESM §1.3, Table S1) have `T_T = F_T`, so `g_T ≡ 0`, and a DOI capped at
`S_T,max = F_T − 1`. Table 1 of the paper writes this case as `F_T = ∞`; we follow the ESM.

### Numerics

As in ESM §5, [`LoucaDOIProblem`](@ref) splits the DOIs `0 … F_s` into `m = 100` bins of width
`F_s/m` and divides every transition rate by the bin width. A host dies when its DOI would reach
`F_s`. Each run starts from a healthy population at `(K_A, K_T)`; after 1.5 years one adult with
DOI 10³ is added; the run ends after 10 years, or earlier when hosts fall below 2 (extinction) or
free plus within-host zoospores fall below 10 (clearance). The paper does not say when an environmental
reservoir starts. By default (`reservoir_onset = :exposure`) its supply `μ_Z Z_o` starts at
exposure, so the population is healthy until the second year, as §2b describes;
`reservoir_onset = :start` puts `Z(0) = Z_o`.

`solve(prob; method = :ode)` integrates the method-of-lines [`to_ode_problem`](@ref) with
Tsit5. The default `method = :semi_implicit` is a fixed-step scheme: backward Euler in the DOI
transitions and losses, density dependence and recruitment from the previous step, and an exact
update of `Z`. It stays positive and stable when the reservoir makes the infection rates huge
(`Z_o` up to 10¹⁵), where an explicit solver needs very small steps. The authors used an explicit
second-order Runge–Kutta scheme with `dt = 10⁻³–10⁻²` d. For the Table S1 parameter sets the two
methods agree to within 1%; in the Monte Carlo study the outcome with `dt = 0.05` d matched
`dt = 0.01` d in all 800 runs checked.

With `stochastic = true` the base mortality becomes random disturbances (ESM §5.3) at rate
`μ_s(0)/f_s`, each killing `Binomial(⌊N_s⌋, f_s)` hosts of the stage.

## Invasion: the reproduction quotient

Linearising a mean-field version of the model about the disease-free state gives (main eq. 3.1;
ESM eqs 18–19)

```math
Q_o = \frac{1}{\bar\mu_Z}\sum_s \frac{\eta_s \bar\beta_s \bar N_s}{\bar R_s - \lambda_s},
```

the zoospores released into the water that descend from one free zoospore. `R̄_s` is the
per-capita recruitment (turnover) rate of stage `s`. [`reproduction_quotient`](@ref) uses the
paper's approximations `N̄_s ≈ K_s`, `R̄_A ≈ K_T/(K_A τ_m)`, `R̄_T ≈ (r/Y) K_A/K_T` (ESM §5.2)
and returns `Inf` when `λ_s ≥ R̄_s` for some stage. [`louca_jacobian`](@ref) is the averaged
Jacobian `J̄_o` of ESM eq. (18), and [`louca_growth_rate`](@ref) its spectral bound or the
approximations of eqs (21) and (23). Over sampled parameter sets the spectral bound is
positive exactly when `Q_o > 1`.

`Q_o` ignores the on-host loss rate `μ_S`, because the mean-field model linearises the
within-host dynamics as `I(S) ≈ λS` with `λ = I′(0)`. In the binned model a newly infected host
sits at DOI 100 (one bin), where `I(n)` is far below `λn` when `μ_S` is large.

## Monte Carlo assessment of extinction risk

[`louca_sample`](@ref) draws a parameter set uniformly on each parameter's scale over the Table 1
ranges ([`louca_ranges`](@ref); log-uniform for `r, Z_o, μ_Z, μ_S,s, K_s, β_s, T_A, T_T, T_m,
T_br`) and rejects sets with a non-positive competition coefficient. [`louca_monte_carlo`](@ref)
runs the exposure experiment for many draws, and [`louca_risk_effects`](@ref) computes the risk
effect of each parameter: the slope of a least-squares line of the extinction indicator on the
parameter (on its sampling scale), times the parameter's range. [`louca_outcome_probabilities`](@ref)
bins the outcomes by `log10 Q_o` (Fig. 5) or by a parameter.

`scripts/louca2014_montecarlo.jl` runs 2000 draws per scenario (the paper ran 10⁴) and stores the
summary in `data/louca2014/reproduction.toml` ([`load_louca2014`](@ref)).

### Paper against ours

Our values are from 2000 runs per scenario. Paper columns: reservoir (`Z_o` from 10 to 10¹⁵) with
tadpoles that suffer or tolerate, then no reservoir. A dash is `|effect| < 0.1`.

| Quantity | Paper | Ours |
|---|---|---|
| Extinctions, `Z_o = 0`, suffer / tolerate | 25% / 21% | 21.5% / 18.4% |
| Risk effect of `Q_o` | +0.66, +0.83, +0.81, +0.71 | +0.90, +0.66, +1.13, +1.02 |
| `Z_o` | +0.66, +0.78 | +1.15, +1.12 |
| `T_A` | −0.38, −0.47, −0.17, −0.28 | −0.11, −0.05, +0.02, 0.00 |
| `μ_Z` | −0.11, −0.20, −0.44, −0.46 | −0.21, −0.15, −0.56, −0.50 |
| `β̄_A` | +0.32, +0.34, +0.38, +0.39 | +0.51, +0.55, +0.54, +0.54 |
| `K_T` (`Z_o = 0`) | −0.27, −0.35 | 0.00, −0.07 |
| `β̄_T` | +0.17, +0.27, +0.24, +0.18 | +0.29, +0.14, +0.41, +0.31 |
| `T_T` (`Z_o > 0`, suffer) | −0.24 | +0.01 |
| `K_A` (`Z_o = 0`) | −0.22, −0.12 | +0.06, +0.04 |
| `λ_A` | +0.17, +0.19, +0.18, +0.12 | +0.12, +0.11, +0.08, +0.09 |
| `μ_S,A`, `η_A`, `μ_S,T`, `T_br` | mostly —, −0.10 to +0.13 | within 0.1 |

What reproduces:

- More extinctions occur when tadpoles suffer than when they tolerate.
- `Q_o` is the best predictor without a reservoir, followed by `β̄_A` and `μ_Z`.
- Without a reservoir, Bd clears below `Q_o ≈ 0.1`, persists most often near `Q_o ≈ 1`, and
  extinction dominates above 10 (Fig. 5c,d).
- Extinction concentrates at high `β̄_A` and low `μ_Z`, and clearance at the opposite corner
  (Fig. 4c).

What does not reproduce:

- The paper finds that host vulnerability (`T_A`, `T_T`) and the carrying capacities matter. In
  our runs they hardly do: extinctions are mostly fast collapses that run hosts to the lethal DOI
  whatever their tolerance.

None of these explain the gap:

- **Step size.** 800/800 outcomes are equal at `dt = 0.01` d.
- **Number of DOI bins.** With `m = 1000`, extinction is 22.3% against 21.0%, and the risk
  effects are unchanged (`scripts/louca2014_hypothesis_bins.jl`).
- **Which parameter is log-uniform.** The ESM lists `r` twice; whether the second is `T_br`,
  `τ_m` or `τ_br` does not matter (`scripts/louca2014_hypothesis_scale.jl`).
- **When the reservoir starts.** It makes no difference whether the reservoir starts at
  exposure or at `t = 0`.

No other explanation could be tested, because the authors' C++ code was not deposited. The
registry lists each comparison (`louca2014_*` in `data/validation/targets.toml`).

## Parameter estimates (ESM §4)

- [`louca_within_host_loss`](@ref): μ_S ≈ 0.25 d⁻¹ from a logistic fit (λ = 0.38 d⁻¹,
  K_S = 14.9) to infected *Dendropsophus meridensis* (unpublished data).
- [`louca_tolerance_from_slope`](@ref): `T ≈ 1.73` from the fitted mortality slope
  `g′(0) = 0.0011`; [`louca_survival_curve`](@ref) is the survival curve of ESM eq. (31).
- [`louca_transmission_from_growth`](@ref): `β̃` between 7 × 10⁻⁶ and 2 × 10⁻³ d⁻¹ from the
  growth of *R. muscosa* epizootics (ρ ≈ 0.15 d⁻¹).
- [`louca_endemic_transmission`](@ref): `β_A` between 2 × 10⁻⁷ and 3 × 10⁻⁵ d⁻¹ for the endemic
  *Atelopus cruciger* population (ESM eq. 25).

All four, and the zoospore release and loss rates the ESM takes from Woodhams et al. (2008),
reproduce (validation targets `louca2014_*`).

## Relation to the other load models

- **Briggs et al. (2010)** ([`BriggsLoadProblem`](@ref)) follow the mean sporangial load of each
  host class. Loads grow by self-reinfection, without a within-host ceiling, until a lethal load
  `Smax = 10⁴` **sporangia**. Louca's DOI is zoospores (read as ZE per swab), the lethal DOI is
  also 10⁴, but on a different basis: `parameter_load(:louca2014, :FA)` equals Vredenburg's 10⁴ ZE
  crash load ([`parameter_load`](@ref)) and refuses comparison with Briggs' `Smax`. Louca's
  within-host growth is sub-exponential (`√n`, colonial growth on the perimeter of a localised
  infection) and levels off at `K_S`. Briggs' pool loses zoospores when they infect
  (`−γ N Z`); Louca's pool has no depletion term, so `β` plays the part of Briggs' `γ ν` without
  the loss.
- **Drawert et al. (2017)** ([`DrawertIBMProblem`](@ref)) make the Briggs model individual-based:
  a whole-number load per frog, like Louca's DOI, but with Briggs' exponential growth and
  demographic stochasticity in every event. Louca's stochastic version keeps the deterministic
  DOI distribution and makes only host mortality random.
- **Wilber et al. (2016)** ([`BdIPMProblem`](@ref)) use a discrete-time integral projection
  model over **log** load with fitted growth, survival and loss kernels. Louca's binned master
  equation is its continuous-time analogue on a linear DOI grid: bins of width 100 zoospores
  resolve high loads, where death happens, but merge every DOI below 100 into one bin.
- **The zoospore pool** ([`ZoosporePool`](@ref)) is the same linear `dZ/dt = shedding − μ_Z Z` in
  all of these; Louca adds an external reservoir `Z_o`.

## Parameters and problem

```@docs
LoucaParameters
LoucaDOIProblem
to_ode_problem(::LoucaDOIProblem)
louca_outcome
tadpoles_tolerate
```

## Functional responses

```@docs
louca_breeding_shape
louca_breeding_rate
louca_transmission
louca_fecundity_factor
louca_metamorphosis_success
louca_base_mortality
louca_disease_mortality
louca_load_rates
louca_load_speed
louca_within_host_capacity
louca_within_host_loss
louca_competition
louca_turnover
```

## Invasion and estimates

```@docs
reproduction_quotient
louca_jacobian
louca_growth_rate
louca_tolerance_from_slope
louca_survival_curve
louca_transmission_from_growth
louca_endemic_transmission
```

## Monte Carlo

```@docs
louca_ranges
louca_sample
louca_monte_carlo
louca_risk_effects
louca_outcome_probabilities
louca2014_path
load_louca2014
```
