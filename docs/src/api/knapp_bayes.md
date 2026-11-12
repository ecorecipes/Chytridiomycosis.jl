# Bayesian fit to field data

```@meta
CurrentModule = Chytridiomycosis
```

!!! note "New analysis"
    The load equation, its priors and these posteriors are new work in this package, not a
    published analysis. The data are from Knapp et al. (2022), CC BY 4.0.

The maximum likelihood fit of [Fitting to field data](fitting.md) gives point estimates, AIC
and profile likelihoods. The Bayesian fit puts a weakly informative prior on every parameter
of a [`KnappLoadModel`](@ref), samples the posterior with NUTS and compares mechanisms by
their leave-one-out predictive accuracy (PSIS-LOO).

- [`KnappPosterior`](@ref) adds the prior to the likelihood. It implements the
  LogDensityProblems interface on the log parameters.
- [`turing_model`](@ref) turns it into a Turing.jl model. It is defined in the package
  extension `ChytridiomycosisTuringExt`, loaded by `using Turing`. The likelihood uses a
  fixed-step integrator, so it is differentiable and NUTS can use ForwardDiff gradients.
- [`pointwise_loglikelihood`](@ref) gives one log-likelihood per swab and draw.
  [`psis_loo`](@ref), [`waic`](@ref) and [`elpd_difference`](@ref) compare mechanisms with it.

```julia
using Chytridiomycosis
using Turing: Turing, NUTS, MCMCThreads, sample, InitFromParams

m = KnappLoadModel(:applied)
post = KnappPosterior(m)                         # default prior: knapp_prior(m)
x0 = Chytridiomycosis.maximum_likelihood(m).x    # Turing also exports maximum_likelihood
chn = sample(turing_model(post), NUTS(1000, 0.99), MCMCThreads(), 1000, 4;
    initial_params = fill(InitFromParams((; x = x0)), 4))
X = reshape(Array(chn), :, length(x0))           # draws × log parameters
loo = psis_loo(pointwise_loglikelihood(m, X))
```

## Priors

Each parameter has a lognormal prior on its natural scale, so its logarithm, which the sampler
works with, has a normal prior. The medians are orders of magnitude that are plausible before
seeing the swabs. The spreads cover several orders of magnitude for the rates, loads and
*J. lividum* levels, and a factor of about 20 either way for half-lives and standard
deviations. The within-frog rate `r0 = ηνf − σ = −0.156` per day stays fixed at the Briggs et
al. (2010) compendium value, as in the maximum likelihood fit, so that the two fits stay
comparable. Untreated frogs sit at equilibrium and fix only `ι/(−r0)`. The rate of return
itself would be informed only by the treated frogs, whose return is also shaped by the
extra clearance, so a free `r0` would trade off against `Δc`, `h` and `τ` or `K`. The table
is in the docstring of [`knapp_prior`](@ref).

## Results

The run below is `scripts/knapp_bayes.jl`: 4 chains of 1000 draws after 1000 adaptation
steps, started near the maximum likelihood estimate with jitter (standard deviation 0.3 on
the log scale), on 4 threads. Target acceptance was 0.99 for `applied`, which had 16
divergent transitions at 0.9, and 0.9 for the others.

**Convergence.** R-hat and effective sample sizes are from MCMCDiagnosticTools on the log
parameters, worst over parameters.

| Mechanism | NUTS time | Divergences | Max R-hat | Min bulk ESS | Min tail ESS |
|---|---|---|---|---|---|
| `none` | 25 s | 0 | 1.004 | 3771 | 2692 |
| `decay` | 161 s | 0 | 1.004 | 1447 | 1588 |
| `applied` | 553 s | 0 | 1.002 | 1399 | 1614 |
| `bacteria` | 74 s | 84 | 1.012 | 458 | 181 |

`bacteria` does not sample cleanly, at 0.99 either (154 divergences). Its posterior turns the
effect off (`Δc` near 0), so `h` and `K` are free and form a funnel. Its predictions, and so
its LOO, are the same as those of `none` (elpd −5755.7 and −5755.9 in two runs).

**Comparing mechanisms.** PSIS-LOO leaves out one swab at a time (371 swabs; every Pareto
`k` is below 0.5, so the estimates are reliable). WAIC agrees to 0.1.

| Mechanism | ΔAIC (ML) | elpd_loo | p_loo | Δelpd vs `applied` (SE) |
|---|---|---|---|---|
| `none` | 178.8 | −5754.9 | 9.9 | −89.0 (11.0) |
| `decay` | 6.2 | −5668.9 | 13.2 | −2.9 (2.7) |
| `applied` | 0 | −5666.0 | 13.6 | 0 |
| `bacteria` | 184.8 | −5755.7 | 10.4 | about −90 (11) |

The improvement of `applied` over `decay` is about one standard error. Most of it comes from
treated frogs on day 37 (1.5 of the 2.9), where `applied` holds loads lower for longer.

**Posteriors.** Medians and 95% intervals on the natural scale. The *J. lividum* parameters
are the same under every mechanism. The well-determined parameters have posterior medians
close to the maximum likelihood estimates. `Δc` and `h` do not: under `applied` the maximum
likelihood values are 1.8 × 10⁵ and 1.1 × 10³, against posterior medians of 3.6 × 10⁵ and
4.6 × 10³, because the posterior spreads along the ridge between them.

| Parameter | `decay` | `applied` | Prior 95% |
|---|---|---|---|
| `ι` (copies per day) | 1.15 × 10⁵ (0.73–1.8 × 10⁵) | 1.08 × 10⁵ (0.71–1.7 × 10⁵) | 28 – 3.6 × 10⁶ |
| `L0` (copies) | 10.4 (4.5–23) | 10.4 (4.4–24) | 0.28 – 3.6 × 10⁴ |
| `σ` | 2.99 (2.77–3.24) | 2.96 (2.74–3.21) | 0.28 – 14 |
| `Δc` (copies per day) | 4.1 × 10⁵ (1.6 × 10⁵ – 2.7 × 10⁶) | 3.6 × 10⁵ (1.1 × 10⁵ – 4.4 × 10⁶) | 28 – 3.6 × 10⁶ |
| `h` (copies) | 2.0 × 10³ (420 – 1.4 × 10⁴) | 4.6 × 10³ (270 – 6.4 × 10⁴) | 2.8 – 3.6 × 10⁵ |
| `Δc/h` (per day) | 219 (88–553) | 90 (37–509) | |
| `τ` (days) | 39 (18–108) | | 0.74 – 265 |
| `K` (genomic equivalents) | | 0.34 (0.040–2.5) | 0.028 – 3.6 × 10³ |

What the posterior adds to the maximum likelihood fit:

- **`h` is weakly, not un-, identified.** It lies along a ridge with `Δc` (correlation of the
  logs 0.86 for `decay`, 0.91 for `applied`). The data fix the low-load clearance rate
  `Δc/h`, and they rule out `h` below a few hundred copies. The posterior standard deviation
  of `log h` is 0.30 (`decay`) and 0.45 (`applied`) of the prior's, so the prior still
  shapes its upper tail.
- **`K` is below a few genomic equivalents.** P(`K` < 1) = 0.85 and P(`K` < 3) = 0.99. The
  lower end of the interval is near the prior's 2.5% point, so the data set only an upper
  bound.
- **`τ` is bounded on both sides**, 18–108 days, narrower above than the profile interval
  (about 17–200 days). P(`τ` > 56 days, the length of the trial) = 0.24.
- **Treated frogs' median load** under `applied`: 1.9 × 10³ (980–3.3 × 10³) on day 12,
  2.7 × 10³ (1.5–5.0 × 10³) on day 37 and 1.0 × 10⁵ (3.3 × 10⁴ – 2.4 × 10⁵) on day 56.
- **The mechanisms are not separable.** ΔAIC of 6.2 reads as a clear preference. The LOO
  difference, 2.9 ± 2.7, shows it is about one standard error.

Swabs are treated as independent, as in the maximum likelihood fit, so all of these
intervals are too narrow if swabs from the same frogs are correlated.

## API

```@docs
KnappPosterior
knapp_prior
turing_model
observed_swabs
pointwise_loglikelihood
psis_loo
waic
elpd_difference
```
