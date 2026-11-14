# Extinction risk with a zoospore pool (Wilber 2017)

```@meta
CurrentModule = Chytridiomycosis
```

Wilber, Knapp, Toothman and Briggs (2017, *Ecology Letters* 20:1169–1181,
[10.1111/ele.12814](https://doi.org/10.1111/ele.12814)) extend the host–parasite integral
projection model of [Wilber et al. (2016)](wilber2016.md) to ask whether extinction of
*Rana muscosa* populations depends more on transmission or on host resistance and tolerance.
They add an environmental zoospore pool fed by infected adults and tadpoles, three tadpole
stages, a seasonal temperature cycle, a yearly demographic update with density-dependent
recruitment, and demographic stochasticity. The implementation here was checked against the
paper, the authors' code ([github.com/mqwilber/host_extinction](https://github.com/mqwilber/host_extinction);
their March 2017 preprint and appendix in that repository stand in for Appendices S1–S4,
which returned HTTP 403) and their data (Dryad
[10.5061/dryad.p614h](https://doi.org/10.5061/dryad.p614h), CC0, bundled).

## Additions to the IPM machinery

Two generic pieces are added to the [IPM machinery](wilber2016.md):

- `R0` when infected hosts also shed into an environmental pool. The pool is one more
  infected state, so the next-generation matrix of Klepac and Caswell (2011) has rank one and
  `R0 = [t; β]ᵀ (I − U)⁻¹ [s0 G0; 0]` (Appendix S2): the host term plus a pool term.
- One step of a matrix model with demographic stochasticity (Caswell 2001; Schreiber and
  Ross 2016): individuals move by multinomial draws on the columns of the transition matrix,
  with death as the remainder, and births or shed pathogens are Poisson.

```@docs
reservoir_reproduction_number
demographic_projection
demographic_projection!
```

## The hybrid model

The state is `u = [T1, T2, T3, S, n_1, …, n_m, Z]`: tadpoles of years 1–3, uninfected adults,
infected adults per load bin (30 bins on `[−5, 18]` ln ZE per swab) and zoospores. Every
3-day step (eqns 1–3; Appendix S3):

```math
\begin{aligned}
T_i(t+1) &= s_{T_i} T_i(t), \qquad
S(t+1) = s_0 (1-\phi) S + \textstyle\sum_x s(x)\, l(x,T)\, n(x), \\
n(x',t+1) &= \textstyle\sum_x s(x)(1-l(x,T)) G(x',x,T)\, n(x) + s_0 \phi S\, G_0(x',T), \\
Z(t+1) &= \nu(T) Z + \mu_A \textstyle\sum_x e^x n(x) + \mu_T \sum_i T_i - \psi ,
\end{aligned}
```

with `φ = 1 − exp(−Λ Δt)` from one of the transmission functions of Table 2 and the
zoospore survival `ν(T)` from Woodhams et al. (2008). Temperature follows a cosine from
4 °C in winter to 20 °C in summer (period 366 days). Once a year, in the step at which the
rising temperature reaches 11.79 °C, adults reproduce (`pA λ` year-1 tadpoles per
surviving adult), year-1 and year-2 tadpoles age or metamorphose, and metamorphs join the
uninfected adults with Ricker survival `exp(−A/K)`. With demographic stochasticity, every
frog and zoospore moves by a multinomial draw on the discretised kernel, and tadpoles and
shed zoospores are Poisson. The stochastic simulator draws from a seeded `rng`, like
[`DrawertIBMProblem`](@ref); the deterministic map lowers to a `SciMLBase.DiscreteProblem`.

```@docs
Wilber2017Parameters
zoospore_survival
infection_probability(::Wilber2017Parameters)
annual_demography_matrix
wilber2017_temperature
HybridIPMProblem
HybridIPMSolution
hybrid_summary
to_discrete_problem(::HybridIPMProblem)
extinction_curve
```

[`IPMKernel`](@ref) has a method for `Wilber2017Parameters`: the within-season adult kernel.

## R0 and invasion (Fig. 2)

```@docs
basic_reproduction_number(::Wilber2017Parameters)
```

## Sensitivity of extinction (Fig. 3)

```@docs
stasis_load
extinction_sensitivity
sensitivity_weights
```

## Data and refits (Tables 1–2)

The within-season vital rates are refitted to the Wilber 2016 temperature experiment
(`data/wilber2016`) and the constant-pool transmission functions to the mesocosm experiment
(`data/wilber2017`), with the authors' Stan priors. Posterior means and WAIC come from
quadrature (up to three parameters) or importance sampling, without MCMC.

```@docs
wilber2017_path
load_wilber2017
Wilber2017Transmission
pointwise_loglikelihood(::Wilber2017Transmission, ::AbstractVector)
Wilber2017Regression
posterior_summary
posterior_draws
```

## What reproduces

Paper and ours are compared in `test/wilber2017.jl`, the validation targets `wilber2017_*`
and `scripts/wilber2017_fig2.jl` and `scripts/wilber2017_sensitivity.jl`. Fig. 2 and Fig. 3
were decoded from the vector graphics of the published PDF; the extinction curves of
Appendix Fig. S5 from Fig. 3A of the authors' preprint.

| Result | Paper | Ours |
|:--|:--|:--|
| survival `b00`, `b10` (posterior means, z-scored load) | 5.295, −2.595 | 5.273, −2.567 (Stan Monte Carlo error ≈ 0.02–0.03) |
| growth, loss, initial-load rows (Table 1) | as printed | the Wilber 2016 REML/ML fits, to every digit |
| transitions in the transmission fits | 333 | 333 (124 infections) |
| constant pool `β0`; WAIC | 0.0807; 401.1 | 0.08068; 401.11 |
| density-dependent, constant pool `β0`, `β1`; WAIC | 0.0418, 0.0525; 336.7 | 0.04178, 0.05264; 336.80 |
| frequency-dependent, constant pool `β0`, `β1`; WAIC | 0.0428, 0.551; 336.6 | 0.04275, 0.5496; 336.53 |
| disease-free annual growth rate | ≈ 1.46 | 1.4633 |
| initial tadpoles for 10 adults | 85, 12, 3 | 85, 12, 3 |
| disease-free adults | 10–15 | 12.7–14.1 over the year |
| Fig. 2c ln R0 with the pool, S = 4, at 4 and 15 °C; without it at 20 °C | 0.291, 2.537; 1.478 | median over 1 000 draws 0.294, 2.536; 1.511 (plug-in 0.290, 2.612; 1.586) |
| Fig. 2b ln R0 at 15 °C, S = 1, with and without the pool | 1.167, 0.090 | 1.154, 0.091 (plug-in 1.231, 0.167) |
| extinction within 8 years | c. 50% | 0.45–0.54 for the three transmission functions (500 runs each) |
| extinction by year 10: pool only, frequency-, density-dependent; high zoospore decay | 0.76, 0.76, 0.74; 0.43 | 0.75, 0.72, 0.71; 0.42 (500 runs each, 25-year runs) |
| logistic-regression weights (Fig. 3c, density-dependent): `b11`, `b21`, `b00`, `b10` | 6.6, 1.8, −1.5, −0.5 | 6.1, 1.4, −1.0, −0.7 (1 000 runs); `b11` first for all three functions |

The dynamic-pool transmission fits (Table 2, rows 4–6) need the tank-level series of the raw
file `DensExpt.csv`, which is only in the authors' repository without a licence; they are
not refitted, and the published values are used. The random forest of Fig. 3 is not
reproduced.

### Printed and coded details that differ

- **Table 1 mixes two fits.** The survival row is a Bayesian posterior mean on the z-scored
  load (whose mean 7.186 and sd 2.961 ln ZE are not printed). The growth, loss and
  initial-load rows are the Wilber (2016) REML and ML fits, but the simulations used Bayesian
  posterior means (`bayesian_parameter_estimates.rds`): for example `b11` 0.796 rather than
  0.799, `v01` 5.84 rather than 5.92 and `v03` 0.71 rather than 0.59.
- **R0 (Fig. 2).** The figure is the median over 500 posterior draws, which is lower than
  `R0` at the point estimates. The authors' R0 code omits `Δt` (it uses `β` per day as a
  probability per step), sheds `ln(eˣ + 1)` zoospores per infected adult (Appendix S2 says
  `f x`; eqn 3 and the simulation use `eˣ`) and sets `s0 = 1`. `as_coded = false` gives `R0`
  of the model as simulated, which is far larger.
- **Zoospore removal.** The text says the uptake of zoospores by newly infected frogs (`ψ`)
  "had very little effect … we do not consider it further"; the code includes it
  (`zoospore_removal = true`).
- **Births.** Appendix S4 says Poisson births with mean `A pA λ`; the code uses `pA λ` times
  the adults that survive the step.
- **Transmission priors.** The latent-pool fits give zoospores a death rate with prior
  Normal(0.3, 0.03) per day "from Woodhams et al. 2008"; the simulation uses Woodhams's
  measured rates, 0.065–0.19 per day over 4–20 °C.
- **Loss-of-infection priors.** The Table 1 footnote gives Normal(0, 5) priors to all `b`;
  the Stan model of `l(x)` has its priors commented out (flat).
- The code's extinction simulations use 25 years (the YAML file says 20).

## How it relates to Wilber 2016 and the Briggs model

The within-season kernel is the 2016 IPM with the survival function refitted on the z-scored
load (on the raw scale, intercept 11.59 and slope −0.876, close to the 2016 deposit) and
environmental infection replaced by the zoospore pool. Like the Briggs et al. (2010) model, it
has tadpoles that are always infected and feed a pool, and density-dependent recruitment; it
differs in using fitted load regressions instead of sporangial growth to a lethal threshold.
The pool makes transmission nearly independent of host density, which is why the three
transmission functions give similar extinction curves.
