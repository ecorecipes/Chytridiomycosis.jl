# Several hosts, one zoospore pool (Daversa et al. 2022)

```@meta
CurrentModule = Chytridiomycosis
```

Daversa et al. (2022, *Am Nat* 200(4), [doi:10.1086/720638](https://doi.org/10.1086/720638))
ran mesocosms of fire salamander, spiny toad and Iberian tree frog larvae, which they inoculated
with Bd. They fitted a model in which every species draws on, and sheds into, one zoospore
pool. For tank `j` (time in days, numbers per 80 L tank):

```math
\begin{aligned}
\frac{dH_i}{dt} &= -\mu_i H_i,\\
\frac{dI_i}{dt} &= \beta_i (H_i - I_i) Z - I_i(\mu_i + \sigma_i),\\
\frac{dZ}{dt} &= \sum_k \lambda_k I_k + \Phi(t) - Z\Big(\gamma + \sum_k H_k \beta_k\Big).
\end{aligned}
```

The equations were checked against the PDF and against the authors' R code
([Zenodo 5068208](https://doi.org/10.5281/zenodo.5068208)).

[`MultiHostPoolProblem`](@ref) takes any number of species. The inoculations Φ(t) are pulses
given by `doses`.

## Reproduction numbers

The contribution of species `i` is

```math
R_{0,i} = \frac{\beta_i \lambda_i H_i}{(\mu_i + \sigma_i)\,(\gamma + \sum_k H_k \beta_k)},
\qquad R_{0,\mathrm{TOT}} = \sum_i R_{0,i}.
```

The paper (p. 6) prints `μ_i` alone in the denominator. The authors' code uses `μ_i + σ_i`, and
so does the next-generation matrix. Only `μ_i + σ_i` reproduces Fig. 4b. The total is the sum
of the contributions, and it is also the spectral radius: all transmission passes through the
one pool, so the host-to-host next-generation matrix has rank one.
[`next_generation_matrix`](@ref) computes the NGM numerically from the right-hand side, and the
tests check that it equals the sum. The method of
[`basic_reproduction_number`](@ref) for `DaversaParameters` returns that sum.

## Parameters

The paper reports its estimates only as a figure (Fig. 4a), and the supplement could not be
obtained. `DaversaParameters(:daversa2022)` therefore holds posterior medians digitised from
that figure, which the compendium marks as unverified. At these values:

- the fit reproduces the paper's r = 0.70 and R² = 0.50;
- per-tank R0 tracks the Fig. 4c–d medians (r = 0.9996), lying about 7% below them.

See `data/validation/targets.toml` for the targets.

```@docs
DaversaParameters
MultiHostPoolProblem
daversa_doses
species_reproduction_numbers
next_generation_matrix
tank_reproduction_numbers
```

## The mesocosm data

The Dryad deposit ([10.5061/dryad.2v6wwpznk](https://doi.org/10.5061/dryad.2v6wwpznk), CC0)
ships in `data/daversa2022/`. It is prepared as in the authors' code.

```@docs
daversa2022_path
MesocosmTank
load_daversa2022
observed_infected
predicted_infected
daversa_loglikelihood
```
