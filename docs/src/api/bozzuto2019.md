# Seasonal Bsal and a Floquet R0 (Bozzuto & Canessa 2019)

```@meta
CurrentModule = Chytridiomycosis
```

Bozzuto & Canessa (2019, *Glob Ecol Conserv* 17:e00551,
[doi:10.1016/j.gecco.2019.e00551](https://doi.org/10.1016/j.gecco.2019.e00551)) asked how
seasonal cycles change the invasion threshold of *B. salamandrivorans* in fire salamanders,
and how much management is needed to keep a population free of it. Time `t` is in decimal
years from the start of January, and densities are per hectare:

```math
\begin{aligned}
\frac{dS}{dt} &= v(t) - dS - c(t)p(t)SI,\\
\frac{dI}{dt} &= c(t)p(t)SI - (\delta(t) + d)I.
\end{aligned}
```

Each seasonal coefficient has the form `θ(t) = θ₀(1 − θ₁ cos 2π(t − ℓθ))`. The proportion of
contacts that transmit is `p(t) = p₀(1 − p₁ cos 2πt)` (the paper sets `p₁ = 0`). Bsal-induced
mortality is `δ(t) = δ₀(1 − δ₁ cos 2πt)`, and natural mortality `d` is constant. The contact
rate has two peaks, in spring and in autumn:

```math
c(t) = \frac{c_0}{m'}\,(1 - c_1\cos 2\pi t)\,(1 - c_2\cos 2\pi(t - \ell_c)).
```

Larval inflow `v(t)` is a step, `v` between `t₁` and `t₂` (late summer). The equations and
values were checked against the PDF and the online Appendix A.

[`SeasonalBsalProblem`](@ref) lowers the model to an `ODEProblem`. `seasonal = false` gives
the non-seasonal model.

## Two things the PDF does not make obvious

- **The contact rate is rescaled.** Table 1 prints `c(t)` without the factor `1/m'`.
  Appendix A1 divides by `m'`, half the sum of the minimum and maximum of `c(t)/c₀`, so that
  `c₀` lies midway between the lowest and highest contact rate. Eq. (2) holds only for the
  rescaled form. `scaled_contact = false` gives the Table 1 formula as printed.
- **Table 1 prints the wrong `c₂` for the hibernating population.** It prints 0.90. The
  Fig. A4 caption gives 0.95, and only 0.95 reproduces the Floquet-R0 of 8.27, the numerical
  8.33 and the elasticities of Fig. A5. The package uses 0.95. `as_printed = true` uses 0.90,
  which gives 8.67.

## Floquet-R0

The infected equation, linearised at the disease-free state, is `dI/dt = J(t) I` with
`J(t) = c(t)p(t)S(t) − d − δ(t)`. Its Floquet exponent `ξ` is the growth of `log I` over one
year, and Bsal cannot invade when `ξ < 0`. Rearranging `ξ < 0` as a ratio gives the paper's
Floquet-R0, and with `S = k₀` and `p₁ = 0` it is eq. (2):

```math
\text{Floquet-}R_0 = R_0\,\frac{2 + c_1c_2\cos 2\pi\ell_c}{m},
\qquad R_0 = \frac{c_0p_0k_0}{d + \delta_0},\qquad m = 2m'.
```

Here `R₀` is the paper's non-seasonal reproduction number. The infected subsystem has only
one compartment, so the Floquet exponent is simply a time average. Three consequences
follow:

- The paper's numerical Floquet-R0, `∫ c p S dt / ∫ (d + δ) dt`, is exactly the periodic R0
  of Bacaër and Guernaoui (the spectral radius of the next-generation operator).
  `method = :spectral` computes it the general way, from Floquet multipliers, and agrees to
  within 1e-7.
- The analytic Floquet-R0 equals the R0 of the model with every coefficient replaced by its
  annual mean ([`averaged_reproduction_number`](@ref)). The roughly 7% gap between
  Floquet-R0 and the paper's `R₀` exists because `c₀` is the midrange of the contact rate,
  not its mean. Seasonality as such, the covariance of contacts with host density, adds
  less than 1% (8.33 against 8.26 for the hibernating population).
- `δ₁` has no effect on the result, as the paper notes.

| | Paper | Package |
|---|---|---|
| Hibernation: R0, Floquet-R0 (eq. 2), numerical | 8.88, 8.27, 8.33 | 8.880, 8.272, 8.328 |
| Hibernation, high activity | 8.22 | 8.221 |
| Aestivation: R0, Floquet-R0, numerical | 2.13, 1.92, 1.92 | 2.095, 1.888, 1.892 |
| Aestivation, high activity | 1.80 | 1.776 |
| Largest increase and decrease (Fig. 1) | +33%, −33% | +33.3%, −33.3% |
| Target proportions ϱF, ϱ (hibernation) | ≈12%, 11% | 12.1%, 11.3% |
| 1/max R0(t), hibernation and aestivation | ≈5%, 20% | 4.8%, 20.0% |
| Elasticity to ℓc (Fig. A5) | ≈2.37, 1.58 | 2.367, 1.581 |

Every aestivation value in the package is about 1.6% below the paper's, while the seasonal
ratios agree (Floquet-R0/R0 is 0.901 in both). One of the aestivation inputs in Table 1
(`c₀`, `k₀` or `δ₀`) must differ from the value the authors used, but the paper does not
say which. These targets are recorded as discrepant in `data/validation/targets.toml`.

```@docs
BozzutoParameters
SeasonalBsalProblem
floquet_exponent
floquet_reproduction_number
averaged_reproduction_number
seasonal_factor
target_proportion
floquet_elasticities
instantaneous_reproduction_number
```

`basic_reproduction_number(p::BozzutoParameters)` returns the non-seasonal `R₀`.

## Seasonal coefficients

```@docs
contact_rate(::BozzutoParameters, ::Any)
contact_scaling
transmission_probability
bsal_mortality
larval_inflow
disease_free_susceptibles
```

## Floquet tools for any periodic model

These functions take matrix-valued functions of time. Appendix A4 of the paper uses the same
construction for a model with an environmental zoospore reservoir.

```@docs
floquet_multipliers
periodic_reproduction_number
```
