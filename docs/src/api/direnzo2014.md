# Load and zoospore output of *Atelopus zeteki* (DiRenzo et al. 2014)

```@meta
CurrentModule = Chytridiomycosis
```

This is a reference implementation of the statistical models of DiRenzo et al. (2014),
*PLOS ONE* 9:e93356
([doi:10.1371/journal.pone.0093356](https://doi.org/10.1371/journal.pone.0093356)), checked
against the published article. Eight *Atelopus zeteki* were inoculated with Bd JEL 423. Five
were naive and three had previously been exposed to an attenuated strain. Each frog was
swabbed (infection intensity, zoospore genomic equivalents) and soaked (zoospore output, per
minute) from day 7 until death.

## Model

Each response was fitted with a linear mixed model of its natural log against time, with a
first-order autoregressive (AR(1)) correlation:

```
ln y_ij = β0 + b_i + r t_ij + ε_ij,   b_i ~ N(0, σ_id²),   cor(ε_ij, ε_ik) = ϕ^|j−k|
```

The mean trajectory is therefore `y(t) = y0 e^{rt}`. The Methods write it as
"dN/dt = y0 e^{rt}", which is the integrated form labelled as a derivative.

!!! warning "What the paper prints"
    The paper prints only the slopes, as a "daily percent increase": 35.4% (SE 0.05) for
    intensity and 13.1% (SE 0.04) for output.

    - **These are 100 × r on the ln scale.** `r/SE` gives 7.08 and 3.28, which match the
      printed t50 = 7.018 and t46 = 3.164. The actual increase per day is
      `100(eʳ − 1)`: 42.5% and 14.0%.
    - **Not printed:** the intercepts, the random-effect and residual SDs, and the AR(1)
      coefficient. No time series were deposited either, so the fitted trajectories need a
      user-supplied `y0`.
    - **Our choices.** The random effect is a random intercept per frog. The AR(1) runs over
      the order of observations, as in nlme's default `corAR1()`.

```@docs
DiRenzoParameters
direnzo_trajectory
daily_percent_increase
direnzo_t_statistic
DiRenzoSeries
DiRenzoLoadModel
direnzo_covariance
```

[`maximum_likelihood`](@ref) and [`parameter_names`](@ref) also have methods for a
[`DiRenzoLoadModel`](@ref).

## Table 1 and the paper's inconsistencies

Table 1 gives the survival time and the terminal values of each frog. It is bundled in
`data/direnzo2014` under CC BY 4.0.

- **Reproduced:** the mean intensity at death, 4,334,422, and its SE in the Results,
  1,156,576.
- **Not reproduced:**
  - **Survival.** The abstract says all frogs died "between 22 and 33 days" (mean 28.88).
    The table's times run from 18 to 33 days, with mean 25.625.
  - **The abstract's SE.** The abstract gives the SE of the intensity at death as
    1,236,431.
  - **Zoospore output.** The text gives a mean output of 23.55 with a maximum of 172.61.
    The table gives 23.29 with a maximum of 170.6.
  - **Log-rank test.** On the table's survival times it gives χ² = 1.30, not 0.7.

```@docs
DiRenzoFrog
direnzo2014_path
load_direnzo2014
direnzo_summary
direnzo_logrank
```
