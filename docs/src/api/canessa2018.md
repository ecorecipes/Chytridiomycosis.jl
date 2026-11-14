# Bsal management decisions with an IPM (Canessa et al. 2018)

```@meta
CurrentModule = Chytridiomycosis
```

Canessa et al. (2018, *J Appl Ecol* 55:1987–1996,
[doi:10.1111/1365-2664.13089](https://doi.org/10.1111/1365-2664.13089)) asked how, if at all,
a *B. salamandrivorans* outbreak in a fire salamander population could be managed. They
compared eight actions over three months, using an integral projection model (IPM) whose
functions were elicited from ten experts. Uncertainty was carried through 1,000 random
parameter sets. The equations were checked against the review copy of the paper, the only
version available (main text, Table 1, Figs 1–4, Supporting Information eqs S1–S12 and Table
S1). The 1,000 sets are the authors' Dryad deposit
([doi:10.5061/dryad.0jh5p](https://doi.org/10.5061/dryad.0jh5p), CC0), bundled as
`data/canessa2018/parameter_draws.csv`.

## The model

There are susceptible adults `S(t)` and infected adults `I(z, t)` indexed by log-load
`z = ln(GE)`. Each step is three days, and there are 30 steps. There is no recruitment, no
density dependence and no recovery (eq. S2):

```math
\begin{aligned}
S(t+1) &= S(t)\,s_0\,\bigl(1-\phi(t)\bigr),\\
I(z',t+1) &= S(t)\,s_0\,\phi(t)\,G_0(z') + \int_l^u I(z,t)\,s(z)\,G(z',z)\,dz,\\
\phi(t) &= 1-\exp\Bigl(-\int_l^u \beta(z)\,I(z,t)\,dz\Bigr).
\end{aligned}
```

- **Survival:** `s(z) = logistic(b₀ + b₁z)` (eq. S4).
- **Transmission:** `β(z) = β₀ logistic(b₂ + b₃z)` (eq. S5).
- **Initial load:** `G₀ = N(ln 200, 0.5)` (eq. S7).
- **Growth:** `z'` is normal with mean `b₄ + b₅z + b₆z²` (eq. S9) and variance
  `σ²_G exp(b₇z)` (eq. S10, borrowed from Bd).

The integrals use the midpoint rule on 300 points in `[-10, 15]`.
[`CanessaIPMProblem`](@ref) iterates the model. Its reproduction number is
`R₀ = S₀ s₀ ∫ β(z) [(I − P)⁻¹ G₀](z) dz`, the infections an adult causes over its infected
life (Wilber et al. 2016); `basic_reproduction_number(p::CanessaParameters)` computes it.

Each of the 1,000 runs refits the functions to one archived draw ([`CanessaDraw`](@ref)):

- logistic regressions of five elicited survival and transmission probabilities (loads 0,
  10, …, 10⁴ GE, with load 0 at `z = −1`);
- a quadratic regression of a Beverton–Holt load series for eq. S9
  ([`canessa_growth_coefficients`](@ref)).

[`canessa_ensemble`](@ref) runs all the draws. [`canessa_summary`](@ref) reduces the runs to
the medians and 95% intervals of Fig. 2, and [`canessa_dispersal`](@ref) gives the movement
distances of Fig. 3.

## What the PDF does not say, and what does not agree

- **Table S1 is not the ensemble.** Its coefficients are fits to the experts' best
  estimates. The draws come from the pooled elicitation and lie far from them. For example,
  survival at 100 GE is 0.62 in Table S1 against a draw median of 0.895. The Table S1 point
  model has `R₀ ≈ 5.1`, while the ensemble median is 9.4.
- **Growth coefficients.** The regression behind eq. S9 is only described. Ten
  Beverton–Holt steps from `L₀ = 10`, with `Lmax = 10⁴` GE and 25 zoospores per zoospore in
  5 days, reproduce Table S1's `b₄`–`b₆` to four decimals. These best estimates are recorded
  in the compendium as inferred (`verified = false`).
- **The cap at `Lmax`.** The SI formula only sets the expected next load to `ln Lmax` once
  the cap is reached. Read literally, half of the loads at the cap then exceed it, where
  survival is lower, and every R0 quantile is about 14% below the paper's (median 8.2).
  The default, `cap = :truncate`, does not let loads exceed `Lmax`. This reproduces the R0
  values of Fig. 2: median 9.43 (95% interval 5.60–14.35) against 9.6 (5.8–14.3). Thinning
  by 50% gives 4.72 against 4.8, and thinning by 90% gives 0.94 against 0.93.
  `cap = :printed` gives the literal reading.
- **Units of `β₀` and the final densities.** Eq. S6 gives the printed `β₀ = 0.007` only
  when it is solved with the number of adults on Stegen's 475 m × 10 m transect
  (503 × 0.475 = 238.9) rather than the density 503 ha⁻¹, which gives 0.0039. The archived
  draws agree: read this way, every `(β₀, S₀)` pair implies a valid `s̄` (mean 0.544,
  printed 0.542), whereas per hectare 13% would need `s̄ < 0`. So `β₀` is per adult on the
  transect. Simulating in those units (`units = :transect`, the default: `S₀ × 0.475`
  susceptibles and one infected adult) reproduces Fig. 2's final densities: 0.951 for
  action e (plotted 1.0) and 0.081 for f (0.083). Per hectare (`units = :density`,
  `I₀ = 0.01` ha⁻¹) they are 0.69 and 0.006, because transmission per susceptible is then
  2.1 times higher. The paper's R0 values keep `S₀` per hectare, and so does
  `basic_reproduction_number`; the simulated outbreaks therefore spread as if R0 were 0.475
  times the reported value. This is an inference, not printed.
- **Printed inconsistencies.**
  - Table S1's `s₀ = 0.9906` is not `logistic(b₀) = 0.9875`.
  - `s₀` is said to match an 8-year life span, but that span gives 0.99897 per 3 days.
  - The simplified model's `R₀` is 7.62, not 7.9.
  - "Reducing transmission and host survival by 87% each (Fig 4a)" should read 54%. The
    87% is the transmission-only intercept of Fig. 4a, whose axes are transmission and
    density.
  - The text says the probiotic increases R0, as here (+2.6%), but Fig. 2 plots it below
    no action.

All of these are validation targets (`canessa2018_*` in `data/validation/targets.toml`).

## Results

The medians and 95% intervals over the 1,000 draws (`cap = :truncate`,
`units = :transect`):

| Action | R0, paper (Fig. 2) | R0, here | P(R0 < 1) | final density ratio, here |
|:--|:--|:--|:--|:--|
| a no action | 9.6 (5.8–14.3) | 9.43 (5.60–14.35) | 0 | <0.001 |
| b body condition | ≈ 9.6 | 8.98 | 0 | <0.001 |
| c probiotic | ≈ 8.3 (text: increases) | 9.68 | 0 | <0.001 |
| d antifungal, 100% | ≈ 0.2 | 0.21 | 1 | 1.00 |
| e antifungal, 80% | ≈ 2.4; P < 2.5% | 2.26 (1.33–3.44) | 0.003 | 0.69 (paper 1.0) |
| f thinning 50% | ≈ 4.8 | 4.72 (2.80–7.17) | 0 | 0.006 (paper 0.083) |
| g thinning 90% | ≈ 0.93; P ≈ 50% | 0.94 (0.56–1.44) | 0.59 | 0.10 |
| h removal 90%/step | ≈ 0.2–0.3 | 0.34 | 1 | 0 |

The paper's ranking is reproduced:

- only the full-coverage antifungal (d) and removal (h) give `R₀ < 1` reliably;
- 90% thinning (g) is a coin toss;
- actions that slow Bsal growth (c–e) lengthen the distances infected adults travel (Fig. 3:
  37 m to 53 m before cohort survival falls below 1%);
- removal keeps them under 10 m.

## The simplified model and Fig. 4

Collapsing the load structure gives a discrete-time SI model (eq. S11):
`S(t+1) = S s₀(1 − φ)` and `I(t+1) = S s₀ φ + I s̄`, with `φ = 1 − exp(−βI)`. Its
`R₀ = β(1−m_β) S₀(1−m_{S₀}) s₀(1−m_s)/(1 − s̄(1−m_s))` (eq. S12) is
`basic_reproduction_number(p::CanessaSIParameters; mβ, mS0, msurv)`. [`required_reduction`](@ref)
solves eq. S12 for `R₀ = 1`. It reproduces Fig. 4a–c:

- 87% for transmission or density alone;
- 75% for survival alone;
- 42% for all three together.

## API

```@docs
CanessaDraw
canessa2018_path
load_canessa2018
CanessaAction
canessa_actions
canessa_action
CanessaParameters
canessa_growth_coefficients
canessa_survival
canessa_transmission
canessa_initial_load
canessa_growth_mean
canessa_growth_variance
CanessaIPMProblem
canessa_density_ratio
canessa_cohort_survival
canessa_dispersal
canessa_ensemble
canessa_summary
CanessaSIParameters
CanessaSIProblem
required_reduction
final_size_transmission
final_density_ratio
```

The `basic_reproduction_number` methods for [`CanessaParameters`](@ref) and
[`CanessaSIParameters`](@ref) are documented with the other methods of
[`basic_reproduction_number`](@ref).
