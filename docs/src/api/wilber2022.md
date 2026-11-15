# Outbreak-state HMM (Wilber 2022)

```@meta
CurrentModule = Chytridiomycosis
```

Wilber, Knapp, Smith and Briggs (2022, *Journal of Animal Ecology* 91:2451–2464,
[10.1111/1365-2656.13823](https://doi.org/10.1111/1365-2656.13823)) asked whether host density
shapes what happens after *Bd* reaches a *Rana sierrae*/*R. muscosa* population, at three branch
points: whether *Bd* invades, how far the population declines, and how heavily frogs become
infected. They used 79 to 111 lakes of the Sierra Nevada surveyed from 2004 to 2019. The
implementation here follows the accepted manuscript and the authors' code and data (Zenodo
[10.5281/zenodo.7114634](https://doi.org/10.5281/zenodo.7114634)). Appendices S1–S4 were not
obtained. The deposit has no open licence, so it is not bundled; see `data/wilber2022/README.md`
and [`wilber2022_path`](@ref).

## Hidden Markov models

A small, model-agnostic layer: scaled forward and backward recursions, smoothed state
probabilities, Viterbi paths, covariate-dependent transition matrices, categorical emissions
and a LogDensityProblems wrapper for several independent sequences sharing parameters.
Emission "likelihoods" are a `T×K` matrix, so an observation model need not be a normalised
distribution (Wilber 2022's is not). The tests check every function against brute-force
enumeration of all paths.

```@docs
hmm_forward
hmm_loglikelihood
hmm_backward
hmm_posterior
hmm_viterbi
hmm_sample
CategoricalEmission
emission_likelihoods
CovariateTransitions
HMMLogDensity
```

## Hypothesis 1: failed and successful invasion

Each lake is in one of four states each June: U (*Bd*-free), LP (*Bd* at prevalence at most
`ρ`), HP (above `ρ`) or E (extirpated). Rows are the state in year `t`, columns in `t + 1`
(Fig. 2A):

```math
P = \begin{pmatrix}
1-a & a\omega & a(1-\omega) & 0\\
\varphi & \omega(1-\varphi) & (1-\omega)(1-\varphi) & 0\\
0 & 0 & 1-e & e\\
0 & 0 & 0 & 1
\end{pmatrix},
\qquad
\operatorname{probit}(\cdot) = \beta_0 + \beta_1 D_t + \beta_2 T_t + \beta_3 S_{t+1}
+ \beta_4 \mathrm{Tad}_t + \beta_5 T_t S_{t+1}
```

for `ω` and `φ` (eqn 1), with `D` the z-scored `log10(count + 1)` of adults and subadults,
`T` the z-scored mean June–September maximum air temperature, `S` the z-scored 1 April snow
water equivalent of the following spring, and `Tad` tadpole presence. `φ` large means failed
invasions; `ω` large means *Bd* stays at low prevalence. Where the count is missing, `φ` and
`ω` are averaged over `D ~ Normal(filled value, σ)`. Observed states are classified from the
swabs (`y` positives of `n`), and the emission weights account for sampling error
([`wilber2022_emission`](@ref)).

Differences between the printed text and the authors' code (the code is followed, with
keywords for two quirks):

- eqn 1 prints `log10(D)`; the code uses `log10(D + 1)`, with `D` the larger of the count and
  the number swabbed;
- every lake starts in U in 2004 and the 2004 observation is not used;
- the emission weights are not normalised over observations;
- the LP/HP split of swab counts is `n ÷ ⌊1/ρ⌋` (`cutoff_rule = :exact` uses `⌊nρ⌋`);
- the adult + subadult count is treated as missing when all subadult entries are missing
  (`subadult_rule = :sum` fixes this).

The posterior is sampled by adaptive random-walk Metropolis (the package has no MCMC
dependency; the model has 13–15 parameters and a cheap likelihood). Leave-one-out uses
[`psis_loo`](@ref) with one observation per lake, as the authors did.

```julia
fit = wilber2022_fit(cutoff = 0.25)                 # full model, abundance, cached
posterior_summary(fit)                              # medians and 95% intervals
wilber2022_effect(fit, :phi_temperature)            # Fig. 3: positive, interval above zero
looic(wilber2022_fit(cutoff = 0.25, model = :linear_temp)) - looic(fit)   # Table 1: 2.30
```

PSIS-LOO over 79 lakes has a Monte Carlo error of 0.2–0.8 LOOIC units per run (a lake or two
with Pareto `k` near 0.7), so `scripts/wilber2022_looic.jl` stores four-seed means in
`data/wilber2022/looic.toml`; the Table 1 targets use those.

Reproduced (validation targets `wilber2022_*`): all ten Fig. 3 effect sizes of the abundance
model at the three cutoffs (nine within 0.03, the long-tailed `φ`-tadpole within 0.15), the
LOOIC values behind Table 1 (absolute values within 0.8 of the authors', differences within
2.0), and the headline: abundance has a weak effect on `ω` in the direction opposite to the
density-threshold prediction and none on `φ`, while warm summers raise `φ` (less so before
severe winters) and tadpoles lower it. The density models of Table 1 are not reproducible
because the lake areas are not in the deposit.

```@docs
wilber2022_path
Wilber2022Data
wilber2022_data
wilber2022_emission
OutbreakHMMProblem
Wilber2022Parameters
wilber2022_parameters
outbreak_transition_matrix
outbreak_rates
outbreak_sequences
outbreak_state_probabilities
outbreak_viterbi
OutbreakHMMFit
looic
wilber2022_effect
wilber2022_fit
```

## Hypothesis 2: magnitude of decline

A beta-binomial regression of the decline from the maximum to the minimum count on the
z-scored `log10` of the maximum count (or density), for lakes where *Bd* prevalence peaked
between the two (eqn 2; 82 lakes counting subadults, 52 adults only).

```@docs
Wilber2022DeclineProblem
Wilber2022DeclineFit
decline_effect
predicted_decline
wilber2022_decline
```

## Hypothesis 3: density and infection intensity

The enzootic path analysis (Yosemite, 111 lake-year pairs from 40 lakes). The translocation
comparison (Fig. 6) needs data that are not in the deposit.

```@docs
wilber2022_path_analysis
```
