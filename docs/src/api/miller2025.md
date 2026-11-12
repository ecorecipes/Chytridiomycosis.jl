# Refugia CTMC (Miller et al. 2025)

```@meta
CurrentModule = Chytridiomycosis
```

This is a reference implementation of the continuous-time Markov chain of Miller et al.
(2025, arXiv:2503.06846 v3, not peer reviewed). It models Bd in outdoor mesocosms of the
green and golden bell frog *Litoria aurea*, given either shaded or sunlight-heated
("unshaded") refugia. It was checked against the paper and the authors' R code (Zenodo
[10.5281/zenodo.15179344](https://doi.org/10.5281/zenodo.15179344)). Time is in weeks, and
states are frog counts per mesocosm.

## The model

There are eight states:

- `S_U`, `I_U1`, `I_U2`, `I_U3`: frogs with no prior infection;
- `S_V`, `I_V1`, `I_V2`, `I_V3`: frogs that have recovered from an infection
  ("vaccinated").

In environment `h`, a naive susceptible frog is infected at rate `β_h F/(N − 1)`, where
`F = Σ m_kj I_kj` weights each infected state by its relative infectiousness. A vaccinated
susceptible frog is infected at rate `α β_h F/(N − 1)`.

- Infected frogs move `I_k1 → I_k2 → I_k3` at rates `γ1 = 1/2.5` and `γ2 = 1/4.5` per week.
- They recover to `S_V` at rate `ω`.
- Frogs leave every state at the loss rate `μ`.

Each mesocosm starts with `(S_U, I_U1, S_V, I_V1) = (5, 5, 5, 5)` and is surveyed in weeks
0, 1, 2, 4, 6, 8, 10 and 15. The simulator follows the authors' GillespieSSA code exactly:

- the 0/0 propensity at `N = 1` is set to zero;
- each survey reads the state at that week.

```@docs
RefugiaParameters
transmissibility
RefugiaCTMCProblem
RefugiaODEProblem
refugia_initial_state
refugia_rates
refugia_stoichiometry
refugia_observation_weeks
```

## Data

The authors' state counts (CC BY 4.0) are bundled in `data/miller2025/`. `:exp` is the eight
mesocosms of Waddle et al. (2024), and `:syn` is the synthetic data of Supplement S2.

```@docs
RefugiaMesocosmData
miller2025_path
load_miller2025
loss_rate
```

## ABC rejection

The paper's ABC fixes `μ` from Eq. (3) and draws `(β_sh, β_un, α, ω)` from uniform priors.

1. For each draw it simulates one shaded and one unshaded mesocosm.
2. It scores the draw by the squared differences from every observed mesocosm in the same
   environment (Eq. 2).
3. It accepts the best 0.2% of 10⁶ draws.

A run takes a few seconds on four threads.

The paper's priors (Table 3) differ from those in the Zenodo code:

| | Table 3 | Zenodo code |
|:--|:--|:--|
| `β_sh`, `β_un` | U(0, 2) | U(0, 1) |
| `ω` | U(0, 2) | U(0, 10) |

[`RefugiaPrior`](@ref) provides both.

```@docs
RefugiaPrior
refugia_distance
refugia_abc
refugia_abc_summary
```

## Reproduction status

`scripts/miller2025_vs_R.jl` (with its R driver `scripts/miller2025_vs_R.R`) compares
10,000 runs of this simulator with 10,000 runs of the authors' R code at three parameter
sets. `scripts/miller2025_abc.jl` repeats the ABC and the identifiability checks. Its
output is stored in `data/miller2025/abc_reproduction.toml`, which is our reproduction and
not a published result.

### Against the authors' R code

We compared 10,000 runs per environment at each of three parameter sets: the Table 3 means,
the Table S2 values, and a high-transmission set.

- **Mean counts.** 313 state × week × environment comparisons. The z statistics have mean
  −0.04 and sd 0.97, and the largest |z| is 3.5. Rerunning that case with 160,000 Julia runs
  gives |z| ≈ 2.
- **Distributions.** 1 of 313 Kolmogorov–Smirnov tests has p < 0.01, against 3 expected by
  chance.
- **ABC summary statistic.** For the total statistic, KS p is 0.53, 0.03 and 0.78 for the
  three sets. At the Table 3 means, the mean total is 3372 in R and 3367 in Julia.

### ABC on the experimental data (paper priors; 20 runs of 10⁶ draws)

The means over seeds below are posterior means; ± is the sd over seeds.

| | Table 3 | ours: mean over seeds ± sd |
|:--|:--|:--|
| `β_sh` | 0.675 | 0.694 ± 0.013 |
| `β_un` | 0.401 | 0.414 ± 0.009 |
| `α` | 0.072 | 0.069 ± 0.002 |
| `αβ_sh` | 0.018 | 0.0185 ± 0.0005 |
| `αβ_un` | 0.009 | 0.0093 ± 0.0002 |
| `ω` | 0.569 | 0.565 ± 0.010 |
| refugia effect (Eq. 5) | 0.406 | 0.404 ± 0.014 |
| vaccination effect, shaded / unshaded (Eq. 4) | 0.973 / 0.978 | 0.9733 / 0.9776 |

Every Table 3 value lies within 3 sd of our mean. The code's priors do not give Table 3:
they give `β_sh` 0.33 and `ω` 1.79. So the published run used the priors stated in the
paper.

In our run, the posterior probability that `β_un < β_sh` is 0.74, and the posterior median
of `β_un/β_sh` is 0.50.

### Identifiability

- **α and β.** The posterior SD of `αβ_sh` is 0.08 of its prior SD. For `α` the ratio is 0.64
  and for `β_sh` it is 0.98, and `log α` and `log β_sh` have correlation −0.70.
  - In a simulation–estimation study with 20 data sets at the Table S2 values, the 95%
    intervals of `β` span 84–97% of the prior's. That of `αβ_sh` spans 7%.
  - All intervals cover the truth. The posterior mean of `αβ_sh` is still about 60% above the
    truth.
- **Loss and recovery.** With the loss rate free, μ ~ U(0, 1), μ is estimated above Eq. (3):
  posterior mean 0.06 per week, against 0.021. `ω` stays close to its prior. The posterior
  correlation of μ and ω is weak (0.15).
