# Multistate survival (Hollanders 2023, Palomar 2023)

```@meta
CurrentModule = Chytridiomycosis
```

Does infection *status* or infection *intensity* predict mortality in enzootic populations?
Two published capture–recapture analyses answer it differently, and both are implemented
here on one multistate layer built on the [hidden Markov model functions](wilber2022.md).

- **Palomar, Fernández-Chacón and Bosch (2023, *Biodiversity and Conservation* 32:793–809,
  [10.1007/s10531-022-02525-3](https://doi.org/10.1007/s10531-022-02525-3))** fitted a
  multi-event model to 13 years of male *Bufo spinosus* capture histories and found lower
  annual survival in *Bd*-positive toads (0.754 against 0.883). The histories were not
  deposited, so the model is implemented and tested on simulated data, and its estimates are
  pending validation targets.
- **Hollanders, Grogan, Nock, McCallum and Newell (2023, *Ecological Applications*
  33:e2724, [10.1002/eap.2724](https://doi.org/10.1002/eap.2724))** fitted a continuous-time
  robust-design multistate model to *Mixophyes fleayi* at three streams. Infection status had
  no effect on mortality, and swab load did. Their data are on Dryad under CC0 and are
  bundled (`data/hollanders2023/`); their code is on Zenodo; their posterior draws come from
  their repository and are not bundled.

## The multistate layer

Every individual's capture history is one HMM sequence, conditioned on its first capture.
There are `K` alive states and an absorbing dead state. Events are coded `0` (not seen) and
`1…E` (seen, with whatever the field record says); `missing` marks an occasion without a
survey.

- **Discrete time.** Survival `S` comes first, then movement between states `Ψ`
  ([`survival_transition`](@ref)). Detection `p` comes next, then state assignment `δ`
  ([`multievent_emission`](@ref)), as in E-SURGE's S, T, B1 and B2 steps.
- **Continuous time.** The generator `Q` holds the hazards. Each interval of length `τ` has
  the transition matrix `exp(Qτ)` ([`ctmc_transition`](@ref)), so unequal intervals need no
  special treatment. For two infection states and death ([`infection_generator`](@ref)),
  the matrix exponential has a closed form.
- **Robust design.** The state is constant within a primary occasion, and its secondary
  surveys are independent draws ([`robust_design_likelihoods`](@ref)).

```@docs
ctmc_transition
infection_generator
survival_transition
multievent_emission
robust_design_likelihoods
quasistationary_distribution
MultiEventStructure
MultiEventProblem
multievent_parameters
simulate_capture_histories
MultiEventFit
qaic
```

The tests check the multi-event likelihood against brute-force summation over every latent
path. They check the closed-form exponential against `LinearAlgebra.exp`. They also check
that simulated histories return the parameters that generated them.

## Palomar et al. 2023: multi-event survival of *Bufo spinosus*

The model has states I (*Bd*+), U (*Bd*−) and dead. Each year a toad survives, then may
change state. A toad alive at the pond is recaptured with probability `p`. A recaptured
toad's status is recorded with probability `δ` (otherwise "unknown").

The main text's event coding contradicts its own matrix B2. The text calls event 1
"Bd-negative", but B2 sends I to event 1. Online Appendix A codes 1 as infected and 2 as
uninfected, which agrees with B2, so that coding is used here.

```@docs
palomar2023_structure
palomar2023_models
palomar2023_parameters
palomar2023_implied_prevalence
palomar2023_model_selection
palomar2023_swabbing_experiment
```

Two checks use only the printed numbers.

1. **Implied prevalence.** The mean survival and transition probabilities imply that about
   71% of living males are infected in the long run. The observed swab prevalence was 27.8%.
2. **State assignment.** Swabbing was a random subsample of captures (588 of 1369), so the
   probability of recording a toad's status cannot depend on its state. The published
   estimates are nonetheless `δ_U = 1.000` (at the boundary) and `δ_I = 0.119`.

Both checks fit one explanation: the fitted model treats captures of unknown status as
infected toads. Without the data this cannot be tested on the real histories.
[`palomar2023_swabbing_experiment`](@ref) tests it on simulated histories with the same
design. Survival is set equal in both states, and the years without swabbing (2006, 2007,
2017 and 2018) are included. That experiment is our analysis, not a published result.

In one simulated data set (seed 2023), the published structure gave `δ_U = 1.000`,
`δ_I = 0.112`, `ψ_U→I = 0.66` and `ψ_I→U = 0.25`, close to the published values. The true
values were `δ = 0.43` in both states, `ψ_U→I = 0.2` and `ψ_I→U = 0.6`. Survival was not
biased. Declaring the unswabbed years unassessed recovered the transitions.

## Hollanders et al. 2023: load-dependent mortality in continuous time

The states are alive *Bd*−, alive *Bd*+ and dead. All rates are hazards per six weeks. The
infected mortality hazard is `φ₂ = φ₁ exp(β_status + β_load x)`, where `x` is the frog's
log10 load (ITS1 copies per swab) at the start of the interval. The load is centred at
3.174 and divided by 2 × 0.911. Recapture is modelled per secondary survey. Frogs were not
swabbed on the last primary occasion, so an infected frog seen then is recorded as seen and
*Bd*−.

The paper's printed text and the authors' code differ in several places, and the code is
followed throughout. The differences are listed in `src/models/hollanders2023.jl`. Two
matter for comparisons.

- The printed slopes are posterior summaries *conditional on RJMCMC inclusion*.
- The annual hazard is the six-weekly hazard × 52/6.

```@docs
hollanders2023_path
Hollanders2023Data
load_hollanders2023
hollanders2023_summary_data
load_hazard_ratio
hollanders2023_draws_path
hollanders2023_posterior
CTMultistateProblem
hollanders2023_load_model
CTMultistateFit
hollanders2023_start
hollanders2023_summary
hollanders2023_fit
hollanders2023_refit_path
```

## What was reproduced

- **From the bundled data.** The sampling summary (686 frogs, 1609 primary detections, 248
  of them infected, and the load range) is reproduced. So are the counts of observed gains
  and clearances.
- **From the authors' posterior draws.** Every printed median and HDI that was checked is
  reproduced: the status and load effects, the annual hazard and survival, and the gain and
  clearance rates and probabilities. The hazard ratios printed in the Discussion are not:
  4.4 at 10⁵ and 11 at 10⁶ ITS1 copies per swab, against about 6.4 and 19 from the draws.
- **From the reduced refit** (maximum likelihood, 26 parameters, about 10 minutes on two
  threads; 593 s recorded in `data/hollanders2023/refit.toml`).
  - Matched within one posterior sd: the clearance hazard (2.62 against 2.623), the gain
    hazard (0.23 against 0.283), the intensity effect (1.29 against 1.953) and temperature on
    clearing (1.11 against 1.101).
  - Outside one posterior sd: the status effect (1.22 against −0.167), the annual hazard
    (0.67 against 0.783) and rainfall on gaining (1.93 against 1.141), although all three
    published values lie inside our Wald intervals.
  - The intensity effect, though within one sd, is not 'present' by the paper's rule: Wald
    P(β ≤ 0) = 0.12 against the paper's pd of 97.5%.
  - The two effects are confounded (sampling correlation −0.89). Fixing status at −0.167
    lowers the log-likelihood by only 0.11, and the profiled fit gives intensity 2.58 and
    annual hazard 0.81.

See `data/validation/targets.toml`.
