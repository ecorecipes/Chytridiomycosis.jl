# Two-trait IPM: Bd and bacteria

```@meta
CurrentModule = Chytridiomycosis
```

!!! note "New model"
    This model is new work in this package, not a published result. The Bd side and the
    zoospore pool are published models (Wilber et al. 2016; Wilber et al. 2017) with values from
    the compendium. The bacterial decline, its variability and the bacterial load just after a
    bath come from this package's fit to the Knapp et al. (2022) swabs (below). The protective
    effects are illustrative, and that fit found no support for an effect on Bd growth.

Each frog carries two continuous traits: its Bd load `x` (ln zoospore equivalents per swab,
as in Wilber et al. 2016) and the load `y` of an applied antifungal bacterium,
*Janthinobacterium lividum* (ln genomic equivalents per swab). The treated/untreated models
([Microbiome treatment](treatment.md)) give every treated frog the same protection, and the
abundance model gives each frog a bacterial abundance but a single Bd load equation. Here the
whole joint distribution of the two loads is projected, so frogs with many bacteria and frogs
with few, and their Bd loads, are followed together.

## The model

Frogs are uninfected or infected, and carry the bacteria or not. A frog loses its bacteria
when `y` falls below the detection floor (1 genomic equivalent per swab, the lower edge of the
`y` mesh). Every 3 days, with `h(y)` the protection given by the bacteria,

```math
h(y) = \frac{1}{1 + e^{-(y - y_{50})/w}},
```

an infected frog survives with `logit s = b00 + b10 x + κs h(y)` (tolerance), clears Bd with
`logit l = b02 + b12 x + b22 T + κl h(y)`, and otherwise moves to

```math
x' \sim \mathcal{N}\!\left(b_{01} + b_{11} x + b_{21} T - κ_g h(y),\; v_{01} e^{2 c_{01} x}\right),
\qquad
y' \mid x' \sim \mathcal{N}\!\left(y - δΔt - ρ x + r\frac{σ_y}{σ_x}(x' - \mathrm{E}x'),\; σ_y^2 (1 - r^2)\right).
```

The bacteria act on Bd through `h(y)` (survival, growth, clearance and, optionally, the
infection probability), and Bd acts on the bacteria through `ρ` and the correlation `r`.
Frogs without bacteria acquire them from the pond and from other frogs, and a probiotic bath
gives a fraction of all frogs a fresh load. Bd spreads through the zoospore pool of Wilber et
al. (2017).

The protective threshold has empirical support: Becker et al. (2009) report that "a threshold
violacein concentration of about 18 μM on a salamander's skin prevents mortality", and cite a
violacein MIC against Bd of 1.82 μM (Brucker et al. 2008). They did not measure Bd load, and
*J. lividum* load predicted violacein only weakly (r = 0.33), so where the threshold lies on
the load axis is not known: `y50` is illustrative.

With every `κ`, `ρ` and `r` zero the two loads are independent and the Bd side is exactly the
Wilber et al. (2016) IPM (tested). Because the Bd axis uses the midpoint rule, as Wilber et al.
did, some frogs are lost where the load distribution is narrow; the bacterial axis uses exact
bin integrals.

Tolerance and resistance pull `R0` in opposite directions: bacteria that let infected frogs
survive longer at a given load (`κ_survival`) raise `R0`, because those frogs shed for longer,
while bacteria that slow Bd growth (`κ_growth`) lower it.

## Fit to the Knapp et al. (2022) cohort

The frogs of the Dusy Basin trial were group-marked, so a frog's two loads cannot be followed
from visit to visit. [`KnappJointModel`](@ref) fits the model to the joint distribution of Bd
and *J. lividum* loads in the treated group on days 12, 37 and 56, starting from its day-0
swabs and a bath of every frog (a **new fit**, not a published analysis;
`scripts/bacteria_ipm_knapp_fit.jl`). Infection is a constant probability per step, the
temperature is fixed at 15 °C, and ITS copies are converted to zoospore equivalents with a
fitted copy number within the range of Longo et al. (2013).

Without coupling (`coupled = false`) the fit gives a *J. lividum* half-life of 4.07 days, a
post-bath load of ln GE 6.5 (SD 2.1), `σy = 0.29` per step, an infection probability of 0.50
per step and 64 ITS copies per zoospore equivalent; these are the defaults of
[`BacteriaIPMParameters`](@ref). The half-life and `σy` trade off along a nearly flat ridge
(4.68 days and `σy = 0.15` fit almost as well, log-likelihood −540.40 against −540.38): with
group-marked frogs the cross-sections do not separate a slower, quieter decline from a faster,
noisier one. Letting the bacteria act on Bd growth finds no protective
effect: with `κ_growth` fixed anywhere between −2 and 4 (and the threshold re-fitted) the
log-likelihood is within 0.6 of the uncoupled fit. The unconstrained maximum (13 log-likelihood
units better) puts an implausibly large effect in the opposite direction (`κ_growth ≈ −107`
above ln GE 10.4), carried by a few swabs. Jani & Briggs (2014) likewise found
*Janthinobacterium* positively associated with Bd in an infection experiment, and not
correlated with it in the field.

## API

```@docs
IPMMesh2D
BacteriaIPMParameters
BacteriaIPMProblem
BacteriaIPMSolution
bacteria_ipm_summary
joint_load_distribution
basic_reproduction_number(::BacteriaIPMProblem)
KnappJointModel
knapp_joint_params
knapp_joint_start
pointwise_loglikelihood(::KnappJointModel, ::AbstractVector)
maximum_likelihood(::KnappJointModel)
```

## References

- Becker MH, Brucker RM, Schwantes CR, Harris RN, Minbiole KPC (2009) The bacterially produced
  metabolite violacein is associated with survival of amphibians infected with a lethal
  fungus. *Appl Environ Microbiol* 75:6635–6638.
  [doi:10.1128/AEM.01294-09](https://doi.org/10.1128/AEM.01294-09)
- Brucker RM et al. (2008) *J Chem Ecol* 34:1422–1429 (the violacein MIC, as cited by Becker
  et al. 2009).
- Jani AJ, Briggs CJ (2014) The pathogen *Batrachochytrium dendrobatidis* disturbs the frog
  skin microbiome during a natural epidemic and experimental infection. *PNAS*
  111:E5049–E5058. [doi:10.1073/pnas.1412752111](https://doi.org/10.1073/pnas.1412752111)
- Knapp RA, Joseph MB, Smith TC, et al. (2022) Effectiveness of antifungal treatments during
  chytridiomycosis epizootics in populations of an endangered frog. *PeerJ* 10:e12712.
  [doi:10.7717/peerj.12712](https://doi.org/10.7717/peerj.12712)
- Wilber MQ, Langwig KE, Kilpatrick AM, McCallum HI, Briggs CJ (2016) Integral projection
  models for host–parasite systems with an application to amphibian chytrid fungus. *Methods
  Ecol Evol* 7:1182–1194. [doi:10.1111/2041-210X.12561](https://doi.org/10.1111/2041-210X.12561)
- Wilber MQ, Knapp RA, Toothman M, Briggs CJ (2017) Resistance, tolerance and environmental
  transmission dynamics determine host extinction risk in a load-dependent amphibian disease.
  *Ecol Lett* 20:1169–1181. [doi:10.1111/ele.12814](https://doi.org/10.1111/ele.12814)
