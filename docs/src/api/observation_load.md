# Swab loads: units, likelihoods and observation operators

```@meta
CurrentModule = Chytridiomycosis
```

A field swab passes through a chain before it becomes a number: the swab picks up part of
the infection, the DNA is extracted and part of the extract is run in a qPCR, the cycle
threshold is read off a standard curve, and the result is scaled to the whole swab. This page
covers the end of that chain, which the package can state explicitly, and the link from model
states to what swabs report. It builds on the load types of [Units and loads](units.md).

| step | function | source |
|---|---|---|
| Ct → quantity in the reaction | [`QPCRStandardCurve`](@ref), [`standard_quantity`](@ref) | Longo et al. 2013, Table 1 |
| ZE of one strain ↔ ITS1 copies ↔ ZE of another | [`restandardize`](@ref), [`strain_copy_number`](@ref) | Longo et al. 2013, Fig. 4A |
| DNA mass ↔ copies, genomes | [`dsdna_mass`](@ref), [`genome_equivalents`](@ref) | Longo et al. 2013, Methods |
| qPCR count → live zoospores | [`shed_viability`](@ref) | Maguire et al. 2016, Table 2 (one lab study) |
| swabs → likelihood | [`LoadHurdleModel`](@ref), [`compare_load_families`](@ref) | standard two-part model; Schrock et al. 2025 for the comparison |
| groups of swabs → aggregation | [`TaylorPowerLawProblem`](@ref) | Schrock et al. 2025 |
| model state → swab distribution | [`swab_distribution`](@ref) | per model; see below |

**Conversions the sources do not give, and the package refuses:** ZE per swab to sporangia
per host (the Briggs et al. 2010 state), ZE per swab to zoospores on the whole host (needs a
swab recovery fraction nobody has measured for these systems), and qPCR DNA to live
zoospores outside the one lab study of Maguire et al. (2016). Each needs an explicit keyword;
the sporangia conversion also warns.

## The units map

ZE (zoospore equivalents) are a qPCR quantity scaled against standards made from counted
zoospores of one strain. The eight strains of Longo et al. (2013) carry 10 to 144 ITS1 copies
per zoospore, so the same swab reads up to 14.4-fold differently depending on which strain the
lab used for its standards:

```julia
strain_copy_number("JEL404")                          # 39.0 copies per zoospore
restandardize(ZEPerSwab(1000), "JEL404" => "CLFT023") # ZEPerSwab(270.8)
restandardize(ZEPerSwab(1000), "JEL404" => "MexMkt")  # ZEPerSwab(3900.0)

c = QPCRStandardCurve(:longo2013, "MexMkt")           # Ct = -3.31 log10(q) + 33.9
qpcr_efficiency(c)                                    # 1.005 (Table 1: 101%)
standard_quantity(c, 30.0)                            # 15.1 zoospores in the reaction
```

Chestnut et al. (2014) divided their water-sample qPCR quantities by 169 copies per zoospore,
above every strain of Longo et al.; `restandardize(x, 169 => "CLFT023")` re-expresses their
densities on a CLFT023 basis (1.17-fold higher).

What `restandardize` cannot fix: a ZE value assumes the infecting strain has the standard's
copy number. When it does not, ZE count ITS1 copies, not zoospores, and only a strain-specific
copy number for the field strain would turn them into zoospores.

```@docs
QPCRStandardCurve
qpcr_efficiency
standard_quantity
cycle_threshold
longo_strains
strain_copy_number
restandardize
dsdna_mass
genome_equivalents
shed_viability
```

## Load likelihoods

[`LoadHurdleModel`](@ref) is the shared swab likelihood: a fraction `π` of swabbed hosts is
infected, an infected host's load follows a lognormal, gamma, Weibull or negative binomial
distribution, and loads below a limit of detection are left-censored:

```math
P(\text{non-detect}) = (1-\pi) + \pi\,P(Y < L), \qquad f(y \ge L) = \pi\, f_D(y).
```

With `lod = 0` it is the classic hurdle (two-part) model. On `scale = :count` every family
gives the probability of the rounded-up load, so the negative binomial can be compared with
the continuous families by AIC; on the continuous scale only the continuous families are
comparable. The model implements LogDensityProblems (flat prior on the unconstrained
parameters), so it can be sampled with any LogDensityProblems sampler or fitted with
`maximum_likelihood(m)`, which also has a method for it (see [Fitting to field data](fitting.md));
`pointwise_loglikelihood(m, x)` feeds [`psis_loo`](@ref).

```julia
y = [s.bd for s in load_knapp2022() if s.group === :wild && !isnan(s.bd)]  # ITS1 copies per swab
f = maximum_likelihood(LoadHurdleModel(y; family = :lognormal))
f.prevalence, f.distribution, f.aic
compare_load_families(y)        # positive loads only, as Schrock et al. fitted each group
```

```@docs
LoadHurdleModel
hurdle_distribution
fit_positive_loads
compare_load_families
```

## Taylor's power law

Schrock et al. (2025; preprint 2024) related the log variance to the log mean of Bd loads
across groups of swabs (species × life stage × site × season × year), zeros included, and
found slopes of 1.90 (Sierra Nevada) to 2.06 (Brazil), steeper than the macroparasite 1.55.
[`TaylorPowerLawProblem`](@ref) follows their code: OLS of log10 variance on log10 mean,
groups with at least 10 swabs and 3 positives, R's `confint` intervals.

A lognormal load with the same log-scale spread in every group gives a slope of exactly 2
when every swab is positive. Their groups include the uninfected swabs, so with prevalence
`π`, `Var/mean² = (e^{σ²} − π)/π` ([`load_moments`](@ref)): groups whose prevalence rises
with their mean load pull the slope below 2.

Their group summaries are on Dryad (CC0) but are not bundled: the download was refused to the
package's automated client. [`load_schrock2024`](@ref) reads the file once it is placed in
`data/schrock2024/` (see the README there), and [`schrock_taylor_slopes`](@ref) reproduces
their Fig. 1 regressions; the validation targets for the published slopes are pending until
then.

```@docs
TaylorPowerLawProblem
load_moments
SchrockGroup
schrock2024_path
load_schrock2024
schrock_taylor_slopes
```

## Observation operators

[`swab_distribution`](@ref) turns a model state into a [`SwabLoadDistribution`](@ref): a
prevalence, a mixture over log10 ZE per swab, optional Gaussian noise `σ` on the log10
scale and an optional limit of detection. Its [`swab_logpdf`](@ref) is a density on the
natural scale, so it can be compared with the [`LoadHurdleModel`](@ref) likelihoods, and
[`load_moments`](@ref) gives the mean and variance that Taylor's power law relates.

| model | state | operator | what is published, what is ours |
|---|---|---|---|
| Wilber et al. 2016 ([`BdIPMSolution`](@ref)) | ln ZE per swab on a mesh | each bin a uniform component | x **is** the swab load (the kernel was fitted to swabs); with `σ = 0`, `lod = 0` nothing is added |
| Wilber et al. 2017 ([`HybridIPMSolution`](@ref)) | adults: ln ZE per swab on a mesh | as above, adults only | as above |
| Louca et al. 2014 | DOI bins, read as ZE per swab (ESM §5.1) | each bin a point | the reading is the paper's; `σ > 0` (needed for a density) is ours |
| Briggs et al. 2010 | sporangia per host | **refused** without `ze_per_sporangium`; warns with it | no published conversion: the result rests on that number |
| any mean-load model | mean ZE per swab | lognormal with stated `sdlog10` | the spread among hosts is ours |

```julia
sol = solve(BdIPMProblem(WilberParameters(:wilber2016; set = "deposited"); temperature = 16))
d = swab_distribution(sol)                    # the published state, as swabs
detection_probability(swab_distribution(sol; lod = 10.0))   # with a 10 ZE LOD (our choice)
load_moments(d)                               # a point on the Taylor plot

# Briggs: the sporangia → ZE step must be stated, and is warned about
swab_distribution(SporangiaPerHost.(loads); ze_per_sporangium = 1.0, σ = 0.3)
```

Field swabs from Knapp et al. (2022) are ITS1 copies, not ZE: compare them with a ZE
prediction only after `convert_load(x, ZEPerSwab; copies_per_zoospore = …)` with a stated
copy number.

```@docs
SwabLoadDistribution
swab_distribution
detection_probability
swab_logpdf
swab_loglikelihood
```

## Validation

The registry block `# --- observation load` of `data/validation/targets.toml` holds the
reproduced numbers (Longo et al.'s efficiencies, copy-number range and amplicon mass; the
slope of 2 for lognormal loads; Maguire et al.'s 80% viability) and two discrepancies: Longo
et al.'s text gives 6.0 × 10⁶ copies in 1.5 ng where Fig. 4B's highest bar reads 6.4 × 10⁶,
and Chestnut et al.'s "approximately four orders of magnitude" between their largest filter
and any previous observation is 2.81 orders against the Nevada filter they name.
