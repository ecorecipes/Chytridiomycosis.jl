# Chytridiomycosis.jl

```@meta
CurrentModule = Chytridiomycosis
```

```@docs
Chytridiomycosis
```

Models of amphibian chytridiomycosis, caused by the fungi *Batrachochytrium dendrobatidis*
(Bd) and *B. salamandrivorans* (Bsal), built on SciML. The package contains reference
implementations of published models, a parameter compendium that records where every
published value comes from, and new models of treating amphibians with protective skin
bacteria.

!!! warning "Published and new work"
    The models under *Published models* reproduce their papers. Everything else is new work
    in this package and **is not a published result**. Apart from values read from the
    compendium, the parameters of the new models are illustrative.

## Published models

| System | Source | Entry point |
|:-------|:-------|:------------|
| Bsal in fire salamanders | Schmidt et al. 2017, [doi:10.1098/rsos.160801](https://doi.org/10.1098/rsos.160801) | [`SchmidtSLIProblem`](@ref) |
| Bd loads in frogs | Briggs et al. 2010, [doi:10.1073/pnas.0912886107](https://doi.org/10.1073/pnas.0912886107) | [`BriggsLoadProblem`](@ref), [`BriggsSeasonalProblem`](@ref) |
| Bd in common toads at a pond | Mitchell et al. 2008, [doi:10.1098/rspb.2007.1356](https://doi.org/10.1098/rspb.2007.1356) | [`MitchellPondProblem`](@ref) |

## New models

| Model | Entry point |
|:------|:------------|
| Stochastic load on one host | [`LoadBDI`](@ref) |
| Load tipping, with host differences and recruitment | [`LoadTippingProblem`](@ref), [`LoadTippingRecruitmentProblem`](@ref) |
| Microbiome treatment on the three host models | [`MicrobiomeTreatment`](@ref), [`ProbioticSLIProblem`](@ref), [`ProbioticPondProblem`](@ref), [`ProbioticLoadProblem`](@ref) |
| Bacterial abundance on each frog, repeated baths | [`ProbioticAbundanceProblem`](@ref) |
| Bsal with an infectable larval stage | [`BsalLarvalProblem`](@ref), [`ProbioticLarvalProblem`](@ref) |
| Toads with infection after metamorphosis | [`ToadStageProblem`](@ref) |
| One nested model with each piece of biology as a switch | [`NestedChytridProblem`](@ref) |
| Stochastic versions | [`StochasticNestedProblem`](@ref), [`StochasticFrogProblem`](@ref), `solve(::ProbioticPondProblem; rng)` |
| Frog load equation fitted to field swabs (Knapp et al. 2022) | [`KnappLoadModel`](@ref), [`maximum_likelihood`](@ref), [`profile_likelihood`](@ref) |

## Getting started

```julia
using Chytridiomycosis
p = SchmidtParameters(:schmidt2017)      # values from the compendium
basic_reproduction_number(p)             # 24.8
sol = solve(SchmidtSLIProblem(p))        # an ODE solution
```

Every problem is a struct describing the model and its set-up, solved with `solve`. Most
return a named tuple of time series and summaries; see each problem's documentation.

## Tutorials

- [Building biological realism](tutorials/01_building_realism.md): for three systems, the
  treatment model is built from the simplest model one piece of biology at a time, ending
  with chance.
- [One model, many switches](tutorials/02_nested_model.md): a single nested model in which
  every combination of biology is compared like for like, and a stochastic version.
- [Fitting to a field trial](tutorials/03_knapp_fit.md): the frog load equation fitted to
  field swabs of treated, wild and control frogs, comparing four explanations for why treated
  frogs stayed low.

The tutorials are Quarto documents in `vignettes/`, rendered to HTML, PDF and Markdown.

## Licence

MIT. Data files in `data/` keep the licences given in their READMEs.
