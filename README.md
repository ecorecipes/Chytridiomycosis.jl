# Chytridiomycosis.jl

Models of amphibian chytridiomycosis (*Batrachochytrium dendrobatidis*, Bd, and
*B. salamandrivorans*, Bsal) in Julia, built on SciML. Work in progress.

## What is here

Reference implementations of published models (`src/models/`), checked against the papers:

| Model | Source | Entry point |
|---|---|---|
| Bsal S–L–I in fire salamanders | Schmidt et al. 2017, doi:10.1098/rsos.160801 | `SchmidtSLIProblem` |
| Bd load model in frogs | Briggs et al. 2010, doi:10.1073/pnas.0912886107 | `BriggsLoadProblem`, `BriggsSeasonalProblem` |
| Bd in common toads at a pond | Mitchell et al. 2008, doi:10.1098/rspb.2007.1356 | `MitchellPondProblem` |

New models (`src/novel/`). **These are not published results**, and most of their added
parameter values are illustrative:

| Model | Entry point |
|---|---|
| Stochastic load birth–death–immigration | `LoadBDI` |
| Load tipping with saturating clearance, host differences and recruitment | `LoadTippingProblem`, `LoadTippingRecruitmentProblem` |
| Microbiome (probiotic) treatment on the three host models | `MicrobiomeTreatment`, `ProbioticSLIProblem`, `ProbioticPondProblem`, `ProbioticLoadProblem` |
| Microbiome treatment on Bsal with larvae | `ProbioticLarvalProblem` |
| Bacterial abundance on each frog, with repeated baths | `ProbioticAbundanceProblem` |
| One nested model with each piece of biology as a switch | `ChytridBiology`, `NestedChytridProblem` |
| Stochastic version of the nested model (whole animals) | `StochasticNestedProblem` |
| Stochastic toad pond (whole animals and zoospores) | `solve(::ProbioticPondProblem; rng)` |
| Stochastic frog loads (whole sporangia and zoospores) | `StochasticFrogProblem` |
| Bsal with an infectable larval stage | `BsalLarvalStage`, `BsalLarvalProblem` |
| Toads with infection after metamorphosis | `PostMetamorphicInfection`, `ToadStageProblem` |
| Frog load equation fitted to field swabs (Knapp et al. 2022), LogDensityProblems interface | `load_knapp2022`, `KnappLoadModel`, `maximum_likelihood`, `profile_likelihood` |

## Vignettes

[Does a protective skin microbiome save the population?](vignettes/01_building_realism/01_building_realism.md)
([PDF](vignettes/01_building_realism/01_building_realism.pdf),
[HTML](vignettes/01_building_realism/01_building_realism.html),
[source](vignettes/01_building_realism/01_building_realism.qmd))

For three systems (Bsal in fire salamanders, Bd in common toads, Bd in frogs) it starts from
the simplest model, adds the biology of the host and pathogen one piece at a time, then the
biology of the treatment, and records after every step how the population fares with and
without treatment. It ends by measuring how far each step moved the answer and, with a
factorial design, whether the effect of one piece depends on the others.

[One model, many switches](vignettes/02_nested_model/02_nested_model.md)
([PDF](vignettes/02_nested_model/02_nested_model.pdf),
[HTML](vignettes/02_nested_model/02_nested_model.html),
[source](vignettes/02_nested_model/02_nested_model.qmd))

Uses a single nested model in which each piece of biology is a switch, runs all 48
combinations at three host densities, and measures the effect of each piece like for like.
It ends with a stochastic version, which overturns one of the deterministic results.

[Fitting the frog load model to a field trial](vignettes/03_knapp_fit/03_knapp_fit.md)
([PDF](vignettes/03_knapp_fit/03_knapp_fit.pdf),
[HTML](vignettes/03_knapp_fit/03_knapp_fit.html),
[source](vignettes/03_knapp_fit/03_knapp_fit.qmd))

Fits the frog load equation to field swabs from frogs treated with itraconazole and then
*Janthinobacterium lividum* (Knapp et al. 2022), building up from no treatment effect to four
explanations. Something held treated frogs' loads down for weeks. A waning drug and the applied
bacteria fit about equally well, but resident *J. lividum* at the measured abundances did not
protect.

Each vignette writes out the equations of every model it runs, with a symbol table and a
diagram of the model.

Vignettes are Quarto documents with executable Julia code, rendered to HTML, PDF and
Markdown: run `quarto render` in `vignettes/`. Drawing the diagrams needs
[Graphviz](https://graphviz.org) (the `dot` program) on the `PATH`, and the PDFs need
`rsvg-convert` (librsvg) to convert them. Render one vignette and one format at a time, since
Quarto's Julia engine shares one notebook server.

## Parameters

Published values are read from a TOML compendium (`data/compendium/`, one file per source
paper) with units, provenance and a `verified` flag saying whether the value was checked
against the primary paper. Nothing published is hard-coded in the models.

```julia
using Chytridiomycosis
p = SchmidtParameters(:schmidt2017)
basic_reproduction_number(p)          # 24.8
sol = solve(SchmidtSLIProblem(p))
```

## Tests

```julia
using Pkg; Pkg.test()
```

Requires Julia 1.12.

## Documentation

[ecorecipes.github.io/Chytridiomycosis.jl](https://ecorecipes.github.io/Chytridiomycosis.jl)
(built by CI from `docs/`; the tutorials are the rendered vignettes, copied with
`julia scripts/sync_vignettes.jl`).

## Licence

MIT (see `LICENSE`). Data files in `data/` keep the licences given in their READMEs.
