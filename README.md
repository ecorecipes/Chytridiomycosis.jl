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
| Bsal with an infectable larval stage | `BsalLarvalStage`, `BsalLarvalProblem` |
| Toads with infection after metamorphosis | `PostMetamorphicInfection`, `ToadStageProblem` |

## Vignette

[Does a protective skin microbiome save the population?](vignettes/01_building_realism/01_building_realism.md)
builds the treatment model one piece of biology at a time, for Bsal in fire salamanders and Bd
in common toads, and shows how each added piece of realism changes the answer. Render it with
`quarto render` in `vignettes/` (Quarto with the Julia engine).

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
