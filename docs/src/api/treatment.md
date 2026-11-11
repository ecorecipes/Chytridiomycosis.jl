# Microbiome treatment

```@meta
CurrentModule = Chytridiomycosis
```

!!! note "New models"
    These models are new work in this package, not published results. Apart from values read
    from the compendium, treatment parameters are illustrative.

A treated class of hosts carrying protective skin bacteria, put on each host–pathogen model.

## The treatment

```@docs
MicrobiomeTreatment
load_reducing_treatment
protection
protection_duration
reversion_rate
microbiome_reproduction_number
equilibrium_coverage
treatment_efficacy
mortality_reduction
```

## On the host models

```@docs
ProbioticSLIProblem
effective_reproduction_number
critical_coverage
ProbioticPondProblem
ProbioticLoadProblem
ProbioticAbundanceProblem
ProbioticLarvalProblem
```
