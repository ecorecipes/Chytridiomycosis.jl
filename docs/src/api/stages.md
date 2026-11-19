# Life stages

```@meta
CurrentModule = Chytridiomycosis
```

!!! note "Modified published models"
    These models add a life stage to a published model: an infectable larval stage to
    Schmidt et al. (2017), and infection after metamorphosis to Mitchell et al. (2008). They
    are new work in this package, not published results.

## Bsal with an infectable larval stage

The Schmidt et al. (2017) model with larvae that can carry infection through metamorphosis
(carry-over from Wang et al. 2021).

```@docs
BsalLarvalStage
BsalLarvalProblem
disease_free_larvae
larval_reproduction_number
```

## Toads with infection after metamorphosis

The Mitchell et al. (2008) model with surviving infected metamorphs and adult carriers.

```@docs
PostMetamorphicInfection
ToadStageProblem
```
