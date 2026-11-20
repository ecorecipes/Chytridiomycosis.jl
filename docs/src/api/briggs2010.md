# Bd loads in frogs (Briggs et al. 2010)

```@meta
CurrentModule = Chytridiomycosis
```

This is a reference implementation of the Bd load model of Briggs et al. (2010), *PNAS*
107:9695–9700 ([doi:10.1073/pnas.0912886107](https://doi.org/10.1073/pnas.0912886107)): the
load on each frog grows by reinfection from its own sporangia and from a shared zoospore pool,
and frogs die when their load passes a threshold. It is checked against the paper.
Parameters are read from the [compendium](compendium.md): `BriggsParameters(:briggs2010)`.

```@docs
BriggsParameters
BriggsLoadProblem
transmission_rate
load_growth_rate(::BriggsParameters, ::Any)
host_density_threshold
```

## With seasons and life stages

```@docs
BriggsStage
BriggsSeasonalParameters
BriggsSeasonalProblem
persistence_probability
```

The six-stage *R. muscosa*/*R. sierrae* model of Fig. 4 (three tadpole and two subadult year
classes plus adults):

```@docs
BriggsLifeStages
BriggsStagedProblem
disease_free_stages
```
