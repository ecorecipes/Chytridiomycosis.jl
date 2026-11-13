# Published models

```@meta
CurrentModule = Chytridiomycosis
```

Reference implementations of published models, checked against the papers. Parameters are
read from the [compendium](compendium.md), for example `SchmidtParameters(:schmidt2017)`.

## Bsal in fire salamanders (Schmidt et al. 2017)

```@docs
SchmidtParameters
SchmidtSLIProblem
invasion_threshold
removal_fraction
basic_reproduction_number
endemic_equilibrium
hopf_threshold
wave_speed
scattered_colony_area
population_crash
```

### Spatial spread and long-distance jumps

The reaction–diffusion version (§2.2) is discretised by finite volumes and lowered to an
`ODEProblem`; its front speed converges to [`wave_speed`](@ref). The scattered-colony model
behind eq. (2.5) is simulated as a jump process with colonies placed in a square habitat, so
that overlaps can be counted.

```@docs
SchmidtPDEProblem
front_position
front_speed
SchmidtColonyProblem
```

## Bd loads in frogs (Briggs et al. 2010)

```@docs
BriggsParameters
BriggsLoadProblem
transmission_rate
load_growth_rate(::BriggsParameters, ::Any)
host_density_threshold
```

### With seasons and life stages

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

## Bd in common toads at a pond (Mitchell et al. 2008)

```@docs
MitchellParameters
MitchellPondProblem
```

## Solving

Every problem is solved with `solve` (from CommonSolve.jl). Those that lower to one ODE can
also be converted directly:

```@docs
to_ode_problem
```
