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

## Bd loads in frogs (Briggs et al. 2010)

```@docs
BriggsParameters
BriggsLoadProblem
transmission_rate
load_growth_rate
host_density_threshold
```

### With seasons and life stages

```@docs
BriggsStage
BriggsSeasonalParameters
BriggsSeasonalProblem
persistence_probability
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
