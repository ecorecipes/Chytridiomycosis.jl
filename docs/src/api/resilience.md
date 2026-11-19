# Resilience and thresholds

```@meta
CurrentModule = Chytridiomycosis
```

The resilience machinery is in
[DynamicalResilience.jl](https://ecorecipes.github.io/DynamicalResilience.jl/dev/) (steady states,
R₀ by the next-generation method, Floquet multipliers, following resilience along a parameter)
and [ResilienceCore.jl](https://ecorecipes.github.io/ResilienceCore.jl/dev/) (what a linearised
operator means and which metrics apply to it). This package only says where its models'
equilibria are and which of their terms are new infections. The tutorial
[Thresholds and resilience](../tutorials/06_resilience.md) uses them on Schmidt's Bsal model and
on the load-tipping model.

```@docs
equilibrium_state
ResilienceCore.linearise(::SchmidtSLIProblem)
DynamicalResilience.reproduction_number(::SchmidtSLIProblem)
to_ode_problem(::LoadTipping)
```

The individual-based seasonal models (for example [`BriggsStagedProblem`](@ref)) have no ODE to
linearise; their resilience is a question for simulation.
