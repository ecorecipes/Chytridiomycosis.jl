# Load models

```@meta
CurrentModule = Chytridiomycosis
```

!!! note "New models"
    These models are new work in this package, not published results.

## Stochastic load on one host

A birth–death–immigration process for the number of sporangia on a host.

```@docs
LoadBDI
load_distribution
mean_load
prevalence
aggregation
infection_gain_rate
infection_loss_rate
enzootic_pool
enzootic_load_distribution
swab
```

## Load tipping

The Briggs load model with saturating clearance by the host, which gives a low-load state and
collective tipping when enough hosts share a pool.

```@docs
LoadTipping
LoadTippingProblem
LoadTippingRecruitmentProblem
pool_feedback
net_growth_rate
critical_immigration
tipping_regime
load_equilibria
pooled_load_equilibria
critical_host_number
mean_time_to_runaway
low_state_distribution
```

## Stochastic frog loads

```@docs
StochasticFrogProblem
```
