# Nested and stochastic models

```@meta
CurrentModule = Chytridiomycosis
```

!!! note "New models"
    These models are new work in this package, not published results.

One host–pathogen–microbiome model in which each piece of biology is a switch, anchored on
Schmidt et al. (2017), and its stochastic version. Both draw on one list of transitions, so
they cannot disagree.

```@docs
ChytridBiology
NestedChytridProblem
invasion_growth_rate
StochasticNestedProblem
```
