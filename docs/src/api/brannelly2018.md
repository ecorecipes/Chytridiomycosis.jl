# Multistate CMR and POPAN (Brannelly 2018)

```@meta
CurrentModule = Chytridiomycosis
```

Brannelly, Chatfield, Sonn, Robak and Richards-Zawacki (2018, *BMC Ecology* 18:34,
[doi:10.1186/s12898-018-0189-5](https://doi.org/10.1186/s12898-018-0189-5)) marked northern
cricket frogs (*Acris crepitans*) for a year in a lowland subtropical forest in Louisiana and
swabbed them at every capture. They fitted **conditional Arnason–Schwarz** multistate models
(M-SURGE) with two states, *Bd*-negative and *Bd*-positive, and a **POPAN** model (Program
MARK) for the super-population size. Infection made frogs easier to catch but did not lower
monthly survival (0.931 uninfected, 0.922 infected, model-averaged).

This is a reimplementation of a published model. The capture histories are the paper's
Additional file 3 (CC0), bundled in `data/brannelly2018`. The refit:

- matches the capture counts of Table 1, except the number of individuals (315 animals in
  the file, 266 printed);
- matches every parameter count of Table 2 and the ranking of its three best models;
- comes within 0.4% of the printed deviances (488.87 against 487.15 for the best model);
- gives model-averaged survival of 0.928 and 0.923, against the printed 0.931 and 0.922;
- does not reproduce the super-population size (1394 against 1854). The first capture
  probability goes to its bound, so `N` is weakly identified.

## Open-population capture–recapture layer

These functions extend the [multistate layer](multistate.md).

[`ArnasonSchwarzProblem`](@ref) is the conditional Arnason–Schwarz model. It conditions on
the state at first capture and estimates no initial-state parameters. Survival is per unit
time, raised to the length of each interval. Its designs are the
[`MultiEventStructure`](@ref) ones.

[`PopanProblem`](@ref) is the Schwarz–Arnason (1996) Jolly–Seber formulation. It has entry
probabilities `pent` and a super-population size `N = u + exp(θ_N)`.

```@docs
ArnasonSchwarzProblem
arnason_schwarz_parameters
ArnasonSchwarzFit
PopanProblem
popan_parameters
PopanFit
akaike_weights
```

## Brannelly et al. 2018

### Occasions and intervals

The paper combined each month's two surveys into one occasion. This gives 10 occasions,
which is what Table 2's parameter counts require. January, May and July were not surveyed,
so the intervals are 1, 2, 1, 1, 2, 2, 1, 1, 1 months.

Two choices were not stated in the paper:

- **Late September and early October.** 28 Sep and 1 Oct 2011 can be taken as September's
  survey pair (the default) or split by calendar month.
- **Interval lengths.** Survival can be raised to each interval's length in months (the
  default) or every interval can count as one step.

The defaults come closest to the printed deviances (see [`load_brannelly2018`](@ref)).

### Table 2 models

In this layer's terms, "f * to" (one parameter for each kind of transition) is `:state`, and
"+ t" adds an additive time effect, `:state_plus_time`.

```@docs
brannelly2018_path
Brannelly2018Data
load_brannelly2018
brannelly2018_summary
brannelly2018_structure
brannelly2018_models
brannelly2018_popan_models
brannelly2018_parameters
brannelly2018_fit
brannelly2018_model_selection
brannelly2018_survival
brannelly2018_popan
```
