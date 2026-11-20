# Overlapping habitats and dilution (Espira et al. 2022)

```@meta
CurrentModule = Chytridiomycosis
```

This is a reference implementation of the two-species environmental transmission model of
Espira et al. (2022), *Am Nat* 199:E43–E56
([doi:10.1086/717413](https://doi.org/10.1086/717413)). It was checked against the author
manuscript (PMC9136953). The online supplement could not be obtained, so the
next-generation matrix was derived again from eq. 1. The paper is analytical: it prints no
parameter values, so there is no compendium entry, and its results are checked as statements
on random parameter sets.

Each species has susceptible, infectious and recovered hosts and sheds pathogens into its own
habitat, of area `ε_i`. A fraction `θ_i` of habitat `i` is shared with the other species.
There are no births, deaths or pathogen growth.

```
dS_i/dt = −κ_i ρ_i π_i S_i (E_i + θ_i E_j)
dI_i/dt =  κ_i ρ_i π_i S_i (E_i + θ_i E_j) − γ_i I_i
dR_i/dt =  γ_i I_i
dE_i/dt =  α_i I_i/ε_i − E_i (κ_i ρ_i N_i/ε_i + θ_j κ_j ρ_j N_j/ε_i + μ)
```

The pick-up efficiency `η_i = a_i/(a_i + μ)`, with `a_i = κ_i ρ_i N_i/ε_i`, places each
species between density-dependent-like (`η → 0`) and frequency-dependent-like (`η → 1`)
transmission. Adding species 2 **dilutes** when the community reproduction number `R0,C` falls
below species 1's own `R0,1`.

!!! warning "Eq. 6 as printed"
    The printed closed form of `R0,C` (eq. 6, and eq. 7 in the same form) has no habitat
    areas, and `θ₁` where eq. 1 gives `θ₂` in the first denominator. It equals the spectral
    radius of the next-generation matrix of eq. 1 when `ε₁ = ε₂ = 1` (then `θ₁ = θ₂`).
    [`community_reproduction_number`](@ref) uses the next-generation matrix;
    `as_printed = true` evaluates eq. 6 as printed. The dilution criterion of eq. 8 agrees
    with the next-generation matrix exactly when `θ₁ε₁ = θ₂ε₂`, as the definition of `θ`
    implies.

Reproduced: `R0,C ≥ min(R0,1, R0,2)`; `R0,C = max(R0,1, R0,2)` without overlap (eq. 10);
the complete-overlap form (eq. 11); the dilution criteria of eqs 8 and 9; and the
attenuation of dilution and amplification as species 1's pick-up efficiency rises (Fig. 3)
and as the overlap shrinks (Fig. 4).

```@docs
EspiraParameters
pickup_efficiency
species_reproduction_numbers(::EspiraParameters)
next_generation_matrix(::EspiraParameters)
community_reproduction_number
dilution_threshold
dilution_effect
espira_parameters_from_efficiencies
OverlappingHabitatProblem
```
