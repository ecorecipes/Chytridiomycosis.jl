# Overlapping habitats and dilution (Espira et al. 2022)

```@meta
CurrentModule = Chytridiomycosis
```

This is a reference implementation of the two-species environmental transmission model of
Espira et al. (2022), *Am Nat* 199:E43–E56
([doi:10.1086/717413](https://doi.org/10.1086/717413)). It was checked against the published
version, the author manuscript (PMC9136953; same equations) and the online supplement
(eqs S1–S18, Figs S1–S2). The paper is analytical: its figures use no printed parameter
values, and its results are checked as statements on random parameter sets. The only numbers
are the supplement's worked example (Table S1), which is in the compendium.

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

!!! warning "Corrections to the paper"
    - **Eq. 6 is a misprint.** The printed closed form of `R0,C` (eq. 6, and eq. 7 in the
      same form, in both the published version and the manuscript) has `μ` without `ε_i`, and
      `θ₁` where eq. 1 gives `θ₂` in the first denominator. The supplement's own S12, S15 and
      S16 have `μ ε_i` and `θ_j`, which is the next-generation matrix of eq. 1.
      [`community_reproduction_number`](@ref) computes S15; `as_printed = true` evaluates
      eq. 6 as printed. The two agree only when `ε₁ = ε₂ = 1`.
    - **Eq. 8 needs consistent overlaps.** The dilution criterion of eq. 8 (and S17) agrees
      with `R0,C < R0,1` exactly when `θ₁ε₁ = θ₂ε₂`, as the definition of `θ` implies.
    - **Figs S1 and S2** (habitat 1 inside habitat 2, and the reverse) are printed with the
      same heat-map raster as Fig. 4. On the figure axes `R0,C/R0,1` depends on `θ₁` and `θ₂`
      separately, so the three cannot coincide: S1 differs from Fig. 4 by up to 0.52 and S2
      by up to 0.17. With habitat 1 inside habitat 2, dilution grows stronger, not weaker, as
      the overlap shrinks.
    - **Worked example.** Its text uses `κ₂ρ₂ = 0.0002` (Table S1: 0.0001) and computes
      `4 × 2 = 4` for the decay term.

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
