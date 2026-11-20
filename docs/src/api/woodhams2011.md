# Adaptive immunity and treatments (Woodhams et al. 2011)

```@meta
CurrentModule = Chytridiomycosis
```

This is a reference implementation of the models in the Appendix of Woodhams et al. (2011),
*Frontiers in Zoology* 8:8
([doi:10.1186/1742-9994-8-8](https://doi.org/10.1186/1742-9994-8-8)). The equations and values
were read from the PMC full text (PMC3098159), including its Appendix "Mathematical Models"
and the extended legend of Figure 6. The models are the [Briggs et al.
(2010)](@ref BriggsLoadProblem) load model in three settings: one frog in a constant
zoospore pool, the same frog with a hypothetical mucosal antibody, and a population of frogs
sharing a dynamic pool. The paper calls the antibody model "purely speculative". Time is in
days and loads are in sporangia per frog.

## One frog in a constant pool

```
dS/dt = ν(A) γZ + ν(A) ηf S − σS,   ν(A) = p exp(−θA)
dA/dt = rS − cA
```

`γZ` is the paper's "γ*Z", the flux `(γ/V)Z` of zoospores encountering the frog, and `ηf` its
"η*f". With no antibody (`θ = 0`) the load has the closed form
`S(t) = S0 e^{at} + b(e^{at} − 1)/a`, with `b = pγZ` and `a = pηf − σ`.

With the Fig. 6B values the package reproduces the following:

- The untreated frog reaches `Smax = 10000` on day 32.2, "in 33 days".
- Removing the frog from the pool on day 20 delays its death to day 36.9.
- Lowering `ν` to 0.15 on day 20 makes `a` negative. The load then falls, but toward
  `−b/a = 1500` sporangia, not to zero. The figure's curve agrees with this, but the legend
  says the frog "eventually clear[s] the infection" (target
  `woodhams2011_fig6b_probiotic_clears`, discrepant).
- With both treatments together, the load decays to zero.

!!! warning "Figures 6C and 6D cannot be rebuilt"
    The paper does not print the antibody production rate `r` or the antibody decay rate
    `c`, and for Fig. 6D it does not print the zoospore death rate `μ` either. The keyed
    constructors therefore require these as keywords. The stable load of 916 sporangia in
    Fig. 6C depends only on `r/c`. With `r/c = 0.01` the equilibrium is 917.4. Our inferred
    values `r = 1e-4` and `c = 0.01` reproduce the figure's curves, but they are not
    published, so the targets stay `pending`.

!!! note "Typesetting slips in the paper"
    - As rendered in PMC, the closed-form solution lacks its division by `a`.
    - The population equation has `λ` where it should have `γ`, and `σS` where it should
      have `σS_i`.
    - The legend of panel (d) writes "n" for `ν`.

```@docs
ImmuneFrogParameters
infection_fraction
constant_pool_load
time_to_lethal_load
constant_pool_equilibria
FrogTreatment
ConstantPoolFrogProblem
```

## A population with a dynamic pool

```
dS_i/dt = ν(A_i) γ Z + ν(A_i) η f S_i − σ S_i
dA_i/dt = r S_i − c A_i
dZ/dt   = η(1 − f) Σ S_i − γ Z N − μ Z
```

Frog `i` dies when `S_i` reaches `Smax`. With `θ = 0` this is the Briggs et al. (2010) model,
and `ImmunePopulationParameters(::BriggsParameters)` converts between the two. In Fig. 6D the
number of frogs is cut from 100 to 10 on day 30 ([`DensityReduction`](@ref)).

```@docs
ImmunePopulationParameters
DensityReduction
ImmunePopulationProblem
```
