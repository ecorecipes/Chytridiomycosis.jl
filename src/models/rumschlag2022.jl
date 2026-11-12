# Rumschlag, Roth, McMahon, Rohr & Civitello (2022) J Anim Ecol 91:170-181,
# doi:10.1111/1365-2656.13612 (preprint doi:10.1101/2021.03.17.435818). Bd zoospores in four
# waters at 21 °C (artificial spring water and three Florida ponds), with Osteopilus
# septentrionalis tadpoles as hosts. Equations checked against the published PDF and the
# preprint, which agree. Wiki: lab/wiki/models/bd-rumschlag2022-zoospore-persistence.md
#
#   dZL/dt = −m ZL           (1)   live zoospores, per capita death rate m
#   dZD/dt =  m ZL − d ZD    (2)   dead zoospores, per capita decomposition rate d
#   dS/dt  = −β S ZL         (3)   susceptible tadpoles
#   dI/dt  =  β S ZL         (4)   infected tadpoles
#
# Units: time in days (the paper's lifespans are 24 h/m and 24 h/d), zoospores per µl (the
# data and Figs 3-4). The paper gives no units for β. Its printed value (7.1e-6) is the
# refitted β·ZL(0) divided by 1e5, so the paper's eq. 10 appears to use a fixed ZL(0) = 1e5
# for the infection assay. The compendium β therefore pairs with ZL0 = 1e5 in
# infection_probability. A refitted β (src/fitting/rumschlag2022.jl) is in
# µl zoospore⁻¹ day⁻¹ and pairs with the fitted ZL(0) in zoospores µl⁻¹.
#
# The paper assumes no host mortality and that infection removes a negligible number of
# zoospores, so eq. (1) has no −βSZL term. Eq. (6) as printed has its bracket reversed; the
# pool code (zoospore_pool.jl) uses the correct solution.

const RUMSCHLAG_WATERS = (:ASW, :TroutCreek, :Flatwoods, :GreenSwamp)

"""
    RumschlagParameters(; m, d, β)
    RumschlagParameters(:rumschlag2022; water = :ASW)

Parameters of the Rumschlag et al. (2022) zoospore persistence model for one water source:

- `m`: per capita death rate of live zoospores (day⁻¹);
- `d`: per capita decomposition rate of dead zoospores (day⁻¹);
- `β`: per capita transmission rate (day⁻¹ per unit of the zoospore concentration used for
  `ZL0`).

The keyed constructor reads the paper's selected model (four water-specific `m`, one `d`,
one `β`) from the compendium. `water` is one of `:ASW` (artificial spring water),
`:TroutCreek`, `:Flatwoods` or `:GreenSwamp`.

- **Rates.** The paper reports lifespans in hours. The constructor converts them with the
  units layer: `m = 1/lifespan` and `d = 1/persistence`, per day.
- **β.** `β` is the printed value, 7.1 × 10⁻⁶, which has no stated units. A refit to the
  deposited data shows that it equals `β ZL0 / 10⁵`, so it reproduces the paper's prevalence
  (Fig. 4) only with `ZL0 = 1e5` in [`infection_probability`](@ref). That convention is
  inferred, not printed.
- **Refitted β.** A β refitted by [`maximum_likelihood`](@ref) is in µl zoospore⁻¹ day⁻¹. It
  goes with the fitted initial concentration (about 143 zoospores µl⁻¹).
"""
Base.@kwdef struct RumschlagParameters
    m::Float64
    d::Float64
    β::Float64
end

const RUMSCHLAG_SETS = Dict(:ASW => "ASW", :TroutCreek => "TroutCreek",
    :Flatwoods => "Flatwoods", :GreenSwamp => "GreenSwamp")

function RumschlagParameters(key::Symbol; water::Symbol = :ASW)
    haskey(RUMSCHLAG_SETS, water) ||
        throw(ArgumentError("water must be one of $(RUMSCHLAG_WATERS); got :$water"))
    set = RUMSCHLAG_SETS[water]
    recs = parameter_records(key)
    function rec(sym)   # the water's own record, else the base record
        i = findfirst(r -> r.symbol === sym && r.set == set, recs)
        i === nothing && (i = findfirst(r -> r.symbol === sym && isempty(r.set), recs))
        i === nothing && throw(ArgumentError("$key has no $sym for water $water"))
        return recs[i]
    end
    return RumschlagParameters(; m = 1 / in_time_unit(rec(:lifespan), :day),
        d = 1 / in_time_unit(rec(:persistence), :day), β = rec(:β).value)
end

time_unit(::Type{<:RumschlagParameters}) = :day

"""
    ZoosporePool(p::RumschlagParameters)

The live/dead [`ZoosporePool`](@ref) of one water source (rates per day).
"""
ZoosporePool(p::RumschlagParameters) = ZoosporePool(p.m, p.d)

"""
    zoospore_lifespan(p::RumschlagParameters) -> Float64

Expected lifespan of a live zoospore, `1/m`, in **hours** (the unit of the paper's results:
"dividing 24 hr by the mortality rate").
"""
zoospore_lifespan(p::RumschlagParameters) = 24 / p.m

"""
    persistence_time(p::RumschlagParameters) -> Float64

Expected time from the death of a zoospore to its complete decomposition, `1/d`, in
**hours**.
"""
persistence_time(p::RumschlagParameters) = 24 / p.d

function rumschlag_rhs!(du, u, p::RumschlagParameters, t)
    ZL, ZD, S, I = u
    du[1] = -p.m * ZL
    du[2] = p.m * ZL - p.d * ZD
    du[3] = -p.β * S * ZL
    du[4] = p.β * S * ZL
    return nothing
end

"""
    ZoosporePersistenceProblem(p::RumschlagParameters; ZL0, ZD0 = 0, S0 = 1, I0 = 0,
                               tspan = (0.0, 2.0))

The four-state model of Rumschlag et al. (2022), eqs 1–4, in days:

```
dZL/dt = −m ZL
dZD/dt =  m ZL − d ZD
dS/dt  = −β S ZL
dI/dt  =  β S ZL
```

`ZL0` has no default. It is the initial live zoospore concentration, on the scale that `β`
is per (see [`RumschlagParameters`](@ref)). For the survival assay it is in zoospores µl⁻¹.
The paper fits it but prints no value; its Fig. 3a and the refit put it near 143 µl⁻¹. `solve(prob)` returns the `ODESolution` of the state
`[ZL, ZD, S, I]`. Use [`to_ode_problem`](@ref) for the `ODEProblem`. To follow the
infection assay, in which a tadpole is added at lag `t1` and removed one day later, start
from `ZL0 = ZL(0) exp(−m t1)` with `tspan = (0, 1)`, or use
[`infection_probability`](@ref).
"""
struct ZoosporePersistenceProblem
    p::RumschlagParameters
    u0::Vector{Float64}
    tspan::Tuple{Float64,Float64}
end
function ZoosporePersistenceProblem(p::RumschlagParameters; ZL0::Real, ZD0::Real = 0.0,
        S0::Real = 1.0, I0::Real = 0.0, tspan = (0.0, 2.0))
    all(>=(0), (ZL0, ZD0, S0, I0)) || throw(ArgumentError("initial states must be non-negative"))
    return ZoosporePersistenceProblem(p, Float64[ZL0, ZD0, S0, I0],
        (Float64(tspan[1]), Float64(tspan[2])))
end

to_ode_problem(prob::ZoosporePersistenceProblem) =
    ODEProblem(rumschlag_rhs!, prob.u0, prob.tspan, prob.p)

CommonSolve.solve(prob::ZoosporePersistenceProblem, alg = Tsit5(); kwargs...) =
    solve(to_ode_problem(prob), alg; kwargs...)

"""
    infection_probability(p::RumschlagParameters, ZL0, t1; τ = 1.0)

Probability that a tadpole is infected when it is added `t1` days after `ZL0` live
zoospores and removed `τ` days later (Rumschlag et al. 2022, eq. 10; `τ = 1` day in the
paper):

```
p = 1 − exp(β ZL0 exp(−m t1) (exp(−m τ) − 1)/m)
```

`ZL0` must be on the scale that `β` is per. With the compendium's β use `ZL0 = 1e5`, the
paper's inferred convention. For example, `infection_probability(RumschlagParameters(:rumschlag2022), 1e5, 0.0)`
gives the paper's Fig. 4a prevalence in artificial spring water at lag 0 (about 0.29). With
a refitted β, use the fitted `ZL0` (zoospores µl⁻¹).
"""
infection_probability(p::RumschlagParameters, ZL0, t1; τ = 1.0) =
    -expm1(-p.β * ZL0 * exp(-p.m * t1) * _expm1_ratio(p.m, τ))
