# Woodhams, Bosch, Briggs et al. (2011) Mitigating amphibian disease: strategies to maintain
# wild populations and control chytridiomycosis. Frontiers in Zoology 8:8,
# doi:10.1186/1742-9994-8-8. Appendix "Mathematical Models" and the extended legend of Fig. 6.
# Equations read from the PMC rendering (PMC3098159) and its equation images; the publisher
# PDF could not be downloaded.
#
# The Briggs et al. (2010) load model (src/models/briggs2010.jl) in three settings. Time in
# days, loads in sporangia per frog.
#
# 1. One frog in a constant zoospore pool, no adaptive immunity:
#      dS/dt = ν(γ/V)Z + ηνf S − σS = b + aS,   b = ν(γ/V)Z,  a = ηνf − σ
#      S(t) = (aS(0)e^{at} + b e^{at} − b)/a
#    The solution as rendered in PMC lacks the division by a; the ODE requires it.
# 2. The same frog with a hypothetical mucosal antibody A ("purely speculative" in the paper):
#      dS/dt = ν(A)(γ/V)Z + ν(A)ηf S − σS,   dA/dt = rS − cA,   ν(A) = p exp(−θA)
# 3. A population of frogs i = 1..N sharing a dynamic pool Z:
#      dS_i/dt = ν_i(A_i)(γ/V)Z + ν_i(A_i)ηf S_i − σS_i
#      dA_i/dt = rS_i − cA_i
#      dZ/dt   = η(1−f)S_tot − (γ/V)Z N − μZ
#    Frog i dies when S_i reaches Smax. As printed, the S_i equation has λ for γ and σS for
#    σS_i (typesetting slips; the text defines only γ). With θ = 0 this is exactly the
#    BriggsLoadProblem.
#
# Unpublished choices:
# - The Fig. 6 legend gives the single-frog parameters only as the products "γ*Z" = (γ/V)Z
#   and "η*f" = ηf, so ConstantPoolFrogProblem is parametrised by those products.
# - r and c (antibody production and decay) are not printed for Fig. 6C or 6D, nor μ for 6D.
#   The keyed constructors therefore require them as keywords; no defaults are invented.
# - A dead frog's load and antibody are frozen at 0 and it stops shedding and encountering
#   zoospores (as in BriggsLoadProblem). A treatment that sets the load (the antifungal of
#   Fig. 6C) leaves A unchanged. A density reduction removes the highest-numbered live frogs
#   (all frogs are identical in the paper's deterministic runs, so the choice is immaterial
#   there).

"""
    ImmuneFrogParameters(; p, γZ, ηf, σ, Smax, θ = 0, r = 0, c = 0)
    ImmuneFrogParameters(:woodhams2011; set = "fig6b")
    ImmuneFrogParameters(:woodhams2011; set = "fig6c", r, c)

One frog in a constant zoospore pool (Woodhams et al. 2011, Appendix), per day:

- `p`: fraction of encounters that infect the frog with no antibodies (`ν` of the model
  without immunity);
- `γZ`: zoospores from the pool encountering the frog, `(γ/V)Z` (zoospores day⁻¹);
- `ηf`: zoospores per sporangium that immediately re-encounter the same frog (day⁻¹);
- `σ`: sporangium loss rate (day⁻¹); `Smax`: lethal load (sporangia);
- `θ`: fall of the infection fraction with antibody, `ν(A) = p exp(−θA)`;
- `r`, `c`: antibody production per sporangium and antibody decay (day⁻¹).

With `θ = 0` (or `r = 0`) there is no adaptive immunity and the load follows the closed form
[`constant_pool_load`](@ref).

The keyed constructor reads the compendium. Set `"fig6b"` is the frog of Fig. 6B, without
immunity. Set `"fig6c"` is the frog of Fig. 6C with immunity. **The paper does not print `r`
or `c`**, so they must be given; no published values exist.
"""
Base.@kwdef struct ImmuneFrogParameters
    p::Float64
    γZ::Float64
    ηf::Float64
    σ::Float64
    Smax::Float64
    θ::Float64 = 0.0
    r::Float64 = 0.0
    c::Float64 = 0.0
end

function ImmuneFrogParameters(key::Symbol; set::AbstractString,
        r::Union{Nothing,Real} = nothing, c::Union{Nothing,Real} = nothing)
    q = parameter_set(key; set)
    θ = hasproperty(q, :θ) ? q.θ : 0.0
    if θ != 0 && (r === nothing || c === nothing)
        throw(ArgumentError("$key set \"$set\" has adaptive immunity, but the paper prints " *
            "no antibody production r or decay c: pass both as keywords"))
    end
    return ImmuneFrogParameters(; p = hasproperty(q, :p) ? q.p : q.ν, γZ = q.γZ, ηf = q.ηf,
        σ = q.σ, Smax = q.Smax, θ, r = something(r, 0.0), c = something(c, 0.0))
end

time_unit(::Type{ImmuneFrogParameters}) = :day

_immune_with(p::ImmuneFrogParameters; kw...) = ImmuneFrogParameters(; (f => getfield(p, f)
    for f in fieldnames(ImmuneFrogParameters))..., kw...)

"""
    infection_fraction(p, A)

Fraction of zoospore encounters that infect a frog with antibody level `A`,
`ν(A) = p exp(−θA)` (Woodhams et al. 2011). Methods for [`ImmuneFrogParameters`](@ref) and
[`ImmunePopulationParameters`](@ref).
"""
infection_fraction(p, A) = p.p * exp(-p.θ * A)

_has_immunity(p) = p.θ != 0 && p.r != 0

"""
    constant_pool_load(p::ImmuneFrogParameters, t; S0 = 0)

Load on a frog in a constant pool without adaptive immunity (Woodhams et al. 2011,
Appendix): `S(t) = S0 e^{at} + b(e^{at} − 1)/a`, with `b = p γZ` and `a = p ηf − σ`
(`S0 + bt` when `a = 0`). Throws when `p` has adaptive immunity.
"""
function constant_pool_load(p::ImmuneFrogParameters, t::Real; S0::Real = 0.0)
    _has_immunity(p) && throw(ArgumentError("no closed form with adaptive immunity (θ, r ≠ 0)"))
    a, b = p.p * p.ηf - p.σ, p.p * p.γZ
    return a == 0 ? S0 + b * t : S0 * exp(a * t) + b * expm1(a * t) / a
end

"""
    time_to_lethal_load(p::ImmuneFrogParameters; S0 = 0) -> Float64

Day on which the load of [`constant_pool_load`](@ref) reaches `Smax`, or `Inf` if it never
does (the load then tends to `−b/a < Smax`). No adaptive immunity.
"""
function time_to_lethal_load(p::ImmuneFrogParameters; S0::Real = 0.0)
    _has_immunity(p) && throw(ArgumentError("no closed form with adaptive immunity (θ, r ≠ 0)"))
    S0 >= p.Smax && return 0.0
    a, b = p.p * p.ηf - p.σ, p.p * p.γZ
    a == 0 && return b > 0 ? (p.Smax - S0) / b : Inf
    ratio = (p.Smax + b / a) / (S0 + b / a)
    return ratio > 0 && log(ratio) / a > 0 ? log(ratio) / a : Inf
end

"""
    constant_pool_equilibria(p::ImmuneFrogParameters; Smax_factor = 10) -> Vector{Float64}

Equilibrium loads of a frog in a constant pool, sorted.

- Without adaptive immunity: `[−b/a]` when `a = p ηf − σ < 0`, otherwise empty (the load
  grows without bound).
- With adaptive immunity (needs `c > 0`): the roots of
  `p exp(−θ r S/c)(γZ + ηf S) = σS`, found by bracketing on `(0, Smax_factor·Smax]`; the
  antibody level is `A* = rS*/c`. Fig. 6C reports 916 sporangia, but `r` and `c` are not
  printed (target `woodhams2011_fig6c_stable_load` is pending).
"""
function constant_pool_equilibria(p::ImmuneFrogParameters; Smax_factor::Real = 10)
    if !_has_immunity(p)
        a, b = p.p * p.ηf - p.σ, p.p * p.γZ
        return a < 0 ? [-b / a] : Float64[]
    end
    p.c > 0 || throw(ArgumentError("an equilibrium with immunity needs c > 0"))
    k = p.θ * p.r / p.c
    g(S) = p.p * exp(-k * S) * (p.γZ + p.ηf * S) - p.σ * S
    grid = vcat(0.0, exp10.(range(-6, log10(Smax_factor * p.Smax); length = 4000)))
    roots = Float64[]
    for i in 1:(length(grid) - 1)
        lo, hi = grid[i], grid[i + 1]
        glo, ghi = g(lo), g(hi)
        glo == 0 && (push!(roots, lo); continue)
        sign(glo) == sign(ghi) && continue
        for _ in 1:200
            mid = (lo + hi) / 2
            (sign(g(mid)) == sign(glo)) ? (lo = mid) : (hi = mid)
        end
        push!(roots, (lo + hi) / 2)
    end
    return roots
end

"""
    FrogTreatment(time; γZ = nothing, p = nothing, S = nothing)

A treatment applied to a frog in a constant pool at `time` (days) (Woodhams et al. 2011,
Fig. 6B, C). Each given field replaces the current value:

- `γZ`: the encounter flux; `γZ = 0` removes the frog from the zoospore pool;
- `p`: the infection fraction with no antibodies, as after a probiotic-like constitutive
  defence (Fig. 6B(c): 0.15);
- `S`: the load; `S = 0` is an antifungal that clears the frog (Fig. 6C).
"""
struct FrogTreatment
    time::Float64
    γZ::Union{Nothing,Float64}
    p::Union{Nothing,Float64}
    S::Union{Nothing,Float64}
end
FrogTreatment(time::Real; γZ = nothing, p = nothing, S = nothing) = FrogTreatment(Float64(time),
    γZ === nothing ? nothing : Float64(γZ), p === nothing ? nothing : Float64(p),
    S === nothing ? nothing : Float64(S))

"""
    ConstantPoolFrogProblem(p::ImmuneFrogParameters; S0 = 0, A0 = 0, tspan = (0.0, 60.0),
                            treatments = FrogTreatment[])

One frog in a constant zoospore pool, with or without adaptive immunity (Woodhams et al.
2011, Appendix; Fig. 6B, C), state `[S, A]`:

```
dS/dt = ν(A) γZ + ν(A) ηf S − σS,   ν(A) = p exp(−θA)
dA/dt = rS − cA
```

`solve(prob)` returns `(sol, death_time)`. The frog dies, and the integration stops, when `S`
reaches `Smax`; `death_time` is `NaN` for a frog that survives. `treatments` are applied at
their times.
"""
struct ConstantPoolFrogProblem
    p::ImmuneFrogParameters
    u0::Vector{Float64}
    tspan::Tuple{Float64,Float64}
    treatments::Vector{FrogTreatment}
end
function ConstantPoolFrogProblem(p::ImmuneFrogParameters; S0::Real = 0.0, A0::Real = 0.0,
        tspan = (0.0, 60.0), treatments = FrogTreatment[])
    S0 >= 0 && A0 >= 0 || throw(ArgumentError("S0 and A0 must be non-negative"))
    return ConstantPoolFrogProblem(p, Float64[S0, A0], (Float64(tspan[1]), Float64(tspan[2])),
        sort(collect(FrogTreatment, treatments); by = t -> t.time))
end

function immune_frog_rhs!(du, u, q::Base.RefValue{ImmuneFrogParameters}, t)
    p = q[]
    S, A = u
    ν = infection_fraction(p, A)
    du[1] = ν * p.γZ + ν * p.ηf * S - p.σ * S
    du[2] = p.r * S - p.c * A
    return nothing
end

to_ode_problem(prob::ConstantPoolFrogProblem) =
    ODEProblem(immune_frog_rhs!, copy(prob.u0), prob.tspan, Ref(prob.p))

function CommonSolve.solve(prob::ConstantPoolFrogProblem, alg = Tsit5(); kwargs...)
    ode = to_ode_problem(prob)
    death = Ref(NaN)
    lethal = ContinuousCallback((u, t, integ) -> u[1] - integ.p[].Smax,
        integ -> (death[] = integ.t; SciMLBase.terminate!(integ)), nothing)
    times = [tr.time for tr in prob.treatments]
    function treat!(integ)
        for tr in prob.treatments
            tr.time == integ.t || continue
            pp = integ.p[]
            integ.p[] = _immune_with(pp; γZ = something(tr.γZ, pp.γZ), p = something(tr.p, pp.p))
            tr.S === nothing || (integ.u[1] = tr.S)
        end
        return nothing
    end
    treatcb = DiscreteCallback((u, t, integ) -> t in times, treat!)
    sol = solve(ode, alg; callback = CallbackSet(lethal, treatcb), tstops = times, kwargs...)
    return (sol = sol, death_time = death[])
end

"""
    ImmunePopulationParameters(; γ, η, f, σ, μ, Smax, p, θ = 0, r = 0, c = 0)
    ImmunePopulationParameters(:woodhams2011; set = "fig6d", μ, r, c)
    ImmunePopulationParameters(b::BriggsParameters; θ = 0, r = 0, c = 0)

Frogs sharing a dynamic zoospore pool, each with a load and an antibody level (Woodhams et
al. 2011, Appendix, "Population-level model of infection with adaptive immunity"), per day:

- `γ`: encounter rate between one zoospore and one frog, `γ/V` of the Briggs model;
- `η`: zoospores released per sporangium; `f`: fraction that immediately re-encounter the
  same frog; `σ`: sporangium loss; `μ`: zoospore death in the pool; `Smax`: lethal load;
- `p`, `θ`, `r`, `c`: as in [`ImmuneFrogParameters`](@ref).

The keyed constructor reads the Fig. 6D values; **`μ`, `r` and `c` are not printed** and must
be given. The third constructor turns a [`BriggsParameters`](@ref) into this form
(`γ = γ/V`, `p = ν`); with `θ = 0` the model is the Briggs et al. (2010) model.
"""
Base.@kwdef struct ImmunePopulationParameters
    γ::Float64
    η::Float64
    f::Float64
    σ::Float64
    μ::Float64
    Smax::Float64
    p::Float64
    θ::Float64 = 0.0
    r::Float64 = 0.0
    c::Float64 = 0.0
end

function ImmunePopulationParameters(key::Symbol; set::AbstractString, μ::Real, r::Real, c::Real)
    q = parameter_set(key; set)
    return ImmunePopulationParameters(; γ = q.γ, η = q.η, f = q.f, σ = q.σ, μ, Smax = q.Smax,
        p = q.p, θ = q.θ, r, c)
end

ImmunePopulationParameters(b::BriggsParameters; θ::Real = 0.0, r::Real = 0.0, c::Real = 0.0) =
    ImmunePopulationParameters(; γ = b.γ / b.V, η = b.η, f = b.f, σ = b.σ, μ = b.μ,
        Smax = b.Smax, p = b.ν, θ, r, c)

time_unit(::Type{ImmunePopulationParameters}) = :day

"""
    DensityReduction(time, N)

Reduce the number of live frogs to `N` at `time` (days), as in Woodhams et al. (2011) Fig.
6D (100 frogs reduced to 10 on day 30). The highest-numbered live frogs are removed.
"""
struct DensityReduction
    time::Float64
    N::Int
end

mutable struct ImmuneHostState
    alive::Vector{Bool}
    death_time::Vector{Float64}
    removal_time::Vector{Float64}
end

struct ImmunePopODEParams
    p::ImmunePopulationParameters
    hosts::ImmuneHostState
end

function immune_population_rhs!(du, u, q::ImmunePopODEParams, t)
    p, alive = q.p, q.hosts.alive
    N = length(alive)
    Z = u[end]
    load = zero(eltype(u))
    nalive = 0
    @inbounds for i in 1:N
        if alive[i]
            S, A = u[i], u[N + i]
            ν = infection_fraction(p, A)
            du[i] = ν * p.γ * Z + ν * p.η * p.f * S - p.σ * S
            du[N + i] = p.r * S - p.c * A
            load += S
            nalive += 1
        else
            du[i] = zero(eltype(u))
            du[N + i] = zero(eltype(u))
        end
    end
    du[end] = p.η * (1 - p.f) * load - p.γ * Z * nalive - p.μ * Z
    return nothing
end

"""
    ImmunePopulationProblem(p::ImmunePopulationParameters; N, S0 = zeros(N), A0 = zeros(N),
                            Z0 = 0, tspan = (0.0, 100.0), reductions = DensityReduction[])

Woodhams et al. (2011) population model, state `[S_1..S_N, A_1..A_N, Z]`:

```
dS_i/dt = ν(A_i) γ Z + ν(A_i) η f S_i − σ S_i
dA_i/dt = r S_i − c A_i
dZ/dt   = η(1 − f) Σ S_i − γ Z N_alive − μ Z
```

Frog `i` dies when `S_i` reaches `Smax`. `solve(prob)` returns
`(sol, death_time, removal_time, alive)`; times are `NaN` for frogs that did not die or were
not removed. With `θ = 0` it reproduces [`BriggsLoadProblem`](@ref).
"""
struct ImmunePopulationProblem
    p::ImmunePopulationParameters
    S0::Vector{Float64}
    A0::Vector{Float64}
    Z0::Float64
    tspan::Tuple{Float64,Float64}
    reductions::Vector{DensityReduction}
end
function ImmunePopulationProblem(p::ImmunePopulationParameters; N::Integer, S0 = zeros(N),
        A0 = zeros(N), Z0::Real = 0.0, tspan = (0.0, 100.0), reductions = DensityReduction[])
    length(S0) == N && length(A0) == N ||
        throw(ArgumentError("S0 and A0 must have length N"))
    return ImmunePopulationProblem(p, collect(float.(S0)), collect(float.(A0)), float(Z0),
        (Float64(tspan[1]), Float64(tspan[2])),
        sort(collect(DensityReduction, reductions); by = r -> r.time))
end

function to_ode_problem(prob::ImmunePopulationProblem)
    N = length(prob.S0)
    hosts = ImmuneHostState(fill(true, N), fill(NaN, N), fill(NaN, N))
    return ODEProblem(immune_population_rhs!, vcat(prob.S0, prob.A0, prob.Z0), prob.tspan,
        ImmunePopODEParams(prob.p, hosts))
end

function CommonSolve.solve(prob::ImmunePopulationProblem, alg = Tsit5(); kwargs...)
    ode = to_ode_problem(prob)
    N = length(prob.S0)
    condition(out, u, t, integ) = (alive = integ.p.hosts.alive;
        @inbounds for i in 1:N
            out[i] = alive[i] ? u[i] - integ.p.p.Smax : -one(eltype(u))
        end)
    function die!(integ, events)
        hosts = integ.p.hosts
        for i in eachindex(events)
            iszero(events[i]) && continue
            (hosts.alive[i] && integ.u[i] >= integ.p.p.Smax / 2) || continue
            hosts.alive[i] = false
            hosts.death_time[i] = integ.t
            integ.u[i] = 0.0
            integ.u[N + i] = 0.0
        end
        return nothing
    end
    times = [r.time for r in prob.reductions]
    function reduce!(integ)
        hosts = integ.p.hosts
        for red in prob.reductions
            red.time == integ.t || continue
            live = findall(hosts.alive)
            for i in reverse(live)
                count(hosts.alive) <= red.N && break
                hosts.alive[i] = false
                hosts.removal_time[i] = integ.t
                integ.u[i] = 0.0
                integ.u[N + i] = 0.0
            end
        end
        return nothing
    end
    cb = CallbackSet(VectorContinuousCallback(condition, die!, N),
        DiscreteCallback((u, t, integ) -> t in times, reduce!))
    sol = solve(ode, alg; callback = cb, tstops = times, kwargs...)
    h = ode.p.hosts
    return (sol = sol, death_time = copy(h.death_time), removal_time = copy(h.removal_time),
        alive = copy(h.alive))
end
