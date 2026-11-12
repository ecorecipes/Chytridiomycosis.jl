# Briggs, Knapp & Vredenburg (2010) PNAS 107:9695. Bd sporangial load on each frog plus a
# shared zoospore pool; a frog dies when its load reaches Smax.
# Units: time in days, loads in sporangia per frog, zoospores per pool of volume V.
# Wiki: lab/wiki/models/bd-briggs2010-load.md
#
#   dS_i/dt = (νγ/V) Z + ηνf S_i − σ S_i                        (main text; β = νγ)
#   dZ/dt   = η(1−f) Σ_alive S_i − μZ − (γ/V) Z N_alive          (checked against the PNAS SI;
#                                                                 the review's β is γ)

"""
    BriggsParameters(; V, γ, ν, f, μ, η, σ, Smax)
    BriggsParameters(:briggs2010)

Pool volume `V`; zoospore–frog encounter rate `γ` (volume day⁻¹); establishment probability
`ν` (so transmission is `β = νγ`); self-reinfection fraction `f`; zoospore death rate `μ`,
zoospore release per sporangium `η` and sporangium loss `σ` (day⁻¹); lethal load `Smax`
(sporangia).
"""
Base.@kwdef struct BriggsParameters{T<:Real}
    V::T
    γ::T
    ν::T
    f::T
    μ::T
    η::T
    σ::T
    Smax::T
end

function BriggsParameters(key::Symbol)
    p = parameter_set(key)
    return BriggsParameters(; V = p.V, γ = p.γ, ν = p.ν, f = p.f, μ = p.μ, η = p.η, σ = p.σ,
        Smax = p.Smax)
end

"""Successful transmission rate from the pool, `β = νγ`."""
transmission_rate(p::BriggsParameters) = p.ν * p.γ

"""Jacobian of the linearised (S, Z) system for `N` identical frogs (the symmetric mode)."""
function briggs_linear_jacobian(p::BriggsParameters, N)
    a = p.η * p.ν * p.f - p.σ
    return [a                         transmission_rate(p)/p.V
            p.η*(1-p.f)*N             -p.μ-p.γ*N/p.V]
end

"""
    load_growth_rate(p, N)

Asymptotic per-day growth rate of Bd load on `N` identical frogs (the dominant eigenvalue of
the linear system). Loads grow toward `Smax` when it is positive and decline toward clearance
when it is negative.
"""
load_growth_rate(p::BriggsParameters, N) = maximum(real, eigvals(briggs_linear_jacobian(p, N)))

"""
    host_density_threshold(p) -> N*

Number of frogs per pool above which Bd load grows (`det J = 0`):
`N* = -aμ / ((aγ + νγη(1-f))/V)`, with `a = ηνf - σ`. Returns `0` when `a ≥ 0` (loads grow
from self-reinfection alone) and `Inf` when no density suffices.
"""
function host_density_threshold(p::BriggsParameters)
    a = p.η * p.ν * p.f - p.σ
    a >= 0 && return zero(a)
    denom = (a * p.γ + transmission_rate(p) * p.η * (1 - p.f)) / p.V
    return denom > 0 ? -a * p.μ / denom : oftype(a, Inf)
end

mutable struct BriggsHostState
    alive::Vector{Bool}
    death_time::Vector{Float64}
end

struct BriggsODEParams{P<:BriggsParameters}
    p::P
    hosts::BriggsHostState
end

function briggs_rhs!(du, u, q::BriggsODEParams, t)
    p, alive = q.p, q.hosts.alive
    N = length(u) - 1
    Z = u[end]
    infection = transmission_rate(p) / p.V * Z
    selfgrowth = p.η * p.ν * p.f - p.σ
    load = zero(eltype(u))
    nalive = 0
    @inbounds for i in 1:N
        if alive[i]
            du[i] = infection + selfgrowth * u[i]
            load += u[i]
            nalive += 1
        else
            du[i] = zero(eltype(u))
        end
    end
    du[end] = p.η * (1 - p.f) * load - p.μ * Z - p.γ / p.V * Z * nalive
    return nothing
end

"""
    BriggsLoadProblem(p; N, S0 = zeros(N), Z0 = 0.0, tspan = (0.0, 365.0))

Deterministic Briggs et al. (2010) model for `N` frogs with initial loads `S0` and `Z0`
zoospores in the pool. Each frog dies (its load is removed and frozen at 0) when its load
reaches `Smax`; `solve` returns `(sol, death_time, alive)`, with `death_time[i] = NaN` for
survivors.
"""
struct BriggsLoadProblem{P<:BriggsParameters,T}
    p::P
    S0::Vector{Float64}
    Z0::Float64
    tspan::T
end
function BriggsLoadProblem(p::BriggsParameters; N::Integer, S0 = zeros(N), Z0 = 0.0,
                           tspan = (0.0, 365.0))
    length(S0) == N || throw(ArgumentError("length(S0) must equal N"))
    return BriggsLoadProblem(p, collect(float.(S0)), float(Z0), float.(tspan))
end

function to_ode_problem(prob::BriggsLoadProblem)
    N = length(prob.S0)
    hosts = BriggsHostState(fill(true, N), fill(NaN, N))
    return ODEProblem(briggs_rhs!, vcat(prob.S0, prob.Z0), prob.tspan,
        BriggsODEParams(prob.p, hosts))
end

function briggs_death_callback(N)
    condition(out, u, t, integ) = (alive = integ.p.hosts.alive;
        @inbounds for i in 1:N
            out[i] = alive[i] ? u[i] - integ.p.p.Smax : -one(eltype(u))
        end)
    # SciMLBase 3 has no separate affect_neg!: act only on an upward crossing by a live frog.
    # DiffEqBase passes an Int8 vector with one entry per frog: 0 = no event, otherwise the
    # sign encodes the crossing direction. Identical frogs cross simultaneously.
    function affect!(integ, events)
        hosts = integ.p.hosts
        for i in eachindex(events)
            iszero(events[i]) && continue
            (hosts.alive[i] && integ.u[i] >= integ.p.p.Smax / 2) || continue
            hosts.alive[i] = false
            hosts.death_time[i] = integ.t
            integ.u[i] = 0.0
        end
        return nothing
    end
    return VectorContinuousCallback(condition, affect!, N)
end

function CommonSolve.solve(prob::BriggsLoadProblem, alg = Tsit5(); kwargs...)
    ode = to_ode_problem(prob)
    cb = briggs_death_callback(length(prob.S0))
    sol = solve(ode, alg; callback = cb, kwargs...)
    hosts = ode.p.hosts
    return (sol = sol, death_time = copy(hosts.death_time), alive = copy(hosts.alive))
end
