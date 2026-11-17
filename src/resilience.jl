# Glue to DynamicalResilience.jl and ResilienceCore.jl: which equilibria each model has and
# which of its terms are new infections. The resilience machinery itself lives there; this file
# only says where the models' equilibria are.

"""
    equilibrium_state(prob, state::Symbol) -> Vector{Float64}

A named equilibrium of a model, as a state vector in the order of its ODE:

| problem | states |
|---|---|
| [`SchmidtSLIProblem`](@ref) | `:endemic` (eq. 2.3; it exists only above the invasion threshold), `:disease_free` (`[K, 0, 0]`) |
| [`BsalLarvalProblem`](@ref) | `:disease_free` (`[K, 0, 0, W*, 0]` with `W* =` [`disease_free_larvae`](@ref)) |
| [`LoadTipping`](@ref) | `:low_load` (the stable low load), `:threshold` (the unstable load above which load runs away) |

Throws an `ArgumentError` if the model has no such state at its parameters (for example, no
endemic state below the invasion threshold, or no low load in the `:runaway` regime).
"""
function equilibrium_state end

function equilibrium_state(prob::SchmidtSLIProblem, state::Symbol)
    p = prob.p
    if state === :endemic
        eq = endemic_equilibrium(p)
        eq === nothing && throw(ArgumentError("no endemic state: K = $(p.K) is at or below the invasion " *
                                              "threshold $(invasion_threshold(p))"))
        return [eq.S, eq.L, eq.I]
    elseif state === :disease_free
        return [float(p.K), 0.0, 0.0]
    end
    throw(ArgumentError("a SchmidtSLIProblem has the states :endemic and :disease_free, not :$state"))
end

function equilibrium_state(prob::BsalLarvalProblem, state::Symbol)
    state === :disease_free ||
        throw(ArgumentError("a BsalLarvalProblem has the state :disease_free, not :$state"))
    return [float(prob.p.K), 0.0, 0.0, disease_free_larvae(prob.p, prob.stage), 0.0]
end

function equilibrium_state(m::LoadTipping, state::Symbol)
    eq = load_equilibria(m)
    S = state === :low_load ? eq.low : state === :threshold ? eq.threshold :
        throw(ArgumentError("a LoadTipping model has the states :low_load and :threshold, not :$state"))
    isfinite(S) || throw(ArgumentError("no :$state state in the $(tipping_regime(m)) regime"))
    return [S]
end

# The reduced one-host load equation dS/dt = ι + R S − c S/(h + S), with the model as parameters
load_tipping_reduced(u, m::LoadTipping, t) = [m.ι + net_growth_rate(m) * u[1] - m.c * u[1] / (m.h + u[1])]

"""
    to_ode_problem(m::LoadTipping; S0 = 0.0, tspan = (0.0, 365.0))

The one-host load equation `dS/dt = ι + (λ − σ + κ) S − c S/(h + S)` as an out-of-place
`ODEProblem` whose parameters are the `LoadTipping` model itself. Use it with DynamicalResilience's
`resilience_along`, varying a field with, for example,
`(m, ι) -> LoadTipping(ι, m.λ, m.σ, m.c, m.h, m.κ)`.
"""
to_ode_problem(m::LoadTipping; S0 = 0.0, tspan = (0.0, 365.0)) =
    ODEProblem(load_tipping_reduced, [float(S0)], float.(tspan), m)

"""
    linearise(prob::SchmidtSLIProblem; state = :endemic)  -> ContinuousJacobian
    linearise(prob::BsalLarvalProblem; state = :disease_free) -> ContinuousJacobian
    linearise(m::LoadTipping; state = :low_load)          -> ContinuousJacobian

ResilienceCore's `linearise` at a named equilibrium of a model (see [`equilibrium_state`](@ref)):
the Jacobian of its ODE there, with the equilibrium attached, so that every ResilienceCore metric
(`resilience`, `reactivity`, …) applies. The Schmidt and load-tipping Jacobians are analytic; the
larval model's is computed by automatic differentiation of its right-hand side. Time is in years
for the Bsal models and in days for the load model.
"""
function ResilienceCore.linearise(prob::SchmidtSLIProblem; state::Symbol = :endemic)
    u = equilibrium_state(prob, state)
    return ContinuousJacobian(schmidt_jacobian(prob.p, u), u)
end

function ResilienceCore.linearise(prob::BsalLarvalProblem; state::Symbol = :disease_free)
    u = equilibrium_state(prob, state)
    return ResilienceCore.linearise(to_ode_problem(prob), u)
end

function ResilienceCore.linearise(m::LoadTipping; state::Symbol = :low_load)
    u = equilibrium_state(m, state)
    S = u[1]
    return ContinuousJacobian(fill(net_growth_rate(m) - m.c * m.h / (m.h + S)^2, 1, 1), u)
end

"""
    reproduction_number(prob::SchmidtSLIProblem)  -> (R0, F, V, K)
    reproduction_number(prob::BsalLarvalProblem)  -> (R0, F, V, K)

R₀ by DynamicalResilience's general next-generation method, with the infection structure of each
model: for Schmidt's model the infected states are `L` and `I` and new infections `βSI` enter
`L`; with the larval stage they are `W_I`, `L` and `I`, new larval infections `λ_W W_S` enter
`W_I` and new adult infections `βSI + β_AW S W_I` enter `L` (carry-over through metamorphosis is a
transition). Equal to [`basic_reproduction_number`](@ref), which gives the same quantity in closed
form; this method also returns the next-generation matrix `K`.
"""
function DynamicalResilience.reproduction_number(prob::SchmidtSLIProblem)
    newinf(u, p, t) = [p.β * u[1] * u[3], zero(eltype(u))]
    return DynamicalResilience.reproduction_number(to_ode_problem(prob), equilibrium_state(prob, :disease_free);
                                                   infected = [2, 3], new_infections = newinf)
end

function DynamicalResilience.reproduction_number(prob::BsalLarvalProblem)
    function newinf(u, q, t)
        p, s = q
        S, L, I, WS, WI = u
        λW = s.adult_to_larva * I + s.larva_to_larva * WI
        return [λW * WS, p.β * S * I + s.larva_to_adult * S * WI, zero(eltype(u))]
    end
    return DynamicalResilience.reproduction_number(to_ode_problem(prob), equilibrium_state(prob, :disease_free);
                                                   infected = [5, 2, 3], new_infections = newinf)
end
