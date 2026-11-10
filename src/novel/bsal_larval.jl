# NEW MODEL (not from a published paper): the Schmidt et al. (2017) Bsal model with an aquatic
# larval stage that can be infected and can carry infection through metamorphosis.
# Wiki: lab/wiki/models/new-bsal-larval-carryover.md
#
# Why: Schmidt's model has adult females only. Wang et al. (2021) infected fire salamander
# larvae experimentally; the larvae showed no disease, and half carried the infection through
# metamorphosis and then died of it.
#
# State (S, L, I, W_S, W_I): susceptible, latent and infectious adults, and uninfected and
# infected larvae, per ha; years. With m = 1/larval_period, μ_W = larval_death and
# N = S + L + I:
#   dW_S/dt = b' S − (m + μ_W) W_S − λ_W W_S
#   dW_I/dt = λ_W W_S − (m + μ_W) W_I               λ_W = β_WA I + β_WW W_I
#   dS/dt   = m W_S + (1 − c) m W_I − (d + γN) S − β S I − β_AW S W_I
#   dL/dt   = c m W_I + β S I + β_AW S W_I − (d + e + γN) L
#   dI/dt   = e L − (d + d_I + γN) I
# b' = b (m + μ_W)/m, so that without Bsal recruitment to adults is b S and the equilibrium is
# Schmidt's (S = K). c is the carry-over probability. Infected larvae do not die of Bsal.
# As in Schmidt's model only susceptible adults breed.
# The larval period, larval death rate and all three larval transmission rates are
# ILLUSTRATIVE; none has a published value here. The years between metamorphosis and
# maturity are not modelled.

"""
    BsalLarvalStage(; carry_over, larval_period = 0.33, larval_death = 1.0,
                    adult_to_larva = 0, larva_to_larva = 0, larva_to_adult = 0)
    BsalLarvalStage(:wang2021; kwargs...)

A larval stage for the Schmidt Bsal model. `carry_over` is the probability that an infected
larva is still infected after metamorphosis (the keyed constructor takes it from the
compendium: 5 of 10 larvae in Wang et al. 2021). `larval_period` (yr) and `larval_death`
(yr⁻¹) describe the stage. The transmission rates (ha yr⁻¹, like Schmidt's `β`) are from
infectious adults to larvae, between larvae, and from infected larvae to adults.
"""
Base.@kwdef struct BsalLarvalStage
    carry_over::Float64
    larval_period::Float64 = 0.33
    larval_death::Float64 = 1.0
    adult_to_larva::Float64 = 0.0
    larva_to_larva::Float64 = 0.0
    larva_to_adult::Float64 = 0.0
    function BsalLarvalStage(carry_over, larval_period, larval_death, adult_to_larva,
            larva_to_larva, larva_to_adult)
        (0 <= carry_over <= 1 && larval_period > 0 && larval_death >= 0 &&
         adult_to_larva >= 0 && larva_to_larva >= 0 && larva_to_adult >= 0) ||
            throw(ArgumentError("need 0 ≤ carry_over ≤ 1, larval_period > 0 and rates ≥ 0"))
        return new(carry_over, larval_period, larval_death, adult_to_larva, larva_to_larva,
            larva_to_adult)
    end
end
function BsalLarvalStage(key::Symbol; kwargs...)
    v = parameter_set(key)
    return BsalLarvalStage(; carry_over = v.carried_over / v.n_larvae, kwargs...)
end

larval_exit(s::BsalLarvalStage) = 1 / s.larval_period + s.larval_death

"""Larvae per ha at the disease-free equilibrium, `b K × larval_period`."""
disease_free_larvae(p::SchmidtParameters, s::BsalLarvalStage) = p.b * p.K * s.larval_period

function bsal_larval_rhs!(du, u, q, t)
    p, s = q
    S, L, I, WS, WI = (max(x, 0.0) for x in u)
    m, out = 1 / s.larval_period, larval_exit(s)
    N = S + L + I
    dd = p.d + density_dependence(p) * N
    λW = s.adult_to_larva * I + s.larva_to_larva * WI
    fromW = s.larva_to_adult * S * WI
    du[1] = m * WS + (1 - s.carry_over) * m * WI - dd * S - p.β * S * I - fromW
    du[2] = s.carry_over * m * WI + p.β * S * I + fromW - (dd + p.e) * L
    du[3] = p.e * L - (dd + p.d_I) * I
    du[4] = p.b * out / m * S - out * WS - λW * WS
    du[5] = λW * WS - out * WI
    return nothing
end

"""
    basic_reproduction_number(p::SchmidtParameters, s::BsalLarvalStage)

R₀ of Bsal with the larval stage: the largest eigenvalue of the next-generation matrix for
the infected classes (infected larvae, latent adults, infectious adults) at the disease-free
state. Carry-over through metamorphosis counts as a transition, not a new infection. With no
larval transmission this is Schmidt's `K / K_threshold`.
"""
function basic_reproduction_number(p::SchmidtParameters, s::BsalLarvalStage)
    m, out, W = 1 / s.larval_period, larval_exit(s), disease_free_larvae(p, s)
    F = [s.larva_to_larva*W  0.0  s.adult_to_larva*W
         s.larva_to_adult*p.K  0.0  p.β*p.K
         0.0  0.0  0.0]
    V = [out  0.0  0.0
         -s.carry_over*m  p.b+p.e  0.0
         0.0  -p.e  p.b+p.d_I]
    return maximum(real, eigvals(F / V))
end

"""
    larval_reproduction_number(p, s)

Infected larvae produced by one infected larva among uninfected larvae, by larva-to-larva
transmission alone: `larva_to_larva × W* / (m + μ_W)`. Above 1, Bsal can persist in the
larvae without any adult.
"""
larval_reproduction_number(p::SchmidtParameters, s::BsalLarvalStage) =
    s.larva_to_larva * disease_free_larvae(p, s) / larval_exit(s)

"""
    BsalLarvalProblem(p::SchmidtParameters, stage; I0 = 0.01, tspan = (0.0, 10.0))

Schmidt's Bsal model with the larval `stage`, started at the disease-free equilibrium with
`I0` infectious adults per ha. Time in years.

`solve` returns series `t`, `adults`, `larvae`, `infectious` (adults), `larval_prevalence`,
and `min_adults`, `mean_adults` and `final_adults` as fractions of `K`.
"""
struct BsalLarvalProblem{P<:SchmidtParameters}
    p::P
    stage::BsalLarvalStage
    I0::Float64
    tspan::Tuple{Float64,Float64}
end
BsalLarvalProblem(p::SchmidtParameters, stage::BsalLarvalStage; I0 = 0.01,
    tspan = (0.0, 10.0)) =
    BsalLarvalProblem(p, stage, float(I0), (float(tspan[1]), float(tspan[2])))

to_ode_problem(prob::BsalLarvalProblem) = ODEProblem(bsal_larval_rhs!,
    [float(prob.p.K), 0.0, prob.I0, disease_free_larvae(prob.p, prob.stage), 0.0], prob.tspan,
    (prob.p, prob.stage))

function CommonSolve.solve(prob::BsalLarvalProblem, alg = Tsit5(); dt = 0.01, abstol = 1e-10,
        reltol = 1e-8, kwargs...)
    sol = solve(to_ode_problem(prob), alg; saveat = dt, abstol, reltol, kwargs...)
    SciMLBase.successful_retcode(sol) || error("Bsal larval model: solver failed")
    K = prob.p.K
    adults = [max(u[1], 0.0) + max(u[2], 0.0) + max(u[3], 0.0) for u in sol.u]
    larvae = [max(u[4], 0.0) + max(u[5], 0.0) for u in sol.u]
    return (t = sol.t, adults = adults, larvae = larvae,
        infectious = [max(u[3], 0.0) for u in sol.u],
        larval_prevalence = [l > 0 ? max(u[5], 0.0) / l : 0.0 for (u, l) in zip(sol.u, larvae)],
        min_adults = minimum(adults) / K, mean_adults = sum(adults) / length(adults) / K,
        final_adults = adults[end] / K)
end
