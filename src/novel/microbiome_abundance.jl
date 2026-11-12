# NEW MODEL (not from a published paper): bacterial abundance as a state of each frog, on the
# load tipping model (load_tipping.jl). Wiki: lab/wiki/models/new-microbiome-abundance.md
#
# The treated/untreated models give every treated animal the same clock. Here each of N frogs
# sharing a pool carries its own abundance B_i of the applied bacteria, in dose units (B = 1
# just after a bath). State (S_1..S_N, B_1..B_N, Z); days.
#   dS_i/dt = infection from the pool + (ηνf − σ) S_i − c_i S_i/(h + S_i)        (load tipping)
#   c_i     = c + (c_treated − c) B_i/(B_i + 1/F)       F = protective_fold
#   dB_i/dt = −δ B_i + x (B̄ − B_i)                      δ = log 2 / abundance_half_life
#   dZ/dt   = as in the load tipping model
# B̄ is the mean abundance over live frogs and x (`direct`, per day) is the rate at which frogs
# exchange bacteria through the pond. Exchange moves bacteria between frogs without creating
# any, so the mean declines at δ whatever x is: the bacteria spread but do not sustain
# themselves. That is what the one field trial showed (Knapp et al. 2022). A bath sets B_i = 1
# for the frogs in the bathed group (the first `coverage × N` frogs).
# c_treated, F and x are ILLUSTRATIVE. The half-life can be taken from our fits to open data
# (4 to 5 days in the field), which are not published estimates.

struct AbundanceODEParams{P<:BriggsParameters}
    p::P
    hosts::BriggsHostState
    c::Vector{Float64}
    c_treated::Vector{Float64}
    h::Float64
    ε::Float64
    decay::Float64
    exchange::Float64
    B50::Float64
end

function abundance_rhs!(du, u, q::AbundanceODEParams, t)
    p, alive = q.p, q.hosts.alive
    N = (length(u) - 1) ÷ 2
    Z = u[end]
    infection = transmission_rate(p) / p.V * Z
    selfgrowth = p.η * p.ν * p.f - p.σ
    load, total, nalive = 0.0, 0.0, 0
    @inbounds for i in 1:N
        alive[i] || continue
        total += max(u[N + i], 0.0)
        nalive += 1
    end
    Bbar = nalive > 0 ? total / nalive : 0.0
    @inbounds for i in 1:N
        if alive[i]
            S, B = max(u[i], 0.0), max(u[N + i], 0.0)
            c = q.c[i] + (q.c_treated[i] - q.c[i]) * (q.B50 > 0 ? B / (B + q.B50) : (B > 0 ? 1.0 : 0.0))
            du[i] = infection + selfgrowth * S - c * S / (q.h + S)
            du[N + i] = -q.decay * B + q.exchange * (Bbar - B)
            load += S
        else
            du[i] = du[N + i] = 0.0
        end
    end
    du[end] = q.ε + p.η * (1 - p.f) * load - p.μ * Z - p.γ / p.V * Z * nalive
    return nothing
end

"""
    ProbioticAbundanceProblem(p::BriggsParameters, treatment; c, c_treated, h, N, ε = 0,
                              coverage = 0, dose_days = [0], days = 365)

Load tipping model for `N` frogs sharing a pool, in which each frog carries its own abundance
of the applied bacteria. A bath on each of `dose_days` sets the abundance of the bathed group
(the first `coverage × N` frogs) to 1. Abundance declines with
`treatment.abundance_half_life` and is exchanged between frogs at rate `treatment.direct`;
clearance rises from `c` toward `c_treated` as `B/(B + 1/treatment.protective_fold)`. `c` and
`c_treated` are one value for every frog or a vector with one value per frog. Time in days.
The other fields of `treatment` are not used.

`solve` returns daily `t`, `alive`, `mean_load`, `abundance_bathed` and `abundance_unbathed`
(means over live frogs in each group), and `death_day` for each frog (`0` for survivors).
"""
struct ProbioticAbundanceProblem{P<:BriggsParameters}
    p::P
    treatment::MicrobiomeTreatment
    c::Vector{Float64}
    c_treated::Vector{Float64}
    h::Float64
    N::Int
    ε::Float64
    coverage::Float64
    dose_days::Vector{Float64}
    days::Int
end
function ProbioticAbundanceProblem(p::BriggsParameters, treatment::MicrobiomeTreatment; c,
        c_treated, h, N::Integer, ε = 0.0, coverage = 0.0, dose_days = [0], days::Integer = 365)
    N > 0 || throw(ArgumentError("need N > 0"))
    cs, cts = per_host(c, N, "c"), per_host(c_treated, N, "c_treated")
    (all(>=(0), cs) && all(>=(0), cts) && h > 0 && ε >= 0 && 0 <= coverage <= 1) ||
        throw(ArgumentError("need c, c_treated, ε ≥ 0, h > 0, N > 0, coverage in [0, 1]"))
    return ProbioticAbundanceProblem(p, treatment, cs, cts, float(h), Int(N),
        float(ε), float(coverage), sort!(collect(Float64, dose_days)), Int(days))
end

function CommonSolve.solve(prob::ProbioticAbundanceProblem, alg = Tsit5(); abstol = 1e-6,
        reltol = 1e-6, kwargs...)
    p, tr, N = prob.p, prob.treatment, prob.N
    hosts = BriggsHostState(fill(true, N), fill(NaN, N))
    q = AbundanceODEParams(p, hosts, prob.c, prob.c_treated, prob.h, prob.ε,
        log(2) / tr.abundance_half_life, tr.direct, 1 / tr.protective_fold)
    bathed = 1:round(Int, prob.coverage * N)
    u = zeros(2N + 1)
    u[end] = prob.ε / (p.μ + p.γ * N / p.V)
    stops = sort!(unique!(vcat(0.0, [d for d in prob.dose_days if 0 <= d < prob.days],
        float(prob.days))))
    ts, us = Float64[], Vector{Float64}[]
    for k in 1:(length(stops) - 1)
        t0, t1 = stops[k], stops[k + 1]
        if t0 in prob.dose_days
            for i in bathed
                hosts.alive[i] && (u[N + i] = 1.0)
            end
        end
        sol = solve(ODEProblem(abundance_rhs!, copy(u), (t0, t1), q), alg;
            callback = briggs_death_callback(N), saveat = 1.0, abstol, reltol, kwargs...)
        SciMLBase.successful_retcode(sol) || error("abundance model: solver failed after day $t0")
        first = k == 1 ? 1 : 2                               # the stop itself is already saved
        append!(ts, sol.t[first:end])
        append!(us, sol.u[first:end])
        u .= max.(sol.u[end], 0.0)
    end
    death = hosts.death_time
    living(t) = [i for i in 1:N if isnan(death[i]) || death[i] > t]
    avg(v, idx) = isempty(idx) ? NaN : sum(v[i] for i in idx) / length(idx)
    unbathed = (length(bathed) + 1):N
    alive = [living(t) for t in ts]
    out = (t = ts, alive = length.(alive),
        mean_load = [avg(v, a) for (v, a) in zip(us, alive)],
        abundance_bathed = [avg(v, [N + i for i in a if i in bathed]) for (v, a) in zip(us, alive)],
        abundance_unbathed = [avg(v, [N + i for i in a if i in unbathed]) for (v, a) in zip(us, alive)])
    death_day = [isnan(d) ? 0 : ceil(Int, d) for d in death]
    return merge(out, (death_day = death_day,))
end
