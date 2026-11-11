# NEW MODEL (not from a published paper): a stochastic version of the frog load models
# (load_tipping.jl, microbiome_abundance.jl). Wiki: lab/wiki/models/new-load-tipping.md
#
# Sporangia on each frog and zoospores in the shared pool are whole numbers. In each step of
# `dt` days, with β = νγ (transmission_rate):
#   new sporangia on frog i   Poisson((β/V) Z dt + ηνf S_i dt)    infection from the pool and
#                                                                  self-reinfection
#   sporangia lost            Binomial(S_i, 1 − exp(−(σ + c_i/(h + S_i)) dt))
#   zoospores released        Poisson(η(1 − f) S_i dt)
#   zoospores lost            Binomial(Z, 1 − exp(−(μ + γN/V) dt)),  N the frogs alive
#   zoospores from outside    Poisson(ε dt)
# These are the rates of the load tipping ODE, so the means agree. A frog dies when its load
# reaches Smax. Zoospores that infect are not removed one by one (in the ODE, encounters
# remove zoospores whether or not they infect, at γN/V).
#
# Applied bacteria are carried per frog as in ProbioticAbundanceProblem: abundance B_i
# (continuous; a bath sets it to 1) declines with `abundance_half_life`, frogs exchange it at
# rate `direct`, and clearance rises from c_i toward c_treated as B/(B + 1/F). In addition a
# frog loses its bacteria entirely at random with `half_life` (the treated class of
# ProbioticLoadProblem). Frogs may differ in c (a vector).

"""
    StochasticFrogProblem(p::BriggsParameters, treatment = MicrobiomeTreatment(); c,
        c_treated = c, h, N, ε = 0, coverage = 0, dose_days = [0], days = 365, dt = 0.05)

`c_treated` is likewise one value or a vector; by default it equals `c` (no effect).

Stochastic load tipping model for `N` frogs sharing a pool: whole sporangia and zoospores,
steps of `dt` days. `c` is one value or a vector (frogs that differ). On each of `dose_days`
the bathed group (the first `coverage × N` frogs) gets bacterial abundance 1, which raises its
clearance toward `c_treated` as described for [`ProbioticAbundanceProblem`](@ref), using
`treatment.abundance_half_life`, `treatment.protective_fold`, `treatment.direct` (exchange)
and `treatment.half_life` (random loss of all bacteria). With the defaults the bacteria last
for life.

`solve(prob; rng)` returns daily `t`, `alive` and `mean_load` (over live frogs), and
`death_day` for each frog (`0` for survivors).
"""
struct StochasticFrogProblem{P<:BriggsParameters}
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
    dt::Float64
end
function StochasticFrogProblem(p::BriggsParameters, treatment::MicrobiomeTreatment = MicrobiomeTreatment();
        c, c_treated = nothing, h, N::Integer, ε = 0.0, coverage = 0.0, dose_days = [0],
        days::Integer = 365, dt = 0.05)
    cs = per_host(c, N, "c")
    ct = c_treated === nothing ? copy(cs) : per_host(c_treated, N, "c_treated")
    (all(>=(0), cs) && all(>=(0), ct) && h > 0 && ε >= 0 && 0 <= coverage <= 1 && N > 0 && dt > 0) ||
        throw(ArgumentError("need c, c_treated, ε ≥ 0, h, dt > 0, N > 0, coverage in [0, 1]"))
    return StochasticFrogProblem(p, treatment, cs, ct, float(h), Int(N), float(ε),
        float(coverage), collect(Float64, dose_days), Int(days), float(dt))
end

function CommonSolve.solve(prob::StochasticFrogProblem; rng::AbstractRNG = Random.default_rng())
    p, tr, N, dt = prob.p, prob.treatment, prob.N, prob.dt
    β = transmission_rate(p)
    S = zeros(Int, N)
    alive = trues(N)
    B = zeros(N)
    Z = round(Int, prob.ε / (p.μ + p.γ * N / p.V))
    bathed = 1:round(Int, prob.coverage * N)
    decay = isinf(tr.abundance_half_life) ? 1.0 : exp(-log(2) / tr.abundance_half_life * dt)
    mixing = exp(-tr.direct * dt)
    lose = -expm1(-reversion_rate(tr) * dt)
    B50 = isinf(tr.protective_fold) ? 0.0 : 1 / tr.protective_fold
    death_day = zeros(Int, N)
    out = (t = Float64[], alive = Int[], mean_load = Float64[])
    steps_per_day = max(1, round(Int, 1 / dt))
    h = 1 / steps_per_day
    doses = Set(round.(Int, prob.dose_days))
    for day in 0:(prob.days - 1)
        if day in doses
            for i in bathed
                alive[i] && (B[i] = 1.0)
            end
        end
        for _ in 1:steps_per_day
            n = count(alive)
            n == 0 && break
            released = 0
            infect = β / p.V * Z * h
            for i in 1:N
                alive[i] || continue
                b = B[i]
                c = prob.c[i] + (prob.c_treated[i] - prob.c[i]) * (B50 > 0 ? b / (b + B50) : (b > 0 ? 1.0 : 0.0))
                s = S[i]
                gain = rand(rng, Poisson(infect + p.η * p.ν * p.f * s * h))
                loss = s > 0 ? rand(rng, Binomial(s, -expm1(-(p.σ + c / (prob.h + s)) * h))) : 0
                released += s > 0 ? rand(rng, Poisson(p.η * (1 - p.f) * s * h)) : 0
                S[i] = s + gain - loss
                if S[i] >= p.Smax
                    alive[i] = false
                    S[i] = 0
                    B[i] = 0.0
                    death_day[i] = day + 1
                end
            end
            lost = Z > 0 ? rand(rng, Binomial(Z, -expm1(-(p.μ + p.γ * n / p.V) * h))) : 0
            Z += released - lost + rand(rng, Poisson(prob.ε * h))
            # bacteria: exchange toward the mean, decline, and random loss of all of them
            if any(>(0), B)
                live = findall(alive)
                Bbar = isempty(live) ? 0.0 : sum(B[live]) / length(live)
                for i in live
                    B[i] = (Bbar + (B[i] - Bbar) * mixing) * decay
                    lose > 0 && rand(rng) < lose && (B[i] = 0.0)
                end
            end
        end
        n = count(alive)
        push!(out.t, day + 1); push!(out.alive, n)
        push!(out.mean_load, n > 0 ? sum(S[alive]) / n : NaN)
    end
    return merge(out, (death_day = death_day,))
end
