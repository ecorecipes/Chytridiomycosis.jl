# NEW MODEL (not from a published paper): a stochastic version of the nested model
# (nested.jl). Wiki: lab/wiki/models/new-nested-model.md
#
# Animals are whole numbers in an area of `area` hectares. Time advances in short steps of
# `dt` years. In each step every animal in a class leaves it with probability 1 − exp(−R dt),
# where R is the total per-animal rate out of that class, and goes to one of the destinations
# in proportion to their rates (a binomial draw then a multinomial split). Births are Poisson.
# The rates come from `nested_flows`, the same list of transitions the ODE uses, evaluated at
# densities (counts / area), so the two models agree in the limit of a large area.
#
# The zoospore pool, which holds very large numbers, is kept continuous: it is updated
# exactly over each step given the shedding at the start of the step. The pathogen counts as
# gone when no animal is infected and the pool would be expected to cause fewer than 0.05
# further infections; the pool is then emptied.
# The step dt is an approximation; dt = 1/1460 (six hours) is short against every rate in
# the model except the latent period (1 day), where it is a quarter of the mean.

"""
    StochasticNestedProblem(p::SchmidtParameters, biology::ChytridBiology,
        treatment = MicrobiomeTreatment(); area = 20.0, introduced = 1, dt = 1/1460,
        saveat = 1/52, kwargs...)

Stochastic version of [`NestedChytridProblem`](@ref): whole animals in `area` hectares,
starting at the disease-free state (rounded), with `introduced` infectious adults arriving at
`t_intro`. Other keywords (`treat`, `coverage`, `treat_times`, `retention`, `t_intro`,
`tspan`) are as for `NestedChytridProblem`; at a treatment each untreated animal is treated
with probability `coverage`. `biology.fadeout_density` is ignored: extinction happens by
itself.

`solve(prob; rng)` returns one run: series `t`, `adults`, `infected` (animals of any stage in
an infected class), and `mean_adults` and `final_adults` (fractions of the starting number of
adults, after `t_intro`), `peak_infected`, `pathogen_extinction` and `host_extinction` (times,
or `NaN`).
"""
struct StochasticNestedProblem{P<:SchmidtParameters}
    model::NestedChytridProblem{P}
    area::Float64
    introduced::Int
    dt::Float64
    saveat::Float64
end
function StochasticNestedProblem(p::SchmidtParameters, biology::ChytridBiology,
        treatment::MicrobiomeTreatment = MicrobiomeTreatment(); area = 20.0,
        introduced::Integer = 1, dt = 1 / 1460, saveat = 1 / 52, kwargs...)
    (area > 0 && introduced >= 0 && dt > 0 && saveat >= dt) ||
        throw(ArgumentError("need area > 0, introduced ≥ 0 and saveat ≥ dt > 0"))
    model = NestedChytridProblem(p, biology, treatment; kwargs...)
    return StochasticNestedProblem(model, float(area), Int(introduced), float(dt), float(saveat))
end

function CommonSolve.solve(prob::StochasticNestedProblem; rng::AbstractRNG = Random.default_rng())
    m = prob.model
    p, bio, A, dt = m.p, m.biology, prob.area, prob.dt
    (t0, t1) = m.tspan
    clock = Ref(t0)
    q = (p, bio, m.treatment, m.retention, clock)
    n = [round(Int, x * A) for x in nested_disease_free(p, bio)[1:16]]
    adults0 = sum(n[nested_index(2, τ, k)] for τ in 1:2, k in 1:4)
    Z = 0.0
    μZ = log(2) / bio.pool_half_life
    φ = bio.environment ? bio.environmental_fraction : 0.0
    infected_idx = [nested_index(s, τ, k) for s in 1:2 for τ in 1:2 for k in 2:4]
    susceptible_idx = [nested_index(s, τ, 1) for s in 1:2 for τ in 1:2]
    stages = m.treat === :adults ? (2,) : m.treat === :larvae ? (1,) : (1, 2)

    events = Tuple{Float64,Symbol}[(t, :treat) for t in m.treat_times if t0 <= t < t1]
    t0 <= m.t_intro < t1 && push!(events, (m.t_intro, :introduce))
    sort!(events; by = e -> (e[1], e[2] === :introduce))

    u = zeros(NESTED_Z)
    rates = [Float64[] for _ in 1:16]
    dests = [Int[] for _ in 1:16]
    births = Tuple{Int,Float64}[]
    Δ = zeros(Int, 16)
    ts, adults, infected = Float64[], Int[], Int[]
    pathogen_gone, host_gone = NaN, NaN
    introduced = false
    peak = 0
    adult_count() = sum(n[nested_index(2, τ, k)] for τ in 1:2, k in 1:4)
    infected_count() = sum(n[i] for i in infected_idx)
    function record!(t)
        push!(ts, t); push!(adults, adult_count()); push!(infected, infected_count())
    end

    t, next_save, e = t0, t0, 1
    while t < t1 - dt / 2
        while e <= length(events) && events[e][1] <= t + dt / 2
            te, kind = events[e]
            if kind === :treat
                clock[] = te
                for s in stages, k in 1:4
                    i, j = nested_index(s, 1, k), nested_index(s, 2, k)
                    moved = rand(rng, Binomial(n[i], m.coverage))
                    n[i] -= moved; n[j] += moved
                end
            elseif isnan(host_gone)
                n[nested_index(2, 1, 4)] += prob.introduced
                introduced = prob.introduced > 0
            end
            e += 1
        end
        if t >= next_save - dt / 2
            record!(t); next_save += prob.saveat
        end
        # rates at the start of the step
        @inbounds for i in 1:16
            u[i] = n[i] / A
            empty!(rates[i]); empty!(dests[i])
        end
        u[NESTED_Z] = Z
        empty!(births)
        dZ = nested_flows(u, q, t) do from, to, f
            f > 0 || return nothing
            if from == 0
                push!(births, (to, f))
            elseif n[from] > 0
                push!(rates[from], f / u[from]); push!(dests[from], to)
            end
            return nothing
        end
        # draw the events of the step
        fill!(Δ, 0)
        @inbounds for i in 1:16
            (n[i] > 0 && !isempty(rates[i])) || continue
            R = sum(rates[i])
            leave = rand(rng, Binomial(n[i], -expm1(-R * dt)))
            Δ[i] -= leave
            rest, Rrest = leave, R
            for j in eachindex(rates[i])
                rest == 0 && break
                k = j == lastindex(rates[i]) ? rest :
                    rand(rng, Binomial(rest, clamp(rates[i][j] / Rrest, 0.0, 1.0)))
                dests[i][j] > 0 && (Δ[dests[i][j]] += k)
                rest -= k
                Rrest -= rates[i][j]
            end
        end
        for (to, f) in births
            Δ[to] += rand(rng, Poisson(f * A * dt))
        end
        n .+= Δ
        if bio.environment                       # exact over the step, shedding held fixed
            HA = dZ + μZ * Z
            decay = exp(-μZ * dt)
            Z = Z * decay + HA / μZ * (1 - decay)
        end
        t += dt
        peak = max(peak, infected_count())
        # extinctions
        if isnan(host_gone) && sum(n) == 0
            host_gone = t
            isnan(pathogen_gone) && introduced && (pathogen_gone = t)
            Z = 0.0
        end
        if introduced && isnan(pathogen_gone) && infected_count() == 0
            S = sum(n[i] for i in susceptible_idx)
            if p.β * φ * Z * S < 0.05             # expected infections still to come from the pool
                pathogen_gone, Z = t, 0.0
            end
        end
    end
    record!(t)
    after = findall(>=(m.t_intro), ts)
    rel = adults[after] ./ max(adults0, 1)
    return (t = ts, adults = adults, infected = infected,
        mean_adults = sum(rel) / length(rel), final_adults = rel[end], peak_infected = peak,
        pathogen_extinction = pathogen_gone, host_extinction = host_gone)
end
