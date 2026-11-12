# Briggs, Knapp & Vredenburg (2010) PNAS 107:9695, the stochastic seasonal model with the
# full R. muscosa / R. sierrae stage structure (SI Methods, "Model with full R. muscosa/R.
# sierrae stage structure"; Fig. 4). Six host stages: three tadpole year classes (tad1-3), two
# subadult year classes (sub1-2) and adults. Each stage has its own γ, ν, f, Smax and
# overwinter survival θ. Tadpoles are capped at Smax but never killed by Bd; subadults and
# adults die when their load reaches Smax.
# Units: days, sporangia per host, zoospores per pool. Wiki: lab/wiki/models/bd-briggs2010-load.md
#
# Summer (100 days): the SI Eqs. S1-S2 step of briggs2010_seasonal.jl, summed over stages.
# Spring, in this order (the SI lists the events but not their order; see below):
#   1. winter survival: each host of stage x with probability θx, zoospores with θZ;
#   2. reproduction: each surviving adult reproduces with probability pF, producing
#      Poisson(R) tadpoles (held back, they enter tad1 at step 4);
#   3. maturation, carrying loads: open adult slots (K minus adults) are filled by sub2 chosen
#      at random, other sub2 are lost; sub1 -> sub2; each tad2 enters metamorphosis with
#      probability m (otherwise -> tad3), and each tad3 metamorphoses; a metamorphosing
#      tadpole reaches sub1 with probability ω_metamorph (otherwise dies) and loses its
#      infection; tad1 -> tad2;
#   4. the newborns enter tad1, uninfected.
# Choices the paper leaves open (each a keyword or stated here):
#   order of spring events               -> as above; adults recruited this spring breed next year
#   tad3 metamorphosis                   -> all tad3 attempt it, succeeding with ω_metamorph
#   how Bd invades in year 20            -> one random adult (else the oldest stage present)
#                                           gets S0 = 10 sporangia at the start of that summer
#   initial population                   -> disease_free_stages(p), rounded
#   Δt, Poisson cut-off                  -> 1 day, 100 (as in briggs2010_seasonal.jl)

const BRIGGS_STAGE_NAMES = (:tad1, :tad2, :tad3, :sub1, :sub2, :adult)

"""
    BriggsLifeStages(; V, μ, η, σ, tad1, tad2, tad3, sub1, sub2, adult, θZ, m, ω, pF, R, K,
                     εZ = 0.0, summer = 100.0)
    BriggsLifeStages(:briggs2010; set = "fig4a" | "fig4b" | "fig4c")

Parameters of the six-stage seasonal model of Briggs et al. (2010, Fig. 4): pool volume `V`,
zoospore death `μ`, zoospore release `η` and sporangium loss `σ` (day⁻¹); one
[`BriggsStage`](@ref) for each of the stages `tad1`, `tad2`, `tad3` (tadpole year classes,
non-lethal load cap), `sub1`, `sub2` (subadult year classes) and `adult`; zoospore overwinter
survival `θZ`; probability `m` that a second-year tadpole enters metamorphosis; probability
`ω` (ω_metamorph) that a metamorphosing tadpole survives it; probability `pF` that an adult
breeds; mean clutch `R` (tadpoles per breeding adult, Poisson); adult ceiling `K`; external
zoospores `εZ` per day; summer length `summer` (days).
"""
Base.@kwdef struct BriggsLifeStages
    V::Float64
    μ::Float64
    η::Float64
    σ::Float64
    tad1::BriggsStage
    tad2::BriggsStage
    tad3::BriggsStage
    sub1::BriggsStage
    sub2::BriggsStage
    adult::BriggsStage
    θZ::Float64
    m::Float64
    ω::Float64
    pF::Float64
    R::Float64
    K::Int
    εZ::Float64 = 0.0
    summer::Float64 = 100.0
end

function BriggsLifeStages(key::Symbol; set::AbstractString)
    p = parameter_set(key; set)
    tad = BriggsStage(; γ = p.γ_tadpole, ν = p.ν, f = p.f, Smax = p.Smax_tadpole,
        θ = p.θ_tadpole, lethal = false)
    sub = BriggsStage(; γ = p.γ_subadult, ν = p.ν, f = p.f, Smax = p.Smax_subadult,
        θ = p.θ_subadult)
    adult = BriggsStage(; γ = p.γ_adult, ν = p.ν, f = p.f, Smax = p.Smax_adult, θ = p.θ_adult)
    return BriggsLifeStages(; V = p.V, μ = p.μ, η = p.η, σ = p.σ, tad1 = tad, tad2 = tad,
        tad3 = tad, sub1 = sub, sub2 = sub, adult, θZ = p.θZ, m = p.m, ω = p.ω_metamorph,
        pF = p.pF, R = p.R, K = Int(p.K), εZ = get(p, :εZ, 0.0), summer = p.summer)
end

life_stages(p::BriggsLifeStages) = (p.tad1, p.tad2, p.tad3, p.sub1, p.sub2, p.adult)

"""
    disease_free_stages(p::BriggsLifeStages) -> NamedTuple

Expected end-of-summer numbers in each stage (`tad1`, `tad2`, `tad3`, `sub1`, `sub2`,
`adult`) of the Bd-free population at its steady state, with the adult class full (`K`).
Throws if subadult recruitment cannot replace adult winter mortality (no Bd-free steady
state with adults at `K`).
"""
function disease_free_stages(p::BriggsLifeStages)
    a = float(p.K)
    tad1 = p.pF * p.R * p.adult.θ * a
    tad2 = p.tad1.θ * tad1
    tad3 = (1 - p.m) * p.tad2.θ * tad2
    sub1 = p.ω * (p.m * p.tad2.θ * tad2 + p.tad3.θ * tad3)
    sub2 = p.sub1.θ * sub1
    p.sub2.θ * sub2 >= (1 - p.adult.θ) * a ||
        throw(ArgumentError("recruitment cannot keep the adult class at K"))
    return (tad1 = tad1, tad2 = tad2, tad3 = tad3, sub1 = sub1, sub2 = sub2, adult = a)
end

"""
    BriggsStagedProblem(p::BriggsLifeStages; years = 50, invade_year = 20, S0 = 10.0,
                        Δt = 1.0, poisson_below = 100.0, init = disease_free_stages(p),
                        Z0 = 0.0, record_year = nothing)

The six-stage stochastic seasonal model of Briggs et al. (2010, Fig. 4). The population starts
Bd-free with `init` hosts per stage (rounded); at the start of summer `invade_year` one random
adult (or, if there are none, a host of the oldest stage present) acquires `S0` sporangia.

`solve(prob; rng)` returns a NamedTuple with end-of-summer yearly series: `counts` and
`mean_load` (each a NamedTuple of vectors over the six stages; the mean is over all hosts of
the stage, as in Fig. 4's "Ave. Sᵢ"), `deaths` (Bd deaths per stage that summer), `Z`, plus
`extinct_year` (first year with no hosts at the end of summer, or `nothing`), `persisted`
(hosts and Bd both present at the end) and, if `record_year` is set, `within`: daily `Z`,
`count` and `mean_load` (days × 6 matrices) and `max_load` during that summer.
"""
struct BriggsStagedProblem
    p::BriggsLifeStages
    years::Int
    invade_year::Int
    S0::Float64
    Δt::Float64
    poisson_below::Float64
    init::NTuple{6,Int}
    Z0::Float64
    record_year::Union{Nothing,Int}
end

function BriggsStagedProblem(p::BriggsLifeStages; years::Integer = 50,
        invade_year::Integer = 20, S0 = 10.0, Δt = 1.0, poisson_below = 100.0,
        init = disease_free_stages(p), Z0 = 0.0, record_year = nothing)
    Δt > 0 || throw(ArgumentError("Δt must be positive"))
    n0 = Tuple(round(Int, init[s]) for s in BRIGGS_STAGE_NAMES)
    all(>=(0), n0) || throw(ArgumentError("initial numbers must be non-negative"))
    return BriggsStagedProblem(p, years, invade_year, S0, Δt, poisson_below, n0, Z0,
        record_year === nothing ? nothing : Int(record_year))
end

# One Δt step for one stage (SI Eq. S1); returns the number of Bd deaths (loads set to NaN).
function briggs_life_stage_step!(rng, S, st::BriggsStage, p::BriggsLifeStages, Z, Δt, below)
    c1 = p.η * st.ν * st.f - p.σ
    grow = exp(c1 * Δt)
    gain = st.ν * st.γ / p.V * Z * growth_integral(c1, Δt)
    deaths = 0
    @inbounds for i in eachindex(S)
        s = maybe_poisson(rng, grow * S[i] + gain, below)
        if s >= st.Smax
            if st.lethal
                s = NaN
                deaths += 1
            else
                s = st.Smax
            end
        end
        S[i] = s
    end
    return deaths
end

function briggs_staged_summer!(rng, hosts, Z, prob::BriggsStagedProblem, rec)
    p, Δt, below = prob.p, prob.Δt, prob.poisson_below
    stages = life_stages(p)
    nsteps = round(Int, p.summer / Δt)
    deaths = zeros(Int, 6)
    t = 0.0
    nextday = 1.0
    day = 0
    for _ in 1:nsteps
        c2 = 0.0
        c3 = p.μ
        for k in 1:6
            c2 += p.η * (1 - stages[k].f) * sum(hosts[k]; init = 0.0)
            c3 += stages[k].γ / p.V * length(hosts[k])
        end
        Znew = exp(-c3 * Δt) * Z + c2 * growth_integral(-c3, Δt)
        for k in 1:6
            d = briggs_life_stage_step!(rng, hosts[k], stages[k], p, Z, Δt, below)
            if d > 0
                deaths[k] += d
                filter!(!isnan, hosts[k])
            end
        end
        Z = maybe_poisson(rng, Znew, below)
        t += Δt
        while t >= nextday - 1e-9
            p.εZ > 0 && (Z += rand(rng, Poisson(p.εZ)))
            nextday += 1.0
            day += 1
            if rec !== nothing && day <= size(rec.count, 1)
                rec.Z[day] = Z
                for k in 1:6
                    n = length(hosts[k])
                    rec.count[day, k] = n
                    rec.mean_load[day, k] = n == 0 ? 0.0 : sum(hosts[k]) / n
                    rec.max_load[day, k] = n == 0 ? 0.0 : maximum(hosts[k])
                end
            end
        end
    end
    return Z, deaths
end

function briggs_staged_spring!(rng, hosts, Z, p::BriggsLifeStages)
    tad1, tad2, tad3, sub1, sub2, adult = hosts
    stages = life_stages(p)
    # 1. winter survival
    for k in 1:6
        θ = stages[k].θ
        filter!(_ -> rand(rng) < θ, hosts[k])
    end
    Z = Float64(rand(rng, Binomial(round(Int, Z), p.θZ)))
    # 2. reproduction by the surviving adults
    breeders = rand(rng, Binomial(length(adult), p.pF))
    newborns = rand(rng, Poisson(p.R * breeders))
    # 3. maturation (loads carried; metamorphs lose their infection)
    Random.shuffle!(rng, sub2)
    open = max(p.K - length(adult), 0)
    append!(adult, view(sub2, 1:min(open, length(sub2))))
    empty!(sub2)
    append!(sub2, sub1)
    empty!(sub1)
    metamorphs = rand(rng, Binomial(length(tad3), p.ω))
    empty!(tad3)
    for s in tad2
        if rand(rng) < p.m
            rand(rng) < p.ω && (metamorphs += 1)
        else
            push!(tad3, s)
        end
    end
    append!(sub1, zeros(metamorphs))
    empty!(tad2)
    append!(tad2, tad1)
    empty!(tad1)
    # 4. newborn tadpoles
    append!(tad1, zeros(newborns))
    return Z
end

function CommonSolve.solve(prob::BriggsStagedProblem; rng::AbstractRNG = Random.default_rng())
    p = prob.p
    hosts = ntuple(k -> zeros(prob.init[k]), 6)
    Z = prob.Z0
    n = prob.years
    mk() = NamedTuple{BRIGGS_STAGE_NAMES}(ntuple(_ -> zeros(n), 6))
    counts = NamedTuple{BRIGGS_STAGE_NAMES}(ntuple(_ -> zeros(Int, n), 6))
    deaths = NamedTuple{BRIGGS_STAGE_NAMES}(ntuple(_ -> zeros(Int, n), 6))
    mean_load = mk()
    Zs = zeros(n)
    ndays = round(Int, p.summer)
    rec = prob.record_year === nothing ? nothing :
          (Z = zeros(ndays), count = zeros(Int, ndays, 6), mean_load = zeros(ndays, 6),
           max_load = zeros(ndays, 6))
    extinct_year = nothing
    for y in 1:n
        if y == prob.invade_year
            k = findlast(!isempty, hosts)
            k === nothing || (hosts[k][rand(rng, eachindex(hosts[k]))] = prob.S0)
        end
        Z, d = briggs_staged_summer!(rng, hosts, Z, prob, y == prob.record_year ? rec : nothing)
        for k in 1:6
            counts[k][y] = length(hosts[k])
            mean_load[k][y] = isempty(hosts[k]) ? 0.0 : sum(hosts[k]) / length(hosts[k])
            deaths[k][y] = d[k]
        end
        Zs[y] = Z
        if extinct_year === nothing && all(isempty, hosts)
            extinct_year = y
        end
        Z = briggs_staged_spring!(rng, hosts, Z, p)
    end
    nhosts = sum(length, hosts)
    bd = Z + sum(h -> sum(h; init = 0.0), hosts)
    return (counts = counts, mean_load = mean_load, deaths = deaths, Z = Zs,
        extinct_year = extinct_year, persisted = nhosts > 0 && bd > 0, within = rec)
end
