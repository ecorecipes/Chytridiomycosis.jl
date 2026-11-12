# Briggs, Knapp & Vredenburg (2010) PNAS 107:9695, stochastic seasonal model (SI Methods).
# Three of the paper's variants: unstructured host, external zoospore source (εZ > 0), and a
# long-lived tadpole stage that is infected but not killed. The full six-stage R. muscosa
# model of Fig. 4 is in briggs2010_stages.jl.
# Units: days, sporangia per host, zoospores per pool. Wiki: lab/wiki/models/bd-briggs2010-load.md
#
# Within a summer, each step of length Δt applies the exact solutions of the two linear
# equations with the other variable frozen at its start-of-step value (SI Eqs. S1, S2):
#   S_i(t+Δt) = e^{c1 Δt} S_i + (β/V) Z (e^{c1 Δt} − 1)/c1,     c1 = ηνf − σ, β = νγ
#   Z(t+Δt)   = e^{−c3 Δt} Z + c2 (1 − e^{−c3 Δt})/c3,          c2 = Σ η(1−f) S_i, c3 = Σγ/V + μ
# Small values are replaced by Poisson draws with that mean.
#
# Choices the paper leaves open (each is a keyword or documented below):
#   Δt ("short")                       -> 1 day
#   "small numbers" for Poisson draws  -> below 100
#   load of the first infected frog    -> 10 sporangia
#   initial tadpoles (tadpole variant) -> R K
#   spring order (tadpole variant)     -> survive winter, adults breed, tadpoles metamorphose

"""
    BriggsStage(; γ, ν, f, Smax, θ, lethal = true)

Stage-specific values: encounter rate `γ`, establishment probability `ν`, self-reinfection
fraction `f`, load limit `Smax`, overwinter survival `θ`. If `lethal`, a host dies at `Smax`;
otherwise its load is capped there.
"""
Base.@kwdef struct BriggsStage
    γ::Float64
    ν::Float64
    f::Float64
    Smax::Float64
    θ::Float64
    lethal::Bool = true
end

"""
    BriggsSeasonalParameters(; V, μ, η, σ, adult, tadpole = nothing, R, K, θZ, m = 0.0,
                             εZ = 0.0, summer = 100.0)
    BriggsSeasonalParameters(:briggs2010; set = "fig3b" | "fig3d" | "fig3f")

Parameters of the stochastic seasonal model. Without a `tadpole` stage, each frog leaves
Poisson(`R`) recruits each spring, up to the ceiling `K`. With one, each adult leaves
Poisson(`R`) tadpoles, and each tadpole metamorphoses with probability `m` per spring (losing
its infection), up to `K` adults. `εZ` is the mean number of external zoospores added per day.
"""
Base.@kwdef struct BriggsSeasonalParameters
    V::Float64
    μ::Float64
    η::Float64
    σ::Float64
    adult::BriggsStage
    tadpole::Union{Nothing,BriggsStage} = nothing
    R::Float64
    K::Int
    θZ::Float64
    m::Float64 = 0.0
    εZ::Float64 = 0.0
    summer::Float64 = 100.0
end

function BriggsSeasonalParameters(key::Symbol; set::AbstractString)
    p = parameter_set(key; set)
    adult = BriggsStage(; γ = p.γ, ν = p.ν, f = p.f, Smax = p.Smax, θ = p.θF)
    tadpole = haskey(p, :γ_tadpole) ?
        BriggsStage(; γ = p.γ_tadpole, ν = p.ν, f = p.f, Smax = p.Smax_tadpole,
            θ = p.θ_tadpole, lethal = false) : nothing
    return BriggsSeasonalParameters(; V = p.V, μ = p.μ, η = p.η, σ = p.σ, adult, tadpole,
        R = p.R, K = Int(p.K), θZ = p.θZ, m = get(p, :m, 0.0), εZ = get(p, :εZ, 0.0),
        summer = p.summer)
end

"""
    BriggsSeasonalProblem(p; years = 10, Δt = 1.0, poisson_below = 100.0, S0 = 10.0,
                          adults0 = p.K, tadpoles0 = round(Int, p.R * p.K), Z0 = 0.0)

Stochastic seasonal Briggs et al. (2010) model, started as in their Fig. 3: `adults0` frogs of
which one carries `S0` sporangia, and `Z0` zoospores. `solve(prob; rng)` returns a NamedTuple
with per-year end-of-summer summaries (`adults`, `tadpoles`, `infected_adults`,
`mean_adult_load`, `Z`, `deaths`) and `persisted` (frogs and Bd both present after the final
winter).
"""
struct BriggsSeasonalProblem
    p::BriggsSeasonalParameters
    years::Int
    Δt::Float64
    poisson_below::Float64
    S0::Float64
    adults0::Int
    tadpoles0::Int
    Z0::Float64
end
function BriggsSeasonalProblem(p::BriggsSeasonalParameters; years::Integer = 10, Δt = 1.0,
        poisson_below = 100.0, S0 = 10.0, adults0::Integer = p.K,
        tadpoles0::Integer = p.tadpole === nothing ? 0 : round(Int, p.R * p.K), Z0 = 0.0)
    Δt > 0 || throw(ArgumentError("Δt must be positive"))
    adults0 >= 1 || throw(ArgumentError("need at least one adult to carry the infection"))
    return BriggsSeasonalProblem(p, years, Δt, poisson_below, S0, adults0, tadpoles0, Z0)
end

# (e^{c Δt} − 1)/c, continuous through c = 0
growth_integral(c, Δt) = abs(c * Δt) < 1e-12 ? Δt : expm1(c * Δt) / c

maybe_poisson(rng, x, below) = (0 < x < below) ? Float64(rand(rng, Poisson(x))) : max(x, 0.0)

# One Δt step for the hosts of one stage; returns the number killed (loads set to NaN).
function briggs_stage_step!(rng, S, st::BriggsStage, p::BriggsSeasonalParameters, Z, Δt, below)
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

function briggs_summer!(rng, adults, tadpoles, Z, prob::BriggsSeasonalProblem)
    p, Δt, below = prob.p, prob.Δt, prob.poisson_below
    nsteps = round(Int, p.summer / Δt)
    deaths = 0
    t = 0.0
    nextday = 1.0
    for _ in 1:nsteps
        c2 = p.η * (1 - p.adult.f) * sum(adults; init = 0.0)
        c3 = p.adult.γ / p.V * length(adults) + p.μ
        if p.tadpole !== nothing
            c2 += p.η * (1 - p.tadpole.f) * sum(tadpoles; init = 0.0)
            c3 += p.tadpole.γ / p.V * length(tadpoles)
        end
        Znew = exp(-c3 * Δt) * Z + c2 * growth_integral(-c3, Δt)
        deaths += briggs_stage_step!(rng, adults, p.adult, p, Z, Δt, below)
        filter!(!isnan, adults)
        p.tadpole === nothing ||
            briggs_stage_step!(rng, tadpoles, p.tadpole, p, Z, Δt, below)
        Z = maybe_poisson(rng, Znew, below)
        t += Δt
        while t >= nextday - 1e-9          # external source arrives once a day
            p.εZ > 0 && (Z += rand(rng, Poisson(p.εZ)))
            nextday += 1.0
        end
    end
    return Z, deaths
end

function briggs_winter!(rng, adults, tadpoles, Z, p::BriggsSeasonalParameters)
    filter!(_ -> rand(rng) < p.adult.θ, adults)
    Z = Float64(rand(rng, Binomial(round(Int, Z), p.θZ)))
    if p.tadpole === nothing
        recruits = rand(rng, Poisson(p.R * length(adults)))      # sum of Poisson(R) per frog
        append!(adults, zeros(clamp(recruits, 0, max(p.K - length(adults), 0))))
    else
        filter!(_ -> rand(rng) < p.tadpole.θ, tadpoles)
        newtads = rand(rng, Poisson(p.R * length(adults)))
        nmeta = rand(rng, Binomial(length(tadpoles), p.m))
        # Metamorphs lose their infection, so only their number matters to the adults, but
        # the tadpoles that leave must be chosen at random with respect to load.
        for _ in 1:nmeta
            deleteat!(tadpoles, rand(rng, 1:length(tadpoles)))
        end
        append!(adults, zeros(clamp(nmeta, 0, max(p.K - length(adults), 0))))
        append!(tadpoles, zeros(newtads))
    end
    return Z
end

function CommonSolve.solve(prob::BriggsSeasonalProblem; rng::AbstractRNG = Random.default_rng())
    p = prob.p
    adults = zeros(prob.adults0)
    adults[1] = prob.S0
    tadpoles = zeros(prob.tadpoles0)
    Z = prob.Z0
    n = prob.years
    out = (adults = zeros(Int, n), tadpoles = zeros(Int, n), infected_adults = zeros(Int, n),
        mean_adult_load = zeros(n), Z = zeros(n), deaths = zeros(Int, n))
    for y in 1:n
        Z, deaths = briggs_summer!(rng, adults, tadpoles, Z, prob)
        out.adults[y] = length(adults)
        out.tadpoles[y] = length(tadpoles)
        out.infected_adults[y] = count(>(0), adults)
        out.mean_adult_load[y] = isempty(adults) ? 0.0 : sum(adults) / length(adults)
        out.Z[y] = Z
        out.deaths[y] = deaths
        Z = briggs_winter!(rng, adults, tadpoles, Z, p)
    end
    hosts = length(adults) + length(tadpoles)
    bd = Z + sum(adults; init = 0.0) + sum(tadpoles; init = 0.0)
    return merge(out, (persisted = hosts > 0 && bd > 0, hosts_final = hosts, bd_final = bd))
end

"""
    persistence_probability(prob; nruns = 100, rng) -> fraction

Fraction of `nruns` runs in which hosts and Bd are both present at the end (the quantity
mapped in Fig. 3 A, C and E of Briggs et al. 2010, for `years = 10`).
"""
function persistence_probability(prob::BriggsSeasonalProblem; nruns::Integer = 100,
        rng::AbstractRNG = Random.default_rng())
    return count(_ -> solve(prob; rng).persisted, 1:nruns) / nruns
end
