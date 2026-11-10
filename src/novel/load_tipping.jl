# NEW MODEL (not from a published chytrid paper): the Briggs et al. (2010) within-host
# equation with one extra term, a saturating clearance of sporangia by constitutive skin
# defences. Wiki: lab/wiki/models/new-load-tipping.md
#
# One host, pool held fixed:
#   dS/dt = ι + (λ − σ) S − c S/(h + S)
#     ι = νγZ/V  sporangia per day from the pool
#     λ = ηνf    self-reinfection per sporangium per day
#     σ          sporangium loss per day
#     c          maximum clearance (sporangia per day)            NEW, no published value
#     h          load at which clearance is half saturated        NEW, no published value
#
# N identical hosts sharing a pool with an outside input of ε zoospores per day: the pool is
# fast, and at its quasi-steady state Z = (ε + η(1−f) N S)/(μ + γN/V) each host obeys the
# same equation with
#   ι → ι₀ = νγ ε / (V (μ + γN/V))
#   λ − σ → R = λ − σ + κ(N),   κ(N) = νγ η(1−f) N / (V (μ + γN/V))   ("pool feedback")
# The equilibria of this reduced equation are exactly those of the full (S, Z) system.
# κ(N) = σ − λ is the Briggs host-density threshold N*.
#
# Regimes, by R and ι:
#   R ≤ 0                      one stable load, no runaway          (:stable)
#   0 < R < c/h and ι ≤ ι_crit a stable low load and an unstable threshold (:bistable)
#   otherwise                  load grows without limit             (:runaway)
# with ι_crit = (√c − √(R h))².
# Units: days, sporangia per host.

"""
    LoadTipping(; ι, λ, σ, c, h, κ = 0)
    LoadTipping(p::BriggsParameters; c, h, Z)        # one host, pool fixed at Z
    LoadTipping(p::BriggsParameters; c, h, N, ε = 0) # N hosts, pool at quasi-steady state

Sporangial load with saturating clearance, `dS/dt = ι + (λ − σ + κ) S − c S/(h + S)`.
`c` (sporangia day⁻¹) and `h` (sporangia) are new parameters with no published values. `κ` is
the pool feedback from `N` identical hosts (see [`pool_feedback`](@ref)); it is zero for a
single host in a fixed pool.
"""
struct LoadTipping{T<:Real}
    ι::T
    λ::T
    σ::T
    c::T
    h::T
    κ::T
    function LoadTipping(ι, λ, σ, c, h, κ)
        ι, λ, σ, c, h, κ = promote(float(ι), float(λ), float(σ), float(c), float(h), float(κ))
        (ι >= 0 && λ >= 0 && σ >= 0 && c >= 0 && h > 0 && κ >= 0) ||
            throw(ArgumentError("need ι, λ, σ, c, κ ≥ 0 and h > 0"))
        return new{typeof(ι)}(ι, λ, σ, c, h, κ)
    end
end
LoadTipping(; ι, λ, σ, c, h, κ = 0.0) = LoadTipping(ι, λ, σ, c, h, κ)

"""
    pool_feedback(p::BriggsParameters, N) -> κ

Per-sporangium growth rate (day⁻¹) that `N` identical hosts gain through the shared pool,
`κ = νγη(1−f)N / (V(μ + γN/V))`. The Briggs threshold `N*` is where `κ = σ − ηνf`.
"""
pool_feedback(p::BriggsParameters, N::Real) =
    transmission_rate(p) * p.η * (1 - p.f) * N / (p.V * (p.μ + p.γ * N / p.V))

function LoadTipping(p::BriggsParameters; c, h, Z = nothing, N = nothing, ε = 0.0)
    λ = p.η * p.ν * p.f
    if Z !== nothing
        N === nothing || throw(ArgumentError("give either Z (one host) or N (a population)"))
        return LoadTipping(transmission_rate(p) * Z / p.V, λ, p.σ, c, h, 0.0)
    end
    N === nothing && throw(ArgumentError("give either Z (one host) or N (a population)"))
    ι0 = transmission_rate(p) * ε / (p.V * (p.μ + p.γ * N / p.V))
    return LoadTipping(ι0, λ, p.σ, c, h, pool_feedback(p, N))
end

"""Net per-sporangium growth rate at high load, `R = λ − σ + κ` (day⁻¹)."""
net_growth_rate(m::LoadTipping) = m.λ - m.σ + m.κ

"""
    critical_immigration(m::LoadTipping) -> ι_crit

Infection rate from the pool above which no low-load state exists: `(√c − √(R h))²` for
`0 ≤ R < c/h`, `0` for `R ≥ c/h`, and `Inf` for `R < 0` (loads never run away).
"""
function critical_immigration(m::LoadTipping)
    R = net_growth_rate(m)
    R < 0 && return oftype(R, Inf)
    return max(sqrt(m.c) - sqrt(R * m.h), zero(R))^2
end

"""
    tipping_regime(m::LoadTipping) -> :stable, :bistable or :runaway

`:stable`: one stable load and no runaway (`R ≤ 0`). `:bistable`: a stable low load and an
unstable threshold above which load runs away. `:runaway`: no low-load state.
"""
function tipping_regime(m::LoadTipping)
    R = net_growth_rate(m)
    R < 0 && return :stable
    R == 0 && return m.ι < m.c ? :stable : :runaway
    return (R * m.h < m.c && m.ι <= critical_immigration(m)) ? :bistable : :runaway
end

"""
    load_equilibria(m::LoadTipping) -> (low, threshold)

The stable low load and the unstable threshold load (sporangia per host): the roots of
`R S² + (ι + R h − c) S + ι h = 0`. In the `:stable` regime `threshold = Inf`; in the
`:runaway` regime both are `NaN`.
"""
function load_equilibria(m::LoadTipping)
    R, ι, c, h = net_growth_rate(m), m.ι, m.c, m.h
    regime = tipping_regime(m)
    regime === :runaway && return (low = oftype(R, NaN), threshold = oftype(R, NaN))
    R == 0 && return (low = ι * h / (c - ι), threshold = oftype(R, Inf))
    B = ι + R * h - c
    root = sqrt(max(B^2 - 4 * R * ι * h, zero(R)))
    if regime === :stable                      # R < 0: the single positive root
        low = iszero(ι) ? zero(R) : 2 * ι * h / (-B + root)
        return (low = low, threshold = oftype(R, Inf))
    end
    return (low = 2 * ι * h / (-B + root), threshold = (-B + root) / (2 * R))
end

"""
    critical_host_number(p::BriggsParameters; c, h, ε = 0) -> N

Smallest number of identical hosts per pool at which the low-load state is lost (the fold).
It is at or above the Briggs threshold `host_density_threshold(p)`, and `Inf` if clearance
holds at every density. Found by a scan and bisection on `N`.
"""
function critical_host_number(p::BriggsParameters; c, h, ε = 0.0, Nmax = 1e6)
    runaway(N) = tipping_regime(LoadTipping(p; c, h, N, ε)) === :runaway
    lo = max(host_density_threshold(p), 0.0)
    isfinite(lo) || return Inf
    runaway(lo) && return lo
    runaway(Nmax) || return Inf
    grid = exp.(range(log(max(lo, 1e-6)), log(Nmax); length = 4001))
    i = findfirst(runaway, grid)
    a, b = i == 1 ? lo : grid[i - 1], grid[i]
    for _ in 1:80
        mid = (a + b) / 2
        runaway(mid) ? (b = mid) : (a = mid)
    end
    return b
end

# ---------------------------------------------------------------------------------------
# Population simulation: N hosts, each with its own load, sharing one zoospore pool. This is
# the Briggs ODE plus the clearance term and the outside input ε; hosts die at Smax.

struct LoadTippingODEParams{P<:BriggsParameters}
    p::P
    hosts::BriggsHostState
    c::Vector{Float64}
    h::Vector{Float64}
    ε::Float64
end

function load_tipping_rhs!(du, u, q::LoadTippingODEParams, t)
    p, alive = q.p, q.hosts.alive
    N = length(u) - 1
    Z = u[end]
    infection = transmission_rate(p) / p.V * Z
    selfgrowth = p.η * p.ν * p.f - p.σ
    load = zero(eltype(u))
    nalive = 0
    @inbounds for i in 1:N
        if alive[i]
            S = max(u[i], zero(eltype(u)))
            du[i] = infection + selfgrowth * S - q.c[i] * S / (q.h[i] + S)
            load += S
            nalive += 1
        else
            du[i] = zero(eltype(u))
        end
    end
    du[end] = q.ε + p.η * (1 - p.f) * load - p.μ * Z - p.γ / p.V * Z * nalive
    return nothing
end

per_host(x::Real, N, name) = fill(float(x), N)
function per_host(x::AbstractVector, N, name)
    length(x) == N || throw(ArgumentError("length($name) must equal N"))
    return collect(float.(x))
end

"""
    pooled_load_equilibria(p::BriggsParameters; c, h, ε = 0, N = length(c))
        -> (regime, ι_low, low, ι_threshold, threshold)

Equilibria of `N` hosts that differ in their clearance parameters (`c`, `h`: a value or a
vector per host) and share one pool. Each host's load is set by the infection rate from the
pool, `ι`, and the pool closes the loop: `ι = ι₀ + (κ/N) Σᵢ Sᵢ(ι)`. `low` and `threshold` are
vectors of per-host loads at the stable low state and at the unstable threshold state
(empty when that state does not exist); `regime` is as in [`tipping_regime`](@ref). Needs
`ηνf < σ` (no host runs away on its own). Found by a scan and bisection on `ι`.
"""
function pooled_load_equilibria(p::BriggsParameters; c, h, ε = 0.0,
        N::Integer = max(length(c), length(h)))
    cs, hs = per_host(c, N, "c"), per_host(h, N, "h")
    λ = p.η * p.ν * p.f
    λ < p.σ || throw(ArgumentError("needs ηνf < σ: otherwise a host can run away on its own"))
    ι0 = transmission_rate(p) * ε / (p.V * (p.μ + p.γ * N / p.V))
    κ = pool_feedback(p, N)
    loads(ι) = [load_equilibria(LoadTipping(ι, λ, p.σ, cs[i], hs[i], 0.0)).low for i in 1:N]
    F(ι) = ι0 + κ / N * sum(loads(ι)) - ι            # > 0: the pool pushes ι up
    function refine(a, b)                            # root of F in (a, b), F changes sign
        fa = F(a)
        for _ in 1:100
            mid = (a + b) / 2
            (F(mid) > 0) == (fa > 0) ? (a = mid) : (b = mid)
        end
        return (a + b) / 2
    end
    grid = ι0 .+ exp10.(range(-9, 9; length = 3601))
    vals = F.(grid)
    none = Float64[]
    down = findfirst(<(0), vals)                     # first place the pool pulls ι down
    down === nothing && return (regime = :runaway, ι_low = NaN, low = none,
        ι_threshold = NaN, threshold = none)
    ιlow = down > 1 ? refine(grid[down - 1], grid[down]) :
           iszero(ι0) ? 0.0 : refine(ι0, grid[1])
    up = findnext(>(0), vals, down)
    up === nothing && return (regime = :stable, ι_low = ιlow, low = loads(ιlow),
        ι_threshold = Inf, threshold = none)
    ιthr = refine(grid[up - 1], grid[up])
    return (regime = :bistable, ι_low = ιlow, low = loads(ιlow), ι_threshold = ιthr,
        threshold = loads(ιthr))
end

"""
    LoadTippingProblem(p::BriggsParameters; c, h, N, ε = 0, S0 = zeros(N), Z0 = nothing,
                       tspan = (0.0, 365.0))

`N` hosts with loads `S0` sharing a zoospore pool that also receives `ε` zoospores per day.
Each host follows the Briggs equation plus clearance `c S/(h + S)`, and dies when its load
reaches `p.Smax`. `c` and `h` are either one value for all hosts or a vector with one value
per host (hosts that differ in their defences). `Z0` defaults to the pool's quasi-steady
state for `S0`. `solve` returns `(sol, death_time, alive)`, as for
[`BriggsLoadProblem`](@ref).
"""
struct LoadTippingProblem{P<:BriggsParameters,T}
    p::P
    c::Vector{Float64}
    h::Vector{Float64}
    ε::Float64
    S0::Vector{Float64}
    Z0::Float64
    tspan::T
end
function LoadTippingProblem(p::BriggsParameters; c, h, N::Integer, ε = 0.0, S0 = zeros(N),
        Z0 = nothing, tspan = (0.0, 365.0))
    length(S0) == N || throw(ArgumentError("length(S0) must equal N"))
    cs, hs = per_host(c, N, "c"), per_host(h, N, "h")
    (all(>=(0), cs) && all(>(0), hs) && ε >= 0) ||
        throw(ArgumentError("need c ≥ 0, h > 0, ε ≥ 0"))
    Z = Z0 === nothing ? (ε + p.η * (1 - p.f) * sum(S0)) / (p.μ + p.γ * N / p.V) : Z0
    return LoadTippingProblem(p, cs, hs, float(ε), collect(float.(S0)), float(Z),
        float.(tspan))
end

function to_ode_problem(prob::LoadTippingProblem)
    N = length(prob.S0)
    hosts = BriggsHostState(fill(true, N), fill(NaN, N))
    return ODEProblem(load_tipping_rhs!, vcat(prob.S0, prob.Z0), prob.tspan,
        LoadTippingODEParams(prob.p, hosts, prob.c, prob.h, prob.ε))
end

function CommonSolve.solve(prob::LoadTippingProblem, alg = Tsit5(); kwargs...)
    ode = to_ode_problem(prob)
    cb = briggs_death_callback(length(prob.S0))
    sol = solve(ode, alg; callback = cb, kwargs...)
    hosts = ode.p.hosts
    return (sol = sol, death_time = copy(hosts.death_time), alive = copy(hosts.alive))
end

# ---------------------------------------------------------------------------------------
# Stochastic form for one host in a fixed pool: the load BDI process (load_bdi.jl) with the
# extra, nonlinear death rate. Counting sporangia n:
#   up    b(n) = ι + λ n
#   down  d(n) = σ n + c n/(h + n)

tipping_up(m::LoadTipping, n) = m.ι + m.λ * n
tipping_down(m::LoadTipping, n) = m.σ * n + m.c * n / (m.h + n)

function require_single_host(m::LoadTipping)
    iszero(m.κ) || throw(ArgumentError(
        "the stochastic results are for one host in a fixed pool (κ = 0)"))
    return nothing
end

logaddexp(a, b) = max(a, b) + log1p(exp(-abs(a - b)))

"""
    mean_time_to_runaway(m::LoadTipping; Smax, from = 0) -> days

Mean time for the stochastic load of one host (pool fixed, `κ = 0`) to first reach `Smax`
sporangia, starting from `from`. Exact for the birth–death chain, by the recursion
`τ(n) = 1/b(n) + d(n)/b(n) τ(n−1)` for the mean time from `n` to `n + 1`. `Inf` when `ι = 0`
(the load can be lost for good). Its inverse is a hazard of death.
"""
function mean_time_to_runaway(m::LoadTipping; Smax::Integer, from::Integer = 0)
    require_single_host(m)
    0 <= from <= Smax || throw(ArgumentError("need 0 ≤ from ≤ Smax"))
    from == Smax && return 0.0
    m.ι > 0 || return Inf
    logτ = -log(tipping_up(m, 0))
    total = from == 0 ? logτ : -Inf
    for n in 1:(Smax - 1)
        b = tipping_up(m, n)
        logτ = logaddexp(-log(b), log(tipping_down(m, n)) - log(b) + logτ)
        n >= from && (total = logaddexp(total, logτ))
    end
    return exp(total)
end

"""
    low_state_distribution(m::LoadTipping; nmax) -> DiscreteNonParametric

Approximate distribution of sporangia on a host in the low-load state (pool fixed, `κ = 0`):
the stationary distribution of the birth–death chain confined to `0:nmax`. `nmax` defaults to
the unstable threshold in the `:bistable` regime and must be given otherwise. It is a good
approximation while runaway is rare. With `c = 0` it is the negative binomial of
[`LoadBDI`](@ref).
"""
function low_state_distribution(m::LoadTipping; nmax::Union{Nothing,Integer} = nothing)
    require_single_host(m)
    m.ι > 0 || throw(ArgumentError("ι = 0: the load stays at zero"))
    if nmax === nothing
        tipping_regime(m) === :bistable ||
            throw(ArgumentError("give nmax: there is no finite threshold in this regime"))
        nmax = floor(Int, load_equilibria(m).threshold)
    end
    logw = zeros(nmax + 1)
    for n in 1:nmax
        logw[n + 1] = logw[n] + log(tipping_up(m, n - 1)) - log(tipping_down(m, n))
    end
    w = exp.(logw .- maximum(logw))
    return DiscreteNonParametric(0:nmax, w ./ sum(w))
end

# ---------------------------------------------------------------------------------------
# Recruitment and inheritance: the N-host model run year after year. Within a year the loads
# follow the ODE above (deterministic). Between years hosts survive with probability
# `survival`, and each survivor leaves Poisson(`recruits`) recruits, up to the ceiling K
# (stochastic). Hosts differ in c, and a recruit's c depends on its one parent's:
#   log c_recruit = h² log c_parent + (1 − h²) log c_median + Normal(0, sdlog √(1 − h⁴))
# so that with no selection the population keeps median c_median and log-sd sdlog, and
# h² = `heritability` is the parent–offspring regression of log c. h² = 0 redraws every
# recruit from the founding distribution.
# Simplifications: recruits start with no load; survivors keep their loads and the pool keeps
# its zoospores from one year to the next; there are no stages and no winter.

"""
    LoadTippingRecruitmentProblem(p::BriggsParameters; c_median, sdlog, h, K, ε = 0,
        heritability = 0, recruits = 0.5, survival = 0.9, years = 50, season = 365.0,
        N0 = K)

The load tipping model over `years` years with host turnover. `N0` founders have clearance
`c` drawn from `LogNormal(log(c_median), sdlog)`; each year the loads run for `season` days,
then hosts survive with probability `survival` and each leaves `Poisson(recruits)` recruits,
up to `K` hosts. A recruit's `c` is inherited from one parent with parent–offspring
regression `heritability` (on the log scale).

`solve(prob; rng)` returns per-year vectors `hosts` (at the start of the year), `deaths` (from
Bd), `median_c` (at the start of the year) and `mean_load` (survivors, end of season), and
the final `c` of each host.
"""
struct LoadTippingRecruitmentProblem{P<:BriggsParameters}
    p::P
    c_median::Float64
    sdlog::Float64
    h::Float64
    K::Int
    ε::Float64
    heritability::Float64
    recruits::Float64
    survival::Float64
    years::Int
    season::Float64
    N0::Int
end
function LoadTippingRecruitmentProblem(p::BriggsParameters; c_median, sdlog, h, K::Integer,
        ε = 0.0, heritability = 0.0, recruits = 0.5, survival = 0.9, years::Integer = 50,
        season = 365.0, N0::Integer = K)
    (c_median > 0 && sdlog >= 0 && h > 0 && ε >= 0 && recruits >= 0 && season > 0) ||
        throw(ArgumentError("need c_median, h, season > 0 and sdlog, ε, recruits ≥ 0"))
    (0 <= heritability <= 1 && 0 <= survival <= 1) ||
        throw(ArgumentError("heritability and survival must be in [0, 1]"))
    0 < N0 <= K || throw(ArgumentError("need 0 < N0 ≤ K"))
    return LoadTippingRecruitmentProblem(p, float(c_median), float(sdlog), float(h), Int(K),
        float(ε), float(heritability), float(recruits), float(survival), Int(years),
        float(season), Int(N0))
end

function recruit_clearance(prob::LoadTippingRecruitmentProblem, c_parent, rng)
    h2 = prob.heritability
    noise = prob.sdlog * sqrt(1 - h2^2)
    return exp(h2 * log(c_parent) + (1 - h2) * log(prob.c_median) + noise * randn(rng))
end

function CommonSolve.solve(prob::LoadTippingRecruitmentProblem, alg = Tsit5();
        rng::AbstractRNG = Random.default_rng(), abstol = 1e-6, reltol = 1e-6, kwargs...)
    p = prob.p
    c = [exp(log(prob.c_median) + prob.sdlog * randn(rng)) for _ in 1:prob.N0]
    S = zeros(prob.N0)
    Z = 0.0
    out = (hosts = zeros(Int, prob.years), deaths = zeros(Int, prob.years),
        median_c = fill(NaN, prob.years), mean_load = fill(NaN, prob.years))
    for y in 1:prob.years
        n = length(c)
        n == 0 && break
        out.hosts[y] = n
        out.median_c[y] = sort(c)[(n + 1) ÷ 2]
        r = solve(LoadTippingProblem(p; c, h = prob.h, N = n, ε = prob.ε, S0 = S, Z0 = Z,
            tspan = (0.0, prob.season)), alg; abstol, reltol, save_everystep = false,
            kwargs...)
        SciMLBase.successful_retcode(r.sol) || error("load tipping: solver failed in year $y")
        out.deaths[y] = n - count(r.alive)
        final = r.sol.u[end]
        Z = max(final[end], 0.0)
        keep = [r.alive[i] && rand(rng) < prob.survival for i in 1:n]
        any(r.alive) && (out.mean_load[y] = sum(final[1:n][r.alive]) / count(r.alive))
        c, S = c[keep], max.(final[1:n][keep], 0.0)
        parents = length(c)
        parents == 0 && continue
        births = min(prob.K - parents, rand(rng, Poisson(prob.recruits * parents)))
        for _ in 1:births
            push!(c, recruit_clearance(prob, c[rand(rng, 1:parents)], rng))
            push!(S, 0.0)
        end
    end
    return merge(out, (c = c,))
end
