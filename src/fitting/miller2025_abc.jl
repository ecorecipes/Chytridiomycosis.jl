# ABC rejection for the Miller et al. (2025) refugia CTMC (§2.3; the authors' code
# simulate_samples.R, calc_summary_stats.R and plot_abc_results.R, Zenodo v1.0).
#
# For each prior draw θ = (β_sh, β_un, α, ω), with μ fixed, the authors simulate one shaded
# and one unshaded mesocosm from x0. The summary statistic (Eq. 2) is the sum of squared
# differences between the simulated counts and the counts of every mesocosm in the same
# environment, over the eight states and the survey weeks. The 0.2% of draws with the
# smallest statistic are accepted.
# - The code includes week 0, which the paper's Eq. (2) leaves out. The simulated week-0
#   state is always x0, so this adds the same constant to every draw and leaves the ranking
#   unchanged. `include_t0` controls it.
# - The paper's priors (Table 3) are β ~ U(0, 2) and ω ~ U(0, 2). The Zenodo code has
#   β ~ U(0, 1) and ω ~ U(0, 10). Both are available: RefugiaPrior(:miller2025; set =
#   "code_priors").

"""
    RefugiaPrior(; β_max, α_max, ω_max, μ_max = nothing)
    RefugiaPrior(:miller2025; set = "")

Independent uniform priors `β_sh, β_un ~ U(0, β_max)`, `α ~ U(0, α_max)` and
`ω ~ U(0, ω_max)`. With `μ_max === nothing` the loss rate is fixed, as in the paper. A
number makes it `μ ~ U(0, μ_max)` too, as in the coded `generate_priors` before the code
overwrites `μ`. That case is used for the μ–ω identifiability check.

The compendium constructor gives the paper's priors (Table 3). `set = "code_priors"` gives
the priors in the authors' Zenodo code.
"""
Base.@kwdef struct RefugiaPrior
    β_max::Float64
    α_max::Float64
    ω_max::Float64
    μ_max::Union{Nothing,Float64} = nothing
end

function RefugiaPrior(key::Symbol; set::AbstractString = "", μ_max = nothing)
    q = parameter_set(key; set)
    return RefugiaPrior(; β_max = q.β_prior_max, α_max = q.α_prior_max, ω_max = q.ω_prior_max,
                        μ_max)
end

# Squared error of one simulated environment (8 × weeks) against every mesocosm in it.
function _refugia_sqerr(d::RefugiaMesocosmData, shaded::Bool, sim::AbstractMatrix, k0::Int)
    s = 0
    @inbounds for m in eachindex(d.mesocosm)
        d.shaded[m] == shaded || continue
        for w in k0:length(d.weeks), i in 1:8
            s += (sim[i, w] - d.counts[i, w, m])^2
        end
    end
    return s
end

"""
    refugia_distance(data, sim_shaded, sim_unshaded; include_t0 = true)
        -> (; shaded, unshaded, total)

The ABC summary statistic of Miller et al. (2025), Eq. (2). Each simulation is an
8 × weeks matrix of counts at `data.weeks`, for example `solve(prob).u`. Each is compared
with every mesocosm of its own environment, and the squared differences are summed. With
`include_t0 = true` week 0 is included, as in the authors' `calc_summary_stats.R`.
"""
function refugia_distance(d::RefugiaMesocosmData, sim_sh::AbstractMatrix,
        sim_un::AbstractMatrix; include_t0::Bool = true)
    size(sim_sh) == size(sim_un) == (8, length(d.weeks)) ||
        throw(DimensionMismatch("simulations must be 8 × $(length(d.weeks))"))
    k0 = include_t0 ? 1 : 2
    sh = _refugia_sqerr(d, true, sim_sh, k0)
    un = _refugia_sqerr(d, false, sim_un, k0)
    return (; shaded = sh, unshaded = un, total = sh + un)
end

const REFUGIA_ABC_NAMES = (:β_sh, :β_un, :α, :ω, :μ)

"""
    refugia_abc(data = load_miller2025(); prior = RefugiaPrior(:miller2025),
                μ = parameter_set(:miller2025).μ, nsamples = 10^6, acceptance = 0.002,
                seed = 47, chunk = 10_000, include_t0 = true, base = RefugiaParameters(:miller2025))

ABC rejection as in Miller et al. (2025) §2.3. It runs on `Threads.nthreads()` threads, and
the result does not depend on the number of threads.

1. Draw `nsamples` parameter sets from `prior`.
2. For each, simulate one shaded and one unshaded mesocosm with [`RefugiaCTMCProblem`](@ref)
   dynamics (fixed constants from `base`, loss rate `μ` unless the prior varies it).
3. Score each by [`refugia_distance`](@ref) against `data`.
4. Accept the `round(acceptance × nsamples)` smallest. Ties are broken by draw order, as
   R's `order` does.

Each block of `chunk` draws has its own `Xoshiro` stream, seeded from `seed` and the block
number.

Returns `(; names, accepted, distance, threshold, nsamples)`:

- `accepted` is an `naccept × 5` matrix with columns `names = (:β_sh, :β_un, :α, :ω, :μ)`;
- `distance` is the accepted draws' statistics.

See [`refugia_abc_summary`](@ref).
"""
function refugia_abc(d::RefugiaMesocosmData = load_miller2025();
        prior::RefugiaPrior = RefugiaPrior(:miller2025), μ::Real = parameter_set(:miller2025).μ,
        nsamples::Integer = 10^6, acceptance::Real = 0.002, seed::Integer = 47,
        chunk::Integer = 10_000, include_t0::Bool = true,
        base::RefugiaParameters = RefugiaParameters(:miller2025))
    naccept = round(Int, acceptance * nsamples)
    1 <= naccept <= nsamples || throw(ArgumentError("acceptance × nsamples must be ≥ 1"))
    θ = zeros(nsamples, 5)
    dist = zeros(Int, nsamples)
    nchunks = cld(nsamples, chunk)
    tspan = (0.0, d.weeks[end])
    k0 = include_t0 ? 1 : 2
    tasks = map(1:nchunks) do c
        Threads.@spawn begin
            rng = Random.Xoshiro(UInt64(seed) << 32 + UInt64(c))
            out_sh = zeros(Int, 8, length(d.weeks))
            out_un = similar(out_sh)
            x = zeros(Int, 8)
            a = zeros(16)
            u0 = refugia_initial_state()
            for i in ((c - 1) * chunk + 1):min(c * chunk, nsamples)
                β_sh = prior.β_max * rand(rng)
                β_un = prior.β_max * rand(rng)
                α = prior.α_max * rand(rng)
                ω = prior.ω_max * rand(rng)
                μi = prior.μ_max === nothing ? Float64(μ) : prior.μ_max * rand(rng)
                p = RefugiaParameters(; β_sh, β_un, α, ω, μ = μi, γ1 = base.γ1, γ2 = base.γ2,
                                      m = base.m)
                _refugia_ssa!(out_sh, x, a, rng, p, β_sh, u0, tspan, d.weeks)
                _refugia_ssa!(out_un, x, a, rng, p, β_un, u0, tspan, d.weeks)
                θ[i, 1], θ[i, 2], θ[i, 3], θ[i, 4], θ[i, 5] = β_sh, β_un, α, ω, μi
                dist[i] = _refugia_sqerr(d, true, out_sh, k0) + _refugia_sqerr(d, false, out_un, k0)
            end
        end
    end
    foreach(wait, tasks)
    idx = sortperm(dist; alg = Base.Sort.DEFAULT_STABLE)[1:naccept]
    return (; names = REFUGIA_ABC_NAMES, accepted = θ[idx, :], distance = dist[idx],
            threshold = dist[idx[end]], nsamples = Int(nsamples))
end

"""
    refugia_abc_summary(res) -> NamedTuple

Posterior summaries of an ABC result from [`refugia_abc`](@ref), laid out like Table 3 of
Miller et al. (2025). For each of `β_sh`, `β_un`, `α`, `αβ_sh`, `αβ_un`, `ω` and `μ`, it
gives `(; mean, median, lo, hi)`, where `lo` and `hi` are the 2.5% and 97.5% quantiles
(R's default type 7).

It also gives the effects of §3.2, which are ratios of posterior means:

- `refugia_effect = 1 − mean(β_un)/mean(β_sh)` (Eq. 5);
- `vaccination_effect_sh = 1 − mean(αβ_sh)/mean(β_sh)` and `vaccination_effect_un`
  likewise (Eq. 4).
"""
function refugia_abc_summary(res)
    A = res.accepted
    cols = (β_sh = A[:, 1], β_un = A[:, 2], α = A[:, 3], αβ_sh = A[:, 3] .* A[:, 1],
            αβ_un = A[:, 3] .* A[:, 2], ω = A[:, 4], μ = A[:, 5])
    summ(v) = (; mean = sum(v) / length(v), median = median(v), lo = quantile(v, 0.025),
               hi = quantile(v, 0.975))
    s = map(summ, cols)
    effects = (; refugia_effect = 1 - s.β_un.mean / s.β_sh.mean,
               vaccination_effect_sh = 1 - s.αβ_sh.mean / s.β_sh.mean,
               vaccination_effect_un = 1 - s.αβ_un.mean / s.β_un.mean)
    return merge(s, effects)
end
