# Catch-curve estimate of annual adult survival before and after Bd, Heard et al. (2024)
# J Anim Ecol 93:1670-1683 (doi:10.1111/1365-2656.14177), following Age_structure_analysis.R
# in the authors' Zenodo archive (doi:10.5281/zenodo.10457229, CC BY 4.0). The likelihood,
# link and priors were checked against the paper (Methods 2.1, read 2026-11-20), which prints
# the posterior means 0.688, 0.352 (L. aurea pre/post-Bd), 0.655, 0.318 (L. raniformis).
#
# Data: skeletochronology ages (lines of arrested growth, LAG) of 403 frogs, Litoria aurea
# (55 historical museum specimens, 49 contemporary) and L. raniformis (45 historical, 254
# contemporary), from 15 sites. Following Scroggie (2012), age − 1 is geometric:
#
#   P(a_i) = (1 − S_i) S_i^{a_i},  a_i = LAG_i − 1
#   logit S_i = α + β₁ species_i + β₂ epoch_i + u_site(i),   u ~ Normal(0, σ_site)
#   α, β ~ Normal(0, precision 0.001), σ_site ~ Uniform(0, 10)
#
# species = 1 for L. raniformis, epoch = 1 for contemporary (post-Bd). The NIMBLE code uses
# the "ones trick" with C = 10,000, which adds the constant −n log C; we drop it. The derived
# survivals are S_h_La = logistic(α), S_c_La = logistic(α + β₂), S_h_Lr = logistic(α + β₁),
# S_c_Lr = logistic(α + β₁ + β₂).
#
# Sampling: the authors ran NIMBLE (3 chains × 50,000, 20,000 burn-in). We use the package's
# adaptive random-walk Metropolis sampler (src/fitting/wilber2022.jl) on the unconstrained
# parameters [α, β₁, β₂, logit(σ_site/10), u_1..u_15] with the Jacobian of the σ transform.

"""
    HeardAgeRecord

One frog of the Heard et al. (2024) age data: `species` (`:aurea` or `:raniformis`), `id`,
`year` of capture, `epoch` (`:historical` = museum specimen before Bd, `:contemporary`),
`site`, `region`, `sul` (snout–urostyle length, mm), `sex` (`'M'`/`'F'`), `age` (lines of
arrested growth) and `reliability` of the age reading (`"Clear"` or `"Satisfactory"`).
"""
struct HeardAgeRecord
    species::Symbol
    id::String
    year::Int
    epoch::Symbol
    site::String
    region::String
    sul::Float64
    sex::Char
    age::Int
    reliability::String
end

"""
    load_heard2024() -> Vector{HeardAgeRecord}

Read the bundled age-structure data of Heard et al. (2024) (`Age_structure_data.csv`, 403
frogs; see [`heard2024_path`](@ref)).
"""
function load_heard2024()
    lines = filter(!isempty ∘ strip, readlines(heard2024_path(:age_structure)))
    hdr = split(lines[1], ',')
    col(n) = (i = findfirst(==(n), hdr); i === nothing ? throw(ArgumentError("no column $n")) : i)
    isp, iid, iy, ie, isi, ir = col("Species"), col("ID"), col("Year of capture"),
        col("Epoch"), col("Site"), col("Region")
    isu, isx, ia, irl = col("SUL"), col("Sex"), col("LAG"), col("Reliability")
    out = HeardAgeRecord[]
    for l in lines[2:end]
        f = split(l, ',')
        sp = f[isp] == "Laur" ? :aurea : f[isp] == "Lran" ? :raniformis :
             throw(ArgumentError("unknown species $(f[isp])"))
        push!(out, HeardAgeRecord(sp, f[iid], parse(Int, f[iy]),
            Symbol(lowercase(f[ie])), f[isi], f[ir], parse(Float64, f[isu]),
            only(f[isx]), parse(Int, f[ia]), f[irl]))
    end
    return out
end

"""
    HeardCatchCurve(data = load_heard2024(); key = :heard2024)

Posterior of the Bayesian catch-curve model of Heard et al. (2024) (Scroggie 2012): the age
`a = LAG − 1` of each frog is geometric with annual survival `S`,

```
P(a) = (1 − S) S^a,   logit S = α + β₁ [L. raniformis] + β₂ [contemporary] + u_site
```

with site effects `u ~ Normal(0, σ_site)` and the priors of the authors' NIMBLE code
(`α, β ~ Normal(0, SD 31.6)`, `σ_site ~ Uniform(0, 10)`, read from the compendium).

It implements the `LogDensityProblems` interface on the unconstrained vector
`[α, β₁, β₂, logit(σ_site/σ_max), u_1, …, u_K]` (see [`parameter_names`](@ref)).
`solve(m; rng, chains = 3, warmup = 5000, iterations = 20_000, thin = 10)` samples it by
adaptive random-walk Metropolis and returns a [`HeardCatchCurveFit`](@ref).
"""
struct HeardCatchCurve
    age::Vector{Int}
    species::Vector{Int}
    epoch::Vector{Int}
    site::Vector{Int}
    sites::Vector{String}
    prior_sd::Float64
    sigma_max::Float64
end

function HeardCatchCurve(data::AbstractVector{HeardAgeRecord} = load_heard2024();
        key::Symbol = :heard2024)
    v = parameter_set(key)
    sites = unique(r.site for r in data)
    idx = Dict(s => i for (i, s) in enumerate(sites))
    return HeardCatchCurve([r.age - 1 for r in data],
        [r.species === :raniformis ? 1 : 0 for r in data],
        [r.epoch === :contemporary ? 1 : 0 for r in data],
        [idx[r.site] for r in data], sites, 1 / sqrt(v.prior_precision), v.sigma_site_max)
end

parameter_names(m::HeardCatchCurve) =
    vcat([:α, :β_species, :β_epoch, :logit_σ_site], [Symbol("u_", s) for s in m.sites])
LogDensityProblems.dimension(m::HeardCatchCurve) = 4 + length(m.sites)
LogDensityProblems.capabilities(::Type{HeardCatchCurve}) =
    LogDensityProblems.LogDensityOrder{0}()

_heard_log1pexp(x) = x > 0 ? x + log1p(exp(-x)) : log1p(exp(x))

"""
    heard_loglikelihood(m::HeardCatchCurve, θ) -> Float64

Geometric (catch-curve) log-likelihood `Σ log(1 − S_i) + a_i log S_i` at the unconstrained
parameters `θ` (no priors).
"""
function heard_loglikelihood(m::HeardCatchCurve, θ::AbstractVector)
    α, b1, b2 = θ[1], θ[2], θ[3]
    ll = 0.0
    for i in eachindex(m.age)
        η = α + b1 * m.species[i] + b2 * m.epoch[i] + θ[4+m.site[i]]
        # log S = −log(1 + e^{−η}),  log(1 − S) = −log(1 + e^{η})
        ll += -_heard_log1pexp(η) - m.age[i] * _heard_log1pexp(-η)
    end
    return ll
end

function LogDensityProblems.logdensity(m::HeardCatchCurve, θ::AbstractVector)
    length(θ) == LogDensityProblems.dimension(m) ||
        throw(DimensionMismatch("expected $(LogDensityProblems.dimension(m)) parameters"))
    all(isfinite, θ) || return -Inf
    z = θ[4]
    σ = m.sigma_max / (1 + exp(-z))
    σ > 0 || return -Inf
    lp = sum(x -> logpdf(Normal(0, m.prior_sd), x), view(θ, 1:3))
    # σ ~ Uniform(0, σmax), σ = σmax·logistic(z): density 1/σmax times |dσ/dz| =
    # σmax·logistic(z)(1 − logistic(z)), so log p(z) = −log1pexp(−z) − log1pexp(z)
    lp += -_heard_log1pexp(-z) - _heard_log1pexp(z)
    for k in eachindex(m.sites)
        lp += logpdf(Normal(0, σ), θ[4+k])
    end
    return lp + heard_loglikelihood(m, θ)
end

"""
    heard_survival(θ) -> NamedTuple

The four derived survivals of the catch-curve model at parameters `θ` (`α, β₁, β₂` first):
`S_h_La`, `S_c_La` (*L. aurea*, historical and contemporary), `S_h_Lr`, `S_c_Lr`
(*L. raniformis*).
"""
function heard_survival(θ::AbstractVector)
    L(x) = 1 / (1 + exp(-x))
    α, b1, b2 = θ[1], θ[2], θ[3]
    return (S_h_La = L(α), S_c_La = L(α + b2), S_h_Lr = L(α + b1), S_c_Lr = L(α + b1 + b2))
end

"""
    HeardCatchCurveFit

Posterior sample from `solve(::HeardCatchCurve)`: `model`, `names`, `draws` (one row per
draw, unconstrained parameters), `rhat` (split-R̂ per parameter) and `acceptance` (per
chain).
"""
struct HeardCatchCurveFit
    model::HeardCatchCurve
    names::Vector{Symbol}
    draws::Matrix{Float64}
    rhat::Vector{Float64}
    acceptance::Vector{Float64}
end

function CommonSolve.solve(m::HeardCatchCurve; rng::AbstractRNG = Random.default_rng(),
        chains::Integer = 3, warmup::Integer = 5000, iterations::Integer = 20_000,
        thin::Integer = 10)
    f = θ -> LogDensityProblems.logdensity(m, θ)
    d = LogDensityProblems.dimension(m)
    x0 = zeros(d)
    x0[1] = 0.5
    x0[4] = -2.0           # σ_site ≈ 1.2
    s = _w22_metropolis(f, x0, Matrix(0.02 * LinearAlgebra.I(d)), rng; chains, warmup,
        iterations, thin)
    return HeardCatchCurveFit(m, parameter_names(m), s.draws, s.rhat, s.acceptance)
end

"""
    heard_survival_summary(fit::HeardCatchCurveFit; level = 0.95)

Posterior `mean`, `median` and central `level` interval (`lower`, `upper`) of each derived
survival of [`heard_survival`](@ref), and of `σ_site`.
"""
function heard_survival_summary(fit::HeardCatchCurveFit; level::Real = 0.95)
    S = [heard_survival(view(fit.draws, i, :)) for i in axes(fit.draws, 1)]
    q(v) = (mean = mean(v), median = median(v), lower = quantile(v, (1 - level) / 2),
            upper = quantile(v, 1 - (1 - level) / 2))
    σ = [fit.model.sigma_max / (1 + exp(-x)) for x in fit.draws[:, 4]]
    return (S_h_La = q(getfield.(S, :S_h_La)), S_c_La = q(getfield.(S, :S_c_La)),
            S_h_Lr = q(getfield.(S, :S_h_Lr)), S_c_Lr = q(getfield.(S, :S_c_Lr)),
            σ_site = q(σ))
end

const _HEARD2024_FIT = Ref{Union{Nothing,HeardCatchCurveFit}}(nothing)

# One posterior sample per session with a fixed seed, for the validation targets.
function _heard2024_reference_fit()
    _HEARD2024_FIT[] === nothing || return _HEARD2024_FIT[]
    _HEARD2024_FIT[] = solve(HeardCatchCurve(); rng = Random.Xoshiro(2024))
    return _HEARD2024_FIT[]
end
