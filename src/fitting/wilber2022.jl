# Refits for Wilber et al. (2022) J Anim Ecol (doi:10.1111/1365-2656.13823) from the authors'
# deposit (Zenodo doi:10.5281/zenodo.7114634; see src/models/wilber2022.jl and
# data/wilber2022/README.md). This reproduces published analyses; it is not a new model.
#
# The authors fitted every model in Stan (3 chains × 1000 post-warmup draws). There is no
# MCMC library among the package's dependencies, so posteriors here come from an adaptive
# random-walk Metropolis sampler (Haario et al. 2001; proposal covariance learned during
# warm-up, scale tuned to an acceptance rate near 0.234), several chains from overdispersed
# starts, with split-R̂ reported. The HMM has 13–15 parameters and a
# likelihood that costs microseconds, so long chains are cheap. LOOIC uses the package's
# PSIS-LOO (psis_loo) on the per-lake log-likelihoods, as loo::loo did on log_lik.
#
# Hypothesis 2 (eqn 2; analysisII_decline.R, betabinomial_decline.stan): beta-binomial
# regression of the decline from the maximum to the minimum count. Code: the decline count
# n_max − n_min ~ BetaBinomial(n_max, pφ, (1 − p)φ), logit p = β0 + β1 z, z the z-scored
# log10 of the maximum count (or density); β0 flat, β1 ~ Normal(0, 5), φ ~ half-Cauchy(0, 2).
# Eqn 2 prints the minimum count as the response with the decline probability; the code's
# orientation is the one implemented (p is the expected proportion lost).
#
# Hypothesis 3, enzootic part (analysisIII_infection_load.Rmd): two least-squares regressions
# on the 111 Yosemite lake-year pairs (location == 1): z(Bd_t) ~ z(Bd_{t−1}) + z(log10 N_t)
# and z(log10 N_t) ~ z(Bd_{t−1}). The translocation part needs
# data/archival/humphreys_translocation_2010.csv, which is not in the deposit.

# --- adaptive random-walk Metropolis ---------------------------------------------------------

function _w22_cov(X::AbstractMatrix)
    μ = sum(X; dims = 1) ./ size(X, 1)
    Z = X .- μ
    return (Z' * Z) ./ (size(X, 1) - 1)
end

function _w22_chol(Σ)
    d = size(Σ, 1)
    S = Symmetric(Matrix(Σ) + 1e-10 * LinearAlgebra.I(d))
    C = LinearAlgebra.cholesky(S; check = false)
    LinearAlgebra.issuccess(C) && return C.L
    return LinearAlgebra.cholesky(Symmetric(Diagonal(max.(diag(Σ), 1e-6)))).L
end

# Split-R̂ (Gelman et al. 2013) of one parameter from a draws × chains matrix.
function _w22_rhat(x::AbstractMatrix)
    n = size(x, 1) ÷ 2
    halves = hcat(x[1:n, :], x[n+1:2n, :])
    m = size(halves, 2)
    means = vec(sum(halves; dims = 1)) ./ n
    vars = [sum(abs2, halves[:, j] .- means[j]) / (n - 1) for j in 1:m]
    W = sum(vars) / m
    B = n * sum(abs2, means .- sum(means) / m) / (m - 1)
    return sqrt(((n - 1) / n * W + B / n) / W)
end

# Returns (draws = (chains·iterations/thin) × d, rhat, acceptance).
function _w22_metropolis(f, x0::AbstractVector, Σ0::AbstractMatrix, rng::AbstractRNG;
                         chains::Integer = 4, warmup::Integer = 4000,
                         iterations::Integer = 12_000, thin::Integer = 4)
    d = length(x0)
    keep = iterations ÷ thin
    out = Array{Float64,3}(undef, keep, d, chains)
    acc = zeros(chains)
    L0 = _w22_chol(Σ0)
    for c in 1:chains
        x = collect(Float64, x0) .+ L0 * randn(rng, d)
        fx = f(x)
        tries = 0
        while !isfinite(fx) && tries < 100
            x = collect(Float64, x0) .+ 0.1 .* (L0 * randn(rng, d))
            fx = f(x)
            tries += 1
        end
        isfinite(fx) || (x = collect(Float64, x0); fx = f(x))
        λ = 2.38^2 / d
        Σ = Matrix{Float64}(Σ0)
        L = _w22_chol(λ .* Σ)
        hist = Matrix{Float64}(undef, warmup, d)
        nacc = 0
        for it in 1:warmup
            y = x .+ L * randn(rng, d)
            fy = f(y)
            if log(rand(rng)) < fy - fx
                x, fx = y, fy
                nacc += 1
            end
            hist[it, :] .= x
            if it % 100 == 0
                λ *= exp(2 * (nacc / 100 - 0.234))
                nacc = 0
                it >= 1000 && (Σ = _w22_cov(view(hist, (it ÷ 2):it, :)))
                L = _w22_chol(λ .* Σ)
            end
        end
        nacc = 0
        for it in 1:iterations
            y = x .+ L * randn(rng, d)
            fy = f(y)
            if log(rand(rng)) < fy - fx
                x, fx = y, fy
                nacc += 1
            end
            it % thin == 0 && (out[it ÷ thin, :, c] .= x)
        end
        acc[c] = nacc / iterations
    end
    rhat = [_w22_rhat(out[:, j, :]) for j in 1:d]
    draws = reduce(vcat, [out[:, :, c] for c in 1:chains])
    return (draws = draws, rhat = rhat, acceptance = acc)
end

function _w22_start(f, x0)
    try
        mode, Σ = _mode(f, x0; iters = 60)
        all(isfinite, Σ) && LinearAlgebra.isposdef(Symmetric(Σ)) && return mode, Σ
        return mode, Matrix(0.05 * LinearAlgebra.I(length(x0)))
    catch
        return collect(Float64, x0), Matrix(0.05 * LinearAlgebra.I(length(x0)))
    end
end

_w22_quantiles(v, level) = (median = median(v), lower = quantile(v, (1 - level) / 2),
                            upper = quantile(v, 1 - (1 - level) / 2))

# --- Hypothesis 1: the outbreak-state HMM -----------------------------------------------------

"""
    OutbreakHMMFit

Posterior sample of an [`OutbreakHMMProblem`](@ref) from `solve`: `problem`, `names`,
`draws` (one row per draw, unconstrained parameters in the order of `names`), `loglik` (per
draw and lake, for [`psis_loo`](@ref)), `rhat` (split-R̂ per parameter) and `acceptance`
(per chain).
"""
struct OutbreakHMMFit
    problem::OutbreakHMMProblem
    names::Vector{Symbol}
    draws::Matrix{Float64}
    loglik::Matrix{Float64}
    rhat::Vector{Float64}
    acceptance::Vector{Float64}
end

# solve(prob::OutbreakHMMProblem; chains = 4, warmup = 5000, iterations = 12_000, thin = 4,
#       rng = Random.Xoshiro(2022)) -> OutbreakHMMFit
# Adaptive random-walk Metropolis: `chains` chains from overdispersed points around zero (the
# log density is not concave there, so no Laplace start), each with `warmup` adaptive
# iterations, then `iterations` of which every `thin`-th is kept: 12 000 draws by default (the
# authors kept 3 000 from Stan's NUTS). Documented in OutbreakHMMProblem.
function CommonSolve.solve(prob::OutbreakHMMProblem; chains::Integer = 4,
                           warmup::Integer = 5000, iterations::Integer = 12_000,
                           thin::Integer = 4, rng::AbstractRNG = Random.Xoshiro(2022))
    f = θ -> LogDensityProblems.logdensity(prob, θ)
    d = length(prob.names)
    s = _w22_metropolis(f, zeros(d), Matrix(0.05 * LinearAlgebra.I(d)), rng;
                        chains, warmup, iterations, thin)
    ll = Matrix{Float64}(undef, size(s.draws, 1), size(prob.data.obs, 1))
    for (k, r) in enumerate(eachrow(s.draws))
        ll[k, :] .= pointwise_loglikelihood(prob, r)
    end
    return OutbreakHMMFit(prob, copy(prob.names), s.draws, ll, s.rhat, s.acceptance)
end

psis_loo(fit::OutbreakHMMFit) = psis_loo(fit.loglik)

"""
    looic(fit) -> Float64

The leave-one-out information criterion `−2 elpd_loo` (PSIS-LOO, [`psis_loo`](@ref)) of an
[`OutbreakHMMFit`](@ref), with one observation per lake as in the authors' `log_lik`.
"""
looic(fit::OutbreakHMMFit) = -2 * psis_loo(fit).elpd

"""
    posterior_summary(fit::OutbreakHMMFit; level = 0.95) -> NamedTuple

Posterior `names`, `mean`, `median`, `lower` and `upper` (central `level` interval) and
`rhat` of every unconstrained parameter of an [`OutbreakHMMFit`](@ref). The `omega_*` and
`phi_*` slopes are the effect sizes of Fig. 3 (per sd of each covariate; tadpole presence
unscaled).
"""
function posterior_summary(fit::OutbreakHMMFit; level::Real = 0.95)
    q = [_w22_quantiles(c, level) for c in eachcol(fit.draws)]
    return (names = fit.names, mean = vec(sum(fit.draws; dims = 1)) ./ size(fit.draws, 1),
            median = [x.median for x in q], lower = [x.lower for x in q],
            upper = [x.upper for x in q], rhat = fit.rhat)
end

"""
    wilber2022_effect(fit::OutbreakHMMFit, name::Symbol; level = 0.95) -> NamedTuple

Posterior `median`, `lower` and `upper` of one parameter (for example `:phi_temperature`).
"""
function wilber2022_effect(fit::OutbreakHMMFit, name::Symbol; level::Real = 0.95)
    j = findfirst(==(name), fit.names)
    j === nothing && throw(ArgumentError("no parameter $name in this model"))
    return _w22_quantiles(view(fit.draws, :, j), level)
end

const _WILBER2022_FITS = Dict{Tuple{Float64,Symbol},OutbreakHMMFit}()

"""
    wilber2022_fit(; cutoff = 0.25, model = :full) -> OutbreakHMMFit

The posterior of one of the Table 1 models (abundance covariate) at prevalence cutoff
`cutoff`, from `solve` with the default settings and seed; cached per `(cutoff, model)`.
Used by the validation targets.
"""
function wilber2022_fit(; cutoff::Real = 0.25, model::Symbol = :full)
    return get!(_WILBER2022_FITS, (Float64(cutoff), model)) do
        solve(OutbreakHMMProblem(; cutoff, model))
    end
end

const WILBER2022_CUTOFFS = (0.25, 0.3333, 0.5)
_w22_rho_label(ρ) = ρ == 0.25 ? "0.25" : ρ == 0.3333 ? "0.33" : ρ == 0.5 ? "0.5" :
                    throw(ArgumentError("cutoff $ρ is not one of $WILBER2022_CUTOFFS"))

# Validation helper: largest |posterior median − Fig. 3 median| over the nine slopes of the
# full abundance model other than φ-tadpole (whose posterior is long-tailed; checked
# separately), against the digitised compendium values.
function _wilber2022_fig3_deviation(ρ::Real)
    fit = wilber2022_fit(; cutoff = ρ)
    pub = parameter_set(:wilber2022; set = "abundance_rho" * _w22_rho_label(ρ))
    names = (:omega_abundance, :omega_temperature, :omega_snow, :omega_tadpole,
             :omega_temp_snow, :phi_abundance, :phi_temperature, :phi_snow, :phi_temp_snow)
    return maximum(abs(wilber2022_effect(fit, n).median - getfield(pub, n)) for n in names)
end

# Stored LOOIC of the Table 1 models (data/wilber2022/looic.toml, written by
# scripts/wilber2022_looic.jl: four sampler seeds per model and cutoff). These are our refits;
# the file holds no deposit data. Returns (mean, sd, looic, max_pareto_k).
function _wilber2022_looic(model::Symbol, ρ::Real)
    t = TOML.parsefile(joinpath(pkgdir(@__MODULE__), "data", "wilber2022", "looic.toml"))
    r = t[string(model)]["rho"*_w22_rho_label(ρ)]
    return (mean = r["mean"], sd = r["sd"], looic = r["looic"], max_pareto_k = r["max_pareto_k"])
end

# Validation helper: does `test(effect)` hold for parameter `name` at all three cutoffs?
_wilber2022_all_cutoffs(test, name::Symbol) =
    all(ρ -> test(wilber2022_effect(wilber2022_fit(; cutoff = ρ), name)), WILBER2022_CUTOFFS)

# --- Hypothesis 2: magnitude of decline ---------------------------------------------------------

const WILBER2022_DECLINE_EXCLUDED = (10055, 10276, 10314, 10315, 10316, 12102, 21111)

"""
    Wilber2022DeclineProblem(; subadults = true, covariate = :abundance,
                             path = wilber2022_path())

The beta-binomial regression of Hypothesis 2 of Wilber et al. (2022) (eqn 2): lakes where the
maximum count came before the minimum and Bd prevalence peaked (≥ 0.25) in between, minus the
seven lakes the authors excluded (82 lakes with subadults, 52 adults only). The decline
`n_max − n_min ~ BetaBinomial(n_max, pφ, (1 − p)φ)`, `logit p = β0 + β1 z`, `z` the z-scored
`log10` of the maximum count (`covariate = :abundance`) or count per metre of shoreline
(`:density`). Parameters `[β0, β1, log φ]`; priors as in the Stan code (β0 flat,
β1 ~ Normal(0, 5), φ ~ half-Cauchy(0, 2)). Fields `start`, `stop`, `z`, `mean`, `sd`
(of the log10 covariate), `sites`. `solve(prob; chains = 4, warmup = 4000,
iterations = 40_000, thin = 10, rng)` samples the posterior by adaptive random-walk
Metropolis from the Laplace mode ([`Wilber2022DeclineFit`](@ref)).
"""
struct Wilber2022DeclineProblem
    start::Vector{Int}
    stop::Vector{Int}
    z::Vector{Float64}
    mean::Float64
    sd::Float64
    sites::Vector{Int}
    subadults::Bool
    covariate::Symbol
end

function Wilber2022DeclineProblem(; subadults::Bool = true, covariate::Symbol = :abundance,
                                  path::AbstractString = wilber2022_path())
    covariate in (:abundance, :density) ||
        throw(ArgumentError("covariate must be :abundance or :density"))
    file = subadults ? "analysisII_magnitudeofdecline_withsubadults.csv" :
           "analysisII_magnitudeofdecline_withoutsubadults.csv"
    hdr, rows = _w22_csv(joinpath(path, "data", "candidate_datasets", file))
    col(n) = findfirst(==(n), hdr)
    keep = filter(rows) do r
        v(n) = _w22num(r[col(n)])
        v("year_max") < v("year_min") && v("prev_overall_max") >= 0.25 &&
            v("year_prev_max") >= v("year_max") && v("year_prev_max") <= v("year_min") &&
            !(parse(Int, r[col("site_id")]) in WILBER2022_DECLINE_EXCLUDED) &&
            v("abund_max") > 0
    end
    start = [round(Int, _w22num(r[col("abund_max")])) for r in keep]
    stop = [round(Int, _w22num(r[col("abund_min")])) for r in keep]
    x = log10.([_w22num(r[col(covariate === :abundance ? "abund_max" : "density")]) for r in keep])
    μ = sum(x) / length(x)
    s = sqrt(sum(abs2, x .- μ) / (length(x) - 1))
    sites = [parse(Int, r[col("site_id")]) for r in keep]
    return Wilber2022DeclineProblem(start, stop, (x .- μ) ./ s, μ, s, sites, subadults, covariate)
end

LogDensityProblems.dimension(::Wilber2022DeclineProblem) = 3
LogDensityProblems.capabilities(::Type{Wilber2022DeclineProblem}) =
    LogDensityProblems.LogDensityOrder{0}()

function pointwise_loglikelihood(m::Wilber2022DeclineProblem, θ::AbstractVector)
    β0, β1, φ = θ[1], θ[2], exp(θ[3])
    return map(eachindex(m.start)) do i
        p = 1 / (1 + exp(-(β0 + β1 * m.z[i])))
        α = (p * φ, (1 - p) * φ)
        x = (m.start[i] - m.stop[i], m.stop[i])
        n = m.start[i]
        _lgamma(n + 1) + _lgamma(φ) - _lgamma(n + φ) +
            sum(_lgamma(x[j] + α[j]) - _lgamma(x[j] + 1) - _lgamma(α[j]) for j in 1:2)
    end
end

function LogDensityProblems.logdensity(m::Wilber2022DeclineProblem, θ::AbstractVector)
    all(isfinite, θ) || return -Inf
    φ = exp(θ[3])
    lp = logpdf(Normal(0, 5), θ[2]) + log(2) + logpdf(Cauchy(0, 2), φ) + θ[3]
    ll = sum(pointwise_loglikelihood(m, θ))
    return isfinite(ll) ? lp + ll : -Inf
end

"""
    Wilber2022DeclineFit

Posterior sample of a [`Wilber2022DeclineProblem`](@ref): `problem`, `draws` (rows
`[β0, β1, log φ]`), `rhat`, `acceptance`.
"""
struct Wilber2022DeclineFit
    problem::Wilber2022DeclineProblem
    draws::Matrix{Float64}
    rhat::Vector{Float64}
    acceptance::Vector{Float64}
end

# solve(prob::Wilber2022DeclineProblem; chains = 4, warmup = 4000, iterations = 40_000,
#       thin = 10, rng = Random.Xoshiro(2022)) -> Wilber2022DeclineFit
# Adaptive random-walk Metropolis from the Laplace mode. Documented in Wilber2022DeclineProblem.
function CommonSolve.solve(prob::Wilber2022DeclineProblem; chains::Integer = 4,
                           warmup::Integer = 4000, iterations::Integer = 40_000,
                           thin::Integer = 10, rng::AbstractRNG = Random.Xoshiro(2022))
    f = θ -> LogDensityProblems.logdensity(prob, θ)
    mode, Σ = _w22_start(f, [1.0, 0.0, 1.0])
    s = _w22_metropolis(f, mode, Σ, rng; chains, warmup, iterations, thin)
    return Wilber2022DeclineFit(prob, s.draws, s.rhat, s.acceptance)
end

"""
    decline_effect(fit::Wilber2022DeclineFit; level = 0.95) -> NamedTuple

Posterior `median`, `lower` and `upper` of the slope β1 (per sd of log10 abundance or
density): the "Density/abundance effect" intervals of Fig. 5.
"""
decline_effect(fit::Wilber2022DeclineFit; level::Real = 0.95) =
    _w22_quantiles(view(fit.draws, :, 2), level)

"""
    predicted_decline(fit::Wilber2022DeclineFit, x; level = 0.95) -> NamedTuple

Posterior `median`, `lower` and `upper` of the expected proportional decline `p` at a
pre-decline count (or density) `x`, on the natural scale (the paper: 7 frogs, 84%).
"""
function predicted_decline(fit::Wilber2022DeclineFit, x::Real; level::Real = 0.95)
    z = (log10(x) - fit.problem.mean) / fit.problem.sd
    p = [1 / (1 + exp(-(r[1] + r[2] * z))) for r in eachrow(fit.draws)]
    return _w22_quantiles(p, level)
end

const _WILBER2022_DECLINE = Dict{Tuple{Bool,Symbol},Wilber2022DeclineFit}()

"""
    wilber2022_decline(; subadults = true, covariate = :abundance) -> Wilber2022DeclineFit

The cached posterior of one of the four decline regressions of Fig. 5.
"""
wilber2022_decline(; subadults::Bool = true, covariate::Symbol = :abundance) =
    get!(_WILBER2022_DECLINE, (subadults, covariate)) do
        solve(Wilber2022DeclineProblem(; subadults, covariate))
    end

# --- Hypothesis 3: host density and load in enzootic populations ------------------------------

"""
    wilber2022_path_analysis(; covariate = :abundance, path = wilber2022_path()) -> NamedTuple

The enzootic path analysis of Wilber et al. (2022) (Hypothesis 3; Fig. S3) by least squares
on the Yosemite lake-year pairs of `analysisIII_loadanddensity.csv`: `host_to_load` is the
standardised coefficient of `log10` abundance (or density) in year `t` on mean log10 Bd load
in year `t` given the load in `t − 1`; `load_to_host` the coefficient of load in `t − 1` on
abundance in `t`. Each has `estimate`, `se`, and the 95% t interval `lower`, `upper`. Also
`n` (pairs) and `populations`.
"""
function wilber2022_path_analysis(; covariate::Symbol = :abundance,
                                  path::AbstractString = wilber2022_path())
    covariate in (:abundance, :density) ||
        throw(ArgumentError("covariate must be :abundance or :density"))
    hdr, rows = _w22_csv(joinpath(path, "data", "candidate_datasets", "analysisIII_loadanddensity.csv"))
    col(n) = findfirst(==(n), hdr)
    yose = filter(r -> _w22num(r[col("location")]) == 1, rows)
    bdp = [_w22num(r[col("mean_bd_present")]) for r in yose]
    bdq = [_w22num(r[col("mean_bd_past")]) for r in yose]
    host = log10.([_w22num(r[col(covariate === :abundance ? "frog_abund_present" : "density")])
                   for r in yose])
    zs(v) = (v .- sum(v) / length(v)) ./ sqrt(sum(abs2, v .- sum(v) / length(v)) / (length(v) - 1))
    function ols(y, X)
        n, k = size(X)
        β = X \ y
        s2 = sum(abs2, y .- X * β) / (n - k)
        se = sqrt.(diag(s2 .* inv(X' * X)))
        t = quantile(TDist(n - k), 0.975)
        return (estimate = β[end], se = se[end], lower = β[end] - t * se[end],
                upper = β[end] + t * se[end])
    end
    n = length(yose)
    h2l = ols(zs(bdp), hcat(ones(n), zs(bdq), zs(host)))
    l2h = ols(zs(host), hcat(ones(n), zs(bdq)))
    return (host_to_load = h2l, load_to_host = l2h, n = n,
            populations = length(unique(r[col("site_id")] for r in yose)))
end
