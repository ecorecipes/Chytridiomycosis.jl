# Load observation layer, part 2: likelihoods for swab loads.
#
# Two-part (hurdle) model with left-censoring. A swab comes from an infected host with
# probability π (the infected fraction among swabbed hosts); an infected host's load Y has
# distribution D (lognormal, gamma, Weibull or negative binomial). A swab with load below the
# limit of detection L reads as a non-detect:
#
#   P(non-detect) = (1 − π) + π P(Y < L)        detected y ≥ L:  π f_D(y)
#
# With L = 0 this is the classic hurdle (two-part) model: zeros are the uninfected, positives
# follow D. With L > 0 the infected hosts below L are censored, not dropped. On the count
# scale (scale = :count) loads are rounded up to integers k = ⌈y⌉ (as Schrock et al. 2025 did
# for their negative binomial fits) and every family gives a probability P(K = k), so the
# negative binomial and the continuous families can be compared by AIC. On the continuous
# scale the densities are of y itself, so AIC compares the continuous families with each
# other only. The structure is standard; it is not taken from one paper. Schrock et al. (2025)
# compared lognormal, gamma, Weibull and exponential fits to the positive loads only (their
# fitdistrplus MLE): fit_positive_loads and compare_load_families reproduce that procedure.

const HURDLE_FAMILIES = (:lognormal, :gamma, :weibull, :negbinomial)
const POSITIVE_FAMILIES = (:lognormal, :gamma, :weibull, :exponential)

"""
    LoadHurdleModel(loads; family = :lognormal, lod = 0, scale = :continuous)

Likelihood of swab loads under a two-part (hurdle) model with left-censoring at a limit of
detection `lod`:

- each swab comes from an infected host with probability `π`;
- an infected host's load follows `family`: `:lognormal`, `:gamma`, `:weibull` or
  `:negbinomial`;
- a load below `lod` (and any zero) is a non-detect, with probability
  `(1 − π) + π P(Y < lod)`; a detected load `y` contributes `π f(y)`.

With `lod = 0` zeros are the uninfected and the model is the classic hurdle model. `loads`
are numbers on one basis (state it: ITS1 copies or ZE per swab) or a vector of
[`PathogenLoad`](@ref)s of one type.

`scale = :count` rounds loads up to integers and uses the probability of each integer for
every family, so that the negative binomial (which requires it) can be compared with the
continuous families by AIC. On the default continuous scale only the continuous families
are comparable with each other.

It implements the LogDensityProblems interface: the log density is the log-likelihood at
unconstrained parameters (a flat prior) in the order of [`parameter_names`](@ref). Fit it
with [`maximum_likelihood`](@ref); [`pointwise_loglikelihood`](@ref) gives one term per
swab. [`hurdle_distribution`](@ref) returns the fitted `π` and distribution.
"""
struct LoadHurdleModel
    y::Vector{Float64}
    family::Symbol
    lod::Float64
    scale::Symbol
end

function LoadHurdleModel(loads::AbstractVector{<:Real}; family::Symbol = :lognormal,
        lod::Real = 0.0, scale::Symbol = family === :negbinomial ? :count : :continuous)
    family in HURDLE_FAMILIES || throw(ArgumentError("family must be one of $HURDLE_FAMILIES"))
    scale in (:continuous, :count) || throw(ArgumentError("scale must be :continuous or :count"))
    family === :negbinomial && scale !== :count &&
        throw(ArgumentError("the negative binomial needs scale = :count"))
    lod >= 0 || throw(ArgumentError("lod must be non-negative"))
    all(y -> y >= 0, loads) || throw(ArgumentError("loads must be non-negative"))
    isempty(loads) && throw(ArgumentError("no loads"))
    return LoadHurdleModel(Float64.(loads), family, Float64(lod), scale)
end

function LoadHurdleModel(loads::AbstractVector{L}; lod::Union{Real,L} = 0.0,
        kwargs...) where {L<:PathogenLoad}
    l = lod isa PathogenLoad ? load_value(lod) : lod
    return LoadHurdleModel(load_value.(loads); lod = l, kwargs...)
end

"""
    parameter_names(m::LoadHurdleModel)

Names of the unconstrained parameters of `m`:

- `logit_prevalence`: logit of `π`;
- `:lognormal`: `meanlog` (mean of ln load) and `log_sdlog`;
- `:gamma`: `log_mean` and `log_shape`;
- `:weibull`: `log_scale` and `log_shape`;
- `:negbinomial`: `log_mean` and `log_size` (`log k`; variance = mean + mean²/k).
"""
function parameter_names(m::LoadHurdleModel)
    m.family === :lognormal && return (:logit_prevalence, :meanlog, :log_sdlog)
    m.family === :gamma && return (:logit_prevalence, :log_mean, :log_shape)
    m.family === :weibull && return (:logit_prevalence, :log_scale, :log_shape)
    return (:logit_prevalence, :log_mean, :log_size)
end

LogDensityProblems.dimension(::LoadHurdleModel) = 3
LogDensityProblems.capabilities(::Type{LoadHurdleModel}) = LogDensityProblems.LogDensityOrder{0}()

function _load_dist(family::Symbol, a, b)
    family === :lognormal && return LogNormal(a, exp(b))
    family === :gamma && (k = exp(b); return Gamma(k, exp(a) / k))
    family === :weibull && return Weibull(exp(b), exp(a))
    family === :exponential && return Exponential(exp(a))
    r, μ = exp(b), exp(a)
    return NegativeBinomial(r, r / (r + μ))
end

"""
    hurdle_distribution(m::LoadHurdleModel, x) -> (; prevalence, distribution)

The infected fraction `π` and the load distribution of infected hosts at unconstrained
parameters `x` (for example the `x` of [`maximum_likelihood`](@ref)).
"""
function hurdle_distribution(m::LoadHurdleModel, x::AbstractVector)
    _check_dim(m, x)
    return (; prevalence = 1 / (1 + exp(-x[1])), distribution = _load_dist(m.family, x[2], x[3]))
end

_check_dim(m, x) = length(x) == LogDensityProblems.dimension(m) ||
    throw(DimensionMismatch("expected $(LogDensityProblems.dimension(m)) parameters"))

_log1pexp(z) = z > 0 ? z + log1p(exp(-z)) : log1p(exp(z))

# log P(K < L) and log P(K = k) on the count scale, K = ⌈Y⌉.
_count_logcdf(d, L) = logcdf(d, L - 1)
_count_logpmf(d::DiscreteUnivariateDistribution, k) = logpdf(d, k)
_count_logpmf(d, k) = k <= 0 ? -Inf : logdiffcdf(d, k, k - 1)

# One swab's log-likelihood. lπ = log π, l1π = log(1 − π).
function _swab_term(m::LoadHurdleModel, d, lπ, l1π, y)
    if m.scale === :continuous
        if y <= 0 || y < m.lod
            m.lod > 0 || return l1π
            lF = logcdf(d, m.lod)
            return _logaddexp(l1π, lπ + lF)
        end
        return lπ + logpdf(d, y)
    end
    L = max(ceil(m.lod), 1.0)
    k = ceil(y)
    if k < L
        return _logaddexp(l1π, lπ + _count_logcdf(d, L))
    end
    return lπ + _count_logpmf(d, k)
end

_logaddexp(a, b) = a == -Inf ? b : b == -Inf ? a : max(a, b) + log1p(exp(-abs(a - b)))

function pointwise_loglikelihood(m::LoadHurdleModel, x::AbstractVector)
    _check_dim(m, x)
    all(isfinite, x) || return fill(-Inf, length(m.y))
    d = _load_dist(m.family, x[2], x[3])
    lπ, l1π = -_log1pexp(-x[1]), -_log1pexp(x[1])
    return [_swab_term(m, d, lπ, l1π, y) for y in m.y]
end

function LogDensityProblems.logdensity(m::LoadHurdleModel, x::AbstractVector)
    ll = sum(pointwise_loglikelihood(m, x))
    return isnan(ll) ? -Inf : ll
end

# Moment-based start: π from the detected fraction, the family from the detected loads.
function _hurdle_start(m::LoadHurdleModel)
    det = m.scale === :count ? filter(y -> ceil(y) >= max(ceil(m.lod), 1.0), m.y) :
          filter(y -> y > 0 && y >= m.lod, m.y)
    p = clamp(length(det) / length(m.y), 0.05, 0.95)
    isempty(det) && return [log(p / (1 - p)), 0.0, 0.0]
    ly = log.(det)
    μl = sum(ly) / length(ly)
    sl = length(ly) > 1 ? sqrt(max(sum(abs2, ly .- μl) / (length(ly) - 1), 0.01)) : 1.0
    mn = sum(det) / length(det)
    v = length(det) > 1 ? max(sum(abs2, det .- mn) / (length(det) - 1), 1e-8) : mn^2
    a, b = m.family === :lognormal ? (μl, log(sl)) :
           m.family === :gamma ? (log(mn), log(clamp(mn^2 / v, 1e-3, 1e3))) :
           m.family === :weibull ? (μl + 0.5772 * sl * sqrt(6) / π, log(clamp(π / (sl * sqrt(6)), 1e-2, 1e2))) :
           (log(mn), log(clamp(mn^2 / max(v - mn, 1e-8), 1e-3, 1e3)))
    return [log(p / (1 - p)), a, b]
end

"""
    maximum_likelihood(m::LoadHurdleModel; starts = 4, iters = 2000, spread = 0.5,
        rng = Random.Xoshiro(1))

Maximum likelihood fit of a [`LoadHurdleModel`](@ref) by Nelder–Mead on the unconstrained
parameters, from a moment-based start and `starts − 1` jittered copies of it. Returns a
`NamedTuple`: `model`, `x`, `prevalence`, `distribution` (of infected hosts' loads),
`loglik`, `aic` (3 parameters) and `logliks` (best of each start).
"""
function maximum_likelihood(m::LoadHurdleModel; starts::Integer = 4, iters::Integer = 2000,
        spread::Real = 0.5, rng::AbstractRNG = Random.Xoshiro(1))
    f(x) = -LogDensityProblems.logdensity(m, x)
    x0 = _hurdle_start(m)
    best, bestf = x0, Inf
    logliks = Float64[]
    for s in 1:starts
        x = s == 1 ? copy(x0) : x0 .+ spread .* randn(rng, 3)
        x, fx = nelder_mead(f, x; iters = iters, step = 0.5)
        x, fx = nelder_mead(f, x; iters = iters ÷ 2, step = 0.05)
        push!(logliks, -fx)
        fx < bestf && ((best, bestf) = (x, fx))
    end
    ld = hurdle_distribution(m, best)
    return (; model = m, x = best, prevalence = ld.prevalence, distribution = ld.distribution,
        loglik = -bestf, aic = 6 + 2bestf, logliks)
end

"""
    fit_positive_loads(loads, family) -> (; family, distribution, loglik, aic)

Maximum likelihood fit of `family` (`:lognormal`, `:gamma`, `:weibull` or `:exponential`) to
the positive loads in `loads` (zeros dropped), as Schrock et al. (2025) fitted the
Bd-positive swabs of each group. Lognormal and exponential are closed-form; gamma and
Weibull are fitted by Nelder–Mead on log parameters. `aic = 2k − 2 log L`, with `k = 1` for
the exponential and 2 otherwise.
"""
function fit_positive_loads(loads::AbstractVector{<:Real}, family::Symbol)
    family in POSITIVE_FAMILIES || throw(ArgumentError("family must be one of $POSITIVE_FAMILIES"))
    y = Float64[v for v in loads if v > 0]
    length(y) >= 2 || throw(ArgumentError("need at least two positive loads"))
    ly = log.(y)
    n = length(y)
    if family === :lognormal
        μ = sum(ly) / n
        d = LogNormal(μ, sqrt(sum(abs2, ly .- μ) / n))
    elseif family === :exponential
        d = Exponential(sum(y) / n)
    else
        nll(x) = (dd = _load_dist(family, x[1], x[2]); s = -sum(logpdf(dd, v) for v in y);
                  isnan(s) ? Inf : s)
        μl = sum(ly) / n
        sl = max(sqrt(sum(abs2, ly .- μl) / (n - 1)), 0.05)
        x0 = family === :gamma ? [log(sum(y) / n), log(clamp(1 / sl^2, 1e-3, 1e3))] :
             [μl + 0.5772 * sl * sqrt(6) / π, log(π / (sl * sqrt(6)))]
        x, _ = nelder_mead(nll, x0; iters = 2000, step = 0.5)
        x, _ = nelder_mead(nll, x; iters = 1000, step = 0.05)
        d = _load_dist(family, x[1], x[2])
    end
    ll = sum(logpdf(d, v) for v in y)
    k = family === :exponential ? 1 : 2
    return (; family, distribution = d, loglik = ll, aic = 2k - 2ll)
end

"""
    compare_load_families(loads; families = (:lognormal, :gamma, :weibull, :exponential))
        -> Vector{NamedTuple}

Fit each family to the positive loads with [`fit_positive_loads`](@ref) and return the fits
sorted by AIC, each with `ΔAIC` from the best. Schrock et al. (2025) counted a family as
fitting a group as well as the best when `ΔAIC ≤ 2`.
"""
function compare_load_families(loads::AbstractVector{<:Real}; families = POSITIVE_FAMILIES)
    fits = [fit_positive_loads(loads, f) for f in families]
    sort!(fits; by = f -> f.aic)
    best = first(fits).aic
    return [merge(f, (; ΔAIC = f.aic - best)) for f in fits]
end
