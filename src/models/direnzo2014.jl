# DiRenzo, Langhammer, Zamudio & Lips (2014) Fungal infection intensity and zoospore output of
# Atelopus zeteki, a potential acute chytrid supershedder. PLOS ONE 9:e93356,
# doi:10.1371/journal.pone.0093356. Read from the published PDF.
#
# The paper fits, separately for Bd infection intensity (swab qPCR, zoospore genomic
# equivalents) and zoospore output (zoospores per minute in a 15-min soak), "linear mixed
# models with a first order autoregressive correlation term to ln transformed response
# variables" (nlme), in eight inoculated frogs sampled from day 7, then every 3-4 days, and at
# death:
#
#   ln y_ij = β0 + b_i + r t_ij + ε_ij,   b_i ~ N(0, σ_id²),   cor(ε_ij, ε_ik) = ϕ^|j−k|
#
# so that the mean trajectory is y(t) = y0 e^{rt} with y0 = e^{β0}. The Methods print the
# trajectory as "dN/dt = y0 e^{rt}", which is the integrated form N = y0 e^{rt} mislabelled as a
# derivative.
#
# What is printed: the slopes only, as "daily percent increase" 35.4% (SE 0.05) and 13.1%
# (SE 0.04), with t50 = 7.018 and t46 = 3.164. These are 100 × r on the ln scale
# (r/SE = 7.08 and 3.28 match the t values; 100(e^r − 1) would not). The implied percent
# increase per day is 100(e^r − 1) = 42.5% and 14.0%. The intercepts, σ_id, σ and ϕ are not
# printed, and no time series were deposited, so the fitted trajectories cannot be drawn on
# the data scale. Prior exposure was tested as a fixed effect and dropped by AIC.
#
# Unpublished choices:
# - The random effect is a random intercept per frog (the paper says "linear mixed models"
#   without naming the random effects).
# - The AR(1) is on the order of the observations within a frog (nlme's corAR1() default);
#   `correlation = :time` gives the continuous-time version ϕ^|t_j − t_k| (corCAR1).
# - The paper does not say how zero counts (3 swabs, 17 soaks) were log-transformed; the log
#   density takes the already transformed responses.

const DIRENZO_RESPONSES = (:intensity, :output)

"""
    DiRenzoParameters(; r_intensity, r_output, se_intensity = NaN, se_output = NaN)
    DiRenzoParameters(:direnzo2014)

Fitted daily growth rates (ln scale, day⁻¹) of Bd infection intensity and zoospore output on
*Atelopus zeteki* (DiRenzo et al. 2014), with their standard errors. The keyed constructor
reads the compendium (`r_intensity = 0.354`, `r_output = 0.131`). These are the only
parameters of the fitted mixed models that the paper prints.
"""
Base.@kwdef struct DiRenzoParameters
    r_intensity::Float64
    r_output::Float64
    se_intensity::Float64 = NaN
    se_output::Float64 = NaN
end

function DiRenzoParameters(key::Symbol)
    recs = parameter_records(key)
    rec(s) = recs[findfirst(r -> r.symbol === s, recs)]
    se(s) = TOML.parsefile(joinpath(COMPENDIUM_DIR, string(key) * ".toml"))["parameter"][
        findfirst(r -> r.symbol === s, recs)]["uncertainty"]["se"]
    return DiRenzoParameters(; r_intensity = rec(:r_intensity).value,
        r_output = rec(:r_output).value, se_intensity = se(:r_intensity),
        se_output = se(:r_output))
end

time_unit(::Type{DiRenzoParameters}) = :day

function _direnzo_rate(p::DiRenzoParameters, response::Symbol)
    response in DIRENZO_RESPONSES ||
        throw(ArgumentError("response must be :intensity or :output; got :$response"))
    return response === :intensity ? (p.r_intensity, p.se_intensity) : (p.r_output, p.se_output)
end

"""
    direnzo_trajectory(p::DiRenzoParameters, t; y0, response = :intensity)

Mean trajectory `y0 exp(r t)` of the log-linear model of DiRenzo et al. (2014): infection
intensity (zoospore genomic equivalents per swab) or zoospore output (zoospores per minute),
`t` days after inoculation. `y0` (the back-transformed intercept) is not printed in the
paper and must be given.
"""
direnzo_trajectory(p::DiRenzoParameters, t; y0::Real, response::Symbol = :intensity) =
    y0 * exp(first(_direnzo_rate(p, response)) * t)

"""
    daily_percent_increase(p::DiRenzoParameters; response = :intensity, as_printed = true)

The paper's "daily percent increase". With `as_printed = true` it is `100 r` (35.4% and
13.1%), which is what the paper reports; with `as_printed = false` it is the actual percent
increase per day, `100(eʳ − 1)` (42.5% and 14.0%).
"""
function daily_percent_increase(p::DiRenzoParameters; response::Symbol = :intensity,
        as_printed::Bool = true)
    r = first(_direnzo_rate(p, response))
    return as_printed ? 100r : 100expm1(r)
end

"""
    direnzo_t_statistic(p::DiRenzoParameters; response = :intensity) -> Float64

`r / SE(r)`, to compare with the paper's t statistics (t50 = 7.018 for intensity, t46 = 3.164
for output).
"""
direnzo_t_statistic(p::DiRenzoParameters; response::Symbol = :intensity) =
    /(_direnzo_rate(p, response)...)

"""
    DiRenzoSeries(id, prior, t, y)

Repeated measurements on one frog for [`DiRenzoLoadModel`](@ref): `t` days post-inoculation
and `y` the **ln-transformed** response, in time order. `prior` is `true` for a frog
previously exposed to the attenuated strain.
"""
struct DiRenzoSeries
    id::String
    prior::Bool
    t::Vector{Float64}
    y::Vector{Float64}
    function DiRenzoSeries(id, prior, t, y)
        length(t) == length(y) || throw(ArgumentError("t and y must have the same length"))
        issorted(t) || throw(ArgumentError("t must be sorted"))
        return new(string(id), prior, collect(float.(t)), collect(float.(y)))
    end
end

"""
    DiRenzoLoadModel(series; prior_effect = false, correlation = :order)

The log-linear mixed model with AR(1) errors of DiRenzo et al. (2014) for one response,

```
ln y_ij = β0 + δ·prior_i + b_i + r t_ij + ε_ij
b_i ~ N(0, σ_id²),  cor(ε_ij, ε_ik) = ϕ^|j − k|
```

(`ϕ^|t_j − t_k|` with `correlation = :time`). It implements the `LogDensityProblems`
interface; the log density is the marginal (ML, not REML) log-likelihood of `series` at
`x = [β0, r, (δ,) log σ_id, log σ, atanh ϕ]`. `prior_effect = true` adds `δ`; the paper dropped
it by AIC. See [`maximum_likelihood`](@ref) and [`direnzo_covariance`](@ref).
"""
struct DiRenzoLoadModel
    series::Vector{DiRenzoSeries}
    prior_effect::Bool
    correlation::Symbol
    function DiRenzoLoadModel(series::AbstractVector{DiRenzoSeries}; prior_effect::Bool = false,
            correlation::Symbol = :order)
        correlation in (:order, :time) ||
            throw(ArgumentError("correlation must be :order or :time"))
        return new(collect(series), prior_effect, correlation)
    end
end

parameter_names(m::DiRenzoLoadModel) = m.prior_effect ?
    (:β0, :r, :δ, :logσ_id, :logσ, :atanhϕ) : (:β0, :r, :logσ_id, :logσ, :atanhϕ)

LogDensityProblems.dimension(m::DiRenzoLoadModel) = m.prior_effect ? 6 : 5
LogDensityProblems.capabilities(::Type{DiRenzoLoadModel}) = LogDensityProblems.LogDensityOrder{0}()

"""
    direnzo_covariance(t; σ_id, σ, ϕ, correlation = :order) -> Matrix

Covariance of the ln responses of one frog sampled at times `t`: `σ_id² 11ᵀ + σ² R`, with
`R_jk = ϕ^|j − k|` (`:order`) or `ϕ^|t_j − t_k|` (`:time`).
"""
function direnzo_covariance(t::AbstractVector; σ_id::Real, σ::Real, ϕ::Real,
        correlation::Symbol = :order)
    n = length(t)
    lag(j, k) = correlation === :order ? abs(j - k) : abs(t[j] - t[k])
    return [σ_id^2 + σ^2 * ϕ^lag(j, k) for j in 1:n, k in 1:n]
end

function LogDensityProblems.logdensity(m::DiRenzoLoadModel, x::AbstractVector)
    length(x) == LogDensityProblems.dimension(m) ||
        throw(DimensionMismatch("expected $(LogDensityProblems.dimension(m)) parameters"))
    all(isfinite, x) || return -Inf
    β0, r = x[1], x[2]
    δ = m.prior_effect ? x[3] : 0.0
    σ_id, σ, ϕ = exp(x[end - 2]), exp(x[end - 1]), tanh(x[end])
    ll = 0.0
    for s in m.series
        μ = β0 .+ δ * s.prior .+ r .* s.t
        Σ = direnzo_covariance(s.t; σ_id, σ, ϕ, correlation = m.correlation)
        ll += logpdf(MvNormal(μ, Symmetric(Σ)), s.y)
    end
    return isfinite(ll) ? ll : -Inf
end

"""
    maximum_likelihood(m::DiRenzoLoadModel; x0 = nothing, iters = 4000)

Maximum likelihood fit of a [`DiRenzoLoadModel`](@ref) by Nelder–Mead. Returns
`(model, x, estimates, loglik, aic)`; `estimates` are on the natural scale
(`σ_id`, `σ`, `ϕ`).
"""
function maximum_likelihood(m::DiRenzoLoadModel; x0 = nothing, iters::Integer = 4000)
    f(x) = -LogDensityProblems.logdensity(m, x)
    if x0 === nothing
        ts = reduce(vcat, (s.t for s in m.series))
        ys = reduce(vcat, (s.y for s in m.series))
        r0 = sum((ts .- mean(ts)) .* (ys .- mean(ys))) / sum(abs2, ts .- mean(ts))
        x0 = vcat(mean(ys) - r0 * mean(ts), r0, m.prior_effect ? [0.0] : Float64[],
            log(std(ys) / 2), log(std(ys) / 2), 0.0)
    end
    x, fx = nelder_mead(f, collect(float.(x0)); iters, step = 0.5)
    x, fx = nelder_mead(f, x; iters = iters ÷ 2, step = 0.05)
    nat = vcat(x[1:(end - 3)], exp(x[end - 2]), exp(x[end - 1]), tanh(x[end]))
    names = m.prior_effect ? (:β0, :r, :δ, :σ_id, :σ, :ϕ) : (:β0, :r, :σ_id, :σ, :ϕ)
    return (model = m, x = x, estimates = NamedTuple{names}(Tuple(nat)), loglik = -fx,
        aic = 2fx + 2length(x))
end
