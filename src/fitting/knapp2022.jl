# NEW MODEL FIT (not a published analysis): the frog load equation fitted to the Dusy Basin
# 2012 field swabs of Knapp et al. 2022 (doi:10.7717/peerj.12712), data in data/knapp2022.
# Wiki: lab/wiki/models/new-knapp-load-fit.md
#
# Each group of frogs is described by its median load L (Bd ITS copies per swab):
#   dL/dt = ι + r0 L − s(t) L/(h + L)
# ι is infection from outside the frog (copies per day), r0 = ηνf − σ the within-frog rate of
# the Briggs et al. 2010 model (from the compendium; negative, so a frog's own reinfection
# cannot sustain its load), and s(t) L/(h + L) extra clearance that saturates at high load.
# Four mechanisms for s(t) are compared:
#   :none      s = 0 for every group.
#   :decay     treated frogs only: s = Δc exp(−t/τ), suppression that wanes (as leftover
#              antifungal drug would).
#   :applied   treated frogs only: s = Δc J/(J + K), J the J. lividum on treated frogs.
#   :bacteria  every group: s = Δc J/(J + K) with its own J, so resident J. lividum on wild
#              and control frogs protects as much as the applied bacteria do.
# J. lividum (genomic equivalents per swab) declines exponentially from its day-12 level, with
# one half-life for treated frogs and one for wild and control frogs (held at the day-12 level
# before day 12). Untreated groups start at their equilibrium. Treated frogs start at a fitted
# load L0 on day 0 (the end of itraconazole); their day −7 swabs are predicted by the
# untreated equilibrium (with control J. lividum for :bacteria).
#
# Observations: log load ~ Normal(log L, σ), and a swab with no copies (or genomic
# equivalents) is censored below 1. Swabs are from group-marked frogs, so they are treated as
# independent; this ignores repeated sampling of the same animals.

"""
    KnappSwab(day, group, bd, jl)

One skin swab from Knapp et al. 2022: `day` relative to 19 July 2012 (the first *J. lividum*
bath), `group` (`:treated`, `:control` or `:wild`), Bd load `bd` (ITS copies per swab) and
*Janthinobacterium lividum* `jl` (genomic equivalents per swab); `NaN` when not measured.
"""
struct KnappSwab
    day::Float64
    group::Symbol
    bd::Float64
    jl::Float64
end

"""
    knapp2022_path()

Path of the bundled Dusy Basin 2012 swab file (Knapp et al. 2022, CC BY 4.0); see
`data/knapp2022/README.md`.
"""
knapp2022_path() = joinpath(pkgdir(@__MODULE__), "data", "knapp2022", "dusy-2012-swabs.csv")

"""
    load_knapp2022(path = knapp2022_path()) -> Vector{KnappSwab}

Read the Dusy Basin 2012 swabs of Knapp et al. 2022 (doi:10.7717/peerj.12712). Days are
counted from 19 July 2012, the first *J. lividum* bath: treated frogs received itraconazole
on days −7 to −1.
"""
function load_knapp2022(path::AbstractString = knapp2022_path())
    lines = readlines(path)
    header = split(lines[1], ',')
    function col(name)
        i = findfirst(==(name), header)
        i === nothing && throw(ArgumentError("column $name missing from $path"))
        return i
    end
    idate, iwild, itrt, ibd, ijl = col("visit_date"), col("expt_or_wild"), col("expt_trt"),
        col("bd_load"), col("jliv_ge")
    start = Date(2012, 7, 19)
    value(s) = s == "NA" ? NaN : parse(Float64, s)
    swabs = KnappSwab[]
    for line in lines[2:end]
        isempty(strip(line)) && continue
        r = split(line, ',')
        group = r[iwild] == "wild" ? :wild : Symbol(r[itrt])
        group in (:treated, :control, :wild) || throw(ArgumentError("unknown group $(r[itrt])"))
        day = (Date(r[idate], dateformat"yyyy-mm-dd") - start).value
        push!(swabs, KnappSwab(day, group, value(r[ibd]), value(r[ijl])))
    end
    return swabs
end

const KNAPP_MECHANISMS = (:none, :decay, :applied, :bacteria)
const KNAPP_COMMON = (:ι, :L0, :σ, :aT, :hT, :aW, :hW, :aC, :σJ)

"""
    KnappLoadModel(mechanism = :bacteria; swabs = load_knapp2022(),
        p = BriggsParameters(:briggs2010), dt = 0.05)

Likelihood of the Knapp et al. 2022 swabs under the frog load equation
`dL/dt = ι + r0 L − s(t) L/(h + L)` (loads in ITS copies per swab, time in days), with
`r0 = ηνf − σ` from `p`. `mechanism` sets the extra clearance `s(t)`:

- `:none`: none.
- `:decay`: treated frogs only, `Δc exp(−t/τ)`.
- `:applied`: treated frogs only, `Δc J/(J + K)`, with `J` the *J. lividum* on treated frogs.
- `:bacteria`: every group, `Δc J/(J + K)` with its own *J. lividum*.

The model implements the LogDensityProblems interface: the log density is the log-likelihood
of the parameters on the log scale (a flat prior), in the order of
[`parameter_names`](@ref). Fit it with [`maximum_likelihood`](@ref), or add a prior with
[`KnappPosterior`](@ref) for a Bayesian fit. `dt` is the step of the exponential integrator.
"""
struct KnappLoadModel
    swabs::Vector{KnappSwab}
    mechanism::Symbol
    r0::Float64
    dt::Float64
end
function KnappLoadModel(mechanism::Symbol = :bacteria; swabs = load_knapp2022(),
        p::BriggsParameters = BriggsParameters(:briggs2010), dt::Real = 0.05)
    mechanism in KNAPP_MECHANISMS ||
        throw(ArgumentError("mechanism must be one of $(KNAPP_MECHANISMS)"))
    r0 = p.η * p.ν * p.f - p.σ
    r0 < 0 || throw(ArgumentError("the fit needs r0 = ηνf − σ < 0; got $r0"))
    dt > 0 || throw(ArgumentError("dt must be positive"))
    return KnappLoadModel(collect(swabs), mechanism, r0, float(dt))
end

"""
    parameter_names(m::KnappLoadModel)

Names of the parameters, all fitted on the log scale:

- `ι`: infection from outside the frog (copies per swab per day);
- `L0`: load of treated frogs on day 0;
- `σ`: standard deviation of log Bd load;
- `aT`, `aW`, `aC`: *J. lividum* on day 12 on treated, wild and control frogs (genomic
  equivalents);
- `hT`, `hW`: its half-life (days) on treated frogs and on wild and control frogs;
- `σJ`: standard deviation of log *J. lividum*;
- `Δc`, `h`: the most extra clearance (copies per day) and the load at which it halves;
- `τ` (`:decay`): the time scale of waning (days);
- `K` (`:applied`, `:bacteria`): *J. lividum* giving half the extra clearance.
"""
parameter_names(m::KnappLoadModel) =
    m.mechanism === :none ? KNAPP_COMMON :
    m.mechanism === :decay ? (KNAPP_COMMON..., :Δc, :h, :τ) : (KNAPP_COMMON..., :Δc, :h, :K)

LogDensityProblems.dimension(m::KnappLoadModel) = length(parameter_names(m))
LogDensityProblems.capabilities(::Type{KnappLoadModel}) = LogDensityProblems.LogDensityOrder{0}()

knapp_params(m::KnappLoadModel, x) = NamedTuple{parameter_names(m)}(Tuple(exp.(x)))

# J. lividum (genomic equivalents per swab) on a group at day t
function knapp_jl(q, group::Symbol, t)
    group === :treated && return q.aT * 2.0^(-(t - 12) / q.hT)
    a = group === :wild ? q.aW : q.aC
    return a * 2.0^(-max(t - 12, 0.0) / q.hW)
end

# extra clearance on a group at day t (treated frogs before day 0 are untreated frogs)
function knapp_suppression(m::KnappLoadModel, q, group::Symbol, t)
    mech = m.mechanism
    mech === :none && return 0.0
    if group === :treated && t >= 0
        mech === :decay && return q.Δc * exp(-t / q.τ)
        J = knapp_jl(q, :treated, t)
        return q.Δc * J / (J + q.K)
    end
    mech === :bacteria || return 0.0
    J = knapp_jl(q, group === :treated ? :control : group, t)
    return q.Δc * J / (J + q.K)
end

# positive root of ι + r0 L − s L/(h + L) = 0 (unique when r0 < 0)
function knapp_equilibrium(ι, r0, s, h)
    s == 0 && return ι / -r0
    a, b, c = r0, ι + r0 * h - s, ι * h
    return (-b - sqrt(b^2 - 4a * c)) / (2a)
end

# one step of dL/dt = ι − k(t, L) L, k = −r0 + s/(h + L) ≥ 0, exactly for frozen k, with a
# midpoint correction; positive and stable for any dt
function knapp_step(L, ι, r0, s, h, dt)
    k = -r0 + s / (h + L)
    Lp = ι / k + (L - ι / k) * exp(-k * dt)
    k = -r0 + s / (h + (L + Lp) / 2)
    return ι / k + (L - ι / k) * exp(-k * dt)
end

# median load of `group` on each of `days` (sorted, ≥ t0), starting from L at t0
function knapp_path(m::KnappLoadModel, q, group::Symbol, L, t0, days)
    out = Vector{promote_type(typeof(L), typeof(q.ι))}(undef, length(days))
    dt, t, i = m.dt, float(t0), 1
    h = m.mechanism === :none ? 1.0 : q.h
    while i <= length(days)
        while t < days[i] - dt / 2
            L = knapp_step(L, q.ι, m.r0, knapp_suppression(m, q, group, t + dt / 2), h, dt)
            t += dt
        end
        out[i] = L
        i += 1
    end
    return out
end

const KNAPP_DAYS = [-7.0, 0.0, 12.0, 37.0, 56.0]

"""
    predicted_loads(m::KnappLoadModel, x; days = [-7, 0, 12, 37, 56])

Median Bd load (ITS copies per swab) and *J. lividum* (genomic equivalents) predicted for each
group on `days` (sorted), at log parameters `x`. Returns a `NamedTuple` with fields `days` and
`treated`, `wild`, `control`, each a `NamedTuple` of vectors `bd` and `jl`. Treated frogs before
day 0 are at their untreated equilibrium; their *J. lividum* is reported as 0 before day 0.
"""
function predicted_loads(m::KnappLoadModel, x::AbstractVector; days = KNAPP_DAYS)
    q = knapp_params(m, x)
    days = collect(Float64, days)
    issorted(days) || throw(ArgumentError("days must be sorted"))
    h = m.mechanism === :none ? 1.0 : q.h
    t0 = min(first(days), -7.0)
    group(g) = begin
        L = knapp_equilibrium(q.ι, m.r0, knapp_suppression(m, q, g, t0), h)
        (bd = knapp_path(m, q, g, L, t0, days), jl = [knapp_jl(q, g, d) for d in days])
    end
    wild, control = group(:wild), group(:control)
    pre = days .< 0
    Lpre = knapp_equilibrium(q.ι, m.r0, knapp_suppression(m, q, :treated, -1.0), h)
    post = knapp_path(m, q, :treated, q.L0, 0.0, days[.!pre])
    treated = (bd = vcat(fill(Lpre, count(pre)), post),
        jl = [d < 0 ? zero(q.aT) : knapp_jl(q, :treated, d) for d in days])
    return (days = days, treated = treated, wild = wild, control = control)
end

# log density of an observation y with log-median μ: lognormal, censored below 1
knapp_obs(y, μ, σ) = y >= 1 ? logpdf(Normal(μ, σ), log(y)) - log(y) : logcdf(Normal(μ, σ), 0.0)

function knapp_check_dimension(m::KnappLoadModel, x::AbstractVector)
    length(x) == LogDensityProblems.dimension(m) ||
        throw(DimensionMismatch("expected $(LogDensityProblems.dimension(m)) parameters"))
    return nothing
end

# a swab contributes when it was taken on an observation day and has a Bd load, or a
# J. lividum count after the first bath (J. lividum on treated frogs is modelled from then on)
knapp_bd_observed(s::KnappSwab) = s.day in KNAPP_DAYS && !isnan(s.bd)
knapp_jl_observed(s::KnappSwab) = s.day in KNAPP_DAYS && !isnan(s.jl) && s.day > 0

"""
    observed_swabs(m::KnappLoadModel)

Indices into `m.swabs` of the swabs that enter the likelihood, in the order of
[`pointwise_loglikelihood`](@ref): those taken on day −7, 0, 12, 37 or 56 with a Bd load, or
with a *J. lividum* count after the first bath.
"""
observed_swabs(m::KnappLoadModel) =
    findall(s -> knapp_bd_observed(s) || knapp_jl_observed(s), m.swabs)

"""
    pointwise_loglikelihood(m::KnappLoadModel, x::AbstractVector) -> Vector
    pointwise_loglikelihood(m::KnappLoadModel, X::AbstractMatrix) -> Matrix

Log-likelihood of each swab that enters the fit (Bd and *J. lividum* together), at log
parameters `x`, in the order of [`observed_swabs`](@ref). The sum is
`LogDensityProblems.logdensity(m, x)`. Every entry is `-Inf` where the model cannot be
evaluated. For a matrix `X` of draws (one draw of the log parameters per row) the result has
one row per draw and one column per swab, as [`psis_loo`](@ref) and [`waic`](@ref) take it.
"""
function pointwise_loglikelihood(m::KnappLoadModel, x::AbstractVector)
    knapp_check_dimension(m, x)
    T = float(eltype(x))
    keep = observed_swabs(m)
    out = fill(T(-Inf), length(keep))
    all(isfinite, x) || return out
    q = knapp_params(m, x)
    pred = predicted_loads(m, x)
    all(g -> all(isfinite, g.bd) && all(>(0), g.bd), (pred.treated, pred.wild, pred.control)) ||
        return out
    for (k, j) in enumerate(keep)
        s = m.swabs[j]
        i = findfirst(==(s.day), pred.days)
        g = getfield(pred, s.group)
        v = zero(T)
        knapp_bd_observed(s) && (v += knapp_obs(s.bd, log(g.bd[i]), q.σ))
        knapp_jl_observed(s) && (v += knapp_obs(s.jl, log(g.jl[i]), q.σJ))
        out[k] = isfinite(v) ? v : T(-Inf)
    end
    return out
end

function pointwise_loglikelihood(m::KnappLoadModel, X::AbstractMatrix)
    size(X, 2) == LogDensityProblems.dimension(m) ||
        throw(DimensionMismatch("X needs one column per parameter"))
    rows = [pointwise_loglikelihood(m, collect(r)) for r in eachrow(X)]
    return permutedims(reduce(hcat, rows))
end

function LogDensityProblems.logdensity(m::KnappLoadModel, x::AbstractVector)
    ll = sum(pointwise_loglikelihood(m, x))
    return isfinite(ll) ? ll : oftype(ll, -Inf)
end

const KNAPP_START = (ι = 3e5, L0 = 50.0, σ = 2.0, aT = 200.0, hT = 4.0, aW = 50.0, hW = 5.0,
    aC = 1500.0, σJ = 1.5, Δc = 1e5, h = 1e4, τ = 30.0, K = 50.0)

"""
    knapp_start(m::KnappLoadModel)

Default starting point (log scale) for [`maximum_likelihood`](@ref).
"""
knapp_start(m::KnappLoadModel) = [log(KNAPP_START[n]) for n in parameter_names(m)]

# Nelder–Mead minimisation of f from x0
function nelder_mead(f, x0::AbstractVector; iters::Integer = 2000, step::Real = 1.0)
    n = length(x0)
    xs = [float.(copy(x0)) for _ in 1:(n + 1)]
    for i in 1:n
        xs[i + 1][i] += step
    end
    fs = map(f, xs)
    for _ in 1:iters
        o = sortperm(fs)
        xs, fs = xs[o], fs[o]
        c = sum(xs[1:n]) / n
        xr = c + (c - xs[end]); fr = f(xr)
        if fr < fs[1]
            xe = c + 2(c - xs[end]); fe = f(xe)
            xs[end], fs[end] = fe < fr ? (xe, fe) : (xr, fr)
        elseif fr < fs[n]
            xs[end], fs[end] = xr, fr
        else
            xc = c + (xs[end] - c) / 2; fc = f(xc)
            if fc < fs[end]
                xs[end], fs[end] = xc, fc
            else
                for i in 2:(n + 1)
                    xs[i] = xs[1] + (xs[i] - xs[1]) / 2
                    fs[i] = f(xs[i])
                end
            end
        end
    end
    i = argmin(fs)
    return xs[i], fs[i]
end

"""
    maximum_likelihood(m::KnappLoadModel; starts = 8, iters = 2000, spread = 1.0,
        x0 = knapp_start(m), rng = Random.Xoshiro(1))

Maximum likelihood fit by Nelder–Mead on the log parameters, restarted from `x0` and from
`starts − 1` points drawn around it (normal, standard deviation `spread`; draws where the
likelihood cannot be evaluated are redrawn). Each start is polished by a second, smaller
simplex. Returns a `NamedTuple`: `model`, `x` (log parameters), `estimates` (natural scale,
by name), `loglik`, `aic`, and `logliks` (best of each start, to check agreement).
"""
function maximum_likelihood(m::KnappLoadModel; starts::Integer = 8, iters::Integer = 2000,
        spread::Real = 1.0, x0::AbstractVector = knapp_start(m),
        rng::AbstractRNG = Random.Xoshiro(1))
    f(x) = -LogDensityProblems.logdensity(m, x)
    best, bestf = copy(x0), Inf
    logliks = Float64[]
    for s in 1:starts
        x = copy(x0)
        if s > 1
            for _ in 1:1000
                x = x0 .+ spread .* randn(rng, length(x0))
                isfinite(f(x)) && break
            end
        end
        x, fx = nelder_mead(f, x; iters = iters)
        x, fx = nelder_mead(f, x; iters = iters ÷ 2, step = 0.2)
        push!(logliks, -fx)
        fx < bestf && ((best, bestf) = (x, fx))
    end
    k = length(best)
    return (model = m, x = best, estimates = knapp_params(m, best), loglik = -bestf,
        aic = 2k + 2bestf, logliks = logliks)
end

"""
    profile_likelihood(fit, name::Symbol, values; iters = 1500)

Profile log-likelihood of parameter `name` at each of `values` (natural scale), for a result of
[`maximum_likelihood`](@ref): the other parameters are re-optimised at each value, moving
outward from the estimate. Returns a vector of log-likelihoods. Values within about 1.92 of
`fit.loglik` form an approximate 95% confidence interval.
"""
function profile_likelihood(fit::NamedTuple, name::Symbol, values; iters::Integer = 1500)
    m = fit.model
    j = findfirst(==(name), parameter_names(m))
    j === nothing && throw(ArgumentError("$name is not a parameter of this model"))
    others = setdiff(1:length(fit.x), j)
    v = collect(Float64, values)
    out = fill(-Inf, length(v))
    centre = log(fit.estimates[name])
    below, above = fit.x[others], fit.x[others]   # warm starts on each side
    for i in sortperm(abs.(log.(v) .- centre))
        lv = log(v[i])
        f(y) = (x = Vector{Float64}(undef, length(fit.x)); x[j] = lv; x[others] = y;
                -LogDensityProblems.logdensity(m, x))
        y, fy = nelder_mead(f, lv < centre ? below : above; iters = iters, step = 0.3)
        y, fy = nelder_mead(f, y; iters = iters ÷ 2, step = 0.1)
        out[i] = -fy
        lv < centre ? (below = y) : (above = y)
    end
    return out
end

# ---------------------------------------------------------------------------------------------
# Bayesian fit. Priors are lognormal on the natural scale of each parameter, i.e. normal on the
# log parameters that the model works with. They are weakly informative: each is centred on an
# order of magnitude that is plausible before seeing these swabs and spans several orders of
# magnitude (scales) or a factor of about 20 (standard deviations and half-lives). r0 stays
# fixed at its compendium value, as in the maximum likelihood fit. Wiki:
# lab/wiki/models/new-knapp-load-fit.md, section "Bayesian fit".

# (median, standard deviation of the log) for each parameter
const KNAPP_PRIOR = (ι = (1e4, 3.0), L0 = (100.0, 3.0), σ = (2.0, 1.0), aT = (100.0, 3.0),
    hT = (7.0, 1.5), aW = (100.0, 3.0), hW = (7.0, 1.5), aC = (100.0, 3.0), σJ = (2.0, 1.0),
    Δc = (1e4, 3.0), h = (1e3, 3.0), τ = (14.0, 1.5), K = (10.0, 3.0))

"""
    knapp_prior(m::KnappLoadModel) -> NamedTuple

Default prior of the Bayesian fit: one `LogNormal(log(median), s)` for each parameter of
[`parameter_names`](@ref), on its natural scale (so normal on the log parameters). Medians
and `s` (95% of the prior within a factor `exp(1.96 s)` of the median):

| Parameter | Unit | Median | `s` | Reason |
|---|---|---|---|---|
| `ι` | ITS copies per swab per day | 10⁴ | 3 | equilibrium `ι/(−r0)` ≈ 6 × 10⁴, with 95% from 2 × 10² to 2 × 10⁷, the range of field swab loads |
| `L0` | ITS copies per swab | 100 | 3 | itraconazole lowers loads by orders of magnitude |
| `σ`, `σJ` | natural log | 2 | 1 | swab loads vary by factors of 10 to 100 |
| `aT`, `aW`, `aC` | genomic equivalents per swab | 100 | 3 | from undetectable to 10⁵ |
| `hT`, `hW` | days | 7 | 1.5 | from hours to months |
| `Δc` | ITS copies per swab per day | 10⁴ | 3 | the same scale as `ι` |
| `h` | ITS copies per swab | 10³ | 3 | anywhere in the range of observed loads |
| `τ` | days | 14 | 1.5 | drug residue lasting from a day to most of a year |
| `K` | genomic equivalents per swab | 10 | 3 | from a fraction of a cell to 10⁴ |

`r0` is fixed from the compendium (`p` of [`KnappLoadModel`](@ref)), as in the maximum
likelihood fit. These are new choices for this package, not published priors.
"""
knapp_prior(m::KnappLoadModel) = NamedTuple{parameter_names(m)}(
    Tuple(LogNormal(log(KNAPP_PRIOR[n][1]), KNAPP_PRIOR[n][2]) for n in parameter_names(m)))

"""
    KnappPosterior(m::KnappLoadModel; prior = knapp_prior(m))

Posterior of the log parameters of `m` (in the order of [`parameter_names`](@ref)) under
`prior`, a `NamedTuple` of `LogNormal` distributions on the natural scale of each parameter.
It implements the LogDensityProblems interface: the log density at log parameters `x` is the
log prior density of `x` (normal, `Normal(μ, s)` for `LogNormal(μ, s)`) plus
`LogDensityProblems.logdensity(m, x)`. Sample it with NUTS through
[`turing_model`](@ref) (`using Turing`), and compare mechanisms with
[`pointwise_loglikelihood`](@ref) and [`psis_loo`](@ref).
"""
struct KnappPosterior
    model::KnappLoadModel
    μ::Vector{Float64}
    s::Vector{Float64}
end
function KnappPosterior(m::KnappLoadModel; prior::NamedTuple = knapp_prior(m))
    names = parameter_names(m)
    Set(keys(prior)) == Set(names) ||
        throw(ArgumentError("prior must have one entry for each of $(names)"))
    all(d -> d isa LogNormal, values(prior)) ||
        throw(ArgumentError("each prior must be a LogNormal distribution"))
    return KnappPosterior(m, [prior[n].μ for n in names], [prior[n].σ for n in names])
end

parameter_names(post::KnappPosterior) = parameter_names(post.model)
LogDensityProblems.dimension(post::KnappPosterior) = length(post.μ)
LogDensityProblems.capabilities(::Type{KnappPosterior}) = LogDensityProblems.LogDensityOrder{0}()
pointwise_loglikelihood(post::KnappPosterior, x) = pointwise_loglikelihood(post.model, x)

function knapp_logprior(post::KnappPosterior, x::AbstractVector)
    knapp_check_dimension(post.model, x)
    return sum(logpdf(Normal(post.μ[i], post.s[i]), x[i]) for i in eachindex(x))
end

function LogDensityProblems.logdensity(post::KnappPosterior, x::AbstractVector)
    lp = knapp_logprior(post, x)
    isfinite(lp) || return oftype(lp, -Inf)
    return lp + LogDensityProblems.logdensity(post.model, x)
end

"""
    turing_model(post::KnappPosterior)

The posterior `post` as a Turing.jl model: the log parameters `x` (a vector in the order of
[`parameter_names`](@ref)) get the normal prior of `post`, and the log-likelihood of the
swabs is added with `@addlogprob!`. Its log joint at `x` equals
`LogDensityProblems.logdensity(post, x)`. Requires `using Turing`. The likelihood is
differentiable (a fixed-step integrator), so NUTS with ForwardDiff works:

```julia
using Turing
chn = sample(turing_model(post), NUTS(0.8), MCMCThreads(), 1000, 4)
X = Array(chn)      # draws × parameters, log scale
```
"""
function turing_model end

# ---------------------------------------------------------------------------------------------
# Comparing models by their predictive accuracy, from pointwise log-likelihoods of posterior
# draws: Pareto smoothed importance sampling leave-one-out cross-validation (PSIS-LOO; Vehtari,
# Gelman & Gabry 2017, doi:10.1007/s11222-016-9696-4; Vehtari et al. 2024,
# doi:10.48550/arXiv.1507.02646) and WAIC (Watanabe 2010). The PSIS steps follow the R
# package loo (tail of min(0.2S, 3√S) draws, Zhang & Stephens 2009 Pareto fit with a weak
# prior on k, weights truncated at the largest raw weight).

function logsumexp_(v)
    mx = maximum(v)
    isfinite(mx) || return mx
    return mx + log(sum(exp(a - mx) for a in v))
end
mean_(v) = sum(v) / length(v)
var_(v) = (mv = mean_(v); sum(abs2(a - mv) for a in v) / (length(v) - 1))

# Zhang & Stephens (2009) fit of a generalized Pareto distribution to sorted exceedances x
function gpd_fit(x::AbstractVector)
    n = length(x)
    m = 30 + floor(Int, sqrt(n))
    xstar = x[max(1, floor(Int, n / 4 + 0.5))]
    θ = [1 / x[n] + (1 - sqrt(m / (j - 0.5))) / 3 / xstar for j in 1:m]
    function prof(t)
        k = mean_(log1p(-t * a) for a in x)
        return log(-t / k) - k - 1
    end
    lθ = [n * prof(t) for t in θ]
    w = exp.(lθ .- logsumexp_(lθ))
    θhat = sum(θ .* w)
    k = mean_(log1p(-θhat * a) for a in x)
    σ = -k / θhat
    k = (k * n + 10 * 0.5) / (n + 10)              # weakly informative prior on k
    return isnan(k) ? (Inf, σ) : (k, σ)
end

# Pareto smoothed log weights (normalised so the largest raw weight is 1) and k̂
function psis_smooth(logratios::AbstractVector)
    S = length(logratios)
    lw = logratios .- maximum(logratios)
    M = min(ceil(Int, 0.2S), ceil(Int, 3sqrt(S)))
    k = Inf
    if M >= 5 && S > M
        o = sortperm(lw)
        tail = lw[o[(S - M + 1):S]]
        if tail[end] - tail[1] > eps() / 100
            cutoff = lw[o[S - M]]
            kf, σf = gpd_fit(exp.(tail) .- exp(cutoff))
            k = kf
            if isfinite(kf) && σf > 0
                p = ((1:M) .- 0.5) ./ M
                lw[o[(S - M + 1):S]] = log.(σf .* expm1.(-kf .* log1p.(-p)) ./ kf .+ exp(cutoff))
            end
        end
    end
    lw .= min.(lw, 0.0)
    return lw, k
end

"""
    psis_loo(ll::AbstractMatrix) -> NamedTuple

Leave-one-out cross-validation by Pareto smoothed importance sampling (Vehtari, Gelman &
Gabry 2017), from a matrix of pointwise log-likelihoods with one row per posterior draw and
one column per observation (as from [`pointwise_loglikelihood`](@ref)). Returns `elpd` (the
expected log predictive density, summed over observations), its standard error `se`, the
effective number of parameters `p` and, per observation, `pointwise` (elpd) and `k` (the
Pareto shape: estimates are reliable when `k < 0.7`). Draws from several chains can be
stacked; the tail length assumes a relative efficiency of 1.
"""
function psis_loo(ll::AbstractMatrix)
    S, n = size(ll)
    S >= 10 || throw(ArgumentError("need at least 10 draws"))
    pointwise, ks, lpd = zeros(n), zeros(n), zeros(n)
    for i in 1:n
        col = ll[:, i]
        lw, ks[i] = psis_smooth(-col)
        pointwise[i] = logsumexp_(lw .+ col) - logsumexp_(lw)
        lpd[i] = logsumexp_(col) - log(S)
    end
    elpd = sum(pointwise)
    return (elpd = elpd, se = sqrt(n * var_(pointwise)), p = sum(lpd) - elpd,
        pointwise = pointwise, k = ks)
end

"""
    waic(ll::AbstractMatrix) -> NamedTuple

Widely applicable information criterion (Watanabe 2010) on the elpd scale, from pointwise
log-likelihoods (one row per draw, one column per observation): `elpd`, its standard error
`se`, the effective number of parameters `p` (summed posterior variances of the pointwise
log-likelihood) and `pointwise`. `−2 elpd` is WAIC on the deviance scale.
"""
function waic(ll::AbstractMatrix)
    S, n = size(ll)
    S >= 2 || throw(ArgumentError("need at least 2 draws"))
    lpd = [logsumexp_(view(ll, :, i)) - log(S) for i in 1:n]
    pv = [var_(view(ll, :, i)) for i in 1:n]
    pointwise = lpd .- pv
    return (elpd = sum(pointwise), se = sqrt(n * var_(pointwise)), p = sum(pv),
        pointwise = pointwise)
end

"""
    elpd_difference(a, b) -> NamedTuple

Difference in expected log predictive density `a.elpd − b.elpd` between two results of
[`psis_loo`](@ref) or [`waic`](@ref) on the same observations, with the standard error of
the difference from the paired pointwise values (`se = √(n var(a_i − b_i))`). A positive
`difference` favours `a`.
"""
function elpd_difference(a::NamedTuple, b::NamedTuple)
    length(a.pointwise) == length(b.pointwise) ||
        throw(DimensionMismatch("a and b must be on the same observations"))
    d = a.pointwise .- b.pointwise
    return (difference = sum(d), se = sqrt(length(d) * var_(d)))
end
