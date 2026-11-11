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
[`parameter_names`](@ref). Fit it with [`maximum_likelihood`](@ref). `dt` is the step of the
exponential integrator.
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
    out = similar(days, Float64)
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
        jl = [d < 0 ? 0.0 : knapp_jl(q, :treated, d) for d in days])
    return (days = days, treated = treated, wild = wild, control = control)
end

# log density of an observation y with log-median μ: lognormal, censored below 1
knapp_obs(y, μ, σ) = y >= 1 ? logpdf(Normal(μ, σ), log(y)) - log(y) : logcdf(Normal(μ, σ), 0.0)

function LogDensityProblems.logdensity(m::KnappLoadModel, x::AbstractVector)
    length(x) == LogDensityProblems.dimension(m) ||
        throw(DimensionMismatch("expected $(LogDensityProblems.dimension(m)) parameters"))
    all(isfinite, x) || return -Inf
    q = knapp_params(m, x)
    pred = predicted_loads(m, x)
    all(g -> all(isfinite, g.bd) && all(>(0), g.bd), (pred.treated, pred.wild, pred.control)) ||
        return -Inf
    ll = 0.0
    for s in m.swabs
        i = findfirst(==(s.day), pred.days)
        i === nothing && continue
        g = getfield(pred, s.group)
        isnan(s.bd) || (ll += knapp_obs(s.bd, log(g.bd[i]), q.σ))
        # J. lividum on treated frogs is modelled from the first bath on
        (!isnan(s.jl) && s.day > 0) && (ll += knapp_obs(s.jl, log(g.jl[i]), q.σJ))
    end
    return isfinite(ll) ? ll : -Inf
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
