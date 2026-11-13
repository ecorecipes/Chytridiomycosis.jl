# Refit of Tompros et al. (2022) Transbound Emerg Dis 69:731-741 (doi:10.1111/tbed.14043) to
# the authors' TRACE deposit (doi:10.7290/xw8yy2RBeu; data/tompros2022). This reproduces a
# published analysis; it is not a new model. Wiki: lab/wiki/models/bsal-transmission-functions.md
#
# Likelihood (appendix eqs 9-10): for each tank and each interval between swab days, the
# probability of the observed number of susceptibles at the end of the interval given the
# numbers susceptible and infected at its start, under the pure-death master equation (eq. 2)
# with no deaths or recoveries. The interval likelihoods are multiplied (a Markov chain).
# Tanks: 20, swab days 0 (initial conditions), 2 and 5. AIC = 2 NLL + 2k.
#
# The authors integrate eq. 2 with deSolve (lsoda, partly on the log scale) and minimise with
# optim(L-BFGS-B) on the natural scale with every parameter ≥ 1e-6. Here the probabilities
# come from a matrix exponential and the minimisation is Nelder–Mead on the log scale, with
# the same lower bound.

const TOMPROS_FILES = Dict(:fit => "counts_days0to5.csv", :validation => "counts_day8.csv",
    :swabs => "swabs.csv")

"""
    tompros2022_path(name::Symbol)

Path of a bundled data file of Tompros et al. (2022) (TRACE deposit
doi:10.7290/xw8yy2RBeu; see `data/tompros2022/README.md`). `name` is one of:

- `:fit`: numbers susceptible and infected in each tank on days 0, 2 and 5
  (`counts_days0to5.csv`, the deposit's `DataFile2_Days0to5.RData`);
- `:validation`: the same on days 5 and 8 for the 10 tanks with no deaths by day 8
  (`counts_day8.csv`, the deposit's `DataFile3_Day8_NoDeaths.RData`);
- `:swabs`: the swab results of every newt (`swabs.csv`, columns of the deposit's
  `DataFile1.csv`).
"""
function tompros2022_path(name::Symbol)
    haskey(TOMPROS_FILES, name) ||
        throw(ArgumentError("name must be one of $(sort!(collect(keys(TOMPROS_FILES))))"))
    return joinpath(pkgdir(@__MODULE__), "data", "tompros2022", TOMPROS_FILES[name])
end

"""
    TomprosTank(tank, days, S, I)

The epidemic in one 1-m² mesocosm of Tompros et al. (2022): swab days `days` and the numbers
of susceptible (`S`) and infected (`I`) newts on those days, after the authors'
false-negative and false-positive rules.
"""
struct TomprosTank
    tank::Int
    days::Vector{Float64}
    S::Vector{Int}
    I::Vector{Int}
end

function _tanks_from_rows(rows)
    out = TomprosTank[]
    for k in sort!(unique(r.tank for r in rows))
        rs = sort!(filter(r -> r.tank == k, rows); by = r -> r.day)
        push!(out, TomprosTank(k, [r.day for r in rs], [r.S for r in rs], [r.I for r in rs]))
    end
    return out
end

"""
    load_tompros2022(; data = :fit) -> Vector{TomprosTank}

The tank epidemics of Tompros et al. (2022). `data = :fit` gives the 20 tanks on days 0, 2
and 5, which the paper fits. `data = :validation` gives the 10 tanks with no deaths or
recoveries by day 8, on days 5 and 8, which the paper holds out (section 2.6). Both are the
authors' processed counts. [`tompros_counts`](@ref) rebuilds them from the swab results.
"""
function load_tompros2022(; data::Symbol = :fit)
    data in (:fit, :validation) || throw(ArgumentError("data must be :fit or :validation"))
    rows, col = _csv_rows(tompros2022_path(data))
    it, id, iS, iI = col("Tanks"), col("Days"), col("S"), col("I")
    return _tanks_from_rows([(tank = parse(Int, r[it]), day = parse(Float64, r[id]),
        S = parse(Int, r[iS]), I = parse(Int, r[iI])) for r in rows])
end

"""
    load_tompros_swabs() -> Vector{NamedTuple}

The swab results of Tompros et al. (2022), one row per newt and swab: `sample`, `tank`,
`infected_start` (`true` for the newts exposed to Bsal before the experiment), `prevalence`
(initial infection prevalence of the tank), `day` (days since the start, or `nothing` for a
necropsy), `positive`, `died` and `survival` (days survived).
"""
function load_tompros_swabs()
    rows, col = _csv_rows(tompros2022_path(:swabs))
    is, itr, itk, ip, id, ipos, idde, isv = col("Sample"), col("Treatment"), col("Tank"),
        col("InfectionPrev"), col("Day"), col("Pos"), col("DDE"), col("DaysSurvival")
    return [(sample = parse(Int, replace(r[is], "HDT " => "")), tank = parse(Int, r[itk]),
             infected_start = r[itr] != "Susceptible", prevalence = parse(Float64, r[ip]),
             day = r[id] == "Necropsy" ? nothing : parse(Float64, r[id]),
             positive = r[ipos] == "1", died = r[idde] == "1",
             survival = parse(Float64, r[isv])) for r in rows]
end

"""
    clean_infection_series(positive, necropsy) -> Vector{Bool}

The decision rules of Tompros et al. (2022, section 2.3), as in the authors'
`CodeFile1_PreProcess.R`, for one newt's swab series in time order (`necropsy` marks a swab
taken at death):

1. a single negative between two positives is a false negative, and so is a negative
   necropsy straight after a positive swab;
2. then a single positive between two negatives is a false positive, and so is a positive
   necropsy straight after a negative swab.

Each rule is applied to the whole series at once (vectorised), as in the R code.
"""
function clean_infection_series(positive::AbstractVector{Bool}, necropsy::AbstractVector{Bool})
    n = length(positive)
    pos = collect(positive)
    n < 2 && return pos
    flip = falses(n)
    for j in 2:(n - 1)
        flip[j] = pos[j - 1] && !pos[j] && pos[j + 1]
    end
    flip[n] = necropsy[n] && pos[n - 1]
    pos[flip] .= true
    flip .= false
    for j in 2:(n - 1)
        flip[j] = !pos[j - 1] && pos[j] && !pos[j + 1]
    end
    flip[n] = necropsy[n] && !pos[n - 1]
    pos[flip] .= false
    return pos
end

"""
    tompros_status(swabs = load_tompros_swabs()) -> Dict{Int,NamedTuple}

Cleaned infection histories of each newt, keyed by sample number: `tank`,
`infected_start`, `days` (with day 0 added, positive for the newts exposed beforehand) and
`positive` after [`clean_infection_series`](@ref), plus `died` and `survival`. Necropsies
are placed at the day of death, as in the authors' code.
"""
function tompros_status(swabs = load_tompros_swabs())
    out = Dict{Int,NamedTuple}()
    for s in unique(r.sample for r in swabs)
        rs = filter(r -> r.sample == s, swabs)
        r1 = first(rs)
        days = [r.day === nothing ? r.survival : r.day for r in rs]
        nec = [r.day === nothing for r in rs]
        pos = [r.positive for r in rs]
        push!(days, 0.0); push!(nec, false); push!(pos, r1.infected_start)
        o = sortperm(days; alg = Base.Sort.DEFAULT_STABLE)
        days, nec, pos = days[o], nec[o], pos[o]
        out[s] = (tank = r1.tank, infected_start = r1.infected_start, days = days,
            necropsy = nec, positive = clean_infection_series(pos, nec), died = r1.died,
            survival = r1.survival)
    end
    return out
end

"""
    tompros_counts(status = tompros_status(); days = (0, 2, 5), tanks = 1:20)
        -> Vector{TomprosTank}

Numbers of susceptible (negative) and infected (positive) newts per tank on `days`, rebuilt
from the cleaned swab histories as in `CodeFile1_PreProcess.R`. With the defaults it
reproduces the authors' processed data, [`load_tompros2022`](@ref)`()`.
"""
function tompros_counts(status = tompros_status(); days = (0, 2, 5), tanks = 1:20)
    rows = NamedTuple[]
    for k in tanks, d in days
        S = I = 0
        for (_, s) in status
            s.tank == k || continue
            for (j, t) in enumerate(s.days)
                t == d || continue
                s.positive[j] ? (I += 1) : (S += 1)
            end
        end
        push!(rows, (tank = k, day = float(d), S = S, I = I))
    end
    return _tanks_from_rows(rows)
end

const TOMPROS_FORMS = (
    frequency = (FrequencyDependent, (:β,), [0.2]),
    density = (DensityDependent, (:β,), [0.02]),
    quadratic = (QuadraticSaturation, (:κ₀, :κ₁, :κ₂), [1.0, 1.0, 1.0]),
    structural = (StructuralSaturation, (:κ₀, :κ_S, :κ_I, :κ₂), [1.0, 1.0, 1.0, 1.0]),
    power = (PowerLaw, (:β, :a, :b), [0.2, 1.0, 1.0]),
    prevalence = (PrevalenceDependent, (:β,), [3.0]),
    prevalence_power = (PrevalencePowerLaw, (:β, :a, :b), [3.0, 1.0, 1.0]),
    prevalence_quadratic = (PrevalenceSaturation, (:κ₀, :κ₁, :κ₂), [0.3, 1.0, 1.0]))

"""
    TomprosModel(form::Symbol; data = load_tompros2022(), lower = 1e-6)

Likelihood of the tank epidemics of Tompros et al. (2022) under one of the eight candidate
transmission functions of their Table 2. `form` is one of:

| `form` | transmission function | parameters |
|---|---|---|
| `:frequency` | [`FrequencyDependent`](@ref) | `β` |
| `:density` | [`DensityDependent`](@ref) | `β` |
| `:quadratic` | [`QuadraticSaturation`](@ref) | `κ₀, κ₁, κ₂` |
| `:structural` | [`StructuralSaturation`](@ref) | `κ₀, κ_S, κ_I, κ₂` |
| `:power` | [`PowerLaw`](@ref) | `β, a, b` |
| `:prevalence` | [`PrevalenceDependent`](@ref) | `β` |
| `:prevalence_power` | [`PrevalencePowerLaw`](@ref) | `β, a, b` |
| `:prevalence_quadratic` | [`PrevalenceSaturation`](@ref) | `κ₀, κ₁, κ₂` |

Time is in days. Each tank contributes the log probability of the number susceptible at
each swab day given the numbers on the previous swab day
([`susceptible_distribution`](@ref)).

The model implements the LogDensityProblems interface. The log density is the
log-likelihood with a flat prior on the **log** parameters; parameters below `lower` (the
paper's bound for `optim`) are set to `lower`. Fit it with [`maximum_likelihood`](@ref).
"""
struct TomprosModel
    form::Symbol
    data::Vector{TomprosTank}
    lower::Float64
end
function TomprosModel(form::Symbol; data = load_tompros2022(), lower::Real = 1e-6)
    haskey(TOMPROS_FORMS, form) ||
        throw(ArgumentError("form must be one of $(keys(TOMPROS_FORMS)); got :$form"))
    return TomprosModel(form, collect(data), float(lower))
end

parameter_names(model::TomprosModel) = TOMPROS_FORMS[model.form][2]
LogDensityProblems.dimension(model::TomprosModel) = length(parameter_names(model))
LogDensityProblems.capabilities(::Type{TomprosModel}) = LogDensityProblems.LogDensityOrder{0}()

"""
    transmission_function(model::TomprosModel, θ) -> TransmissionFunction

The transmission function of `model` with natural-scale parameters `θ`, in the order of
[`parameter_names`](@ref).
"""
function transmission_function(model::TomprosModel, θ::AbstractVector)
    length(θ) == LogDensityProblems.dimension(model) ||
        throw(DimensionMismatch("expected $(LogDensityProblems.dimension(model)) parameters"))
    return TOMPROS_FORMS[model.form][1](θ...)
end

# log-likelihood of the tanks at natural-scale parameters θ
function _tompros_loglik(model::TomprosModel, θ)
    tf = transmission_function(model, θ)
    ll = 0.0
    for tank in model.data, j in 1:(length(tank.days) - 1)
        S0, I0, S1 = tank.S[j], tank.I[j], tank.S[j + 1]
        S1 <= S0 || return -Inf
        S0 == 0 && continue
        p = susceptible_distribution(tf, S0, I0, tank.days[j + 1] - tank.days[j])
        ll += log(p[S1 + 1])
    end
    return ll
end

function LogDensityProblems.logdensity(model::TomprosModel, x::AbstractVector)
    length(x) == LogDensityProblems.dimension(model) ||
        throw(DimensionMismatch("expected $(LogDensityProblems.dimension(model)) parameters"))
    all(isfinite, x) || return -Inf
    θ = max.(exp.(x), model.lower)
    ll = _tompros_loglik(model, θ)
    return isfinite(ll) ? ll : -Inf
end

"""
    maximum_likelihood(model::TomprosModel; starts = 6, iters = 1500, rng = Random.Xoshiro(1))

Maximum likelihood fit of a [`TomprosModel`](@ref) by Nelder–Mead on the log parameters.
The first start is a default point; the others are drawn as in the authors'
`CodeFile3_Optimisation.R` (rates and κ's uniform on 0.1–2 times the default, exponents
uniform on 0.5–1.5). Each start is polished by a second, smaller simplex. Parameters are
held at or above `model.lower`.

Returns a `NamedTuple` with `model`, `x` (log parameters), `estimates` (natural scale, by
name; values at the lower bound are reported as `model.lower`), `nll` (negative
log-likelihood), `k`, `aic = 2 nll + 2k` and `nlls` (best value from each start).
"""
function maximum_likelihood(model::TomprosModel; starts::Integer = 6, iters::Integer = 1500,
        rng::AbstractRNG = Random.Xoshiro(1))
    lo = log(model.lower)
    f(x) = -LogDensityProblems.logdensity(model, x)
    names = parameter_names(model)
    θ0 = TOMPROS_FORMS[model.form][3]
    best, bestf = log.(θ0), Inf
    nlls = Float64[]
    for s in 1:starts
        θ = s == 1 ? θ0 : [n in (:a, :b) ? 0.5 + rand(rng) : θ0[i] * (0.1 + 1.9rand(rng))
                           for (i, n) in enumerate(names)]
        x, fx = _nelder_mead_tol(f, log.(θ); iters = iters, step = 0.5, lower = lo)
        x, fx = _nelder_mead_tol(f, x; iters = iters ÷ 2, step = 0.05, lower = lo)
        push!(nlls, fx)
        fx < bestf && ((best, bestf) = (x, fx))
    end
    best = max.(best, lo)
    k = length(best)
    θ = exp.(best)
    θ[best .<= lo + 1e-6] .= model.lower      # report values at the bound as the bound
    est = NamedTuple{names}(Tuple(θ))
    return (model = model, x = best, estimates = est, nll = bestf, k = k,
        aic = 2bestf + 2k, nlls = nlls)
end

"""
    tompros_model_selection(; data = load_tompros2022(), forms = keys(TOMPROS_FORMS), kwargs...)
        -> Vector{NamedTuple}

Fit the candidate transmission functions of Tompros et al. (2022, Table 2) and rank them by
AIC. Each entry has the fields of [`maximum_likelihood`](@ref), plus `form` and `Δaic`.
Keywords are passed to `maximum_likelihood`.
"""
function tompros_model_selection(; data = load_tompros2022(), forms = keys(TOMPROS_FORMS),
        kwargs...)
    fits = [maximum_likelihood(TomprosModel(f; data = data); kwargs...) for f in forms]
    sort!(fits; by = f -> f.aic)
    a0 = first(fits).aic
    return [merge(f, (form = f.model.form, Δaic = f.aic - a0)) for f in fits]
end

"""
    likelihood_interval(model::TomprosModel, values; level = 0.95) -> (; lower, upper, mode)

Interval of a one-parameter [`TomprosModel`](@ref) as Tompros et al. (2022) compute it
(`CodeFile4_FreqDependentProfile.R`, Fig. 3): the likelihood is evaluated at the parameter
`values` (natural scale, increasing and equally spaced), normalised to unit area by the
trapezium rule, and the interval is cut where the cumulative area reaches
`(1 − level)/2` and `(1 + level)/2`, by linear interpolation. This is a credible interval
under a flat prior on `values`, not a likelihood-ratio interval. `mode` is the grid value
with the largest likelihood.
"""
function likelihood_interval(model::TomprosModel, values::AbstractVector; level::Real = 0.95)
    LogDensityProblems.dimension(model) == 1 ||
        throw(ArgumentError("likelihood_interval needs a one-parameter model"))
    v = collect(float.(values))
    ll = [_tompros_loglik(model, [max(b, model.lower)]) for b in v]
    L = exp.(ll .- maximum(ll))
    h = (v[end] - v[1]) / (length(v) - 1)
    cum = zeros(length(v))
    for i in 2:length(v)
        cum[i] = cum[i - 1] + h * (L[i - 1] + L[i]) / 2
    end
    cum ./= cum[end]
    q(p) = (j = findfirst(>=(p), cum); j == 1 ? v[1] :
        v[j - 1] + (p - cum[j - 1]) / (cum[j] - cum[j - 1]) * (v[j] - v[j - 1]))
    return (lower = q((1 - level) / 2), upper = q((1 + level) / 2), mode = v[argmax(L)])
end

"""
    tompros_negloglik(form::Symbol, θ; data = load_tompros2022()) -> Float64

Negative log-likelihood of the tank data under transmission function `form` with
natural-scale parameters `θ` (in the order of [`parameter_names`](@ref)), for example the
estimates printed in Table 2 of Tompros et al. (2022).
"""
tompros_negloglik(form::Symbol, θ::AbstractVector; data = load_tompros2022()) =
    -_tompros_loglik(TomprosModel(form; data = data), float.(θ))

# Nelder–Mead (as nelder_mead in knapp2022.jl) with every vertex projected onto x ≥ lower,
# stopping once the simplex has collapsed: the spread of the function values is below ftol
# and every vertex is within xtol of the best.
function _nelder_mead_tol(f0, x0::AbstractVector; iters::Integer = 1500, step::Real = 0.5,
        ftol::Real = 1e-10, xtol::Real = 1e-8, lower::Real = -Inf)
    n = length(x0)
    f(x) = f0(max.(x, lower))
    xs = [max.(float.(copy(x0)), lower) for _ in 1:(n + 1)]
    for i in 1:n
        xs[i + 1][i] += step
    end
    fs = map(f, xs)
    for _ in 1:iters
        o = sortperm(fs)
        xs, fs = xs[o], fs[o]
        if isfinite(fs[1]) && fs[end] - fs[1] < ftol &&
           maximum(maximum(abs, x - xs[1]) for x in xs) < xtol
            break
        end
        c = sum(xs[1:n]) / n
        xr = max.(c + (c - xs[end]), lower); fr = f(xr)
        if fr < fs[1]
            xe = max.(c + 2(c - xs[end]), lower); fe = f(xe)
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

# One fit per form per session, for the validation registry (data/validation/targets.toml).
const _TOMPROS_REFERENCE = Dict{Symbol,Any}()
_tompros_reference_fit(form::Symbol) =
    get!(() -> maximum_likelihood(TomprosModel(form)), _TOMPROS_REFERENCE, form)
