# Fits of Gajewski et al. 2021 (Ecol Evol 11:17920, doi:10.1002/ece3.8379) to the authors'
# data (Zenodo doi:10.5281/zenodo.4525954; not bundled, no licence file). This reproduces a
# published analysis; it is not a new model.
#
# Posterior of one TPC fit (Methods 2.5 and the authors' HierarchicalModel_<form>.bug), all
# constant-temperature wells pooled, with one Y0, K and d shared by every temperature:
#   μ_i = Eq. 1 with r = f(T_i, θ)
#   AS CODED  Briere 2, Ratkowsky: Y_i ~ LogNormal(log(μ_i + 1e-7), σ)
#             Ikemoto, Logan 10, Stinner: Y_i ~ Normal(μ_i, σ) truncated to (0, ∞)
#             Ikemoto only: the delay in μ is exp(d)
#   AS PRINTED (Methods 2.5): Y_i ~ LogNormal(μ_i, σ) for every form, delay d.
# The paper writes "log N(μ_i, τ_i)"; the code's lognormal is on log(μ + 1e-7), which we
# read as the intended meaning of the printed form. Priors are those of the code (and of
# Appendix S1), read from the compendium sets "prior" and "prior_<form>"; JAGS normals are
# (mean, precision) and gammas (shape, rate). Temperatures enter the kelvin curves as
# T + 273.15, as in RunModels.R.

"""
    gajewski2021_path() -> String

Folder with the data of Gajewski et al. (2021) (`optical_QLD.csv`, `Vary_OD.csv`): the
authors' Zenodo release [10.5281/zenodo.4525954](https://doi.org/10.5281/zenodo.4525954),
extracted. It is not bundled (the release has no licence file). The path is
`ENV["CHYTRIDIOMYCOSIS_GAJEWSKI2021"]` if set, else
`papers/gajewski2021/zenodo/gzach93-Predict_Bd_Growth_VaryEnvir-c4ade0e` under the first
ancestor of the package directory that has it. Throws an `ArgumentError` if neither exists.
"""
function gajewski2021_path()
    p = _gajewski2021_find()
    p === nothing && throw(ArgumentError(
        "Gajewski et al. 2021 data not found: download doi:10.5281/zenodo.4525954, extract it " *
        "and set ENV[\"CHYTRIDIOMYCOSIS_GAJEWSKI2021\"] to the folder holding optical_QLD.csv"))
    return p
end

const GAJEWSKI2021_DEPOSIT = joinpath("gajewski2021", "zenodo",
    "gzach93-Predict_Bd_Growth_VaryEnvir-c4ade0e")

function _gajewski2021_find()
    e = get(ENV, "CHYTRIDIOMYCOSIS_GAJEWSKI2021", "")
    isempty(e) || return isfile(joinpath(e, "optical_QLD.csv")) ? e : nothing
    d = pkgdir(@__MODULE__)
    while true
        c = joinpath(d, "papers", GAJEWSKI2021_DEPOSIT)
        isfile(joinpath(c, "optical_QLD.csv")) && return c
        parent = dirname(d)
        parent == d && return nothing
        d = parent
    end
end

_gj_unquote(s) = strip(strip(s), '"')

function _gj_csv(path)
    lines = filter(!isempty ∘ strip, readlines(path))
    hdr = _gj_unquote.(split(lines[1], ','))
    rows = [_gj_unquote.(split(l, ',')) for l in lines[2:end]]
    col(n) = (i = findfirst(==(n), hdr); i === nothing ? throw(ArgumentError("no $n in $path")) : i)
    return rows, col
end

const GAJEWSKI_REGIME_CODES = Dict("Dry_Low_Air" => :DLA, "Wet_High_Frog" => :WHF,
    "Wet_Low_Frog" => :WLF)

"""
    load_gajewski2021(; path = gajewski2021_path()) -> (constant, varying)

The optical density data of Gajewski et al. (2021):

- `constant`: Stevenson et al.'s (2013) wells at ten constant temperatures (4693 readings),
  as vectors `T` (°C), `t` (days) and `y` (the `adjusted` OD the authors fitted);
- `varying`: a `Dict` from regime (`:DLA`, `:WHF`, `:WLF`) to vectors `t` and `y` (Stevenson
  et al. 2020).
"""
function load_gajewski2021(; path::AbstractString = gajewski2021_path())
    rows, col = _gj_csv(joinpath(path, "optical_QLD.csv"))
    it, id, iy = col("Treatment"), col("Day"), col("adjusted")
    T = [parse(Float64, chop(r[it])) for r in rows]
    constant = (T = T, t = [parse(Float64, r[id]) for r in rows],
                y = [parse(Float64, r[iy]) for r in rows])
    rows, col = _gj_csv(joinpath(path, "Vary_OD.csv"))
    ir, id, iy = col("trt"), col("day"), col("od")
    varying = Dict(g => (t = [parse(Float64, r[id]) for r in rows if GAJEWSKI_REGIME_CODES[r[ir]] === g],
                         y = [parse(Float64, r[iy]) for r in rows if GAJEWSKI_REGIME_CODES[r[ir]] === g])
                   for g in GAJEWSKI_REGIMES)
    return (constant = constant, varying = varying)
end

const GAJEWSKI_TPC_NAMES = Dict(:briere2 => (:c, :Tmin, :Tmax, :b),
    :ratkowsky => (:c, :k, :Tmin, :Tmax),
    :ikemoto => (:ρ, :TΦ, :ΔHA, :ΔHL, :TL, :ΔHH, :TH),
    :logan10 => (:α, :c, :b, :Tmax, :ΔT),
    :stinner => (:C, :k1, :k2, :Topt))

"""
    GajewskiModel(; form = :stinner, likelihood = :as_coded, data = load_gajewski2021().constant)

The posterior of one hierarchical fit of Gajewski et al. (2021): Eq. 1 for every
constant-temperature well, with `r` given by TPC `form` (see [`gajewski_curve`](@ref)) and
one `Y0`, `K` and `d` shared by all temperatures.

- `likelihood = :as_coded` follows the authors' JAGS models: lognormal for `:briere2` and
  `:ratkowsky`, normal truncated at zero for `:ikemoto`, `:logan10` and `:stinner`, and the
  delay `exp(d)` for `:ikemoto`. The paper (Methods 2.5) prints a lognormal likelihood for
  all five; `likelihood = :lognormal` gives that, with the delay `d`.
- Priors are the authors' (compendium sets `"prior"` and `"prior_<form>"`).

Parameters are on their natural scale, in the order of [`parameter_names`](@ref):
`Y0, K, d, σ`, then the curve's own. The model implements the LogDensityProblems interface;
the log density is the log posterior up to a constant (`-Inf` outside the prior support).
Use [`posterior_mode`](@ref) for its mode and `scripts/gajewski2021_mcmc.jl` to sample it.
"""
struct GajewskiModel
    form::Symbol
    likelihood::Symbol
    T::Vector{Float64}
    t::Vector{Float64}
    y::Vector{Float64}
    temps::Vector{Float64}
    index::Vector{Int}
    prior::NamedTuple
    R::Float64
end

function GajewskiModel(; form::Symbol = :stinner, likelihood::Symbol = :as_coded,
        data = load_gajewski2021().constant)
    _gajewski_check_form(form)
    likelihood in (:as_coded, :lognormal) ||
        throw(ArgumentError("likelihood must be :as_coded or :lognormal"))
    temps = sort(unique(data.T))
    index = [searchsortedfirst(temps, T) for T in data.T]
    pr = merge(parameter_set(:gajewski2021; set = "prior"),
               parameter_set(:gajewski2021; set = "prior_" * string(form)))
    R = parameter_set(:gajewski2021; set = "prior_ikemoto").R
    return GajewskiModel(form, likelihood, collect(Float64, data.T), collect(Float64, data.t),
        collect(Float64, data.y), temps, index, pr, R)
end

parameter_names(m::GajewskiModel) = (:Y0, :K, :d, :σ, GAJEWSKI_TPC_NAMES[m.form]...)
LogDensityProblems.dimension(m::GajewskiModel) = length(parameter_names(m))
LogDensityProblems.capabilities(::Type{GajewskiModel}) = LogDensityProblems.LogDensityOrder{0}()

_gj_lognormal(x, lo, hi) = lo <= x <= hi ? -log(hi - lo) : -Inf
_gj_exp(x, rate) = x >= 0 ? log(rate) - rate * x : -Inf
_gj_gamma(x, shape, rate) = x > 0 ? logpdf(Gamma(shape, 1 / rate), x) : -Inf
_gj_norm(x, mean, prec) = logpdf(Normal(mean, 1 / sqrt(prec)), x)
_gj_tnorm(x, mean, prec, lo, hi) = lo <= x <= hi ?
    logpdf(Normal(mean, 1 / sqrt(prec)), x) -
    log(cdf(Normal(mean, 1 / sqrt(prec)), hi) - cdf(Normal(mean, 1 / sqrt(prec)), lo)) : -Inf

function _gajewski_logprior(form::Symbol, θ, p)
    lp = _gj_exp(θ.σ, p.σ_rate) + _gj_lognormal(θ.Y0, p.Y0_lo, p.Y0_hi) +
         _gj_exp(θ.d, p.d_rate) + _gj_lognormal(θ.K, p.K_lo, p.K_hi)
    if form === :briere2
        lp += _gj_gamma(θ.c, p.c_shape, p.c_rate) + _gj_norm(θ.Tmin, p.Tmin_mean, p.Tmin_precision) +
              _gj_lognormal(θ.Tmax, p.Tmax_lo, p.Tmax_hi) + _gj_gamma(θ.b, p.b_shape, p.b_rate)
    elseif form === :ratkowsky
        lp += _gj_gamma(θ.c, p.c_shape, p.c_rate) + _gj_gamma(θ.k, p.k_shape, p.k_rate) +
              _gj_norm(θ.Tmin, p.Tmin_mean, p.Tmin_precision) +
              _gj_lognormal(θ.Tmax, p.Tmax_lo, p.Tmax_hi)
    elseif form === :ikemoto
        lp += _gj_lognormal(θ.ρ, p.ρ_lo, p.ρ_hi) + _gj_lognormal(θ.TΦ, p.TΦ_lo, p.TΦ_hi) +
              _gj_gamma(θ.TH, p.TH_shape, p.TH_rate) + _gj_gamma(θ.TL, p.TL_shape, p.TL_rate) +
              _gj_norm(θ.ΔHA, p.ΔHA_mean, p.ΔHA_precision) +
              _gj_norm(θ.ΔHH, p.ΔHH_mean, p.ΔHH_precision) +
              _gj_norm(θ.ΔHL, p.ΔHL_mean, p.ΔHL_precision)
    elseif form === :logan10
        lp += _gj_lognormal(θ.α, p.α_lo, p.α_hi) + _gj_norm(θ.c, p.c_mean, p.c_precision) +
              _gj_gamma(θ.b, p.b_shape, p.b_rate) + _gj_gamma(θ.ΔT, p.ΔT_shape, p.ΔT_rate) +
              _gj_lognormal(θ.Tmax, p.Tmax_lo, p.Tmax_hi)
    else
        lp += _gj_gamma(θ.C, p.C_shape, p.C_rate) +
              _gj_tnorm(θ.k1, p.k1_mean, p.k1_precision, p.k1_lo, p.k1_hi) +
              _gj_tnorm(θ.k2, p.k2_mean, p.k2_precision, p.k2_lo, p.k2_hi) +
              _gj_lognormal(θ.Topt, p.Topt_lo, p.Topt_hi)
    end
    return lp
end

_gj_named(m, x) = NamedTuple{parameter_names(m)}(Tuple(x))

# log-likelihood of the data at named parameters θ (prior support already checked)
function _gajewski_loglik(m::GajewskiModel, θ)
    curve = gajewski_curve(m.form, merge(θ, (R = m.R,)))
    rs = [curve(T) for T in m.temps]
    coded = m.likelihood === :as_coded
    lognormal = !coded || m.form in (:briere2, :ratkowsky)
    delay = coded && m.form === :ikemoto ? exp(θ.d) : θ.d
    ll = 0.0
    for i in eachindex(m.y)
        μ = delayed_logistic(m.t[i], θ.Y0, θ.K, rs[m.index[i]], delay)
        if lognormal
            ll += logpdf(LogNormal(log(μ + 1e-7), θ.σ), m.y[i])
        else
            D = Normal(μ, θ.σ)
            ll += logpdf(D, m.y[i]) - logccdf(D, 0.0)
        end
    end
    return ll
end

function _gajewski_valid(m::GajewskiModel, θ)
    θ.σ > 0 || return false
    m.form === :briere2 && return θ.Tmin < θ.Tmax && θ.b > 0
    m.form === :ratkowsky && return θ.Tmin < θ.Tmax
    return true
end

function LogDensityProblems.logdensity(m::GajewskiModel, x::AbstractVector)
    length(x) == LogDensityProblems.dimension(m) ||
        throw(DimensionMismatch("expected $(LogDensityProblems.dimension(m)) parameters"))
    all(isfinite, x) || return -Inf
    θ = _gj_named(m, x)
    _gajewski_valid(m, θ) || return -Inf
    lp = _gajewski_logprior(m.form, θ, m.prior)
    isfinite(lp) || return -Inf
    ll = _gajewski_loglik(m, θ)
    return isfinite(ll) ? lp + ll : -Inf
end

"""
    gajewski_start(m::GajewskiModel) -> Vector{Float64}

A starting point for [`posterior_mode`](@ref) or a sampler: the printed posterior medians of
the fit (compendium set `"<form>"`), with two exceptions. The Logan-10 `α` starts at the
midpoint of its printed HDI, because the printed median is inconsistent with it; and `σ`,
which the paper does not report, is the root mean square residual at that point on the scale
of the likelihood. For `:ikemoto` with the as-coded likelihood, `d` is the printed value,
which the code exponentiates.
"""
function gajewski_start(m::GajewskiModel)
    p = parameter_set(:gajewski2021; set = string(m.form))
    θ = Dict{Symbol,Float64}(n => (n === :σ ? 1.0 : Float64(p[n])) for n in parameter_names(m))
    if m.form === :logan10
        lo, hi = _gajewski_interval("logan10", :α)
        θ[:α] = (lo + hi) / 2
    end
    x = [θ[n] for n in parameter_names(m)]
    nt = _gj_named(m, x)
    curve = gajewski_curve(m.form, merge(nt, (R = m.R,)))
    coded = m.likelihood === :as_coded
    lognormal = !coded || m.form in (:briere2, :ratkowsky)
    delay = coded && m.form === :ikemoto ? exp(nt.d) : nt.d
    res = map(eachindex(m.y)) do i
        μ = delayed_logistic(m.t[i], nt.Y0, nt.K, curve(m.T[i]), delay)
        lognormal ? log(m.y[i]) - log(μ + 1e-7) : m.y[i] - μ
    end
    x[4] = sqrt(sum(abs2, res) / length(res))
    return x
end

# The printed interval (uncertainty.lo, uncertainty.hi) of a compendium record; the reader
# (src/compendium.jl) keeps only the value, so this reads the TOML file.
function _gajewski_interval(set::AbstractString, symbol::Symbol)
    raw = TOML.parsefile(joinpath(COMPENDIUM_DIR, "gajewski2021.toml"))
    r = only(r for r in raw["parameter"] if get(r, "set", "") == set && r["symbol"] == string(symbol))
    return Float64(r["uncertainty"]["lo"]), Float64(r["uncertainty"]["hi"])
end

"""
    GajewskiLogisticModel(regime; data = load_gajewski2021().varying)

The logistic model (Eq. 1) without a TPC, fitted by Gajewski et al. (2021) to each
varying-temperature regime (Methods 2.7; Appendix S1 Table B.1; the authors'
`jags-logistic.bug`): `Y_i ~ Normal(μ_i, σ)`, not truncated, with priors `r ~ Gamma(1, 1)`,
`K ~ Uniform(0.08, 0.2)` (compendium set `"prior_logistic"`) and the `Y0`, `d` and `σ` priors
of the hierarchical fits. Parameters `Y0, K, d, r, σ` on their natural scale. Implements the
LogDensityProblems interface (the log posterior up to a constant).
"""
struct GajewskiLogisticModel
    regime::Symbol
    t::Vector{Float64}
    y::Vector{Float64}
    prior::NamedTuple
end

function GajewskiLogisticModel(regime::Symbol; data = load_gajewski2021().varying)
    regime in GAJEWSKI_REGIMES ||
        throw(ArgumentError("regime must be one of $(GAJEWSKI_REGIMES); got :$regime"))
    pr = merge(parameter_set(:gajewski2021; set = "prior"),
               parameter_set(:gajewski2021; set = "prior_logistic"))
    return GajewskiLogisticModel(regime, collect(Float64, data[regime].t),
        collect(Float64, data[regime].y), pr)
end

parameter_names(::GajewskiLogisticModel) = (:Y0, :K, :d, :r, :σ)
LogDensityProblems.dimension(::GajewskiLogisticModel) = 5
LogDensityProblems.capabilities(::Type{GajewskiLogisticModel}) =
    LogDensityProblems.LogDensityOrder{0}()

function LogDensityProblems.logdensity(m::GajewskiLogisticModel, x::AbstractVector)
    length(x) == 5 || throw(DimensionMismatch("expected 5 parameters"))
    all(isfinite, x) || return -Inf
    Y0, K, d, r, σ = x
    p = m.prior
    σ > 0 || return -Inf
    lp = _gj_exp(σ, p.σ_rate) + _gj_lognormal(Y0, p.Y0_lo, p.Y0_hi) + _gj_exp(d, p.d_rate) +
         _gj_lognormal(K, p.K_lo, p.K_hi) + _gj_gamma(r, p.r_shape, p.r_rate)
    isfinite(lp) || return -Inf
    ll = sum(logpdf(Normal(delayed_logistic(m.t[i], Y0, K, r, d), σ), m.y[i])
             for i in eachindex(m.y))
    return isfinite(ll) ? lp + ll : -Inf
end

"""
    gajewski_start(m::GajewskiLogisticModel)

The printed posterior medians of Table B.1 for the regime, with `σ` the root mean square
residual there.
"""
function gajewski_start(m::GajewskiLogisticModel)
    p = parameter_set(:gajewski2021; set = "logistic_" * string(m.regime))
    σ = sqrt(sum(abs2, m.y[i] - delayed_logistic(m.t[i], p.Y0, p.K, p.r, p.d)
                 for i in eachindex(m.y)) / length(m.y))
    return [p.Y0, p.K, p.d, p.r, σ]
end

"""
    posterior_mode(m; x0 = gajewski_start(m), rounds = 3, iters = 4000)

The mode of the log posterior of a [`GajewskiModel`](@ref) or
[`GajewskiLogisticModel`](@ref), by Nelder–Mead from `x0` in coordinates scaled by
`max(|x0|, 1e-3)`, restarted `rounds` times from the best point with a smaller simplex.
With thousands of readings the mode is close to the posterior medians that the paper
reports. Returns a `NamedTuple` with `model`, `x`, `estimates` (by name) and `logdensity`.
"""
function posterior_mode(m::Union{GajewskiModel,GajewskiLogisticModel};
        x0::AbstractVector = gajewski_start(m), rounds::Integer = 3, iters::Integer = 4000)
    s = max.(abs.(x0), 1e-3)
    f(z) = -LogDensityProblems.logdensity(m, x0 .+ s .* z)
    z, fz = zeros(length(x0)), f(zeros(length(x0)))
    isfinite(fz) || throw(ArgumentError("the log density is not finite at x0"))
    for k in 1:rounds
        z, fz = nelder_mead(f, z; iters, step = 0.1 / 4^(k - 1))
    end
    x = x0 .+ s .* z
    return (model = m, x = x, estimates = _gj_named(m, x), logdensity = -fz)
end


# One posterior mode per fit per session, for the validation registry (pending targets that
# need the authors' unbundled data)
const _GAJEWSKI_REFERENCE = Dict{Any,Any}()
_gajewski_reference_mode(form::Symbol) =
    get!(() -> posterior_mode(GajewskiModel(; form)), _GAJEWSKI_REFERENCE, form)
_gajewski_reference_mode(m::Symbol, regime::Symbol) =
    get!(() -> posterior_mode(GajewskiLogisticModel(regime)), _GAJEWSKI_REFERENCE, (m, regime))

# DIC ≈ D(mode) + 2k, for comparison with the penalized deviances of Table 1
function _gajewski_dic(form::Symbol)
    fit = _gajewski_reference_mode(form)
    return -2 * _gajewski_loglik(fit.model, fit.estimates) + 2 * length(fit.x)
end
