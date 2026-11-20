# Gajewski, Stevenson, Pike, Roznik, Alford & Johnson (2021) Ecol Evol 11:17920-17931,
# doi:10.1002/ece3.8379. Bd optical density (OD) in culture as a delayed logistic curve whose
# growth rate follows a thermal performance curve (TPC). Equations and values checked against
# the full text (PMC8717292), Appendix S1, and the authors' code (Zenodo
# doi:10.5281/zenodo.4525954: functions.R, TPCurves.R, RunModels.R, HierarchicalModel_*.bug).
#
#   Y(t) = Y0                                         t < d       (Eq. 1)
#        = K Y0 / (Y0 + (K − Y0) exp(−r(T) (t − d)))  t ≥ d
#   dx/dt = r(T(t)) x (1 − x/K) 1(t > d),  x(0) = Y0              (Eq. 3)
#
# Time in days, OD in absorbance units at 492 nm (minus heat-killed controls), T in °C.
# Eq. 3 with constant T is Eq. 1. The five TPCs are Brière-2 ([`Briere`](@ref)), Ratkowsky
# ([`Ratkowsky`](@ref)), Ikemoto ([`SharpeSchoolfieldIkemoto`](@ref)), Logan-10
# ([`Logan10`](@ref)) and Stinner ([`Stinner`](@ref)), in the forms of TPCurves.R.
#
# Where the paper and the code differ (all documented, none tuned):
# - Kelvin. The Ratkowsky and Ikemoto fits use T + 273.15 (RunModels.R); the predictions use
#   T + 273.16 (functions.R). `kelvin_offset` selects either.
# - Wet Low Frog regime. The last segment of temperature.wlf() returns 25 without the kelvin
#   offset, so the Ratkowsky and Ikemoto predictions see 25 K (r = 0) for 3.3 h a day.
#   `gajewski_prediction(...; as_coded = true)` reproduces this.
# - Ikemoto delay. HierarchicalModel_Ikemoto.bug uses exp(d) as the delay, but Table 2 prints d
#   and Predictions.R uses d. `ikemoto_delay = :fit` gives exp(d).
# - Logan-10 α. Appendix S1 prints a median of 0.134 with an HDI of (0.885, 0.957); see the
#   compendium notes. The keyed constructor uses the printed value unless `α` is given.
#
# Choice where the paper is silent: Table B.2 summarises r(T(t)) over one day at the posterior
# medians; we take the median and the 95% HDI of r on a grid of `step` days over [0, 1).

const GAJEWSKI_FORMS = (:briere2, :ratkowsky, :ikemoto, :logan10, :stinner)
const GAJEWSKI_REGIMES = (:DLA, :WHF, :WLF)
const GAJEWSKI_KELVIN_FORMS = (:ratkowsky, :ikemoto)
const GAJEWSKI_ISOLATE = "Paluma-Lgenimaculata #2-2011-CO (QLD)"
const GAJEWSKI_PREDICTION_KELVIN = 273.16     # functions.R, temperature.*(…, Kelvin = TRUE)

"""
    delayed_logistic(t, Y0, K, r, d = 0)

The logistic curve with a delay (Gajewski et al. 2021, Eq. 1; Voyles et al. 2017): `Y0` for
`t < d`, then `K Y0 / (Y0 + (K − Y0) exp(−r (t − d)))`. With `d = 0` it is the standard
logistic curve that Sheets et al. (2021) fit to normalised MTT readings.
"""
delayed_logistic(t, Y0, K, r, d = 0) =
    t < d ? Y0 * one(r * t) : K * Y0 / (Y0 + (K - Y0) * exp(-r * (t - d)))

function _gajewski_check_form(form::Symbol)
    form in GAJEWSKI_FORMS ||
        throw(ArgumentError("form must be one of $(GAJEWSKI_FORMS); got :$form"))
    return form
end

"""
    gajewski_curve(form, θ; kelvin_offset = 273.15) -> ThermalPerformanceCurve

The thermal performance curve `form` of Gajewski et al. (2021) with parameters `θ` (any
property collection, such as a `NamedTuple` or the compendium set of that form), called with
temperature in °C:

| `form` | Curve | Parameters |
|---|---|---|
| `:briere2` | [`Briere`](@ref), `m = b` | `c, Tmin, Tmax, b` (°C) |
| `:ratkowsky` | [`Ratkowsky`](@ref) | `c, k, Tmin, Tmax` (kelvin) |
| `:ikemoto` | [`SharpeSchoolfieldIkemoto`](@ref) | `ρ, TΦ, ΔHA, ΔHL, TL, ΔHH, TH` (K, cal mol⁻¹), `R` |
| `:logan10` | [`Logan10`](@ref) | `α, c, b, Tmax, ΔT` (°C) |
| `:stinner` | [`Stinner`](@ref) | `C, k1, k2, Topt` (°C) |

The kelvin temperatures of Ratkowsky and Ikemoto are converted with `kelvin_offset`:
273.15 in the authors' fits and 273.16 in their predictions.
"""
function gajewski_curve(form::Symbol, θ; kelvin_offset::Real = KELVIN_OFFSET)
    _gajewski_check_form(form)
    if form === :briere2
        return Briere(; a = θ.c, Tmin = θ.Tmin, Tmax = θ.Tmax, m = θ.b)
    elseif form === :ratkowsky
        return Ratkowsky(; b = θ.c, Tmin = θ.Tmin - kelvin_offset, c = θ.k,
                         Tmax = θ.Tmax - kelvin_offset)
    elseif form === :ikemoto
        return SharpeSchoolfieldIkemoto(; ρ = θ.ρ, TΦ = θ.TΦ, ΔHA = θ.ΔHA, ΔHL = θ.ΔHL,
            TL = θ.TL, ΔHH = θ.ΔHH, TH = θ.TH, R = θ.R, kelvin_offset)
    elseif form === :logan10
        return Logan10(; α = θ.α, c = θ.c, b = θ.b, Tmax = θ.Tmax, ΔT = θ.ΔT)
    else
        return Stinner(; C = θ.C, k1 = θ.k1, k2 = θ.k2, Topt = θ.Topt)
    end
end

"""
    GajewskiParameters(; form, Y0, K, d, tpc)
    GajewskiParameters(:gajewski2021; form = :stinner, kelvin_offset = 273.15,
                       ikemoto_delay = :printed, α = nothing)

The delayed logistic growth of Bd optical density with a temperature-dependent rate
(Gajewski et al. 2021, Eqs 1 and 3), for one of the five TPC fits `form` (see
[`gajewski_curve`](@ref)):

- `Y0`: initial optical density;
- `K`: maximum optical density;
- `d`: delay before growth (days);
- `tpc`: the growth rate `r(T)` (day⁻¹, `T` in °C), an [`IsolateTPC`](@ref).

The keyed constructor reads the posterior medians of that fit from the compendium (Appendix
S1 Tables B.3–B.7, Table 2). The medians of the parameters do not give the posterior median
curve, so summaries of the curve built from them (its optimum, peak, adjusted limits) only
approximate the paper's Table 1 and 2 values.

- `kelvin_offset`: 273.15 (the fits) or 273.16 (the authors' predictions), for Ratkowsky
  and Ikemoto.
- `ikemoto_delay`: `:printed` uses `d` = 0.007 as in Table 2 and Predictions.R; `:fit` uses
  `exp(d)`, the delay in the authors' JAGS model for the Ikemoto fit.
- `α`: replaces the Logan-10 `α`, whose printed median (0.134) contradicts its own HDI
  (0.885, 0.957).
"""
struct GajewskiParameters{C}
    form::Symbol
    Y0::Float64
    K::Float64
    d::Float64
    tpc::IsolateTPC{C}
end

function GajewskiParameters(; form::Symbol, Y0::Real, K::Real, d::Real, tpc)
    _gajewski_check_form(form)
    c = tpc isa IsolateTPC ? tpc : IsolateTPC("Bd", GAJEWSKI_ISOLATE, "",
        "logistic growth rate r of optical density (day⁻¹)", tpc, :gajewski2021)
    return GajewskiParameters(form, Float64(Y0), Float64(K), Float64(d), c)
end

function GajewskiParameters(key::Symbol; form::Symbol = :stinner,
        kelvin_offset::Real = KELVIN_OFFSET, ikemoto_delay::Symbol = :printed,
        α::Union{Nothing,Real} = nothing)
    _gajewski_check_form(form)
    ikemoto_delay in (:printed, :fit) ||
        throw(ArgumentError("ikemoto_delay must be :printed or :fit"))
    θ = parameter_set(key; set = string(form))
    if form === :ikemoto
        θ = merge(θ, (R = parameter_set(key; set = "prior_ikemoto").R,))
    end
    α === nothing || form !== :logan10 || (θ = merge(θ, (α = α,)))
    d = form === :ikemoto && ikemoto_delay === :fit ? exp(θ.d) : θ.d
    return GajewskiParameters(; form, Y0 = θ.Y0, K = θ.K, d,
        tpc = gajewski_curve(form, θ; kelvin_offset))
end

time_unit(::Type{<:GajewskiParameters}) = :day

"""
    optical_density(p::GajewskiParameters, t, T)

Optical density at time `t` (days) at constant temperature `T` (°C): Eq. 1 with
`r = p.tpc(T)`.
"""
optical_density(p::GajewskiParameters, t, T) = delayed_logistic(t, p.Y0, p.K, p.tpc(T), p.d)

"""
    GajewskiRegime(name; as_coded_kelvin = false)

One of the three varying-temperature regimes of Gajewski et al. (2021), `:DLA` (Dry Low
Air), `:WHF` (Wet High Frog) or `:WLF` (Wet Low Frog): the authors' piecewise-linear fit to
the incubator logger (functions.R), repeated every day. Callable as `T(t)`, `t` in days, `T`
in °C. The segments are read from the compendium (sets `"regime_<name>"`).

With `as_coded_kelvin = true`, the last segment of `:WLF` (from 0.861 d) returns
`25 − 273.16` °C, which is what the authors' code feeds the kelvin-based curves (Ratkowsky
and Ikemoto) there: it omits the offset in that segment, so they see 25 K.
"""
struct GajewskiRegime
    name::Symbol
    starts::Vector{Float64}
    T0::Vector{Float64}
    slopes::Vector{Float64}
    as_coded_kelvin::Bool
end

function GajewskiRegime(name::Symbol; as_coded_kelvin::Bool = false)
    name in GAJEWSKI_REGIMES ||
        throw(ArgumentError("regime must be one of $(GAJEWSKI_REGIMES); got :$name"))
    p = parameter_set(:gajewski2021; set = "regime_" * string(name))
    n = count(k -> startswith(string(k), "start_"), keys(p))
    seg(k) = [p[Symbol(k, "_", i)] for i in 1:n]
    starts = seg("start")
    issorted(starts) && starts[1] == 0 || throw(ArgumentError("bad segments for $name"))
    return GajewskiRegime(name, starts, seg("T0"), seg("slope"), as_coded_kelvin)
end

function (g::GajewskiRegime)(t)
    s = t - floor(t)
    i = searchsortedlast(g.starts, s)
    g.as_coded_kelvin && g.name === :WLF && i == length(g.starts) &&
        return g.T0[i] - GAJEWSKI_PREDICTION_KELVIN
    return g.T0[i] + g.slopes[i] * s
end

"""
    GajewskiGrowthProblem(p::GajewskiParameters; temperature, tspan = (0.0, 14.0))

Logistic growth of optical density under a time-varying temperature (Gajewski et al. 2021,
Eq. 3), in days:

```
dx/dt = r(T(t)) x (1 − x/K) 1(t > d),   x(0) = Y0
```

`temperature` is a number (°C, constant) or any callable `T(t)`, such as a
[`GajewskiRegime`](@ref) or a [`SeasonalTemperature`](@ref). `solve(prob)` returns the
`ODESolution`; [`to_ode_problem`](@ref) gives the `ODEProblem`. At constant temperature the
solution is [`delayed_logistic`](@ref) (Eq. 1).
"""
struct GajewskiGrowthProblem{P<:GajewskiParameters,F}
    p::P
    temperature::F
    tspan::Tuple{Float64,Float64}
end

function GajewskiGrowthProblem(p::GajewskiParameters; temperature, tspan = (0.0, 14.0))
    T = temperature isa Real ? (let c = Float64(temperature); _ -> c; end) : temperature
    return GajewskiGrowthProblem(p, T, (Float64(tspan[1]), Float64(tspan[2])))
end

function _gajewski_rhs(u, prob::GajewskiGrowthProblem, t)
    p = prob.p
    t > p.d || return zero(u)
    return p.tpc(prob.temperature(t)) .* u .* (1 .- u ./ p.K)
end

to_ode_problem(prob::GajewskiGrowthProblem) =
    ODEProblem(_gajewski_rhs, [prob.p.Y0], prob.tspan, prob)

CommonSolve.solve(prob::GajewskiGrowthProblem, alg = Tsit5(); kwargs...) =
    solve(to_ode_problem(prob), alg; tstops = filter(t -> prob.tspan[1] < t < prob.tspan[2],
        [prob.p.d]), kwargs...)

"""
    gajewski_prediction(form, regime; as_coded = true, tspan = (0.0, 14.0), kwargs...)

The optical density predicted by the posterior medians of TPC fit `form` under varying
regime `regime` (`:DLA`, `:WHF` or `:WLF`), as in Gajewski et al. (2021) Fig. 4 but at the
medians rather than over posterior draws. Returns the `ODESolution` of
[`GajewskiGrowthProblem`](@ref).

With `as_coded = true` (the default), the kelvin-based curves use the offset 273.16 and the
`:WLF` regime keeps the authors' slip in its last segment ([`GajewskiRegime`](@ref)), as in
Predictions.R. With `as_coded = false` the offset is 273.15 and the slip is removed. Other
keywords go to [`GajewskiParameters`](@ref) (`ikemoto_delay`, `α`); solver tolerances are
`reltol = abstol = 1e-8`.
"""
function gajewski_prediction(form::Symbol, regime::Symbol; as_coded::Bool = true,
        tspan = (0.0, 14.0), kwargs...)
    kelvin = as_coded && form in GAJEWSKI_KELVIN_FORMS
    p = GajewskiParameters(:gajewski2021; form,
        kelvin_offset = kelvin ? GAJEWSKI_PREDICTION_KELVIN : KELVIN_OFFSET, kwargs...)
    T = GajewskiRegime(regime; as_coded_kelvin = kelvin)
    return solve(GajewskiGrowthProblem(p; temperature = T, tspan); reltol = 1e-8,
        abstol = 1e-8)
end

"""
    regime_growth_rate(p::GajewskiParameters, regime; step = 1e-3, mass = 0.95)
        -> (median, lo, hi)

The logistic growth rate `r(T(t))` over one day of a varying regime (a
[`GajewskiRegime`](@ref) or any callable `T(t)`): its median and its `mass` highest density
interval over a grid of `step` days on [0, 1). Gajewski et al. (2021) report this summary in
Appendix S1 Table B.2; they do not say how it was computed, and this is our reading of it.
"""
function regime_growth_rate(p::GajewskiParameters, regime; step::Real = 1e-3,
        mass::Real = 0.95)
    rs = [p.tpc(regime(t)) for t in 0:step:(1 - step)]
    lo, hi = _hdi(rs; mass)
    return (median = median(rs), lo = lo, hi = hi)
end

"""
    adjusted_limits(c; threshold = 0.01, range = (-10.0, 50.0), step = 0.01) -> (lo, hi)

The adjusted thermal limits of Gajewski et al. (2021): the lowest and highest temperatures
(°C) on a grid of `step` over `range` at which the curve `c` exceeds `threshold`. Their
code takes the grid point next to the first and last values above the threshold
(Predictions.R); this returns the crossing found by bisection instead.
"""
function adjusted_limits(c; threshold::Real = 0.01, range = (-10.0, 50.0), step::Real = 0.01)
    Ts = range[1]:step:range[2]
    above = findall(T -> c(T) > threshold, Ts)
    isempty(above) && throw(ArgumentError("the curve never exceeds $threshold"))
    function cross(a, b)          # c(a) ≤ threshold < c(b) or the reverse
        fa = c(a) > threshold
        for _ in 1:60
            m = (a + b) / 2
            (c(m) > threshold) == fa ? (a = m) : (b = m)
        end
        return (a + b) / 2
    end
    i, j = first(above), last(above)
    lo = i == 1 ? Ts[1] : cross(Ts[i - 1], Ts[i])
    hi = j == length(Ts) ? Ts[end] : cross(Ts[j], Ts[j + 1])
    return (lo, hi)
end

"""
    gajewski_summary(form; as_coded = true, kwargs...) -> NamedTuple

The curve summaries of Gajewski et al. (2021) Tables 1, 2 and B.2, computed from the curve
at the posterior medians of fit `form` ([`GajewskiParameters`](@ref); keywords go to it):
`Topt` and `r_max` ([`thermal_optimum`](@ref) on 0–45 °C), `Tmin_adj` and `Tmax_adj`
([`adjusted_limits`](@ref)), and `r_DLA`, `r_WHF`, `r_WLF`, the median rates over each regime
([`regime_growth_rate`](@ref)). With `as_coded = true` the regime rates of the kelvin-based
curves use the authors' prediction offset (273.16) and the `:WLF` slip, as in
[`gajewski_prediction`](@ref). The paper's values are posterior medians of each summary,
which the summary of the median curve only approximates.
"""
function gajewski_summary(form::Symbol; as_coded::Bool = true, kwargs...)
    c = GajewskiParameters(:gajewski2021; form, kwargs...).tpc.curve
    To = thermal_optimum(c; range = (0.0, 45.0))
    lo, hi = adjusted_limits(c)
    kelvin = as_coded && form in GAJEWSKI_KELVIN_FORMS
    pr = GajewskiParameters(:gajewski2021; form,
        kelvin_offset = kelvin ? GAJEWSKI_PREDICTION_KELVIN : KELVIN_OFFSET, kwargs...)
    rate(g) = regime_growth_rate(pr, GajewskiRegime(g; as_coded_kelvin = kelvin)).median
    return (Topt = To, r_max = c(To), Tmin_adj = lo, Tmax_adj = hi,
            r_DLA = rate(:DLA), r_WHF = rate(:WHF), r_WLF = rate(:WLF))
end
