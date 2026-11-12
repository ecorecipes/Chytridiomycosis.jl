# Thermal performance curves: rate as a function of temperature (°C), for pathogen isolates
# and hosts, plus seasonal temperature forcing. Spec: lab/wiki/models/thermal-performance.md
#
# The review gives no functional form (§9.1), so these are the standard families of thermal
# biology. The published primary data are mostly cardinal temperatures (limits and optima)
# and a few fitted curves (Kelly et al. 2024 polynomials). Those live in the compendium.

const BOLTZMANN_EV = 8.617333262e-5          # eV K⁻¹ (CODATA 2018)
const KELVIN_OFFSET = 273.15

"""
    ThermalPerformanceCurve

Supertype of the thermal performance curves. Every curve `c` is callable as `c(T)` with `T` in
°C, and returns a rate in the units of whatever it was built from. The curves are generic in
the number type, so they work with automatic differentiation.

Families: [`Briere`](@ref), [`CardinalTPC`](@ref), [`GaussianQuadratic`](@ref),
[`SharpeSchoolfield`](@ref), [`Ratkowsky`](@ref), [`PolynomialTPC`](@ref) and
[`InterpolatedTPC`](@ref). Use [`thermal_optimum`](@ref) to locate a curve's peak.
"""
abstract type ThermalPerformanceCurve end

"""
    Briere(; a, Tmin, Tmax, m = 2)
    Briere(; Tmin, Topt, Tmax, rmax)

The Brière curve (Brière et al. 1999), `r(T) = a T (T − Tmin) (Tmax − T)^(1/m)`, for
`Tmin < T < Tmax` and zero outside. `m = 2` is Brière-1, and a free `m` is Brière-2. It is
left-skewed, with an abrupt fall near `Tmax`.

The second form takes the cardinal temperatures and the peak rate, which is what most chytrid
papers report. The optimum satisfies `1/m = (Tmax − Topt)(1/Topt + 1/(Topt − Tmin))`, which
fixes `m`, and then `a` gives `r(Topt) = rmax`. It needs `0 ≤ Tmin < Topt < Tmax`.
"""
struct Briere{R<:Real} <: ThermalPerformanceCurve
    a::R
    Tmin::R
    Tmax::R
    m::R
end

function Briere(; Tmin, Tmax, a = nothing, m = nothing, Topt = nothing, rmax = nothing)
    Tmin < Tmax || throw(ArgumentError("need Tmin < Tmax"))
    if Topt === nothing
        a === nothing && throw(ArgumentError("give either a (and m) or Topt and rmax"))
        m = m === nothing ? 2 : m
        m > 0 || throw(ArgumentError("need m > 0"))
        return Briere(promote(a, Tmin, Tmax, m)...)
    end
    (a === nothing && m === nothing && rmax !== nothing) ||
        throw(ArgumentError("with Topt give rmax, and not a or m"))
    0 <= Tmin < Topt < Tmax || throw(ArgumentError("need 0 ≤ Tmin < Topt < Tmax"))
    m = 1 / ((Tmax - Topt) * (1 / Topt + 1 / (Topt - Tmin)))
    a = rmax / (Topt * (Topt - Tmin) * (Tmax - Topt)^(1 / m))
    return Briere(promote(a, Tmin, Tmax, m)...)
end

(c::Briere)(T) = c.Tmin < T < c.Tmax ? c.a * T * (T - c.Tmin) * (c.Tmax - T)^(1 / c.m) :
                 zero(c.a * T)

"""
    CardinalTPC(; Tmin, Topt, Tmax, rmax)

The cardinal temperature model with inflection (Rosso et al. 1993). It is set directly by the
three cardinal temperatures and the rate at the optimum:

    r(T) = rmax (T − Tmax)(T − Tmin)² /
           ((Topt − Tmin)[(Topt − Tmin)(T − Topt) − (Topt − Tmax)(Topt + Tmin − 2T)])

for `Tmin < T < Tmax`, and zero outside. Rosso et al. require `Topt > (Tmin + Tmax)/2`, and
the constructor enforces it. Bsal's reported cardinal temperatures (optimum 10–15 °C,
growth from 5 to 22 °C; Martel et al. 2013) fall outside that domain; use [`Briere`](@ref)
for them.
"""
struct CardinalTPC{R<:Real} <: ThermalPerformanceCurve
    Tmin::R
    Topt::R
    Tmax::R
    rmax::R
end

function CardinalTPC(; Tmin, Topt, Tmax, rmax)
    Tmin < Topt < Tmax || throw(ArgumentError("need Tmin < Topt < Tmax"))
    Topt > (Tmin + Tmax) / 2 ||
        throw(ArgumentError("the cardinal temperature model needs Topt > (Tmin + Tmax)/2"))
    return CardinalTPC(promote(Tmin, Topt, Tmax, rmax)...)
end

function (c::CardinalTPC)(T)
    c.Tmin < T < c.Tmax || return zero(c.rmax * T)
    d = (c.Topt - c.Tmin) * ((c.Topt - c.Tmin) * (T - c.Topt) -
                             (c.Topt - c.Tmax) * (c.Topt + c.Tmin - 2T))
    return c.rmax * (T - c.Tmax) * (T - c.Tmin)^2 / d
end

"""
    GaussianQuadratic(; rmax, Topt, σ, CTmax)

A Gaussian rise and a quadratic fall (Deutsch et al. 2008):
`rmax exp(−((T − Topt)/(2σ))²)` for `T ≤ Topt`, and
`rmax (1 − ((T − Topt)/(Topt − CTmax))²)` for `Topt < T < CTmax`. It is zero from `CTmax`
up. The curve is asymmetric, with a long cool tail and an abrupt upper decline.
"""
struct GaussianQuadratic{R<:Real} <: ThermalPerformanceCurve
    rmax::R
    Topt::R
    σ::R
    CTmax::R
end

function GaussianQuadratic(; rmax, Topt, σ, CTmax)
    Topt < CTmax || throw(ArgumentError("need Topt < CTmax"))
    σ > 0 || throw(ArgumentError("need σ > 0"))
    return GaussianQuadratic(promote(rmax, Topt, σ, CTmax)...)
end

function (c::GaussianQuadratic)(T)
    T <= c.Topt && return c.rmax * exp(-((T - c.Topt) / (2c.σ))^2)
    T < c.CTmax && return c.rmax * (1 - ((T - c.Topt) / (c.Topt - c.CTmax))^2)
    return zero(c.rmax * T)
end

"""
    SharpeSchoolfield(; r_ref, E, Eh, Th, El = nothing, Tl = nothing, Tref = 20)

The Sharpe–Schoolfield curve (Schoolfield, Sharpe & Magnuson 1981). It is an Arrhenius rise
with activation energy `E` (eV), divided by terms for enzyme inactivation at high, and
optionally low, temperatures:

    r(T) = r_ref exp(E/k (1/Tref − 1/T)) /
           (1 + exp(Eh/k (1/Th − 1/T)) + exp(El/k (1/T − 1/Tl)))

Temperatures are given in °C and converted to kelvin, and `k` is Boltzmann's constant in
eV K⁻¹. `Eh, El > 0`, and each inactivation term equals 1 at its half-inactivation
temperature `Th` or `Tl`. Leave out `El` and `Tl` to drop low-temperature inactivation.
`r_ref` is the uninhibited rate at `Tref`.
"""
struct SharpeSchoolfield{R<:Real} <: ThermalPerformanceCurve
    r_ref::R
    E::R
    Eh::R
    Th::R
    El::R
    Tl::R
    Tref::R
    low::Bool
end

function SharpeSchoolfield(; r_ref, E, Eh, Th, El = nothing, Tl = nothing, Tref = 20)
    (El === nothing) == (Tl === nothing) || throw(ArgumentError("give both El and Tl, or neither"))
    low = El !== nothing
    vals = promote(r_ref, E, Eh, Th, low ? El : 0, low ? Tl : 0, Tref)
    return SharpeSchoolfield(vals..., low)
end

function (c::SharpeSchoolfield)(T)
    K(x) = x + KELVIN_OFFSET
    k, TK = BOLTZMANN_EV, K(T)
    num = c.r_ref * exp(c.E / k * (1 / K(c.Tref) - 1 / TK))
    den = 1 + exp(c.Eh / k * (1 / K(c.Th) - 1 / TK))
    c.low && (den += exp(c.El / k * (1 / TK - 1 / K(c.Tl))))
    return num / den
end

"""
    Ratkowsky(; b, Tmin, c, Tmax)

The Ratkowsky model for microbial growth (Ratkowsky et al. 1983):
`√r = b (T − Tmin)(1 − exp(c (T − Tmax)))` for `Tmin < T < Tmax`, and zero outside.
"""
struct Ratkowsky{R<:Real} <: ThermalPerformanceCurve
    b::R
    Tmin::R
    c::R
    Tmax::R
end

function Ratkowsky(; b, Tmin, c, Tmax)
    Tmin < Tmax || throw(ArgumentError("need Tmin < Tmax"))
    return Ratkowsky(promote(b, Tmin, c, Tmax)...)
end

(c::Ratkowsky)(T) = c.Tmin < T < c.Tmax ?
    (c.b * (T - c.Tmin) * (1 - exp(c.c * (T - c.Tmax))))^2 : zero(c.b * T)

"""
    PolynomialTPC(coefficients; Tmin = -Inf, Tmax = Inf, nonnegative = true)

A polynomial in temperature, `Σ coefficients[k+1] T^k`. This is the form Kelly et al. (2024)
fitted to the growth rates of each Bsal isolate. Outside `[Tmin, Tmax]`, the range the
polynomial was fitted on, it returns zero rather than extrapolating. With
`nonnegative = true`, negative values are set to zero.
"""
struct PolynomialTPC{R<:Real} <: ThermalPerformanceCurve
    coefficients::Vector{R}
    Tmin::Float64
    Tmax::Float64
    nonnegative::Bool
end

PolynomialTPC(coefficients::AbstractVector{<:Real}; Tmin::Real = -Inf, Tmax::Real = Inf,
              nonnegative::Bool = true) =
    PolynomialTPC(collect(float.(coefficients)), Float64(Tmin), Float64(Tmax), nonnegative)

function (c::PolynomialTPC)(T)
    c.Tmin <= T <= c.Tmax || return zero(first(c.coefficients) * T)
    y = evalpoly(T, c.coefficients)
    return c.nonnegative ? max(y, zero(y)) : y
end

"""
    InterpolatedTPC(temperatures, values; extrapolation = :zero)

A rate measured at a few temperatures, interpolated linearly between them. Outside the
measured range it is zero (`extrapolation = :zero`) or the nearest measured value
(`:flat`). Use it when the data are too sparse to choose a family, as with Woodhams et al.'s
(2008) zoospore rates at four temperatures (`thermal_data(:woodhams2008, :mZ; set_prefix =
"table1")`).
"""
struct InterpolatedTPC{R<:Real} <: ThermalPerformanceCurve
    temperatures::Vector{Float64}
    values::Vector{R}
    extrapolation::Symbol
end

function InterpolatedTPC(temperatures::AbstractVector{<:Real}, values::AbstractVector{<:Real};
                         extrapolation::Symbol = :zero)
    length(temperatures) == length(values) >= 2 ||
        throw(ArgumentError("need at least two temperatures, one value each"))
    extrapolation in (:zero, :flat) || throw(ArgumentError("extrapolation is :zero or :flat"))
    i = sortperm(temperatures)
    allunique(temperatures) || throw(ArgumentError("temperatures must be distinct"))
    return InterpolatedTPC(Float64.(temperatures[i]), collect(values[i]), extrapolation)
end

function (c::InterpolatedTPC)(T)
    ts, vs = c.temperatures, c.values
    if T < ts[1] || T > ts[end]
        c.extrapolation === :zero && return zero(vs[1] * T)
        return T < ts[1] ? vs[1] * one(T) : vs[end] * one(T)
    end
    j = clamp(searchsortedlast(ts, T), 1, length(ts) - 1)
    w = (T - ts[j]) / (ts[j + 1] - ts[j])
    return (1 - w) * vs[j] + w * vs[j + 1]
end

"""
    thermal_optimum(c::ThermalPerformanceCurve; range = (-5, 45), step = 0.01) -> Float64

The temperature (°C) at which `c` peaks within `range`. It searches a grid of spacing `step`
and then refines with a golden-section search around the best grid point.
"""
function thermal_optimum(c::ThermalPerformanceCurve; range = (-5.0, 45.0), step::Real = 0.01)
    lo, hi = float.(range)
    Ts = lo:step:hi
    i = argmax(c.(Ts))
    a, b = max(lo, Ts[i] - step), min(hi, Ts[i] + step)
    φ = (sqrt(5) - 1) / 2
    for _ in 1:60
        x1, x2 = b - φ * (b - a), a + φ * (b - a)
        c(x1) >= c(x2) ? (b = x2) : (a = x1)
    end
    return (a + b) / 2
end

"""
    IsolateTPC(pathogen, isolate, lineage, trait, curve, source)
    IsolateTPC(key::Symbol, isolate::AbstractString)

A thermal performance curve for one trait of one pathogen isolate, with its provenance. It is
callable like its `curve`.

- `pathogen`: `"Bd"` or `"Bsal"`.
- `trait`: for example `"growth rate r (day⁻¹)"`.
- `source`: the compendium key.

The keyed form builds the curve from the compendium. So far this means
`IsolateTPC(:kelly2024, "Eifel2019")`: Kelly et al.'s (2024) fitted polynomial for that
isolate, on the 0–30 °C range they tested.
"""
struct IsolateTPC{C<:ThermalPerformanceCurve}
    pathogen::String
    isolate::String
    lineage::String
    trait::String
    curve::C
    source::Symbol
end

(c::IsolateTPC)(T) = c.curve(T)

function IsolateTPC(key::Symbol, isolate::AbstractString)
    key === :kelly2024 ||
        throw(ArgumentError("no isolate curves for $key; the compendium has :kelly2024"))
    recs = parameter_records(key)
    any(r -> r.set == isolate, recs) ||
        throw(ArgumentError("$key has no isolate \"$isolate\""))
    p = parameter_set(key; set = isolate)
    coefs = [getfield(p, Symbol("c$k")) for k in 0:4 if haskey(p, Symbol("c$k"))]
    sys = recs[findfirst(r -> r.set == isolate, recs)].system
    return IsolateTPC(get(sys, "pathogen", "Bsal"), String(isolate), get(sys, "lineage", ""),
        "growth rate r (day⁻¹), total counts (GrowthCurver)",
        PolynomialTPC(coefs; Tmin = p.T_tested_min, Tmax = p.T_tested_max), key)
end

"""
    HostTPC(species, trait, curve, source)

A thermal performance curve for a host trait, such as performance, immune function or
survival of infection. It is kept separate from pathogen curves ([`IsolateTPC`](@ref)) so
that thermal mismatch can be built from the two. It is callable like its `curve`. No host
curve has yet been verified for the compendium.
"""
struct HostTPC{C<:ThermalPerformanceCurve}
    species::String
    trait::String
    curve::C
    source::Symbol
end

(c::HostTPC)(T) = c.curve(T)

"""
    thermal_data(key::Symbol, symbol::Symbol; set_prefix = "") -> (temperatures, values)

The values of compendium symbol `symbol` measured at several temperatures. It reads the
records whose `system.temperature_C` is a number and whose parameter set starts with
`set_prefix`, sorted by temperature. For example, `thermal_data(:woodhams2008, :mZ;
set_prefix = "table1")` gives the encystment rates of Table 1. Pass them to
[`InterpolatedTPC`](@ref) or to a fit.
"""
function thermal_data(key::Symbol, symbol::Symbol; set_prefix::AbstractString = "")
    recs = [r for r in parameter_records(key) if r.symbol === symbol &&
            startswith(r.set, set_prefix) && get(r.system, "temperature_C", "") isa Real]
    isempty(recs) && throw(ArgumentError("$key has no $symbol records at stated temperatures"))
    Ts = [Float64(r.system["temperature_C"]) for r in recs]
    allunique(Ts) || throw(ArgumentError("$key.$symbol: several records at one temperature; " *
        "narrow them with set_prefix"))
    i = sortperm(Ts)
    return Ts[i], [r.value for r in recs][i]
end

"""
    SeasonalTemperature(; mean, amplitude, peak, period = 365)

A sinusoidal temperature (°C) through the year:
`T(t) = mean + amplitude cos(2π (t − peak)/period)`. The phase is stated explicitly as the
time `peak` of the warmest point, in the same time unit as `period`. Callable as `T(t)`.
"""
Base.@kwdef struct SeasonalTemperature{R<:Real}
    mean::R
    amplitude::R
    peak::R
    period::R = 365.0
end

SeasonalTemperature(mean, amplitude, peak, period) =
    SeasonalTemperature(promote(mean, amplitude, peak, period)...)

(s::SeasonalTemperature)(t) = s.mean + s.amplitude * cos(2π * (t - s.peak) / s.period)

"""
    temperature_dependent(rate, T) -> t -> rate(T(t))

Compose a rate function of temperature (a [`ThermalPerformanceCurve`](@ref), an
[`IsolateTPC`](@ref) or any callable) with a temperature series `T(t)`, such as a
[`SeasonalTemperature`](@ref) or an interpolant of observed temperatures. The result is a
time-varying rate to use inside an ODE, SDE or jump right-hand side.
"""
temperature_dependent(rate, T) = t -> rate(T(t))
