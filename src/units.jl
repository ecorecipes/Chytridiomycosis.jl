# Units layer: typed pathogen loads with explicit conversions, and explicit time units for
# rates. Spec: lab/wiki/parameters/units.md
#
# A light internal layer rather than Unitful.jl. Most of the package's units are domain
# counts (sporangia, ZE, ITS1 copies, frogs, an arbitrary pool volume) that Unitful does not
# know and that should not convert into one another silently. Time is the only physical
# dimension the models share, and model structs hold plain floats in a declared time unit.

"""
    PathogenLoad

Supertype of the four load types: [`ZEPerSwab`](@ref), [`ITS1CopiesPerSwab`](@ref),
[`Zoospores`](@ref) and [`SporangiaPerHost`](@ref). Each stores log10 of the value, so loads
spanning many orders of magnitude keep their precision.

- Construct a load from the natural value, for example `ZEPerSwab(1e4)`, or from log10 with
  `ZEPerSwab(; log10 = 4)`.
- Read it back with [`load_value`](@ref) or `log10(x)`.
- Loads of the same type compare (`<`, `==`, `isapprox`). Comparing loads of different types
  throws an `ArgumentError`, because there is no basis to compare them on.
- Change basis only with [`convert_load`](@ref), which takes the assumptions as keyword
  arguments.
"""
abstract type PathogenLoad end

_log10load(x::Real) = x >= 0 ? log10(float(x)) :
                      throw(DomainError(x, "a pathogen load cannot be negative"))

"""
    ZEPerSwab(value)
    ZEPerSwab(; log10)

Zoospore equivalents per swab. This is the qPCR estimate, scaled against standards made from
counted zoospores of one Bd strain. It is the unit of most field data and of the
"10,000 ZE" crash load.
"""
struct ZEPerSwab <: PathogenLoad
    log10value::Float64
    ZEPerSwab(x::Real) = new(_log10load(x))
    ZEPerSwab(; log10::Real) = new(Float64(log10))
end

"""
    ITS1CopiesPerSwab(value)
    ITS1CopiesPerSwab(; log10)

ITS1 gene copies per swab, from qPCR scaled against amplicon standards. The number of copies
per zoospore varies among Bd strains (Longo et al. 2013), so converting to
[`ZEPerSwab`](@ref) needs `copies_per_zoospore`.
"""
struct ITS1CopiesPerSwab <: PathogenLoad
    log10value::Float64
    ITS1CopiesPerSwab(x::Real) = new(_log10load(x))
    ITS1CopiesPerSwab(; log10::Real) = new(Float64(log10))
end

"""
    Zoospores(value)
    Zoospores(; log10)

A count of actual zoospores, such as zoospores on one host, in a culture well or in the pool
of the Briggs model. The compartment it counts is the model's to state. It is not the
qPCR-based [`ZEPerSwab`](@ref).
"""
struct Zoospores <: PathogenLoad
    log10value::Float64
    Zoospores(x::Real) = new(_log10load(x))
    Zoospores(; log10::Real) = new(Float64(log10))
end

"""
    SporangiaPerHost(value)
    SporangiaPerHost(; log10)

Sporangia (zoosporangia) on one host. This is the state variable of the Briggs et al. (2010)
load model, whose lethal load `Smax` is in these units. Converting it to or from
[`ZEPerSwab`](@ref) has no published basis.
"""
struct SporangiaPerHost <: PathogenLoad
    log10value::Float64
    SporangiaPerHost(x::Real) = new(_log10load(x))
    SporangiaPerHost(; log10::Real) = new(Float64(log10))
end

"""
    load_value(x::PathogenLoad) -> Float64

The natural-scale value of load `x`, in the units of its type. `log10(x)` gives log10 of the
value.
"""
load_value(x::PathogenLoad) = exp10(x.log10value)
Base.log10(x::PathogenLoad) = x.log10value

Base.show(io::IO, x::PathogenLoad) = print(io, nameof(typeof(x)), "(", load_value(x), ")")

_incomparable(a, b) = throw(ArgumentError("cannot compare a $(nameof(typeof(a))) with a " *
    "$(nameof(typeof(b))): convert one with convert_load, stating the assumptions"))

Base.isless(a::L, b::L) where {L<:PathogenLoad} = isless(a.log10value, b.log10value)
Base.:(==)(a::L, b::L) where {L<:PathogenLoad} = a.log10value == b.log10value
Base.isapprox(a::L, b::L; kw...) where {L<:PathogenLoad} =
    isapprox(load_value(a), load_value(b); kw...)
Base.isless(a::PathogenLoad, b::PathogenLoad) = _incomparable(a, b)
Base.:(==)(a::PathogenLoad, b::PathogenLoad) = _incomparable(a, b)
Base.isapprox(a::PathogenLoad, b::PathogenLoad; kw...) = _incomparable(a, b)
Base.hash(x::L, h::UInt) where {L<:PathogenLoad} = hash(x.log10value, hash(L, h))

"""
    convert_load(x::PathogenLoad, ::Type{<:PathogenLoad}; assumptions...)

Convert load `x` to another basis. Every assumption the conversion rests on is a required
keyword argument, so it is visible at the call site. Conversions are exact inverses of each
other.

| from → to | keywords | rule |
|---|---|---|
| `ITS1CopiesPerSwab` ↔ `ZEPerSwab` | `copies_per_zoospore`, `extraction_factor = 1` | ZE = copies × `extraction_factor` / `copies_per_zoospore` |
| `ZEPerSwab` ↔ `Zoospores` | `swab_fraction` | zoospores on the host = ZE / `swab_fraction` |
| `ZEPerSwab` ↔ `SporangiaPerHost` | `ze_per_sporangium` | sporangia = ZE / `ze_per_sporangium`; **warns** |

- `copies_per_zoospore` is the ITS1 copy number of the strain used for the standards. It
  ranges from 10 to 144 among strains (`parameter_set(:longo2013)`).
- `extraction_factor` carries any protocol rescaling between labs, such as dilution of the
  extract.
- `swab_fraction` is the fraction of the zoospores on a host that one swab recovers.
- ZE to sporangia has no published basis (review §2.4). That conversion exists so a model can
  state its assumption, and it warns each time it is called.

Converting to the same type returns `x`. Other pairs throw an `ArgumentError`; go through
`ZEPerSwab` in two explicit steps.
"""
convert_load(x::L, ::Type{L}) where {L<:PathogenLoad} = x
convert_load(x::PathogenLoad, ::Type{L}; kw...) where {L<:PathogenLoad} =
    throw(ArgumentError("no direct conversion from $(nameof(typeof(x))) to $(nameof(L)); " *
        "convert through ZEPerSwab in two explicit steps"))

_positive(name, v) = v > 0 ? log10(float(v)) :
                     throw(ArgumentError("$name must be positive, got $v"))

convert_load(x::ITS1CopiesPerSwab, ::Type{ZEPerSwab}; copies_per_zoospore::Real,
             extraction_factor::Real = 1) =
    ZEPerSwab(; log10 = x.log10value - _positive("copies_per_zoospore", copies_per_zoospore) +
        _positive("extraction_factor", extraction_factor))
convert_load(x::ZEPerSwab, ::Type{ITS1CopiesPerSwab}; copies_per_zoospore::Real,
             extraction_factor::Real = 1) =
    ITS1CopiesPerSwab(; log10 = x.log10value + _positive("copies_per_zoospore", copies_per_zoospore) -
        _positive("extraction_factor", extraction_factor))

function _fraction(v)
    0 < v <= 1 || throw(ArgumentError("swab_fraction must be in (0, 1], got $v"))
    return log10(float(v))
end
convert_load(x::ZEPerSwab, ::Type{Zoospores}; swab_fraction::Real) =
    Zoospores(; log10 = x.log10value - _fraction(swab_fraction))
convert_load(x::Zoospores, ::Type{ZEPerSwab}; swab_fraction::Real) =
    ZEPerSwab(; log10 = x.log10value + _fraction(swab_fraction))

const _SPORANGIA_WARNING = "ZE per swab and sporangia per host have no published conversion " *
    "(review §2.4); the result rests entirely on the assumed ze_per_sporangium"
function convert_load(x::ZEPerSwab, ::Type{SporangiaPerHost}; ze_per_sporangium::Real)
    @warn _SPORANGIA_WARNING ze_per_sporangium
    return SporangiaPerHost(; log10 = x.log10value - _positive("ze_per_sporangium", ze_per_sporangium))
end
function convert_load(x::SporangiaPerHost, ::Type{ZEPerSwab}; ze_per_sporangium::Real)
    @warn _SPORANGIA_WARNING ze_per_sporangium
    return ZEPerSwab(; log10 = x.log10value + _positive("ze_per_sporangium", ze_per_sporangium))
end

# Units strings in the compendium that name a load basis (whitespace-normalised).
const LOAD_UNITS = Dict("ZE swab^-1" => ZEPerSwab, "ITS1 copies swab^-1" => ITS1CopiesPerSwab,
    "zoospores" => Zoospores, "sporangia host^-1" => SporangiaPerHost)

"""
    parameter_load(key::Symbol, symbol::Symbol; set = "") -> PathogenLoad

The compendium value `symbol` of source paper `key`, typed by the load basis its `units`
name. Recognised units are `"ZE swab^-1"`, `"ITS1 copies swab^-1"`, `"zoospores"` and
`"sporangia host^-1"`. Throws an `ArgumentError` if the record is not a load.

For example, `parameter_load(:briggs2010, :Smax)` is a [`SporangiaPerHost`](@ref) and
`parameter_load(:kinney2011, :crash_load)` is a [`ZEPerSwab`](@ref). They cannot be compared
without an explicit conversion.
"""
function parameter_load(key::Symbol, symbol::Symbol; set::AbstractString = "")
    recs = parameter_records(key)
    i = findfirst(r -> r.symbol === symbol && r.set == set, recs)
    i === nothing && throw(ArgumentError("$key has no record $symbol in set \"$set\""))
    r = recs[i]
    u = join(split(r.units), " ")
    haskey(LOAD_UNITS, u) ||
        throw(ArgumentError("$key.$symbol has units \"$(r.units)\", which name no load basis"))
    return LOAD_UNITS[u](r.value)
end

# ---------------------------------------------------------------------------------------
# Time units

const TIME_UNIT_NAMES = Dict("h" => :hour, "hr" => :hour, "hour" => :hour, "hours" => :hour,
    "d" => :day, "day" => :day, "days" => :day, "wk" => :week, "week" => :week,
    "weeks" => :week, "yr" => :year, "year" => :year, "years" => :year)

function _days(u::Symbol, days_per_year::Real)
    u === :hour && return 1 / 24
    u === :day && return 1.0
    u === :week && return 7.0
    u === :year && return float(days_per_year)
    throw(ArgumentError("unknown time unit :$u; use :hour, :day, :week or :year"))
end

"""
    convert_rate(r, (from => to)::Pair{Symbol,Symbol}; days_per_year = 365)

Convert a rate `r` per `from` into a rate per `to`. Time units are `:hour`, `:day`, `:week`
and `:year`. A year is `days_per_year` days.

The default year is 365 days. Pass the source paper's own convention when it differs. For
example, Schmidt et al. (2017) use 52 weeks (364 days): their `e = 364 yr⁻¹` is one day.

```julia
convert_rate(52.0, :year => :day; days_per_year = 364)   # 1/7 per day
```
"""
convert_rate(r, units::Pair{Symbol,Symbol}; days_per_year::Real = 365) =
    r * _days(units.second, days_per_year) / _days(units.first, days_per_year)

# The time unit in a compendium units string and its exponent, or nothing.
function _time_part(units::AbstractString)
    found = nothing
    for tok in split(units)
        m = match(r"^([A-Za-z]+)(?:\^(-?\d+))?$", tok)
        m === nothing && continue
        haskey(TIME_UNIT_NAMES, m.captures[1]) || continue
        found === nothing ||
            throw(ArgumentError("units \"$units\" contain more than one time unit"))
        found = (TIME_UNIT_NAMES[m.captures[1]],
                 m.captures[2] === nothing ? 1 : parse(Int, m.captures[2]))
    end
    return found
end

"""
    in_time_unit(r::ParameterRecord, to::Symbol; days_per_year = 365) -> Float64

The value of compendium record `r` with its time component expressed in `to` (`:hour`,
`:day`, `:week` or `:year`). It reads the time unit from the record's `units`:

- a rate (`"yr^-1"`, `"ha yr^-1"`, `"zoospore^-1 yr^-1"`) is converted as a rate;
- a duration (`"day"`) is converted as a time.

Other unit factors are left as they are. Throws an `ArgumentError` when the units contain no
time unit. Use it at the boundary where a model reads a published value into its own time
unit ([`time_unit`](@ref)).
"""
function in_time_unit(r::ParameterRecord, to::Symbol; days_per_year::Real = 365)
    tp = _time_part(r.units)
    tp === nothing &&
        throw(ArgumentError("$(r.key).$(r.symbol): units \"$(r.units)\" contain no time unit"))
    from, k = tp
    return r.value * (_days(from, days_per_year) / _days(to, days_per_year))^k
end

"""
    survival_to_hazard(p, Δt)

The constant hazard rate (per unit of `Δt`) that gives survival probability `p` over an
interval `Δt`: `−log(p)/Δt`. A hazard computed this way from a published survival is a
derived (`D`) value, and its compendium record should say so in `notes`.
"""
function survival_to_hazard(p::Real, Δt::Real)
    0 < p <= 1 || throw(DomainError(p, "survival probability must be in (0, 1]"))
    Δt > 0 || throw(DomainError(Δt, "the interval must be positive"))
    return -log(p) / Δt
end

"""
    hazard_to_survival(h, Δt)

The probability of surviving an interval `Δt` at constant hazard `h`: `exp(−hΔt)`. It is the
inverse of [`survival_to_hazard`](@ref).
"""
hazard_to_survival(h::Real, Δt::Real) = exp(-h * Δt)

"""
    time_unit(p) -> Symbol
    time_unit(::Type)

The time unit in which a model's parameter struct holds its rates: `:hour`, `:day`, `:week`
or `:year`. Declared so far:

- [`SchmidtParameters`](@ref): `:year`
- [`BriggsParameters`](@ref) and [`BriggsSeasonalParameters`](@ref): `:day`
- [`MitchellParameters`](@ref): `:year`
- [`WoodhamsParameters`](@ref): `:hour`
- [`RumschlagParameters`](@ref): `:day`
- [`RefugiaParameters`](@ref): `:week`
- [`BozzutoParameters`](@ref): `:year`

A new model declares its own with
`Chytridiomycosis.time_unit(::Type{<:MyParameters}) = :day`.
"""
time_unit(p) = time_unit(typeof(p))
time_unit(T::Type) = throw(ArgumentError("$T declares no time unit"))
time_unit(::Type{<:SchmidtParameters}) = :year
time_unit(::Type{<:BriggsParameters}) = :day
time_unit(::Type{<:BriggsSeasonalParameters}) = :day
time_unit(::Type{<:MitchellParameters}) = :year
