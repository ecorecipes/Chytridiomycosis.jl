# Units and loads

```@meta
CurrentModule = Chytridiomycosis
```

The field reports pathogen load on several bases that are easy to conflate: zoospore
equivalents per swab, ITS1 copies per swab, actual zoospores, and sporangia per host. The
"10,000 ZE" crash load and the Briggs model's lethal load of 10⁴ sporangia are different
quantities. The package gives each basis its own type and converts between them only when
the assumptions are stated:

```julia
smax  = parameter_load(:briggs2010, :Smax)        # SporangiaPerHost(10000.0)
crash = parameter_load(:kinney2011, :crash_load)  # ZEPerSwab(10000.0)
smax == crash                                     # ArgumentError: different bases

q = parameter_set(:longo2013)                     # 10 to 144 ITS1 copies per zoospore
convert_load(ITS1CopiesPerSwab(5000), ZEPerSwab; copies_per_zoospore = q.copies_per_zoospore_max)
```

## How models use the layer

Model structs keep plain floats, so the ODE right-hand sides are unchanged. The layer is used
at the boundaries:

- **Declare the model's time unit.** Use `Chytridiomycosis.time_unit(::Type{<:MyParameters}) = :day`.
  The published models declare theirs.
- **Read published values into that unit.** Use
  `in_time_unit(record, time_unit(MyParameters))` instead of hand-written factors of 365 or
  24. State `days_per_year` where the paper's year is not 365 days (Schmidt et al. 2017 use
  52 weeks).
- **Type every load a model exposes.** Return loads as [`SporangiaPerHost`](@ref),
  [`ZEPerSwab`](@ref) and so on, and compare model loads with field data only after an
  explicit [`convert_load`](@ref).
- **Turn survival probabilities into hazards** with [`survival_to_hazard`](@ref), and record
  the result as a derived (`D`) value.

Why not Unitful.jl: most of the package's units are domain counts (ZE, ITS1 copies, sporangia,
frogs, an arbitrary pool volume) that Unitful does not define and that should never convert
into one another implicitly. Time is the only physical dimension the models share. A small
internal layer covers both, with no new dependency.

## Loads

```@docs
PathogenLoad
ZEPerSwab
ITS1CopiesPerSwab
Zoospores
SporangiaPerHost
load_value
convert_load
parameter_load
```

## Time units

```@docs
convert_rate
in_time_unit
survival_to_hazard
hazard_to_survival
time_unit
```
