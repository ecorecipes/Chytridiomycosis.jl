# Parameter compendium: one TOML file per source paper in data/compendium/.
# Schema: lab/wiki/parameters/compendium-schema.md

const COMPENDIUM_DIR = normpath(joinpath(@__DIR__, "..", "data", "compendium"))
const PROVENANCE_CODES = ("M", "F", "E", "A", "D")

"""
    ParameterRecord

One compendium entry: a value with its units, provenance code (`M` measured, `F` fitted,
`E` elicited, `A` assumed, `D` derived), role, system, source location and whether it has
been checked against the primary source.
"""
struct ParameterRecord
    key::Symbol
    symbol::Symbol
    quantity::String
    value::Float64
    units::String
    provenance::String
    role::String
    system::Dict{String,Any}
    location::String
    verified::Bool
    notes::String
    set::String
end

"""
    compendium_sources() -> Vector{Symbol}

Keys of all source papers in the compendium.
"""
compendium_sources() =
    sort!([Symbol(splitext(f)[1]) for f in readdir(COMPENDIUM_DIR) if endswith(f, ".toml")])

"""
    parameter_records(key::Symbol) -> Vector{ParameterRecord}

All compendium records for source paper `key` (e.g. `:schmidt2017`).
"""
function parameter_records(key::Symbol)
    path = joinpath(COMPENDIUM_DIR, string(key) * ".toml")
    isfile(path) || throw(ArgumentError("no compendium file for $key in $COMPENDIUM_DIR"))
    raw = TOML.parsefile(path)
    recs = ParameterRecord[]
    for r in get(raw, "parameter", Any[])
        prov = r["provenance"]
        prov in PROVENANCE_CODES ||
            throw(ArgumentError("$key.$(r["symbol"]): provenance $prov not in $PROVENANCE_CODES"))
        push!(recs, ParameterRecord(key, Symbol(r["symbol"]), r["quantity"], Float64(r["value"]),
            r["units"], prov, get(r, "role", ""), get(r, "system", Dict{String,Any}()),
            get(r, "location", ""), get(r, "verified", false), get(r, "notes", ""),
            get(r, "set", "")))
    end
    length(unique((r.set, r.symbol) for r in recs)) == length(recs) ||
        throw(ArgumentError("$key: duplicate symbols within a parameter set"))
    return recs
end

"""
    parameter_set(key::Symbol; set = "") -> NamedTuple

The values for source paper `key`, keyed by symbol, in the units recorded in the
compendium (each model documents the units it expects). Records without a `set` field form
the base set; `set = "fig3b"` returns the base set with that named set's records added or
substituted (for papers that use different values in different figures).
"""
function parameter_set(key::Symbol; set::AbstractString = "")
    recs = parameter_records(key)
    isempty(set) || any(r -> r.set == set, recs) ||
        throw(ArgumentError("$key has no parameter set \"$set\""))
    vals = Dict{Symbol,Float64}()
    order = Symbol[]
    for pass in ("", set), r in recs
        r.set == pass || continue
        haskey(vals, r.symbol) || push!(order, r.symbol)
        vals[r.symbol] = r.value
    end
    return NamedTuple{Tuple(order)}(Tuple(vals[s] for s in order))
end
