# Validation-target registry: published numbers and outcomes that the implementations must
# reproduce, each with its source, the call that reproduces it, a tolerance and a status.
# Registry: data/validation/targets.toml. Spec: lab/wiki/validation/targets.md

const VALIDATION_FILE = normpath(joinpath(@__DIR__, "..", "data", "validation", "targets.toml"))
const TARGET_STATUSES = (:reproduced, :qualitative, :discrepant, :pending)

"""
    ValidationTarget

One published reproduction target, read from `data/validation/targets.toml`. Fields:

- `id`, `model` (the wiki model page or model name), `source` (a `references.tsv` key) and
  `location` (where in the paper);
- `quantity` and `published`, the value or outcome the paper reports;
- `call`, the package code that reproduces it, evaluated in the package module;
- `tolerance`, a `Dict` with one of `digits` (equal after rounding), `atol`, `rtol`, or
  `lower` and `upper` (strict interval); an empty `Dict` means exact equality;
- `status`: one of
  - `:reproduced`: the published number is matched within tolerance;
  - `:qualitative`: a published statement is matched as a condition;
  - `:discrepant`: the paper's number is not reproduced; `regression` holds ours;
  - `:pending`: no implementation yet;
- `regression`, our value and its tolerance as a `Dict`, or `nothing`;
- `notes`.

See [`validation_targets`](@ref) and [`check_target`](@ref).
"""
struct ValidationTarget
    id::String
    model::String
    source::Symbol
    location::String
    quantity::String
    published::Any
    call::String
    tolerance::Dict{String,Any}
    status::Symbol
    regression::Union{Nothing,Dict{String,Any}}
    notes::String
end

function _target(d::Dict{String,Any})
    st = Symbol(d["status"])
    st in TARGET_STATUSES ||
        throw(ArgumentError("target $(d["id"]): status $st not in $TARGET_STATUSES"))
    call = strip(get(d, "call", ""))
    st === :pending || !isempty(call) ||
        throw(ArgumentError("target $(d["id"]): only pending targets may have no call"))
    reg = get(d, "regression", nothing)
    st === :discrepant && reg === nothing &&
        throw(ArgumentError("target $(d["id"]): a discrepant target needs a regression value"))
    src = d["source"]
    return ValidationTarget(d["id"], get(d, "model", ""), Symbol(src["ref"]),
        get(src, "location", ""), d["quantity"], get(d, "published", nothing), call,
        get(d, "tolerance", Dict{String,Any}()), st, reg, get(d, "notes", ""))
end

"""
    validation_targets(; status = nothing, source = nothing) -> Vector{ValidationTarget}

All registered reproduction targets, or those with the given `status` (for example
`:discrepant`) or `source` key (for example `:schmidt2017`).
"""
function validation_targets(; status::Union{Nothing,Symbol} = nothing,
                            source::Union{Nothing,Symbol} = nothing)
    ts = [_target(d) for d in TOML.parsefile(VALIDATION_FILE)["target"]]
    allunique(t.id for t in ts) || throw(ArgumentError("duplicate target ids in $VALIDATION_FILE"))
    status === nothing || filter!(t -> t.status === status, ts)
    source === nothing || filter!(t -> t.source === source, ts)
    return ts
end

"""
    validation_target(id::AbstractString) -> ValidationTarget

The registered target with this `id`.
"""
function validation_target(id::AbstractString)
    i = findfirst(t -> t.id == id, validation_targets())
    i === nothing && throw(ArgumentError("no validation target \"$id\""))
    return validation_targets()[i]
end

"""
    reproduce(t::ValidationTarget)

Run the target's `call` in the package module and return what it computes. A `:pending`
target has no call and throws an `ArgumentError`.
"""
function reproduce(t::ValidationTarget)
    t.status === :pending && throw(ArgumentError("target $(t.id) is pending: no implementation"))
    return Core.eval(@__MODULE__, Meta.parse("begin\n" * t.call * "\nend"))
end

function _within(value, target, tol::AbstractDict)
    haskey(tol, "lower") && return tol["lower"] < value < tol["upper"]
    haskey(tol, "digits") && return round(value; digits = tol["digits"]) == target
    haskey(tol, "atol") && return abs(value - target) <= tol["atol"]
    haskey(tol, "rtol") && return isapprox(value, target; rtol = tol["rtol"])
    return value == target
end

"""
    check_target(t::ValidationTarget) -> (; value, published, regression)

Reproduce target `t` and compare the result.

- `published` is `true` when the value matches the paper within `t.tolerance`.
- `regression` is `true` when it matches our recorded value, and `nothing` when there is no
  regression value.

For a `:discrepant` target `published` is expected to be `false`.
"""
function check_target(t::ValidationTarget)
    v = reproduce(t)
    pub = _within(v, t.published, t.tolerance)
    reg = t.regression === nothing ? nothing :
          _within(v, t.regression["value"], filter(p -> p.first != "value", t.regression))
    return (; value = v, published = pub, regression = reg)
end
