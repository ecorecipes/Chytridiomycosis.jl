# Drawert, Griesemer, Petzold & Briggs (2017) J R Soc Interface 14:20170480,
# "Using stochastic epidemiological models to evaluate conservation strategies for endangered
# amphibians". doi:10.1098/rsif.2017.0480
# The Briggs et al. (2010) load model made individual-based and demographically stochastic:
# integer sporangial loads on each frog, integer zoospores in one pool, and one-time
# interventions (cleaning frogs, culling frogs, cleaning the environment).
# Checked against the paper (§3, Tables 1–2), the SI text (figshare 10.6084/m9.figshare.5328160;
# §1, Table S3, eqs 10–12, Table S9) and the authors' C code and parameter table (figshare
# 10.6084/m9.figshare.5328127; WithinHostModel.c, FrogFunctions.c, FrogAdditions.c,
# issa/core.c, determ_frogs.c). Units: days; sporangia per frog; zoospores per pool.
# Wiki: lab/wiki/models/bd-drawert2017-stochastic-ibm.md
#
# Deterministic model (eqs 3.1–3.2; ESM eqs 1–2), for live frogs i:
#   dS_i/dt = ν(S_i) γ Z + ν(S_i) η f S_i − σ(S_i) S_i,   frog i dies when S_i > Smax
#   dZ/dt   = Σ_i {η(1 − f) S_i − γ Z} − μ Z
# with ν(S) = ν0 (Baseline, Sigma) or ν0 exp(−ν1 S) (Nu), σ(S) = σ0 (Baseline, Nu) or
# σ0 + σ1 S (Sigma). Pool volume V = 1, so γ stands for γ/V. With ν1 = σ1 = 0 these are the
# Briggs et al. (2010) equations (src/models/briggs2010.jl).
#
# Stochastic model (ESM Table S3, as coded in WithinHostModel.c). Channels and rates:
#   encounter      γ Z N          Z − 1; a live frog j chosen uniformly; S_j + 1 with prob ν(S_j)
#   release        η ΣS           with prob f: a frog j chosen ∝ S_j, S_j + 1 with prob ν(S_j)
#                                 (the zoospore is lost otherwise); with prob 1 − f: Z + 1
#   sporangium     Σ σ(S_i) S_i   S_j − 1, frog j chosen ∝ σ(S_j) S_j
#   zoospore death μ Z            Z − 1
# A frog dies when an infection would take its load above Smax; its sporangia go with it.
# Here the release channel is split, without changing the process, into self-reinfection at
# rate η f ν0 ΣS (accepted with prob ν(S_j)/ν0) and release to the pool at η(1 − f) ΣS, which
# drops the null events. A run ends when all frogs are dead (frog extinction), when no
# sporangia or zoospores are left (fungal extinction) or at tfinal (persistence).
#
# Interventions (§3.3; FrogAdditions.c): one pulse at a whole day.
#   cleaning frogs   Binomial(N, fraction) live frogs, chosen uniformly, lose all sporangia
#   culling frogs    Binomial(N, fraction) live frogs, chosen uniformly, are removed
#   cleaning the env Binomial(Z, fraction) zoospores are removed. The authors' clean_env()
#                    sets Z = 0 before drawing the binomial, so their code removes every
#                    zoospore whatever the fraction; `as_coded = true` does the same.
# An intervention is not applied once the run's outcome is settled (the C code stops then).
#
# Two printed formulas in the ESM disagree with the rest of the paper:
# - Eq. (11) for ν1 drops the logarithm and the ηf term of the S* expression in ESM §1.1.2.
#   drawert_variant(p, :nu) inverts the §1.1.2 expression; `as_printed = true` uses eq. (11).
# - Eq. (10) divides log-differences by log(c_range); the authors' 2132 points are spaced as
#   if divided by the range of log c (nearest neighbours ≥ 0.0901 > ε = 0.09). The default of
#   nondimensional_distance uses the range of log c; `as_printed = true` the printed form.

const DRAWERT_VARIANTS = (:baseline, :nu, :sigma)
const DRAWERT_OUTCOMES = (:frogs_extinct, :fungus_extinct, :persistence)
const DRAWERT_INTERVENTION_KINDS = (:clean_frogs, :cull_frogs, :clean_environment)

"""
    DrawertParameters(; N, γ, η, ν0, f, σ0, μ, ν1 = 0, σ1 = 0, Smax, S0, Z0 = 0, V = 1)
    DrawertParameters(:drawert2017; set = "", point = nothing, kwargs...)

Parameters of the Drawert et al. (2017) model. Rates are per day.

- `N`: initial number of frogs;
- `γ`: encounter rate between a zoospore and a frog (frog⁻¹ day⁻¹, for a pool of volume `V`);
- `η`: zoospores released per sporangium per day;
- `ν0`: probability that a zoospore encysts when it meets an uninfected frog;
- `f`: fraction of released zoospores that meet their own frog again;
- `σ0`: sporangium loss rate;
- `μ`: zoospore death rate in the pool;
- `ν1`: decline of encystment with load, `ν(S) = ν0 exp(−ν1 S)` (the Nu model);
- `σ1`: increase of sporangium loss with load, `σ(S) = σ0 + σ1 S` (the Sigma model);
- `Smax`: lethal load (sporangia); a frog dies when its load would exceed it;
- `S0`: sporangia on the one frog infected at time 0; `Z0`: zoospores in the pool then;
- `V`: pool volume (1 in the paper).

`ν1 = σ1 = 0` is the Baseline model. See [`drawert_variant`](@ref) for the Nu and Sigma
models. `BriggsParameters(p)` gives the [`BriggsParameters`](@ref) with the same Baseline
dynamics (`ν = ν0`, `σ = σ0`, the same `V`, `γ`, `f`, `μ`, `η` and `Smax`): with
`ν1 = σ1 = 0` the deterministic model is the Briggs et al. (2010) model.

The compendium form reads `Smax`, `S0`, `Z0` and `V` from `data/compendium/drawert2017.toml`:

- with `set = "groupA"` … `"groupE"`, the biological parameters are the representative
  point of that group (ESM Table S9, as printed);
- with `point = id`, they are point `id` of the authors' table of 2132 sampled points
  ([`load_drawert2017`](@ref)), at full precision.

Keyword arguments override single values.
"""
Base.@kwdef struct DrawertParameters
    N::Int
    γ::Float64
    η::Float64
    ν0::Float64
    f::Float64
    σ0::Float64
    μ::Float64
    ν1::Float64 = 0.0
    σ1::Float64 = 0.0
    Smax::Int
    S0::Int
    Z0::Int = 0
    V::Float64 = 1.0
end

function DrawertParameters(key::Symbol; set::AbstractString = "",
        point::Union{Nothing,Integer} = nothing, kwargs...)
    q = parameter_set(key; set)
    base = (; Smax = round(Int, q.Smax), S0 = round(Int, q.S0), Z0 = round(Int, q.Z0), V = q.V)
    if point !== nothing
        isempty(set) || throw(ArgumentError("give either `set` or `point`, not both"))
        r = _drawert_point(point)
        base = merge(base, (; N = r.N, γ = r.γ, η = r.η, ν0 = r.ν0, f = r.f, σ0 = r.σ0, μ = r.μ))
    elseif !isempty(set)
        base = merge(base, (; N = round(Int, q.N), γ = q.γ, η = q.η, ν0 = q.ν0, f = q.f,
                            σ0 = q.σ0, μ = q.μ))
    end
    return DrawertParameters(; merge(base, values(kwargs))...)
end

time_unit(::Type{<:DrawertParameters}) = :day

function _check(p::DrawertParameters)
    (p.N >= 1 && p.S0 >= 0 && p.Z0 >= 0 && p.Smax >= 1 && p.S0 <= p.Smax) ||
        throw(ArgumentError("need N ≥ 1, Smax ≥ 1 and 0 ≤ S0 ≤ Smax, Z0 ≥ 0"))
    (p.γ >= 0 && p.η >= 0 && p.σ0 >= 0 && p.μ >= 0 && p.ν1 >= 0 && p.σ1 >= 0 && p.V > 0) ||
        throw(ArgumentError("rates must be non-negative and V positive"))
    (0 <= p.ν0 <= 1 && 0 <= p.f <= 1) || throw(ArgumentError("ν0 and f must lie in [0, 1]"))
    return p
end

"""
    model_variant(p::DrawertParameters) -> Symbol

`:baseline` (`ν1 = σ1 = 0`), `:nu` (`ν1 > 0`) or `:sigma` (`σ1 > 0`). The paper never uses
both; a parameter set with both is reported as `:both`.
"""
function model_variant(p::DrawertParameters)
    p.ν1 > 0 && p.σ1 > 0 && return :both
    p.ν1 > 0 && return :nu
    p.σ1 > 0 && return :sigma
    return :baseline
end

"""
    encystment_probability(p::DrawertParameters, S) -> ν(S)

Probability that a zoospore meeting a frog with load `S` encysts: `ν0 exp(−ν1 S)`.
"""
encystment_probability(p::DrawertParameters, S) = p.ν0 * exp(-p.ν1 * S)

"""
    sporangium_loss_rate(p::DrawertParameters, S) -> σ(S)

Per-sporangium loss rate (day⁻¹) on a frog with load `S`: `σ0 + σ1 S`.
"""
sporangium_loss_rate(p::DrawertParameters, S) = p.σ0 + p.σ1 * S

# The Briggs et al. (2010) parameters with the same Baseline dynamics: ν = ν0, σ = σ0 and the
# same V, γ, f, μ, η and Smax (documented under DrawertParameters).
BriggsParameters(p::DrawertParameters) =
    BriggsParameters(; V = p.V, γ = p.γ, ν = p.ν0, f = p.f, μ = p.μ, η = p.η, σ = p.σ0,
                     Smax = Float64(p.Smax))

"""
    load_growth_rate(p::DrawertParameters) -> λ

The linearised Bd growth rate of the Baseline model with `p.N` frogs: the dominant
eigenvalue of the linear system (ESM §1.1.1, §1.4), the same as
`load_growth_rate(BriggsParameters(p), p.N)`. The paper kept points with
`0.035 ≤ λ ≤ 0.35` day⁻¹.
"""
load_growth_rate(p::DrawertParameters) = load_growth_rate(BriggsParameters(p), p.N)

"""
    nondimensional_parameters(p::DrawertParameters) -> (; c1, c2, c3, c4)

The non-dimensional groups of Table 2 (ESM eqs 3–6), with `γ` read as `γ/V`:

- `c1 = η ν0 f / σ0`: self-reinfection;
- `c2 = γ η (1 − f) ν0 N / σ0²`: infectiousness of the frogs;
- `c3 = γ N / σ0`: zoospore loss by meeting frogs;
- `c4 = μ / σ0`: background zoospore loss.
"""
function nondimensional_parameters(p::DrawertParameters)
    γ = p.γ / p.V
    return (; c1 = p.η * p.ν0 * p.f / p.σ0,
            c2 = γ * p.η * (1 - p.f) * p.ν0 * p.N / p.σ0^2,
            c3 = γ * p.N / p.σ0,
            c4 = p.μ / p.σ0)
end

"""
    nondimensional_ranges(cs) -> (; c1, c2, c3, c4)

The ranges used to scale [`nondimensional_distance`](@ref), from a collection `cs` of
`(c1, c2, c3, c4)` named tuples: `max − min` of `c1`, `log c2`, `log c3` and `c4`. With
`as_printed = true` the `c2` and `c3` entries are `max − min` of `c2` and `c3` themselves.
"""
function nondimensional_ranges(cs; as_printed::Bool = false)
    rng(f) = (v = [f(c) for c in cs]; maximum(v) - minimum(v))
    return as_printed ?
        (; c1 = rng(c -> c.c1), c2 = rng(c -> c.c2), c3 = rng(c -> c.c3), c4 = rng(c -> c.c4)) :
        (; c1 = rng(c -> c.c1), c2 = rng(c -> log(c.c2)), c3 = rng(c -> log(c.c3)),
         c4 = rng(c -> c.c4))
end

"""
    nondimensional_distance(a, b, ranges; as_printed = false)

The distance between two points of non-dimensional parameter space (ESM eq. 10), each a
`(c1, c2, c3, c4)` named tuple:

`d = √(Δc1²/r1² + Δlog c2²/r2² + Δlog c3²/r3² + Δc4²/r4²)`

`ranges` comes from [`nondimensional_ranges`](@ref). By default `r2` and `r3` are the ranges
of `log c2` and `log c3`, which is how the authors' points are spaced: no two of the 2132 are
closer than 0.0901, just above the threshold ε = 0.09. Eq. (10) prints `log(c2range)`, the
log of the range of `c2`; `as_printed = true` uses that (pass ranges made with
`as_printed = true`).
"""
function nondimensional_distance(a, b, r; as_printed::Bool = false)
    d2 = r.c2
    d3 = r.c3
    if as_printed
        d2 = abs(log(r.c2))
        d3 = abs(log(r.c3))
    end
    return sqrt(((a.c1 - b.c1) / r.c1)^2 + ((log(a.c2) - log(b.c2)) / d2)^2 +
                ((log(a.c3) - log(b.c3)) / d3)^2 + ((a.c4 - b.c4) / r.c4)^2)
end

# The quasi-steady per-frog gain K: with Z at γ/V Z = γη(1−f)N S/(γN + μV) for N frogs that
# all carry S, dS/dt = S [ν(S) K − σ(S)] with K = γη(1−f)N/(γN + μV) + ηf.
_drawert_K(p::DrawertParameters) = p.γ * p.η * (1 - p.f) * p.N / (p.γ * p.N + p.μ * p.V) + p.η * p.f

"""
    equilibrium_load_ratio(p::DrawertParameters) -> S*

The load at which Bd growth on each frog stops, as a multiple of `Smax`, for `N` frogs with
equal loads and the pool at quasi-steady state (ESM §1.1.2, §1.1.3). With
`K = γη(1 − f)N/(γN + μV) + ηf`:

- Nu model: `S* = ln(ν0 K/σ0)/(ν1 Smax)`;
- Sigma model: `S* = (ν0 K − σ0)/(σ1 Smax)`;
- Baseline: `Inf`.

Frogs die in the deterministic model when `S* > 1`.
"""
function equilibrium_load_ratio(p::DrawertParameters)
    v = model_variant(p)
    v === :both && throw(ArgumentError("S* is defined for one load-dependence at a time"))
    v === :baseline && return Inf
    K = _drawert_K(p)
    v === :nu && return log(p.ν0 * K / p.σ0) / (p.ν1 * p.Smax)
    return (p.ν0 * K - p.σ0) / (p.σ1 * p.Smax)
end

"""
    drawert_variant(p, variant; S_star, as_printed = false) -> DrawertParameters

The Baseline point `p` turned into the `:nu` or `:sigma` model (or back to `:baseline`), with
the load dependence chosen so that loads level off at `S_star × Smax` (ESM §1.4; the paper
uses `S_star = 1.2` and `2`, which the compendium holds as `S_star_low`, `S_star_high`):

- `:nu`: `ν1 = ln(ν0 K/σ0)/(S* Smax)`, the inverse of the S* formula of ESM §1.1.2;
- `:sigma`: `σ1 = (ν0 K − σ0)/(S* Smax)`, ESM eq. (12);

with `K` as in [`equilibrium_load_ratio`](@ref). ESM eq. (11) prints
`ν1 = ν0 γη(1 − f)N/(σ0 (γN + μ))/(S* Smax)`, without the logarithm and the `ηf` term;
`as_printed = true` uses it. The printed ν1 does not give loads that level off at
`S* × Smax`.
"""
function drawert_variant(p::DrawertParameters, variant::Symbol; S_star::Real = 1.2,
        as_printed::Bool = false)
    variant in DRAWERT_VARIANTS || throw(ArgumentError("variant must be one of $DRAWERT_VARIANTS"))
    b = DrawertParameters(; N = p.N, γ = p.γ, η = p.η, ν0 = p.ν0, f = p.f, σ0 = p.σ0, μ = p.μ,
        Smax = p.Smax, S0 = p.S0, Z0 = p.Z0, V = p.V)
    variant === :baseline && return b
    S_star > 0 || throw(ArgumentError("S_star must be positive"))
    K = _drawert_K(b)
    p.ν0 * K > p.σ0 ||
        throw(ArgumentError("loads do not grow (ν0 K ≤ σ0), so no S* can be reached"))
    scale = S_star * p.Smax
    if variant === :nu
        ν1 = as_printed ? p.ν0 * p.γ * p.η * (1 - p.f) * p.N / (p.σ0 * (p.γ * p.N + p.μ * p.V)) / scale :
             log(p.ν0 * K / p.σ0) / scale
        return DrawertParameters(; N = b.N, γ = b.γ, η = b.η, ν0 = b.ν0, f = b.f, σ0 = b.σ0,
            μ = b.μ, ν1, Smax = b.Smax, S0 = b.S0, Z0 = b.Z0, V = b.V)
    end
    return DrawertParameters(; N = b.N, γ = b.γ, η = b.η, ν0 = b.ν0, f = b.f, σ0 = b.σ0,
        μ = b.μ, σ1 = (p.ν0 * K - p.σ0) / scale, Smax = b.Smax, S0 = b.S0, Z0 = b.Z0, V = b.V)
end

# --------------------------------------------------------------------------------------------
# Interventions

"""
    DrawertIntervention(kind, day, fraction)

A one-time intervention (§3.3), applied at time `day`:

- `:clean_frogs`: `Binomial(N, fraction)` live frogs, chosen at random, are cleared of all
  sporangia (antifungal treatment);
- `:cull_frogs`: `Binomial(N, fraction)` live frogs, chosen at random, are removed;
- `:clean_environment`: `Binomial(Z, fraction)` zoospores are removed from the pool (all of
  them with `as_coded = true` in [`DrawertIBMProblem`](@ref)).

`N` and `Z` are the live frogs and pool zoospores at `day`.
"""
struct DrawertIntervention
    kind::Symbol
    day::Float64
    fraction::Float64
    function DrawertIntervention(kind::Symbol, day::Real, fraction::Real)
        kind in DRAWERT_INTERVENTION_KINDS ||
            throw(ArgumentError("kind must be one of $DRAWERT_INTERVENTION_KINDS"))
        day > 0 || throw(ArgumentError("the intervention day must be positive"))
        0 <= fraction <= 1 || throw(ArgumentError("fraction must lie in [0, 1]"))
        return new(kind, Float64(day), Float64(fraction))
    end
end

"""
    drawert_interventions(key = :drawert2017) -> Vector{DrawertIntervention}

The 27 one-time strategies of the paper (§3.3): each of `:clean_frogs`, `:cull_frogs` and
`:clean_environment`, on day 14, 28 or 56, treating 25%, 50% or 75%. Days and fractions are
read from the compendium. Ordered by kind, then day, then fraction.
"""
function drawert_interventions(key::Symbol = :drawert2017)
    q = parameter_set(key)
    days = (q.day_1, q.day_2, q.day_3)
    fracs = (q.fraction_1, q.fraction_2, q.fraction_3)
    return [DrawertIntervention(k, d, x) for k in DRAWERT_INTERVENTION_KINDS for d in days
            for x in fracs]
end

# --------------------------------------------------------------------------------------------
# The authors' table of sampled points (bundled, CC BY 4.0)

"""
    drawert2017_path()

Path of the bundled table of the 2132 Baseline parameter points, with the authors' own
30-run summaries for each (figshare 10.6084/m9.figshare.5328127, CC BY 4.0). See
`data/drawert2017/README.md`.
"""
drawert2017_path() = joinpath(pkgdir(@__MODULE__), "data", "drawert2017", "baseline_points.csv")

const _DRAWERT_COLUMNS = (
    id = "ID", N = "param_Nfrogs", γ = "param_gamma", η = "param_eta", ν0 = "param_nu",
    f = "param_f", σ0 = "param_sigma", μ = "param_mu",
    mean_days_to_first_death = "data_meanDaysToSmax", std_days_to_first_death = "data_stdDaysToSmax",
    mean_surviving = "data_meanFracSurviving", std_surviving = "data_stdFracSurviving",
    mean_days_to_Z0 = "data_meanDaysToZ0", std_days_to_Z0 = "data_stdDaysToZ0",
    mean_days_to_F0 = "data_meanDaysToF0", std_days_to_F0 = "data_stdDaysToF0",
    frac_persistence = "data_fracPersistence", frac_frogs_extinct = "data_fracFrogsExtinct",
    frac_fungus_extinct = "data_fracFungusExtinct", λ = "data_lambda",
    simulation_time = "simulation_time")

const _DRAWERT_TABLE = Ref{Any}(nothing)

"""
    load_drawert2017() -> Vector{NamedTuple}

The authors' 2132 Baseline parameter points ([`drawert2017_path`](@ref)), one named tuple per
point, ordered by `id`:

- parameters `id`, `N`, `γ`, `η`, `ν0`, `f`, `σ0`, `μ`;
- the authors' 30-run results: `mean_surviving`, `std_surviving` (population sd),
  `frac_frogs_extinct`, `frac_fungus_extinct`, `frac_persistence`,
  `mean_days_to_first_death` / `std_…` (the C code's `daysToSmax`: the first daily output
  with fewer than `N` frogs), `mean_days_to_Z0`, `mean_days_to_F0` and their sds;
- `λ`, the authors' linearised growth rate, and `simulation_time`, their CPU seconds for the
  30 runs.
"""
function load_drawert2017()
    _DRAWERT_TABLE[] === nothing || return _DRAWERT_TABLE[]::Vector
    lines = filter(!isempty ∘ strip, readlines(drawert2017_path()))
    header = split(lines[1], ',')
    idx = map(c -> (i = findfirst(==(c), header);
                    i === nothing ? throw(ArgumentError("column $c missing")) : i), _DRAWERT_COLUMNS)
    rows = map(lines[2:end]) do l
        v = split(l, ',')
        get(k) = parse(Float64, v[idx[k]])
        NamedTuple{keys(_DRAWERT_COLUMNS)}(Tuple(k in (:id, :N) ? round(Int, get(k)) : get(k)
                                                 for k in keys(_DRAWERT_COLUMNS)))
    end
    sort!(rows; by = r -> r.id)
    _DRAWERT_TABLE[] = rows
    return rows
end

function _drawert_point(id::Integer)
    rows = load_drawert2017()
    i = searchsortedfirst(rows, id; by = r -> r isa Integer ? r : r.id)
    (i <= length(rows) && rows[i].id == id) || throw(ArgumentError("no point with id $id"))
    return rows[i]
end

# Design constants read once from the compendium.
const _DRAWERT_DESIGN = Ref{Any}(nothing)
function _drawert_design()
    if _DRAWERT_DESIGN[] === nothing
        q = parameter_set(:drawert2017)
        _DRAWERT_DESIGN[] = (; tfinal = Float64(q.tfinal), n_runs = round(Int, q.n_runs),
            group_sd = q.group_sd, group_mean_low = q.group_mean_low,
            group_mean_high = q.group_mean_high, effects = (q.effect_1, q.effect_2, q.effect_3))
    end
    return _DRAWERT_DESIGN[]::NamedTuple
end

# --------------------------------------------------------------------------------------------
# Stochastic individual-based model

"""
    DrawertIBMProblem(p; intervention = nothing, tfinal = 1000, saveat = Float64[],
                      save_loads = false, as_coded = false)

The stochastic individual-based model of Drawert et al. (2017) (ESM §1.3, Table S3): an
integer load on each frog and an integer number of zoospores in the pool, simulated exactly
with the direct method (Gillespie). One frog starts with `p.S0` sporangia, the others with
none, and the pool with `p.Z0` zoospores. `tfinal` (1000 days) comes from the compendium.

`intervention` is a [`DrawertIntervention`](@ref) or `nothing`. With `as_coded = true`,
cleaning the environment removes every zoospore, as the authors' C code does (the paper
describes removing the given fraction).

`solve(prob; rng)` returns a named tuple:

- `outcome`: `:frogs_extinct`, `:fungus_extinct` (no sporangia or zoospores left, some frogs
  alive) or `:persistence` (both present at `tfinal`);
- `survivors`, `surviving_fraction` (survivors / `N`; culled frogs do not survive);
- `death_time`: per frog, the time it died of chytridiomycosis (`NaN` otherwise);
  `culled`: per frog, whether it was culled;
- `first_loss`: the first time a frog died or was culled (`Inf` if none);
- `t_end`: when the outcome was settled (or `tfinal`); `events`: the number of reactions;
- `t`, `frogs`, `zoospores`, `sporangia`: the state at the `saveat` times, and, with
  `save_loads = true`, `loads`, an `N × length(saveat)` matrix (`-1` for a dead frog).

The run stops when its outcome is settled unless later `saveat` times remain.
"""
struct DrawertIBMProblem
    p::DrawertParameters
    intervention::Union{Nothing,DrawertIntervention}
    tfinal::Float64
    saveat::Vector{Float64}
    save_loads::Bool
    as_coded::Bool
end

function DrawertIBMProblem(p::DrawertParameters;
        intervention::Union{Nothing,DrawertIntervention} = nothing,
        tfinal::Real = _drawert_design().tfinal, saveat = Float64[], save_loads::Bool = false,
        as_coded::Bool = false)
    _check(p)
    model_variant(p) === :both && throw(ArgumentError("set ν1 or σ1, not both"))
    tfinal > 0 || throw(ArgumentError("tfinal must be positive"))
    ts = sort!(Float64.(collect(saveat)))
    all(t -> 0 <= t <= tfinal, ts) || throw(ArgumentError("saveat must lie in [0, tfinal]"))
    return DrawertIBMProblem(p, intervention, Float64(tfinal), ts, save_loads, as_coded)
end

# Fenwick (binary indexed) tree over non-negative integer weights.
@inline function _fw_add!(tree::Vector{Int}, i::Int, δ::Int)
    n = length(tree)
    @inbounds while i <= n
        tree[i] += δ
        i += i & (-i)
    end
    return nothing
end

# Smallest index i with prefix sum ≥ r, for 1 ≤ r ≤ total.
@inline function _fw_find(tree::Vector{Int}, r::Int)
    n = length(tree)
    pos = 0
    step = n == 0 ? 0 : prevpow(2, n)
    @inbounds while step > 0
        nxt = pos + step
        if nxt <= n && tree[nxt] < r
            pos = nxt
            r -= tree[nxt]
        end
        step >>= 1
    end
    return pos + 1
end

mutable struct _DrawertState
    S::Vector{Int}          # load per frog; -1 once dead or culled
    alive::Vector{Int}      # indices of live frogs
    slot::Vector{Int}       # position of each frog in `alive` (0 if not alive)
    fw1::Vector{Int}        # Fenwick tree of S
    fw2::Vector{Int}        # Fenwick tree of S² (Sigma model)
    ΣS::Int
    ΣS2::Int                # kept only for the Sigma model (sq = true)
    Z::Int
    sq::Bool
end

function _drawert_state(p::DrawertParameters)
    N = p.N
    st = _DrawertState(zeros(Int, N), collect(1:N), collect(1:N), zeros(Int, N), zeros(Int, N),
        0, 0, p.Z0, p.σ1 > 0)
    _set_load!(st, 1, p.S0)
    return st
end

@inline function _set_load!(st::_DrawertState, j::Int, s::Int)
    @inbounds old = st.S[j]
    @inbounds st.S[j] = s
    st.ΣS += s - old
    _fw_add!(st.fw1, j, s - old)
    if st.sq
        st.ΣS2 += s * s - old * old
        _fw_add!(st.fw2, j, s * s - old * old)
    end
    return nothing
end

# Remove live frog j (death or cull) at time t.
@inline function _remove_frog!(st::_DrawertState, j::Int)
    _set_load!(st, j, 0)
    @inbounds begin
        st.S[j] = -1
        k = st.slot[j]
        last = st.alive[end]
        st.alive[k] = last
        st.slot[last] = k
        pop!(st.alive)
        st.slot[j] = 0
    end
    return nothing
end

# One sporangium more on frog j; the frog dies if that would take it above Smax.
@inline function _infect!(st::_DrawertState, j::Int, Smax::Int, t::Float64, death_time)
    @inbounds s = st.S[j]
    if s + 1 > Smax
        _remove_frog!(st, j)
        @inbounds death_time[j] = t
    else
        _set_load!(st, j, s + 1)
    end
    return nothing
end

function _apply_intervention!(st::_DrawertState, iv::DrawertIntervention, rng::AbstractRNG,
        as_coded::Bool, culled, t)
    n = length(st.alive)
    if iv.kind === :clean_environment
        st.Z = as_coded ? 0 : st.Z - rand(rng, Binomial(st.Z, iv.fraction))
        return nothing
    end
    n == 0 && return nothing
    k = rand(rng, Binomial(n, iv.fraction))
    chosen = copy(st.alive)                   # k distinct live frogs, uniformly
    for i in 1:k
        r = rand(rng, i:n)
        chosen[i], chosen[r] = chosen[r], chosen[i]
    end
    for i in 1:k
        j = chosen[i]
        if iv.kind === :clean_frogs
            _set_load!(st, j, 0)
        else
            _remove_frog!(st, j)
            culled[j] = true
        end
    end
    return nothing
end

function CommonSolve.solve(prob::DrawertIBMProblem; rng::AbstractRNG = Random.default_rng())
    p = prob.p
    N, Smax = p.N, p.Smax
    st = _drawert_state(p)
    γ = p.γ / p.V
    a_self = p.η * p.f * p.ν0           # self-reinfection attempts per sporangium
    a_pool = p.η * (1 - p.f)            # releases to the pool per sporangium
    nu_dep = p.ν1 > 0
    sig_dep = p.σ1 > 0
    death_time = fill(NaN, N)
    culled = falses(N)
    ns = length(prob.saveat)
    frogs = zeros(Int, ns)
    zoos = zeros(Int, ns)
    spor = zeros(Int, ns)
    loads = prob.save_loads ? fill(-1, N, ns) : zeros(Int, 0, 0)
    k = 1
    function record!(k)
        frogs[k] = length(st.alive)
        zoos[k] = st.Z
        spor[k] = st.ΣS
        if prob.save_loads
            @inbounds for j in 1:N
                loads[j, k] = st.S[j]
            end
        end
    end
    t = 0.0
    tf = prob.tfinal
    ti = prob.intervention === nothing ? Inf : prob.intervention.day
    outcome = :none
    t_end = tf
    events = 0
    while true
        n = length(st.alive)
        if outcome === :none
            if n == 0
                outcome, t_end, ti = :frogs_extinct, t, Inf
            elseif st.ΣS == 0 && st.Z == 0
                outcome, t_end, ti = :fungus_extinct, t, Inf
            end
        end
        outcome !== :none && k > ns && break
        a1 = γ * st.Z * n
        a2 = a_self * st.ΣS
        a3 = a_pool * st.ΣS
        a4 = sig_dep ? p.σ0 * st.ΣS + p.σ1 * st.ΣS2 : p.σ0 * st.ΣS
        a5 = p.μ * st.Z
        a0 = a1 + a2 + a3 + a4 + a5
        tnext = a0 > 0 ? t + Random.randexp(rng) / a0 : Inf
        tstop = min(ti, tf)
        if tnext >= tstop
            # nothing happens before tstop (the drawn reaction is discarded: memoryless)
            while k <= ns && prob.saveat[k] < tstop
                record!(k); k += 1
            end
            t = tstop
            if tstop == ti
                _apply_intervention!(st, prob.intervention, rng, prob.as_coded, culled, t)
                ti = Inf
                continue
            end
            while k <= ns
                record!(k); k += 1
            end
            outcome === :none && (outcome = :persistence; t_end = tf)
            break
        end
        while k <= ns && prob.saveat[k] < tnext
            record!(k); k += 1
        end
        t = tnext
        events += 1
        r = rand(rng) * a0
        if r < a1                                   # encounter: zoospore meets a live frog
            st.Z -= 1
            j = @inbounds st.alive[rand(rng, 1:n)]
            if rand(rng) < (nu_dep ? encystment_probability(p, @inbounds st.S[j]) : p.ν0)
                _infect!(st, j, Smax, t, death_time)
            end
        elseif r < a1 + a2                          # zoospore re-encounters its own frog
            j = _fw_find(st.fw1, rand(rng, 1:st.ΣS))
            if !nu_dep || rand(rng) < exp(-p.ν1 * @inbounds st.S[j])
                _infect!(st, j, Smax, t, death_time)
            end
        elseif r < a1 + a2 + a3                     # zoospore released to the pool
            st.Z += 1
        elseif r < a1 + a2 + a3 + a4                # sporangium lost
            j = if sig_dep && rand(rng) * a4 >= p.σ0 * st.ΣS
                _fw_find(st.fw2, rand(rng, 1:st.ΣS2))
            else
                _fw_find(st.fw1, rand(rng, 1:st.ΣS))
            end
            _set_load!(st, j, @inbounds(st.S[j]) - 1)
        else                                        # zoospore dies
            st.Z -= 1
        end
    end
    survivors = length(st.alive)
    first_loss = minimum((isnan(death_time[j]) ? (culled[j] ? prob.intervention.day : Inf) : death_time[j])
                         for j in 1:N)
    res = (; outcome, survivors, surviving_fraction = survivors / N, death_time,
           culled = collect(culled), first_loss, t_end, events,
           t = copy(prob.saveat), frogs, zoospores = zoos, sporangia = spor)
    return prob.save_loads ? merge(res, (; loads)) : res
end

"""
    drawert_summary(prob::DrawertIBMProblem; nruns = 30, rng) -> NamedTuple

Run `nruns` realisations (30 by default, from the compendium) and summarise them as the
authors' code does (ESM Table S4; `print_stats` in WithinHostModel.c):

- `mean_surviving`, `std_surviving`: mean and population sd of the surviving fraction;
- `frac_frogs_extinct`, `frac_fungus_extinct`, `frac_persistence`;
- `mean_days_to_first_death`, `std_days_to_first_death`: the first whole day on which fewer
  than `N` frogs are alive (deaths and culls; `tfinal` if none), as the C code's
  `daysToSmax`;
- `outcome_class` ([`outcome_class`](@ref)) and `group` ([`survival_group`](@ref));
- `surviving`, the surviving fraction of each run.
"""
function drawert_summary(prob::DrawertIBMProblem; nruns::Integer = _drawert_design().n_runs,
        rng::AbstractRNG = Random.default_rng())
    nruns >= 1 || throw(ArgumentError("nruns must be at least 1"))
    surv = zeros(nruns)
    first = zeros(nruns)
    counts = Dict(o => 0 for o in DRAWERT_OUTCOMES)
    for i in 1:nruns
        r = solve(prob; rng)
        surv[i] = r.surviving_fraction
        first[i] = isfinite(r.first_loss) ? ceil(r.first_loss) : prob.tfinal
        counts[r.outcome] += 1
    end
    m = sum(surv) / nruns
    sd = sqrt(max(sum(abs2, surv) / nruns - m^2, 0.0))
    mf = sum(first) / nruns
    sf = sqrt(max(sum(abs2, first) / nruns - mf^2, 0.0))
    fe = counts[:frogs_extinct] / nruns
    fu = counts[:fungus_extinct] / nruns
    pe = counts[:persistence] / nruns
    return (; mean_surviving = m, std_surviving = sd, frac_frogs_extinct = fe,
            frac_fungus_extinct = fu, frac_persistence = pe, mean_days_to_first_death = mf,
            std_days_to_first_death = sf, outcome_class = outcome_class(fe, fu, pe),
            group = survival_group(m, sd), surviving = surv, nruns = Int(nruns))
end

"""
    outcome_class(frac_frogs_extinct, frac_fungus_extinct, frac_persistence) -> Symbol

The outcome category of a parameter point in Fig. 2a: `:frog_extinction`,
`:fungal_extinction` or `:coexistence` when every run ends that way, and `:variable`
otherwise.
"""
function outcome_class(fe::Real, fu::Real, pe::Real)
    fe == 1 && return :frog_extinction
    fu == 1 && return :fungal_extinction
    pe == 1 && return :coexistence
    return :variable
end

"""
    survival_group(mean, sd) -> Symbol

The group of a parameter point (Fig. 3a) from the mean and sd of the fraction of frogs that
survive its runs. With the thresholds of the compendium (0.1 and 0.5):

- `:A`: mean and sd both zero (always extinction);
- `:C`: sd above 0.1;
- `:B`: otherwise, mean below 0.1;
- `:D`: otherwise, mean above 0.5;
- `:E`: otherwise (mean in [0.1, 0.5]; partial die-offs).
"""
function survival_group(m::Real, sd::Real)
    d = _drawert_design()
    m == 0 && sd == 0 && return :A
    sd > d.group_sd && return :C
    m < d.group_mean_low && return :B
    m > d.group_mean_high && return :D
    return :E
end

"""
    intervention_effect(Δ) -> Float64

The shading class of Fig. 4 for a change `Δ` in the mean surviving fraction (with minus
without the intervention): `±0.5`, `±0.3` or `±0.1` for the largest threshold that `|Δ|`
reaches, else `0`. The thresholds come from the compendium.
"""
function intervention_effect(Δ::Real)
    for e in reverse(_drawert_design().effects)
        abs(Δ) >= e - 1e-12 && return sign(Δ) * e
    end
    return 0.0
end

# --------------------------------------------------------------------------------------------
# Deterministic model

"""
    DrawertODEProblem(p; intervention = nothing, tfinal = 1000)

The deterministic version of the model (eqs 3.1–3.2; ESM eqs 1–2): one load per frog and
the pool, from the same initial state as [`DrawertIBMProblem`](@ref). A frog dies (its load
is removed) when it reaches `Smax`.

Interventions follow the authors' `determ_frogs.c`, which the paper does not describe:

- cleaning and culling act on `floor(fraction × N)` live frogs, at least one, taken in index
  order (frog 1 is the one infected first);
- culling does nothing when one frog or none is left;
- cleaning the environment multiplies `Z` by `1 − fraction`.

`solve(prob; kwargs...)` integrates with `Tsit5` and returns
`(; sol, outcome, survivors, surviving_fraction, death_time, t_end)`. The outcome is as for
the stochastic model: `:fungus_extinct` when, on some whole day, the total load and `Z` are
both below `1e-6`, as in the authors' code. `to_ode_problem(prob)` gives the `ODEProblem`.
"""
struct DrawertODEProblem
    p::DrawertParameters
    intervention::Union{Nothing,DrawertIntervention}
    tfinal::Float64
end

function DrawertODEProblem(p::DrawertParameters;
        intervention::Union{Nothing,DrawertIntervention} = nothing,
        tfinal::Real = _drawert_design().tfinal)
    _check(p)
    model_variant(p) === :both && throw(ArgumentError("set ν1 or σ1, not both"))
    return DrawertODEProblem(p, intervention, Float64(tfinal))
end

struct _DrawertODEParams
    p::DrawertParameters
    alive::Vector{Bool}
    death_time::Vector{Float64}
end

function _drawert_rhs!(du, u, q::_DrawertODEParams, t)
    p, alive = q.p, q.alive
    N = length(u) - 1
    Z = u[end]
    γZ = p.γ / p.V * Z
    pool = zero(eltype(u))
    @inbounds for i in 1:N
        if alive[i]
            S = u[i]
            ν = encystment_probability(p, S)
            du[i] = ν * (γZ + p.η * p.f * S) - sporangium_loss_rate(p, S) * S
            pool += p.η * (1 - p.f) * S - γZ
        else
            du[i] = zero(eltype(u))
        end
    end
    du[end] = pool - p.μ * Z
    return nothing
end

function to_ode_problem(prob::DrawertODEProblem)
    p = prob.p
    u0 = zeros(p.N + 1)
    u0[1] = p.S0
    u0[end] = p.Z0
    return ODEProblem(_drawert_rhs!, u0, (0.0, prob.tfinal),
        _DrawertODEParams(p, fill(true, p.N), fill(NaN, p.N)))
end

function _drawert_death_callback(N)
    condition(out, u, t, integ) = (q = integ.p;
        @inbounds for i in 1:N
            out[i] = q.alive[i] ? u[i] - q.p.Smax : -one(eltype(u))
        end)
    function affect!(integ, events)
        q = integ.p
        for i in eachindex(events)
            iszero(events[i]) && continue
            (q.alive[i] && integ.u[i] >= q.p.Smax / 2) || continue
            q.alive[i] = false
            q.death_time[i] = integ.t
            integ.u[i] = 0.0
        end
        any(q.alive) || SciMLBase.terminate!(integ)
        return nothing
    end
    return VectorContinuousCallback(condition, affect!, N)
end

function _drawert_ode_intervention(iv::DrawertIntervention)
    condition(u, t, integ) = t == iv.day
    function affect!(integ)
        q = integ.p
        live = findall(q.alive)
        n = length(live)
        if iv.kind === :clean_environment
            integ.u[end] *= 1 - iv.fraction
        elseif n > 0 && !(iv.kind === :cull_frogs && n <= 1)
            m = max(1, floor(Int, iv.fraction * n))
            for i in live[1:min(m, n)]
                integ.u[i] = 0.0
                iv.kind === :cull_frogs && (q.alive[i] = false)
            end
            any(q.alive) || SciMLBase.terminate!(integ)
        end
        return nothing
    end
    return DiscreteCallback(condition, affect!)
end

function CommonSolve.solve(prob::DrawertODEProblem; abstol = 1e-8, reltol = 1e-8, kwargs...)
    ode = to_ode_problem(prob)
    q = ode.p
    N = prob.p.N
    cbs = Any[_drawert_death_callback(N)]
    tstops = Float64[]
    if prob.intervention !== nothing
        push!(cbs, _drawert_ode_intervention(prob.intervention))
        push!(tstops, prob.intervention.day)
    end
    sol = solve(ode, Tsit5(); callback = SciMLBase.CallbackSet(cbs...), tstops,
        saveat = 0:1.0:prob.tfinal, abstol, reltol, kwargs...)
    survivors = count(q.alive)
    outcome = :persistence
    t_end = prob.tfinal
    if survivors == 0
        outcome = :frogs_extinct
        t_end = sol.t[end]
    else
        for (t, u) in zip(sol.t, sol.u)
            t >= 1 || continue
            if sum(@view u[1:N]) <= 1e-6 && u[end] <= 1e-6
                outcome, t_end = :fungus_extinct, t
                break
            end
        end
    end
    return (; sol, outcome, survivors, surviving_fraction = survivors / N,
            death_time = copy(q.death_time), t_end)
end
