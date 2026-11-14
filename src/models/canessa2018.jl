# Canessa S, Bozzuto C, Campbell Grant EH, et al. (2018) Decision-making for mitigating wildlife
# diseases: from theory to practice for an emerging fungal pathogen of amphibians. J Appl Ecol
# 55:1987-1996. doi:10.1111/1365-2664.13089
# Wiki: lab/wiki/models/bsal-canessa2018-ipm-decision.md
#
# Fire salamanders and Bsal over three months (30 steps of 3 days), with no recruitment, no
# density dependence and no recovery. Checked against the review copy of the paper (main text,
# Table 1, Figs 1-4, Online Supporting Information eqs S1-S12 and Table S1) and the archived
# parameter draws (Dryad doi:10.5061/dryad.0jh5p, bundled as data/canessa2018/).
#
# Integral projection model (eqs S2-S10), z = ln(load in GE):
#   S(t+1)     = S(t) s0 (1 - φ(t))
#   I(z',t+1)  = S(t) s0 φ(t) G0(z') + ∫ I(z,t) s(z) G(z',z) dz
#   φ(t)       = 1 - exp(-∫ β(z) I(z,t) dz)
#   s(z) = logistic(b0 + b1 z),  β(z) = β0 logistic(b2 + b3 z),  G0 = N(μG0, σ2G0),
#   G(z',z) = N(b4 + b5 z + b6 z², σ2G exp(b7 z)), capped at z ≥ ln Lmax.
# R0 = S0 s0 ∫ β(z) [(I - P)⁻¹ G0](z) dz (Wilber et al. 2016; Klepac & Caswell 2011).
# Simplified SI model (eqs S11-S12): S(t+1) = S s0 (1 - φ), I(t+1) = S s0 φ + I s̄,
#   φ = 1 - exp(-β I), R0 = β(1 - mβ) S0(1 - mS0) s0(1 - ms)/(1 - s̄(1 - ms)).
#
# Points the PDF does not make obvious (details in the compendium notes and the targets):
# - Table S1 gives the fits to the experts' best estimates, but every ensemble run refits the
#   functions to one row of the archived draws, whose distributions lie far from Table S1
#   (the Table S1 point model has R0 ≈ 5.1; the draws give a median ≈ 9.4).
# - The growth regression of eq. S9 is only described; 10 Beverton-Holt steps from L0 = 10 with
#   Lmax = 10⁴ and 25 zoospores per zoospore in 5 days reproduce Table S1's b4-b6 exactly, so
#   each draw is refitted that way.
# - The cap at Lmax: the SI formula (mean ln Lmax above the cap, normal spread) lets loads
#   exceed Lmax and gives R0 = 8.2 (95% interval 5.0-12.2), uniformly ~14% below the paper's
#   9.6 (5.8-14.3). Not letting loads exceed Lmax (cap = :truncate, the default) gives 9.4
#   (5.6-14.3) and reproduces the other R0 values of Fig. 2.
# - Units of β0 (inferred): eq. S6 gives the printed β0 = 0.007 only with S0 as the count on
#   the 0.475-ha Stegen transect (238.9 adults), not 503 ha⁻¹; the archived (β0, S0) pairs
#   agree. Simulating in the same units (units = :transect: S0 × 0.475 susceptibles, one
#   infected adult introduced) reproduces Fig. 2's final densities (e 0.95, f 0.081 against
#   1.0 and 0.083); per hectare (units = :density, I0 = 0.01 ha⁻¹) they are 0.69 and 0.006.
#   R0 keeps S0 per ha, as the paper's R0 values do.

const CANESSA_LOADS = (0.0, 10.0, 100.0, 1000.0, 10000.0)
const CANESSA_ACTION_NAMES = (:a, :b, :c, :d, :e, :f, :g, :h)

_canessa_logistic(x) = 1 / (1 + exp(-x))

# --- Archived parameter draws -------------------------------------------------------------------

"""
    CanessaDraw

One of the 1,000 archived parameter sets of Canessa et al. (2018) (Dryad
doi:10.5061/dryad.0jh5p, file `Canessa_JAPPL_2017_00978_SI2.xls`). Fields:

- `id`: row number, 1-1000;
- `survival`: 3-day survival probabilities at loads 0, 10, 100, 1000 and 10000 GE;
- `transmission`: per-encounter transmission probabilities at the same loads;
- `Lmax`: maximum load per animal (GE); `growth`: zoospores produced by one zoospore in
  5 days;
- `β0`: encounter rate (ha per 3-day step), derived from Stegen et al. (2017);
- `S0`: initial density of susceptible adults (ha⁻¹), paired with `β0`.

The survival, transmission, `Lmax` and `growth` columns are draws from the pooled expert
elicitation (beta-PERT distributions). See [`load_canessa2018`](@ref).
"""
struct CanessaDraw
    id::Int
    survival::NTuple{5,Float64}
    transmission::NTuple{5,Float64}
    Lmax::Float64
    growth::Float64
    β0::Float64
    S0::Float64
end

"""
    canessa2018_path()

Path of the bundled table of the 1,000 archived parameter draws (Dryad
doi:10.5061/dryad.0jh5p, CC0). See `data/canessa2018/README.md`.
"""
canessa2018_path() = joinpath(pkgdir(@__MODULE__), "data", "canessa2018", "parameter_draws.csv")

const _CANESSA_DRAWS = Ref{Union{Nothing,Vector{CanessaDraw}}}(nothing)

"""
    load_canessa2018() -> Vector{CanessaDraw}

The 1,000 archived parameter draws ([`canessa2018_path`](@ref)), in row order, as
[`CanessaDraw`](@ref)s. The result is cached.
"""
function load_canessa2018()
    _CANESSA_DRAWS[] === nothing || return _CANESSA_DRAWS[]::Vector{CanessaDraw}
    lines = filter(!isempty ∘ strip, readlines(canessa2018_path()))
    header = split(lines[1], ',')
    cols = ("draw", "S_0", "S_10", "S_100", "S_1000", "S_10000", "T_0", "T_10", "T_100",
            "T_1000", "T_10000", "L_max", "L_growth", "beta0", "S0")
    idx = map(c -> (i = findfirst(==(c), header);
                    i === nothing ? throw(ArgumentError("column $c missing")) : i), cols)
    draws = map(lines[2:end]) do l
        v = split(l, ',')
        x(k) = parse(Float64, v[idx[k]])
        CanessaDraw(round(Int, x(1)), ntuple(j -> x(1 + j), 5), ntuple(j -> x(6 + j), 5),
            x(12), x(13), x(14), x(15))
    end
    _CANESSA_DRAWS[] = draws
    return draws
end

# --- Management actions (Table 1) -----------------------------------------------------------

"""
    CanessaAction(name; label = string(name), β_factor = 1, β_low_factor = 1,
                  low_load = 100, growth_factor = 1, density_factor = 1, removal = 0)

A management action of Canessa et al. (2018, Table 1), as changes to the IPM:

- `β_factor`: multiplies the transmission function `β(z)` at every load;
- `β_low_factor`: multiplies `β(z)` for loads at or below `low_load` (GE);
- `growth_factor`: multiplies the elicited zoospores-per-zoospore value before the growth
  function is refitted (Table 1's `λ = (⁵√(0.2·data))³` is `growth_factor = 0.2`);
- `density_factor`: multiplies the initial density `S0` (removal before Bsal arrives);
- `removal`: the proportion of all adults removed at each step after Bsal arrives, so that
  `s0` and `s(z)` are multiplied by `1 - removal`.

[`canessa_actions`](@ref) builds the eight actions of Table 1 from the compendium.
"""
struct CanessaAction
    name::Symbol
    label::String
    β_factor::Float64
    β_low_factor::Float64
    low_load::Float64
    growth_factor::Float64
    density_factor::Float64
    removal::Float64
    function CanessaAction(name, label, β_factor, β_low_factor, low_load, growth_factor,
            density_factor, removal)
        all(>=(0), (β_factor, β_low_factor, low_load, density_factor)) ||
            throw(ArgumentError("action factors must be non-negative"))
        growth_factor > 0 || throw(ArgumentError("growth_factor must be positive"))
        0 <= removal < 1 || throw(ArgumentError("removal must lie in [0, 1)"))
        return new(name, label, β_factor, β_low_factor, low_load, growth_factor,
            density_factor, removal)
    end
end

CanessaAction(name::Symbol; label::AbstractString = string(name), β_factor = 1.0,
    β_low_factor = 1.0, low_load = 100.0, growth_factor = 1.0, density_factor = 1.0,
    removal = 0.0) = CanessaAction(name, String(label), β_factor, β_low_factor, low_load,
    growth_factor, density_factor, removal)

"""
    canessa_actions(key = :canessa2018) -> Vector{CanessaAction}

The eight exploratory actions of Table 1, with values from the compendium:

- `a` no action;
- `b` improve body condition: `0.5 β(z)` for `z ≤ 100 GE`;
- `c` probiotic: as `b`, and Bsal growth slowed by 80% (`growth_factor = 0.2`);
- `d` antifungal, full coverage: `(1 - 0.98) β(z)` and growth slowed by 80%;
- `e` antifungal, 80% coverage: `(1 - 0.8 · 0.98) β(z)` and growth slowed by 80%;
- `f`, `g` pre-emptive thinning: `0.5 S0` and `0.1 S0`;
- `h` post-detection removal: 90% of all adults per step (`0.1 s(z)`, and `0.1 s0`).
"""
function canessa_actions(key::Symbol = :canessa2018)
    p = parameter_set(key)
    slow, eff, cov = p.growth_factor, p.antifungal_efficacy, p.antifungal_coverage
    low, half = p.low_load_threshold, p.body_condition_factor
    return [CanessaAction(:a; label = "No action"),
        CanessaAction(:b; label = "Improve body condition", β_low_factor = half, low_load = low),
        CanessaAction(:c; label = "Probiotic treatment", β_low_factor = half, low_load = low,
            growth_factor = slow),
        CanessaAction(:d; label = "Antifungal, 100% coverage", β_factor = 1 - eff,
            growth_factor = slow),
        CanessaAction(:e; label = "Antifungal, 80% coverage", β_factor = 1 - cov * eff,
            growth_factor = slow),
        CanessaAction(:f; label = "Pre-emptive thinning 50%", density_factor = 1 - p.thinning_low),
        CanessaAction(:g; label = "Pre-emptive thinning 90%", density_factor = 1 - p.thinning_high),
        CanessaAction(:h; label = "Post-detection removal", removal = p.removal_rate)]
end

"""
    canessa_action(name::Symbol; key = :canessa2018) -> CanessaAction

The Table 1 action `name` (`:a` to `:h`); see [`canessa_actions`](@ref).
"""
function canessa_action(name::Symbol; key::Symbol = :canessa2018)
    name in CANESSA_ACTION_NAMES ||
        throw(ArgumentError("action must be one of $CANESSA_ACTION_NAMES"))
    return canessa_actions(key)[findfirst(==(name), CANESSA_ACTION_NAMES)]
end

const _NO_ACTION = CanessaAction(:a; label = "No action")

# --- Fits ------------------------------------------------------------------------------------

# Logistic regression of proportions y on z (binomial likelihood with fractional responses,
# as Matlab's glmfit with 'binomial'), by iteratively reweighted least squares.
function _canessa_logistic_fit(z::AbstractVector, y::AbstractVector; maxiter = 100,
        tol = 1e-12)
    all(p -> 0 <= p <= 1, y) || throw(ArgumentError("proportions must lie in [0, 1]"))
    b = [0.0, 0.0]
    X = hcat(ones(length(z)), z)
    for _ in 1:maxiter
        η = X * b
        μ = _canessa_logistic.(η)
        w = μ .* (1 .- μ)
        W = Diagonal(w)
        bnew = (X' * W * X) \ (X' * (W * η .+ (y .- μ)))
        done = maximum(abs, bnew .- b) < tol
        b = bnew
        done && break
    end
    return (b[1], b[2])
end

_canessa_load_points(z_zero::Real) = [z_zero; log.(collect(CANESSA_LOADS[2:end]))]

"""
    canessa_growth_coefficients(Lmax, growth; steps = 10, L0 = 10, step_days = 3,
                                growth_days = 5) -> (b4, b5, b6)

The expected next log-load `z' = b4 + b5 z + b6 z²` of eq. S9. A Beverton-Holt series
`L(t+1) = λ L(t)/(1 + (λ - 1) L(t)/Lmax)` (eq. S8) with `λ = (growth^(1/growth_days))^step_days`
is run for `steps` steps from `L0`, and `ln L(t+1)` is regressed on `ln L(t)` (quadratic,
least squares). The paper does not give `steps`; 10 steps with `Lmax = 10⁴` and
`growth = 25` reproduce Table S1 (1.0194, 1.5288, -0.0692).
"""
function canessa_growth_coefficients(Lmax::Real, growth::Real; steps::Integer = 10,
        L0::Real = 10.0, step_days::Real = 3, growth_days::Real = 5)
    λ = (growth^(1 / growth_days))^step_days
    λ == 1 && throw(ArgumentError("growth = 1 gives a constant series; nothing to fit"))
    steps >= 3 || throw(ArgumentError("need at least 3 steps to fit a quadratic"))
    L = zeros(steps + 1)
    L[1] = L0
    for t in 1:steps
        L[t+1] = λ * L[t] / (1 + (λ - 1) * L[t] / Lmax)
    end
    z = log.(L)
    x, y = z[1:end-1], z[2:end]
    b = hcat(ones(steps), x, x .^ 2) \ y
    return (b[1], b[2], b[3])
end

# --- Parameters ----------------------------------------------------------------------------

"""
    CanessaParameters(; b0, b1, s0 = logistic(b0), b2, b3, β0, μG0, σ2G0, b4, b5, b6, Lmax,
                      σ2G, b7, S0, β_low_factor = 1, low_load = 100, removal = 0,
                      step_days = 3)
    CanessaParameters(:canessa2018; action = canessa_action(:a))
    CanessaParameters(d::CanessaDraw; action = canessa_action(:a), s0 = :fit, key = :canessa2018)

The functions of the Canessa et al. (2018) integral projection model for one parameter set,
on the log-load scale `z = ln(load in GE)`. Probabilities and `β0` are per step of
`step_days = 3` days, and densities per hectare.

- Survival (eq. S4): `s(z) = (1 - removal) logistic(b0 + b1 z)`; susceptibles survive with
  `s0` (already multiplied by `1 - removal`).
- Transmission (eqs. S3, S5): `β(z) = β0 logistic(b2 + b3 z)`, times `β_low_factor` for
  `z ≤ ln(low_load)`; `β0` is in ha per step.
- Initial load (eq. S7): `G0 = N(μG0, σ2G0)`.
- Growth (eqs. S9, S10): `z' ~ N(b4 + b5 z + b6 z², σ2G exp(b7 z))` for `z < ln Lmax`, and
  `N(ln Lmax, σ2G exp(b7 ln Lmax))` above.
- `S0`: initial density of susceptible adults (ha⁻¹).

The keyed constructor gives Table S1, the fits to the experts' best estimates (`s0 = 0.9906`
as printed). An action with `growth_factor ≠ 1` refits `b4`-`b6` from the best estimates
`Lmax = 10⁴` and 25 zoospores per zoospore, inferred from Table S1.

The `CanessaDraw` constructor refits every function to one archived draw, as each of the
paper's 1,000 runs did: logistic regressions of the elicited survival and transmission
probabilities on `z` (load 0 at `z = -1`), and eq. S9 from the draw's `Lmax` and growth
([`canessa_growth_coefficients`](@ref)). `s0 = :fit` (default) takes `s(0)` of the survival
fit, as Table S1 describes; `s0 = :elicited` takes the elicited survival at load 0. The
action's factors are applied on top. `G0` and the variance come from the compendium.
"""
struct CanessaParameters
    b0::Float64
    b1::Float64
    s0::Float64
    b2::Float64
    b3::Float64
    β0::Float64
    μG0::Float64
    σ2G0::Float64
    b4::Float64
    b5::Float64
    b6::Float64
    Lmax::Float64
    σ2G::Float64
    b7::Float64
    S0::Float64
    β_low_factor::Float64
    low_load::Float64
    removal::Float64
    step_days::Float64
    function CanessaParameters(b0, b1, s0, b2, b3, β0, μG0, σ2G0, b4, b5, b6, Lmax, σ2G, b7,
            S0, β_low_factor, low_load, removal, step_days)
        0 <= s0 <= 1 || throw(ArgumentError("s0 must lie in [0, 1]"))
        all(>=(0), (β0, S0, β_low_factor, low_load)) ||
            throw(ArgumentError("β0, S0, β_low_factor and low_load must be non-negative"))
        all(>(0), (σ2G0, σ2G, Lmax, step_days)) ||
            throw(ArgumentError("variances, Lmax and step_days must be positive"))
        0 <= removal < 1 || throw(ArgumentError("removal must lie in [0, 1)"))
        return new(b0, b1, s0, b2, b3, β0, μG0, σ2G0, b4, b5, b6, Lmax, σ2G, b7, S0,
            β_low_factor, low_load, removal, step_days)
    end
end

CanessaParameters(; b0, b1, s0 = _canessa_logistic(b0), b2, b3, β0, μG0, σ2G0, b4, b5, b6, Lmax, σ2G,
    b7, S0, β_low_factor = 1.0, low_load = 100.0, removal = 0.0, step_days = 3.0) =
    CanessaParameters(b0, b1, s0, b2, b3, β0, μG0, σ2G0, b4, b5, b6, Lmax, σ2G, b7, S0,
        β_low_factor, low_load, removal, step_days)

function _canessa_apply(p::NamedTuple, a::CanessaAction)
    keep = 1 - a.removal
    return CanessaParameters(; p..., s0 = p.s0 * keep, β0 = p.β0 * a.β_factor,
        S0 = p.S0 * a.density_factor, β_low_factor = a.β_low_factor, low_load = a.low_load,
        removal = a.removal)
end

function CanessaParameters(key::Symbol; action::CanessaAction = _NO_ACTION)
    c = parameter_set(key)
    b4, b5, b6 = action.growth_factor == 1 ? (c.b4, c.b5, c.b6) :
        canessa_growth_coefficients(c.Lmax_best, c.growth_best * action.growth_factor;
            steps = round(Int, c.growth_series_steps), L0 = c.L0_growth,
            step_days = c.step_length, growth_days = c.growth_period)
    base = (; b0 = c.b0, b1 = c.b1, s0 = c.s0, b2 = c.b2, b3 = c.b3, β0 = c.β0, μG0 = c.μG0,
        σ2G0 = c.σ2G0, b4, b5, b6, Lmax = c.Lmax_best, σ2G = c.σ2G, b7 = c.b7, S0 = c.S0,
        step_days = c.step_length)
    return _canessa_apply(base, action)
end

function CanessaParameters(d::CanessaDraw; action::CanessaAction = _NO_ACTION,
        s0::Symbol = :fit, key::Symbol = :canessa2018)
    s0 in (:fit, :elicited) || throw(ArgumentError("s0 must be :fit or :elicited"))
    c = parameter_set(key)
    z = _canessa_load_points(c.z_zero_load)
    b0, b1 = _canessa_logistic_fit(z, collect(d.survival))
    b2, b3 = _canessa_logistic_fit(z, collect(d.transmission))
    b4, b5, b6 = canessa_growth_coefficients(d.Lmax, d.growth * action.growth_factor;
        steps = round(Int, c.growth_series_steps), L0 = c.L0_growth, step_days = c.step_length,
        growth_days = c.growth_period)
    base = (; b0, b1, s0 = s0 === :fit ? _canessa_logistic(b0) : d.survival[1], b2, b3, β0 = d.β0,
        μG0 = c.μG0, σ2G0 = c.σ2G0, b4, b5, b6, Lmax = d.Lmax, σ2G = c.σ2G, b7 = c.b7,
        S0 = d.S0, step_days = c.step_length)
    return _canessa_apply(base, action)
end

time_unit(::Type{<:CanessaParameters}) = :day

# --- Vital functions -------------------------------------------------------------------------

"""
    canessa_survival(p::CanessaParameters, z)

The 3-day survival probability of an infected adult with log-load `z` (eq. S4), including
any per-step removal.
"""
canessa_survival(p::CanessaParameters, z::Real) = (1 - p.removal) * _canessa_logistic(p.b0 + p.b1 * z)

"""
    canessa_transmission(p::CanessaParameters, z)

The load-dependent transmission function `β(z) = β0 logistic(b2 + b3 z)` of eq. S5 (ha per
step), with the factor of actions b and c for `z ≤ ln(low_load)`. The probability that a
susceptible is infected in a step is `1 - exp(-∫ β(z) I(z) dz)` (eq. S3).
"""
canessa_transmission(p::CanessaParameters, z::Real) =
    p.β0 * _canessa_logistic(p.b2 + p.b3 * z) * (z <= log(p.low_load) ? p.β_low_factor : 1.0)

"""
    canessa_initial_load(p::CanessaParameters, z)

The density `G0(z)` of the log-load of a newly infected adult (eq. S7), normal with mean
`μG0` and variance `σ2G0`.
"""
canessa_initial_load(p::CanessaParameters, z::Real) =
    exp(-(z - p.μG0)^2 / (2p.σ2G0)) / sqrt(2π * p.σ2G0)

"""
    canessa_growth_mean(p::CanessaParameters, z)

The expected log-load one step later, `b4 + b5 z + b6 z²` (eq. S9), or `ln Lmax` once
`z ≥ ln Lmax` (the cap of SI lines 116-118).
"""
canessa_growth_mean(p::CanessaParameters, z::Real) =
    z >= log(p.Lmax) ? log(p.Lmax) : p.b4 + p.b5 * z + p.b6 * z^2

"""
    canessa_growth_variance(p::CanessaParameters, z)

The variance of the log-load one step later, `σ2G exp(b7 z)` (eq. S10, from Wilber et al.
2016), evaluated at `min(z, ln Lmax)`.
"""
canessa_growth_variance(p::CanessaParameters, z::Real) = p.σ2G * exp(p.b7 * min(z, log(p.Lmax)))

# --- Discretised kernel ------------------------------------------------------------------------

function _canessa_mesh(mesh)
    lo, hi, n = mesh
    n >= 2 && hi > lo || throw(ArgumentError("mesh must be (lower, upper, n) with n ≥ 2"))
    h = (hi - lo) / n
    return [lo + (i - 0.5) * h for i in 1:n], h
end

const _CANESSA_DEFAULT_MESH = (-10.0, 15.0, 300)
const CANESSA_CAPS = (:printed, :truncate)

"""Midpoint-rule discretisation: mesh `z`, width `h`, `P[i,j] = h s(z_j) G(z_i, z_j)`,
`g0[i] = h G0(z_i)` and `β[i] = β(z_i)`. With `cap = :truncate` no load exceeds `ln Lmax`:
`G(·, z_j)` is cut at the last mesh point below `ln Lmax` and rescaled to its mass on the
mesh."""
function _canessa_kernel(p::CanessaParameters, mesh, cap::Symbol = :truncate)
    cap in CANESSA_CAPS || throw(ArgumentError("cap must be one of $CANESSA_CAPS"))
    z, h = _canessa_mesh(mesh)
    n = length(z)
    top = cap === :truncate ? something(findlast(<=(log(p.Lmax)), z), 0) : n
    top >= 1 || throw(ArgumentError("ln Lmax lies below the mesh"))
    P = zeros(n, n)
    @inbounds for j in 1:n
        μ = canessa_growth_mean(p, z[j])
        v = canessa_growth_variance(p, z[j])
        c = h / sqrt(2π * v)
        full = 0.0
        for i in 1:n
            g = c * exp(-(z[i] - μ)^2 / (2v))
            full += g
            i <= top && (P[i, j] = g)
        end
        kept = cap === :truncate ? sum(@view P[1:top, j]) : full
        f = kept > 0 ? canessa_survival(p, z[j]) * full / kept : 0.0
        for i in 1:top
            P[i, j] *= f
        end
    end
    g0 = h .* canessa_initial_load.(Ref(p), z)
    β = canessa_transmission.(Ref(p), z)
    return (; z, h, P, g0, β)
end

# y = P x, in plain Julia (multithreaded BLAS on 300 × 300 matrices is slow and erratic on a
# loaded machine, and these are called thousands of times).
function _canessa_mul(P::Matrix{Float64}, x::Vector{Float64})
    n = length(x)
    y = zeros(n)
    @inbounds for j in 1:n
        xj = x[j]
        xj == 0 && continue
        @simd for i in 1:n
            y[i] += P[i, j] * xj
        end
    end
    return y
end

# Solve (I - P) x = b by Gaussian elimination with partial pivoting, in plain Julia.
function _canessa_fundamental(P::Matrix{Float64}, b::Vector{Float64})
    n = length(b)
    A = -P
    @inbounds for i in 1:n
        A[i, i] += 1
    end
    x = copy(b)
    l = zeros(n)
    @inbounds for k in 1:n
        piv, r = abs(A[k, k]), k
        for i in k+1:n
            abs(A[i, k]) > piv && ((piv, r) = (abs(A[i, k]), i))
        end
        piv > 0 || throw(ErrorException("I - P is singular"))
        if r != k
            for j in k:n
                A[k, j], A[r, j] = A[r, j], A[k, j]
            end
            x[k], x[r] = x[r], x[k]
        end
        for i in k+1:n
            l[i] = A[i, k] / A[k, k]
        end
        for j in k+1:n
            a = A[k, j]
            a == 0 && continue
            @simd for i in k+1:n
                A[i, j] -= l[i] * a
            end
        end
        for i in k+1:n
            x[i] -= l[i] * x[k]
        end
    end
    @inbounds for k in n:-1:1
        s = x[k]
        for j in k+1:n
            s -= A[k, j] * x[j]
        end
        x[k] = s / A[k, k]
    end
    return x
end

_canessa_R0(p::CanessaParameters, k) = p.S0 * p.s0 * dot(k.β, _canessa_fundamental(k.P, k.g0))

"""
    basic_reproduction_number(p::CanessaParameters; mesh = (-10, 15, 300), cap = :truncate)

The basic reproduction number of the IPM, `R0 = S0 s0 ∫ β(z) [(I - P)⁻¹ G0](z) dz`: the
infections caused over its life by an adult that enters the infected class with a load from
`G0`, in a population of `S0` susceptibles (Wilber et al. 2016, Appendix S1; the next-
generation operator has rank one). `P` is the survival-growth kernel `s(z) G(z', z)`,
discretised by the midpoint rule on `mesh = (lower, upper, n)` (the paper's -10, 15 and 300).

`cap` sets how the maximum load `Lmax` acts:

- `:truncate` (default): no load exceeds `ln Lmax`. The growth kernel `G(·, z)` is cut at the
  last mesh point below `ln Lmax` and rescaled to its mass on the mesh, so survivors at the
  cap stay there. The SI does not print this, but it takes literally "the hypothesis of a
  maximum possible load per animal", and it reproduces Fig. 2 (median R0 9.4 against 9.6,
  thinning 50% 4.7 against 4.8, thinning 90% 0.94 against 0.93).
- `:printed`: the formula of SI lines 116-118 alone. Above `ln Lmax` the expected next load
  is `ln Lmax`, but the normal kernel lets loads exceed it, where survival is lower. The
  median R0 is then 8.2, uniformly about 14% below the paper's.
"""
function basic_reproduction_number(p::CanessaParameters; mesh = _CANESSA_DEFAULT_MESH,
        cap::Symbol = :truncate)
    return _canessa_R0(p, _canessa_kernel(p, mesh, cap))
end

# --- The IPM problem -----------------------------------------------------------------------------

const CANESSA_UNITS = (:transect, :density)

# Area (ha) of the Stegen et al. (2017) transect, 475 m × 10 m.
function _canessa_transect_area(key::Symbol = :canessa2018)
    c = parameter_set(key)
    return c.transect_length * c.search_width / 10_000
end

# Initial susceptibles and default introduction for the simulation units (see
# CanessaIPMProblem).
function _canessa_units(p::CanessaParameters, units::Symbol, I0)
    units in CANESSA_UNITS || throw(ArgumentError("units must be one of $CANESSA_UNITS"))
    c = parameter_set(:canessa2018)
    if units === :transect
        return (S0 = p.S0 * _canessa_transect_area(), I0 = I0 === nothing ? c.I0_count : I0)
    else
        return (S0 = p.S0, I0 = I0 === nothing ? c.I0 / 100 : I0)   # km⁻² to ha⁻¹
    end
end

"""
    CanessaIPMProblem(p::CanessaParameters; units = :transect, I0 = nothing, nsteps = 30,
                      mesh = (-10, 15, 300), cap = :truncate)

The Canessa et al. (2018) integral projection model (eq. S2) for parameters `p`, starting
from `S0` susceptible and `I0` infected adults with loads distributed as `G0`, for `nsteps`
steps of 3 days. The integrals use the midpoint rule on `mesh = (lower, upper, n)`; `cap` is
as in [`basic_reproduction_number`](@ref).

`units` sets what the transmission term `∫ β(z) I(z) dz` counts:

- `:transect` (default): numbers of adults on the 0.475-ha transect of Stegen et al. (2017)
  (475 m × 10 m), from which `S0` and `β0` were derived. The population starts with
  `S0 × 0.475` susceptibles and, by default, `I0 = 1` infected adult. This is an inference
  (not printed): eq. S6 gives the printed `β0 = 0.007` only when it is solved with the
  transect count (see [`final_size_transmission`](@ref)), and simulating in the same units
  reproduces the final densities of Fig. 2 (actions e and f).
- `:density`: adults per hectare, as Table S1 labels `S0` (503 ha⁻¹) and `I0` (one per km²,
  `I0 = 0.01` ha⁻¹ by default). With `β0` per hectare, every outbreak runs 1/0.475 = 2.1
  times faster than with `:transect`.

[`basic_reproduction_number`](@ref) keeps `S0` per hectare, as the paper's R0 values do; with
`units = :transect` the outbreak therefore spreads as if `R0` were 0.475 times the reported
value.

`solve(prob)` iterates the model and returns a named tuple: `t` (days), `S` (susceptibles),
`I` (infected), `φ` (infection probability in each step), `z` (mesh) and `n`, the
`length(z) × (nsteps + 1)` matrix of infected adults in each mesh cell, all in `units`.
"""
struct CanessaIPMProblem
    p::CanessaParameters
    I0::Float64
    nsteps::Int
    mesh::Tuple{Float64,Float64,Int}
    cap::Symbol
    units::Symbol
end

function CanessaIPMProblem(p::CanessaParameters; units::Symbol = :transect,
        I0::Union{Nothing,Real} = nothing, nsteps::Integer = 30,
        mesh = _CANESSA_DEFAULT_MESH, cap::Symbol = :truncate)
    u = _canessa_units(p, units, I0)
    u.I0 >= 0 || throw(ArgumentError("I0 must be non-negative"))
    nsteps >= 0 || throw(ArgumentError("nsteps must be non-negative"))
    cap in CANESSA_CAPS || throw(ArgumentError("cap must be one of $CANESSA_CAPS"))
    return CanessaIPMProblem(p, Float64(u.I0), Int(nsteps),
        (Float64(mesh[1]), Float64(mesh[2]), Int(mesh[3])), cap, units)
end

function _canessa_iterate(p::CanessaParameters, k, I0, nsteps; save::Bool = true,
        S0::Real = p.S0)
    n = I0 .* k.g0
    S = Float64(S0)
    Ss, Is, φs = [S], [sum(n)], Float64[]
    ns = [copy(n)]
    for _ in 1:nsteps
        φ = 1 - exp(-dot(k.β, n))
        new = S * p.s0 * φ
        n = _canessa_mul(k.P, n) .+ new .* k.g0
        S = S * p.s0 * (1 - φ)
        push!(Ss, S); push!(Is, sum(n)); push!(φs, φ)
        save && push!(ns, n)
    end
    return (; S = Ss, I = Is, φ = φs, n = save ? reduce(hcat, ns) : nothing)
end

function CommonSolve.solve(prob::CanessaIPMProblem)
    p = prob.p
    k = _canessa_kernel(p, prob.mesh, prob.cap)
    S0 = _canessa_units(p, prob.units, prob.I0).S0
    r = _canessa_iterate(p, k, prob.I0, prob.nsteps; S0)
    return (; t = p.step_days .* (0:prob.nsteps), S = r.S, I = r.I, φ = r.φ, z = k.z, n = r.n)
end

"""
    canessa_density_ratio(p::CanessaParameters, reference::CanessaParameters;
                          units = :transect, I0 = nothing, nsteps = 30,
                          mesh = (-10, 15, 300), cap = :truncate)

The outcome on the y-axis of Fig. 2: the susceptibles after `nsteps` under `p` (an action,
with Bsal introduced as `I0` infected adults), divided by the susceptibles with neither
infection nor action, `S0 s0^nsteps` of `reference` (the same draw without action). `units`
and `I0` are as in [`CanessaIPMProblem`](@ref); the ratio itself has no units.
"""
function canessa_density_ratio(p::CanessaParameters, reference::CanessaParameters;
        units::Symbol = :transect, I0::Union{Nothing,Real} = nothing, nsteps::Integer = 30,
        mesh = _CANESSA_DEFAULT_MESH, cap::Symbol = :truncate)
    k = _canessa_kernel(p, mesh, cap)
    u = _canessa_units(p, units, I0)
    S = _canessa_iterate(p, k, u.I0, nsteps; save = false, S0 = u.S0).S[end]
    return S / (_canessa_units(reference, units, I0).S0 * reference.s0^nsteps)
end

# --- Movement of infected adults (Fig. 3) --------------------------------------------------------

"""
    canessa_cohort_survival(p::CanessaParameters; threshold = 0.01, mesh = (-10, 15, 300),
                            cap = :truncate, maxsteps = 10_000) -> Vector{Float64}

The proportion of a cohort of newly infected adults (loads from `G0`) still alive after
0, 1, 2, ... steps, using the survival-growth kernel, up to and including the first step at
which it falls below `threshold`. The paper's life span is the time at which the cumulative
survival falls below 1% (main text lines 171-174).
"""
function canessa_cohort_survival(p::CanessaParameters; threshold::Real = 0.01,
        mesh = _CANESSA_DEFAULT_MESH, cap::Symbol = :truncate, maxsteps::Integer = 10_000)
    k = _canessa_kernel(p, mesh, cap)
    return _canessa_survivorship(k, threshold, maxsteps)
end

function _canessa_survivorship(k, threshold, maxsteps)
    n = k.g0 ./ sum(k.g0)
    surv = [1.0]
    while surv[end] >= threshold
        length(surv) > maxsteps && throw(ErrorException("cohort survives beyond $maxsteps steps"))
        n = _canessa_mul(k.P, n)
        push!(surv, sum(n))
    end
    return surv
end

"""
    canessa_dispersal(p::CanessaParameters; v = 0.5, threshold = 0.01,
                      mesh = (-10, 15, 300), cap = :truncate)
        -> (; distance, proportion, lifespan, max_distance, mean_distance)

How far infected adults move before dying (Fig. 3): `proportion[k]` of a newly infected
cohort moves at least `distance[k]` metres. Adults move `v` km yr⁻¹ (Table S1, 0.5; a 365-day
year), so `3 v/365` km per step. `lifespan` is the number of steps until cumulative survival
falls below `threshold`, `max_distance` the distance covered in that time, and
`mean_distance` the expected distance (`v` times the expected steps survived).
"""
function canessa_dispersal(p::CanessaParameters; v::Real = 0.5, threshold::Real = 0.01,
        mesh = _CANESSA_DEFAULT_MESH, cap::Symbol = :truncate)
    surv = canessa_cohort_survival(p; threshold, mesh, cap)
    return _canessa_dispersal(surv, p.step_days, v)
end

function _canessa_dispersal(surv, step_days, v)
    d = 1000 * v * step_days / 365           # metres per step
    L = length(surv) - 1
    return (; distance = d .* (0:L), proportion = surv, lifespan = L, max_distance = d * L,
        mean_distance = d * sum(surv[2:end]))
end

# --- Ensemble over the archived draws (Fig. 2, Fig. 3) ---------------------------------------

"""
    canessa_ensemble(action = canessa_action(:a); draws = load_canessa2018(),
                     units = :transect, I0 = nothing, nsteps = 30, mesh = (-10, 15, 300),
                     cap = :truncate, v = 0.5, s0 = :fit)

Runs the IPM once per archived draw, as the paper's 1,000 simulations did, and returns one
named tuple per draw: `id`, `R0`
([`basic_reproduction_number`](@ref)),
`ratio` (final susceptibles over the uninfected, unmanaged population,
[`canessa_density_ratio`](@ref); `units` and `I0` as in [`CanessaIPMProblem`](@ref)), and
from [`canessa_dispersal`](@ref) `lifespan` (steps), `max_distance` and `mean_distance` (m).
Summarise with [`canessa_summary`](@ref).
"""
function canessa_ensemble(action::CanessaAction = _NO_ACTION;
        draws::AbstractVector{CanessaDraw} = load_canessa2018(), units::Symbol = :transect,
        I0::Union{Nothing,Real} = nothing,
        nsteps::Integer = 30, mesh = _CANESSA_DEFAULT_MESH, cap::Symbol = :truncate,
        v::Real = 0.5, s0::Symbol = :fit)
    units in CANESSA_UNITS || throw(ArgumentError("units must be one of $CANESSA_UNITS"))
    out = Vector{NamedTuple{(:id, :R0, :ratio, :lifespan, :max_distance, :mean_distance),
        NTuple{6,Float64}}}(undef, length(draws))
    Threads.@threads for i in eachindex(draws)
        d = draws[i]
        ref = CanessaParameters(d; s0)
        p = action == _NO_ACTION ? ref : CanessaParameters(d; action, s0)
        k = _canessa_kernel(p, mesh, cap)
        u = _canessa_units(p, units, I0)
        S = _canessa_iterate(p, k, u.I0, nsteps; save = false, S0 = u.S0).S[end]
        S_free = _canessa_units(ref, units, I0).S0 * ref.s0^nsteps
        mv = _canessa_dispersal(_canessa_survivorship(k, 0.01, 10_000), p.step_days, v)
        out[i] = (; id = Float64(d.id), R0 = _canessa_R0(p, k),
            ratio = S / S_free, lifespan = Float64(mv.lifespan),
            max_distance = mv.max_distance, mean_distance = mv.mean_distance)
    end
    return out
end

"""
    canessa_summary(results; level = 0.95)

Summary of a [`canessa_ensemble`](@ref) in the form of Fig. 2: for `R0`, `ratio`,
`max_distance` and `mean_distance` a named tuple `(median, lower, upper)` of quantiles over the
draws. Also `P_R0_below_1`, the proportion of draws with `R0 < 1` (the chance of
eradication), and `ratio_at_R0_bounds`, the ratios of the draws whose `R0` lies nearest the
lower and upper quantiles (Fig. 2 draws its ratio bars this way).
"""
function canessa_summary(results::AbstractVector; level::Real = 0.95)
    α = (1 - level) / 2
    q(x) = (; median = median(x), lower = quantile(x, α), upper = quantile(x, 1 - α))
    R0 = [r.R0 for r in results]
    ratio = [r.ratio for r in results]
    near(v) = ratio[argmin(abs.(R0 .- v))]
    qR = q(R0)
    return (; R0 = qR, ratio = q(ratio), max_distance = q([r.max_distance for r in results]),
        mean_distance = q([r.mean_distance for r in results]),
        P_R0_below_1 = count(<(1), R0) / length(R0),
        ratio_at_R0_bounds = (lower_R0 = near(qR.lower), upper_R0 = near(qR.upper)))
end

const _CANESSA_ENSEMBLES = Dict{Tuple{Symbol,Symbol,Symbol},Any}()

# The ensembles of the eight Table 1 actions over all 1,000 draws, computed once per session
# for the validation targets and the tests.
function _canessa_reference(name::Symbol; cap::Symbol = :truncate, units::Symbol = :transect)
    key = (name, cap, units)
    haskey(_CANESSA_ENSEMBLES, key) && return _CANESSA_ENSEMBLES[key]
    r = canessa_ensemble(canessa_action(name); cap, units)
    _CANESSA_ENSEMBLES[key] = (; results = r, summary = canessa_summary(r))
    return _CANESSA_ENSEMBLES[key]
end

# --- Simplified SI model (eqs S11, S12; Fig. 4) --------------------------------------------------

"""
    CanessaSIParameters(; β, S0, s0, s̄, step_days = 3)
    CanessaSIParameters(:canessa2018)

The simplified discrete-time SI model of Canessa et al. (2018, eq. S11), the IPM collapsed
over load: transmission `β` (ha per 3-day step), initial susceptibles `S0` (ha⁻¹), survival of
susceptibles `s0` and the mean survival of infected adults `s̄` per step. The keyed
constructor takes Table S1: `β = 0.007`, `S0 = 503`, `s0 = 0.9906` and `s̄ = 0.542`.
"""
struct CanessaSIParameters
    β::Float64
    S0::Float64
    s0::Float64
    s̄::Float64
    step_days::Float64
    function CanessaSIParameters(β, S0, s0, s̄, step_days)
        β >= 0 && S0 >= 0 || throw(ArgumentError("β and S0 must be non-negative"))
        0 <= s0 <= 1 && 0 <= s̄ < 1 || throw(ArgumentError("need 0 ≤ s0 ≤ 1 and 0 ≤ s̄ < 1"))
        return new(β, S0, s0, s̄, step_days)
    end
end

CanessaSIParameters(; β, S0, s0, s̄, step_days = 3.0) =
    CanessaSIParameters(β, S0, s0, s̄, step_days)

function CanessaSIParameters(key::Symbol)
    c = parameter_set(key)
    return CanessaSIParameters(; β = c.β0, S0 = c.S0, s0 = c.s0, s̄ = c.s_bar,
        step_days = c.step_length)
end

time_unit(::Type{<:CanessaSIParameters}) = :day

function _canessa_reductions(mβ, mS0, msurv)
    all(m -> 0 <= m <= 1, (mβ, mS0, msurv)) ||
        throw(ArgumentError("reductions must lie in [0, 1]"))
    return nothing
end

"""
    basic_reproduction_number(p::CanessaSIParameters; mβ = 0, mS0 = 0, msurv = 0)

Eq. S12: `R0 = β(1 - mβ) S0(1 - mS0) s0(1 - msurv)/(1 - s̄(1 - msurv))`, the reproduction
number of the simplified model after management reduces transmission by `mβ`, the initial
density by `mS0` and the survival of all adults (`s0` and `s̄`) by `msurv`. With no management
it is eq. S11's `S0 s0 β/(1 - s̄)`; SI line 70 prints the same quantity without `s0`.
"""
function basic_reproduction_number(p::CanessaSIParameters; mβ::Real = 0, mS0::Real = 0,
        msurv::Real = 0)
    _canessa_reductions(mβ, mS0, msurv)
    return p.β * (1 - mβ) * p.S0 * (1 - mS0) * p.s0 * (1 - msurv) / (1 - p.s̄ * (1 - msurv))
end

"""
    required_reduction(p::CanessaSIParameters, which; mβ = 0, mS0 = 0, msurv = 0)

The proportional reduction that brings eq. S12 to `R0 = 1` (Fig. 4a-c). `which` names what is
solved for: `:β`, `:S0` or `:survival` (the others held at the given values), or a tuple of
these, all reduced by the same proportion (for example `(:β, :S0, :survival)`). Returns 0 if
`R0 ≤ 1` already, and `NaN` if no reduction in `[0, 1]` suffices.
"""
function required_reduction(p::CanessaSIParameters, which; mβ::Real = 0, mS0::Real = 0,
        msurv::Real = 0)
    targets = which isa Symbol ? (which,) : Tuple(which)
    all(w -> w in (:β, :S0, :survival), targets) ||
        throw(ArgumentError("which must be :β, :S0, :survival or a tuple of them"))
    R(m) = basic_reproduction_number(p; mβ = :β in targets ? m : mβ,
        mS0 = :S0 in targets ? m : mS0, msurv = :survival in targets ? m : msurv)
    R(0.0) <= 1 && return 0.0
    R(1.0) > 1 && return NaN
    lo, hi = 0.0, 1.0
    for _ in 1:200
        mid = (lo + hi) / 2
        R(mid) > 1 ? (lo = mid) : (hi = mid)
        hi - lo < 1e-14 && break
    end
    return (lo + hi) / 2
end

"""
    final_size_transmission(S0, S∞, s̄) -> β

Eq. S6, the final-size relation `ln(S0/S∞) = R0 (1 - S∞/S0)` with `R0 = S0 β/(1 - s̄)`
(SI line 70), solved for `β`: `β = (1 - s̄) ln(S0/S∞)/(S0 - S∞)`. The paper used it with
data from Stegen et al. (2017) to derive the encounter rate `β0`.

`β` is per unit of whatever `S0` and `S∞` count. With Table S1's `S0 = 503` ha⁻¹, `S∞ = 7`
and `s̄ = 0.542` it is 0.0039, not the printed 0.007. With `S0` as the count on the 0.475-ha
transect, `503 × 0.475 = 238.9` adults, and `S∞ = 7`, it is 0.00697. The archived draws agree
with that reading (validation targets `canessa2018_beta0_*`).
"""
function final_size_transmission(S0::Real, S∞::Real, s̄::Real)
    0 < S∞ < S0 || throw(ArgumentError("need 0 < S∞ < S0"))
    0 <= s̄ < 1 || throw(ArgumentError("need 0 ≤ s̄ < 1"))
    return (1 - s̄) * log(S0 / S∞) / (S0 - S∞)
end

"""
    CanessaSIProblem(p::CanessaSIParameters; mβ = 0, mS0 = 0, msurv = 0, I0 = 0.01,
                     nsteps = 30)

The simplified SI model (eq. S11) with the management reductions of eq. S12, from
`S0(1 - mS0)` susceptibles and `I0` infected adults per hectare, for `nsteps` steps of 3 days.
`solve(prob)` returns `(; t, S, I)` with `t` in days; see also [`final_density_ratio`](@ref).
"""
struct CanessaSIProblem
    p::CanessaSIParameters
    mβ::Float64
    mS0::Float64
    msurv::Float64
    I0::Float64
    nsteps::Int
end

function CanessaSIProblem(p::CanessaSIParameters; mβ::Real = 0, mS0::Real = 0,
        msurv::Real = 0, I0::Real = 0.01, nsteps::Integer = 30)
    _canessa_reductions(mβ, mS0, msurv)
    I0 >= 0 && nsteps >= 0 || throw(ArgumentError("I0 and nsteps must be non-negative"))
    return CanessaSIProblem(p, mβ, mS0, msurv, I0, nsteps)
end

function CommonSolve.solve(prob::CanessaSIProblem)
    p = prob.p
    β = p.β * (1 - prob.mβ)
    s0 = p.s0 * (1 - prob.msurv)
    s̄ = p.s̄ * (1 - prob.msurv)
    S, I = [p.S0 * (1 - prob.mS0)], [prob.I0]
    for _ in 1:prob.nsteps
        φ = 1 - exp(-β * I[end])
        push!(I, S[end] * s0 * φ + I[end] * s̄)
        push!(S, S[end] * s0 * (1 - φ))
    end
    return (; t = p.step_days .* (0:prob.nsteps), S, I)
end

"""
    final_density_ratio(prob::CanessaSIProblem; quasi_extinction = 0.0)

Fig. 4d: susceptibles after `prob.nsteps` with management and infection, divided by
susceptibles after the same time with neither (`S0 s0^nsteps`). Ratios below
`quasi_extinction` (the paper uses 0.01) are returned as 0.
"""
function final_density_ratio(prob::CanessaSIProblem; quasi_extinction::Real = 0.0)
    S = solve(prob).S[end]
    r = S / (prob.p.S0 * prob.p.s0^prob.nsteps)
    return r < quasi_extinction ? 0.0 : r
end
