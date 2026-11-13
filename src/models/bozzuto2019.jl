# Bozzuto & Canessa (2019) Glob Ecol Conserv 17:e00551. doi:10.1016/j.gecco.2019.e00551
# A two-compartment Bsal model for fire salamanders with seasonal larval inflow, contact rate,
# transmission proportion and Bsal-induced mortality, analysed with a "Floquet-R0".
# Units: time in decimal years (t = 0 is the start of January); densities per hectare.
# Wiki: lab/wiki/models/bsal-bozzuto2019-seasonal.md
#
# Model (eq. 1, checked against the PDF and the online Appendix A):
#   dS/dt = v(t) - d S - c(t) p(t) S I
#   dI/dt = c(t) p(t) S I - (δ(t) + d) I
# with θ(t) = θ0 (1 - θ1 cos 2π(t - ℓθ)) for p and δ (no lag), a two-peaked contact rate
#   c(t) = (c0/m') (1 - c1 cos 2πt)(1 - c2 cos 2π(t - ℓc)),
# and a step inflow v(t) = v on [t1, t2] (numerical analysis) or the sinusoid
# v0 (1 - cos 2π(t - 1/4)) (analytical approximation, Appendix A2; v0 = d k0 here).
#
# Two points the PDF does not make obvious:
# - Table 1 prints c(t) without the factor 1/m'. Appendix A1 rescales it so that c0 is the
#   midrange of c(t) (m' = (min + max)/2 of c(t)/c0), and eq. (2) is only consistent with the
#   rescaled form. `scaled_contact = false` gives the Table 1 formula literally.
# - Table 1 prints c2 = 0.90 for the hibernating population. The Fig. A4 caption gives 0.95,
#   which alone reproduces Floquet-R0 = 8.27, the numerical 8.33 and the Fig. A5 elasticities.
#   `as_printed = true` uses 0.90.
#
# For this model the infected subsystem is scalar, so the Floquet exponent is the time average
# of c p S - d - δ, and the paper's "Floquet-R0" (∫ c p S dt / ∫ (d + δ) dt) is exactly the
# Bacaër-Guernaoui periodic R0. Its departure from R0 = c0 p0 k0/(d + δ0) comes almost
# entirely from c0 being the midrange rather than the mean of c(t).

const BOZZUTO_POPULATIONS = (:hibernation, :aestivation)
const BOZZUTO_INFLOWS = (:step, :sinusoidal, :constant)

"""
    BozzutoParameters(; c0, c1, c2, ℓc, p0 = 1, p1 = 0, k0, v, t1, t2, d, δ0, δ1)
    BozzutoParameters(:bozzuto2019; population = :hibernation, activity = :table1,
                      as_printed = false)

Parameters of the Bozzuto & Canessa (2019) seasonal Bsal model (eq. 1, Table 1). Time is in
decimal years and densities are per hectare.

- `c0` (ha yr⁻¹), `c1`, `c2`, `ℓc` (yr): the contact rate
  `c(t) = (c0/m') (1 - c1 cos 2πt)(1 - c2 cos 2π(t - ℓc))`. `c1` sets the winter trough, `c2`
  the trough lagged by `ℓc`, and `m'` ([`contact_scaling`](@ref)) makes `c0` the midrange.
- `p0`, `p1`: the proportion of contacts that transmit, `p(t) = p0 (1 - p1 cos 2πt)`.
- `k0` (ha⁻¹): average carrying capacity, used where the disease-free density is held fixed.
- `v` (larvae ha⁻¹ yr⁻¹), `t1`, `t2` (yr): the step larval inflow, `v` on `[t1, t2]`.
- `d` (yr⁻¹): natural mortality, constant.
- `δ0` (yr⁻¹), `δ1`: Bsal-induced mortality `δ(t) = δ0 (1 - δ1 cos 2πt)`.

The keyed constructor reads the compendium. `population` is `:hibernation` (north-western
Europe) or `:aestivation` (coastal NE Spain). `activity = :high` substitutes the `c1`, `c2`
of Fig. A4 (~25% activity in winter and summer). `as_printed = true` uses Table 1's
`c2 = 0.90` for the hibernating population instead of the 0.95 that reproduces the paper's
results (see the compendium notes).
"""
struct BozzutoParameters
    c0::Float64
    c1::Float64
    c2::Float64
    ℓc::Float64
    p0::Float64
    p1::Float64
    k0::Float64
    v::Float64
    t1::Float64
    t2::Float64
    d::Float64
    δ0::Float64
    δ1::Float64
    function BozzutoParameters(c0, c1, c2, ℓc, p0, p1, k0, v, t1, t2, d, δ0, δ1)
        all(>=(0), (c0, p0, k0, v, d, δ0)) ||
            throw(ArgumentError("rates and densities must be non-negative"))
        all(x -> 0 <= x <= 1, (c1, c2, p1, δ1)) ||
            throw(ArgumentError("amplitudes c1, c2, p1 and δ1 must lie in [0, 1]"))
        0 <= ℓc < 1 || throw(ArgumentError("phase lag ℓc must lie in [0, 1)"))
        0 <= t1 <= t2 <= 1 || throw(ArgumentError("need 0 ≤ t1 ≤ t2 ≤ 1"))
        return new(c0, c1, c2, ℓc, p0, p1, k0, v, t1, t2, d, δ0, δ1)
    end
end

BozzutoParameters(; c0, c1, c2, ℓc, p0 = 1.0, p1 = 0.0, k0, v, t1, t2, d, δ0, δ1) =
    BozzutoParameters(c0, c1, c2, ℓc, p0, p1, k0, v, t1, t2, d, δ0, δ1)

function BozzutoParameters(key::Symbol; population::Symbol = :hibernation,
        activity::Symbol = :table1, as_printed::Bool = false)
    population in BOZZUTO_POPULATIONS ||
        throw(ArgumentError("population must be one of $BOZZUTO_POPULATIONS"))
    activity in (:table1, :high) || throw(ArgumentError("activity must be :table1 or :high"))
    sets = ["", string(population)]
    activity === :high && push!(sets, "$(population)_high_activity")
    as_printed && population === :hibernation && activity === :table1 &&
        push!(sets, "hibernation_table1")
    recs = parameter_records(key)
    vals = Dict{Symbol,Float64}()
    for s in sets, r in recs
        r.set == s || continue
        vals[r.symbol] = _time_part(r.units) === nothing ? r.value : in_time_unit(r, :year)
    end
    val(k) = haskey(vals, k) ? vals[k] : throw(ArgumentError("$key has no $k for $population"))
    return BozzutoParameters(; c0 = val(:c0), c1 = val(:c1), c2 = val(:c2), ℓc = val(:ℓc),
        p0 = val(:p0), p1 = val(:p1), k0 = val(:k0), v = val(:v), t1 = val(:t1),
        t2 = val(:t2), d = 1 / val(:lifespan), δ0 = val(:δ0), δ1 = val(:δ1))
end

time_unit(::Type{<:BozzutoParameters}) = :year

"""Return a copy of `p` with fields replaced, e.g. `with(p; k0 = 250.0)`."""
with(p::BozzutoParameters; kw...) =
    BozzutoParameters(; (f => get(kw, f, getfield(p, f)) for f in fieldnames(BozzutoParameters))...)

# --- Seasonal functions ----------------------------------------------------------------------

_contact_shape(c1, c2, ℓc, t) = (1 - c1 * cospi(2t)) * (1 - c2 * cospi(2 * (t - ℓc)))

# Minimum or maximum of a 1-periodic function on [0, 1): a grid, then golden-section search
# around the best grid point.
function _periodic_extremum(f, sense::Symbol; n::Int = 2048)
    s = sense === :max ? -1.0 : 1.0
    g(t) = s * f(t)
    ts = range(0, 1; length = n + 1)[1:n]
    i = argmin(g.(ts))
    a, b = ts[i] - 1 / n, ts[i] + 1 / n
    φ = (sqrt(5) - 1) / 2
    x1, x2 = b - φ * (b - a), a + φ * (b - a)
    f1, f2 = g(x1), g(x2)
    for _ in 1:80
        if f1 < f2
            b, x2, f2 = x2, x1, f1
            x1 = b - φ * (b - a); f1 = g(x1)
        else
            a, x1, f1 = x1, x2, f2
            x2 = a + φ * (b - a); f2 = g(x2)
        end
    end
    return s * min(f1, f2, g(ts[i]))
end

"""
    contact_scaling(p::BozzutoParameters)
    contact_scaling(c1, c2, ℓc)

The factor `m' = (min c(t)/c0 + max c(t)/c0)/2` of the unscaled two-peaked contact rate
`c0 (1 - c1 cos 2πt)(1 - c2 cos 2π(t - ℓc))`. Appendix A1 divides the contact rate by `m'`
so that `c0` lies midway between its minimum and maximum. The paper found no closed form
and computed it numerically, as here (a grid refined by golden-section search). The `m` of
eq. (2) is `2m'`.
"""
contact_scaling(c1::Real, c2::Real, ℓc::Real) =
    (_periodic_extremum(t -> _contact_shape(c1, c2, ℓc, t), :min) +
     _periodic_extremum(t -> _contact_shape(c1, c2, ℓc, t), :max)) / 2
contact_scaling(p::BozzutoParameters) = contact_scaling(p.c1, p.c2, p.ℓc)

"""
    contact_rate(p::BozzutoParameters, t; scale = contact_scaling(p))

The contact rate `c(t) = (c0/scale)(1 - c1 cos 2πt)(1 - c2 cos 2π(t - ℓc))`, in ha yr⁻¹.
The default `scale` is the `m'` of Appendix A1. `scale = 1` gives the formula as printed in
Table 1. Pass a precomputed `scale` when calling it many times.
"""
contact_rate(p::BozzutoParameters, t; scale::Real = contact_scaling(p)) =
    p.c0 * _contact_shape(p.c1, p.c2, p.ℓc, t) / scale

"""
    transmission_probability(p::BozzutoParameters, t)

The proportion of contacts that transmit, `p(t) = p0 (1 - p1 cos 2πt)` (Appendix A1). The
paper sets `p1 = 0`.
"""
transmission_probability(p::BozzutoParameters, t) = p.p0 * (1 - p.p1 * cospi(2t))

"""
    bsal_mortality(p::BozzutoParameters, t)

The Bsal-induced mortality rate `δ(t) = δ0 (1 - δ1 cos 2πt)`, in yr⁻¹. It is lowest at the
start of the year, when infected salamanders live longest (40 days in the north-west).
"""
bsal_mortality(p::BozzutoParameters, t) = p.δ0 * (1 - p.δ1 * cospi(2t))

"""
    larval_inflow(p::BozzutoParameters, t; inflow = :step)

The inflow of metamorphs, in larvae ha⁻¹ yr⁻¹.

- `:step`: Table 1's step function, `v` while `t1 ≤ t mod 1 ≤ t2` and zero otherwise. The
  paper uses it in its numerical analysis.
- `:sinusoidal`: Appendix A2's approximation `v0 (1 - cos 2π(t - 1/4))`, peaking in
  September. The paper does not give `v0`; here `v0 = d k0`, the value that makes the mean
  disease-free density `k0`.
- `:constant`: `d k0`, the non-seasonal model.
"""
function larval_inflow(p::BozzutoParameters, t; inflow::Symbol = :step)
    if inflow === :step
        τ = mod(t, 1)
        return p.t1 <= τ <= p.t2 ? p.v : zero(p.v)
    elseif inflow === :sinusoidal
        return p.d * p.k0 * (1 - cospi(2 * (t - 0.25)))
    elseif inflow === :constant
        return p.d * p.k0
    end
    throw(ArgumentError("inflow must be one of $BOZZUTO_INFLOWS"))
end

"""
    disease_free_susceptibles(p::BozzutoParameters, t; inflow = :step)

The periodic disease-free density `S(t)` (ha⁻¹), the solution of `dS/dt = v(t) - d S` with
`S(0) = S(1)`.

- `:step`: eq. (A1), with `S(0) = A/(1 - e^{-d})` and `A = (v/d)(e^{-d(1-t2)} - e^{-d(1-t1)})`.
  Its mean is `v (t2 - t1)/d`, which is 499.2 for the hibernating population, against
  `k0 = 500`.
- `:sinusoidal`: the exact periodic solution for the sinusoidal inflow,
  `k0 + k0 d (2π cos 2πt - d sin 2πt)/(d² + 4π²)`. Eq. (A2b), `k0 (1 + (d/2π) cos 2πt)`, is
  this with `d² = 0`.
- `:constant`: `k0`. This is what eq. (2) assumes.
"""
function disease_free_susceptibles(p::BozzutoParameters, t; inflow::Symbol = :step)
    d = p.d
    if inflow === :step
        τ = mod(t, 1)
        A = p.v / d * (exp(-d * (1 - p.t2)) - exp(-d * (1 - p.t1)))
        S0 = A / (1 - exp(-d))
        τ < p.t1 && return S0 * exp(-d * τ)
        S1 = S0 * exp(-d * p.t1)
        τ <= p.t2 && return (p.v - exp(-d * (τ - p.t1)) * (p.v - d * S1)) / d
        S2 = (p.v - exp(-d * (p.t2 - p.t1)) * (p.v - d * S1)) / d
        return S2 * exp(-d * (τ - p.t2))
    elseif inflow === :sinusoidal
        Δ = d^2 + 4π^2
        return p.k0 * (1 + d * (2π * cospi(2t) - d * sinpi(2t)) / Δ)
    elseif inflow === :constant
        return p.k0
    end
    throw(ArgumentError("inflow must be one of $BOZZUTO_INFLOWS"))
end

# --- Reproduction numbers ---------------------------------------------------------------------

"""
    basic_reproduction_number(p::BozzutoParameters)

The non-seasonal `R0 = c0 p0 k0/(d + δ0)` of Bozzuto & Canessa (2019, p. 3). Because `c0`
is the midrange of the contact rate, not its mean, this is not the R0 of the time-averaged
model; see [`averaged_reproduction_number`](@ref).
"""
basic_reproduction_number(p::BozzutoParameters) = p.c0 * p.p0 * p.k0 / (p.d + p.δ0)

"""
    instantaneous_reproduction_number(p::BozzutoParameters, t; inflow = :constant,
                                      scale = contact_scaling(p))

`R0(t) = c(t) p(t) S(t)/(d + δ(t))` (p. 3), the grey curves of Fig. 2. The paper writes `k(t)`
for the density. `inflow` picks it as in [`disease_free_susceptibles`](@ref) (`:constant`
is `k0`). The paper says `R0(t)` does not decide whether Bsal invades.
"""
instantaneous_reproduction_number(p::BozzutoParameters, t; inflow::Symbol = :constant,
        scale::Real = contact_scaling(p)) =
    contact_rate(p, t; scale) * transmission_probability(p, t) *
    disease_free_susceptibles(p, t; inflow) / (p.d + bsal_mortality(p, t))

"""
    averaged_reproduction_number(p::BozzutoParameters; inflow = :constant)

The R0 of the model with every coefficient replaced by its annual mean:
`mean(c) mean(p) mean(S)/(d + δ0)`. With `S = k0` and `p1 = 0` this equals the analytical
Floquet-R0 of eq. (2): the infected subsystem is scalar, so its Floquet exponent is a time
average. With `inflow = :step` it differs from the numerical Floquet-R0 only by the
covariance of `c(t)` with `S(t)`.
"""
function averaged_reproduction_number(p::BozzutoParameters; inflow::Symbol = :constant,
        n::Int = 4096)
    m′ = contact_scaling(p)
    ts = ((1:n) .- 0.5) ./ n
    mean_(f) = sum(f, ts) / n
    return mean_(t -> contact_rate(p, t; scale = m′)) * mean_(t -> transmission_probability(p, t)) *
           mean_(t -> disease_free_susceptibles(p, t; inflow)) / (p.d + p.δ0)
end

"""
    seasonal_factor(c1, c2, ℓc)
    seasonal_factor(p::BozzutoParameters)

The ratio Floquet-R0/R0 of eq. (2), `(2 + c1 c2 cos 2πℓc)/m` with `m = 2m'`
([`contact_scaling`](@ref)). It is the mean of the rescaled contact cycle `c(t)/c0`. The paper
finds that it lies between about 2/3 and 4/3 (±33%).
"""
seasonal_factor(c1::Real, c2::Real, ℓc::Real) =
    (2 + c1 * c2 * cospi(2ℓc)) / (2 * contact_scaling(c1, c2, ℓc))
seasonal_factor(p::BozzutoParameters) = seasonal_factor(p.c1, p.c2, p.ℓc)

# ∫₀¹ J(t) dt for the scalar infected equation dI/dt = J(t) I at the disease-free state, by
# integrating z' = J (log of eq. A3's y). Returns (∫ c p S, ∫ (d + δ)).
function _floquet_integrals(p::BozzutoParameters, inflow::Symbol; abstol = 1e-12,
        reltol = 1e-12)
    m′ = contact_scaling(p)
    function rhs!(dz, z, _, t)
        dz[1] = contact_rate(p, t; scale = m′) * transmission_probability(p, t) *
                disease_free_susceptibles(p, t; inflow)
        dz[2] = p.d + bsal_mortality(p, t)
        return nothing
    end
    stops = inflow === :step ? [p.t1, p.t2] : Float64[]
    sol = solve(ODEProblem(rhs!, zeros(2), (0.0, 1.0)), Tsit5(); abstol, reltol,
        tstops = stops, save_everystep = false)
    return sol.u[end][1], sol.u[end][2]
end

"""
    floquet_exponent(p::BozzutoParameters; method = :analytic, inflow = ...)

The Floquet exponent `ξ` of the infected equation linearised at the disease-free state:
`log y(1)` for `dy/dt = (c(t) p(t) S(t) - d - δ(t)) y`, `y(0) = 1` (eq. A3; the period is one
year). Bsal cannot invade when `ξ < 0`.

- `method = :analytic`: eq. (A4a),
  `ξ = (c0 p0 k0/2m')(2 - c1 k1 + c1 p1 - k1 p1 + cos 2πℓc (c1 c2 - c2 k1 + c2 p1 - 3c1c2k1p1/4)) - δ0 - d`,
  with `k1 = 0` for `inflow = :constant` (the paper's choice, eq. A4c) and `k1 = d/2π` for
  `inflow = :sinusoidal` (eq. A2b).
- `method = :numerical`: integrates eq. (A3) for the disease-free `S(t)` of `inflow`
  (default `:step`, eq. A1).
"""
function floquet_exponent(p::BozzutoParameters; method::Symbol = :analytic,
        inflow::Symbol = method === :analytic ? :constant : :step)
    if method === :analytic
        inflow in (:constant, :sinusoidal) ||
            throw(ArgumentError("the analytic exponent needs inflow = :constant or :sinusoidal"))
        k1 = inflow === :sinusoidal ? p.d / 2π : 0.0
        c1, c2, p1 = p.c1, p.c2, p.p1
        Q = 2 - c1 * k1 + c1 * p1 - k1 * p1 +
            cospi(2p.ℓc) * (c1 * c2 - c2 * k1 + c2 * p1 - 3c1 * c2 * k1 * p1 / 4)
        return p.c0 * p.p0 * p.k0 / (2 * contact_scaling(p)) * Q - p.δ0 - p.d
    elseif method === :numerical
        gain, loss = _floquet_integrals(p, inflow)
        return gain - loss
    end
    throw(ArgumentError("method must be :analytic or :numerical"))
end

"""
    floquet_reproduction_number(p::BozzutoParameters; method = :analytic, inflow = ...)

The paper's Floquet-R0: the rearrangement of the Floquet exponent that is below one exactly
when `ξ < 0`, `Floquet-R0 = 1 + ξ/(d + δ0)`.

- `method = :analytic` (default `inflow = :constant`): eq. (2),
  `R0 (2 + c1 c2 cos 2πℓc)/m` when `p1 = 0`.
- `method = :numerical` (default `inflow = :step`): `∫ c p S dt / ∫ (d + δ) dt` over a year,
  with `S(t)` from eq. (A1), as in the paper's numerical check.
- `method = :spectral` (default `inflow = :step`): the Bacaër–Guernaoui periodic R0, the
  `R` for which the infected subsystem with transmission divided by `R` has dominant Floquet
  multiplier one, from [`periodic_reproduction_number`](@ref). For this scalar subsystem it
  equals `:numerical`; it is here as an independent check.
"""
function floquet_reproduction_number(p::BozzutoParameters; method::Symbol = :analytic,
        inflow::Symbol = method === :analytic ? :constant : :step)
    if method === :spectral
        m′ = contact_scaling(p)
        F(t) = fill(contact_rate(p, t; scale = m′) * transmission_probability(p, t) *
                    disease_free_susceptibles(p, t; inflow), 1, 1)
        V(t) = fill(p.d + bsal_mortality(p, t), 1, 1)
        stops = inflow === :step ? [p.t1, p.t2] : Float64[]
        return periodic_reproduction_number(F, V, 1; tstops = stops)
    elseif method === :numerical
        gain, loss = _floquet_integrals(p, inflow)
        return gain / loss
    end
    return 1 + floquet_exponent(p; method, inflow) / (p.d + p.δ0)
end

"""
    target_proportion(p::BozzutoParameters; basis = :floquet, kwargs...)

The proportion of the susceptible density to keep (by removal, say) so that Bsal cannot
invade, capped at one (§2.3):

- `basis = :floquet`: `ϱF = 1/Floquet-R0` (`kwargs` go to
  [`floquet_reproduction_number`](@ref));
- `basis = :nonseasonal`: `ϱ = 1/R0`;
- `basis = :peak`: `1/max_t R0(t)`, which also prevents temporary outbreaks (`kwargs` go
  to [`instantaneous_reproduction_number`](@ref)).
"""
function target_proportion(p::BozzutoParameters; basis::Symbol = :floquet, kwargs...)
    R = if basis === :floquet
        floquet_reproduction_number(p; kwargs...)
    elseif basis === :nonseasonal
        basic_reproduction_number(p)
    elseif basis === :peak
        m′ = contact_scaling(p)
        _periodic_extremum(t -> instantaneous_reproduction_number(p, t; scale = m′, kwargs...),
            :max; n = 8192)
    else
        throw(ArgumentError("basis must be :floquet, :nonseasonal or :peak"))
    end
    return min(1.0, 1 / R)
end

"""
    floquet_elasticities(p::BozzutoParameters; σ = 1e-6)
        -> (; R0, c1, c2, ℓc)

The elasticities `(∂F/∂θ)(θ/F)` of the analytical Floquet-R0 `F` (eq. 2) with respect to `R0`,
`c1`, `c2` and `ℓc`, by the central difference of eq. (A7b) with step `σ`, as in Fig. A5a.
The elasticities of the target proportion `ϱF = 1/F` (Fig. A5b) are their negatives.
"""
function floquet_elasticities(p::BozzutoParameters; σ::Real = 1e-6)
    θ = (basic_reproduction_number(p), p.c1, p.c2, p.ℓc)
    F(x) = x[1] * seasonal_factor(x[2], x[3], x[4])
    F0 = F(θ)
    e = map(1:4) do i
        hi = ntuple(j -> j == i ? θ[j] + σ / 2 : θ[j], 4)
        lo = ntuple(j -> j == i ? θ[j] - σ / 2 : θ[j], 4)
        (F(hi) - F(lo)) / σ * θ[i] / F0
    end
    return (; R0 = e[1], c1 = e[2], c2 = e[3], ℓc = e[4])
end

# --- Generic Floquet tools -----------------------------------------------------------------------

"""
    floquet_multipliers(A, n; period = 1.0, tstops = Float64[], abstol = 1e-10,
                        reltol = 1e-10)

The Floquet multipliers of the linear periodic system `dx/dt = A(t) x` in `n` dimensions:
the eigenvalues of the monodromy matrix `Φ(period)`, where `dΦ/dt = A(t) Φ`, `Φ(0) = I`
(eqs. A6a–c). `A(t)` returns an `n × n` matrix. The disease-free state is stable when every
multiplier has modulus below one. Pass the times where `A` jumps as `tstops`.
"""
function floquet_multipliers(A, n::Integer; period::Real = 1.0, tstops = Float64[],
        abstol = 1e-10, reltol = 1e-10)
    rhs(Φ, _, t) = A(t) * Φ
    Φ0 = _identity(n)
    sol = solve(ODEProblem(rhs, Φ0, (0.0, float(period))), Tsit5(); abstol, reltol,
        tstops = collect(Float64, tstops), save_everystep = false)
    return eigvals(sol.u[end])
end

_identity(n) = [i == j ? 1.0 : 0.0 for i in 1:n, j in 1:n]

"""
    periodic_reproduction_number(F, V, n; period = 1.0, tstops = Float64[],
                                 rtol = 1e-10, kwargs...)

The basic reproduction number of a periodic model (Bacaër & Guernaoui 2006; Wang & Zhao
2008): the `R` at which `dx/dt = (F(t)/R - V(t)) x` has dominant Floquet multiplier one.
`F(t)` is the `n × n` matrix of new infections and `V(t)` of transitions out of the infected
compartments, at the disease-free state. Found by bisection on `log R`. `kwargs` go to
[`floquet_multipliers`](@ref).
"""
function periodic_reproduction_number(F, V, n::Integer; period::Real = 1.0,
        tstops = Float64[], rtol::Real = 1e-10, kwargs...)
    ρ(R) = maximum(abs, floquet_multipliers(t -> F(t) ./ R .- V(t), n; period, tstops,
        kwargs...))
    ρ(Inf) < 1 || throw(ArgumentError("the infected subsystem without transmission is not stable"))
    lo, hi = 1.0, 1.0
    while ρ(lo) <= 1
        lo /= 2
        lo < 1e-12 && return 0.0
    end
    while ρ(hi) > 1
        hi *= 2
    end
    while hi - lo > rtol * hi
        mid = sqrt(lo * hi)
        ρ(mid) > 1 ? (lo = mid) : (hi = mid)
    end
    return sqrt(lo * hi)
end

# --- The ODE problem ------------------------------------------------------------------------------

"""In-place right-hand side of eq. (1); state `u = [S, I]`."""
function bozzuto_rhs!(du, u, q, t)
    p = q.p
    S, I = u
    if q.seasonal
        β = contact_rate(p, t; scale = q.scale) * transmission_probability(p, t)
        δ = bsal_mortality(p, t)
        v = larval_inflow(p, t; inflow = q.inflow)
    else
        β, δ, v = p.c0 * p.p0, p.δ0, p.d * p.k0
    end
    du[1] = v - p.d * S - β * S * I
    du[2] = β * S * I - (δ + p.d) * I
    return nothing
end

"""
    SeasonalBsalProblem(p::BozzutoParameters; I0 = 0.01, u0 = nothing, tspan = (0.0, 1.0),
                        inflow = :step, scaled_contact = true, seasonal = true)

The Bozzuto & Canessa (2019) model, eq. (1), with state `[S, I]` in salamanders per hectare
and time in years from the start of January.

- `u0` defaults to the disease-free `S(t0)` of `inflow` ([`disease_free_susceptibles`](@ref))
  plus `I0` infected. The paper simulates no outbreak, so the default `I0 = 0.01 ha⁻¹` (one
  per km²) is only a convenient small introduction.
- `inflow`: `:step`, `:sinusoidal` or `:constant` ([`larval_inflow`](@ref)).
- `scaled_contact = false` uses Table 1's contact formula without the `1/m'` of Appendix A1.
- `seasonal = false` gives the non-seasonal model: `c0`, `p0`, `δ0` and inflow `d k0`.

`solve(prob)` returns the `ODESolution`. `to_ode_problem(prob)` gives the `ODEProblem`; the
times where a step inflow switches on and off are passed as `tstops`.
"""
struct SeasonalBsalProblem
    p::BozzutoParameters
    u0::Vector{Float64}
    tspan::Tuple{Float64,Float64}
    inflow::Symbol
    scaled_contact::Bool
    seasonal::Bool
end

function SeasonalBsalProblem(p::BozzutoParameters; I0::Real = 0.01, u0 = nothing,
        tspan = (0.0, 1.0), inflow::Symbol = :step, scaled_contact::Bool = true,
        seasonal::Bool = true)
    inflow in BOZZUTO_INFLOWS || throw(ArgumentError("inflow must be one of $BOZZUTO_INFLOWS"))
    t0 = float(tspan[1])
    S0 = seasonal ? disease_free_susceptibles(p, t0; inflow) : p.k0
    u = u0 === nothing ? [S0, float(I0)] : collect(Float64, u0)
    length(u) == 2 || throw(ArgumentError("u0 must be [S, I]"))
    return SeasonalBsalProblem(p, u, (t0, float(tspan[2])), inflow, scaled_contact, seasonal)
end

function to_ode_problem(prob::SeasonalBsalProblem)
    p = prob.p
    q = (; p, scale = prob.scaled_contact ? contact_scaling(p) : 1.0, inflow = prob.inflow,
        seasonal = prob.seasonal)
    t0, tf = prob.tspan
    stops = Float64[]
    if prob.seasonal && prob.inflow === :step
        for k in floor(Int, t0):ceil(Int, tf), τ in (p.t1, p.t2)
            t0 < k + τ < tf && push!(stops, k + τ)
        end
    end
    return ODEProblem(bozzuto_rhs!, copy(prob.u0), prob.tspan, q; tstops = stops)
end

CommonSolve.solve(prob::SeasonalBsalProblem, alg = Tsit5(); abstol = 1e-8, reltol = 1e-8,
    kwargs...) = solve(to_ode_problem(prob), alg; abstol, reltol, kwargs...)
