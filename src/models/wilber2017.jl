# Wilber MQ, Knapp RA, Toothman M, Briggs CJ (2017) Resistance, tolerance and environmental
# transmission dynamics determine host extinction risk in a load-dependent amphibian disease.
# Ecology Letters 20:1169–1181. doi:10.1111/ele.12814
# Wiki: lab/wiki/models/bd-wilber2017-ipm.md
#
# Checked against the published PDF (eqns 1–3, Tables 1–2, Methods, Figs 1–3) and the authors'
# code (github.com/mqwilber/host_extinction, master archive: IPM_functions_for_R.R,
# model_parameters.yml, R0_environmental.R, stochastic_multiseason_simulation.R,
# sensitivity_analysis.R, stable_stage_structure.R, zdecay_plots.R), with the March 2017
# preprint appendix in the same repository standing in for Appendices S1–S4 (HTTP 403).
#
# The hybrid model. State u = [T1, T2, T3, S, n_1 … n_m, Z] (tadpoles of years 1–3,
# uninfected adults, infected adults per load bin, zoospores in the pool). Every 3-day step
# (Fig. 1b; eqns 1–3; Appendix S3 eqns 20–21), at temperature T(t):
#   T_i(t+1)  = s_Ti T_i(t)
#   S(t+1)    = s0 (1 − φ) S + Σ s(x) l(x,T) n(x)
#   n(x′,t+1) = Σ s(x)(1 − l(x,T)) G(x′,x,T) n(x) + s0 φ S G0(x′,T)
#   Z(t+1)    = ν(T) Z + μA Σ eˣ n(x) + μT Σ T_i − ψ
# with φ = 1 − exp(−Λ Δt) (Table 2) computed from the state at t. Once a year, when the
# rising temperature first reaches T_repro, the demographic update follows the disease step
# in the same step (Appendix S3 eqns 22–27): T1 ← pA λ (adult survivors), T2 ← pT1 T1,
# T3 ← pT2 T2, S ← S + [(1−pT1) mT1 T1 + (1−pT2) mT2 T2 + mT3 T3] exp(−A(t)/K).
#
# As coded, and differing from or missing in the text (see the wiki card for the tests):
# - ψ, the zoospores removed by newly infected frogs, is in the code (nonstochastic_step,
#   stochastic_step) although the text says it "had very little effect … we do not consider
#   it further": min(φ_Z S exp(μ0 + σ0²/2), ν Z), with φ_Z the pool-only infection
#   probability. zoospore_removal = false drops it.
# - Year-1 tadpoles are produced by adults that survive the step (pA λ times the survivors of
#   each adult class), not pA λ A(t) as Appendix S4 says.
# - ν(T) is an R smooth.spline through Woodhams et al. (2008, Table 1) zoospore death rates,
#   which with four points selects df = 4: the interpolating natural cubic spline, used here.
# - The R0 code (R0_environmental.R) uses Λ without the factor Δt, shedding ln(eˣ + 1) per
#   infected frog (Appendix S2 says f x; eqn 3 and the simulation use eˣ) and s0 = 1.

const WILBER2017_TRANSMISSION = (:zpool, :zpool_dd, :zpool_fd, :const, :const_dd, :const_fd)
const WILBER2017_FIELDS = (:b00, :b10, :surv_mean, :surv_sd, :b01, :b11, :b21, :v01, :c01,
    :b02, :b12, :b22, :b03, :b13, :v03, :c03, :s0, :β0, :β1, :μA, :μT, :sT1, :sT2, :sT3, :pT1,
    :pT2, :pT3, :mT1, :mT2, :mT3, :pA, :λ, :K, :ν, :Δt)

"""
    Wilber2017Parameters(; b00, b10, surv_mean, surv_sd, b01, b11, b21, v01, c01, b02, b12,
                         b22, b03, b13, v03, c03, s0, transmission = :zpool_dd, β0, β1 = 0,
                         μA, μT, sT1, sT2, sT3, pT1, pT2, pT3, mT1, mT2, mT3, pA, λ, K,
                         ν = NaN, Δt = 3)
    Wilber2017Parameters(:wilber2017; transmission = :zpool_dd, set = "", kwargs...)

Parameters of the Wilber et al. (2017) hybrid model of Bd in *Rana muscosa* (Table 1),
with `x` the natural log of the load in zoospore equivalents per swab and `T` in °C. All
probabilities are per step of `Δt = 3` days.

| part | form | fields |
|:--|:--|:--|
| infected survival `s(x)` | `logit s = b00 + b10 (x − surv_mean)/surv_sd` | `b00`, `b10`, `surv_mean`, `surv_sd` |
| growth `G(x′, x, T)` | Normal, mean `b01 + b11 x + b21 T`, variance `v01 exp(2 c01 x)` | `b01`, `b11`, `b21`, `v01`, `c01` |
| loss of infection `l(x, T)` | `logit l = b02 + b12 x + b22 T` | `b02`, `b12`, `b22` |
| initial load `G0(x′, T)` | Normal, mean `b03 + b13 T`, variance `v03 exp(2 c03 T)` | `b03`, `b13`, `v03`, `c03` |
| uninfected adult survival | constant | `s0` |
| transmission `φ` | `1 − exp(−Λ Δt)`, Table 2 | `transmission`, `β0`, `β1` |
| zoospore pool | `ν(T)` survival; shedding `μA eˣ` per adult, `μT` per tadpole | `ν`, `μA`, `μT` |
| tadpoles | survival per step; not metamorphosing; surviving metamorphosis | `sT1…3`, `pT1…3`, `mT1…3` |
| reproduction | probability, clutch, Ricker recruitment `exp(−A/K)` | `pA`, `λ`, `K` |

The survival regression was fitted to the z-scored load; `surv_mean` and `surv_sd` (7.186
and 2.961 ln ZE) are not printed and come from the bundled Wilber 2016 data.

`transmission` (Table 2; `Λ` per day, `I` infected adults, `A` adults, `Z` zoospores):

- `:zpool`: `Λ = β0 ln(Z + 1)` (dynamic zoospore pool);
- `:zpool_dd`: `Λ = β0 ln(Z + 1) + β1 I` (the best model, Fig. 2);
- `:zpool_fd`: `Λ = β0 ln(Z + 1) + β1 I/A`;
- `:const`, `:const_dd`, `:const_fd`: `Λ = β0 [Z > 0]` plus `β1 I` or `β1 I/A` (a constant
  pool, which the authors' code switches on only when zoospores are present).

`ν = NaN` (default) uses the temperature-dependent zoospore survival
[`zoospore_survival`](@ref); a number fixes it (the "high rate of decay" scenario uses 0.05).

The keyed constructor reads `data/compendium/wilber2017.toml`: Table 1 as printed (`set =
""`) or with the authors' full-precision `s0 = 0.9^(1/122)` and tadpole survivals
(`set = "coded"`), and `β0`, `β1` of the chosen Table 2 model.
"""
struct Wilber2017Parameters
    b00::Float64
    b10::Float64
    surv_mean::Float64
    surv_sd::Float64
    b01::Float64
    b11::Float64
    b21::Float64
    v01::Float64
    c01::Float64
    b02::Float64
    b12::Float64
    b22::Float64
    b03::Float64
    b13::Float64
    v03::Float64
    c03::Float64
    s0::Float64
    transmission::Symbol
    β0::Float64
    β1::Float64
    μA::Float64
    μT::Float64
    sT1::Float64
    sT2::Float64
    sT3::Float64
    pT1::Float64
    pT2::Float64
    pT3::Float64
    mT1::Float64
    mT2::Float64
    mT3::Float64
    pA::Float64
    λ::Float64
    K::Float64
    ν::Float64
    Δt::Float64
    function Wilber2017Parameters(b00, b10, surv_mean, surv_sd, b01, b11, b21, v01, c01, b02,
            b12, b22, b03, b13, v03, c03, s0, transmission, β0, β1, μA, μT, sT1, sT2, sT3,
            pT1, pT2, pT3, mT1, mT2, mT3, pA, λ, K, ν, Δt)
        transmission in WILBER2017_TRANSMISSION ||
            throw(ArgumentError("transmission must be one of $WILBER2017_TRANSMISSION"))
        (v01 > 0 && v03 > 0 && surv_sd > 0) ||
            throw(ArgumentError("v01, v03 and surv_sd must be positive"))
        prob(v) = 0 <= v <= 1
        all(prob, (s0, sT1, sT2, sT3, pT1, pT2, pT3, mT1, mT2, mT3, pA)) ||
            throw(ArgumentError("survival and transition probabilities must lie in [0, 1]"))
        (isnan(ν) || 0 <= ν < 1) || throw(ArgumentError("ν must be NaN or lie in [0, 1)"))
        (β0 >= 0 && β1 >= 0 && μA >= 0 && μT >= 0 && λ >= 0 && K > 0 && Δt > 0) ||
            throw(ArgumentError("β0, β1, μA, μT, λ must be non-negative and K, Δt positive"))
        return new(b00, b10, surv_mean, surv_sd, b01, b11, b21, v01, c01, b02, b12, b22, b03,
            b13, v03, c03, s0, transmission, β0, β1, μA, μT, sT1, sT2, sT3, pT1, pT2, pT3,
            mT1, mT2, mT3, pA, λ, K, ν, Δt)
    end
end

Wilber2017Parameters(; b00, b10, surv_mean, surv_sd, b01, b11, b21, v01, c01, b02, b12, b22,
        b03, b13, v03, c03, s0, transmission::Symbol = :zpool_dd, β0, β1 = 0.0, μA, μT, sT1,
        sT2, sT3, pT1, pT2, pT3, mT1, mT2, mT3, pA, λ, K, ν = NaN, Δt = 3.0) =
    Wilber2017Parameters(b00, b10, surv_mean, surv_sd, b01, b11, b21, v01, c01, b02, b12, b22,
        b03, b13, v03, c03, s0, transmission, β0, β1, μA, μT, sT1, sT2, sT3, pT1, pT2, pT3,
        mT1, mT2, mT3, pA, λ, K, ν, Δt)

function Wilber2017Parameters(key::Symbol; transmission::Symbol = :zpool_dd,
        set::AbstractString = "", kwargs...)
    transmission in WILBER2017_TRANSMISSION ||
        throw(ArgumentError("transmission must be one of $WILBER2017_TRANSMISSION"))
    q = parameter_set(key; set)
    tr = parameter_set(key; set = string(transmission))
    names = filter(s -> s ∉ (:β0, :β1, :ν, :Δt) && haskey(q, s), WILBER2017_FIELDS)
    base = NamedTuple{names}(Tuple(getfield(q, s) for s in names))
    rec = only(filter(r -> r.symbol === :Δt && r.set == "", parameter_records(key)))
    base = merge(base, (; Δt = in_time_unit(rec, :day), transmission, β0 = tr.β0,
        β1 = haskey(tr, :β1) ? tr.β1 : 0.0))
    return Wilber2017Parameters(; merge(base, values(kwargs))...)
end

time_unit(::Type{<:Wilber2017Parameters}) = :day

_with(p::Wilber2017Parameters; kw...) = Wilber2017Parameters(;
    (f => get(kw, f, getfield(p, f)) for f in fieldnames(Wilber2017Parameters))...)

# IPMKernel(p::Wilber2017Parameters; temperature): the within-season kernel (documented
# under IPMKernel).
function IPMKernel(p::Wilber2017Parameters; temperature::Real)
    T = Float64(temperature)
    survival = x -> _logistic(p.b00 + p.b10 * (x - p.surv_mean) / p.surv_sd)
    growth = x -> Normal(p.b01 + p.b11 * x + p.b21 * T, sqrt(p.v01 * exp(2 * p.c01 * x)))
    loss = x -> _logistic(p.b02 + p.b12 * x + p.b22 * T)
    initial = Normal(p.b03 + p.b13 * T, sqrt(p.v03 * exp(2 * p.c03 * T)))
    return IPMKernel(survival, growth, loss, initial, p.s0)
end

# --- Zoospore survival: natural cubic spline through Woodhams et al. (2008) Table 1 -----------

struct _NaturalSpline
    x::Vector{Float64}
    y::Vector{Float64}
    M::Vector{Float64}           # second derivatives at the knots (0 at the ends)
end

function _NaturalSpline(x::AbstractVector, y::AbstractVector)
    n = length(x)
    n >= 3 && issorted(x) || throw(ArgumentError("need at least three increasing knots"))
    h = diff(x)
    A = zeros(n, n)
    r = zeros(n)
    A[1, 1] = A[n, n] = 1
    for i in 2:n-1
        A[i, i-1], A[i, i], A[i, i+1] = h[i-1], 2(h[i-1] + h[i]), h[i]
        r[i] = 6 * ((y[i+1] - y[i]) / h[i] - (y[i] - y[i-1]) / h[i-1])
    end
    return _NaturalSpline(collect(Float64, x), collect(Float64, y), A \ r)
end

# Cubic inside [x1, xn]; linear beyond, with the end slopes (as R's predict.smooth.spline).
function (s::_NaturalSpline)(t::Real)
    x, y, M = s.x, s.y, s.M
    n = length(x)
    if t < x[1] || t > x[n]
        i = t < x[1] ? 1 : n - 1
        h = x[i+1] - x[i]
        slope = (y[i+1] - y[i]) / h - h * (2M[i] + M[i+1]) / 6
        t > x[n] && (slope = (y[n] - y[n-1]) / h + h * (M[n-1] + 2M[n]) / 6)
        return t < x[1] ? y[1] + slope * (t - x[1]) : y[n] + slope * (t - x[n])
    end
    i = clamp(searchsortedlast(x, t), 1, n - 1)
    h = x[i+1] - x[i]
    a, b = x[i+1] - t, t - x[i]
    return M[i] * a^3 / (6h) + M[i+1] * b^3 / (6h) + (y[i] / h - M[i] * h / 6) * a +
           (y[i+1] / h - M[i+1] * h / 6) * b
end

const _WOODHAMS_DEATH = Ref{Union{Nothing,_NaturalSpline}}(nothing)

# Zoospore death rate per day as a function of temperature, from the compendium.
function _woodhams_death_spline()
    _WOODHAMS_DEATH[] === nothing || return _WOODHAMS_DEATH[]::_NaturalSpline
    recs = filter(r -> r.symbol === :dZ && startswith(r.set, "table1"),
        parameter_records(:woodhams2008))
    sort!(recs; by = r -> Float64(r.system["temperature_C"]))
    s = _NaturalSpline([Float64(r.system["temperature_C"]) for r in recs],
        [in_time_unit(r, :day) for r in recs])
    _WOODHAMS_DEATH[] = s
    return s
end

"""
    zoospore_survival(p::Wilber2017Parameters, temperature) -> Float64

Probability `ν(T)` that a zoospore in the pool survives one step of `p.Δt` days at
`temperature` (°C): `exp(−d(T) Δt)`, with `d(T)` the death rate (per day) interpolated by a
natural cubic spline through the four rates of Woodhams et al. (2008, Table 1; compendium
`woodhams2008`, sets `table1_*`), as the authors' `get_zoospore_surv_prob` does with R's
`smooth.spline` (which on four points selects the interpolating spline). Linear
extrapolation outside 4–28 °C. If `p.ν` is a number it is returned instead.

At 4, 12 and 20 °C, `ν` is 0.823, 0.804 and 0.565 per 3 days.
"""
function zoospore_survival(p::Wilber2017Parameters, temperature::Real)
    isnan(p.ν) || return p.ν
    return exp(-_woodhams_death_spline()(temperature) * p.Δt)
end

"""
    infection_probability(p::Wilber2017Parameters; Z, I, A) -> Float64

Probability `φ = 1 − exp(−Λ Δt)` that an uninfected adult becomes infected in one step,
for `Z` zoospores in the pool, `I` infected adults and `A` adults (Table 2; see
[`Wilber2017Parameters`](@ref) for the forms of `Λ`). For frequency dependence, `I/A = 0`
when `A = 0`, as in the authors' code.
"""
function infection_probability(p::Wilber2017Parameters; Z::Real, I::Real = 0.0, A::Real = 0.0)
    return 1 - exp(-_force(p, Z, I, A) * p.Δt)
end

function _force(p::Wilber2017Parameters, Z, I, A)
    tr = p.transmission
    env = tr in (:zpool, :zpool_dd, :zpool_fd) ? p.β0 * log(Z + 1) : p.β0 * (Z > 0)
    host = tr in (:zpool_dd, :const_dd) ? p.β1 * I :
           tr in (:zpool_fd, :const_fd) ? (A > 0 ? p.β1 * I / A : 0.0) : 0.0
    return env + host
end

# --- The hybrid model ------------------------------------------------------------------------

"""
    annual_demography_matrix(p::Wilber2017Parameters; steps_per_year = 122) -> Matrix{Float64}

The 4 × 4 disease-free, density-independent annual projection matrix on
`[T1, T2, T3, adults]` that the authors use for the initial stage structure
(`stable_stage_structure.R`): per-step survivals raised to `steps_per_year`, `pA λ` year-1
tadpoles per adult, the tadpole transitions `pT`, metamorph survival `mT` and recruitment
without density dependence. Its dominant eigenvalue is the "λ_growth rate ≈ 1.46" of Table 1
and its stable distribution, scaled to 10 adults, the initial 85, 12 and 3 tadpoles.
"""
function annual_demography_matrix(p::Wilber2017Parameters; steps_per_year::Integer = 122)
    k = steps_per_year
    sT = (p.sT1^k, p.sT2^k, p.sT3^k)
    pT = (p.pT1, p.pT2, p.pT3)
    mT = (p.mT1, p.mT2, p.mT3)
    A = zeros(4, 4)
    A[1, 4] = p.λ * p.pA
    A[2, 1] = sT[1] * pT[1]
    A[3, 2] = sT[2] * pT[2]
    for i in 1:3
        A[4, i] = sT[i] * (1 - pT[i]) * mT[i]
    end
    A[4, 4] = p.s0^k
    return A
end

"""
    HybridIPMProblem(p::Wilber2017Parameters; mesh = IPMMesh(:wilber2017), years = 25,
                     temperature = nothing, u0 = nothing, stochastic = true,
                     zoospore_removal = true)

The hybrid model of Wilber et al. (2017, Fig. 1; Appendix S3): the within-season IPM of
adult load with three tadpole stages and a zoospore pool, stepped every `p.Δt = 3` days,
and the yearly demographic update. The state is
`u = [T1, T2, T3, S, n_1, …, n_m, Z]` (numbers per mesh bin for infected adults).

- `mesh`: the load mesh, 30 bins on `[−5, 18]` ln ZE by default (compendium).
- `years`: length in years of `days_per_year = 366` days (122 steps; compendium).
- `temperature`: °C as a number or a function of time in days; by default the authors'
  cosine, `T(k) = 12 − 8 cos(2π k/122 + acos((T_start − 12)/−8))` at step `k`, from 4 °C in
  winter to 20 °C in summer ([`SeasonalTemperature`](@ref)).
- `u0`: initial state; by default 10 uninfected adults, 85, 12 and 3 tadpoles and no
  zoospores (compendium).
- `stochastic = true`: demographic stochasticity (Appendix S4): multinomial transitions of
  individuals and zoospores, Poisson births and shedding ([`demographic_projection`](@ref)).
  `false` iterates the mean map.
- `zoospore_removal = true`: the authors' code removes the zoospores taken up by newly
  infected frogs (`ψ` of eqn 3); `false` drops the term as the text describes.

The step to time `k` uses the temperature at time `k` and the infection probability of the
state at `k − 1`; the yearly update happens in the step at which `round(T, 2)` first equals
`T_repro = 11.79` °C in spring (step 30 of each year), as in the authors' code.
`solve(prob; rng)` returns a [`HybridIPMSolution`](@ref); [`extinction_curve`](@ref)
repeats stochastic runs.
"""
struct HybridIPMProblem{TT}
    p::Wilber2017Parameters
    mesh::IPMMesh
    temperature::TT
    u0::Vector{Float64}
    steps::Int
    reproduce::BitVector
    stochastic::Bool
    zoospore_removal::Bool
    days_per_year::Float64
end

"""
    wilber2017_temperature(; key = :wilber2017) -> SeasonalTemperature

The seasonal lake temperature of the hybrid model, as a [`SeasonalTemperature`](@ref) in
days: mean 12 °C, amplitude 8 °C, period 366 days, starting at `T_start = 4` °C (winter), so
that `T(3k)` is the authors' `get_temp(k, 4)` (compendium).
"""
function wilber2017_temperature(; key::Symbol = :wilber2017)
    q = parameter_set(key)
    period = q.days_per_year
    φ0 = acos(clamp((q.T_start - q.T_mean) / -q.T_amplitude, -1, 1))
    return SeasonalTemperature(; mean = q.T_mean, amplitude = q.T_amplitude,
        peak = (π - φ0) * period / (2π), period)
end

function HybridIPMProblem(p::Wilber2017Parameters; mesh::IPMMesh = IPMMesh(:wilber2017),
        years::Union{Nothing,Real} = nothing, temperature = nothing,
        u0::Union{Nothing,AbstractVector} = nothing, stochastic::Bool = true,
        zoospore_removal::Bool = true)
    q = parameter_set(:wilber2017)
    years = years === nothing ? q.years : years
    temperature = temperature === nothing ? wilber2017_temperature() : temperature
    temperature isa Real || applicable(temperature, 0.0) ||
        throw(ArgumentError("temperature must be a number or a function of time (days)"))
    n = mesh.n
    if u0 === nothing
        u0 = vcat(q.T1_0, q.T2_0, q.T3_0, q.A0, zeros(n), q.Z0)
    end
    length(u0) == n + 5 || throw(DimensionMismatch("u0 must have length mesh.n + 5"))
    all(>=(0), u0) || throw(ArgumentError("initial numbers must be non-negative"))
    stochastic && !all(isinteger, u0) &&
        throw(ArgumentError("a stochastic run needs integer initial numbers"))
    steps = round(Int, years * q.days_per_year / p.Δt)
    reproduce = falses(steps)
    flag = true
    for k in 1:steps
        T = _temperature_at(temperature, k * p.Δt)
        if round(T; digits = 2) == q.T_repro
            reproduce[k] = flag
            flag = !flag
        end
    end
    return HybridIPMProblem(p, mesh, temperature, collect(Float64, u0), steps, reproduce,
        stochastic, zoospore_removal, Float64(q.days_per_year))
end

_temperature_at(T::Real, _) = Float64(T)
_temperature_at(T, t) = Float64(T(t))

"""
    HybridIPMSolution

The result of `solve(::HybridIPMProblem)`:

- `t`: times in days, `0, Δt, …`;
- `u`: the state at each time (`(mesh.n + 5) × length(t)`), rows
  `[T1, T2, T3, S, n_1 … n_m, Z]`;
- `temperature`: °C used in each step (`length(t) − 1`);
- `φ`: the infection probability of each step;
- `mesh`: the [`IPMMesh`](@ref).

`hybrid_summary(sol)` gives tadpoles, adults, prevalence, mean load and extinction over time.
"""
struct HybridIPMSolution
    t::Vector{Float64}
    u::Matrix{Float64}
    temperature::Vector{Float64}
    φ::Vector{Float64}
    mesh::IPMMesh
end

"""
    hybrid_summary(sol::HybridIPMSolution) -> NamedTuple

Time series (one entry per time in `sol.t`, in days; `years = t/366`) of `tadpoles` (all
three stages), `adults`, `uninfected`, `infected`, `prevalence` among adults, `mean_load`
(mean ln ZE of infected adults, `NaN` when none), `zoospores` and `extinct` (no tadpoles
and no adults).
"""
function hybrid_summary(sol::HybridIPMSolution)
    n = sol.mesh.n
    U = sol.u
    ys = meshpoints(sol.mesh)
    tad = vec(sum(U[1:3, :]; dims = 1))
    S = U[4, :]
    Iv = U[5:4+n, :]
    I = vec(sum(Iv; dims = 1))
    A = S .+ I
    ml = [I[k] > 0 ? dot(view(Iv, :, k), ys) / I[k] : NaN for k in eachindex(I)]
    return (t = sol.t, years = sol.t ./ 366, tadpoles = tad, adults = A, uninfected = S,
        infected = I, prevalence = [a > 0 ? i / a : NaN for (i, a) in zip(I, A)],
        mean_load = ml, zoospores = U[end, :], extinct = (tad .+ A) .== 0)
end

# Per-temperature pieces of the within-season matrix.
struct _HybridPieces
    d::DiscretizedIPM
    ν::Float64
    g0μ::Float64
    g0σ::Float64
    W::Matrix{Float64}         # within-season matrix without the S column
end

struct _HybridCache
    prob::HybridIPMProblem
    pieces::Dict{Float64,_HybridPieces}
end
_HybridCache(prob::HybridIPMProblem) = _HybridCache(prob, Dict{Float64,_HybridPieces}())

function _pieces(c::_HybridCache, T::Float64)
    get!(c.pieces, T) do
        p, mesh = c.prob.p, c.prob.mesh
        n = mesh.n
        k = IPMKernel(p; temperature = T)
        d = discretize(k, mesh)
        ν = zoospore_survival(p, T)
        N = n + 5
        W = zeros(N, N)
        W[1, 1], W[2, 2], W[3, 3] = p.sT1, p.sT2, p.sT3
        W[4, 5:4+n] .= d.loss
        W[5:4+n, 5:4+n] .= d.P
        ys = meshpoints(mesh)
        W[N, 1:3] .= p.μT
        W[N, 5:4+n] .= p.μA .* exp.(ys)
        W[N, N] = ν
        _HybridPieces(d, ν, k.initial_load.μ, k.initial_load.σ, W)
    end
end

# The full mean matrix of the step to time k, given the state u at k − 1. Returns φ, the
# pool-only infection probability and the pieces.
function _hybrid_matrix!(A, c::_HybridCache, u, k)
    prob = c.prob
    p, n = prob.p, prob.mesh.n
    T = _temperature_at(prob.temperature, k * p.Δt)
    pc = _pieces(c, T)
    S = u[4]
    Itot = sum(@view u[5:4+n])
    Ad = S + Itot
    Z = u[end]
    φ = infection_probability(p; Z, I = Itot, A = Ad)
    φZ = infection_probability(p; Z, I = 0.0, A = Ad)
    copyto!(A, pc.W)
    A[4, 4] = p.s0 * (1 - φ)
    @views A[5:4+n, 4] .= (p.s0 * φ) .* pc.d.recruitment
    prob.reproduce[k] && _demography!(A, p, n, Ad)
    return T, φ, φZ, pc
end

# Left-multiply the within-season matrix by the yearly demography matrix D (Appendix S3
# eqns 22–27; build_demography_matrix): rows T1, T2, T3 and S change, the rest are identity.
function _demography!(A, p::Wilber2017Parameters, n, adults)
    r = exp(-adults / p.K)
    adult_rows = 4:4+n
    newT1 = p.λ * p.pA .* vec(sum(@view(A[adult_rows, :]); dims = 1))
    newT2 = p.pT1 .* A[1, :]
    newT3 = p.pT2 .* A[2, :]
    recruits = ((1 - p.pT1) * p.mT1 * r) .* A[1, :] .+ ((1 - p.pT2) * p.mT2 * r) .* A[2, :] .+
               ((1 - p.pT3) * p.mT3 * r) .* A[3, :]
    A[4, :] .+= recruits
    A[1, :] .= newT1
    A[2, :] .= newT2
    A[3, :] .= newT3
    return A
end

# One step; writes the state at k into `out`. Returns (T, φ).
function _hybrid_step!(out, buf, c::_HybridCache, u, k, rng)
    prob = c.prob
    p, n = prob.p, prob.mesh.n
    N = n + 5
    A, Tm, F, v = buf
    T, φ, φZ, pc = _hybrid_matrix!(A, c, u, k)
    S, Z = u[4], u[N]
    if !prob.stochastic
        mul!(out, A, u)
        if prob.zoospore_removal
            out[N] -= min(φZ * S * exp(pc.g0μ + pc.g0σ^2 / 2), Z * pc.ν)
        end
        return T, φ
    end
    # T: survival and transitions (row T1 only its own survival, row Z only pool survival);
    # F: births of year-1 tadpoles and shedding into the pool (stochastic_step).
    copyto!(Tm, A)
    fill!(F, 0.0)
    @views F[1, 2:N] .= A[1, 2:N]
    @views F[N, 1:N-1] .= A[N, 1:N-1]
    @views Tm[1, 2:N] .= 0
    @views Tm[N, 1:N-1] .= 0
    Tm[N, N] = 0.0
    copyto!(v, u)
    v[N] = 0.0
    demographic_projection!(rng, out, Tm, F, v)
    zleft = Z > 0 ? Float64(rand(rng, Binomial(Int(Z), pc.ν))) : 0.0
    out[N] += zleft
    if prob.zoospore_removal && S > 0
        infected = rand(rng, Binomial(Int(S), φZ))
        gained = 0.0
        for _ in 1:infected
            gained += exp(pc.g0μ + pc.g0σ * randn(rng))
        end
        out[N] -= min(round(gained), zleft)
    end
    return T, φ
end

_hybrid_buffers(N) = (zeros(N, N), zeros(N, N), zeros(N, N), zeros(N))

function CommonSolve.solve(prob::HybridIPMProblem; rng::AbstractRNG = Random.default_rng())
    c = _HybridCache(prob)
    N = prob.mesh.n + 5
    U = Matrix{Float64}(undef, N, prob.steps + 1)
    U[:, 1] .= prob.u0
    Ts = Vector{Float64}(undef, prob.steps)
    φs = similar(Ts)
    buf = _hybrid_buffers(N)
    for k in 1:prob.steps
        Ts[k], φs[k] = @views _hybrid_step!(U[:, k+1], buf, c, U[:, k], k, rng)
    end
    t = collect(0:prob.steps) .* prob.p.Δt
    return HybridIPMSolution(t, U, Ts, φs, prob.mesh)
end

# The deterministic map for a DiscreteProblem.
struct _HybridMap{B}
    cache::_HybridCache
    buf::B
end

function (m::_HybridMap)(du, u, _, t)
    k = round(Int, t) + 1
    _hybrid_step!(du, m.buf, m.cache, u, k, nothing)
    return nothing
end

"""
    to_discrete_problem(prob::HybridIPMProblem) -> SciMLBase.DiscreteProblem

The deterministic hybrid map (`stochastic = false` required) as an in-place
`SciMLBase.DiscreteProblem` on `[T1, T2, T3, S, n_1 … n_m, Z]`, with time in steps
(`tspan = (0, steps)`; multiply by `Δt` for days). Solve it with `FunctionMap()` from
OrdinaryDiffEq, or call `solve(prob)`.
"""
function to_discrete_problem(prob::HybridIPMProblem)
    prob.stochastic && throw(ArgumentError("to_discrete_problem needs stochastic = false"))
    m = _HybridMap(_HybridCache(prob), _hybrid_buffers(prob.mesh.n + 5))
    return DiscreteProblem{true}(m, copy(prob.u0), (0.0, Float64(prob.steps)))
end

"""
    extinction_curve(prob::HybridIPMProblem; nsims = 500, rng = Random.default_rng())
        -> NamedTuple

Time-dependent extinction probability of the hybrid model (Methods, 'Simulating the hybrid
model'; Appendix Fig. S5): the fraction of `nsims` stochastic runs of `prob` in which no
tadpoles and no adults remain, at each time. Returns `years` (`t/366`), `probability` and
`extinction_step` (the first extinct step of each run, `typemax(Int)` if none). Extinction
is absorbing, so each run stops when it occurs.
"""
function extinction_curve(prob::HybridIPMProblem; nsims::Union{Nothing,Integer} = nothing,
        rng::AbstractRNG = Random.default_rng())
    prob.stochastic || throw(ArgumentError("extinction_curve needs a stochastic problem"))
    nsims = nsims === nothing ? round(Int, parameter_set(:wilber2017).n_sims) : nsims
    c = _HybridCache(prob)
    ext = _extinction_steps(c, nsims, rng)
    prob_ext = [count(<=(k), ext) / nsims for k in 0:prob.steps]
    years = collect(0:prob.steps) .* (prob.p.Δt / prob.days_per_year)
    return (years = years, probability = prob_ext, extinction_step = ext)
end

function _extinction_steps(c::_HybridCache, nsims, rng)
    prob = c.prob
    N = prob.mesh.n + 5
    buf = _hybrid_buffers(N)
    u, w = similar(prob.u0), similar(prob.u0)
    ext = fill(typemax(Int), nsims)
    for s in 1:nsims
        u .= prob.u0
        for k in 1:prob.steps
            _hybrid_step!(w, buf, c, u, k, rng)
            u, w = w, u
            if sum(@view u[1:N-1]) == 0
                ext[s] = k
                break
            end
        end
    end
    return ext
end

"""
    basic_reproduction_number(p::Wilber2017Parameters; temperature, S,
                              mesh = IPMMesh(:wilber2017; n = 50), as_coded = true)
        -> NamedTuple

`R0` of Bd invading `S` susceptible adults at a constant `temperature` (°C), with and without
the environmental zoospore pool (Fig. 2; Appendix S2): [`reservoir_reproduction_number`](@ref)
of the within-season IPM, linearising `Λ` at the disease-free state (`ln(Z + 1) ≈ Z`).
Returns `total` (with the pool), `host` (`β0 = 0`), `environment` and the invasion
probabilities `1 − 1/R0` (0 when `R0 ≤ 1`) as `invasion` and `invasion_host`.

- `as_coded = true` (default) follows `R0_environmental.R`, which made Fig. 2: transmission
  `β0 S` and `β1 S` per step without the factor `Δt`, shedding `ln(eˣ + 1)` zoospores per
  infected adult and `s0 = 1`.
- `as_coded = false` uses the model as simulated (eqns 1–3): `Λ Δt`, shedding `μA eˣ`, and
  `s0` of `p`.

Density dependence (`:zpool_dd`, `:const_dd`) gives each infected adult `β1 S` infections per
day; frequency dependence (`:zpool_fd`, `:const_fd`) gives `β1`. The 50-bin mesh is the
authors'.
"""
function basic_reproduction_number(p::Wilber2017Parameters; temperature::Real, S::Real,
        mesh::Union{Nothing,IPMMesh} = nothing, as_coded::Bool = true)
    mesh = mesh === nothing ?
           IPMMesh(:wilber2017; n = round(Int, parameter_set(:wilber2017).mesh_n_R0)) : mesh
    q = as_coded ? _with(p; s0 = 1.0) : p
    d = discretize(IPMKernel(q; temperature), mesh)
    ys = meshpoints(mesh)
    dt = as_coded ? 1.0 : p.Δt
    tr = p.transmission
    host = tr in (:zpool_dd, :const_dd) ? p.β1 * S : tr in (:zpool_fd, :const_fd) ? p.β1 : 0.0
    shed = as_coded ? log.(exp.(ys) .+ 1) : p.μA .* exp.(ys)
    r = reservoir_reproduction_number(d; host_transmission = host * dt,
        pool_transmission = p.β0 * S * dt, shedding = shed,
        pool_survival = zoospore_survival(p, temperature))
    inv(R) = R > 1 ? 1 - 1 / R : 0.0
    return merge(r, (invasion = inv(r.total), invasion_host = inv(r.host)))
end

"""
    stasis_load(p::Wilber2017Parameters; temperature) -> NamedTuple

The graphical argument of Fig. 4: the load at which the mean growth function crosses the
1:1 line, `x* = (b01 + b21 T)/(1 − b11)` (ln ZE per swab; resistance), and the probability
`s(x*)` of surviving a step there (tolerance). Because `s` is steep near `x*`, small changes
in `b11` or `b00` move `s(x*)` a lot. Returns `x` and `survival`; `x = Inf` and survival
`s(∞)` when `b11 ≥ 1` (loads grow without bound).
"""
function stasis_load(p::Wilber2017Parameters; temperature::Real)
    x = p.b11 < 1 ? (p.b01 + p.b21 * temperature) / (1 - p.b11) : Inf
    s = IPMKernel(p; temperature).survival
    return (x = x, survival = isfinite(x) ? s(x) : s(1e6))
end

# --- Sensitivity analysis (Fig. 3) ------------------------------------------------------------

# The 19 parameters perturbed by sensitivity_analysis.R, in its order, and the ones dropped
# for each transmission function.
const WILBER2017_SENSITIVITY = (:b00, :b10, :b01, :b21, :b11, :v01, :c01, :b03, :b13, :v03,
    :c03, :b02, :b22, :b12, :β0, :β1)

"""
    extinction_sensitivity(p::Wilber2017Parameters; nsims = 1000, σ = 0.3, years = 8,
                           mesh = IPMMesh(:wilber2017), rng = Random.default_rng(),
                           zoospore_removal = true) -> NamedTuple

The global sensitivity analysis of Fig. 3 (Methods; `sensitivity_analysis.R`, after Sobie
2009): in each of `nsims` runs every transmission, resistance and tolerance parameter is
multiplied by an independent lognormal factor with median 1 and `sdlog = σ`, and one
stochastic run of [`HybridIPMProblem`](@ref) records whether the population went extinct
within `years`. The parameters are `b00, b10` (tolerance: survival), `b01, b21, b11, v01,
c01, b03, b13, v03, c03, b02, b22, b12` (resistance: growth, initial load, loss) and `β0`,
`β1` (transmission; `β1` only for density or frequency dependence): 16 for `:zpool_dd`
and `:zpool_fd`, 15 for `:zpool`.

Returns `names`, the perturbed values `X` (`nsims × length(names)`) and `extinct`
(`Bool` per run). [`sensitivity_weights`](@ref) fits the logistic regression of Fig. 3.
"""
function extinction_sensitivity(p::Wilber2017Parameters; nsims::Union{Nothing,Integer} = nothing,
        σ::Union{Nothing,Real} = nothing, years::Union{Nothing,Real} = nothing,
        mesh::IPMMesh = IPMMesh(:wilber2017), rng::AbstractRNG = Random.default_rng(),
        zoospore_removal::Bool = true)
    q = parameter_set(:wilber2017)
    nsims = nsims === nothing ? round(Int, q.n_sensitivity) : nsims
    σ = σ === nothing ? q.σ_sensitivity : σ
    years = years === nothing ? q.years_sensitivity : years
    names = p.transmission in (:zpool, :const) ? WILBER2017_SENSITIVITY[1:end-1] :
            WILBER2017_SENSITIVITY
    X = Matrix{Float64}(undef, nsims, length(names))
    extinct = falses(nsims)
    for s in 1:nsims
        vals = [getfield(p, f) * exp(σ * randn(rng)) for f in names]
        X[s, :] .= vals
        ps = _with(p; (f => v for (f, v) in zip(names, vals))...)
        prob = HybridIPMProblem(ps; mesh, years, stochastic = true, zoospore_removal)
        extinct[s] = _extinction_steps(_HybridCache(prob), 1, rng)[1] <= prob.steps
    end
    return (names = names, X = X, extinct = extinct)
end

"""
    sensitivity_weights(X, y; C = 1.0) -> NamedTuple

The regularised logistic regression of Fig. 3: columns of `X` are z-scored (population
standard deviation, as scikit-learn's `StandardScaler`) and extinction `y` is regressed on
them by L2-penalised logistic regression with inverse penalty `C`, minimising
`wᵀw/2 + C Σ log(1 + exp(−yᵢ wᵀxᵢ))` with `yᵢ = ±1` and the intercept as an extra, penalised
feature (scikit-learn's default `liblinear` solver of 2017). Returns `weights` (one per
column; positive means the parameter raises the extinction probability) and `intercept`.
"""
function sensitivity_weights(X::AbstractMatrix, y::AbstractVector; C::Real = 1.0)
    n, k = size(X)
    length(y) == n || throw(DimensionMismatch("y must have one entry per row of X"))
    μ = vec(sum(X; dims = 1)) ./ n
    sd = sqrt.(vec(sum((X .- μ') .^ 2; dims = 1)) ./ n)
    Z = hcat((X .- μ') ./ ifelse.(sd .> 0, sd, 1.0)', ones(n))
    s = [b ? 1.0 : -1.0 for b in y]
    w = zeros(k + 1)
    for _ in 1:200
        m = s .* (Z * w)
        π_ = _logistic.(-m)                   # 1 − σ(m)
        g = w .- C .* (Z' * (s .* π_))
        h = π_ .* (1 .- π_)
        H = LinearAlgebra.I(k + 1) + C .* (Z' * (h .* Z))
        Δ = H \ g
        w .-= Δ
        maximum(abs, Δ) < 1e-12 && break
    end
    return (weights = w[1:k], intercept = w[end])
end
