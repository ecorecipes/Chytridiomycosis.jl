# Louca, Lampo & Doebeli (2014) Proc R Soc B 281:20132783, "Assessing host extinction risk
# following exposure to Batrachochytrium dendrobatidis". doi:10.1098/rspb.2013.2783
# Checked against the published PDF and the electronic supplementary material (ESM), which the
# first author's site bundles after the article (lab/papers/louca2014_si/). No code or data
# were deposited. Units: days; degree of infection (DOI) = zoospores on a host, read in the
# ESM (§5.1) as the qPCR load of a standardised swab (ZE swab⁻¹).
# Wiki: lab/wiki/models/bd-louca2014-load-distribution.md
#
# State: the numbers of adults P_A(t, n) and tadpoles P_T(t, n) with DOI n = 0, 1, 2, … and the
# free zoospore pool Z(t). Master equations (main eqs 2.1–2.2; ESM eqs 1–3):
#
#   dP_s(n)/dt = [β_s(t)Z + I_s,u(n−1)] P_s(n−1) + I_s,d(n+1) P_s(n+1)
#                − [β_s(t)Z + I_s,u(n) + I_s,d(n)] P_s(n)
#                − [μ_s(N_s) + g_s(n) (+ 1/τ_m for tadpoles)] P_s(n) + δ_n0 (recruitment)
#   recruitment into adults   (1/τ_m)   Σ σ_m(n′)  P_T(n′)
#   recruitment into tadpoles r_br(t)   Σ σ_br(n′) P_A(n′)
#   dZ/dt = Σ_s η_s N_s S_s − μ_Z (Z − Z_o)
#
# Functional responses (ESM §1.4): wrapped-Cauchy breeding r_br(t) (eq. 4), σ_br(n) =
# 1/(1 + n/T_br), σ_m(n) = 1/(1 + n/T_m), μ_s(N) = μ_s(0) + c_s N_s (eq. 5) with
# μ_s(0) = −ln σ_s / Y and c_s from the carrying capacities (eq. 9), g_s(n) =
# (F_s − T_s) ln2/(Y T_s) · n/(F_s − n) (eq. 6), I_u(n) = 2 max(0, μ_S + λ)[√(n+1) − 1],
# I_d(n) = μ_S n − 2 min(0, μ_S + λ)[√(n+1) − 1] (eq. 7), β_s(t) = β̄_s [1 + b cos(2π(t − φ)/Y)].
#
# Numerics as in ESM §5: the DOIs 0 … F_s are split into m = 100 bins of width w = F_s/m (bin k
# holds DOI k w) and every transition rate of the first three rows is divided by w. A host whose
# DOI would reach F_s dies. Tolerant tadpoles (T_T = F_T, so g_T ≡ 0) have their DOI capped at
# S_T,max = F_T − 1, i.e. no upward transition from the top bin. The authors integrated with
# an explicit second-order Runge–Kutta scheme (dt = 10⁻³–10⁻² d). We offer the method-of-lines
# ODEProblem (Tsit5 by default) and a fixed-step semi-implicit scheme (backward Euler in the
# DOI transitions and losses, density dependence and recruitment lagged one step, exact
# exponential update of Z), which stays positive and stable when Z_o is large.
#
# Stochastic model (ESM §5.3): the base mortality μ_s(0) is replaced by random disturbances at
# rate μ_s(0)/f_s, each killing Binomial(⌊N_s⌋, f_s) hosts of stage s, every DOI bin being
# scaled by (1 − k_s/N_s).
#
# Departures and interpretations (see the wiki page):
# - Table 1 writes F_T = ∞ for tolerant tadpoles; the ESM (§1.3, Table S1) uses T_T = F_T with
#   S_T,max = F_T − 1. We follow the ESM with F_T = 10⁴.
# - ESM §5.2 lists r twice among the log-uniform parameters; the second is taken as T_br.
# - Q_o's risk effect is computed on log10 Q_o over the samples where Q_o is finite.

const LOUCA_STAGES = (:adult, :tadpole)
const LOUCA_OUTCOMES = (:extinction, :clearance, :persistence)

"""
    LoucaParameters(; Y = 365, Zo = 0, μZ, μSA, μST, λA, λT, βA, βT, φ = 0, b = 0, ηA, ηT,
                    FA, FT, TA, TT, STmax, Tm, Tbr, r, τm, τbr, σA, σT, KA, KT,
                    fA = 0.2, fT = 0.2, m = 100)
    LoucaParameters(:louca2014; set = "", kwargs...)

Parameters of the Louca et al. (2014) model. Rates are per day, DOIs (loads) in zoospores per
host (ZE per standardised swab in the paper's reading). A suffix `A` is the adult value and
`T` the tadpole value.

- `Y`: days in a year; `Zo`: environmental reservoir (zoospores);
- `μZ`: free zoospore loss rate; `μSA`, `μST`: asymptotic on-host zoospore loss rates;
- `λA`, `λT`: early within-host growth rates `I′(0)` (negative: the infection decays);
- `βA`, `βT`: mean pool-to-host transmission rates (zoospore⁻¹ day⁻¹); `b`, `φ`: relative
  amplitude and phase lag (days) of their annual oscillation;
- `ηA`, `ηT`: zoospores released per zoospore on a host per day;
- `FA`, `FT`: lethal DOIs; `TA`, `TT`: DOI tolerances (yearly survival halved);
  `STmax`: maximum tadpole DOI (`FT − 1` for tolerant tadpoles, which have `TT = FT`, or
  `FT`); `Tm`, `Tbr`: DOIs at which metamorphosis success and fecundity are halved;
- `r`: yearly fecundity; `τm`: expected tadpole lifetime; `τbr`: breeding period (full width
  at half maximum of the breeding function; `τbr ≥ Y/2` is constant breeding);
- `σA`, `σT`: yearly survivorships without disease or crowding; `KA`, `KT`: carrying
  capacities, which set the competition coefficients ([`louca_competition`](@ref));
- `fA`, `fT`: fraction killed by one disturbance (stochastic model);
- `m`: number of DOI bins.

The compendium form reads ESM Table S1 from `data/compendium/louca2014.toml`; `set` is one of
`"fig2ab"`, `"fig2cd"`, `"fig3a"` and `"fig3b"`. Keyword arguments override single values.
"""
Base.@kwdef struct LoucaParameters
    Y::Float64 = 365.0
    Zo::Float64 = 0.0
    μZ::Float64
    μSA::Float64
    μST::Float64
    λA::Float64
    λT::Float64
    βA::Float64
    βT::Float64
    φ::Float64 = 0.0
    b::Float64 = 0.0
    ηA::Float64
    ηT::Float64
    FA::Float64
    FT::Float64
    TA::Float64
    TT::Float64
    STmax::Float64
    Tm::Float64
    Tbr::Float64
    r::Float64
    τm::Float64
    τbr::Float64
    σA::Float64
    σT::Float64
    KA::Float64
    KT::Float64
    fA::Float64 = 0.2
    fT::Float64 = 0.2
    m::Int = 100
end

const _LOUCA_FIELDS = fieldnames(LoucaParameters)

function LoucaParameters(key::Symbol; set::AbstractString = "", kwargs...)
    q = parameter_set(key; set)
    vals = Dict{Symbol,Any}()
    for f in _LOUCA_FIELDS
        haskey(q, f) && (vals[f] = f === :m ? round(Int, q[f]) : q[f])
    end
    for (k, v) in kwargs
        vals[k] = v
    end
    return LoucaParameters(; vals...)
end

time_unit(::Type{<:LoucaParameters}) = :day

"""
    with(p::LoucaParameters; kwargs...) -> LoucaParameters

A copy of `p` with the given fields replaced.
"""
with(p::LoucaParameters; kw...) =
    LoucaParameters(; (f => get(kw, f, getfield(p, f)) for f in _LOUCA_FIELDS)...)

"""
    tadpoles_tolerate(p::LoucaParameters) -> Bool

Whether tadpoles tolerate infection: their DOI is capped below the lethal DOI
(`STmax < FT`), and they then carry no disease mortality (`TT = FT`; ESM §1.3).
"""
tadpoles_tolerate(p::LoucaParameters) = p.STmax < p.FT

function _check(p::LoucaParameters)
    p.m >= 2 || throw(ArgumentError("need at least 2 DOI bins"))
    (p.FA > 0 && p.FT > 0 && p.TA > 0 && p.TT > 0 && p.Tm > 0 && p.Tbr > 0) ||
        throw(ArgumentError("DOI thresholds must be positive"))
    (p.TA <= p.FA && p.TT <= p.FT) || throw(ArgumentError("need T_s ≤ F_s"))
    (p.STmax <= p.FT) || throw(ArgumentError("need STmax ≤ FT"))
    tadpoles_tolerate(p) && p.TT != p.FT &&
        throw(ArgumentError("tolerant tadpoles (STmax < FT) need TT = FT (ESM §1.3)"))
    (0 < p.σA < 1 && 0 < p.σT < 1) || throw(ArgumentError("survivorships must lie in (0, 1)"))
    (0 < p.fA <= 1 && 0 < p.fT <= 1) || throw(ArgumentError("fA, fT must lie in (0, 1]"))
    (p.μZ > 0 && p.μSA > 0 && p.μST > 0 && p.βA >= 0 && p.βT >= 0 && p.ηA >= 0 && p.ηT >= 0) ||
        throw(ArgumentError("rates must be non-negative (μZ, μS positive)"))
    (0 <= p.b <= 1) || throw(ArgumentError("b must lie in [0, 1]"))
    (p.r > 0 && p.τm > 0 && p.τbr > 0 && p.KA > 0 && p.KT > 0 && p.Zo >= 0) ||
        throw(ArgumentError("r, τm, τbr, KA, KT must be positive and Zo non-negative"))
    return p
end

# Per-stage view of the parameters.
function _louca_stage(p::LoucaParameters, s::Symbol)
    s === :adult && return (; μS = p.μSA, λ = p.λA, β = p.βA, η = p.ηA, F = p.FA, T = p.TA,
                            Smax = p.FA, σ = p.σA, K = p.KA, f = p.fA)
    s === :tadpole && return (; μS = p.μST, λ = p.λT, β = p.βT, η = p.ηT, F = p.FT, T = p.TT,
                              Smax = p.STmax, σ = p.σT, K = p.KT, f = p.fT)
    throw(ArgumentError("stage must be :adult or :tadpole"))
end

# --------------------------------------------------------------------------------------------
# Functional responses (ESM §1.4)

"""
    louca_breeding_shape(p::LoucaParameters) -> x

The shape `x = 2πγ/Y` of the wrapped-Cauchy breeding function (ESM eq. 4), set by the
breeding period `τbr` (its full width at half maximum, `r_br(τbr/2) = (max + min)/2`):
`cosh x = 1/cos(π τbr/Y)`. `Inf` (constant breeding) for `τbr ≥ Y/2`.
"""
function louca_breeding_shape(p::LoucaParameters)
    p.τbr >= p.Y / 2 && return Inf
    return acosh(1 / cos(π * p.τbr / p.Y))
end

"""
    louca_breeding_rate(p::LoucaParameters, t) -> r_br(t)

Tadpoles hatched per healthy adult per day (ESM eq. 4), a wrapped Cauchy density centred on
`t = 0` (mod `Y`) that integrates to `r` over a year:

`r_br(t) = (r/Y) sinh x / (cosh x − cos(2πt/Y))`, `x` from [`louca_breeding_shape`](@ref).
"""
function louca_breeding_rate(p::LoucaParameters, t::Real, x::Real = louca_breeding_shape(p))
    isinf(x) && return p.r / p.Y
    # sinh x/(cosh x − cos θ) without overflow for large x
    θ = 2π * t / p.Y
    return p.r / p.Y * sinh(x) / (cosh(x) - cos(θ))
end

"""
    louca_transmission(p::LoucaParameters, stage, t) -> β_s(t)

The pool-to-host transmission rate of `stage` (`:adult` or `:tadpole`) at time `t` (days):
`β̄_s [1 + b cos(2π(t − φ)/Y)]`.
"""
louca_transmission(p::LoucaParameters, s::Symbol, t::Real) =
    _louca_stage(p, s).β * (1 + p.b * cos(2π * (t - p.φ) / p.Y))

"""
    louca_fecundity_factor(p::LoucaParameters, n) -> σ_br(n)

Relative fecundity of an adult with DOI `n`: `1/(1 + n/T_br)`.
"""
louca_fecundity_factor(p::LoucaParameters, n::Real) = 1 / (1 + n / p.Tbr)

"""
    louca_metamorphosis_success(p::LoucaParameters, n) -> σ_m(n)

Probability that a tadpole with DOI `n` survives metamorphosis: `1/(1 + n/T_m)`.
"""
louca_metamorphosis_success(p::LoucaParameters, n::Real) = 1 / (1 + n / p.Tm)

"""
    louca_base_mortality(p::LoucaParameters, stage) -> μ_s(0)

Per-capita mortality (day⁻¹) of healthy hosts without crowding: `−ln σ_s / Y`.
"""
louca_base_mortality(p::LoucaParameters, s::Symbol) = -log(_louca_stage(p, s).σ) / p.Y

"""
    louca_disease_mortality(p::LoucaParameters, stage, n) -> g_s(n)

Bd-induced mortality (day⁻¹) at DOI `n` (ESM eq. 6):
`(F_s − T_s) ln2/(Y T_s) · n/(F_s − n)`, so that yearly survival is halved at `n = T_s`
(`exp(−Y g_s(T_s)) = 1/2`) and death is certain at `n = F_s`. Zero for tolerant tadpoles
(`T_T = F_T`).
"""
function louca_disease_mortality(p::LoucaParameters, s::Symbol, n::Real)
    q = _louca_stage(p, s)
    q.T >= q.F && return 0.0
    n >= q.F && return Inf
    return (q.F - q.T) * log(2) / (p.Y * q.T) * n / (q.F - n)
end

"""
    louca_load_rates(p::LoucaParameters, stage, n) -> (up, down)

Within-host transition rates (ESM eq. 7) from DOI `n` to `n + 1` (`up`, `I_u`) and to `n − 1`
(`down`, `I_d`), with `a = μ_S + λ`:

`I_u(n) = 2 max(0, a)(√(n+1) − 1)`, `I_d(n) = μ_S n − 2 min(0, a)(√(n+1) − 1)`.

Colonial growth on the perimeter of a localised infection gives the `√n` terms; `μ_S n` is the
intrinsic zoospore loss. `I′(0) = λ`.
"""
function louca_load_rates(p::LoucaParameters, s::Symbol, n::Real)
    q = _louca_stage(p, s)
    a = q.μS + q.λ
    sq = sqrt(n + 1) - 1
    return (2 * max(0.0, a) * sq, q.μS * n - 2 * min(0.0, a) * sq)
end

"""
    louca_load_speed(p::LoucaParameters, stage, n) -> I_s(n)

The expected DOI growth speed `I_u(n) − I_d(n)` of an individual with DOI `n` without
reinfection (ESM Fig. S1b).
"""
louca_load_speed(p::LoucaParameters, s::Symbol, n::Real) =
    (r = louca_load_rates(p, s, n); r[1] - r[2])

"""
    louca_within_host_capacity(λ, μS) -> K_S
    louca_within_host_capacity(p::LoucaParameters, stage) -> K_S

The DOI at which `I(n) = 0` for `λ > 0` (ESM §1.4): `K_S = 4λ(λ + μ_S)/μ_S²`. Infections
grow below it and decay above it; if `K_S ≥ F_s` they run to the lethal DOI. `0` for `λ ≤ 0`.
"""
louca_within_host_capacity(λ::Real, μS::Real) = λ > 0 ? 4λ * (λ + μS) / μS^2 : 0.0
louca_within_host_capacity(p::LoucaParameters, s::Symbol) =
    (q = _louca_stage(p, s); louca_within_host_capacity(q.λ, q.μS))

"""
    louca_within_host_loss(λ, KS) -> μ_S

The asymptotic on-host loss rate that gives within-host capacity `KS` at early growth rate
`λ > 0`: the positive root of `KS μ² = 4λ(λ + μ)` (the inverse of
[`louca_within_host_capacity`](@ref)). ESM §4.1 adopts λ = 0.38 d⁻¹ and K_S = 14.9 from a
logistic fit to *D. meridensis* and obtains μ_S ≈ 0.25 d⁻¹.
"""
function louca_within_host_loss(λ::Real, KS::Real)
    (λ > 0 && KS > 0) || throw(ArgumentError("need λ > 0 and KS > 0"))
    return (4λ + sqrt(16λ^2 + 16KS * λ^2)) / (2KS)
end

"""
    louca_competition(p::LoucaParameters) -> (; cA, cT)

The competition coefficients of the density-dependent mortality `μ_s(N) = μ_s(0) + c_s N_s`
that put the disease-free equilibrium with constant breeding at `(K_A, K_T)` (ESM eq. 9):

`c_A = [σ_m(0) K_T/(τ_m K_A) − μ_A(0)]/K_A`,  `c_T = [(r/Y) σ_br(0) K_A/K_T − 1/τ_m − μ_T(0)]/K_T`.

The parameter set is valid only if both are positive (ESM §5.2).
"""
function louca_competition(p::LoucaParameters)
    μA0, μT0 = louca_base_mortality(p, :adult), louca_base_mortality(p, :tadpole)
    cA = (louca_metamorphosis_success(p, 0) * p.KT / (p.τm * p.KA) - μA0) / p.KA
    cT = (p.r / p.Y * louca_fecundity_factor(p, 0) * p.KA / p.KT - 1 / p.τm - μT0) / p.KT
    return (; cA, cT)
end

"""
    louca_turnover(p::LoucaParameters) -> (; RA, RT)

Per-capita recruitment rates of the disease-free community, approximated as in ESM §5.2 from
the carrying capacities: `R̄_A ≈ σ_m(0) K_T/(K_A τ_m)` and `R̄_T ≈ (r/Y) K_A/K_T`.
"""
louca_turnover(p::LoucaParameters) =
    (; RA = louca_metamorphosis_success(p, 0) * p.KT / (p.KA * p.τm), RT = p.r / p.Y * p.KA / p.KT)

# --------------------------------------------------------------------------------------------
# Early invasion (main §3a; ESM §2.3)

"""
    reproduction_quotient(p::LoucaParameters) -> Q_o
    reproduction_quotient(; μZ, η, β, N, R, λ) -> Q_o

The basic reproduction quotient of Bd (main eq. 3.1; ESM eq. 19): the zoospores released into
the water that descend from one free zoospore early in the invasion,

`Q_o = (1/μ̄_Z) Σ_s η_s β̄_s N̄_s/(R̄_s − λ_s)`.

For a `LoucaParameters`, `N̄_s ≈ K_s` and `R̄_s` from [`louca_turnover`](@ref) (ESM §5.2).
The keyword form takes vectors over host types (several species and stages, ESM §2). `Inf` if
`λ_s ≥ R̄_s` for some type: within-host growth outpaces turnover and Bd invades whatever
the pool does.
"""
function reproduction_quotient(; μZ::Real, η, β, N, R, λ)
    any(R .- λ .<= 0) && return Inf
    return sum(η .* β .* N ./ (R .- λ)) / μZ
end
function reproduction_quotient(p::LoucaParameters)
    R = louca_turnover(p)
    return reproduction_quotient(; μZ = p.μZ, η = (p.ηA, p.ηT), β = (p.βA, p.βT),
                                 N = (p.KA, p.KT), R = (R.RA, R.RT), λ = (p.λA, p.λT))
end

"""
    louca_jacobian(p::LoucaParameters) -> Matrix

The time-averaged Jacobian `J̄_o` of the mean-field model at the disease-free state (ESM
eq. 18), in the variables (mean DOI of adults, of tadpoles, free zoospores):

```
[λ_A − R̄_A   0           β̄_A
 0           λ_T − R̄_T   β̄_T
 η_A N̄_A     η_T N̄_T     −μ̄_Z]
```

with `N̄_s ≈ K_s` and `R̄_s` from [`louca_turnover`](@ref). Its spectral bound
([`louca_growth_rate`](@ref)) is positive exactly when `Q_o > 1`.
"""
function louca_jacobian(p::LoucaParameters)
    R = louca_turnover(p)
    return [p.λA-R.RA 0.0 p.βA
            0.0 p.λT-R.RT p.βT
            p.ηA*p.KA p.ηT*p.KT -p.μZ]
end

"""
    louca_growth_rate(p::LoucaParameters; method = :eigen) -> ρ

The early exponential growth rate of Bd after exposure (ESM §2.3):

- `:eigen`: the spectral bound of [`louca_jacobian`](@ref);
- `:leading`: `ρ_o`, the positive root of `(ρ + μ_Z)ρ = Σ η_s β̄_s N̄_s` (ESM eq. 21), which
  ignores within-host dynamics and turnover;
- `:asymptotic`: the two-term approximation of ESM eq. (23),
  `ρ_o + Σ η_s β̄_s (λ_s − R̄_s) N̄_s / (ρ_o² + Σ η_s β̄_s N̄_s)`, for `|λ_s − R̄_s| ≪ ρ`.
"""
function louca_growth_rate(p::LoucaParameters; method::Symbol = :eigen)
    method === :eigen && return maximum(real, eigvals(louca_jacobian(p)))
    R = louca_turnover(p)
    A = p.ηA * p.βA * p.KA + p.ηT * p.βT * p.KT
    ρo = (-p.μZ + sqrt(p.μZ^2 + 4A)) / 2
    method === :leading && return ρo
    method === :asymptotic || throw(ArgumentError("method must be :eigen, :leading or :asymptotic"))
    B = p.ηA * p.βA * (p.λA - R.RA) * p.KA + p.ηT * p.βT * (p.λT - R.RT) * p.KT
    return ρo + B / (ρo^2 + A)
end

# --------------------------------------------------------------------------------------------
# Parameter estimates (ESM §4)

"""
    louca_tolerance_from_slope(g′0; Y = 365) -> T

The DOI tolerance implied by a mortality that rises linearly in DOI with slope `g′(0)`
(day⁻¹ per zoospore) when `F ≫ T` in ESM eq. (6): `T = ln 2/(Y g′(0))`. ESM §4.1 fits
`g′(0) = 0.0011` to *D. meridensis* survival and obtains `T ≈ 1.73`.
"""
louca_tolerance_from_slope(g′0::Real; Y::Real = 365) = log(2) / (Y * g′0)

"""
    louca_survival_curve(T; xo, KS, λ, g′0, μo = 0) -> σ(T)

Survival to `T` days after exposure of hosts whose mean DOI grows logistically from `xo` with
rate `λ` to capacity `KS`, under a mortality `g′0 · DOI` plus base mortality `μo` (ESM eq. 31):

`σ(T) = [1 + (xo/KS)(e^{λT} − 1)]^{−KS g′0/λ} e^{−μo T}`.
"""
louca_survival_curve(T::Real; xo::Real, KS::Real, λ::Real, g′0::Real, μo::Real = 0.0) =
    (1 + xo / KS * (exp(λ * T) - 1))^(-KS * g′0 / λ) * exp(-μo * T)

"""
    louca_transmission_from_growth(ρ; μZ, η, N) -> β̃

The mean transmission rate implied by an epizootic growth rate `ρ` when within-host dynamics
are ignored (ESM §4.2, from eq. 21): `β̃ = (ρ + μ_Z)ρ/(η N)`.
"""
louca_transmission_from_growth(ρ::Real; μZ::Real, η::Real, N::Real) = (ρ + μZ) * ρ / (η * N)

"""
    louca_endemic_transmission(; μZ, R, η, N, λ = 0) -> β

The transmission rate implied by a single host type in a long-term endemic state (ESM eq. 25,
`η β N/(μ_Z (R − λ)) ≈ 1`): `β = μ_Z (R − λ)/(η N)`, with `R` the per-capita recruitment
rate (day⁻¹). ESM §4.2 applies it to *Atelopus cruciger*.
"""
louca_endemic_transmission(; μZ::Real, R::Real, η::Real, N::Real, λ::Real = 0.0) =
    μZ * (R - λ) / (η * N)

# --------------------------------------------------------------------------------------------
# DOI-distribution model: binned master equations

# Bin-level rates of one stage, independent of time and state.
struct _LoucaBins
    w::Float64              # bin width (DOI per bin)
    n::Vector{Float64}      # DOI of each bin, n_k = (k − 1) w
    up::Vector{Float64}     # I_u(n_k)/w
    down::Vector{Float64}   # I_d(n_k)/w
    g::Vector{Float64}      # g_s(n_k)
    rec::Vector{Float64}    # σ_br(n_k) (adults) or σ_m(n_k) (tadpoles)
    capped::Bool            # tolerant tadpoles: no upward move out of the top bin
end

function _louca_bins(p::LoucaParameters, s::Symbol)
    q = _louca_stage(p, s)
    w = q.F / p.m
    n = [(k - 1) * w for k in 1:p.m]
    rates = [louca_load_rates(p, s, x) for x in n]
    up = first.(rates) ./ w
    down = last.(rates) ./ w
    g = [louca_disease_mortality(p, s, x) for x in n]
    rec = s === :adult ? [louca_fecundity_factor(p, x) for x in n] :
          [louca_metamorphosis_success(p, x) for x in n]
    capped = s === :tadpole && tadpoles_tolerate(p)
    capped && (up[end] = 0.0)
    return _LoucaBins(w, n, up, down, g, rec, capped)
end

struct _LoucaCache
    p::LoucaParameters
    A::_LoucaBins
    T::_LoucaBins
    x::Float64      # breeding shape
    μA0::Float64
    μT0::Float64
    cA::Float64
    cT::Float64
    base::Bool      # deterministic base mortality on (false in the stochastic model)
    status::Base.RefValue{Symbol}   # outcome, set by the end-of-run callback
    t_reservoir::Float64            # the reservoir supply μ_Z Z_o starts at this time
end

function _louca_cache(p::LoucaParameters; stochastic::Bool = false, t_reservoir::Real = -Inf)
    c = louca_competition(p)
    return _LoucaCache(p, _louca_bins(p, :adult), _louca_bins(p, :tadpole),
        louca_breeding_shape(p), louca_base_mortality(p, :adult),
        louca_base_mortality(p, :tadpole), c.cA, c.cT, !stochastic, Ref(:persistence),
        Float64(t_reservoir))
end

# u = [P_A (m); P_T (m); Z]
@inline _louca_views(u, m) = (view(u, 1:m), view(u, (m + 1):(2m)), 2m + 1)

# du for one stage given the infection rate per bin ι, extra loss and recruitment into bin 1.
@inline function _louca_stage_rhs!(du, P, B::_LoucaBins, ι, loss, rec_in)
    m = length(P)
    @inbounds for k in 1:m
        upk = (k == m && B.capped) ? 0.0 : ι + B.up[k]
        out = upk + B.down[k] + loss + B.g[k]
        v = -out * P[k]
        k > 1 && (v += (ι + B.up[k-1]) * P[k-1])
        k < m && (v += B.down[k+1] * P[k+1])
        du[k] = v
    end
    @inbounds du[1] += rec_in
    return du
end

function _louca_rhs!(du, u, c::_LoucaCache, t)
    p = c.p
    m = p.m
    PA, PT, iz = _louca_views(u, m)
    dA, dT, _ = _louca_views(du, m)
    Z = max(u[iz], 0.0)
    NA, NT = sum(PA), sum(PT)
    season = 1 + p.b * cos(2π * (t - p.φ) / p.Y)
    ιA = p.βA * season * Z / c.A.w
    ιT = p.βT * season * Z / c.T.w
    lossA = (c.base ? c.μA0 : 0.0) + c.cA * NA
    lossT = (c.base ? c.μT0 : 0.0) + c.cT * NT + 1 / p.τm
    recA = sum(c.T.rec[k] * PT[k] for k in 1:m) / p.τm
    recT = louca_breeding_rate(p, t, c.x) * sum(c.A.rec[k] * PA[k] for k in 1:m)
    _louca_stage_rhs!(dA, PA, c.A, ιA, lossA, recA)
    _louca_stage_rhs!(dT, PT, c.T, ιT, lossT, recT)
    shed = p.ηA * sum(c.A.n[k] * PA[k] for k in 1:m) + p.ηT * sum(c.T.n[k] * PT[k] for k in 1:m)
    Zo = t >= c.t_reservoir ? p.Zo : 0.0
    du[iz] = shed - p.μZ * (u[iz] - Zo)
    return nothing
end

"""
    LoucaDOIProblem(p::LoucaParameters; years = 10, exposure_time = 1.5 Y,
                    exposure_doi = 1000, stochastic = false, saveat = 1.0,
                    save_distributions = false, u0 = nothing,
                    reservoir_onset = :exposure)

The DOI-distribution model of Louca et al. (2014) (main eqs 2.1–2.2, ESM eqs 1–3), binned as
in ESM §5: `p.m` bins of width `F_s/m` per stage, transition rates divided by the bin width, a
host dies when its DOI would reach `F_s`, tolerant tadpoles stay in the top bin. The state is
`[P_A; P_T; Z]` (hosts per bin, then free zoospores).

The run starts from a healthy population at `(K_A, K_T)` (or `u0`) and no free zoospores;
at `exposure_time` one adult with DOI `exposure_doi` (rounded to a bin) is added. With
`reservoir_onset = :exposure` (default) the reservoir supply `μ_Z Z_o` also starts then, so the
population is first exposed in the second year (§2b); `:start` has the reservoir from `t = 0`
with `Z(0) = Z_o`, which infects the population during the burn-in. The paper does not say
which it used. It ends at `years` years, when
tadpoles plus adults fall below 2 (host extinction) or, after exposure, when free plus
within-host zoospores fall below 10 (Bd clearance). Horizon, exposure and thresholds come
from the compendium (ESM §5).

With `stochastic = true`, base mortality is replaced by random disturbances (ESM §5.3): at
rate `μ_s(0)/f_s`, `Binomial(⌊N_s⌋, f_s)` hosts of stage `s` die, every bin scaled alike.

`solve(prob; method = :semi_implicit, dt = 0.05, rng)` integrates with a fixed-step scheme
(backward Euler in the DOI transitions and losses, density dependence and recruitment from the
previous step, exact update of `Z`), which stays positive when `Z` is large; disturbances
happen in a step with probability `dt μ_s(0)/f_s`, as in the authors' code. `method = :ode`
solves [`to_ode_problem`](@ref) with `alg` (default `Tsit5()`); the stochastic version then
draws Poisson disturbance times in advance. Both return a named tuple:

- `outcome` (`:extinction`, `:clearance` or `:persistence`) and `t_end`;
- `t`, `NA`, `NT` (hosts), `SA`, `ST` (mean DOI), `prevalence_A`, `prevalence_T` (fraction
  with DOI > 0), `Z` at the save times;
- `doi`: the DOI of each bin per stage, `PA`, `PT`: the final distributions, and with
  `save_distributions = true`, `PA_t`, `PT_t`: `m × length(t)` matrices.
"""
struct LoucaDOIProblem
    p::LoucaParameters
    tspan::Tuple{Float64,Float64}
    exposure_time::Float64
    exposure_doi::Float64
    stochastic::Bool
    saveat::Float64
    save_distributions::Bool
    extinction_threshold::Float64
    clearance_threshold::Float64
    u0::Vector{Float64}
    reservoir_onset::Symbol
end

const _LOUCA_DESIGN = Ref{Any}(nothing)
function _louca_design()
    if _LOUCA_DESIGN[] === nothing
        q = parameter_set(:louca2014)
        _LOUCA_DESIGN[] = (; Y = q.Y, years = q.years, exposure_time = q.exposure_time,
            exposure_doi = q.exposure_doi, extinction = q.extinction_threshold,
            clearance = q.clearance_threshold, n_samples = round(Int, q.n_samples))
    end
    return _LOUCA_DESIGN[]::NamedTuple
end

function LoucaDOIProblem(p::LoucaParameters; years::Real = _louca_design().years,
        exposure_time::Real = _louca_design().exposure_time * p.Y,
        exposure_doi::Real = _louca_design().exposure_doi, stochastic::Bool = false,
        saveat::Real = 1.0, save_distributions::Bool = false, u0 = nothing,
        reservoir_onset::Symbol = :exposure)
    _check(p)
    reservoir_onset in (:exposure, :start) ||
        throw(ArgumentError("reservoir_onset must be :exposure or :start"))
    c = louca_competition(p)
    (c.cA > 0 && c.cT > 0) ||
        throw(ArgumentError("competition coefficients must be positive (ESM eq. 9): $c"))
    years > 0 && saveat > 0 || throw(ArgumentError("years and saveat must be positive"))
    exposure_doi < p.FA || throw(ArgumentError("exposure_doi must be below FA"))
    m = p.m
    if u0 === nothing
        u = zeros(2m + 1)
        u[1] = p.KA
        u[m+1] = p.KT
        u[end] = reservoir_onset === :start ? p.Zo : 0.0
    else
        length(u0) == 2m + 1 || throw(ArgumentError("u0 must have length 2m + 1"))
        u = Float64.(collect(u0))
    end
    d = _louca_design()
    return LoucaDOIProblem(p, (0.0, years * p.Y), Float64(exposure_time),
        Float64(exposure_doi), stochastic, Float64(saveat), save_distributions,
        d.extinction, d.clearance, u, reservoir_onset)
end

_louca_t_reservoir(prob::LoucaDOIProblem) =
    prob.reservoir_onset === :start ? -Inf : prob.exposure_time

_louca_exposure_bin(prob::LoucaDOIProblem) =
    clamp(round(Int, prob.exposure_doi / (prob.p.FA / prob.p.m)) + 1, 1, prob.p.m)

function _louca_total_zoospores(u, c::_LoucaCache)
    m = c.p.m
    PA, PT, iz = _louca_views(u, m)
    return max(u[iz], 0.0) + sum(c.A.n[k] * PA[k] for k in 1:m) + sum(c.T.n[k] * PT[k] for k in 1:m)
end

function _louca_status(u, c::_LoucaCache, prob::LoucaDOIProblem, t)
    m = c.p.m
    sum(view(u, 1:2m)) < prob.extinction_threshold && return :extinction
    t > prob.exposure_time && _louca_total_zoospores(u, c) < prob.clearance_threshold &&
        return :clearance
    return :persistence
end

# A disturbance of stage s (ESM §5.3).
function _louca_disturb!(u, s::Symbol, p::LoucaParameters, rng::AbstractRNG)
    m = p.m
    idx = s === :adult ? (1:m) : ((m + 1):(2m))
    N = sum(view(u, idx))
    N <= 0 && return u
    f = _louca_stage(p, s).f
    k = rand(rng, Binomial(floor(Int, N), f))
    scale = max(0.0, 1 - k / N)
    @inbounds for i in idx
        u[i] *= scale
    end
    return u
end

"""
    to_ode_problem(prob::LoucaDOIProblem; rng = Random.default_rng()) -> ODEProblem

The method-of-lines `ODEProblem` of the binned master equations, with callbacks for the
exposure, the end of the run (extinction or clearance) and, for the stochastic model, the
disturbances at Poisson times drawn with `rng`. The callbacks keep state (the next
disturbance, the outcome in `ode.p.status`), so lower afresh for each solve.
"""
function to_ode_problem(prob::LoucaDOIProblem; rng::AbstractRNG = Random.default_rng())
    p = prob.p
    c = _louca_cache(p; stochastic = prob.stochastic, t_reservoir = _louca_t_reservoir(prob))
    m = p.m
    kexp = _louca_exposure_bin(prob)
    texp = prob.exposure_time
    status = c.status
    # disturbance times, by stage
    events = Tuple{Float64,Symbol}[]
    if prob.stochastic
        for s in LOUCA_STAGES
            rate = louca_base_mortality(p, s) / _louca_stage(p, s).f
            t = prob.tspan[1]
            while true
                t += Random.randexp(rng) / rate
                t > prob.tspan[2] && break
                push!(events, (t, s))
            end
        end
        sort!(events; by = first)
    end
    next = Ref(1)
    tstops = sort!(vcat([texp], first.(events)))
    cb_exp = DiscreteCallback((u, t, i) -> t == texp, i -> (i.u[kexp] += 1; nothing);
        save_positions = (false, false))
    cb_dist = DiscreteCallback((u, t, i) -> next[] <= length(events) && t == events[next[]][1],
        function (i)
            while next[] <= length(events) && events[next[]][1] == i.t
                _louca_disturb!(i.u, events[next[]][2], p, rng)
                next[] += 1
            end
            return nothing
        end; save_positions = (false, false))
    cb_end = DiscreteCallback((u, t, i) -> (status[] = _louca_status(u, c, prob, t)) !== :persistence,
        i -> SciMLBase.terminate!(i); save_positions = (false, true))
    cbs = SciMLBase.CallbackSet(cb_exp, cb_dist, cb_end)
    return ODEProblem(_louca_rhs!, copy(prob.u0), prob.tspan, c; callback = cbs, tstops)
end

function _louca_summary(c::_LoucaCache, u)
    m = c.p.m
    PA, PT, iz = _louca_views(u, m)
    NA, NT = sum(PA), sum(PT)
    SA = NA > 0 ? sum(c.A.n[k] * PA[k] for k in 1:m) / NA : 0.0
    ST = NT > 0 ? sum(c.T.n[k] * PT[k] for k in 1:m) / NT : 0.0
    prevA = NA > 0 ? 1 - PA[1] / NA : 0.0
    prevT = NT > 0 ? 1 - PT[1] / NT : 0.0
    return (NA, NT, SA, ST, prevA, prevT, u[iz])
end

function _louca_output(prob, c, ts, us, outcome, t_end)
    m = prob.p.m
    s = [_louca_summary(c, u) for u in us]
    out = (; outcome, t_end, t = ts, NA = getindex.(s, 1), NT = getindex.(s, 2),
           SA = getindex.(s, 3), ST = getindex.(s, 4), prevalence_A = getindex.(s, 5),
           prevalence_T = getindex.(s, 6), Z = getindex.(s, 7),
           doi = (; adult = copy(c.A.n), tadpole = copy(c.T.n)),
           PA = isempty(us) ? Float64[] : us[end][1:m],
           PT = isempty(us) ? Float64[] : us[end][(m+1):(2m)])
    if prob.save_distributions
        out = merge(out, (; PA_t = reduce(hcat, [u[1:m] for u in us]),
                          PT_t = reduce(hcat, [u[(m+1):(2m)] for u in us])))
    end
    return out
end

# Thomas algorithm for one stage: (1 + dt out_k) x_k − dt(ι + up_{k−1}) x_{k−1}
# − dt down_{k+1} x_{k+1} = rhs_k, overwriting P with x. cp is scratch.
@inline function _louca_implicit_stage!(P, B::_LoucaBins, ι, loss, rec_in, dt, cp)
    m = length(P)
    @inbounds P[1] += dt * rec_in
    # forward sweep
    @inbounds begin
        upk = (m == 1 && B.capped) ? 0.0 : ι + B.up[1]
        bk = 1 + dt * (upk + B.down[1] + loss + B.g[1])
        ck = m > 1 ? -dt * B.down[2] : 0.0
        cp[1] = ck / bk
        P[1] = P[1] / bk
        for k in 2:m
            upk = (k == m && B.capped) ? 0.0 : ι + B.up[k]
            ak = -dt * (ι + B.up[k-1])
            bk = 1 + dt * (upk + B.down[k] + loss + B.g[k])
            ck = k < m ? -dt * B.down[k+1] : 0.0
            den = bk - ak * cp[k-1]
            cp[k] = ck / den
            P[k] = (P[k] - ak * P[k-1]) / den
        end
        for k in (m-1):-1:1
            P[k] -= cp[k] * P[k+1]
        end
        for k in 1:m
            P[k] = max(P[k], 0.0)
        end
    end
    return P
end

function _louca_semi_implicit(prob::LoucaDOIProblem, dt::Real, rng::AbstractRNG)
    p = prob.p
    c = _louca_cache(p; stochastic = prob.stochastic, t_reservoir = _louca_t_reservoir(prob))
    m = p.m
    u = copy(prob.u0)
    PA, PT, iz = _louca_views(u, m)
    cp = zeros(m)
    t0, tf = prob.tspan
    nsteps = ceil(Int, (tf - t0) / dt)
    dt = (tf - t0) / nsteps
    kexp = _louca_exposure_bin(prob)
    exposed = prob.exposure_time < t0
    ratesA = louca_base_mortality(p, :adult) / p.fA
    ratesT = louca_base_mortality(p, :tadpole) / p.fT
    ts = Float64[t0]
    us = [copy(u)]
    nextsave = t0 + prob.saveat
    outcome = :persistence
    t = t0
    eμ = exp(-p.μZ * dt)
    for i in 1:nsteps
        if !exposed && t >= prob.exposure_time - 1e-9
            u[kexp] += 1
            exposed = true
        end
        NA, NT = sum(PA), sum(PT)
        Z = max(u[iz], 0.0)
        season = 1 + p.b * cos(2π * (t - p.φ) / p.Y)
        ιA = p.βA * season * Z / c.A.w
        ιT = p.βT * season * Z / c.T.w
        lossA = (c.base ? c.μA0 : 0.0) + c.cA * NA
        lossT = (c.base ? c.μT0 : 0.0) + c.cT * NT + 1 / p.τm
        recA = sum(c.T.rec[k] * PT[k] for k in 1:m) / p.τm
        recT = louca_breeding_rate(p, t, c.x) * sum(c.A.rec[k] * PA[k] for k in 1:m)
        _louca_implicit_stage!(PA, c.A, ιA, lossA, recA, dt, cp)
        _louca_implicit_stage!(PT, c.T, ιT, lossT, recT, dt, cp)
        shed = p.ηA * sum(c.A.n[k] * PA[k] for k in 1:m) +
               p.ηT * sum(c.T.n[k] * PT[k] for k in 1:m)
        Zo = t >= c.t_reservoir ? p.Zo : 0.0
        u[iz] = Zo + (Z - Zo) * eμ + shed * (1 - eμ) / p.μZ
        if prob.stochastic
            rand(rng) < dt * ratesA && _louca_disturb!(u, :adult, p, rng)
            rand(rng) < dt * ratesT && _louca_disturb!(u, :tadpole, p, rng)
        end
        t = t0 + i * dt
        if t >= nextsave - 1e-9
            push!(ts, t)
            push!(us, copy(u))
            nextsave += prob.saveat
        end
        outcome = _louca_status(u, c, prob, t)
        if outcome !== :persistence
            ts[end] == t || (push!(ts, t); push!(us, copy(u)))
            break
        end
    end
    return _louca_output(prob, c, ts, us, outcome, t)
end

function CommonSolve.solve(prob::LoucaDOIProblem; method::Symbol = :semi_implicit,
        dt::Real = 0.05, rng::AbstractRNG = Random.default_rng(), alg = Tsit5(),
        reltol = 1e-6, abstol = 1e-6, kwargs...)
    method === :semi_implicit && return _louca_semi_implicit(prob, dt, rng)
    method === :ode || throw(ArgumentError("method must be :semi_implicit or :ode"))
    ode = to_ode_problem(prob; rng)
    sol = solve(ode, alg; saveat = prob.saveat, reltol, abstol, kwargs...)
    c = ode.p
    status = c.status[]
    return _louca_output(prob, c, sol.t, sol.u, status, sol.t[end])
end

"""
    louca_outcome(sol) -> Symbol

The outcome of a solved [`LoucaDOIProblem`](@ref): `:extinction` (hosts below 2),
`:clearance` (free and within-host zoospores below 10 after exposure) or `:persistence`.
"""
louca_outcome(sol::NamedTuple) = sol.outcome

# --------------------------------------------------------------------------------------------
# Monte Carlo assessment of extinction risk (§2b, §3c; ESM §5.2)

# Sampled parameters, with their scale; ranges from the compendium.
const _LOUCA_SAMPLED = (:Zo, :μZ, :μSA, :μST, :λA, :λT, :βA, :βT, :φ, :b, :ηA, :ηT, :TA, :TT,
    :Tm, :Tbr, :r, :τm, :τbr, :σA, :σT, :KA, :KT)
const _LOUCA_LOG = (:Zo, :μZ, :μSA, :μST, :βA, :βT, :TA, :TT, :Tm, :Tbr, :r, :KA, :KT)

"""
    louca_ranges(; tadpoles = :suffer, reservoir = false, log_params) -> NamedTuple

The sampling range `(lo, hi, log)` of each sampled parameter (Table 1; scales from ESM §5.2),
in days, for one of the paper's four scenarios: tadpoles that `:suffer` or `:tolerate`
infection, with (`reservoir = true`, `Z_o` from 10 to 10¹⁵) or without (`Z_o = 0`) an
environmental reservoir. Fixed parameters (`Zo` without a reservoir, `TT` for tolerant
tadpoles) are omitted. `log = true` marks log-uniform sampling, for the parameters in
`log_params`. ESM §5.2 lists `r, Z_o, μ_Z, μ_S,A, μ_S,T, K_A, K_T, β_A, β_T, T_A, T_T, T_m, r`;
the second `r` is read as `T_br` (the default).
"""
function louca_ranges(; tadpoles::Symbol = :suffer, reservoir::Bool = false,
        log_params = _LOUCA_LOG)
    tadpoles in (:suffer, :tolerate) || throw(ArgumentError("tadpoles must be :suffer or :tolerate"))
    q = parameter_set(:louca2014)
    Y = q.Y
    rg = Dict{Symbol,Tuple{Float64,Float64}}(
        :Zo => (q.Zo_min, q.Zo_max), :μZ => (q.μZ_min, q.μZ_max),
        :μSA => (q.μS_min, q.μS_max), :μST => (q.μS_min, q.μS_max),
        :λA => (q.λ_min, q.λ_max), :λT => (q.λ_min, q.λ_max),
        :βA => (q.β_min, q.β_max), :βT => (q.β_min, q.β_max),
        :φ => (q.φ_min * Y, q.φ_max * Y), :b => (q.b_min, q.b_max),
        :ηA => (q.η_min, q.η_max), :ηT => (q.η_min, q.η_max),
        :TA => (q.T_min, q.FA), :TT => (q.T_min, q.FT),
        :Tm => (q.T_min, q.Tmbr_max), :Tbr => (q.T_min, q.Tmbr_max),
        :r => (q.r_min, q.r_max), :τm => (q.τm_min, q.τm_max * Y),
        :τbr => (q.τbr_min, q.τbr_max * Y), :σA => (q.σA_min, q.σA_max),
        :σT => (q.σT_min, q.σT_max), :KA => (q.K_min, q.K_max), :KT => (q.K_min, q.K_max))
    keep = filter(s -> (reservoir || s !== :Zo) && (tadpoles === :suffer || s !== :TT),
                  collect(_LOUCA_SAMPLED))
    return NamedTuple{Tuple(keep)}(Tuple((lo = rg[s][1], hi = rg[s][2], log = s in log_params)
                                         for s in keep))
end

"""
    louca_sample(rng; tadpoles = :suffer, reservoir = false, maxtries = 10_000, log_params)
        -> (; p, u)

One parameter set drawn as in ESM §5.2: each parameter of [`louca_ranges`](@ref) uniform on
its (log or linear) scale, rejected and redrawn until both competition coefficients of
[`louca_competition`](@ref) are positive. The fixed values (`FA`, `FT`, `fA`, `fT`, `m`, `Y`)
come from the compendium. Tolerant tadpoles get `TT = FT` and `STmax = FT − 1`; suffering
ones `STmax = FT`. Returns the parameters `p` and the scaled draws `u` (each in [0, 1], on
the sampling scale).
"""
function louca_sample(rng::AbstractRNG; tadpoles::Symbol = :suffer, reservoir::Bool = false,
        maxtries::Integer = 10_000, log_params = _LOUCA_LOG)
    rg = louca_ranges(; tadpoles, reservoir, log_params)
    q = parameter_set(:louca2014)
    for _ in 1:maxtries
        u = map(_ -> rand(rng), rg)
        v = map(rg, u) do r, x
            r.log ? exp10(log10(r.lo) + x * (log10(r.hi) - log10(r.lo))) : r.lo + x * (r.hi - r.lo)
        end
        TT = tadpoles === :tolerate ? q.FT : v.TT
        STmax = tadpoles === :tolerate ? q.FT - 1 : q.FT
        p = LoucaParameters(; Y = q.Y, Zo = reservoir ? v.Zo : 0.0, μZ = v.μZ, μSA = v.μSA,
            μST = v.μST, λA = v.λA, λT = v.λT, βA = v.βA, βT = v.βT, φ = v.φ, b = v.b,
            ηA = v.ηA, ηT = v.ηT, FA = q.FA, FT = q.FT, TA = v.TA, TT, STmax, Tm = v.Tm,
            Tbr = v.Tbr, r = v.r, τm = v.τm, τbr = v.τbr, σA = v.σA, σT = v.σT, KA = v.KA,
            KT = v.KT, fA = q.fA, fT = q.fT, m = round(Int, q.m))
        c = louca_competition(p)
        (c.cA > 0 && c.cT > 0) && return (; p, u)
    end
    throw(ErrorException("no valid parameter set in $maxtries draws"))
end

"""
    louca_monte_carlo(n; tadpoles = :suffer, reservoir = false, rng, stochastic = false,
                      dt = 0.05, years = 10, reservoir_onset = :exposure) -> Vector{NamedTuple}

`n` runs of the paper's Monte Carlo assessment (§2b; ESM §5.2): each draws a parameter set with
[`louca_sample`](@ref), exposes a healthy population to one infected adult after 1.5 years and
runs the [`LoucaDOIProblem`](@ref) (semi-implicit, step `dt`) for up to `years` years. Each
element holds `p`, the scaled draws `u`, `Qo` ([`reproduction_quotient`](@ref)), `outcome`
and `t_end`. The paper used `n = 10⁴` per scenario.
"""
function louca_monte_carlo(n::Integer; tadpoles::Symbol = :suffer, reservoir::Bool = false,
        rng::AbstractRNG = Random.default_rng(), stochastic::Bool = false, dt::Real = 0.05,
        years::Real = _louca_design().years, reservoir_onset::Symbol = :exposure)
    return map(1:n) do _
        s = louca_sample(rng; tadpoles, reservoir)
        sol = solve(LoucaDOIProblem(s.p; years, stochastic, saveat = years * s.p.Y,
                                    reservoir_onset); dt, rng)
        (; p = s.p, u = s.u, Qo = reproduction_quotient(s.p), outcome = sol.outcome,
         t_end = sol.t_end)
    end
end

"""
    louca_risk_effects(results; outcome = :extinction) -> NamedTuple

The risk effect of every sampled parameter and of `Q_o` (§2b; ESM §5.2): the slope `A` of the
least-squares line of the outcome indicator (1 for `outcome`, else 0) on the parameter, times
the parameter's range `X_max − X_min`. Parameters are taken on their sampling scale (log10 for
log-uniform ones), so the risk effect is the slope on the scaled draw `u ∈ [0, 1]`. For `Q_o`
(not sampled) the regression is on `log10 Q_o` over the runs with finite `Q_o`, scaled by its
observed range.

`results` is the output of [`louca_monte_carlo`](@ref), or any vector of named tuples with
fields `u` (named tuple of scaled draws), `Qo` and `outcome`.
"""
function louca_risk_effects(results; outcome::Symbol = :extinction)
    isempty(results) && throw(ArgumentError("no results"))
    y = [r.outcome === outcome ? 1.0 : 0.0 for r in results]
    slope(x, y) = (xm = sum(x) / length(x); ym = sum(y) / length(y);
                   sum((x .- xm) .* (y .- ym)) / sum((x .- xm) .^ 2))
    names = keys(first(results).u)
    eff = map(names) do s
        x = [Float64(getfield(r.u, s)) for r in results]
        slope(x, y) * (maximum(x) - minimum(x))
    end
    fin = findall(r -> isfinite(r.Qo) && r.Qo > 0, results)
    xq = [log10(results[i].Qo) for i in fin]
    eQ = length(fin) > 2 ? slope(xq, y[fin]) * (maximum(xq) - minimum(xq)) : NaN
    return merge((; Qo = eQ), NamedTuple{names}(Tuple(eff)))
end

"""
    louca_outcome_probabilities(results; edges, variable = :Qo) -> (; centers, extinction,
                                clearance, persistence, n)

The fraction of runs ending in each outcome, binned by `log10` of `variable` (`:Qo`, or a
parameter of the named tuples `r.p`) at the given `edges` (in log10 units), as in Figs 4–5.
Runs with infinite `Q_o` are dropped. `n` is the number of runs per bin.
"""
function louca_outcome_probabilities(results; edges, variable::Symbol = :Qo)
    x = [variable === :Qo ? r.Qo : getfield(r.p, variable) for r in results]
    keep = findall(v -> isfinite(v) && v > 0, x)
    lx = log10.(x[keep])
    out = [results[i].outcome for i in keep]
    nb = length(edges) - 1
    cnt = zeros(Int, nb)
    acc = Dict(o => zeros(nb) for o in LOUCA_OUTCOMES)
    for (v, o) in zip(lx, out)
        j = searchsortedlast(edges, v)
        (1 <= j <= nb) || continue
        cnt[j] += 1
        acc[o][j] += 1
    end
    frac(o) = [cnt[j] > 0 ? acc[o][j] / cnt[j] : NaN for j in 1:nb]
    centers = [(edges[j] + edges[j+1]) / 2 for j in 1:nb]
    return (; centers, extinction = frac(:extinction), clearance = frac(:clearance),
            persistence = frac(:persistence), n = cnt)
end

# --------------------------------------------------------------------------------------------
# Stored reproduction output

"""
    louca2014_path()

Path of the stored output of our Monte Carlo reproduction of Louca et al. (2014)
(`data/louca2014/reproduction.toml`, written by `scripts/louca2014_montecarlo.jl`). These are
our simulations, not published data; see `data/louca2014/README.md`.
"""
louca2014_path() = joinpath(pkgdir(@__MODULE__), "data", "louca2014", "reproduction.toml")

"""
    load_louca2014() -> Dict

The stored Monte Carlo summary ([`louca2014_path`](@ref)): per scenario (`"suffer_noreservoir"`,
`"tolerate_noreservoir"`, `"suffer_reservoir"`, `"tolerate_reservoir"`) the number of runs,
the fraction of each outcome, the risk effects and the outcome probabilities binned by
`log10 Q_o`.
"""
load_louca2014() = TOML.parsefile(louca2014_path())
