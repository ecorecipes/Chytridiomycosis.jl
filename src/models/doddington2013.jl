# Doddington, Bosch, Oliver, Grassly, Garcia, Schmidt, Garner & Fisher (2013) Ecology
# 94:1795-1804, doi:10.1890/12-1270.1. Mallorcan midwife toad (Alytes muletensis) stages and
# free-living Bd zoospores in one pond. Equations B.1-B.8 and Table B1 are in Appendix B; the
# stochastic model that produced the paper's results is the authors' Berkeley Madonna code
# (Supplement 1, model_code.txt). Both are on figshare (doi:10.6084/m9.figshare.c.3306075, CC0)
# and were checked against this file. Time is in YEARS; zoospores are numbers in the pond.
#
# Printed model (Appendix B), with T = t - floor(t) the time of year:
#   dAX/dt = (1 - AX/K)(εJ JX + εJ JY) - μA AX                                    (B.1)
#   dLX/dt = ψ AX - (μL + εL + ν Z + δ(T - tM2)) LX                                (B.2)
#   dOX/dt = δ(T - tM2) LX - (μO + εO + ν Z) OX                                    (B.3)
#   dJX/dt = εL LX + εO OX - (μJ + εJ) JX                                         (B.4)
#   dLY/dt = ν Z LX - (εL + μL + δ(T - tM2)) LY                                    (B.5)
#   dOY/dt = ν Z OX + δ(T - tM2) LY - (εO + μO) OY                                 (B.6)
#   dJY/dt = εL LY + εO OY - (μJ + α + εJ) JY                                     (B.7)
#   dZ/dt  = ρ LY + ρ OY - μZ Z                                                    (B.8)
# ψ is nonzero only in (tB1, tB2), εL only in (tM1, tM2), εO only in (tOM, tM2), and
# "juveniles only leave due to εJ after a year has passed". Infected juveniles mature into
# uninfected adults.
#
# Juvenile age classes. The code implements "after a year" with two classes per juvenile
# type: this year's juveniles (JuvX) join the class that can mature (JMatX) at the time of year
# tAge = 0.001 (MatSwitch). It also keeps juveniles from tadpoles and from overwintering
# tadpoles apart. Both models here use the code's 14 states:
#   u = [A, LX, OX, JX, JMX, OJX, OJMX, LY, OY, JY, JMY, OJY, OJMY, Z]
# (J = this year's juveniles from tadpoles, JM = those that can mature, OJ/OJM the same for
# juveniles from overwintering tadpoles).
#
# Printed versus coded maturation (`maturation` keyword):
#   :as_printed  JM and OJM mature at the constant εJ (Table B1, 12.3 yr⁻¹) all year;
#   :as_coded    the code's JMR = TMR and OWJMR = OWTMR: JM mature at εL only inside
#                (tM1, tM2) and OJM at εO only inside (tOM, tM2). This is the default,
#                because the paper's results come from the code.
#
# Choices where the code and the paper are silent or ambiguous (documented, not published):
# - Switch times. The code fires the overwintering switch (OWSwitch), the juvenile ageing
#   (MatSwitch) and the zoospore introduction (ZspIntroduce) on steps with
#   `threshold < time < threshold + DT`. With DT = 0.0005 every threshold (0.679, 0.001,
#   15.42) is itself a multiple of DT, so no grid time satisfies the strict inequalities
#   exactly and whether a switch fires depends on Berkeley Madonna's rounding of TIME. Here
#   each switch fires once, on the first step whose time is at or after the threshold.
# - Adult recruitment. The code draws binomial(1 - AX/K, maturing juveniles). The probability
#   is clamped to [0, 1] (it is negative if AX > K).
# - The deterministic model takes the code's 14 states and time windows, with the jumps
#   (δ terms and ageing) as callbacks, and the introduction as setting Z to Z_intro.
# - The overwintering death rate: the code uses TadDeath (μL) where Table B1 has μO; they are
#   equal (0.05 yr⁻¹). Both models here use μO for overwintering tadpoles.
#
# Not implemented: the SIR model with temperature-dependent recovery fitted in Fig. 4 (the
# datalogger series are not deposited and the mapping of Eq. 1 to a recovery rate is not
# printed), the fungal-load individual-based model (Fig. B3; no code), and the state-space
# trend model of Appendix A (statistical; counts not deposited).

const DODDINGTON_STATES = (:A, :LX, :OX, :JX, :JMX, :OJX, :OJMX, :LY, :OY, :JY, :JMY, :OJY,
    :OJMY, :Z)

"""
    DoddingtonParameters(; tB1, tB2, tM1, tM2, tOM, tAge, ψ, εL, εO, εJ, μL, μO, μJ, μA, K,
                         ρ, ν, μZ, α)
    DoddingtonParameters(:doddington2013; set = "", timing = :code, kw...)

Parameters of the Doddington et al. (2013) model of Mallorcan midwife toads and Bd zoospores
in one pond (Appendix B, Table B1). Rates are per **year**.

- Timing, as fractions of the year: `tB1`–`tB2` breeding season, `tM1`–`tM2` tadpole
  metamorphosis, `tOM` start of overwintering-tadpole metamorphosis, `tAge` the time of year at
  which last year's juveniles become able to mature.
- Host: `ψ` tadpoles per adult (breeding season), `εL`, `εO`, `εJ` metamorphosis and
  maturation rates, `μL`, `μO`, `μJ`, `μA` natural death rates, `K` the adult population limit
  in the density dependence `1 - A/K`.
- Pathogen: `ρ` zoospores released per infected tadpole, `ν` infection rate per zoospore in
  the pond, `μZ` zoospore death rate, `α` infection-specific death rate of juveniles.

The keyed constructor reads Table B1 from the compendium. `set` selects a figure's values
(`"fig3a"`, `"fig3b"`, `"appB_a"`, `"appB_b"`, `"appB_c"`); keywords override single fields
(for example `ν` and `ρ` for a cell of Fig. 3). `timing = :code` takes the code's fractions of
the year (0.249, …); `timing = :table` takes Table B1's days of the year divided by 365.
`tAge` is in the code only.
"""
Base.@kwdef struct DoddingtonParameters
    tB1::Float64
    tB2::Float64
    tM1::Float64
    tM2::Float64
    tOM::Float64
    tAge::Float64
    ψ::Float64
    εL::Float64
    εO::Float64
    εJ::Float64
    μL::Float64
    μO::Float64
    μJ::Float64
    μA::Float64
    K::Float64
    ρ::Float64
    ν::Float64
    μZ::Float64
    α::Float64
end

with(p::DoddingtonParameters; kw...) = DoddingtonParameters(;
    (f => Float64(get(kw, f, getfield(p, f))) for f in fieldnames(DoddingtonParameters))...)

time_unit(::Type{<:DoddingtonParameters}) = :year

# A compendium record of `key`: from `set` if it has one, else from the base set.
function _doddington_record(key::Symbol, sym::Symbol, set::AbstractString)
    recs = parameter_records(key)
    i = isempty(set) ? nothing : findfirst(r -> r.symbol === sym && r.set == set, recs)
    i === nothing && (i = findfirst(r -> r.symbol === sym && isempty(r.set), recs))
    i === nothing && throw(ArgumentError("$key has no $sym (set \"$set\")"))
    return recs[i]
end

function DoddingtonParameters(key::Symbol; set::AbstractString = "", timing::Symbol = :code,
        kw...)
    timing in (:code, :table) || throw(ArgumentError("timing must be :code or :table"))
    isempty(set) || any(r -> r.set == set, parameter_records(key)) ||
        throw(ArgumentError("$key has no parameter set \"$set\""))
    rec(s) = _doddington_record(key, s, set)
    tset = timing === :code ? "code" : ""
    tim(s) = in_time_unit(_doddington_record(key, s, tset), :year)
    rate(s) = in_time_unit(rec(s), :year)
    p = DoddingtonParameters(; tB1 = tim(:tB1), tB2 = tim(:tB2), tM1 = tim(:tM1),
        tM2 = tim(:tM2), tOM = tim(:tOM),
        tAge = in_time_unit(_doddington_record(key, :tAge, "code"), :year),
        ψ = rate(:ψ), εL = rate(:εL), εO = rate(:εO), εJ = rate(:εJ), μL = rate(:μL),
        μO = rate(:μO), μJ = rate(:μJ), μA = rate(:μA), K = rec(:K).value,
        ρ = rate(:ρ), ν = rate(:ν), μZ = rate(:μZ), α = rate(:α))
    return isempty(kw) ? p : with(p; kw...)
end

# Run settings of the authors' code (Supplement 1), read from the compendium.
function _doddington_settings(key::Symbol = :doddington2013)
    r(s) = _doddington_record(key, s, "")
    return (; A0 = r(:A0).value, Z_intro = r(:Z_intro).value,
        t_intro = in_time_unit(r(:t_intro), :year), t_end = in_time_unit(r(:t_end), :year),
        dt = in_time_unit(r(:dt), :year), runs = Int(r(:n_runs).value))
end

# Seasonal rates at time t: (ψ, εL, εO, maturation of JM, maturation of OJM).
function _doddington_rates(p::DoddingtonParameters, t::Real, maturation::Symbol)
    day = mod(t, 1)
    ψt = p.tB1 < day < p.tB2 ? p.ψ : 0.0
    εLt = p.tM1 < day < p.tM2 ? p.εL : 0.0
    εOt = p.tOM < day < p.tM2 ? p.εO : 0.0
    mL, mO = maturation === :as_coded ? (εLt, εOt) : (p.εJ, p.εJ)
    return ψt, εLt, εOt, mL, mO
end

function _check_maturation(m::Symbol)
    m in (:as_coded, :as_printed) ||
        throw(ArgumentError("maturation must be :as_coded or :as_printed; got :$m"))
    return m
end

function doddington_rhs!(du, u, pm, t)
    p, maturation = pm
    A, LX, OX, JX, JMX, OJX, OJMX, LY, OY, JY, JMY, OJY, OJMY, Z = u
    ψt, εLt, εOt, mL, mO = _doddington_rates(p, t, maturation)
    λ = p.ν * Z
    du[1] = (1 - A / p.K) * (mL * (JMX + JMY) + mO * (OJMX + OJMY)) - p.μA * A
    du[2] = ψt * A - (p.μL + εLt + λ) * LX
    du[3] = -(p.μO + εOt + λ) * OX
    du[4] = εLt * LX - p.μJ * JX
    du[5] = -(p.μJ + mL) * JMX
    du[6] = εOt * OX - p.μJ * OJX
    du[7] = -(p.μJ + mO) * OJMX
    du[8] = λ * LX - (p.μL + εLt) * LY
    du[9] = λ * OX - (p.μO + εOt) * OY
    du[10] = εLt * LY - (p.μJ + p.α) * JY
    du[11] = -(p.μJ + p.α + mL) * JMY
    du[12] = εOt * OY - (p.μJ + p.α) * OJY
    du[13] = -(p.μJ + p.α + mO) * OJMY
    du[14] = p.ρ * (LY + OY) - p.μZ * Z
    return nothing
end

# The δ(T - tM2) jump (tadpoles become overwintering tadpoles).
function _doddington_overwinter!(u)
    u[3] += u[2]; u[2] = 0
    u[9] += u[8]; u[8] = 0
    return u
end
# Ageing: this year's juveniles join the classes that can mature.
function _doddington_age!(u)
    for (j, jm) in ((4, 5), (6, 7), (10, 11), (12, 13))
        u[jm] += u[j]; u[j] = 0
    end
    return u
end

"""
    DoddingtonODEProblem(p::DoddingtonParameters; A0, Z_intro, t_intro, tspan,
                         maturation = :as_coded, u0 = nothing)

The deterministic model of Doddington et al. (2013), Appendix B eqs B.1–B.8, in years, on the
14 states `[A, LX, OX, JX, JMX, OJX, OJMX, LY, OY, JY, JMY, OJY, OJMY, Z]` of the authors' code
(juveniles are split by age and by tadpole type; see the source header). The jumps are
callbacks: at time of year `tM2` tadpoles become overwintering tadpoles (the δ terms); at
`tAge` last year's juveniles become able to mature; at `t_intro` the zoospore count is set to
`Z_intro`.

`A0`, `Z_intro`, `t_intro` and `tspan = (0, t_end)` default to the settings of the authors'
code (150 adults, 20 000 zoospores at 15.42 yr, 50 yr) from the compendium. `maturation`
is `:as_coded` (the code's seasonal maturation) or `:as_printed` (constant `εJ`). The paper's
results come from the stochastic version, [`DoddingtonStochasticProblem`](@ref).
"""
struct DoddingtonODEProblem
    p::DoddingtonParameters
    u0::Vector{Float64}
    Z_intro::Float64
    t_intro::Float64
    tspan::Tuple{Float64,Float64}
    maturation::Symbol
end
function DoddingtonODEProblem(p::DoddingtonParameters; A0 = nothing, Z_intro = nothing,
        t_intro = nothing, tspan = nothing, maturation::Symbol = :as_coded, u0 = nothing)
    s = _doddington_settings()
    A0 = something(A0, s.A0)
    x0 = u0 === nothing ? vcat(Float64(A0), zeros(13)) : collect(Float64, u0)
    length(x0) == 14 && all(>=(0), x0) ||
        throw(ArgumentError("u0 must have 14 non-negative entries"))
    tsp = something(tspan, (0.0, s.t_end))
    return DoddingtonODEProblem(p, x0, Float64(something(Z_intro, s.Z_intro)),
        Float64(something(t_intro, s.t_intro)), (Float64(tsp[1]), Float64(tsp[2])),
        _check_maturation(maturation))
end

# Event times in (t0, t1]: (time, kind) with kind :overwinter, :age or :intro.
function _doddington_events(p::DoddingtonParameters, t0, t1, t_intro)
    ev = Tuple{Float64,Symbol}[]
    for k in floor(Int, t0):ceil(Int, t1)
        for (f, kind) in ((p.tM2, :overwinter), (p.tAge, :age))
            t = k + f
            t0 < t <= t1 && push!(ev, (t, kind))
        end
    end
    t0 < t_intro <= t1 && push!(ev, (t_intro, :intro))
    return sort!(ev)
end

function to_ode_problem(prob::DoddingtonODEProblem)
    p = prob.p
    t0, t1 = prob.tspan
    ev = _doddington_events(p, t0, t1, prob.t_intro)
    times = unique(first.(ev))
    function affect!(integ)
        for (t, kind) in ev
            t == integ.t || continue
            kind === :overwinter && _doddington_overwinter!(integ.u)
            kind === :age && _doddington_age!(integ.u)
            kind === :intro && (integ.u[14] = prob.Z_intro)
        end
        return nothing
    end
    cb = DiscreteCallback((u, t, integ) -> t in times, affect!; save_positions = (true, true))
    # the seasonal windows switch rates on and off: stop there too
    edges = Float64[]
    for k in floor(Int, t0):ceil(Int, t1), f in (p.tB1, p.tB2, p.tM1, p.tOM)
        t0 < k + f < t1 && push!(edges, k + f)
    end
    return ODEProblem(doddington_rhs!, copy(prob.u0), prob.tspan, (p, prob.maturation);
        callback = cb, tstops = sort!(vcat(times, edges)))
end

CommonSolve.solve(prob::DoddingtonODEProblem, alg = Tsit5(); abstol = 1e-8, reltol = 1e-8,
    kwargs...) = solve(to_ode_problem(prob), alg; abstol, reltol, kwargs...)

"""
    DoddingtonStochasticProblem(p::DoddingtonParameters; A0, Z_intro, t_intro, tspan, dt,
                                maturation = :as_coded, saveat = 0.01)

The stochastic model of Doddington et al. (2013): a line-by-line port of the authors'
Berkeley Madonna code (Supplement 1). Each step of length `dt` (0.0005 yr) draws every
event with a binomial (or, for zoospore release, Poisson) distribution, as the code does.
Competing exits (death, metamorphosis, infection) are split with nested binomials in the
code's order. States are integers, in the order of [`DoddingtonODEProblem`](@ref).

`solve(prob; rng)` returns a named tuple `(t, u, outcome)`: save times, the state at each,
and the outcome at the end of the run (see [`doddington_outcome`](@ref)). Defaults come from
the code: 150 adults, 20 000 zoospores introduced at 15.42 yr, a run of 50 yr. See the
source header for how the switch times are handled.
"""
struct DoddingtonStochasticProblem
    p::DoddingtonParameters
    u0::Vector{Int}
    Z_intro::Int
    t_intro::Float64
    tspan::Tuple{Float64,Float64}
    dt::Float64
    maturation::Symbol
    saveat::Float64
end
function DoddingtonStochasticProblem(p::DoddingtonParameters; A0 = nothing, Z_intro = nothing,
        t_intro = nothing, tspan = nothing, dt = nothing, maturation::Symbol = :as_coded,
        saveat::Real = 0.01, u0 = nothing)
    s = _doddington_settings()
    x0 = u0 === nothing ? vcat(round(Int, something(A0, s.A0)), zeros(Int, 13)) :
         collect(Int, u0)
    length(x0) == 14 && all(>=(0), x0) ||
        throw(ArgumentError("u0 must have 14 non-negative entries"))
    tsp = something(tspan, (0.0, s.t_end))
    h = Float64(something(dt, s.dt))
    h > 0 || throw(ArgumentError("dt must be positive"))
    return DoddingtonStochasticProblem(p, x0, round(Int, something(Z_intro, s.Z_intro)),
        Float64(something(t_intro, s.t_intro)), (Float64(tsp[1]), Float64(tsp[2])), h,
        _check_maturation(maturation), Float64(saveat))
end

@inline _bin(rng, n::Integer, q::Real) =
    (n <= 0 || q <= 0) ? 0 : q >= 1 ? Int(n) : rand(rng, Binomial(Int(n), Float64(q)))
@inline _pois(rng, m::Real) = m <= 0 ? 0 : rand(rng, Poisson(m))

# The first step index i (time t0 + i dt) at or after time τ, with a small tolerance.
_first_step(t0, dt, τ) = max(0, ceil(Int, (τ - t0) / dt - 1e-9))

# One step of the code. `ow`, `ag` and `intro` say whether the switches fire on this step.
function _doddington_step!(u, p::DoddingtonParameters, t, dt, rng, maturation, ow, ag,
        intro, Z_intro)
    A, LX, OX, JX, JMX, OJX, OJMX, LY, OY, JY, JMY, OJY, OJMY, Z = u
    ψt, TMR, OWTMR, JMR, OWJMR = _doddington_rates(p, t, maturation)
    λ = Z * p.ν
    # uninfected tadpoles
    nLleave = _bin(rng, LX, dt * (p.μL + TMR + λ))
    nLnotdead = _bin(rng, nLleave, (TMR + λ) / (p.μL + TMR + λ))
    nLmorph = _bin(rng, nLnotdead, TMR + λ == 0 ? 0.0 : TMR / (TMR + λ))
    nLinf = nLnotdead - nLmorph
    nhatch = ψt > 0 ? _bin(rng, A, dt * ψt) : 0
    # uninfected overwintering tadpoles
    nOleave = _bin(rng, OX, dt * (p.μO + OWTMR + λ))
    nOnotdead = _bin(rng, nOleave, (OWTMR + λ) / (p.μO + OWTMR + λ))
    nOmorph = _bin(rng, nOnotdead, OWTMR + λ == 0 ? 0.0 : OWTMR / (OWTMR + λ))
    nOinf = nOnotdead - nOmorph
    # uninfected juveniles
    dJX = _bin(rng, JX, dt * p.μJ)
    nJMleave = _bin(rng, JMX, dt * (p.μJ + JMR))
    nJMmat = _bin(rng, nJMleave, JMR / (p.μJ + JMR))
    dOJX = _bin(rng, OJX, dt * p.μJ)
    nOJMleave = _bin(rng, OJMX, dt * (p.μJ + OWJMR))
    nOJMmat = _bin(rng, nOJMleave, OWJMR / (p.μJ + OWJMR))
    # infected tadpoles
    nLYleave = _bin(rng, LY, dt * (p.μL + TMR))
    nLYmorph = _bin(rng, nLYleave, TMR / (p.μL + TMR))
    nOYleave = _bin(rng, OY, dt * (p.μO + OWTMR))
    nOYmorph = _bin(rng, nOYleave, OWTMR / (p.μO + OWTMR))
    # infected juveniles
    dJY = _bin(rng, JY, dt * (p.μJ + p.α))
    nJMYleave = _bin(rng, JMY, dt * (p.μJ + p.α + JMR))
    nJMYmat = _bin(rng, nJMYleave, JMR / (p.μJ + p.α + JMR))
    dOJY = _bin(rng, OJY, dt * (p.μJ + p.α))
    nOJMYleave = _bin(rng, OJMY, dt * (p.μJ + p.α + OWJMR))
    nOJMYmat = _bin(rng, nOJMYleave, OWJMR / (p.μJ + p.α + OWJMR))
    # adults
    recruits = _bin(rng, nJMmat + nOJMmat + nJMYmat + nOJMYmat, clamp(1 - A / p.K, 0, 1))
    dA = _bin(rng, A, dt * p.μA)
    # zoospores
    Znew = intro ? Z_intro :
           Z + _pois(rng, dt * LY * p.ρ) + _pois(rng, dt * OY * p.ρ) - _bin(rng, Z, dt * p.μZ)
    # updates (all from the old state, as Berkeley Madonna's `next` does)
    u[1] = A + recruits - dA
    u[2] = ow ? 0 : LX + nhatch - nLleave
    u[3] = ow ? OX + LX : OX - nOleave
    u[4] = ag ? 0 : JX + nLmorph - dJX
    u[5] = ag ? JMX + JX : JMX - nJMleave
    u[6] = ag ? 0 : OJX + nOmorph - dOJX
    u[7] = ag ? OJMX + OJX : OJMX - nOJMleave
    u[8] = ow ? 0 : LY + nLinf - nLYleave
    u[9] = ow ? OY + LY : OY + nOinf - nOYleave
    u[10] = ag ? 0 : JY + nLYmorph - dJY
    u[11] = ag ? JMY + JY : JMY - nJMYleave
    u[12] = ag ? 0 : OJY + nOYmorph - dOJY
    u[13] = ag ? OJMY + OJY : OJMY - nOJMYleave
    u[14] = Znew
    for k in eachindex(u)   # the code's `limit X >= 0`
        u[k] = max(u[k], 0)
    end
    return u
end

"""
    doddington_outcome(u) -> Symbol

The outcome of a run of Doddington et al. (2013) with final state `u`: `:host_extinction`
when no toads are left (adults, tadpoles or juveniles), `:Bd_extinction` when toads remain
but there are no zoospores and no infected tadpoles, and `:coexistence` otherwise.
"""
function doddington_outcome(u::AbstractVector)
    sum(@view u[1:13]) == 0 && return :host_extinction
    u[14] + u[8] + u[9] == 0 && return :Bd_extinction
    return :coexistence
end

function CommonSolve.solve(prob::DoddingtonStochasticProblem;
        rng::AbstractRNG = Random.default_rng())
    p, dt = prob.p, prob.dt
    t0, t1 = prob.tspan
    nsteps = round(Int, (t1 - t0) / dt)
    every = isfinite(prob.saveat) ? max(1, round(Int, prob.saveat / dt)) : typemax(Int)
    iintro = _first_step(t0, dt, prob.t_intro)
    # the year of the next overwintering switch and of the next ageing
    k_ow = floor(Int, t0); k_ow + p.tM2 < t0 && (k_ow += 1)
    k_ag = floor(Int, t0); k_ag + p.tAge < t0 && (k_ag += 1)
    next_ow = _first_step(t0, dt, k_ow + p.tM2)
    next_ag = _first_step(t0, dt, k_ag + p.tAge)
    u = copy(prob.u0)
    ts = [t0]
    us = [copy(u)]
    for i in 0:(nsteps - 1)
        t = t0 + i * dt
        ow = i == next_ow
        ag = i == next_ag
        ow && (k_ow += 1; next_ow = _first_step(t0, dt, k_ow + p.tM2))
        ag && (k_ag += 1; next_ag = _first_step(t0, dt, k_ag + p.tAge))
        _doddington_step!(u, p, t, dt, rng, prob.maturation, ow, ag, i == iintro,
            prob.Z_intro)
        if (i + 1) % every == 0 || i == nsteps - 1
            push!(ts, t0 + (i + 1) * dt)
            push!(us, copy(u))
        end
    end
    return (; t = ts, u = us, outcome = doddington_outcome(u))
end

"""
    doddington_outcomes(p::DoddingtonParameters; runs, rng, kwargs...) -> NamedTuple

Run [`DoddingtonStochasticProblem`](@ref) `runs` times (default 500, as in the paper) and
count the outcomes: `(; host_extinction, Bd_extinction, coexistence, runs, host_survival)`,
where `host_survival` is the fraction of runs with toads left at the end (the quantity
mapped in Fig. 3). Other keywords go to the problem constructor.
"""
function doddington_outcomes(p::DoddingtonParameters; runs::Union{Nothing,Integer} = nothing,
        rng::AbstractRNG = Random.default_rng(), kwargs...)
    n = something(runs, _doddington_settings().runs)
    prob = DoddingtonStochasticProblem(p; saveat = Inf, kwargs...)
    c = Dict(:host_extinction => 0, :Bd_extinction => 0, :coexistence => 0)
    for _ in 1:n
        c[solve(prob; rng).outcome] += 1
    end
    return (; host_extinction = c[:host_extinction], Bd_extinction = c[:Bd_extinction],
        coexistence = c[:coexistence], runs = n, host_survival = 1 - c[:host_extinction] / n)
end

"""
    doddington_grid(key = :doddington2013; n = 10) -> (; ν, ρ)

The axes of Fig. 3 of Doddington et al. (2013): `n` log-spaced values of the transmission
parameter `ν` and of the zoospore release rate `ρ` between the ends of their Table B1 ranges
(2.57e-6–7.27e-5 and 1e4–6.6e6). With `n = 10` they are the printed tick labels.
"""
function doddington_grid(key::Symbol = :doddington2013; n::Integer = 10)
    ends(s) = _doddington_range(key, s)
    lν, hν = ends(:ν)
    lρ, hρ = ends(:ρ)
    lg(a, b) = exp.(range(log(a), log(b); length = n))
    return (; ν = lg(lν, hν), ρ = lg(lρ, hρ))
end

# The ends of a Table B1 range, from the compendium file (the reader keeps only values).
function _doddington_range(key::Symbol, sym::Symbol)
    d = TOML.parsefile(joinpath(@__DIR__, "..", "..", "data", "compendium", "$key.toml"))
    for r in d["parameter"]
        r["symbol"] == String(sym) && !haskey(r, "set") && haskey(r, "uncertainty") &&
            return Float64(r["uncertainty"]["lo"]), Float64(r["uncertainty"]["hi"])
    end
    throw(ArgumentError("$key has no range for $sym"))
end

"""
    excess_death_rate(deaths, time, deaths0, time0; z = 1.96) -> (; rate, lo, hi)

The isolate-specific death rate of Doddington et al. (2013), Table 1: the death rate
`deaths/time` of an exposed group minus that of the control group, with a confidence interval
from the standard error of a difference in rates, `sqrt(deaths/time² + deaths0/time0²)`.
"""
function excess_death_rate(deaths::Real, time::Real, deaths0::Real, time0::Real;
        z::Real = 1.96)
    r = deaths / time - deaths0 / time0
    se = sqrt(deaths / time^2 + deaths0 / time0^2)
    return (; rate = r, lo = r - z * se, hi = r + z * se)
end

"""
    doddington_exposure_trial(key = :doddington2013) -> NamedTuple

Table 1 of Doddington et al. (2013) recomputed from its counts in the compendium: group
death rates (`control`, `CAPE`, `GPL`, per year) and the isolate-specific excess death rates
with 95% intervals (`CAPE_excess`, `GPL_excess`, from [`excess_death_rate`](@ref)).
"""
function doddington_exposure_trial(key::Symbol = :doddington2013)
    v(s) = _doddington_record(key, s, "").value
    d0, t0 = v(:deaths_control), v(:time_control)
    return (; control = d0 / t0, CAPE = v(:deaths_CAPE) / v(:time_CAPE),
        GPL = v(:deaths_GPL) / v(:time_GPL),
        CAPE_excess = excess_death_rate(v(:deaths_CAPE), v(:time_CAPE), d0, t0),
        GPL_excess = excess_death_rate(v(:deaths_GPL), v(:time_GPL), d0, t0))
end

# The most frequent outcome (as a String) of `runs` stochastic runs, for the validation targets.
function _doddington_mode(p::DoddingtonParameters; runs::Integer = 100, seed::Integer = 1,
        kwargs...)
    o = doddington_outcomes(p; runs, rng = Random.Xoshiro(seed), kwargs...)
    return String(argmax(s -> getfield(o, s), (:host_extinction, :Bd_extinction, :coexistence)))
end
