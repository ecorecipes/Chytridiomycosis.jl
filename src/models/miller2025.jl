# Miller, Lydeamore, Flegg, Berger, Skerratt, Waddle & Campbell (2025) arXiv:2503.06846 (v3),
# "Sunlight-heated refugia protect frogs from chytridiomycosis". A continuous-time Markov
# chain (CTMC) for Bd in outdoor mesocosms of Litoria aurea with shaded or sunlight-heated
# ("unshaded") refugia. Checked against the paper (§2.1, Table 1, Table S1, eqs S1–S8) and
# the authors' R code (Zenodo 10.5281/zenodo.15179344 v1.0, fns_stochastic_model.R).
# Units: weeks; frogs per mesocosm. Wiki: lab/wiki/models/bd-miller2025-ctmc-refugia.md
#
# States, in the authors' order (their names in brackets):
#   S_U [S], I_U1 [I1], I_U2 [I2], I_U3 [I3], S_V [SV], I_V1 [IV1], I_V2 [IV2], I_V3 [IV3]
# U = no prior infection, V = recovered from a prior infection ("vaccinated").
# The 16 transitions (Table 1; columns of the authors' state-change matrix nu):
#   1  S_U -> I_U1   β_h S_U F / (N − 1)        2  S_U lost   μ S_U
#   3  I_U1 -> I_U2  γ1 I_U1                    4  I_U1 lost  μ I_U1
#   5  I_U2 -> I_U3  γ2 I_U2                    6  I_U2 lost  μ I_U2
#   7  I_U3 -> S_V   ω I_U3                     8  I_U3 lost  μ I_U3
#   9  S_V -> I_V1   α β_h S_V F / (N − 1)      10 S_V lost   μ S_V
#   11 I_V1 -> I_V2  γ1 I_V1                    12 I_V1 lost  μ I_V1
#   13 I_V2 -> I_V3  γ2 I_V2                    14 I_V2 lost  μ I_V2
#   15 I_V3 -> S_V   ω I_V3                     16 I_V3 lost  μ I_V3
# with F = Σ m_kj I_kj (Eq. 1 without β) and N the number of frogs present.
#
# Behaviour of the authors' simulator kept here:
# - The infection propensity is 0/0 when one frog is left. GillespieSSA (0.6.2, ssa()) sets
#   NaN propensities to 0, so infection stops at N = 1.
# - get_timepoints reads the state at the last event strictly before each observation week,
#   which for a CTMC is the state at that week with probability one.
# - Simulations stop at t_f = 15 weeks, or earlier if all frogs are lost.
# - Both environments start from x0 = (5, 5, 0, 0, 5, 5, 0, 0) and run independently.

const REFUGIA_STATES = (:S_U, :I_U1, :I_U2, :I_U3, :S_V, :I_V1, :I_V2, :I_V3)
const REFUGIA_ENVIRONMENTS = (:shaded, :unshaded)
# (from, to) state indices of the 16 transitions; to = 0 means the frog is lost.
const REFUGIA_TRANSITIONS = ((1, 2), (1, 0), (2, 3), (2, 0), (3, 4), (3, 0), (4, 5), (4, 0),
                             (5, 6), (5, 0), (6, 7), (6, 0), (7, 8), (7, 0), (8, 5), (8, 0))

"""
    RefugiaParameters(; β_sh, β_un, α, ω, μ, γ1, γ2, m)
    RefugiaParameters(:miller2025; set = "", kwargs...)

Parameters of the Miller et al. (2025) refugia CTMC, as rates per week:

- `β_sh`, `β_un`: transmissibility with shaded and with sunlight-heated ("unshaded") refugia;
- `α`: susceptibility of previously infected ("vaccinated") frogs relative to naive frogs;
- `ω`: recovery rate from `I_k,3` to `S_V`;
- `μ`: loss rate (death or not found) from every state;
- `γ1`, `γ2`: progression rates `I_k,1 → I_k,2 → I_k,3`;
- `m`: infectiousness of `(I_U1, I_U2, I_U3, I_V1, I_V2, I_V3)` relative to `I_U3`.

The compendium constructor reads every value from `data/compendium/miller2025.toml`. The
base set holds the fixed constants (Table S1) and the posterior means of Table 3. The set
`"synthetic"` holds the Table S2 values used for Supplement S2. Keyword arguments override
single values, for example `RefugiaParameters(:miller2025; β_sh = 0.5)`.
"""
Base.@kwdef struct RefugiaParameters
    β_sh::Float64
    β_un::Float64
    α::Float64
    ω::Float64
    μ::Float64
    γ1::Float64
    γ2::Float64
    m::NTuple{6,Float64}
end

function RefugiaParameters(key::Symbol; set::AbstractString = "", kwargs...)
    q = parameter_set(key; set)
    base = (; β_sh = q.β_sh, β_un = q.β_un, α = q.α, ω = q.ω, μ = q.μ, γ1 = q.γ1, γ2 = q.γ2,
            m = (q.m_U1, q.m_U2, q.m_U3, q.m_V1, q.m_V2, q.m_V3))
    return RefugiaParameters(; merge(base, values(kwargs))...)
end

time_unit(::Type{<:RefugiaParameters}) = :week

"""
    transmissibility(p::RefugiaParameters, environment) -> β

`p.β_sh` for `environment = :shaded`, `p.β_un` for `:unshaded` (sunlight-heated refugia).
"""
function transmissibility(p::RefugiaParameters, environment::Symbol)
    environment === :shaded && return p.β_sh
    environment === :unshaded && return p.β_un
    throw(ArgumentError("environment must be one of $REFUGIA_ENVIRONMENTS"))
end

"""
    refugia_initial_state(key = :miller2025) -> NTuple{8,Int}

Frogs per mesocosm at the start (§2.1.2), in the order of the states `S_U, I_U1, I_U2, I_U3,
S_V, I_V1, I_V2, I_V3`: `(5, 5, 0, 0, 5, 5, 0, 0)`, read from the compendium.
"""
function refugia_initial_state(key::Symbol = :miller2025)
    key === :miller2025 && return _miller_design().u0
    q = parameter_set(key)
    return (round(Int, q.S_U0), round(Int, q.I_U1_0), 0, 0, round(Int, q.S_V0),
            round(Int, q.I_V1_0), 0, 0)
end

# The design constants (initial state, T_f, survey weeks), read once from the compendium and
# the bundled data so that constructing many problems stays cheap.
const _MILLER_DESIGN = Ref{Any}(nothing)
function _miller_design()
    if _MILLER_DESIGN[] === nothing
        q = parameter_set(:miller2025)
        u0 = (round(Int, q.S_U0), round(Int, q.I_U1_0), 0, 0, round(Int, q.S_V0),
              round(Int, q.I_V1_0), 0, 0)
        _MILLER_DESIGN[] = (; u0, T_f = q.T_f, weeks = load_miller2025(:exp).weeks)
    end
    return _MILLER_DESIGN[]::NamedTuple{(:u0, :T_f, :weeks),Tuple{NTuple{8,Int},Float64,Vector{Float64}}}
end

# Force of infection per susceptible frog without the α factor: β F / (N − 1); zero when
# fewer than two frogs are present (GillespieSSA sets the 0/0 propensity to zero).
@inline function _refugia_foi(p::RefugiaParameters, β, x)
    N = x[1] + x[2] + x[3] + x[4] + x[5] + x[6] + x[7] + x[8]
    N > 1 || return zero(float(eltype(x)))
    m = p.m
    F = m[1] * x[2] + m[2] * x[3] + m[3] * x[4] + m[4] * x[6] + m[5] * x[7] + m[6] * x[8]
    return β * F / (N - 1)
end

@inline function _refugia_rates!(a, p::RefugiaParameters, β, x)
    λ = _refugia_foi(p, β, x)
    μ = p.μ
    @inbounds begin
        a[1] = λ * x[1];          a[2] = μ * x[1]
        a[3] = p.γ1 * x[2];       a[4] = μ * x[2]
        a[5] = p.γ2 * x[3];       a[6] = μ * x[3]
        a[7] = p.ω * x[4];        a[8] = μ * x[4]
        a[9] = p.α * λ * x[5];    a[10] = μ * x[5]
        a[11] = p.γ1 * x[6];      a[12] = μ * x[6]
        a[13] = p.γ2 * x[7];      a[14] = μ * x[7]
        a[15] = p.ω * x[8];       a[16] = μ * x[8]
    end
    return a
end

"""
    refugia_rates(p::RefugiaParameters, environment, x) -> Vector{Float64}

The 16 transition rates (per week) in state `x` (eight counts, ordered as `S_U, I_U1, I_U2,
I_U3, S_V, I_V1, I_V2, I_V3`), in the order of Table 1 and of the authors' reaction list:
infection, then loss, for each state in turn. See also [`refugia_stoichiometry`](@ref).
"""
function refugia_rates(p::RefugiaParameters, environment::Symbol, x)
    length(x) == 8 || throw(ArgumentError("state must have 8 entries"))
    return _refugia_rates!(zeros(16), p, transmissibility(p, environment), x)
end

"""
    refugia_stoichiometry() -> Matrix{Int}

The 8 × 16 state-change matrix of the refugia CTMC. Column `j` is the change caused by
transition `j` of [`refugia_rates`](@ref). It is the matrix `nu` of the authors'
`build_stochastic_system()`.
"""
function refugia_stoichiometry()
    ν = zeros(Int, 8, 16)
    for (j, (from, to)) in enumerate(REFUGIA_TRANSITIONS)
        ν[from, j] -= 1
        to == 0 || (ν[to, j] += 1)
    end
    return ν
end

"""
    RefugiaCTMCProblem(p; environment = :shaded, u0 = refugia_initial_state(),
                       tspan = (0.0, T_f), saveat = refugia_observation_weeks())

The Miller et al. (2025) CTMC for one mesocosm in one environment. Defaults:

- `T_f = 15` weeks, from the compendium;
- `saveat` is the survey weeks of the data: 0, 1, 2, 4, 6, 8, 10 and 15.

`solve(prob; rng)` simulates the chain exactly with the direct method (Gillespie). It
returns `(; t, u, events)`:

- `t` is `saveat`;
- `u` is an 8 × `length(t)` matrix of counts. Column `k` is the state at `t[k]`.
- `events` is the number of transitions.

With `trajectory = true` it also returns `jump_times` and `jump_states`, every transition.
"""
struct RefugiaCTMCProblem
    p::RefugiaParameters
    environment::Symbol
    u0::NTuple{8,Int}
    tspan::Tuple{Float64,Float64}
    saveat::Vector{Float64}
end

function RefugiaCTMCProblem(p::RefugiaParameters; environment::Symbol = :shaded,
        u0 = refugia_initial_state(), tspan = (0.0, _miller_design().T_f),
        saveat = refugia_observation_weeks())
    transmissibility(p, environment)             # validates the environment
    length(u0) == 8 && all(>=(0), u0) || throw(ArgumentError("u0 must be 8 non-negative counts"))
    ts = sort!(Float64.(collect(saveat)))
    all(t -> tspan[1] <= t <= tspan[2], ts) || throw(ArgumentError("saveat must lie in tspan"))
    return RefugiaCTMCProblem(p, environment, Tuple(Int.(u0)), Float64.(tspan), ts)
end

# Direct-method SSA writing the state at each `saveat` time into the columns of `out`.
# `x` (length 8) and `a` (length 16) are work vectors. Returns the number of transitions.
function _refugia_ssa!(out::AbstractMatrix{Int}, x::Vector{Int}, a::Vector{Float64},
        rng::AbstractRNG, p::RefugiaParameters, β::Float64, u0, tspan, saveat;
        jumps = nothing)
    @inbounds for i in 1:8
        x[i] = u0[i]
    end
    t, tf = tspan
    k, ns = 1, length(saveat)
    nev = 0
    @inbounds while true
        _refugia_rates!(a, p, β, x)
        a0 = 0.0
        for j in 1:16
            a0 += a[j]
        end
        tnext = a0 > 0 ? t + Random.randexp(rng) / a0 : Inf
        # surveys up to the next transition see the current state (the authors read the
        # state at the last event strictly before each survey week)
        while k <= ns && saveat[k] <= tnext
            for i in 1:8
                out[i, k] = x[i]
            end
            k += 1
        end
        (tnext > tf || k > ns && jumps === nothing) && break
        r = rand(rng) * a0
        j = 1
        c = a[1]
        while c < r && j < 16
            j += 1
            c += a[j]
        end
        while a[j] == 0 && j > 1          # guard against round-off picking a zero rate
            j -= 1
        end
        from, to = REFUGIA_TRANSITIONS[j]
        x[from] -= 1
        to == 0 || (x[to] += 1)
        t = tnext
        nev += 1
        jumps === nothing || (push!(jumps[1], t); push!(jumps[2], Tuple(x)))
    end
    return nev
end

function CommonSolve.solve(prob::RefugiaCTMCProblem; rng::AbstractRNG = Random.default_rng(),
        trajectory::Bool = false)
    out = zeros(Int, 8, length(prob.saveat))
    jumps = trajectory ? (Float64[], NTuple{8,Int}[]) : nothing
    β = transmissibility(prob.p, prob.environment)
    nev = _refugia_ssa!(out, zeros(Int, 8), zeros(16), rng, prob.p, β, prob.u0, prob.tspan,
        prob.saveat; jumps)
    res = (; t = copy(prob.saveat), u = out, events = nev)
    return trajectory ? merge(res, (; jump_times = jumps[1], jump_states = jumps[2])) : res
end

"""
    RefugiaODEProblem(p; environment = :shaded, u0 = refugia_initial_state(),
                      tspan = (0.0, T_f))

The deterministic version of the refugia model (Supplement S3, eqs S1–S8). It has the
CTMC's rates, with `N(t)` the sum of the eight states. The authors fitted it by least
squares to obtain the synthetic-data parameters (Table S2). `solve(prob; saveat, kwargs...)`
integrates with `Tsit5` and returns the `ODESolution`; `to_ode_problem(prob)` gives the
`ODEProblem`. The force of infection is set to zero when `N(t) ≤ 1`, as in the CTMC.
"""
struct RefugiaODEProblem
    p::RefugiaParameters
    environment::Symbol
    u0::Vector{Float64}
    tspan::Tuple{Float64,Float64}
end

function RefugiaODEProblem(p::RefugiaParameters; environment::Symbol = :shaded,
        u0 = refugia_initial_state(), tspan = (0.0, _miller_design().T_f))
    transmissibility(p, environment)
    length(u0) == 8 || throw(ArgumentError("u0 must have 8 entries"))
    return RefugiaODEProblem(p, environment, Float64.(collect(u0)), Float64.(tspan))
end

function _refugia_ode!(du, u, prm, t)
    p, β = prm
    a = _refugia_rates!(zeros(eltype(u), 16), p, β, u)
    fill!(du, zero(eltype(du)))
    @inbounds for (j, (from, to)) in enumerate(REFUGIA_TRANSITIONS)
        du[from] -= a[j]
        to == 0 || (du[to] += a[j])
    end
    return nothing
end

to_ode_problem(prob::RefugiaODEProblem) =
    ODEProblem(_refugia_ode!, copy(prob.u0), prob.tspan,
               (prob.p, transmissibility(prob.p, prob.environment)))

CommonSolve.solve(prob::RefugiaODEProblem; kwargs...) =
    solve(to_ode_problem(prob), Tsit5(); kwargs...)

# --------------------------------------------------------------------------------------------
# Data: the authors' state counts (Zenodo, CC BY 4.0), bundled in data/miller2025/.

"""
    RefugiaMesocosmData

State counts for the mesocosms of Miller et al. (2025), read by [`load_miller2025`](@ref):

- `weeks`: the survey weeks;
- `mesocosm`: the mesocosm ids;
- `shaded`: `true` for shade-cloth refugia, `false` for sunlight-heated refugia;
- `counts`: an 8 × weeks × mesocosms array of frogs per state, in the order `S_U, I_U1, I_U2,
  I_U3, S_V, I_V1, I_V2, I_V3`.
"""
struct RefugiaMesocosmData
    weeks::Vector{Float64}
    mesocosm::Vector{Int}
    shaded::Vector{Bool}
    counts::Array{Int,3}
end

"""
    miller2025_path(dataset = :exp)

Path of the bundled data of Miller et al. (2025, Zenodo 10.5281/zenodo.15179344, CC BY 4.0).

- `:exp` is the eight outdoor mesocosms of Waddle et al. (2024), assigned to states.
- `:syn` is the authors' synthetic data (Supplement S2).

See `data/miller2025/README.md`.
"""
function miller2025_path(dataset::Symbol = :exp)
    dataset in (:exp, :syn) || throw(ArgumentError("dataset must be :exp or :syn"))
    return joinpath(pkgdir(@__MODULE__), "data", "miller2025", string(dataset) * "_data.csv")
end

const _MILLER_COMPARTMENTS = ("S", "I1", "I2", "I3", "SV", "IV1", "IV2", "IV3")

"""
    load_miller2025(dataset = :exp) -> RefugiaMesocosmData

Read the bundled mesocosm state counts (`:exp` or `:syn`; see [`miller2025_path`](@ref)).
Mesocosms are sorted by id, and missing (mesocosm, week, state) rows count as zero, as in
the authors' `pivot_wider(values_fill = 0)`.
"""
function load_miller2025(dataset::Symbol = :exp)
    path = miller2025_path(dataset)
    lines = filter(!isempty ∘ strip, readlines(path))
    header = split(lines[1], ',')
    function col(name)
        i = findfirst(==(name), header)
        i === nothing && throw(ArgumentError("column $name missing from $path"))
        return i
    end
    im, iw, ic, ish, iN = col("Mesocosm"), col("week"), col("compartment"), col("shaded"),
        col("N")
    rows = [split(l, ',') for l in lines[2:end]]
    meso = sort!(unique(parse(Int, r[im]) for r in rows))
    weeks = sort!(unique(parse(Float64, r[iw]) for r in rows))
    counts = zeros(Int, 8, length(weeks), length(meso))
    shaded = Vector{Union{Missing,Bool}}(missing, length(meso))
    for r in rows
        s = findfirst(==(r[ic]), _MILLER_COMPARTMENTS)
        s === nothing && throw(ArgumentError("unknown compartment $(r[ic]) in $path"))
        m = searchsortedfirst(meso, parse(Int, r[im]))
        w = searchsortedfirst(weeks, parse(Float64, r[iw]))
        counts[s, w, m] = parse(Int, r[iN])
        sh = r[ish] == "TRUE"
        ismissing(shaded[m]) ? (shaded[m] = sh) :
            shaded[m] == sh || throw(ArgumentError("mesocosm $(meso[m]) changes treatment"))
    end
    return RefugiaMesocosmData(weeks, meso, Bool.(shaded), counts)
end

"""
    refugia_observation_weeks() -> Vector{Float64}

The survey weeks of the experiment, `[0, 1, 2, 4, 6, 8, 10, 15]`, read from the bundled data
(the authors' `sample_times`).
"""
refugia_observation_weeks() = copy(_miller_design().weeks)

"""
    loss_rate(data::RefugiaMesocosmData; N0 = 20) -> μ

The loss rate of Eq. (3) of Miller et al. (2025): `−log(N̄(T_f)/N0)/T_f`, where `N̄(T_f)` is
the mean number of frogs found per mesocosm at the last survey `T_f`.
"""
function loss_rate(d::RefugiaMesocosmData; N0::Real = 20)
    Tf = d.weeks[end]
    NT = sum(d.counts[:, end, :]) / length(d.mesocosm)
    return -log(NT / N0) / Tf
end
