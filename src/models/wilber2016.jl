# Wilber MQ, Langwig KE, Kilpatrick AM, McCallum HI, Briggs CJ (2016) Integral projection
# models for host–parasite systems with an application to amphibian chytrid fungus.
# Methods Ecol Evol 7:1182–1194. doi:10.1111/2041-210X.12561
# Wiki: lab/wiki/models/bd-wilber2016-ipm.md
#
# Checked against the published PDF (eqns 1–13, Table 1, Figs 2–5) and the authors' R code
# and results (github.com/mqwilber/ipm_for_parasites, commit 240b5e67; IPM_functions_for_R.R,
# ipm_vital_rate_parameter_analysis.R, ipm_bd_analysis_for_manuscript.R,
# lower_level_sensitivity_analysis.R, transmission_simulation.R, R0_calculation.R). Data:
# Dryad doi:10.5061/dryad.06sd9 (CC0), bundled as data/wilber2016/. Appendices S1–S3 were not
# retrievable; the code stands in for them.
#
# Model (eqns 7–8), x = ln(ZE swab⁻¹) for infected frogs, T in °C, one step = 3 days:
#   S(t+1)    = S(t) s0 (1 − φ(T)) + ∫ I(x,t) s(x) l(x,T) dx
#   I(x′,t+1) = ∫ I(x,t) s(x) (1 − l(x,T)) G(x′,x,T) dx + S(t) s0 φ(T) G0(x′,T)
#   logit s(x)    = b00 + b10 x                              (eqn 9; temperature-independent)
#   G(x′,x,T)     = Normal(b01 + b11 x + b21 T, v01 exp(2 c01 x))   (eqns 10–11, variance)
#   logit l(x,T)  = b02 + b12 x + b22 T                      (eqn 12)
#   G0(x′,T)      = Normal(b03 + b13 T, v03 exp(2 c03 T))
#   logit φ(T)    = b04 + b14 T
# Discretised by the midpoint rule on [−5, 18] with 100 bins and no eviction correction.
# With density-dependent transmission (eqn 13), φ = 1 − exp(−(a + β ∫ x I(x,t) dx)) with
# a = −ln(1 − ω). The authors' code (density_dependent_zoospores) instead uses
# ln(Σ e^x n(x) + 1), the log of the total zoospore load plus one: LoadTransmission(...;
# as_coded = true). See the validation targets for what this changes.
#
# Printed and coded details that differ (see also the compendium notes):
# - Table 1's survival coefficients (11.824, −0.8605) are not those of the deposited fit or of
#   a refit of the deposited data (11.597, −0.8873). WilberParameters(:wilber2016) uses Table
#   1; set = "deposited" uses the deposit, which is what the authors' figures were made with.
# - The text says 20 frogs at three temperatures, 5 per temperature: the data have 20 frogs at
#   four temperatures (4, 12, 20, 26 °C); 26 °C was excluded from all fits.
# - The bootstrap (monte_carlo_params) stores the survival slope draw as surv_size, which
#   set_temp_params never reads (it reads surv_slope), so the slope is never resampled.
#   growth_rate_bootstrap(...; as_coded = true) reproduces this.
# - The density-dependent scan used ω in 0.01–0.3 in the deposited code, not 0.01–0.6.

const WILBER_FIELDS = (:b00, :b10, :b01, :b11, :b21, :v01, :c01, :b02, :b12, :b22, :b03,
    :b13, :v03, :c03, :b04, :b14, :s0)

"""
    WilberParameters(; b00, b10, b01, b11, b21, v01, c01, b02, b12, b22, b03, b13, v03, c03,
                     b04, b14, s0 = 1, Δt = 3)
    WilberParameters(:wilber2016; set = "", kwargs...)
    WilberParameters(fit::WilberFit; s0 = 1, Δt = 3)

Regression coefficients of the vital-rate functions of the Wilber et al. (2016) Bd–*Rana
muscosa* IPM (Table 1), with `x` the natural log of the load in zoospore equivalents per
swab and `T` the temperature in °C. Probabilities are per time step of `Δt` days (3).

| function | form | coefficients |
|:--|:--|:--|
| survival `s(x)` | `logit s = b00 + b10 x` | `b00`, `b10` |
| growth `G(x′, x, T)` | Normal, mean `b01 + b11 x + b21 T`, variance `v01 exp(2 c01 x)` | `b01`, `b11`, `b21`, `v01`, `c01` |
| loss of infection `l(x, T)` | `logit l = b02 + b12 x + b22 T` | `b02`, `b12`, `b22` |
| initial load `G0(x′, T)` | Normal, mean `b03 + b13 T`, variance `v03 exp(2 c03 T)` | `b03`, `b13`, `v03`, `c03` |
| infection `φ(T)` | `logit φ = b04 + b14 T` | `b04`, `b14` |
| uninfected survival | constant | `s0` |

The keyed constructor reads `data/compendium/wilber2016.toml`: Table 1 as printed by
default, or with `set = "deposited"` the authors' deposited fit at full precision (identical
to Table 1 to the printed digits, except survival: `b00 = 11.597`, `b10 = −0.8873`, which a
refit of the data reproduces). The authors' figures were made with the deposited values.
`WilberParameters(fit)` takes the coefficients of [`fit_wilber2016`](@ref).
"""
struct WilberParameters
    b00::Float64
    b10::Float64
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
    b04::Float64
    b14::Float64
    s0::Float64
    Δt::Float64
    function WilberParameters(b00, b10, b01, b11, b21, v01, c01, b02, b12, b22, b03, b13,
            v03, c03, b04, b14, s0, Δt)
        (v01 > 0 && v03 > 0) || throw(ArgumentError("variances v01 and v03 must be positive"))
        0 <= s0 <= 1 || throw(ArgumentError("s0 must lie in [0, 1]"))
        Δt > 0 || throw(ArgumentError("Δt must be positive"))
        return new(b00, b10, b01, b11, b21, v01, c01, b02, b12, b22, b03, b13, v03, c03, b04,
            b14, s0, Δt)
    end
end

WilberParameters(; b00, b10, b01, b11, b21, v01, c01, b02, b12, b22, b03, b13, v03, c03, b04,
        b14, s0 = 1.0, Δt = 3.0) =
    WilberParameters(b00, b10, b01, b11, b21, v01, c01, b02, b12, b22, b03, b13, v03, c03,
        b04, b14, s0, Δt)

function WilberParameters(key::Symbol; set::AbstractString = "", kwargs...)
    q = parameter_set(key; set)
    base = NamedTuple{WILBER_FIELDS}(Tuple(getfield(q, f) for f in WILBER_FIELDS))
    rec = only(filter(r -> r.symbol === :Δt && r.set == "", parameter_records(key)))
    base = merge(base, (; Δt = in_time_unit(rec, :day)))
    return WilberParameters(; merge(base, values(kwargs))...)
end

time_unit(::Type{<:WilberParameters}) = :day

# Copy with fields replaced.
_with(p::WilberParameters; kw...) =
    WilberParameters(; (f => get(kw, f, getfield(p, f)) for f in fieldnames(WilberParameters))...)

# IPMMesh(:wilber2016; n = 100): the mesh of Wilber et al. (2016), `n` bins on [−5, 18]
# (natural-log ZE per swab), bounds and default `n` from the compendium (documented under
# IPMMesh). The authors' bootstrap, elasticity and R0 code uses n = 50.
function IPMMesh(key::Symbol; n::Union{Nothing,Integer} = nothing)
    q = parameter_set(key)
    return IPMMesh(q.mesh_lower, q.mesh_upper, n === nothing ? round(Int, q.mesh_n) : n)
end

_logistic(z) = 1 / (1 + exp(-z))

# IPMKernel(p::WilberParameters; temperature): the vital-rate functions of the Wilber et al.
# (2016) IPM at `temperature` (°C) (documented under IPMKernel).
function IPMKernel(p::WilberParameters; temperature::Real)
    T = Float64(temperature)
    survival = x -> _logistic(p.b00 + p.b10 * x)
    growth = x -> Normal(p.b01 + p.b11 * x + p.b21 * T, sqrt(p.v01 * exp(2 * p.c01 * x)))
    loss = x -> _logistic(p.b02 + p.b12 * x + p.b22 * T)
    initial = Normal(p.b03 + p.b13 * T, sqrt(p.v03 * exp(2 * p.c03 * T)))
    return IPMKernel(survival, growth, loss, initial, p.s0)
end

"""
    infection_probability(p::WilberParameters, temperature) -> Float64

Probability `φ(T)` that an uninfected frog becomes infected in one 3-day step at `temperature`
(°C), `logit φ = b04 + b14 T` (Table 1). In the experiment frogs were housed singly, so this
is infection from the environment (the inoculum).
"""
infection_probability(p::WilberParameters, temperature::Real) =
    _logistic(p.b04 + p.b14 * temperature)

"""
    load_growth_rate(p::WilberParameters; temperature, x = 0.0) -> Float64

Expected change of ln load per day of an infected frog that survives and stays infected,
`(b01 + (b11 − 1) x + b21 T)/Δt`, at load `x` (ln ZE per swab) and `temperature` (°C). Unlike
the exponential growth of the Briggs et al. (2010) model
([`load_growth_rate(::BriggsParameters, ::Any)`](@ref)), the fitted growth function
(`b11 < 1`) pulls loads towards `x* = (b01 + b21 T)/(1 − b11)`, where the rate is zero.
"""
load_growth_rate(p::WilberParameters; temperature::Real, x::Real = 0.0) =
    (p.b01 + (p.b11 - 1) * x + p.b21 * temperature) / p.Δt

"""
    LoadTransmission(ω, β; as_coded = false)

Density-dependent transmission with an environmental reservoir (Wilber et al. 2016, eqn 13):
the probability that an uninfected host is infected in one step is

    φ = 1 − exp(−(a + β Z)),   a = −ln(1 − ω),

where `ω` is the probability of infection from the environment alone and `Z` is the
pathogen in the host population.

- `as_coded = false` (default): `Z = ∫ x I(x, t) dx`, as printed, the sum of the (log)
  loads of infected hosts.
- `as_coded = true`: `Z = ln(Σ eˣ n(x) + 1)`, the log of the total zoospore load plus one, as
  in the authors' `density_dependent_zoospores`. This is what produced the paper's
  density-dependent results (Figs S7–S8).

Use it as the `transmission` of a [`BdIPMProblem`](@ref);
`infection_probability(tr, mesh, u)` evaluates it for a state `u = [S; n]`.
"""
struct LoadTransmission
    ω::Float64
    β::Float64
    as_coded::Bool
    function LoadTransmission(ω::Real, β::Real, as_coded::Bool)
        0 <= ω < 1 || throw(ArgumentError("ω must lie in [0, 1)"))
        β >= 0 || throw(ArgumentError("β must be non-negative"))
        return new(ω, β, as_coded)
    end
end
LoadTransmission(ω::Real, β::Real; as_coded::Bool = false) = LoadTransmission(ω, β, as_coded)

"""
    infection_probability(tr::LoadTransmission, mesh::IPMMesh, u) -> Float64

The infection probability of [`LoadTransmission`](@ref) `tr` for the state
`u = [S; n_1, …, n_m]` on `mesh`.
"""
function infection_probability(tr::LoadTransmission, mesh::IPMMesh, u::AbstractVector)
    ys = meshpoints(mesh)
    n = @view u[2:end]
    Z = tr.as_coded ? log(sum(exp.(ys) .* n) + 1) : dot(ys, n)
    return 1 - exp(-(-log1p(-tr.ω) + tr.β * Z))
end

"""
    BdIPMProblem(p::WilberParameters; temperature, mesh = IPMMesh(:wilber2016),
                 S0 = 100, I0 = zeros(mesh.n), days = 120, transmission = nothing,
                 integration = :midpoint, eviction = :none)

The Wilber et al. (2016) host–parasite IPM projected for `days` (in steps of `p.Δt` = 3 days)
from `S0` uninfected frogs and infected frogs `I0` (numbers per mesh bin).

- `temperature`: °C, a number or a function of time in days (the kernel is then rebuilt
  each step, for seasonal runs).
- `transmission = nothing`: the infection probability is `φ(T)` of the regression (eqns 7–8,
  density-independent). A [`LoadTransmission`](@ref) makes it depend on the state (eqn 13).
- `integration`, `eviction`: passed to [`discretize`](@ref).

`S0` and `days` default to the compendium values of the paper's summer simulations (100
frogs, 120 days). `solve(prob)` iterates the map directly and returns a
[`BdIPMSolution`](@ref); [`to_discrete_problem`](@ref) gives the same map as a
`SciMLBase.DiscreteProblem`.
"""
struct BdIPMProblem{TT,TR}
    p::WilberParameters
    temperature::TT
    mesh::IPMMesh
    u0::Vector{Float64}
    steps::Int
    transmission::TR
    integration::Symbol
    eviction::Symbol
end

function BdIPMProblem(p::WilberParameters; temperature, mesh::IPMMesh = IPMMesh(:wilber2016),
        S0::Union{Nothing,Real} = nothing, I0::AbstractVector = zeros(mesh.n),
        days::Union{Nothing,Real} = nothing, transmission = nothing,
        integration::Symbol = :midpoint, eviction::Symbol = :none)
    q = parameter_set(:wilber2016)
    S0 = S0 === nothing ? q.S_init : S0
    days = days === nothing ? q.summer_days : days
    length(I0) == mesh.n || throw(DimensionMismatch("I0 must have one entry per mesh bin"))
    (S0 >= 0 && all(>=(0), I0)) || throw(ArgumentError("initial numbers must be non-negative"))
    steps = round(Int, days / p.Δt)
    isapprox(steps * p.Δt, days) || throw(ArgumentError("days must be a multiple of Δt = $(p.Δt)"))
    (transmission === nothing || transmission isa LoadTransmission) ||
        throw(ArgumentError("transmission must be nothing or a LoadTransmission"))
    temperature isa Real || applicable(temperature, 0.0) ||
        throw(ArgumentError("temperature must be a number or a function of time (days)"))
    return BdIPMProblem(p, temperature, mesh, vcat(Float64(S0), Float64.(I0)), steps,
        transmission, integration, eviction)
end

"""
    BdIPMSolution

The result of `solve(::BdIPMProblem)`:

- `t`: times in days, `0, Δt, …`;
- `S`: uninfected frogs at each time;
- `I`: infected frogs per mesh bin (`n × length(t)`);
- `φ`: the infection probability used in each step (`length(t) − 1`);
- `mesh`: the [`IPMMesh`](@ref).

`load_summary(sol.mesh, [sol.S[k]; sol.I[:, k]])` summarises time `k`.
"""
struct BdIPMSolution
    t::Vector{Float64}
    S::Vector{Float64}
    I::Matrix{Float64}
    φ::Vector{Float64}
    mesh::IPMMesh
end

# The map u(t) → u(t + 1) at step k (time k Δt days), with a kernel cache for a fixed
# temperature.
struct _BdIPMMap{P<:BdIPMProblem}
    prob::P
    cache::Base.RefValue{Union{Nothing,DiscretizedIPM}}
end
_BdIPMMap(prob::BdIPMProblem) = _BdIPMMap(prob, Ref{Union{Nothing,DiscretizedIPM}}(nothing))

function _kernel_at(m::_BdIPMMap, tdays)
    prob = m.prob
    if prob.temperature isa Real
        m.cache[] === nothing && (m.cache[] = discretize(IPMKernel(prob.p;
            temperature = prob.temperature), prob.mesh; integration = prob.integration,
            eviction = prob.eviction))
        return m.cache[]::DiscretizedIPM, Float64(prob.temperature)
    end
    T = Float64(prob.temperature(tdays))
    return discretize(IPMKernel(prob.p; temperature = T), prob.mesh;
        integration = prob.integration, eviction = prob.eviction), T
end

function _step_phi(prob::BdIPMProblem, T, u)
    tr = prob.transmission
    return tr === nothing ? infection_probability(prob.p, T) :
           infection_probability(tr, prob.mesh, u)
end

function (m::_BdIPMMap)(du, u, _, k)
    d, T = _kernel_at(m, k * m.prob.p.Δt)
    _ipm_step!(du, d, u, _step_phi(m.prob, T, u))
    return nothing
end

"""
    to_discrete_problem(prob::BdIPMProblem) -> SciMLBase.DiscreteProblem

The IPM map as an in-place `SciMLBase.DiscreteProblem` on the state `[S; n_1, …, n_m]`, with
time counted in steps (`tspan = (0, steps)`; multiply by `Δt` for days). Solve it with
`FunctionMap()` from OrdinaryDiffEq (not a dependency of this package), or call
`solve(prob)` to iterate directly.
"""
to_discrete_problem(prob::BdIPMProblem) =
    DiscreteProblem{true}(_BdIPMMap(prob), copy(prob.u0), (0.0, Float64(prob.steps)))

function CommonSolve.solve(prob::BdIPMProblem)
    f = _BdIPMMap(prob)
    n = prob.mesh.n
    U = Matrix{Float64}(undef, n + 1, prob.steps + 1)
    U[:, 1] .= prob.u0
    φs = Vector{Float64}(undef, prob.steps)
    for k in 1:prob.steps
        u = @view U[:, k]
        d, T = _kernel_at(f, (k - 1) * prob.p.Δt)
        φs[k] = _step_phi(prob, T, u)
        @views _ipm_step!(U[:, k+1], d, u, φs[k])
    end
    t = collect(0:prob.steps) .* prob.p.Δt
    return BdIPMSolution(t, U[1, :], U[2:end, :], φs, prob.mesh)
end

# --- Analyses of the density-independent IPM -----------------------------------------------

function _wilber_matrix(p::WilberParameters, T::Real, mesh::IPMMesh; φ = nothing,
        integration::Symbol = :midpoint, eviction::Symbol = :none)
    d = discretize(IPMKernel(p; temperature = T), mesh; integration, eviction)
    return projection_matrix(d, φ === nothing ? infection_probability(p, T) : φ), d
end

"""
    population_growth_rate(p::WilberParameters; temperature, mesh = IPMMesh(:wilber2016),
                           integration = :midpoint, eviction = :none) -> Float64

The population growth rate `λ` per 3-day step of the density-independent IPM (eqns 7–8;
Fig. 4b) at `temperature` (°C): the dominant eigenvalue of the full matrix with the
uninfected class. There is no reproduction, so `λ ≤ 1` and `λ^(120/3)` is the fraction
surviving a 120-day summer in the long run.
"""
population_growth_rate(p::WilberParameters; temperature::Real,
        mesh::IPMMesh = IPMMesh(:wilber2016), kwargs...) =
    population_growth_rate(_wilber_matrix(p, temperature, mesh; kwargs...)[1])

"""
    stable_load_distribution(p::WilberParameters; temperature, mesh = IPMMesh(:wilber2016))
        -> NamedTuple

The stable load distribution of infected frogs at `temperature` (Fig. 5): the dominant
eigenvector of the full matrix, conditioned on infection. Returns the midpoints `x`, the
probability in each bin `probability` (summing to 1), its `mean`, `variance` and
variance-to-mean ratio `vm_ratio` (of ln ZE, as in the paper), and `uninfected`, the stable
fraction uninfected.
"""
function stable_load_distribution(p::WilberParameters; temperature::Real,
        mesh::IPMMesh = IPMMesh(:wilber2016), kwargs...)
    A, _ = _wilber_matrix(p, temperature, mesh; kwargs...)
    w = stable_distribution(A)
    c = w[2:end] ./ (1 - w[1])
    ys = meshpoints(mesh)
    m = dot(c, ys)
    v = dot(c, ys .^ 2) - m^2
    return (x = ys, probability = c, mean = m, variance = v, vm_ratio = v / m,
            uninfected = w[1])
end

"""
    mean_time_to_death(p::WilberParameters; temperature, mesh = IPMMesh(:wilber2016))
        -> NamedTuple

Expected time to death in days (Fig. 4a) from each state of the density-independent IPM,
with loss of infection and reinfection included: the absorption times of the full matrix
([`absorption_times`](@ref)) times `Δt`. Returns `uninfected` (from the uninfected class),
`x` (mesh midpoints), `infected` (from each load) and the standard deviations `sd_uninfected`
and `sd_infected`. Hosts evicted from the mesh count as dead.
"""
function mean_time_to_death(p::WilberParameters; temperature::Real,
        mesh::IPMMesh = IPMMesh(:wilber2016), kwargs...)
    A, _ = _wilber_matrix(p, temperature, mesh; kwargs...)
    a = absorption_times(A)
    sd = sqrt.(max.(a.variance, 0.0)) .* p.Δt
    return (uninfected = a.mean[1] * p.Δt, x = meshpoints(mesh), infected = a.mean[2:end] .* p.Δt,
            sd_uninfected = sd[1], sd_infected = sd[2:end])
end

"""
    basic_reproduction_number(p::WilberParameters; temperature, β = 9.82e-4, S = 100,
                              mesh = IPMMesh(:wilber2016), method = :full) -> Float64

`R0` of the Wilber et al. (2016) IPM with density-dependent transmission between frogs, as in
the authors' Appendix S1 code (`R0_calculation.R`): each infected frog infects `β S` of `S`
susceptible frogs per step. `β` and `S` default to the compendium values of that code
(`β_R0`, `S_R0`).

- `method = :full`: the next-generation `R0` of the discretised IPM
  (the method of [`basic_reproduction_number`](@ref) for a [`DiscretizedIPM`](@ref)).
- `method = :collapsed`: the SIS approximation `R0 = β s0 S / (1 − s̄(1 − l̄))`, with `s̄`
  and `l̄` the survival and loss probabilities averaged over the stable load distribution of
  the IPM run with `φ = 1 − exp(−β)`.

The authors' code uses a 50-bin mesh (`IPMMesh(:wilber2016; n = 50)`).
"""
function basic_reproduction_number(p::WilberParameters; temperature::Real,
        β::Union{Nothing,Real} = nothing, S::Union{Nothing,Real} = nothing,
        mesh::IPMMesh = IPMMesh(:wilber2016), method::Symbol = :full)
    q = parameter_set(:wilber2016)
    β = β === nothing ? q.β_R0 : β
    S = S === nothing ? q.S_R0 : S
    k = IPMKernel(p; temperature)
    d = discretize(k, mesh)
    if method === :full
        return basic_reproduction_number(d, β * S)
    elseif method === :collapsed
        w = stable_distribution(projection_matrix(d, 1 - exp(-β)))
        c = w[2:end] ./ sum(w[2:end])
        ys = meshpoints(mesh)
        sbar = sum(c .* k.survival.(ys))
        lbar = sum(c .* k.loss.(ys))
        return β * p.s0 * S / (1 - sbar * (1 - lbar))
    end
    throw(ArgumentError("method must be :full or :collapsed"))
end

"""
    vital_rate_elasticities(p::WilberParameters; temperature, δ = 0.001,
                            mesh = IPMMesh(:wilber2016), of = :λ) -> NamedTuple

Elasticities of the population growth rate (`of = :λ`) or of the variance-to-mean ratio of
the stable load distribution (`of = :vm_ratio`) to each regression coefficient (Wilber et
al. 2016, 'Analysing the IPM', Figs S5–S6): `e_b = (f(b(1 + δ)) − f(b)) / (δ f(b))`, keyed by
coefficient (`b00`, `b10`, …, `s0`). `s0 = 1` is perturbed downwards (`s0(1 − δ)`), since a
survival probability cannot exceed 1. The values are signed; the authors plot absolute values.
"""
function vital_rate_elasticities(p::WilberParameters; temperature::Real,
        δ::Union{Nothing,Real} = nothing, mesh::IPMMesh = IPMMesh(:wilber2016), of::Symbol = :λ)
    δ = δ === nothing ? parameter_set(:wilber2016).δ : δ
    f = of === :λ ? (q -> population_growth_rate(q; temperature, mesh)) :
        of === :vm_ratio ? (q -> stable_load_distribution(q; temperature, mesh).vm_ratio) :
        throw(ArgumentError("of must be :λ or :vm_ratio"))
    base = f(p)
    vals = map(WILBER_FIELDS) do s
        b = getfield(p, s)
        h = s === :s0 ? -δ : δ
        (f(_with(p; s => b * (1 + h))) - base) / (h * base)
    end
    return NamedTuple{WILBER_FIELDS}(vals)
end

# --- Data ------------------------------------------------------------------------------------

"""
    wilber2016_path()

Path of the bundled data of the Wilber et al. (2016) temperature experiment
(`data/wilber2016/converted_temp_data.csv`; Dryad doi:10.5061/dryad.06sd9, CC0). See
`data/wilber2016/README.md`.
"""
wilber2016_path() = joinpath(pkgdir(@__MODULE__), "data", "wilber2016", "converted_temp_data.csv")

"""
    load_wilber2016() -> Vector{NamedTuple}

The 589 swab-to-swab transitions of the 20 adult frogs of the Wilber et al. (2016) experiment
([`wilber2016_path`](@ref)), one named tuple per transition:

- `individual` (e.g. `"F20_3"`), `temperature` (°C: 4, 12, 20 or 26);
- `interval`: days between the swabs (3, or 8 for the first and 120 for a few gaps);
- `load`, `load_next`: qPCR loads in zoospore equivalents per swab (natural scale, 0 when
  uninfected; `load_next` is `NaN` when the frog died);
- `survived`: whether the frog was alive at the second swab.
"""
function load_wilber2016()
    lines = filter(!isempty ∘ strip, readlines(wilber2016_path()))
    header = split(lines[1], ',')
    header == ["size", "sizeNext", "time", "individual", "temp", "surv", "stage"] ||
        throw(ArgumentError("unexpected columns in $(wilber2016_path())"))
    return map(lines[2:end]) do l
        v = split(l, ',')
        (individual = String(v[4]), temperature = parse(Int, v[5]),
         interval = parse(Int, v[3]), load = parse(Float64, v[1]),
         load_next = isempty(v[2]) ? NaN : parse(Float64, v[2]), survived = v[6] == "1")
    end
end
