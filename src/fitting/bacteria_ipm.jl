# NEW MODEL FIT (not a published analysis): the two-trait Bd–bacteria IPM
# (src/novel/bacteria_ipm.jl) fitted to the joint Bd and Janthinobacterium lividum loads of the
# treated cohort of Knapp et al. (2022), Dusy Basin 2012 (data/knapp2022, CC BY 4.0).
#
# Frogs were group-marked, so a frog's two loads cannot be followed from visit to visit; what
# the swabs give is the joint distribution of (Bd, J. lividum) in the treated group on days 0,
# 12, 37 and 56. The model is fitted to those cross-sections:
#   - start: the day-0 swabs (after itraconazole, before the baths), Bd loads binned on the
#     mesh, J. lividum absent; a bath of every frog on day 0 (the day-1 bath of the trial falls in
#     the same 3-day step and is not added again);
#   - infection: a constant probability φ per step from the pond (no zoospore data);
#   - temperature: fixed at `temperature` (15 °C by default; the trial has no water
#     temperatures in the deposit), which with the Wilber et al. (2016) growth regression sets
#     where Bd loads settle;
#   - units: Knapp et al. report ITS copies per swab and the Wilber kernel uses zoospore
#     equivalents, so ln ZE = ln copies − c with c, the log ITS copies per zoospore, fitted
#     within the range of strains in Longo et al. (2013) (compendium longo2013: 10–144);
#   - observations on days 12, 37 and 56 (model times 12, 36 and 57 days, the nearest steps):
#     each swab with both loads measured contributes the model's probability of its class (Bd
#     and J. lividum each present or absent; J. lividum below 1 GE is absent) times the
#     density of its loads per ln unit, among live frogs. A small floor (1e-10) keeps the log
#     finite for loads the model puts no mass on.
# Two versions: `coupled = false` (bacteria do not act on Bd; κ_growth = 0) and `coupled =
# true` (κ_growth and the threshold y50 fitted). Survival effects are held at 0: live
# cross-sections carry little information about survival. Comparing the two asks whether frogs
# with more J. lividum carry less Bd than the uncoupled model predicts.

const _KJ_NAMES_BASE = (:decline, :bath_mean, :bath_sd, :σy, :φ, :copies_per_ze)
const _KJ_NAMES_COUPLED = (_KJ_NAMES_BASE..., :κ_growth, :y50)

"""
    KnappJointModel(; coupled = true, temperature = 15.0, swabs = load_knapp2022(),
                    xmesh = IPMMesh(-5.0, 18.0, 30), ymesh = IPMMesh(0.0, 14.0, 20),
                    base = BacteriaIPMParameters(κ_survival = 0.0, κ_growth = 0.0))

The two-trait Bd–bacteria IPM ([`BacteriaIPMProblem`](@ref)) as a likelihood for the joint Bd
and *J. lividum* loads of the treated cohort of Knapp et al. (2022): a **new fit**, not a
published analysis. The cohort starts from its day-0 swabs, is bathed on day 0, is infected
with a constant probability `φ` per 3-day step, and is compared with the swabs of days 12, 37
and 56 (classes and joint densities among live frogs).

Parameters (natural scale; [`parameter_names`](@ref)): `decline` (per day), `bath_mean`,
`bath_sd`, `σy`, `φ`, `copies_per_ze` (ITS copies per zoospore equivalent, between the lowest
and highest of Longo et al. 2013) and, if `coupled`, `κ_growth` and `y50`. Other values are
those of `base`. The log density takes the unconstrained vector `x` (log, logit and
log-ratio transforms); [`maximum_likelihood`](@ref) fits it.
"""
struct KnappJointModel
    swabs::Vector{KnappSwab}
    temperature::Float64
    xmesh::IPMMesh
    ymesh::IPMMesh
    coupled::Bool
    base::BacteriaIPMParameters
    cz_range::Tuple{Float64,Float64}
end

function KnappJointModel(; coupled::Bool = true, temperature::Real = 15.0,
        swabs::AbstractVector{KnappSwab} = load_knapp2022(),
        xmesh::IPMMesh = IPMMesh(-5.0, 18.0, 30), ymesh::IPMMesh = IPMMesh(0.0, 14.0, 20),
        base::BacteriaIPMParameters = BacteriaIPMParameters(κ_survival = 0.0, κ_growth = 0.0))
    tr = filter(s -> s.group === :treated && !isnan(s.bd), swabs)
    any(s -> s.day == 0, tr) || throw(ArgumentError("no day-0 swabs of treated frogs"))
    q = parameter_set(:longo2013)
    return KnappJointModel(collect(tr), Float64(temperature), xmesh, ymesh, coupled, base,
        (Float64(q.copies_per_zoospore_min), Float64(q.copies_per_zoospore_max)))
end

parameter_names(m::KnappJointModel) = m.coupled ? _KJ_NAMES_COUPLED : _KJ_NAMES_BASE

LogDensityProblems.dimension(m::KnappJointModel) = length(parameter_names(m))
LogDensityProblems.capabilities(::Type{KnappJointModel}) = LogDensityProblems.LogDensityOrder{0}()

"""
    knapp_joint_params(m::KnappJointModel, x) -> NamedTuple

The natural-scale parameters of [`KnappJointModel`](@ref) `m` for the unconstrained vector `x`.
"""
function knapp_joint_params(m::KnappJointModel, x::AbstractVector)
    length(x) == LogDensityProblems.dimension(m) ||
        throw(DimensionMismatch("x must have length $(LogDensityProblems.dimension(m))"))
    lo, hi = log.(m.cz_range)
    vals = (exp(x[1]), x[2], exp(x[3]), exp(x[4]), _logistic(x[5]),
        exp(lo + (hi - lo) * _logistic(x[6])))
    m.coupled && (vals = (vals..., x[7], x[8]))
    return NamedTuple{parameter_names(m)}(vals)
end

"""
    knapp_joint_start(m::KnappJointModel) -> Vector{Float64}

Starting point (unconstrained) for [`maximum_likelihood`](@ref): the decline and bath values
of [`BacteriaIPMParameters`](@ref), `σy = 0.5`, `φ = 0.1`, the geometric middle of the
copy-number range and, if coupled, `κ_growth = 0.5`, `y50 = ln 100`.
"""
function knapp_joint_start(m::KnappJointModel)
    b = m.base
    x = [log(b.decline), b.bath_mean, log(b.bath_sd), log(b.σy), log(0.1 / 0.9), 0.0]
    m.coupled && append!(x, [0.5, log(100.0)])
    return x
end

# Parameters and solution for x.
function _knapp_joint_solve(m::KnappJointModel, x::AbstractVector)
    θ = knapp_joint_params(m, x)
    p = with(m.base; decline = θ.decline, bath_mean = θ.bath_mean, bath_sd = θ.bath_sd,
        σy = θ.σy, κ_growth = m.coupled ? θ.κ_growth : 0.0,
        y50 = m.coupled ? θ.y50 : m.base.y50)
    mesh = IPMMesh2D(m.xmesh, m.ymesh)
    R = _bipm_ranges(mesh)
    u0 = zeros(R.Z)
    c = log(θ.copies_per_ze)
    for s in m.swabs
        s.day == 0 || continue
        if s.bd > 0
            u0[R.I0[_nearest_bin(m.xmesh, log(s.bd) - c)]] += 1
        else
            u0[R.S0] += 1
        end
    end
    u0 ./= sum(u0)
    prob = BacteriaIPMProblem(p; temperature = m.temperature, xmesh = m.xmesh,
        ymesh = m.ymesh, u0, days = 57, baths = [(day = 0, coverage = 1.0)],
        transmission = :constant, φ = θ.φ)
    return θ, solve(prob)
end

_nearest_bin(mesh::IPMMesh, v) = clamp(floor(Int, (v - mesh.lower) / meshwidth(mesh)) + 1, 1, mesh.n)

# The model time index for a swab day (steps of 3 days, nearest step).
_knapp_joint_index(day) = round(Int, day / 3) + 1

"""
    pointwise_loglikelihood(m::KnappJointModel, x) -> Vector{Float64}

Log-likelihood of each swab of days 12, 37 and 56 with both loads measured, in the order of
`m.swabs`.
"""
function pointwise_loglikelihood(m::KnappJointModel, x::AbstractVector)
    θ, sol = _knapp_joint_solve(m, x)
    c = log(θ.copies_per_ze)
    hx, hy = meshwidth(m.xmesh), meshwidth(m.ymesh)
    out = Float64[]
    for s in m.swabs
        (s.day > 0 && !isnan(s.jl)) || continue
        d = joint_load_distribution(sol, _knapp_joint_index(s.day))
        alive = d.uninfected_without + sum(d.uninfected) + sum(d.infected_without) +
                sum(d.infected)
        bdpos, jlpos = s.bd > 0, s.jl >= 1
        ix = bdpos ? _nearest_bin(m.xmesh, log(s.bd) - c) : 0
        iy = jlpos ? _nearest_bin(m.ymesh, log(s.jl)) : 0
        P = !bdpos && !jlpos ? d.uninfected_without :
            bdpos && !jlpos ? d.infected_without[ix] / hx :
            !bdpos && jlpos ? d.uninfected[iy] / hy : d.infected[ix, iy] / (hx * hy)
        push!(out, log(max(P / alive, 0.0) + 1e-10))
    end
    return out
end

function LogDensityProblems.logdensity(m::KnappJointModel, x::AbstractVector)
    all(isfinite, x) || return -Inf
    v = try
        sum(pointwise_loglikelihood(m, x))
    catch e
        e isa ArgumentError || e isa DomainError || rethrow()
        -Inf
    end
    return isfinite(v) ? v : -Inf
end

"""
    maximum_likelihood(m::KnappJointModel; starts = 3, iters = 1500, spread = 0.5,
                       x0 = knapp_joint_start(m), rng = Random.Xoshiro(1))

Maximum likelihood fit of [`KnappJointModel`](@ref) by Nelder–Mead on the unconstrained
parameters, from `x0` and `starts − 1` points around it, each polished by a smaller simplex.
Returns `model`, `x`, `estimates` (natural scale), `loglik`, `aic`, `n` (swabs) and `logliks`
(best of each start).
"""
function maximum_likelihood(m::KnappJointModel; starts::Integer = 3, iters::Integer = 1500,
        spread::Real = 0.5, x0::AbstractVector = knapp_joint_start(m),
        rng::AbstractRNG = Random.Xoshiro(1))
    f(x) = -LogDensityProblems.logdensity(m, x)
    best, bestf = copy(x0), Inf
    logliks = Float64[]
    for s in 1:starts
        x = s == 1 ? copy(x0) : x0 .+ spread .* randn(rng, length(x0))
        x, fx = nelder_mead(f, x; iters)
        x, fx = nelder_mead(f, x; iters = iters ÷ 2, step = 0.2)
        push!(logliks, -fx)
        fx < bestf && ((best, bestf) = (x, fx))
    end
    k = length(best)
    return (model = m, x = best, estimates = knapp_joint_params(m, best), loglik = -bestf,
        aic = 2k + 2bestf, n = length(pointwise_loglikelihood(m, best)), logliks = logliks)
end
