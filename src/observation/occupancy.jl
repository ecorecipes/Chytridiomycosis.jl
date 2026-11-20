# Occupancy layer: likelihoods for detection/non-detection data at sites surveyed repeatedly.
# Model-agnostic; users are Adams et al. (2010) (single-season occupancy of Bd in host
# populations, src/models/adams2010.jl) and Mosher et al. (2018) (dynamic two-species
# host–pathogen occupancy, src/models/mosher2018.jl, src/fitting/mosher2018.jl).
#
# Single season (MacKenzie et al. 2002, Ecology 83:2248). A site is occupied with probability
# ψ and stays so during the season; each survey j detects the species with probability p_j if
# the site is occupied and never if not (no false positives). For one site's history y:
#   L = ψ ∏_j p_j^{y_j} (1 − p_j)^{1 − y_j} + (1 − ψ) 1{no detection}.
#
# Dynamic multistate occupancy (MacKenzie et al. 2009, Ecology 90:823; robust design). A site
# is in one of K states in each season (primary occasion); states change between seasons by a
# row-stochastic transition matrix, and stay fixed over the surveys (secondary occasions)
# within a season. Each survey records one observation column; the likelihood is the HMM
# forward algorithm over seasons with the product of survey emission probabilities as the
# season's emission. This reuses hmm_loglikelihood (observation/hmm.jl) and
# robust_design_likelihoods (observation/multistate.jl).
#
# Two-species conditional parameterisation (Richmond, Hines & Beissinger 2010, Ecol Appl
# 20:2036; dynamic form as in Miller et al. 2012 and Mosher 2017, PhD dissertation, Colorado
# State University, Tables 1.1-1.2). States in the order (A, B, U, AB): species A only,
# species B only, neither, both. Species A is the host and B the pathogen in Mosher et al.
# (2018). Initial state:
#   φ0 = [ψA(1−ψBA), (1−ψA)ψBa, (1−ψA)(1−ψBa), ψA ψBA]
# Transitions (rows: state at t, columns: state at t + 1), each species' colonisation (γ) and
# extinction (ε) conditional on whether the other species is present at t:
#   A  → [(1−εA)(1−γBA),  εA γBA,          εA(1−γBA),       (1−εA)γBA      ]
#   B  → [γAB εB,          (1−γAB)(1−εB),   (1−γAB)εB,       γAB(1−εB)      ]
#   U  → [γA(1−γB),        (1−γA)γB,        (1−γA)(1−γB),    γA γB          ]
#   AB → [(1−εAB)εBA,      εAB(1−εBA),      εAB εBA,         (1−εAB)(1−εBA) ]
# Detection (rows: true state, columns: observed A, B, U, AB in one survey):
#   A  → [pA, 0, 1−pA, 0];  B → [0, pB, 1−pB, 0];  U → [0, 0, 1, 0]
#   AB → [rA(1−rBA), (1−rA)rBa, (1−rA)(1−rBa), rA rBA]
# with pB = rBa = 0 when the pathogen is detected only on swabs of detected hosts. Surveys
# in which one species was not recorded (a toad seen but not swabbed) are coded with the
# partial columns 5-8 below, whose probabilities are sums of the full columns.

# ---------------------------------------------------------------------------------------
# Single season

"""
    occupancy_loglikelihood(y::AbstractVector, ψ, p) -> Real
    occupancy_loglikelihood(Y::AbstractMatrix, ψ, p) -> Real

Single-season occupancy log-likelihood (MacKenzie et al. 2002). `y` is one site's detection
history (`1` detected, `0` not, `missing` no survey); `ψ` is the occupancy probability and
`p` the detection probability, one number or one per survey. For a matrix `Y` (rows are
sites) the site log-likelihoods are summed; `ψ` may then be a number or a vector over sites
and `p` a number, a vector over sites or a matrix the size of `Y`. Assumes no false
positives and closure within the season.
"""
function occupancy_loglikelihood(y::AbstractVector, ψ::Real, p)
    lp = zero(promote_type(typeof(ψ), _occ_eltype(p)))
    detected = false
    for j in eachindex(y)
        yj = y[j]
        yj === missing && continue
        pj = p isa Real ? p : p[j]
        if yj == 1
            detected = true
            lp += log(pj)
        elseif yj == 0
            lp += log1p(-pj)
        else
            throw(ArgumentError("detections must be 0, 1 or missing; got $yj"))
        end
    end
    return detected ? log(ψ) + lp : log(ψ * exp(lp) + (1 - ψ))
end

_occ_eltype(p::Real) = typeof(p)
_occ_eltype(p) = eltype(p)

function occupancy_loglikelihood(Y::AbstractMatrix, ψ, p)
    s = zero(promote_type(_occ_eltype(ψ), _occ_eltype(p)))
    for i in axes(Y, 1)
        ψi = ψ isa Real ? ψ : ψ[i]
        pi = p isa Real ? p : (p isa AbstractMatrix ? view(p, i, :) : p[i])
        s += occupancy_loglikelihood(view(Y, i, :), ψi, pi)
    end
    return s
end

_logit(x) = log(x / (1 - x))
_expit(x) = 1 / (1 + exp(-x))

"""
    occupancy_mle(Y) -> NamedTuple

Maximum likelihood estimates of a single-season occupancy model with constant `ψ` and `p`
from a detection matrix `Y` (rows sites, columns surveys, `missing` allowed). Fitted by BFGS
on the logit scale. Returns `ψ`, `p`, their standard errors `se_ψ`, `se_p` (delta method
from the observed information), `loglik` and `naive`, the fraction of sites with at least one
detection.
"""
function occupancy_mle(Y::AbstractMatrix)
    nll(x) = -occupancy_loglikelihood(Y, _expit(x[1]), _expit(x[2]))
    naive = count(i -> any(==(1), skipmissing(view(Y, i, :))), axes(Y, 1)) / size(Y, 1)
    r = _ms_bfgs(nll, [_logit(clamp(naive, 0.05, 0.95)), 0.0]; gtol = 1e-8)
    V, _ = _ms_vcov(-_ms_hessian(nll, r.x; h = 1e-4))
    ψ, p = _expit(r.x[1]), _expit(r.x[2])
    return (ψ = ψ, p = p, se_ψ = ψ * (1 - ψ) * sqrt(V[1, 1]),
            se_p = p * (1 - p) * sqrt(V[2, 2]), loglik = -r.f, naive = naive)
end

# ---------------------------------------------------------------------------------------
# Dynamic multistate occupancy

"""
    dynamic_occupancy_loglikelihood(initial, transitions, emission, codes) -> Real

Log-likelihood of one site's history under a dynamic multistate occupancy model (robust
design). `initial` is the state distribution in the first season, `transitions` one `K×K`
row-stochastic matrix or a vector of them (season `t` to `t + 1`), `emission` a `K×O`
matrix of observation probabilities (or a function `(t, j) -> matrix` for survey `j` of
season `t`), and `codes` a `T×J` matrix of observed column indices with `missing` for no
survey. Computed with [`robust_design_likelihoods`](@ref) and [`hmm_loglikelihood`](@ref).
"""
dynamic_occupancy_loglikelihood(initial, transitions, emission, codes::AbstractMatrix) =
    hmm_loglikelihood(initial, transitions, robust_design_likelihoods(emission, codes))

"""
    simulate_dynamic_occupancy(rng, initial, transitions, emission, nseasons, nsurveys;
                               observed = nothing) -> NamedTuple

Simulate one site of a dynamic multistate occupancy model. Returns `states` (the latent
state index in each season) and `codes` (`nseasons×nsurveys` observation columns). Survey
observations are drawn from the columns `observed` (default: all) of `emission` (a matrix or a function
`(t, j) -> matrix`), whose rows must sum to one over those columns.
"""
function simulate_dynamic_occupancy(rng::AbstractRNG, initial::AbstractVector, transitions,
                                    emission, nseasons::Integer, nsurveys::Integer;
                                    observed = nothing)
    states = zeros(Int, nseasons)
    codes = Matrix{Union{Missing,Int}}(undef, nseasons, nsurveys)
    states[1] = _occ_draw(rng, initial)
    for t in 2:nseasons
        states[t] = _occ_draw(rng, view(_hmm_P(transitions, t - 1), states[t-1], :))
    end
    for t in 1:nseasons, j in 1:nsurveys
        E = _rd_emission(emission, t, j)
        cols = observed === nothing ? axes(E, 2) : observed
        w = [E[states[t], c] for c in cols]
        codes[t, j] = collect(cols)[_occ_draw(rng, w)]
    end
    return (states = states, codes = codes)
end

function _occ_draw(rng::AbstractRNG, w)
    u = rand(rng) * sum(w)
    c = zero(eltype(w))
    for (k, wk) in enumerate(w)
        c += wk
        u <= c && return k
    end
    return length(w)
end

# ---------------------------------------------------------------------------------------
# Two-species conditional parameterisation

const TWO_SPECIES_STATES = (:A, :B, :U, :AB)

"""
    two_species_initial(ψA, ψBA, ψBa) -> Vector

Initial state distribution of the two-species conditional occupancy model (Richmond et al.
2010) in the state order `(A, B, U, AB)`: `ψA` is the occupancy probability of species A and
`ψBA`, `ψBa` that of species B given A present or absent.
"""
two_species_initial(ψA, ψBA, ψBa) =
    [ψA * (1 - ψBA), (1 - ψA) * ψBa, (1 - ψA) * (1 - ψBa), ψA * ψBA]

"""
    two_species_transition(; γA, γAB, γB, γBA, εA, εAB, εB, εBA) -> Matrix

The 4×4 season-to-season transition matrix of the dynamic two-species occupancy model,
states `(A, B, U, AB)`, rows the state at `t` (Mosher 2017, Table 1.1). Colonisation `γ` and
extinction `ε` of each species are conditional on the other: `γA`, `εA` for species A when
B is absent and `γAB`, `εAB` when B is present; `γB`, `εB` for B when A is absent and `γBA`,
`εBA` when A is present.
"""
function two_species_transition(; γA, γAB, γB, γBA, εA, εAB, εB, εBA)
    return [(1-εA)*(1-γBA)  εA*γBA          εA*(1-γBA)        (1-εA)*γBA
            γAB*εB          (1-γAB)*(1-εB)  (1-γAB)*εB        γAB*(1-εB)
            γA*(1-γB)       (1-γA)*γB       (1-γA)*(1-γB)     γA*γB
            (1-εAB)*εBA     εAB*(1-εBA)     εAB*εBA           (1-εAB)*(1-εBA)]
end

"""
    two_species_detection(; pA, rA, rBA, pB = 0, rBa = 0) -> Matrix

The 4×8 per-survey observation matrix of the two-species occupancy model (rows the true
state `(A, B, U, AB)`). Columns 1-4 are the observed states `A`, `B`, `U`, `AB` (Mosher 2017,
Table 1.2): `pA`, `pB` detect a species alone, `rA` detects A with B present, and `rBA`,
`rBa` detect B with A present and detected or not detected. Columns 5-8 are partial records,
the sums of the matching full columns: A detected and B not recorded, A not detected and B
not recorded, B detected and A not recorded, B not detected and A not recorded. Use
[`two_species_code`](@ref) to map a survey's records to a column.
"""
function two_species_detection(; pA, rA, rBA, pB = 0.0, rBa = 0.0)
    T = promote_type(typeof(pA), typeof(rA), typeof(rBA), typeof(pB), typeof(rBa))
    D = zeros(T, 4, 8)
    D[1, 1], D[1, 3] = pA, 1 - pA
    D[2, 2], D[2, 3] = pB, 1 - pB
    D[3, 3] = 1
    D[4, 1] = rA * (1 - rBA)
    D[4, 2] = (1 - rA) * rBa
    D[4, 3] = (1 - rA) * (1 - rBa)
    D[4, 4] = rA * rBA
    for k in 1:4
        D[k, 5] = D[k, 1] + D[k, 4]     # A = 1, B unrecorded
        D[k, 6] = D[k, 2] + D[k, 3]     # A = 0, B unrecorded
        D[k, 7] = D[k, 2] + D[k, 4]     # B = 1, A unrecorded
        D[k, 8] = D[k, 1] + D[k, 3]     # B = 0, A unrecorded
    end
    return D
end

"""
    two_species_code(a, b) -> Union{Int,Missing}

Observation column of [`two_species_detection`](@ref) for one survey with records `a` (species
A) and `b` (species B), each `1`, `0` or `missing`. Both missing (no survey) gives `missing`.
"""
function two_species_code(a, b)
    if a === missing
        b === missing && return missing
        return b == 1 ? 7 : 8
    elseif b === missing
        return a == 1 ? 5 : 6
    end
    a == 1 && b == 1 && return 4
    a == 1 && return 1
    b == 1 && return 2
    return 3
end

"""
    two_species_records(code) -> Tuple

Inverse of [`two_species_code`](@ref): the records `(a, b)` of an observation column.
"""
function two_species_records(code)
    code === missing && return (missing, missing)
    return ((1, 0), (0, 1), (0, 0), (1, 1), (1, missing), (0, missing), (missing, 1),
            (missing, 0))[code]
end

# ---------------------------------------------------------------------------------------
# Single season with covariates (logit links), as fitted in Program MARK or PRESENCE

"""
    SingleSeasonOccupancyModel(Y, Xψ, Xp; ψnames, pnames)

Single-season occupancy model with logit-linear covariates: `Y` is the `sites × surveys`
detection matrix (`1`, `0`, `missing`), `Xψ` the `sites × kψ` design matrix of occupancy
(include a column of ones for the intercept) and `Xp` the design of detection, either
`sites × kp` (site covariates) or `sites × surveys × kp` (survey or animal covariates, as
in Adams et al. 2010 where each survey is one animal tested). Coefficients are ordered
`[βψ; βp]`. Implements the `LogDensityProblems` interface (the log-likelihood) and
[`parameter_names`](@ref); fit with [`maximum_likelihood`](@ref).
"""
struct SingleSeasonOccupancyModel{TY<:AbstractMatrix,TA<:AbstractMatrix,TB<:AbstractArray}
    Y::TY
    Xψ::TA
    Xp::TB
    names::Vector{String}
end

function SingleSeasonOccupancyModel(Y::AbstractMatrix, Xψ::AbstractMatrix, Xp::AbstractArray;
        ψnames = ["ψ:" * string(k) for k in 1:size(Xψ, 2)],
        pnames = ["p:" * string(k) for k in 1:size(Xp, ndims(Xp))])
    size(Xψ, 1) == size(Y, 1) || throw(DimensionMismatch("Xψ needs one row per site"))
    if ndims(Xp) == 2
        size(Xp, 1) == size(Y, 1) || throw(DimensionMismatch("Xp needs one row per site"))
    elseif ndims(Xp) == 3
        size(Xp)[1:2] == size(Y) || throw(DimensionMismatch("Xp must be sites × surveys × k"))
    else
        throw(ArgumentError("Xp must be a matrix or a 3-d array"))
    end
    return SingleSeasonOccupancyModel(Y, Xψ, Xp, vcat(string.(ψnames), string.(pnames)))
end

parameter_names(m::SingleSeasonOccupancyModel) = m.names
LogDensityProblems.dimension(m::SingleSeasonOccupancyModel) = length(m.names)
LogDensityProblems.capabilities(::Type{<:SingleSeasonOccupancyModel}) =
    LogDensityProblems.LogDensityOrder{0}()

function LogDensityProblems.logdensity(m::SingleSeasonOccupancyModel, θ::AbstractVector)
    kψ = size(m.Xψ, 2)
    βψ, βp = view(θ, 1:kψ), view(θ, kψ+1:length(θ))
    s = zero(eltype(θ))
    J = size(m.Y, 2)
    for i in axes(m.Y, 1)
        ψ = _expit(dot(view(m.Xψ, i, :), βψ))
        p = ndims(m.Xp) == 2 ? _expit(dot(view(m.Xp, i, :), βp)) :
            [_expit(dot(view(m.Xp, i, j, :), βp)) for j in 1:J]
        s += occupancy_loglikelihood(view(m.Y, i, :), ψ, p)
    end
    return s
end

"""
    maximum_likelihood(m::SingleSeasonOccupancyModel; x0 = zeros(dim)) -> NamedTuple

Maximum likelihood fit by BFGS. Returns `θ`, `names`, `se` (Wald), `loglik`, `deviance`
(`−2 loglik`), `k`, `n` (sites, the AICc sample size) and `aicc`.
"""
function maximum_likelihood(m::SingleSeasonOccupancyModel; x0 = nothing, iters::Integer = 2000,
                            gtol::Real = 1e-6)
    k = length(m.names)
    nll(θ) = -LogDensityProblems.logdensity(m, θ)
    r = _ms_bfgs(nll, x0 === nothing ? zeros(k) : collect(float.(x0)); iters, gtol)
    V, _ = _ms_vcov(-_ms_hessian(nll, r.x; h = 1e-4))
    n = size(m.Y, 1)
    ll = -r.f
    return (θ = r.x, names = m.names, se = [v >= 0 ? sqrt(v) : NaN for v in LinearAlgebra.diag(V)],
            loglik = ll, deviance = -2ll, k = k, n = n,
            aicc = -2ll + 2k + 2k * (k + 1) / (n - k - 1))
end
