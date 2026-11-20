# Pradel (1996) temporal-symmetry capture–recapture model, in the robust-design form of
# Program MARK ("Robust Design Pradel" with Huggins closed captures and p = c).
#
# Pradel R (1996) Utilization of capture-mark-recapture for the study of recruitment and
# population growth rate. Biometrics 52:703-709. Williams, Nichols & Conroy (2002) Analysis and
# Management of Animal Populations, section 18.4. Hardy et al. (2023, Appendix S1.2) restate the
# model in the f (per-capita recruitment) parameterization used here.
#
# Primary periods t = 1…T, each with one or more secondary occasions j. Parameters:
#   φ_t  apparent survival from t to t + 1                (t = 1…T−1)
#   f_t  per-capita recruitment from t to t + 1           (t = 1…T−1)
#   p_tj capture probability at secondary occasion j of t (p = c: no behavioural effect)
# Derived: λ_t = φ_t + f_t (eq. 2 of Hardy et al. Appendix S1.2), seniority
# γ_{t+1} = φ_t / λ_t (their eq. 3), p*_t = 1 − ∏_j (1 − p_tj), and N_t = n_t / p*_t.
#
# Likelihood. The history of an animal first caught in period e and last caught in period l
# has probability (conditional on the number of animals ever caught)
#
#   P(h) = N_e ξ_e  ∏_{t=e}^{l−1} φ_t  ∏_{t=e}^{l} c_t(h)  χ_l  /  Σ_t N_t ξ_t p*_t
#
# with N_t ∝ ∏_{s<t} λ_s, c_t(h) = ∏_j p_tj^{x_tj} (1 − p_tj)^{1−x_tj} when the animal was caught
# in period t and 1 − p*_t when not, χ_T = 1, χ_t = 1 − φ_t + φ_t (1 − p*_{t+1}) χ_{t+1} (never
# seen after t) and ξ_1 = 1, ξ_t = 1 − γ_t + γ_t (1 − p*_{t−1}) ξ_{t−1} (never seen before t).
# For a caught period, MARK's robust design factors c_t into p*_t times the Huggins conditional
# probability of the within-period pattern, ∏_j p^x (1 − p)^{1−x} / p*_t; the product is the
# same. The denominator is the expected number of animals ever caught, so the probabilities
# of all histories with at least one capture sum to one (checked in test/hardy2023.jl).
#
# The log-likelihood depends on the histories only through counts (captures and misses per
# occasion, first and last captures per period, animals known alive per interval), which are
# computed once; each evaluation is O(number of occasions). Unsampled occasions (MARK's ".")
# contribute nothing and are left out of p*. The sampling pattern must be the same for all
# animals of a group, as in the data this layer was written for.

"""
    PradelGroup(histories; primary = 1:size(histories, 2))

Capture histories of one group (for example one site) for the Pradel (1996)
temporal-symmetry model.

- `histories` is an `n×J` matrix with one row per animal and one column per secondary
  occasion: `1` caught, `0` not caught, `missing` (or `-1`) for an occasion that was not
  sampled. The sampling pattern must be the same in every row.
- `primary[j]` is the primary period (`1…T`) of column `j`. The default, one occasion per
  period, gives the standard Pradel model; several columns per period give the robust design
  with Huggins closed captures and `p = c`, as in Program MARK.

Rows with no capture are ignored. The constructor reduces the histories to the counts that
the likelihood needs; see [`pradel_loglikelihood`](@ref).
"""
struct PradelGroup
    T::Int
    primary::Vector{Int}
    sampled::Vector{Bool}
    captures::Vector{Int}   # per occasion: animals caught
    zeros::Vector{Int}      # per occasion: animals not caught in a period in which they were
    missed::Vector{Int}     # per period: animals not caught in it, between first and last capture
    alive::Vector{Int}      # per interval t → t + 1: animals with first ≤ t < last
    first::Vector{Int}      # per period: animals first caught in it
    last::Vector{Int}       # per period: animals last caught in it
    caught::Vector{Int}     # per period: animals caught at least once in it (n_t)
    n::Int                  # animals ever caught
end

_pradel_code(x) = x === missing ? -1 : Int(x)

function PradelGroup(histories::AbstractMatrix; primary::AbstractVector{<:Integer} =
                     1:size(histories, 2))
    nrow, J = size(histories)
    length(primary) == J || throw(DimensionMismatch("primary must have one entry per column"))
    issorted(primary) && first(primary) == 1 && all(diff(primary) .<= 1) ||
        throw(ArgumentError("primary must run 1, …, T in order without gaps"))
    T = last(primary)
    T >= 2 || throw(ArgumentError("need at least two primary periods"))
    H = [_pradel_code(histories[i, j]) for i in 1:nrow, j in 1:J]
    all(in((-1, 0, 1)), H) || throw(ArgumentError("histories must be 0, 1 or missing"))
    sampled = nrow == 0 ? fill(true, J) : H[1, :] .!= -1
    all(i -> all((H[i, :] .!= -1) .== sampled), 1:nrow) ||
        throw(ArgumentError("the sampling pattern (missing occasions) must be the same in every row"))
    captures = zeros(Int, J); zs = zeros(Int, J)
    missed = zeros(Int, T); alive = zeros(Int, T - 1)
    fst = zeros(Int, T); lst = zeros(Int, T); caught = zeros(Int, T)
    n = 0
    cols = [findall(==(t), primary) for t in 1:T]
    for i in 1:nrow
        inper = [any(j -> H[i, j] == 1, cols[t]) for t in 1:T]
        e = findfirst(inper)
        e === nothing && continue
        l = findlast(inper)
        n += 1
        fst[e] += 1; lst[l] += 1
        for t in e:l
            if inper[t]
                caught[t] += 1
                for j in cols[t]
                    H[i, j] == 1 && (captures[j] += 1)
                    H[i, j] == 0 && (zs[j] += 1)
                end
            else
                missed[t] += 1
            end
            t < l && (alive[t] += 1)
        end
    end
    return PradelGroup(T, collect(Int, primary), sampled, captures, zs, missed, alive, fst, lst,
                       caught, n)
end

"""
    pradel_capture_probability(g::PradelGroup, p) -> Vector

The probability `p*_t = 1 − ∏_j (1 − p_tj)` of catching an animal present in period `t` at
least once, over the sampled occasions of each period (zero for a period that was not
sampled). `p` has one entry per occasion (column of the histories).
"""
function pradel_capture_probability(g::PradelGroup, p::AbstractVector)
    length(p) == length(g.primary) || throw(DimensionMismatch("p needs one entry per occasion"))
    R = eltype(p)
    q = ones(R, g.T)
    for j in eachindex(p)
        g.sampled[j] && (q[g.primary[j]] *= 1 - p[j])
    end
    return 1 .- q
end

"""
    pradel_loglikelihood(g::PradelGroup, φ, f, p) -> Real

Log-likelihood of the Pradel (1996) temporal-symmetry model with apparent survival `φ` and
per-capita recruitment `f` (both of length `T − 1`, for the intervals between primary
periods) and capture probabilities `p` (one per secondary occasion), conditional on the
number of animals ever caught. With several occasions per period it is the robust-design
Pradel model of Program MARK (Huggins closed captures, `p = c`). The realized population
growth rate is `λ = φ + f`.

On the boreal toad data of Hardy et al. (2023) it reproduces MARK's −2 log L (for example
29334.27 for the constant model; see `test/hardy2023.jl`).
"""
function pradel_loglikelihood(g::PradelGroup, φ::AbstractVector, f::AbstractVector,
                              p::AbstractVector)
    T = g.T
    length(φ) == T - 1 && length(f) == T - 1 ||
        throw(DimensionMismatch("φ and f need T − 1 = $(T - 1) entries"))
    g.n == 0 && return zero(promote_type(eltype(φ), eltype(f), eltype(p)))
    ps = pradel_capture_probability(g, p)
    R = promote_type(eltype(φ), eltype(f), eltype(p), Float64)
    χ = ones(R, T)
    for t in T-1:-1:1
        χ[t] = 1 - φ[t] + φ[t] * (1 - ps[t+1]) * χ[t+1]
    end
    ξ = ones(R, T)
    logN = zeros(R, T)
    for t in 2:T
        λ = φ[t-1] + f[t-1]
        γ = φ[t-1] / λ
        ξ[t] = 1 - γ + γ * (1 - ps[t-1]) * ξ[t-1]
        logN[t] = logN[t-1] + log(λ)
    end
    m = maximum(logN)
    D = sum(exp(logN[t] - m) * ξ[t] * ps[t] for t in 1:T)
    ll = -g.n * (m + log(D))
    _xlog(k, x) = k == 0 ? zero(R) : k * log(x)
    for t in 1:T
        ll += _xlog(g.first[t], ξ[t]) + g.first[t] * logN[t]
        ll += _xlog(g.last[t], χ[t]) + _xlog(g.missed[t], 1 - ps[t])
        t < T && (ll += _xlog(g.alive[t], φ[t]))
    end
    for j in eachindex(p)
        g.sampled[j] || continue
        ll += _xlog(g.captures[j], p[j]) + _xlog(g.zeros[j], 1 - p[j])
    end
    return ll
end

"""
    pradel_abundance(g::PradelGroup, p) -> Vector

Abundance in each primary period estimated as `N_t = n_t / p*_t`, where `n_t` is the number of
animals caught in period `t` (Huggins estimator, as derived by MARK's robust-design Pradel
model). `NaN` for a period that was not sampled.
"""
function pradel_abundance(g::PradelGroup, p::AbstractVector)
    ps = pradel_capture_probability(g, p)
    return [ps[t] > 0 ? g.caught[t] / ps[t] : NaN for t in 1:g.T]
end

"""
    PradelProblem(groups, designs; names = nothing)

Maximum-likelihood problem for a multi-group Pradel model with linear predictors.

- `groups` is a vector of [`PradelGroup`](@ref)s.
- `designs[g]` is a named tuple `(φ = Xφ, f = Xf, p = Xp)` of design matrices for group `g`:
  `Xφ` and `Xf` have `T − 1` rows (intervals), `Xp` one row per occasion. All groups share the
  coefficient vector `θ = [βφ; βf; βp]` (the matrices of one parameter must have the same
  number of columns in every group).
- Links: `logit φ = Xφ βφ`, `log f = Xf βf`, `logit p = Xp βp`, as in Hardy et al. (2023,
  Appendix S1.2). Rows of `Xp` for unsampled occasions are ignored.

`names` optionally labels the coefficients (default `φ_1, …, f_1, …, p_1, …`). The problem
implements the LogDensityProblems interface (log density = log-likelihood, flat prior on `θ`).
`solve(prob)` returns a [`PradelFit`](@ref).
"""
struct PradelProblem
    groups::Vector{PradelGroup}
    designs::Vector{NamedTuple{(:φ, :f, :p),NTuple{3,Matrix{Float64}}}}
    blocks::NamedTuple{(:φ, :f, :p),NTuple{3,UnitRange{Int}}}
    names::Vector{Symbol}
end

function PradelProblem(groups::AbstractVector{PradelGroup}, designs::AbstractVector;
                       names = nothing)
    length(groups) == length(designs) || throw(DimensionMismatch("one design per group"))
    isempty(groups) && throw(ArgumentError("need at least one group"))
    D = [(φ = Matrix{Float64}(d.φ), f = Matrix{Float64}(d.f), p = Matrix{Float64}(d.p))
         for d in designs]
    k = map(s -> size(D[1][s], 2), (φ = :φ, f = :f, p = :p))
    for (g, d) in zip(groups, D)
        size(d.φ, 1) == g.T - 1 && size(d.f, 1) == g.T - 1 ||
            throw(DimensionMismatch("φ and f designs need T − 1 rows"))
        size(d.p, 1) == length(g.primary) || throw(DimensionMismatch("p design needs one row per occasion"))
        all(s -> size(d[s], 2) == k[s], (:φ, :f, :p)) ||
            throw(DimensionMismatch("every group must have the same number of coefficients"))
    end
    blocks = (φ = 1:k.φ, f = k.φ .+ (1:k.f), p = (k.φ + k.f) .+ (1:k.p))
    nm = names === nothing ?
         vcat([Symbol("φ_", i) for i in 1:k.φ], [Symbol("f_", i) for i in 1:k.f],
              [Symbol("p_", i) for i in 1:k.p]) : collect(Symbol, names)
    length(nm) == last(blocks.p) || throw(DimensionMismatch("need $(last(blocks.p)) names"))
    return PradelProblem(collect(groups), D, blocks, nm)
end

LogDensityProblems.dimension(prob::PradelProblem) = length(prob.names)
LogDensityProblems.capabilities(::Type{PradelProblem}) = LogDensityProblems.LogDensityOrder{0}()

_pradel_logistic(x) = 1 / (1 + exp(-x))

function _pradel_rates(prob::PradelProblem, θ::AbstractVector, g::Integer)
    d = prob.designs[g]
    φ = _pradel_logistic.(d.φ * θ[prob.blocks.φ])
    f = exp.(d.f * θ[prob.blocks.f])
    p = _pradel_logistic.(d.p * θ[prob.blocks.p])
    return (φ = φ, f = f, p = p)
end

function LogDensityProblems.logdensity(prob::PradelProblem, θ::AbstractVector)
    length(θ) == LogDensityProblems.dimension(prob) ||
        throw(DimensionMismatch("expected $(LogDensityProblems.dimension(prob)) parameters"))
    ll = 0.0
    for g in eachindex(prob.groups)
        r = _pradel_rates(prob, θ, g)
        ll += pradel_loglikelihood(prob.groups[g], r.φ, r.f, r.p)
    end
    return isfinite(ll) ? ll : -Inf
end

"""
    nobs(prob::PradelProblem)

Sample size used for AICc: the total number of captures (sum of the capture histories). With
this sample size the AICc values of Hardy et al. (2023, Appendices S1.3–S1.5), which come from
Program MARK, are reproduced from their −2 log L and K.
"""
nobs(prob::PradelProblem) = sum(sum(g.captures) for g in prob.groups)

"""
    PradelFit

Maximum-likelihood fit of a [`PradelProblem`](@ref) from `solve(prob)`. Fields:

- `problem`, `θ` (coefficients, in the order of `problem.names`), `vcov` (inverse observed
  information of `θ`; coefficients with no curvature, such as an estimate on a boundary, are
  left out and get `NaN`; all `NaN` if the rest is singular) and `gradient` (at the optimum);
- `loglik`, `deviance` (−2 log L, MARK's "−2log(L)"), `k` (number of coefficients), `n`
  (`nobs(prob)`, the total number of captures) and `aicc`;
- `groups`: for each group a named tuple of `φ`, `f`, `λ = φ + f` (intervals), `p`
  (occasions), `pstar` and `N` (periods; `N = n_t / p*_t`).

Standard errors of derived quantities: [`pradel_se`](@ref).
"""
struct PradelFit
    problem::PradelProblem
    θ::Vector{Float64}
    vcov::Matrix{Float64}
    gradient::Vector{Float64}
    loglik::Float64
    deviance::Float64
    k::Int
    n::Int
    aicc::Float64
    groups::Vector{NamedTuple}
end

function _pradel_summary(prob::PradelProblem, θ::AbstractVector, g::Integer)
    r = _pradel_rates(prob, θ, g)
    G = prob.groups[g]
    return (φ = r.φ, f = r.f, λ = r.φ .+ r.f, p = r.p,
            pstar = pradel_capture_probability(G, r.p), N = pradel_abundance(G, r.p))
end

# solve(prob::PradelProblem; x0, iters) -> PradelFit; documented in PradelProblem.
function CommonSolve.solve(prob::PradelProblem;
                           x0::AbstractVector = _pradel_start(prob), iters::Integer = 3000,
                           gtol::Real = 1e-4, vcov::Bool = true)
    obj(θ) = -LogDensityProblems.logdensity(prob, θ)
    r = _ms_bfgs(obj, x0; iters, gtol)
    V = vcov ? _pradel_vcov(_ms_hessian(obj, r.x)) : fill(NaN, length(r.x), length(r.x))
    ll = -r.f
    k, n = length(r.x), nobs(prob)
    aicc = -2ll + 2k + 2k * (k + 1) / (n - k - 1)
    return PradelFit(prob, r.x, V, r.gradient, ll, -2ll, k, n, aicc,
                     [_pradel_summary(prob, r.x, g) for g in eachindex(prob.groups)])
end

# Inverse of the observed information H (Hessian of −log L). Coefficients with no curvature
# (an estimate on the boundary, such as a survival probability of zero, or a parameter the data
# do not inform) are left out and get NaN, one at a time (the one loading most on the
# flattest direction) until the rest can be inverted, as MARK reports the
# remaining standard errors when an estimate is at a boundary.
function _pradel_vcov(H::AbstractMatrix)
    d = size(H, 1)
    keep = [H[i, i] > 1e-4 for i in 1:d]
    V = fill(NaN, d, d)
    while any(keep)
        F = LinearAlgebra.eigen(Symmetric(H[keep, keep]))
        if F.values[1] > 1e-8 * max(1.0, F.values[end])
            V[keep, keep] .= F.vectors * Diagonal(1 ./ F.values) * F.vectors'
            return V
        end
        # drop the coefficient that loads most on the flattest direction, and try again
        idx = findall(keep)
        keep[idx[argmax(abs.(F.vectors[:, 1]))]] = false
    end
    return V
end

# start: φ = 0.6, f = 0.3, p = 0.3 when each block has an intercept-like first column
function _pradel_start(prob::PradelProblem)
    θ = zeros(LogDensityProblems.dimension(prob))
    θ[first(prob.blocks.f)] = log(0.3)
    θ[first(prob.blocks.p)] = log(0.3 / 0.7)
    θ[first(prob.blocks.φ)] = log(0.6 / 0.4)
    return θ
end

"""
    pradel_se(fit::PradelFit, quantity::Function) -> (value, se)

Value and delta-method standard error of a scalar function `quantity(θ)` of the coefficients
of a [`PradelFit`](@ref), with a central-difference gradient. Coefficients on which the
quantity does not depend (zero gradient) are skipped, so a boundary estimate elsewhere does
not make the result `NaN`. For example the realized growth
rate of group 2 in interval 5: `pradel_se(fit, θ -> pradel_derived(fit.problem, θ, 2).λ[5])`.
"""
function pradel_se(fit::PradelFit, quantity::Function)
    v = quantity(fit.θ)
    g = _ms_grad(quantity, fit.θ; h = 1e-6)
    use = abs.(g) .> 1e-9     # coefficients the quantity does not depend on are skipped
    se = any(use) ? sqrt(max(0.0, g[use]' * fit.vcov[use, use] * g[use])) : 0.0
    return (value = v, se = se)
end

"""
    pradel_derived(prob::PradelProblem, θ, g) -> NamedTuple

Natural-scale parameters of group `g` at coefficients `θ`: `φ`, `f`, `λ = φ + f`, `p`,
`pstar` and `N` (see [`PradelFit`](@ref)).
"""
pradel_derived(prob::PradelProblem, θ::AbstractVector, g::Integer) = _pradel_summary(prob, θ, g)

"""
    simulate_pradel(rng, N1, φ, f, p; primary = 1:length(p)) -> Matrix{Int}

Simulate capture histories under the Pradel model. `N1` animals are present in period 1;
between periods `t` and `t + 1` each survives with probability `φ[t]` and
`Poisson(f[t] N_t)` recruits enter. Each animal present in period `t` is caught at each
secondary occasion `j` (with `primary[j] = t`) with probability `p[j]`. Returns the histories
of the animals caught at least once (rows) by occasion (columns).
"""
function simulate_pradel(rng::AbstractRNG, N1::Integer, φ::AbstractVector, f::AbstractVector,
                         p::AbstractVector; primary::AbstractVector{<:Integer} = 1:length(p))
    T = last(primary)
    length(φ) == T - 1 == length(f) || throw(DimensionMismatch("φ and f need T − 1 entries"))
    alive = collect(1:N1)
    entry = ones(Int, N1); exit = fill(T, N1)
    for t in 1:T-1
        survivors = [i for i in alive if rand(rng) < φ[t]]
        for i in setdiff(alive, survivors)
            exit[i] = t
        end
        nrec = rand(rng, Poisson(f[t] * length(alive)))
        newids = length(entry) .+ (1:nrec)
        append!(entry, fill(t + 1, nrec)); append!(exit, fill(T, nrec))
        alive = vcat(survivors, collect(newids))
    end
    H = zeros(Int, length(entry), length(primary))
    for i in eachindex(entry), j in eachindex(primary)
        entry[i] <= primary[j] <= exit[i] && rand(rng) < p[j] && (H[i, j] = 1)
    end
    return H[vec(any(H .== 1; dims = 2)), :]
end
