# Open-population capture–recapture likelihoods that extend the multistate layer
# (src/observation/multistate.jl): the conditional Arnason–Schwarz model with unequal
# intervals (as M-SURGE and MARK fit it) and the POPAN (Schwarz–Arnason) formulation of the
# Jolly–Seber model, with super-population size. Plus Akaike weights for model averaging.
# User: Brannelly et al. 2018 (src/models/brannelly2018.jl).
#
# Conditional Arnason–Schwarz (Arnason 1973; Schwarz, Schweigert & Arnason 1993, Biometrics
# 49:177; Lebreton et al. 2009, Adv Ecol Res 41:87). States are observed without error when
# an animal is caught (event k = state k), and the likelihood is conditioned on the first
# capture AND on the state at first capture, so no initial-state parameters enter (unlike
# MultiEventProblem, which estimates the state distribution at first capture). Survival is a
# probability per unit time; over an interval of length τ it is S^τ (the MARK/M-SURGE "time
# intervals" convention). Transitions and detection are per interval and per occasion.
#
# POPAN (Schwarz & Arnason 1996, Biometrics 52:860). A super-population of N animals enters
# the study area over the T occasions: a fraction pent_b is present for the first time at
# occasion b (Σ pent_b = 1). Present animals are caught with probability p_t and survive to
# the next occasion with probability φ_t^τ_t. With u distinct animals caught and n_ω animals
# with history ω,
#
#   L = N! / ((N − u)! Π n_ω!) · Π P(ω)^{n_ω} · P(0)^{N − u},
#
#   ψ_1 = pent_1,  ψ_{t+1} = ψ_t (1 − p_t) φ_t + pent_{t+1}   (present, not yet caught),
#   χ_T = 1,       χ_t = (1 − φ_t) + φ_t (1 − p_{t+1}) χ_{t+1}  (never caught after t),
#   P(ω) = ψ_f p_f Π_{t=f}^{l−1} φ_t p_{t+1}^{ω_{t+1}} (1 − p_{t+1})^{1−ω_{t+1}} · χ_l,
#   P(0) = 1 − Σ_t ψ_t p_t,
#
# for first capture f and last capture l. As in MARK, φ and p are logit-linear, pent is a
# multinomial logit with pent_1 the baseline, and N = u + exp(θ_N) (log link on N − u).

const POPAN_STRUCTURES = (:constant, :time)

"""
    ArnasonSchwarzProblem(histories, structure::MultiEventStructure; intervals = ones(T − 1))

The conditional Arnason–Schwarz multistate model (as fitted by M-SURGE and MARK): `histories`
is an `n×T` matrix with `0` (not caught) or the state `k ∈ 1…K` of an animal at capture, and
`structure` gives the survival, transition and detection designs (its `assignment` must be
the `K×K` identity: states are observed without error). The likelihood is conditioned on the
first capture and on the state there, so, unlike [`MultiEventProblem`](@ref), no
initial-state parameter is estimated. Survival is per unit time: over interval `t` of length
`intervals[t]` it is `S^intervals[t]`; transitions and detection are per interval and per
occasion.

A LogDensityProblems target (log-likelihood, flat prior on the unconstrained parameters
`prob.names`); `solve(prob; x0, iters)` maximises it and returns an [`ArnasonSchwarzFit`](@ref).
"""
struct ArnasonSchwarzProblem
    inner::MultiEventProblem
    intervals::Vector{Float64}
    free::Vector{Int}
    names::Vector{Symbol}
end

function ArnasonSchwarzProblem(histories::AbstractMatrix, st::MultiEventStructure;
                               intervals::AbstractVector{<:Real} = ones(size(histories, 2) - 1))
    K = st.nstates
    size(st.assignment) == (K, K) && st.assignment == Matrix{Bool}(LinearAlgebra.I(K)) ||
        throw(ArgumentError("an Arnason–Schwarz model needs the identity state assignment"))
    length(intervals) == size(histories, 2) - 1 ||
        throw(DimensionMismatch("need $(size(histories, 2) - 1) interval lengths"))
    all(>(0), intervals) || throw(ArgumentError("interval lengths must be positive"))
    inner = MultiEventProblem(histories, st)
    free = setdiff(1:length(inner.names), inner.blocks.initial)
    return ArnasonSchwarzProblem(inner, collect(Float64, intervals), free, inner.names[free])
end

LogDensityProblems.dimension(p::ArnasonSchwarzProblem) = length(p.free)
LogDensityProblems.capabilities(::Type{<:ArnasonSchwarzProblem}) =
    LogDensityProblems.LogDensityOrder{0}()

function _as_full(prob::ArnasonSchwarzProblem, θ::AbstractVector)
    x = zeros(eltype(θ), length(prob.inner.names))
    x[prob.free] .= θ
    return x
end

"""
    arnason_schwarz_parameters(prob::ArnasonSchwarzProblem, θ) -> NamedTuple

Natural-scale parameters at `θ`: `survival` (`K×(T − 1)`, per unit time; raise to the
interval length for the interval survival), `transition` (`K×K×(T − 1)`) and `detection`
(`K×(T − 1)`, column `t` is occasion `t + 1`).
"""
function arnason_schwarz_parameters(prob::ArnasonSchwarzProblem, θ::AbstractVector)
    par = multievent_parameters(prob.inner, _as_full(prob, θ))
    return (survival = par.survival, transition = par.transition, detection = par.detection)
end

function _as_logliks(prob::ArnasonSchwarzProblem, θ::AbstractVector)
    par = arnason_schwarz_parameters(prob, θ)
    H = prob.inner.histories
    K = prob.inner.structure.nstates
    T = size(H, 2)
    Φ = [survival_transition(view(par.survival, :, t) .^ prob.intervals[t],
                             view(par.transition, :, :, t)) for t in 1:T-1]
    δ = Matrix{eltype(par.detection)}(LinearAlgebra.I(K))
    B = [multievent_emission(view(par.detection, :, t - 1), δ) for t in 2:T]
    R = eltype(B[1])
    out = Vector{R}(undef, size(H, 1))
    for r in axes(H, 1)
        f = prob.inner.first[r]
        init = zeros(R, K + 1)
        init[H[r, f]] = 1
        L = ones(R, T - f + 1, K + 1)
        for t in f+1:T
            o = H[r, t]
            o === missing && continue
            L[t-f+1, :] .= view(B[t-1], :, o + 1)
        end
        out[r] = hmm_loglikelihood(init, Φ[f:T-1], L)
    end
    return out
end

pointwise_loglikelihood(prob::ArnasonSchwarzProblem, θ::AbstractVector) = _as_logliks(prob, θ)

function LogDensityProblems.logdensity(prob::ArnasonSchwarzProblem, θ::AbstractVector)
    ll = sum(prob.inner.counts .* _as_logliks(prob, θ))
    return isnan(ll) ? oftype(ll, -Inf) : ll
end

"""
    ArnasonSchwarzFit

Maximum-likelihood fit of an [`ArnasonSchwarzProblem`](@ref): `problem`, `θ`, `vcov` (inverse
observed information; `NaN` when singular), `rank_deficiency`, `loglik`, `deviance`
(`−2 loglik`, as M-SURGE and MARK print it for this model), `k` (estimable parameters),
`aic`, `parameters` (from [`arnason_schwarz_parameters`](@ref)) and `gradient`.
"""
struct ArnasonSchwarzFit
    problem::ArnasonSchwarzProblem
    θ::Vector{Float64}
    vcov::Matrix{Float64}
    rank_deficiency::Int
    loglik::Float64
    deviance::Float64
    k::Int
    aic::Float64
    parameters::NamedTuple
    gradient::Vector{Float64}
end

# solve(prob::ArnasonSchwarzProblem; x0 = zeros, iters = 2000). Documented in the problem.
function CommonSolve.solve(prob::ArnasonSchwarzProblem;
                           x0::AbstractVector = zeros(LogDensityProblems.dimension(prob)),
                           iters::Integer = 2000)
    f(θ) = -LogDensityProblems.logdensity(prob, θ)
    r = _ms_bfgs(f, x0; iters)
    V, def = _ms_vcov(-_ms_hessian(f, r.x))
    ll = -r.f
    k = length(r.x) - def
    return ArnasonSchwarzFit(prob, r.x, V, def, ll, -2ll, k, -2ll + 2k,
                             arnason_schwarz_parameters(prob, r.x), r.gradient)
end

# --- POPAN ---------------------------------------------------------------------------------------

"""
    PopanProblem(histories; survival = :constant, detection = :time, entry = :time,
                 intervals = ones(T − 1))

The POPAN (Schwarz & Arnason 1996) open-population model with super-population size `N`.
`histories` is an `n×T` matrix; any positive entry is a capture (states are ignored). Each
of `survival` (`φ`, per unit time, `φ^τ` over an interval of length `τ`), `detection` (`p`
at each of the `T` occasions) and `entry` (`pent_2…pent_T`, multinomial logit with `pent_1`
as baseline) is `:constant` (MARK's "(.)") or `:time` ("(t)"). `N = u + exp(θ_N)`, `u` the
number of distinct animals caught.

The log density is the full multinomial log-likelihood (with its `N!/(N − u)!` and
`1/Π n_ω!` terms) at the unconstrained parameters `prob.names`. `solve(prob)` returns a
[`PopanFit`](@ref); [`popan_parameters`](@ref) gives the natural-scale parameters.
"""
struct PopanProblem
    seen::Matrix{Bool}
    counts::Vector{Int}
    u::Int
    intervals::Vector{Float64}
    survival::Symbol
    detection::Symbol
    entry::Symbol
    blocks::NamedTuple
    names::Vector{Symbol}
end

function PopanProblem(histories::AbstractMatrix; survival::Symbol = :constant,
                      detection::Symbol = :time, entry::Symbol = :time,
                      intervals::AbstractVector{<:Real} = ones(size(histories, 2) - 1))
    n, T = size(histories)
    T >= 3 || throw(ArgumentError("need at least three occasions"))
    for s in (survival, detection, entry)
        s in POPAN_STRUCTURES || throw(ArgumentError("structure :$s not in $POPAN_STRUCTURES"))
    end
    length(intervals) == T - 1 || throw(DimensionMismatch("need $(T - 1) interval lengths"))
    pooled = Dict{Vector{Bool},Int}()
    order = Vector{Vector{Bool}}()
    for i in 1:n
        h = [x !== missing && x > 0 for x in histories[i, :]]
        any(h) || continue
        haskey(pooled, h) || push!(order, h)
        pooled[h] = get(pooled, h, 0) + 1
    end
    S = Matrix{Bool}(undef, length(order), T)
    for (r, h) in enumerate(order)
        S[r, :] .= h
    end
    counts = [pooled[h] for h in order]
    sizes = (survival = survival === :time ? T - 1 : 1, detection = detection === :time ? T : 1,
             entry = entry === :time ? T - 1 : 1, N = 1)
    names = Symbol[]
    blocks = Dict{Symbol,UnitRange{Int}}()
    for (c, m) in pairs(sizes)
        blocks[c] = length(names) .+ (1:m)
        append!(names, [Symbol(c, "_", j) for j in 1:m])
    end
    nb = NamedTuple{keys(sizes)}(Tuple(blocks[c] for c in keys(sizes)))
    return PopanProblem(S, counts, sum(counts), collect(Float64, intervals), survival,
                        detection, entry, nb, names)
end

LogDensityProblems.dimension(p::PopanProblem) = length(p.names)
LogDensityProblems.capabilities(::Type{<:PopanProblem}) = LogDensityProblems.LogDensityOrder{0}()

"""
    popan_parameters(prob::PopanProblem, θ) -> NamedTuple

Natural-scale POPAN parameters at `θ`: `survival` (length `T − 1`, per unit time),
`detection` (length `T`), `entry` (`pent`, length `T`, sums to one) and `N` (the
super-population size).
"""
function popan_parameters(prob::PopanProblem, θ::AbstractVector)
    T = size(prob.seen, 2)
    φ = _logistic.(θ[prob.blocks.survival])
    φ = prob.survival === :time ? φ : fill(φ[1], T - 1)
    p = _logistic.(θ[prob.blocks.detection])
    p = prob.detection === :time ? p : fill(p[1], T)
    η = θ[prob.blocks.entry]
    η = prob.entry === :time ? η : fill(η[1], T - 1)
    w, w0 = _ms_softmax(η)
    pent = vcat(w0, w)
    N = prob.u + exp(θ[prob.blocks.N[1]])
    return (survival = φ, detection = p, entry = pent, N = N)
end

function LogDensityProblems.logdensity(prob::PopanProblem, θ::AbstractVector)
    par = popan_parameters(prob, θ)
    S = prob.seen
    T = size(S, 2)
    φ = par.survival .^ prob.intervals
    p, pent, N = par.detection, par.entry, par.N
    R = eltype(p)
    ψ = Vector{R}(undef, T)
    ψ[1] = pent[1]
    for t in 1:T-1
        ψ[t+1] = ψ[t] * (1 - p[t]) * φ[t] + pent[t+1]
    end
    χ = Vector{R}(undef, T)
    χ[T] = one(R)
    for t in T-1:-1:1
        χ[t] = (1 - φ[t]) + φ[t] * (1 - p[t+1]) * χ[t+1]
    end
    P0 = 1 - sum(ψ .* p)
    lg = Distributions.SpecialFunctions.loggamma
    ll = lg(N + 1) - lg(N - prob.u + 1) + (N - prob.u) * log(max(P0, zero(R)))
    for r in axes(S, 1)
        f = findfirst(view(S, r, :))
        l = findlast(view(S, r, :))
        lp = log(ψ[f] * p[f]) + log(χ[l])
        for t in f:l-1
            lp += log(φ[t]) + (S[r, t+1] ? log(p[t+1]) : log1p(-p[t+1]))
        end
        ll += prob.counts[r] * lp - lg(prob.counts[r] + 1)
    end
    return isnan(ll) ? oftype(ll, -Inf) : ll
end

"""
    PopanFit

Maximum-likelihood fit of a [`PopanProblem`](@ref): `problem`, `θ`, `vcov`,
`rank_deficiency`, `loglik`, `deviance` (`−2 loglik`; MARK's POPAN deviance is on another
basis), `k`, `aic`, `parameters` (from [`popan_parameters`](@ref)), `N` (super-population
size), `N_ci` (95% interval `u + f₀ exp(±1.96 SE(log f₀))`, `f₀ = N − u`, the log-based
interval MARK reports) and `gradient`.
"""
struct PopanFit
    problem::PopanProblem
    θ::Vector{Float64}
    vcov::Matrix{Float64}
    rank_deficiency::Int
    loglik::Float64
    deviance::Float64
    k::Int
    aic::Float64
    parameters::NamedTuple
    N::Float64
    N_ci::Tuple{Float64,Float64}
    gradient::Vector{Float64}
end

# solve(prob::PopanProblem; x0, iters = 3000). Documented in PopanProblem.
function CommonSolve.solve(prob::PopanProblem; x0::Union{Nothing,AbstractVector} = nothing,
                           iters::Integer = 3000)
    if x0 === nothing
        x0 = zeros(LogDensityProblems.dimension(prob))
        x0[prob.blocks.survival] .= 2.0
        x0[prob.blocks.detection] .= -1.0
        x0[prob.blocks.N] .= log(prob.u)
    end
    f(θ) = -LogDensityProblems.logdensity(prob, θ)
    r = _ms_bfgs(f, x0; iters)
    V, def = _ms_vcov(-_ms_hessian(f, r.x))
    ll = -r.f
    k = length(r.x) - def
    par = popan_parameters(prob, r.x)
    iN = prob.blocks.N[1]
    f0 = exp(r.x[iN])
    se = sqrt(V[iN, iN])
    ci = isfinite(se) ? (prob.u + f0 * exp(-1.96se), prob.u + f0 * exp(1.96se)) : (NaN, NaN)
    return PopanFit(prob, r.x, V, def, ll, -2ll, k, -2ll + 2k, par, par.N, ci, r.gradient)
end

"""
    akaike_weights(aic) -> Vector

Akaike weights `exp(−Δᵢ/2) / Σⱼ exp(−Δⱼ/2)` of a set of models with information criteria
`aic` (`Δᵢ = aicᵢ − min aic`; Burnham & Anderson 2002). A model-averaged estimate is
`sum(weights .* estimates)`.
"""
function akaike_weights(aic::AbstractVector{<:Real})
    Δ = aic .- minimum(aic)
    w = exp.(-Δ ./ 2)
    return w ./ sum(w)
end
