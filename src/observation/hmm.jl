# Generic hidden Markov model (HMM) layer: forward (scaled, log-likelihood accumulated in log
# space), backward, smoothed state probabilities, Viterbi, covariate-dependent transition
# matrices, categorical emissions and a LogDensityProblems wrapper. Model-agnostic: the
# outbreak-state HMM of Wilber et al. (2022) (src/models/wilber2022.jl) is one user.
#
# Conventions. K hidden states, T time steps. `initial` is the distribution of the state at
# t = 1 (length K). `transitions` is one K×K row-stochastic matrix (time-homogeneous) or a
# vector of T − 1 matrices, the t-th mapping the state at t to the state at t + 1:
# P[t][i, j] = Pr(X_{t+1} = j | X_t = i). `likelihoods` is a T×K matrix with
# L[t, k] = p(y_t | X_t = k); a missing observation is a row of ones. Emission rows need not
# sum to one over observations (any non-negative weights are accepted), so models whose
# "emission" is a likelihood factor, not a distribution, can use the same code.
# Everything is generic in the element type, so AD through the transitions works.

_hmm_P(P::AbstractMatrix, t) = P
_hmm_P(P::AbstractVector{<:AbstractMatrix}, t) = P[t]

function _hmm_check(initial, transitions, L)
    T, K = size(L)
    T >= 1 || throw(ArgumentError("need at least one time step"))
    length(initial) == K ||
        throw(DimensionMismatch("initial has length $(length(initial)), likelihoods have $K states"))
    if transitions isa AbstractMatrix
        size(transitions) == (K, K) || throw(DimensionMismatch("transition matrix must be $K×$K"))
    else
        length(transitions) >= T - 1 ||
            throw(DimensionMismatch("need $(T - 1) transition matrices, got $(length(transitions))"))
        all(P -> size(P) == (K, K), transitions) ||
            throw(DimensionMismatch("transition matrices must be $K×$K"))
    end
    return T, K
end

"""
    hmm_forward(initial, transitions, likelihoods; emit_first = true) -> NamedTuple

Scaled forward recursion of a hidden Markov model with `K` states and `T` steps.
`initial` is the state distribution at the first step; `transitions` is one `K×K`
row-stochastic matrix or a vector of `T − 1` of them (the `t`-th gives
`Pr(X_{t+1} = j | X_t = i)` in row `i`, column `j`); `likelihoods` is the `T×K` matrix of
emission likelihoods `p(y_t | X_t = k)`, with a row of ones for a missing observation.

The filtered vector is renormalised at every step and the log normalising constants are
summed, so the log-likelihood does not underflow for long sequences. With
`emit_first = false` the first row of `likelihoods` is not used: the first observation only
fixes the start (the convention of the Wilber et al. 2022 Stan code, which starts every
population in a known state).

Returns `loglik`, `filtered` (`T×K`, row `t` is `Pr(X_t | y_1:t)`) and `logc` (the log
normalising constant of each step; `sum(logc) == loglik`). A zero-probability sequence gives
`loglik = -Inf` and `NaN` filtered rows from that step on.
"""
function hmm_forward(initial::AbstractVector, transitions, L::AbstractMatrix;
                     emit_first::Bool = true)
    T, K = _hmm_check(initial, transitions, L)
    R = promote_type(eltype(initial), eltype(L), _hmm_eltype(transitions))
    α = Matrix{R}(undef, T, K)
    logc = Vector{R}(undef, T)
    a = emit_first ? initial .* view(L, 1, :) : collect(R, initial)
    for t in 1:T
        t > 1 && (a = (transpose(_hmm_P(transitions, t - 1)) * α[t - 1, :]) .* view(L, t, :))
        c = sum(a)
        if !(c > 0)                       # zero-probability sequence: stop here
            logc[t] = -Inf
            logc[t+1:T] .= 0
            α[t:T, :] .= NaN
            return (loglik = convert(R, -Inf), filtered = α, logc = logc)
        end
        logc[t] = log(c)
        α[t, :] .= a ./ c
    end
    return (loglik = sum(logc), filtered = α, logc = logc)
end

_hmm_eltype(P::AbstractMatrix) = eltype(P)
_hmm_eltype(P::AbstractVector{<:AbstractMatrix}) = isempty(P) ? Float64 : eltype(first(P))


"""
    hmm_loglikelihood(initial, transitions, likelihoods; emit_first = true)

Log-likelihood of one observed sequence: [`hmm_forward`](@ref)`(…).loglik`.
"""
hmm_loglikelihood(initial, transitions, L; emit_first::Bool = true) =
    hmm_forward(initial, transitions, L; emit_first).loglik

"""
    hmm_backward(transitions, likelihoods, logc) -> Matrix

Scaled backward recursion matching [`hmm_forward`](@ref): with `logc` the forward log
normalising constants, row `t` of the result is `β_t(k) = p(y_{t+1:T} | X_t = k)` divided by
`exp(sum(logc[t+1:T]))`, so `filtered .* backward` gives smoothed probabilities.
"""
function hmm_backward(transitions, L::AbstractMatrix, logc::AbstractVector)
    T, K = size(L)
    R = promote_type(eltype(L), _hmm_eltype(transitions), eltype(logc))
    β = ones(R, T, K)
    for t in (T - 1):-1:1
        β[t, :] .= (_hmm_P(transitions, t) * (view(β, t + 1, :) .* view(L, t + 1, :))) ./
                   exp(logc[t + 1])
    end
    return β
end

"""
    hmm_posterior(initial, transitions, likelihoods; emit_first = true) -> Matrix

Smoothed state probabilities `Pr(X_t = k | y_1:T)` (`T×K`, rows sum to one) by the
forward–backward algorithm.
"""
function hmm_posterior(initial, transitions, L::AbstractMatrix; emit_first::Bool = true)
    f = hmm_forward(initial, transitions, L; emit_first)
    isfinite(f.loglik) || throw(DomainError(f.loglik, "the sequence has zero probability"))
    β = hmm_backward(transitions, L, f.logc)
    γ = f.filtered .* β
    return γ ./ sum(γ; dims = 2)
end

"""
    hmm_viterbi(initial, transitions, likelihoods; emit_first = true) -> NamedTuple

Most probable hidden state sequence (Viterbi, in log space). Returns `path` (a vector of
state indices) and `logprob`, the joint log probability of `path` and the observations.
Ties go to the lowest state index.
"""
function hmm_viterbi(initial, transitions, L::AbstractMatrix; emit_first::Bool = true)
    T, K = _hmm_check(initial, transitions, L)
    δ = Matrix{Float64}(undef, T, K)
    ψ = zeros(Int, T, K)
    δ[1, :] .= log.(initial) .+ (emit_first ? log.(view(L, 1, :)) : 0.0)
    for t in 2:T
        lP = log.(_hmm_P(transitions, t - 1))
        for j in 1:K
            best, arg = -Inf, 1
            for i in 1:K
                v = δ[t - 1, i] + lP[i, j]
                if v > best
                    best, arg = v, i
                end
            end
            δ[t, j] = best + log(L[t, j])
            ψ[t, j] = arg
        end
    end
    lp, last = findmax(view(δ, T, :))
    path = zeros(Int, T)
    path[T] = last
    for t in (T - 1):-1:1
        path[t] = ψ[t + 1, path[t + 1]]
    end
    return (path = path, logprob = lp)
end

"""
    hmm_sample(rng, initial, transitions, T) -> Vector{Int}

Draw a hidden state path of length `T` from the Markov chain (states only; draw observations
from the path with the model's emission).
"""
function hmm_sample(rng::AbstractRNG, initial::AbstractVector, transitions, T::Integer)
    path = zeros(Int, T)
    path[1] = _hmm_draw(rng, initial)
    for t in 2:T
        path[t] = _hmm_draw(rng, view(_hmm_P(transitions, t - 1), path[t - 1], :))
    end
    return path
end

function _hmm_draw(rng, p)
    u = rand(rng) * sum(p)
    c = zero(u)
    for (k, pk) in enumerate(p)
        c += pk
        u <= c && return k
    end
    return findlast(>(0), p)
end

"""
    CategoricalEmission(B)

Emission with a finite observation alphabet: `B[k, o] = Pr(Y = o | X = k)`. Use
[`emission_likelihoods`](@ref) to turn a sequence of observation codes (`missing` for no
observation) into the `T×K` likelihood matrix the HMM functions take.
"""
struct CategoricalEmission{M<:AbstractMatrix}
    B::M
end

"""
    emission_likelihoods(e::CategoricalEmission, obs) -> Matrix

The `T×K` matrix of emission likelihoods for observation codes `obs` (integers indexing the
columns of `e.B`, or `missing`, which gives a row of ones).
"""
function emission_likelihoods(e::CategoricalEmission, obs::AbstractVector)
    K = size(e.B, 1)
    L = ones(eltype(e.B), length(obs), K)
    for (t, o) in enumerate(obs)
        o === missing && continue
        L[t, :] .= view(e.B, :, o)
    end
    return L
end

"""
    CovariateTransitions(f, covariates)

Covariate-dependent transition matrices: `f(θ, x)` returns the `K×K` row-stochastic matrix
for parameters `θ` and one covariate row `x`, and `covariates` holds one row per transition
(`T − 1` rows; row `t` drives the step from `t` to `t + 1`). Call it as `ct(θ)` to get the
vector of matrices for [`hmm_forward`](@ref).
"""
struct CovariateTransitions{F,X<:AbstractMatrix}
    f::F
    covariates::X
end

(ct::CovariateTransitions)(θ) = [ct.f(θ, x) for x in eachrow(ct.covariates)]

"""
    HMMLogDensity(sequences, logprior, dimension; emit_first = true)

A LogDensityProblems target for a set of independent HMM sequences that share parameters.
`sequences(θ)` returns an iterable of `(initial, transitions, likelihoods)` tuples (or named
tuples with those fields), one per sequence; `logprior(θ)` is the log prior density on the
same (unconstrained) parameter vector. `logdensity` is the prior plus the sum of the
[`hmm_loglikelihood`](@ref)s; [`pointwise_loglikelihood`](@ref) returns one value per
sequence.
"""
struct HMMLogDensity{F,G}
    sequences::F
    logprior::G
    dimension::Int
    emit_first::Bool
end

HMMLogDensity(sequences, logprior, dimension::Integer; emit_first::Bool = true) =
    HMMLogDensity(sequences, logprior, Int(dimension), emit_first)

LogDensityProblems.capabilities(::Type{<:HMMLogDensity}) = LogDensityProblems.LogDensityOrder{0}()
LogDensityProblems.dimension(m::HMMLogDensity) = m.dimension

_hmm_triple(s::Tuple) = s
_hmm_triple(s::NamedTuple) = (s.initial, s.transitions, s.likelihoods)

function pointwise_loglikelihood(m::HMMLogDensity, θ::AbstractVector)
    return [hmm_loglikelihood(_hmm_triple(s)...; emit_first = m.emit_first)
            for s in m.sequences(θ)]
end

function LogDensityProblems.logdensity(m::HMMLogDensity, θ::AbstractVector)
    lp = m.logprior(θ)
    isfinite(lp) || return oftype(float(lp), -Inf)
    ll = sum(pointwise_loglikelihood(m, θ))
    return isnan(ll) ? oftype(float(lp + ll), -Inf) : lp + ll
end
