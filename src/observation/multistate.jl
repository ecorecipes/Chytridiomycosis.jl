# Generic multistate capture–recapture layer: Cormack–Jolly–Seber (CJS) and multi-event models
# with state uncertainty, discrete- or continuous-time state processes, and robust-design
# (secondary-occasion) observations. Built on the HMM layer (src/observation/hmm.jl): every
# individual capture history is one HMM sequence, conditioned on first capture.
# Users: Palomar et al. 2023
# (src/models/palomar2023.jl) and Hollanders et al. 2023 (src/models/hollanders2023.jl).
#
# Conventions. K alive states, then one absorbing dead state (index K + 1). Events are coded
# 0 (not seen), 1, …, E − 1 (seen, with whatever the field record says), as in E-SURGE and
# MARK; `missing` marks an occasion with no survey. An individual enters at its first
# capture f: the state there is drawn from `initial` and only the event's state-assignment
# probability enters the likelihood (not detection), so histories are conditioned on first
# capture (Pradel 2005, Biometrics 61:442; Lebreton et al. 2009, Adv Ecol Res 41:87).

# --- state processes ---------------------------------------------------------------------------

"""
    ctmc_transition(Q, τ) -> Matrix

Transition probability matrix `exp(Qτ)` of a continuous-time Markov chain with generator
`Q` (rows sum to zero, off-diagonals are rates) over an interval of length `τ` in the time
unit of the rates. Unequal intervals between occasions are handled by passing each one's
`τ` (Conn et al. 2012; Glennie et al. 2022). A 3×3 generator whose last row is zero (two
transient states and an absorbing state, as in [`infection_generator`](@ref)) uses a closed
form that is exact and differentiable; other generators use `LinearAlgebra.exp`.
"""
function ctmc_transition(Q::AbstractMatrix, τ::Real)
    n = LinearAlgebra.checksquare(Q)
    τ >= 0 || throw(DomainError(τ, "the interval must be non-negative"))
    if n == 3 && iszero(Q[3, 1]) && iszero(Q[3, 2]) && iszero(Q[3, 3])
        return _ctmc_absorbing2(Q[1, 1], Q[1, 2], Q[2, 1], Q[2, 2], τ)
    end
    return exp(Matrix(Q) .* τ)
end

# exp(Aτ) of the 2×2 transient block A = [a11 a12; a21 a22] (a12, a21 ≥ 0), via
# exp(Aτ) = e^{sτ}[cosh(qτ) I + sinh(qτ)/q (A − sI)], s = tr A / 2, q² = ((a11 − a22)/2)² + a12 a21.
function _ctmc_absorbing2(a11, a12, a21, a22, τ)
    s = (a11 + a22) / 2
    q = sqrt(((a11 - a22) / 2)^2 + a12 * a21)
    e1, e2 = exp((s + q) * τ), exp((s - q) * τ)
    C = (e1 + e2) / 2
    S = q > 1e-10 ? (e1 - e2) / (2q) : τ * exp(s * τ)
    m11, m12 = C + S * (a11 - s), S * a12
    m21, m22 = S * a21, C + S * (a22 - s)
    z = zero(m11)
    return [m11 m12 max(1 - m11 - m12, z)
            m21 m22 max(1 - m21 - m22, z)
            z z one(z)]
end

"""
    infection_generator(gain, clearance, mortality_uninfected, mortality_infected) -> Matrix

The 3×3 generator of the infection–survival process of Hollanders et al. (2023), states
(alive uninfected, alive infected, dead), all arguments hazard rates in one time unit:

    [−(gain + μ_U)   gain              μ_U
     clearance       −(clearance + μ_I) μ_I
     0               0                 0  ]

Pass it to [`ctmc_transition`](@ref) with the interval length.
"""
function infection_generator(gain, clearance, mortality_uninfected, mortality_infected)
    z = zero(promote_type(typeof(gain), typeof(clearance), typeof(mortality_uninfected),
                          typeof(mortality_infected)))
    return [-(gain + mortality_uninfected) gain mortality_uninfected
            clearance -(clearance + mortality_infected) mortality_infected
            z z z]
end

"""
    survival_transition(S, Ψ) -> Matrix

Discrete-time multistate transition matrix with `K = length(S)` alive states and a dead
state appended: survive with probability `S[k]`, then move between alive states with the
row-stochastic `K×K` matrix `Ψ` (the two-step S then T pattern of E-SURGE). Entry `[k, l]` is
`S[k] Ψ[k, l]`, `[k, K + 1]` is `1 − S[k]`, and the dead state is absorbing.
"""
function survival_transition(S::AbstractVector, Ψ::AbstractMatrix)
    K = length(S)
    size(Ψ) == (K, K) || throw(DimensionMismatch("Ψ must be $K×$K"))
    R = promote_type(eltype(S), eltype(Ψ))
    Φ = zeros(R, K + 1, K + 1)
    for k in 1:K
        for l in 1:K
            Φ[k, l] = S[k] * Ψ[k, l]
        end
        Φ[k, K+1] = 1 - S[k]
    end
    Φ[K+1, K+1] = one(R)
    return Φ
end

"""
    multievent_emission(p, δ) -> Matrix

Observation matrix of a multi-event model with `K = length(p)` alive states and a dead
state: an animal in state `k` is encountered with probability `p[k]` (matrix B1 of E-SURGE),
and an encountered animal gives seen-event `e` with probability `δ[k, e]` (B2; rows of `δ`
sum to one). Column 1 is event 0 (not seen), column `e + 1` is event `e`. The dead state is
never seen. With `δ` the identity this is a multistate (Arnason–Schwarz) emission; with
off-diagonal or "unknown" columns it carries state uncertainty (Pradel 2005).
"""
function multievent_emission(p::AbstractVector, δ::AbstractMatrix)
    K, E = size(δ)
    length(p) == K || throw(DimensionMismatch("p has length $(length(p)), δ has $K rows"))
    R = promote_type(eltype(p), eltype(δ))
    B = zeros(R, K + 1, E + 1)
    for k in 1:K
        B[k, 1] = 1 - p[k]
        for e in 1:E
            B[k, e+1] = p[k] * δ[k, e]
        end
    end
    B[K+1, 1] = one(R)
    return B
end

"""
    robust_design_likelihoods(emission, codes) -> Matrix

The `T×K` emission-likelihood matrix (for [`hmm_forward`](@ref)) of one individual under the
robust design: the state is constant within a primary occasion and each of its secondary
surveys is an independent observation. `codes` is a `T×J` matrix of observation-column
indices (`missing` where there was no survey); `emission(t, j)` returns the `K×O`
observation matrix of secondary `j` of primary `t` (or pass one matrix for all). Row `t` of
the result is the product over the non-missing secondaries; a primary with none is a row of
ones.
"""
function robust_design_likelihoods(emission, codes::AbstractMatrix)
    T, J = size(codes)
    E1 = _rd_emission(emission, 1, 1)
    K = size(E1, 1)
    L = ones(eltype(E1), T, K)
    for t in 1:T, j in 1:J
        o = codes[t, j]
        o === missing && continue
        E = _rd_emission(emission, t, j)
        for k in 1:K
            L[t, k] *= E[k, o]
        end
    end
    return L
end

_rd_emission(E::AbstractMatrix, t, j) = E
_rd_emission(f, t, j) = f(t, j)

"""
    quasistationary_distribution(Φ) -> Vector

Long-run distribution over the alive states of a transition matrix `Φ` whose last state is
an absorbing death: the normalised left Perron vector of the alive block (the state mix of
the survivors once the initial condition is forgotten). For two states (infected,
uninfected) its first entry is the infection prevalence the transition and survival
estimates imply among living animals.
"""
function quasistationary_distribution(Φ::AbstractMatrix)
    n = LinearAlgebra.checksquare(Φ)
    M = Matrix(Φ[1:n-1, 1:n-1])
    F = eigen(transpose(M))
    i = argmax(real.(F.values))
    v = abs.(real.(F.vectors[:, i]))
    return v ./ sum(v)
end

# --- multi-event model structures (E-SURGE-style design matrices) -------------------------------

const MULTIEVENT_STRUCTURES = (:constant, :state, :time, :state_plus_time, :state_time)

"""
    MultiEventStructure(; nstates, assignment, baseline, survival = :state,
                          transition = :state, detection = :state, state_assignment = :state)

The model pattern of a multi-event capture–recapture model (Pradel 2005) with `nstates`
alive states and a dead state.

- `assignment` is a `K×E` `Bool` matrix: `assignment[k, e]` is true when an encountered animal
  in state `k` can give seen-event `e` (events `1…E`; event 0 is "not seen").
  `baseline[k]` is the column whose probability is one minus the others (E-SURGE's `*`).
- Each parameter class is logit- (or multinomial-logit-) linear with a design given by a
  symbol: `:constant`, `:state`, `:time`, `:state_plus_time` (additive) or `:state_time`
  (interaction). Survival `S[k, t]` and transitions `ψ[k → l, t]` vary over the `T − 1`
  intervals, detection `p[k, t]` over occasions `2…T`; state assignment has no time effect
  (`:constant` or `:state`). Transitions are a multinomial logit on each row with staying as
  the baseline; the initial-state probabilities `π` (at first capture) are free (`K − 1`
  parameters).

The model names of Palomar et al. (2023, Table 1) map directly: "state + time" is
`:state_plus_time`, "state * time" `:state_time`, "(.)" `:constant`.
"""
struct MultiEventStructure
    nstates::Int
    assignment::Matrix{Bool}
    baseline::Vector{Int}
    survival::Symbol
    transition::Symbol
    detection::Symbol
    state_assignment::Symbol
end

function MultiEventStructure(; nstates::Integer, assignment::AbstractMatrix{Bool},
                             baseline::AbstractVector{<:Integer}, survival::Symbol = :state,
                             transition::Symbol = :state, detection::Symbol = :state,
                             state_assignment::Symbol = :state)
    K = Int(nstates)
    size(assignment, 1) == K || throw(DimensionMismatch("assignment must have $K rows"))
    length(baseline) == K || throw(DimensionMismatch("baseline must have length $K"))
    all(k -> assignment[k, baseline[k]], 1:K) ||
        throw(ArgumentError("each baseline cell must be an allowed assignment"))
    for s in (survival, transition, detection)
        s in MULTIEVENT_STRUCTURES ||
            throw(ArgumentError("structure :$s not in $MULTIEVENT_STRUCTURES"))
    end
    state_assignment in (:constant, :state) ||
        throw(ArgumentError("state_assignment must be :constant or :state"))
    return MultiEventStructure(K, Matrix{Bool}(assignment), collect(Int, baseline), survival,
                               transition, detection, state_assignment)
end

# Design matrix for ncells × ntimes linear predictors, rows ordered (cell, time).
function _ms_design(ncells::Int, ntimes::Int, s::Symbol)
    n = ncells * ntimes
    row(c, t) = (c - 1) * ntimes + t
    if s === :constant
        return ones(n, 1)
    elseif s === :state
        X = zeros(n, ncells)
        for c in 1:ncells, t in 1:ntimes
            X[row(c, t), c] = 1
        end
        return X
    elseif s === :time
        X = zeros(n, ntimes)
        for c in 1:ncells, t in 1:ntimes
            X[row(c, t), t] = 1
        end
        return X
    elseif s === :state_plus_time
        X = zeros(n, ncells + ntimes - 1)
        for c in 1:ncells, t in 1:ntimes
            X[row(c, t), c] = 1
            t > 1 && (X[row(c, t), ncells+t-1] = 1)
        end
        return X
    elseif s === :state_time
        return Matrix(1.0 * LinearAlgebra.I(n))
    end
    throw(ArgumentError("unknown structure :$s"))
end

# Free (non-baseline) assignment cells, row by row.
_ms_free_cells(st::MultiEventStructure) =
    [(k, e) for k in 1:st.nstates for e in axes(st.assignment, 2)
     if st.assignment[k, e] && e != st.baseline[k]]

"""
    MultiEventProblem(histories, structure::MultiEventStructure; unassessed = Int[])

Maximum-likelihood problem for a multi-event (or multistate CJS) model. `histories` is an
`n×T` matrix of event codes (`0` not seen, `1…E` seen events; `missing` for an occasion
not surveyed), one row per marked individual. Rows are conditioned on first capture;
individuals never seen, or first seen on the last occasion, carry no information and are
dropped, and identical histories are pooled with counts. `unassessed` lists occasions at
which no animal's state could be assessed (for example, years without swabbing): there an
encountered animal gives its state's `baseline` event with probability one, whatever the
state-assignment parameters.

It is a LogDensityProblems target whose log density is the log-likelihood (a flat prior on
the unconstrained parameters, in the order of `prob.names`). `solve(prob)` maximises it
(BFGS with finite-difference gradients from `x0 = zeros`) and returns a
[`MultiEventFit`](@ref). [`multievent_parameters`](@ref) gives the natural-scale
parameters at a parameter vector; [`simulate_capture_histories`](@ref) draws data.
"""
struct MultiEventProblem
    structure::MultiEventStructure
    histories::Matrix{Union{Missing,Int}}
    counts::Vector{Int}
    first::Vector{Int}
    designs::NamedTuple
    blocks::NamedTuple
    names::Vector{Symbol}
    unassessed::Vector{Bool}
end

function MultiEventProblem(histories::AbstractMatrix, st::MultiEventStructure;
                           unassessed::AbstractVector{<:Integer} = Int[])
    n, T = size(histories)
    T >= 2 || throw(ArgumentError("need at least two occasions"))
    E = size(st.assignment, 2)
    for h in histories
        h === missing || 0 <= h <= E || throw(ArgumentError("event code $h outside 0…$E"))
    end
    pooled = Dict{Vector{Union{Missing,Int}},Int}()
    order = Vector{Vector{Union{Missing,Int}}}()
    for i in 1:n
        h = Vector{Union{Missing,Int}}(histories[i, :])
        f = findfirst(x -> x !== missing && x > 0, h)
        (f === nothing || f == T) && continue
        haskey(pooled, h) || push!(order, h)
        pooled[h] = get(pooled, h, 0) + 1
    end
    H = Matrix{Union{Missing,Int}}(undef, length(order), T)
    for (r, h) in enumerate(order)
        H[r, :] .= h
    end
    first = [findfirst(x -> x !== missing && x > 0, H[r, :]) for r in 1:size(H, 1)]
    K = st.nstates
    designs = (initial = Matrix(1.0 * LinearAlgebra.I(K - 1)),
               survival = _ms_design(K, T - 1, st.survival),
               transition = _ms_design(K * (K - 1), T - 1, st.transition),
               detection = _ms_design(K, T - 1, st.detection),
               assignment = _ms_design(length(_ms_free_cells(st)), 1, st.state_assignment))
    names = Symbol[]
    blocks = Dict{Symbol,UnitRange{Int}}()
    for (cls, X) in pairs(designs)
        r = length(names) .+ (1:size(X, 2))
        blocks[cls] = r
        append!(names, [Symbol(cls, "_", j) for j in 1:size(X, 2)])
    end
    nb = NamedTuple{keys(designs)}(Tuple(blocks[c] for c in keys(designs)))
    all(t -> 1 <= t <= T, unassessed) || throw(ArgumentError("unassessed occasions outside 1…$T"))
    return MultiEventProblem(st, H, [pooled[h] for h in order], first, designs, nb, names,
                             [t in unassessed for t in 1:T])
end

LogDensityProblems.dimension(p::MultiEventProblem) = length(p.names)
LogDensityProblems.capabilities(::Type{<:MultiEventProblem}) = LogDensityProblems.LogDensityOrder{0}()

# softmax with an implicit zero for the baseline
function _ms_softmax(η::AbstractVector)
    m = max(zero(eltype(η)), maximum(η; init = zero(eltype(η))))
    e = exp.(η .- m)
    z = exp(-m) + sum(e)
    return e ./ z, exp(-m) / z
end

"""
    multievent_parameters(prob::MultiEventProblem, θ) -> NamedTuple

Natural-scale parameters of a [`MultiEventProblem`](@ref) at the unconstrained vector `θ`:
`initial` (length `K`, the state distribution at first capture), `survival` (`K×(T − 1)`),
`transition` (`K×K×(T − 1)`, rows sum to one), `detection` (`K×(T − 1)`, column `t` is
occasion `t + 1`) and `assignment` (`K×E`, rows sum to one).
"""
function multievent_parameters(prob::MultiEventProblem, θ::AbstractVector)
    st = prob.structure
    K = st.nstates
    T = size(prob.histories, 2)
    R = eltype(θ)
    lin(c) = prob.designs[c] * view(θ, prob.blocks[c])
    ηπ = lin(:initial)
    w, w0 = _ms_softmax(ηπ)
    initial = vcat(w, w0)
    S = reshape(_logistic.(lin(:survival)), T - 1, K)'
    p = reshape(_logistic.(lin(:detection)), T - 1, K)'
    ηψ = reshape(lin(:transition), T - 1, K * (K - 1))
    Ψ = zeros(R, K, K, T - 1)
    for t in 1:T-1
        c = 0
        for k in 1:K
            others = [l for l in 1:K if l != k]
            η = [ηψ[t, c+i] for i in eachindex(others)]
            c += length(others)
            pr, stay = _ms_softmax(η)
            Ψ[k, k, t] = stay
            for (i, l) in enumerate(others)
                Ψ[k, l, t] = pr[i]
            end
        end
    end
    E = size(st.assignment, 2)
    δ = zeros(R, K, E)
    cells = _ms_free_cells(st)
    ηδ = lin(:assignment)
    for k in 1:K
        idx = [i for (i, c) in enumerate(cells) if c[1] == k]
        pr, base = _ms_softmax([ηδ[i] for i in idx])
        δ[k, st.baseline[k]] = base
        for (j, i) in enumerate(idx)
            δ[k, cells[i][2]] = pr[j]
        end
    end
    return (initial = initial, survival = Matrix(S), transition = Ψ, detection = Matrix(p),
            assignment = δ)
end

# state assignment at an occasion when no state can be assessed: the baseline event
function _ms_unassessed(st::MultiEventStructure, R::Type = Float64)
    δ = zeros(R, st.nstates, size(st.assignment, 2))
    for k in 1:st.nstates
        δ[k, st.baseline[k]] = one(R)
    end
    return δ
end

# log-likelihood of each pooled history (not multiplied by counts)
function _ms_logliks(prob::MultiEventProblem, θ::AbstractVector)
    par = multievent_parameters(prob, θ)
    K = prob.structure.nstates
    T = size(prob.histories, 2)
    Φ = [survival_transition(view(par.survival, :, t), view(par.transition, :, :, t))
         for t in 1:T-1]
    δ0 = _ms_unassessed(prob.structure, eltype(par.assignment))
    B = [multievent_emission(view(par.detection, :, t - 1),
                             prob.unassessed[t] ? δ0 : par.assignment) for t in 2:T]
    init = vcat(par.initial, zero(eltype(par.initial)))
    R = eltype(B[1])
    out = Vector{R}(undef, size(prob.histories, 1))
    for r in axes(prob.histories, 1)
        f = prob.first[r]
        L = ones(R, T - f + 1, K + 1)
        e = prob.histories[r, f]
        δf = prob.unassessed[f] ? δ0 : par.assignment
        for k in 1:K
            L[1, k] = δf[k, e]
        end
        L[1, K+1] = 0
        for t in f+1:T
            o = prob.histories[r, t]
            o === missing && continue
            L[t-f+1, :] .= view(B[t-1], :, o + 1)
        end
        out[r] = hmm_loglikelihood(init, Φ[f:T-1], L)
    end
    return out
end

"""
    pointwise_loglikelihood(prob::MultiEventProblem, θ) -> Vector

Log-likelihood of each distinct capture history at `θ` (multiply by `prob.counts` and sum
for the total).
"""
pointwise_loglikelihood(prob::MultiEventProblem, θ::AbstractVector) = _ms_logliks(prob, θ)

function LogDensityProblems.logdensity(prob::MultiEventProblem, θ::AbstractVector)
    ll = sum(prob.counts .* _ms_logliks(prob, θ))
    return isnan(ll) ? oftype(ll, -Inf) : ll
end

"""
    simulate_capture_histories(rng, structure, parameters, releases; assignment = nothing)
        -> Matrix{Int}

Draw capture histories from a multi-event model. `parameters` is a named tuple as returned by
[`multievent_parameters`](@ref) (`initial`, `survival`, `transition`, `detection`,
`assignment`); `releases[t]` animals are first captured at occasion `t`, in a state drawn
from `initial` and with a first event drawn from the state-assignment matrix. Pass
`assignment`, a vector of `T` matrices, to let state assignment differ between occasions
(for example, no swabbing in some years). Returns one row per animal of event codes (0 not
seen).
"""
function simulate_capture_histories(rng::AbstractRNG, st::MultiEventStructure,
                                    par::NamedTuple, releases::AbstractVector{<:Integer};
                                    assignment = nothing)
    T = length(releases)
    K = st.nstates
    n = sum(releases)
    H = zeros(Int, n, T)
    i = 0
    for f in 1:T, _ in 1:releases[f]
        i += 1
        δ(t) = assignment === nothing ? par.assignment : assignment[t]
        z = _hmm_draw(rng, par.initial)
        H[i, f] = _hmm_draw(rng, view(δ(f), z, :))
        for t in f+1:T
            Φ = survival_transition(view(par.survival, :, t - 1), view(par.transition, :, :, t - 1))
            z = _hmm_draw(rng, view(Φ, z, :))
            z == K + 1 && break
            B = multievent_emission(view(par.detection, :, t - 1), δ(t))
            H[i, t] = _hmm_draw(rng, view(B, z, :)) - 1
        end
    end
    return H
end

# --- maximum likelihood --------------------------------------------------------------------------

function _ms_grad(f, x::AbstractVector; h::Real = 1e-5)
    g = similar(x, Float64)
    xp = collect(Float64, x)
    for i in eachindex(xp)
        hi = h * max(1.0, abs(x[i]))
        xp[i] = x[i] + hi
        fp = f(xp)
        xp[i] = x[i] - hi
        fm = f(xp)
        xp[i] = x[i]
        g[i] = (fp - fm) / (2hi)
    end
    return g
end

# BFGS minimisation of f with central-difference gradients and Armijo backtracking.
function _ms_bfgs(f, x0::AbstractVector; iters::Integer = 1000, gtol::Real = 1e-4,
                  h::Real = 1e-5)
    x = collect(Float64, x0)
    fx = f(x)
    isfinite(fx) || throw(ArgumentError("the objective is not finite at the starting point"))
    g = _ms_grad(f, x; h)
    d = length(x)
    H = Matrix{Float64}(LinearAlgebra.I, d, d)
    it = 0
    for k in 1:iters
        it = k
        maximum(abs, g) < gtol && break
        dir = -(H * g)
        slope = dot(dir, g)
        if !(slope < 0)
            H = Matrix{Float64}(LinearAlgebra.I, d, d)
            dir = -g
            slope = -dot(g, g)
        end
        t = 1.0
        xn, fn = x, fx
        while t > 1e-14
            xn = x .+ t .* dir
            fn = f(xn)
            isfinite(fn) && fn <= fx + 1e-4 * t * slope && break
            t /= 2
        end
        t <= 1e-14 && break
        gn = _ms_grad(f, xn; h)
        s = xn .- x
        y = gn .- g
        sy = dot(s, y)
        if sy > 1e-12
            k == 1 && (H .*= sy / dot(y, y))
            ρ = 1 / sy
            Hy = H * y
            H .+= (ρ^2 * dot(y, Hy) + ρ) .* (s * s') .- ρ .* (Hy * s' .+ s * Hy')
        end
        x, fx, g = xn, fn, gn
    end
    return (x = x, f = fx, gradient = g, iterations = it)
end

# Central-difference Hessian of f (steps h·max(1, |x|), large enough that rounding in a
# log-likelihood of order 10³ does not swamp the second differences).
function _ms_hessian(f, x::AbstractVector; h::Real = 1e-3)
    d = length(x)
    hs = [h * max(1.0, abs(v)) for v in x]
    fmid = f(x)
    H = zeros(d, d)
    e(i) = (v = zeros(d); v[i] = hs[i]; v)
    for i in 1:d
        H[i, i] = (f(x .+ e(i)) - 2 * fmid + f(x .- e(i))) / hs[i]^2
        for j in 1:i-1
            H[i, j] = H[j, i] = (f(x .+ e(i) .+ e(j)) - f(x .+ e(i) .- e(j)) -
                                 f(x .- e(i) .+ e(j)) + f(x .- e(i) .- e(j))) / (4hs[i] * hs[j])
        end
    end
    return H
end

# inverse of a negative-definite Hessian of a log-likelihood, or NaNs when it is singular
function _ms_vcov(H::AbstractMatrix)
    F = LinearAlgebra.eigen(Symmetric(-H))
    if minimum(F.values) <= 1e-8 * max(1.0, maximum(abs, F.values))
        return fill(NaN, size(H)), count(<=(1e-8 * max(1.0, maximum(abs, F.values))), F.values)
    end
    return Matrix(F.vectors * Diagonal(1 ./ F.values) * F.vectors'), 0
end

"""
    MultiEventFit

Maximum-likelihood fit of a [`MultiEventProblem`](@ref) from `solve`: `problem`, `θ`
(unconstrained estimates in the order of `problem.names`), `vcov` (inverse observed
information on that scale; `NaN` when singular), `rank_deficiency` (the number of
non-identifiable directions, as E-SURGE reports), `loglik`, `deviance` (`−2 loglik`), `k` (the
number of estimable parameters, `length(θ) − rank_deficiency`), `aic`, `parameters` (from
[`multievent_parameters`](@ref)) and `gradient` (at the optimum, to check convergence). Use
[`qaic`](@ref) for overdispersed data.
"""
struct MultiEventFit
    problem::MultiEventProblem
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

# solve(prob::MultiEventProblem; x0 = zeros(dimension), iters = 2000) -> MultiEventFit.
# Documented in MultiEventProblem.
function CommonSolve.solve(prob::MultiEventProblem;
                           x0::AbstractVector = zeros(LogDensityProblems.dimension(prob)),
                           iters::Integer = 2000)
    f(θ) = -LogDensityProblems.logdensity(prob, θ)
    r = _ms_bfgs(f, x0; iters)
    V, def = _ms_vcov(-_ms_hessian(f, r.x))
    ll = -r.f
    k = length(r.x) - def
    return MultiEventFit(prob, r.x, V, def, ll, -2ll, k, -2ll + 2k,
                         multievent_parameters(prob, r.x), r.gradient)
end

"""
    qaic(fit, ĉ) -> Float64

Quasi-AIC `deviance / ĉ + 2k` of a [`MultiEventFit`](@ref) with overdispersion `ĉ`
(Burnham & Anderson 2002; Palomar et al. 2023 used ĉ = 1.474 from U-CARE).
"""
qaic(fit::MultiEventFit, ĉ::Real) = fit.deviance / ĉ + 2fit.k
