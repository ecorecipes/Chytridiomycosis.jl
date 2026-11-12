# NEW MODEL (N5, not a published analysis): an observation model linking the Briggs et al.
# (2010) load and zoospore-pool states to skin-swab qPCR (per frog, or pooled) and to eDNA
# in pond water, and an identifiability analysis of the Briggs parameters from these data.
# Wiki: lab/wiki/models/new-swab-edna-identifiability.md
#
# Latent state (one pond of N frogs, linear phase: every load well below Smax). Frogs start
# free of Bd and the pool is seeded with Z0 zoospores. Counting sporangia as a linear
# birth–death–immigration process given the pool (N1, src/novel/load_bdi.jl), with
# ι = νγZ/V, λ = ηνf and r0 = λ − σ:
#   dm/dt = ι + r0 m                       mean sporangia per frog (the Briggs equation)
#   dZ/dt = η(1−f) N m − (μ + γN/V) Z      zoospores in the pool
#   dw/dt = 2λ m + 2 r0 w                  w = Var − mean of sporangia among frogs
# The system is linear, so it is solved exactly with a matrix exponential. The aggregation of
# loads among frogs is k(t) = m²/w (the negative-binomial shape; exactly ι/λ at constant ι).
#
# Observations (units in docstrings):
#   swab tube of n pooled swabs ~ NegativeBinomial(shape n k, mean n ρ m)     (count of copies)
#   eDNA sample                  ~ LogNormal(log(κ Z/V), σe), censored below 1 copy
# ρ (copies per swab per sporangium) and κ (copies per sample per zoospore per volume unit)
# are the assay scales. The NB swab likelihood is exact when ρ ≤ 1 is a recovery fraction
# (binomial thinning keeps k, N1) and an approximation otherwise. qPCR measurement noise is
# not modelled separately, so it would be absorbed into k.

"""
    BriggsMomentProblem(p::BriggsParameters; N, Z0, m0 = 0.0, w0 = 0.0, tspan = (0.0, 60.0))

Moments of the sporangial load on `N` frogs sharing one pool, in the linear (sub-lethal)
phase of the Briggs et al. (2010) model with each frog's load a birth–death–immigration
process given the pool (see [`LoadBDI`](@ref)): the mean `m` (sporangia per frog), the pool
`Z` (zoospores) and `w = Var − m`, the excess variance of load among frogs. `Z0` zoospores
start in the pool; `m0`, `w0` set the initial moments (default: frogs free of Bd).

`solve(prob; saveat)` is exact (a matrix exponential) and returns a `NamedTuple` with
vectors `t`, `m`, `Z`, `w` and the aggregation `k = m²/w`, the negative-binomial shape of
loads among frogs. `Smax` is ignored: use it only while loads stay well below `Smax`.
This is new work (model N5), not a published result.
"""
struct BriggsMomentProblem{P<:BriggsParameters,T}
    p::P
    N::Float64
    u0::Vector{Float64}
    tspan::T
end
function BriggsMomentProblem(p::BriggsParameters; N::Real, Z0::Real, m0::Real = 0.0,
        w0::Real = 0.0, tspan = (0.0, 60.0))
    N > 0 || throw(ArgumentError("N must be positive"))
    (Z0 >= 0 && m0 >= 0) || throw(ArgumentError("Z0 and m0 must be non-negative"))
    return BriggsMomentProblem(p, float(N), [float(m0), float(Z0), float(w0)], float.(tspan))
end

# matrix of the linear moment system, state (m, Z, w). `p` is a BriggsParameters or any
# object with the same fields (a NamedTuple of complex values for complex-step derivatives)
function briggs_moment_matrix(p, N)
    λ = p.η * p.ν * p.f
    r0 = λ - p.σ
    b = p.ν * p.γ / p.V
    z = zero(r0)
    return [r0                 b                       z
            p.η*(1-p.f)*N      -(p.μ + p.γ * N / p.V)  z
            2λ                 z                       2r0]
end

# exact moments (m, Z, w) at each of `days` from u0 at day 0
briggs_moments(p, N, u0, days) =
    [exp(briggs_moment_matrix(p, N) * t) * u0 for t in days]

function moment_rhs!(du, u, prob::BriggsMomentProblem, t)
    du .= briggs_moment_matrix(prob.p, prob.N) * u
    return nothing
end

to_ode_problem(prob::BriggsMomentProblem) = ODEProblem(moment_rhs!, copy(prob.u0), prob.tspan, prob)

function CommonSolve.solve(prob::BriggsMomentProblem;
        saveat = range(prob.tspan[1], prob.tspan[2]; length = 101))
    t = collect(Float64, saveat)
    xs = [exp(briggs_moment_matrix(prob.p, prob.N) * (s - prob.tspan[1])) * prob.u0 for s in t]
    m, Z, w = getindex.(xs, 1), getindex.(xs, 2), getindex.(xs, 3)
    return (t = t, m = m, Z = Z, w = w, k = m .^ 2 ./ w)
end

"""
    SurveyPond(; N, days, swabs = 0, pool = 1, edna = 0)

Sampling design for one pond of `N` frogs: on each of `days` (after the pool is seeded on
day 0), `swabs` qPCR results from skin swabs, each from `pool` frogs' swabs combined in one
tube (`pool = 1`: individual swabs; `swabs * pool ≤ N`), and `edna` water samples.
"""
struct SurveyPond
    N::Int
    days::Vector{Float64}
    swabs::Int
    pool::Int
    edna::Int
end
function SurveyPond(; N::Integer, days, swabs::Integer = 0, pool::Integer = 1,
        edna::Integer = 0)
    N > 0 || throw(ArgumentError("N must be positive"))
    (swabs >= 0 && edna >= 0 && pool >= 1) || throw(ArgumentError("invalid sample numbers"))
    swabs * pool <= N || throw(ArgumentError("swabs * pool must not exceed N"))
    d = collect(Float64, days)
    (issorted(d) && all(>(0), d)) || throw(ArgumentError("days must be sorted and positive"))
    return SurveyPond(N, d, swabs, pool, edna)
end

"""
    SurveySample(pond, day, kind, pool, copies)

One qPCR result: pond index, `day`, `kind` (`:swab` or `:edna`), the number of swabs `pool`
combined in the tube (1 for eDNA), and `copies` (Bd ITS copies per swab tube, or per eDNA
sample; eDNA values below 1 copy are non-detections).
"""
struct SurveySample
    pond::Int
    day::Float64
    kind::Symbol
    pool::Int
    copies::Float64
end

"""
    SurveyObservation(; ρ, κ, σe)

Assay scales and noise: a swab carries on average `ρ` Bd ITS copies per sporangium on the frog
(copies sporangium⁻¹); an eDNA sample carries on average `κ` copies per unit of zoospore
density `Z/V` (copies per sample per zoospore volume-unit⁻¹, i.e. copies per zoospore times the
volume filtered); `σe` is the standard deviation of log eDNA copies. None of these has a
published value for the Briggs model: they must be supplied (or fitted), and stated.
"""
Base.@kwdef struct SurveyObservation{T<:Real}
    ρ::T
    κ::T
    σe::T
end

"""
    simulate_survey(rng, p::BriggsParameters, obs::SurveyObservation, ponds; Z0)
        -> Vector{SurveySample}

Synthetic data for the designs `ponds` (a vector of [`SurveyPond`](@ref)), with `Z0[i]`
zoospores seeded into pond `i` on day 0. Swab tubes are negative binomial with shape `pool·k`
and mean `pool·ρ·m` (the sum of `pool` frogs); eDNA is lognormal around `κZ/V`, with values
below 1 copy recorded as 0. For identifiability studies only: the output of a new model.
"""
function simulate_survey(rng::AbstractRNG, p::BriggsParameters, obs::SurveyObservation,
        ponds::AbstractVector{SurveyPond}; Z0)
    length(Z0) == length(ponds) || throw(ArgumentError("need one Z0 per pond"))
    out = SurveySample[]
    for (i, pd) in enumerate(ponds)
        xs = briggs_moments(p, pd.N, [0.0, float(Z0[i]), 0.0], pd.days)
        for (t, x) in zip(pd.days, xs)
            m, Z, w = x
            for _ in 1:pd.swabs
                y = rand(rng, survey_swab_distribution(m, w, obs.ρ, pd.pool))
                push!(out, SurveySample(i, t, :swab, pd.pool, y))
            end
            for _ in 1:pd.edna
                y = exp(log(obs.κ * Z / p.V) + obs.σe * randn(rng))
                push!(out, SurveySample(i, t, :edna, 1, y >= 1 ? y : 0.0))
            end
        end
    end
    return out
end

# sum of `n` swabs from frogs with load mean m and excess variance w, at assay scale ρ
function survey_swab_distribution(m, w, ρ, n)
    k = m^2 / w
    M = n * ρ * m
    r = n * k
    return NegativeBinomial(r, r / (r + M))
end

const SURVEY_BRIGGS = (:γ, :ν, :f, :μ, :η, :σ)
const SURVEY_OBS = (:ρ, :κ, :σe)

"""
    SurveyModel(samples, ponds; p::BriggsParameters, obs::SurveyObservation, Z0,
        free = (:γ, :ν, :f, :μ, :η, :σ, :ρ, :κ, :σe, :Z0_1, …))

Likelihood of swab and eDNA `samples` from the `ponds` designs (see [`simulate_survey`](@ref))
under the Briggs moment model (new model N5). The parameters are the Briggs rates `γ, ν, f,
μ, η, σ` (`V` and `Smax` are fixed from `p`), the assay scales `ρ, κ, σe` and the seeded pool
`Z0_i` of each pond. Those in `free` are estimated on the log scale; the rest are fixed at
the values in `p`, `obs` and `Z0`, which also give the starting point. Fixing `κ` means an
eDNA assay on an absolute scale, fixing `ρ` a swab assay calibrated in sporangia.

Implements LogDensityProblems (log-likelihood of the log parameters, flat prior, in the order
of [`parameter_names`](@ref)); fit with [`maximum_likelihood`](@ref) and
[`profile_likelihood`](@ref); check identifiability with [`sensitivity_rank`](@ref),
[`expected_information`](@ref) and [`observed_information`](@ref).
"""
struct SurveyModel{P<:BriggsParameters}
    samples::Vector{SurveySample}
    ponds::Vector{SurveyPond}
    p::P
    names::Tuple
    values::Vector{Float64}     # natural-scale values of every parameter, in `names` order
    free::Vector{Int}           # indices into `names` of the estimated parameters
end
function SurveyModel(samples::AbstractVector{SurveySample}, ponds::AbstractVector{SurveyPond};
        p::BriggsParameters, obs::SurveyObservation, Z0, free = nothing)
    length(Z0) == length(ponds) || throw(ArgumentError("need one Z0 per pond"))
    all(s -> 1 <= s.pond <= length(ponds), samples) ||
        throw(ArgumentError("sample refers to a pond that is not in ponds"))
    names = (SURVEY_BRIGGS..., SURVEY_OBS..., (Symbol("Z0_", i) for i in eachindex(ponds))...)
    values = Float64[p.γ, p.ν, p.f, p.μ, p.η, p.σ, obs.ρ, obs.κ, obs.σe, Z0...]
    free = free === nothing ? names : Tuple(free)
    idx = Int[]
    for n in free
        j = findfirst(==(n), names)
        j === nothing && throw(ArgumentError("$n is not a parameter; choose from $names"))
        push!(idx, j)
    end
    allunique(idx) || throw(ArgumentError("free has duplicates"))
    all(>(0), values[idx]) || throw(ArgumentError("free parameters must start positive"))
    return SurveyModel(collect(samples), collect(ponds), p, names, values, idx)
end

"""
    parameter_names(m::SurveyModel)

Names of the estimated (free) parameters of a [`SurveyModel`](@ref), all on the log scale.
"""
parameter_names(m::SurveyModel) = m.names[m.free]

LogDensityProblems.dimension(m::SurveyModel) = length(m.free)
LogDensityProblems.capabilities(::Type{<:SurveyModel}) = LogDensityProblems.LogDensityOrder{0}()

# all parameters (natural scale) with the free ones set from log values x; generic in eltype
function survey_values(m::SurveyModel, x::AbstractVector)
    v = Vector{promote_type(Float64, eltype(x))}(m.values)
    v[m.free] .= exp.(x)
    return v
end

survey_start(m::SurveyModel) = log.(m.values[m.free])

survey_briggs(m::SurveyModel, v) =
    (V = m.p.V, γ = v[1], ν = v[2], f = v[3], μ = v[4], η = v[5], σ = v[6])

function LogDensityProblems.logdensity(m::SurveyModel, x::AbstractVector)
    length(x) == LogDensityProblems.dimension(m) ||
        throw(DimensionMismatch("expected $(LogDensityProblems.dimension(m)) parameters"))
    all(isfinite, x) || return -Inf
    v = survey_values(m, x)
    v[3] < 1 || return -Inf                                   # f is a fraction
    p = survey_briggs(m, v)
    ρ, κ, σe = v[7], v[8], v[9]
    moments = [briggs_moments(p, pd.N, [0.0, v[9 + i], 0.0], pd.days)
               for (i, pd) in enumerate(m.ponds)]
    ll = 0.0
    for s in m.samples
        pd = m.ponds[s.pond]
        j = findfirst(==(s.day), pd.days)
        j === nothing && continue
        mm, Z, w = moments[s.pond][j]
        (mm > 0 && Z > 0 && w > 0) || return -Inf
        if s.kind === :swab
            ll += logpdf(survey_swab_distribution(mm, w, ρ, s.pool), round(Int, s.copies))
        else
            ll += knapp_obs(s.copies, log(κ * Z / p.V), σe)
        end
    end
    return isfinite(ll) ? ll : -Inf
end

"""
    maximum_likelihood(m::SurveyModel; starts = 3, rounds = 3, iters = 2000, spread = 1.0,
        x0 = log of the starting values, rng = Random.Xoshiro(1))

Maximum likelihood fit of a [`SurveyModel`](@ref) by Nelder–Mead on the log parameters, as
for the Knapp fit but preconditioned: the simplex moves in coordinates whitened by the
[`expected_information`](@ref) (one unit = one standard error along each of its
eigenvectors), recomputed at the current best point for each of `rounds`. Without this the
simplex stalls in the narrow curved ridges of this likelihood. Restarts from `x0` and from
`starts − 1` points `spread` standard errors around it. Returns a `NamedTuple` `model`, `x`,
`estimates` (natural scale, by name), `loglik`, `aic` and `logliks` (best of each start).
"""
function maximum_likelihood(m::SurveyModel; starts::Integer = 3, rounds::Integer = 3,
        iters::Integer = 2000, spread::Real = 1.0, x0::AbstractVector = survey_start(m),
        rng::AbstractRNG = Random.Xoshiro(1))
    f(x) = -LogDensityProblems.logdensity(m, x)
    isfinite(f(x0)) || throw(ArgumentError("the likelihood cannot be evaluated at x0"))
    best, bestf = collect(float.(x0)), f(x0)
    logliks = Float64[]
    for s in 1:starts
        x = collect(float.(x0))
        if s > 1
            T = survey_whitening(m, x0)
            for _ in 1:1000
                x = x0 .+ T * (spread .* randn(rng, length(x0)))
                isfinite(f(x)) && break
            end
        end
        fx = f(x)
        for _ in 1:rounds
            T = survey_whitening(m, x)
            z, fz = nelder_mead(z -> f(x .+ T * z), zeros(length(x)); iters = iters, step = 1.0)
            z, fz = nelder_mead(z -> f(x .+ T * z), z; iters = iters ÷ 2, step = 0.1)
            fz < fx && ((x, fx) = (x .+ T * z, fz))
        end
        push!(logliks, -fx)
        fx < bestf && ((best, bestf) = (x, fx))
    end
    names = parameter_names(m)
    return (model = m, x = best, estimates = NamedTuple{names}(Tuple(exp.(best))),
        loglik = -bestf, aic = 2length(best) + 2bestf, logliks = logliks)
end

# x = x0 + T z with T = V Λ^{-1/2} from the expected information at x (eigenvalues floored,
# so unidentified directions get large but finite steps)
function survey_whitening(m::SurveyModel, x::AbstractVector)
    E = eigen(Symmetric(expected_information(m, x)))
    λ = max.(E.values, 1e-8 * maximum(E.values))
    return E.vectors * Diagonal(1 ./ sqrt.(λ))
end

"""
    profile_likelihood(m::SurveyModel, name::Symbol, values; x0 = log of the starting values,
        kwargs...)

Profile log-likelihood of parameter `name` of a [`SurveyModel`](@ref) at each of `values`
(natural scale): the model is refitted by [`maximum_likelihood`](@ref) (keywords passed on)
with `name` fixed at each value and the other free parameters re-estimated, starting from
`x0` (pass `fit.x` of the full fit). Values within 1.92 of the maximum log-likelihood form an
approximate 95% confidence interval.
"""
function profile_likelihood(m::SurveyModel, name::Symbol, values;
        x0::AbstractVector = survey_start(m), kwargs...)
    j = findfirst(==(name), parameter_names(m))
    j === nothing && throw(ArgumentError("$name is not a free parameter of this model"))
    others = setdiff(1:length(m.free), j)
    return map(collect(Float64, values)) do val
        vals = copy(m.values)
        vals[m.free] .= exp.(x0)
        vals[m.free[j]] = val
        mj = SurveyModel(m.samples, m.ponds, m.p, m.names, vals, m.free[others])
        maximum_likelihood(mj; starts = 1, kwargs...).loglik
    end
end

"""
    observed_information(m::SurveyModel, x = log of the starting values; h = 1e-3)

Observed Fisher information (minus the Hessian of the log-likelihood, central differences
with step `h`) at log parameters `x`. Its inverse approximates the covariance of the log
estimates, so `sqrt.(diag(inv(I)))` are standard errors on the log scale (≈ coefficients
of variation); a near-singular matrix flags a practically unidentified direction.
"""
function observed_information(m::SurveyModel, x::AbstractVector = survey_start(m);
        h::Real = 1e-3)
    f(y) = LogDensityProblems.logdensity(m, y)
    n = length(x)
    H = zeros(n, n)
    e(i) = (u = zeros(n); u[i] = h; u)
    for i in 1:n, j in i:n
        ei, ej = e(i), e(j)
        H[i, j] = (f(x + ei + ej) - f(x + ei - ej) - f(x - ei + ej) + f(x - ei - ej)) / (4h^2)
        H[j, i] = H[i, j]
    end
    return -H
end

# log of the three summaries that fix the sample distributions at every pond and day, in
# pond-then-day order: swab mean ρm, aggregation k = m²/w and eDNA level κZ/V (generic eltype)
function survey_summaries(m::SurveyModel, x::AbstractVector)
    v = survey_values(m, x)
    p = survey_briggs(m, v)
    out = NTuple{3,eltype(v)}[]
    for (i, pd) in enumerate(m.ponds)
        for (mm, Z, w) in briggs_moments(p, pd.N, [zero(v[1]), v[9 + i], zero(v[1])], pd.days)
            push!(out, (log(v[7] * mm), log(mm^2 / w), log(v[8] * Z / p.V)))
        end
    end
    return out
end

# Jacobians of the three log summaries with respect to the log free parameters, by complex
# steps (exact to rounding): a 3-vector of matrices (rows: pond and day)
function survey_jacobians(m::SurveyModel, x::AbstractVector)
    n, h = length(x), 1e-30
    cols = map(1:n) do j
        y = complex.(float.(x))
        y[j] += im * h
        [imag.(t) ./ h for t in survey_summaries(m, y)]
    end
    return [reduce(hcat, [getindex.(c, s) for c in cols]) for s in 1:3]
end

# (pond, day) rows in survey_summaries order, and each pond's design
survey_rows(m::SurveyModel) = [(i, pd) for (i, pd) in enumerate(m.ponds) for _ in pd.days]

"""
    sensitivity_rank(m::SurveyModel, x = log of the starting values; aggregation = true,
        rtol = 1e-7)

Local structural identifiability of the free parameters of `m` under its sampling design,
from the rank of the Jacobian of the noise-free outputs (log swab mean `ρm`, log aggregation
`k = m²/w` if `aggregation`, log eDNA level `κZ/V`, at every sampled pond and day) with
respect to the log parameters, computed exactly by complex steps. `aggregation = false`
describes swab totals (pooled samples whose spread is not used). The data are not used.

Returns a `NamedTuple`: `rank`, `nparams`, `singular` (singular values, largest first),
`names`, and `null`, a vector of `NamedTuple`s (parameter name ⇒ log-scale component, the
largest scaled to ±1), one per unidentified direction: moving along it leaves every output
unchanged. The null space is a basis, so with several directions any mixture of them is
also unidentified.
"""
function sensitivity_rank(m::SurveyModel, x::AbstractVector = survey_start(m);
        aggregation::Bool = true, rtol::Real = 1e-7)
    n = length(x)
    Jm, Jk, Je = survey_jacobians(m, x)
    rows = survey_rows(m)
    swab = [pd.swabs > 0 for (_, pd) in rows]
    edna = [pd.edna > 0 for (_, pd) in rows]
    J = vcat(Jm[swab, :], aggregation ? Jk[swab, :] : zeros(0, n), Je[edna, :])
    F = svd(J; full = true)
    s = F.S
    r = count(>(rtol * s[1]), s)
    names = parameter_names(m)
    null = [NamedTuple{names}(Tuple(F.V[:, j] ./ maximum(abs, F.V[:, j]))) for j in (r + 1):n]
    return (rank = r, nparams = n, singular = s, names = names, null = null)
end

# Fisher information for log r of one NegativeBinomial with shape r and mean M (mean and shape
# are orthogonal). The score in r is ψ(y+r) − ψ(r) + log(r/(r+M)) + (M − y)/(r+M), and
# ψ(y+r) − ψ(r) = Σ_{j<y} 1/(r+j), so no special functions are needed.
function nb_log_shape_information(r::Real, M::Real)
    d = NegativeBinomial(r, r / (r + M))
    ymax = ceil(Int, M + 12 * sqrt(M * (1 + M / r)) + 30)
    c = log(r / (r + M)) + M / (r + M)
    acc, ψ = 0.0, 0.0
    for y in 0:ymax
        acc += pdf(d, y) * (ψ + c - y / (r + M))^2
        ψ += 1 / (r + y)
    end
    return r^2 * acc
end

"""
    expected_information(m::SurveyModel, x = log of the starting values; aggregation = true)

Expected Fisher information of the free log parameters of `m` for its sampling design
(`m.ponds`; the samples are not used), at log parameters `x`. Each swab tube contributes
through its negative-binomial mean (information `rM/(r+M)` for log `M`) and, if
`aggregation`, its shape (computed by summation); each eDNA sample through its log level
(`1/σe²`; censoring below one copy is ignored) and log `σe` (2). `aggregation = false` gives
the information in swab totals alone. `sqrt.(diag(inv(I)))` are the smallest attainable
standard errors on the log scale (≈ coefficients of variation): use it to compare designs.
"""
function expected_information(m::SurveyModel, x::AbstractVector = survey_start(m);
        aggregation::Bool = true)
    n = length(x)
    v = survey_values(m, x)
    p = survey_briggs(m, v)
    Jm, Jk, Je = survey_jacobians(m, x)
    I = zeros(n, n)
    row = 0
    for (i, pd) in enumerate(m.ponds)
        for (mm, Z, w) in briggs_moments(p, pd.N, [0.0, v[9 + i], 0.0], pd.days)
            row += 1
            if pd.swabs > 0
                r, M = pd.pool * mm^2 / w, pd.pool * v[7] * mm
                I .+= pd.swabs * (r * M / (r + M)) .* (Jm[row, :] * Jm[row, :]')
                aggregation && (I .+= pd.swabs * nb_log_shape_information(r, M) .*
                                      (Jk[row, :] * Jk[row, :]'))
            end
            pd.edna > 0 && (I .+= (pd.edna / v[9]^2) .* (Je[row, :] * Je[row, :]'))
        end
    end
    j = findfirst(==(9), m.free)                       # log σe
    j === nothing || (I[j, j] += 2 * sum(pd.edna * length(pd.days) for pd in m.ponds))
    return I
end

"""
    identifiable_combinations(p::BriggsParameters)

The combinations of Briggs parameters that the N5 hand derivation shows are identified by
richer and richer data (new work, not a published result):

- `γ_V = γ/V`, `μ` and `r0 = ηνf − σ`: swab means (or eDNA) at two or more host densities,
  from the trace `r0 − μ − γN/V` and determinant `−r0(μ + γN/V) − νγη(1−f)N/V` of the linear
  system; at one density only the trace and determinant are identified;
- `νη_1mf = νη(1−f)`: the same data (slope of the determinant in `N`);
- `ν` and `η(1−f)` separately: add eDNA, with both assays on absolute scales (otherwise only
  `ρν/κ`);
- `ηf`: add the aggregation of loads among frogs, `k = γZ/(Vηf)`, with the pool `Z` on an
  absolute scale; then `η`, `f` and `σ` follow.

See the wiki page `models/new-swab-edna-identifiability.md`.
"""
identifiable_combinations(p::BriggsParameters) =
    (γ_V = p.γ / p.V, μ = p.μ, r0 = p.η * p.ν * p.f - p.σ, νη_1mf = p.ν * p.η * (1 - p.f),
     ηf = p.η * p.f)
