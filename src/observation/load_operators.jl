# Load observation layer, part 4: observation operators. Each maps a model's load state to
# the distribution of what one field swab reports, in ZE per swab:
#
#   prevalence π (infected fraction of the hosts swabbed)
#   × a mixture over components k, weight w_k, of log10 ZE uniform on [lo_k, hi_k]
#     (a point when lo_k = hi_k), plus Gaussian noise of sd σ on the log10 scale
#   → swabs below the limit of detection read as non-detects.
#
# What is published and what is ours:
# - Wilber 2016/2017 IPMs: x = ln(ZE per swab) IS the swab measurement (the growth kernel was
#   fitted to swab loads, so swab noise is inside it), and an uninfected host reads zero. The
#   operator with σ = 0 and lod = 0 only re-expresses the published state. σ > 0 or lod > 0
#   add observation error or censoring that the paper does not have: OUR choice.
# - Louca 2014: the DOI is read in the ESM (§5.1) as the qPCR load of a standardised swab
#   (ZE per swab); the binned DOI n_k is a point component. Same caveat for σ and lod.
# - Briggs 2010: the state is sporangia per host. No source converts sporangia to ZE per swab,
#   so the operator refuses unless ze_per_sporangium is given, and then warns; the result
#   rests entirely on that number. A lognormal spread among hosts around a mean load is OUR
#   choice (the deterministic model has no among-host spread).

"""
    SwabLoadDistribution(prevalence, lower, upper, weights; σ = 0, lod = 0)

The distribution of the load reported by one swab, in ZE per swab:

- with probability `1 − prevalence` the host is uninfected and the swab reads zero;
- otherwise its log10 load is drawn from a mixture: component `k` (weight `weights[k]`) is
  uniform on `[lower[k], upper[k]]` in log10 ZE (a single point when they are equal), plus
  Gaussian noise of standard deviation `σ` on the log10 scale;
- a swab below the limit of detection `lod` (ZE per swab) reads as a non-detect.

Built by [`swab_distribution`](@ref) from a model state. Use
[`detection_probability`](@ref), [`swab_logpdf`](@ref), [`swab_loglikelihood`](@ref),
[`load_moments`](@ref) and `rand(rng, d, n)` on it.
"""
struct SwabLoadDistribution
    prevalence::Float64
    lower::Vector{Float64}
    upper::Vector{Float64}
    weights::Vector{Float64}
    σ::Float64
    lod::Float64
    function SwabLoadDistribution(prevalence::Real, lower::AbstractVector{<:Real},
            upper::AbstractVector{<:Real}, weights::AbstractVector{<:Real}; σ::Real = 0.0,
            lod::Real = 0.0)
        0 <= prevalence <= 1 || throw(ArgumentError("prevalence must lie in [0, 1]"))
        length(lower) == length(upper) == length(weights) ||
            throw(DimensionMismatch("lower, upper and weights differ in length"))
        all(lower .<= upper) || throw(ArgumentError("need lower ≤ upper"))
        all(isfinite, lower) && all(isfinite, upper) ||
            throw(ArgumentError("component bounds must be finite (log10 ZE)"))
        all(>=(0), weights) || throw(ArgumentError("weights must be non-negative"))
        σ >= 0 || throw(ArgumentError("σ must be non-negative"))
        lod >= 0 || throw(ArgumentError("lod must be non-negative"))
        s = sum(weights)
        prevalence > 0 && !(s > 0) && throw(ArgumentError("an infected population needs positive weights"))
        w = s > 0 ? Float64.(weights) ./ s : Float64.(weights)
        return new(Float64(prevalence), Float64.(lower), Float64.(upper), w, Float64(σ), Float64(lod))
    end
end

_lod10(d::SwabLoadDistribution) = d.lod > 0 ? log10(d.lod) : -Inf

# P(z ≥ c) for one component. ∫Φ(u)du = uΦ(u) + φ(u).
function _component_above(lo, hi, σ, c)
    c == -Inf && return 1.0
    if σ == 0
        lo == hi && return lo >= c ? 1.0 : 0.0
        return clamp((hi - max(lo, c)) / (hi - lo), 0.0, 1.0)
    end
    N = Normal()
    lo == hi && return ccdf(N, (c - lo) / σ)
    G(u) = u * cdf(N, u) + pdf(N, u)
    return σ * (G((hi - c) / σ) - G((lo - c) / σ)) / (hi - lo)
end

# Density of z (log10 load) for one component.
function _component_density(lo, hi, σ, z)
    if σ == 0
        lo == hi && throw(ArgumentError("a point component has no density without noise: " *
                                        "give σ > 0"))
        return lo <= z <= hi ? 1 / (hi - lo) : 0.0
    end
    lo == hi && return pdf(Normal(lo, σ), z)
    return (cdf(Normal(), (z - lo) / σ) - cdf(Normal(), (z - hi) / σ)) / (hi - lo)
end

"""
    detection_probability(d::SwabLoadDistribution) -> Float64

Probability that a swab is positive (at or above the limit of detection): the prevalence
times the fraction of infected hosts' swabs at or above `lod`.
"""
detection_probability(d::SwabLoadDistribution) = d.prevalence *
    sum(w * _component_above(lo, hi, d.σ, _lod10(d)) for (lo, hi, w) in zip(d.lower, d.upper, d.weights); init = 0.0)

"""
    swab_logpdf(d::SwabLoadDistribution, y) -> Float64

Log-likelihood of one swab reading `y` (ZE per swab, a number or a [`ZEPerSwab`](@ref)).
A zero or a reading below `d.lod` is a non-detect, with log probability
`log(1 − detection_probability(d))`; a positive reading contributes the density of `y` on the
natural scale (the log10 density divided by `y ln 10`), so likelihoods are comparable with
those of [`LoadHurdleModel`](@ref) on the continuous scale. Point components need `σ > 0`.
"""
function swab_logpdf(d::SwabLoadDistribution, y::Real)
    y >= 0 || throw(DomainError(y, "a swab load cannot be negative"))
    (y == 0 || y < d.lod) && return log1p(-detection_probability(d))
    z = log10(y)
    f = sum(w * _component_density(lo, hi, d.σ, z) for (lo, hi, w) in zip(d.lower, d.upper, d.weights); init = 0.0)
    return log(d.prevalence) + log(f) - log(y * log(10))
end
swab_logpdf(d::SwabLoadDistribution, y::ZEPerSwab) = swab_logpdf(d, load_value(y))

"""
    swab_loglikelihood(d::SwabLoadDistribution, ys) -> Float64

Sum of [`swab_logpdf`](@ref) over the swab readings `ys` (ZE per swab).
"""
swab_loglikelihood(d::SwabLoadDistribution, ys) = sum(y -> swab_logpdf(d, y), ys; init = 0.0)

function Base.rand(rng::AbstractRNG, d::SwabLoadDistribution)
    rand(rng) < d.prevalence || return 0.0
    u, k, acc = rand(rng), length(d.weights), 0.0
    for i in eachindex(d.weights)
        acc += d.weights[i]
        u <= acc && (k = i; break)
    end
    z = d.lower[k] + rand(rng) * (d.upper[k] - d.lower[k]) + d.σ * randn(rng)
    y = exp10(z)
    return y < d.lod ? 0.0 : y
end
Base.rand(rng::AbstractRNG, d::SwabLoadDistribution, n::Integer) = [rand(rng, d) for _ in 1:n]
Base.rand(d::SwabLoadDistribution, n::Integer) = rand(Random.default_rng(), d, n)

# E[10^(j z) 1{z ≥ c}] for a point component at x with noise σ.
function _point_moment(x, σ, c, j)
    κ = log(10)
    σ == 0 && return x >= c ? exp10(j * x) : 0.0
    c == -Inf && return exp10(j * x) * exp((j * σ * κ)^2 / 2)
    return exp10(j * x) * exp((j * σ * κ)^2 / 2) * ccdf(Normal(), (c - x - j * σ^2 * κ) / σ)
end

"""
    load_moments(d::SwabLoadDistribution) -> (; mean, variance)

Mean and variance of the reported swab load (ZE per swab), with uninfected hosts and
non-detects read as zero: the quantities of Taylor's power law ([`TaylorPowerLawProblem`](@ref)).
Uniform components are integrated with a 16-point midpoint rule.
"""
function load_moments(d::SwabLoadDistribution)
    c = _lod10(d)
    m1 = m2 = 0.0
    for (lo, hi, w) in zip(d.lower, d.upper, d.weights)
        if lo == hi
            m1 += w * _point_moment(lo, d.σ, c, 1)
            m2 += w * _point_moment(lo, d.σ, c, 2)
        else
            xs = lo .+ (hi - lo) .* ((1:16) .- 0.5) ./ 16
            m1 += w * sum(_point_moment(x, d.σ, c, 1) for x in xs) / 16
            m2 += w * sum(_point_moment(x, d.σ, c, 2) for x in xs) / 16
        end
    end
    π = d.prevalence
    mean = π * m1
    return (; mean, variance = π * m2 - mean^2)
end

# Operators -----------------------------------------------------------------------------------

"""
    swab_distribution(loads::AbstractVector{ZEPerSwab}; σ = 0, lod = 0)
    swab_distribution(loads::AbstractVector{SporangiaPerHost}; ze_per_sporangium, σ = 0, lod = 0)
    swab_distribution(mean_load::ZEPerSwab; sdlog10, prevalence = 1, lod = 0)
    swab_distribution(sol::BdIPMSolution; step = length(sol.t), σ = 0, lod = 0)
    swab_distribution(sol::HybridIPMSolution; step = length(sol.t), σ = 0, lod = 0)
    swab_distribution(p::LoucaParameters, out; stage = :adult, step = nothing, σ = 0, lod = 0)
        -> SwabLoadDistribution

Observation operators: the distribution of one field swab's reading (ZE per swab) predicted
from a model's load state. `σ` is extra Gaussian noise on the log10 scale and `lod` a limit
of detection (ZE per swab); both default to none.

- **Host loads in ZE per swab** (for example a stochastic or individual-based model): each
  host is a point component; hosts with zero load are uninfected.
- **Sporangia per host** (Briggs et al. 2010): no published conversion exists, so this
  throws an `ArgumentError` unless `ze_per_sporangium` is given, and warns when it is. Pass
  the loads of the living frogs.
- **A mean load** (for a deterministic model with no spread among hosts): a lognormal among
  infected hosts with that arithmetic mean and log10 standard deviation `sdlog10`. The
  spread is an assumption of ours, not of any of the published models.
- **Wilber et al. (2016) IPM** ([`BdIPMSolution`](@ref)): the state at `step` (an index into
  `sol.t`); each mesh bin of ln(ZE per swab) is a uniform component (converted to log10),
  prevalence = infected/(infected + uninfected). The published x is the swab load, so with
  `σ = 0`, `lod = 0` this only re-expresses the published state.
- **Wilber et al. (2017) hybrid IPM** ([`HybridIPMSolution`](@ref)): adults only (tadpoles
  carry no load distribution in that model).
- **Louca et al. (2014)**: `out` is the result of `solve(::LoucaDOIProblem)`. The final
  distribution of `stage` (`:adult` or `:tadpole`), or with `step` a column of `PA_t`/`PT_t`
  (needs `save_distributions = true`). Bin 1 (DOI 0) is uninfected; bin `k` is a point at
  DOI `(k − 1)w`, read as ZE per swab (ESM §5.1). Point components need `σ > 0` for
  [`swab_logpdf`](@ref) of a detected swab.
"""
function swab_distribution(loads::AbstractVector{ZEPerSwab}; σ::Real = 0.0, lod::Real = 0.0)
    isempty(loads) && throw(ArgumentError("no hosts"))
    inf = [log10(x) for x in loads if load_value(x) > 0]
    prev = length(inf) / length(loads)
    return SwabLoadDistribution(prev, inf, inf, ones(length(inf)); σ, lod)
end

function swab_distribution(loads::AbstractVector{SporangiaPerHost};
        ze_per_sporangium::Union{Nothing,Real} = nothing, σ::Real = 0.0, lod::Real = 0.0)
    ze_per_sporangium === nothing && throw(ArgumentError("sporangia per host have no published " *
        "conversion to ZE per swab; pass ze_per_sporangium explicitly to state the assumption"))
    ze_per_sporangium > 0 || throw(ArgumentError("ze_per_sporangium must be positive"))
    @warn _SPORANGIA_WARNING ze_per_sporangium
    ze = [ZEPerSwab(; log10 = log10(x) + log10(ze_per_sporangium)) for x in loads]
    return swab_distribution(ze; σ, lod)
end

function swab_distribution(mean_load::ZEPerSwab; sdlog10::Real, prevalence::Real = 1.0,
        lod::Real = 0.0)
    sdlog10 > 0 || throw(ArgumentError("sdlog10 must be positive"))
    load_value(mean_load) > 0 || throw(ArgumentError("the mean load must be positive"))
    # arithmetic mean of a lognormal: 10^(μ10 + sdlog10² ln10 / 2)
    μ10 = log10(mean_load) - sdlog10^2 * log(10) / 2
    return SwabLoadDistribution(prevalence, [μ10], [μ10], [1.0]; σ = sdlog10, lod)
end

function _ipm_swabs(S::Real, n::AbstractVector{<:Real}, mesh::IPMMesh, σ, lod)
    e = meshedges(mesh) ./ log(10)
    I = sum(n)
    tot = S + I
    tot > 0 || throw(ArgumentError("no hosts in the state"))
    return SwabLoadDistribution(I / tot, e[1:end-1], e[2:end], max.(n, 0.0); σ, lod)
end

function swab_distribution(sol::BdIPMSolution; step::Integer = length(sol.t), σ::Real = 0.0,
        lod::Real = 0.0)
    1 <= step <= length(sol.t) || throw(BoundsError(sol.t, step))
    return _ipm_swabs(sol.S[step], view(sol.I, :, step), sol.mesh, σ, lod)
end

function swab_distribution(sol::HybridIPMSolution; step::Integer = length(sol.t),
        σ::Real = 0.0, lod::Real = 0.0)
    1 <= step <= length(sol.t) || throw(BoundsError(sol.t, step))
    n = sol.mesh.n
    return _ipm_swabs(sol.u[4, step], view(sol.u, 5:4+n, step), sol.mesh, σ, lod)
end

function swab_distribution(p::LoucaParameters, out::NamedTuple; stage::Symbol = :adult,
        step::Union{Nothing,Integer} = nothing, σ::Real = 0.0, lod::Real = 0.0)
    stage in (:adult, :tadpole) || throw(ArgumentError("stage must be :adult or :tadpole"))
    doi = stage === :adult ? out.doi.adult : out.doi.tadpole
    P = if step === nothing
        stage === :adult ? out.PA : out.PT
    else
        haskey(out, :PA_t) || throw(ArgumentError("solve with save_distributions = true to pick a step"))
        stage === :adult ? out.PA_t[:, step] : out.PT_t[:, step]
    end
    length(P) == length(doi) == p.m || throw(DimensionMismatch("distribution and DOI bins differ"))
    tot = sum(max.(P, 0.0))
    tot > 0 || throw(ArgumentError("no hosts of stage $stage in the state"))
    inf = 2:length(doi)
    x = log10.(doi[inf])
    return SwabLoadDistribution(sum(max.(P[inf], 0.0)) / tot, x, x, max.(P[inf], 0.0); σ, lod)
end
