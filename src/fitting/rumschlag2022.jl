# Refit of Rumschlag et al. 2022 (J Anim Ecol 91:170, doi:10.1111/1365-2656.13612) to the
# authors' deposited assay data (figshare doi:10.6084/m9.figshare.16814701, CC BY 4.0;
# data/rumschlag2022). This reproduces a published analysis; it is not a new model.
# Wiki: lab/wiki/models/bd-rumschlag2022-zoospore-persistence.md
#
# Likelihood (paper eqs 7 and 10):
#   survival assay, each unit:  log(Z_obs + 0.1) ~ Normal(log(Z_pred + 0.1), SD_i), i = L, D
#   infection assay, each tadpole: infected ~ Bernoulli(p), p from eq. 10 with τ = 1 day
# with ZL(0) one fitted initial concentration shared by both assays (inferred from the
# degrees of freedom in Table 1: m, d, β, SD_L, SD_D and ZL(0) give df = 6 for the model
# with one estimate of each). Parameters with two estimates split ASW from the three ponds.
# Survival times are the recorded times of the counts (time.hrs); with the nominal times
# (0, 6, 16, 24, 48 h) Table 1 is not reproduced.

"""
    RumschlagCount(water, t, live, dead)

One unit of the Rumschlag et al. (2022) zoospore survival assay: water source `water` (one
of `:ASW`, `:TroutCreek`, `:Flatwoods`, `:GreenSwamp`), sampling time `t` (days since the
zoospores were added), and live and dead zoospores counted by Trypan blue (cells per µl).
"""
struct RumschlagCount
    water::Symbol
    t::Float64
    live::Float64
    dead::Float64
end

"""
    RumschlagTadpole(water, lag, infected)

One tadpole of the Rumschlag et al. (2022) infection assay: water source `water`, the lag
`lag` (days) between adding zoospores and adding the tadpole, and whether it was infected
(qPCR of the mouthparts two weeks after a one-day exposure).
"""
struct RumschlagTadpole
    water::Symbol
    lag::Float64
    infected::Bool
end

const RUMSCHLAG_CODES = Dict("ASW" => :ASW, "A" => :TroutCreek, "B" => :Flatwoods,
    "C" => :GreenSwamp)

"""
    rumschlag2022_path(assay::Symbol)

Path of the bundled data of Rumschlag et al. (2022) for `assay` `:survival` or `:infection`
(figshare doi:10.6084/m9.figshare.16814701, CC BY 4.0; see `data/rumschlag2022/README.md`).
"""
function rumschlag2022_path(assay::Symbol)
    assay in (:survival, :infection) ||
        throw(ArgumentError("assay must be :survival or :infection; got :$assay"))
    return joinpath(pkgdir(@__MODULE__), "data", "rumschlag2022", string(assay) * ".csv")
end

function _read_rows(path)
    lines = filter(!isempty ∘ strip, readlines(path))
    header = split(lines[1], ',')
    function col(name)
        i = findfirst(==(name), header)
        i === nothing && throw(ArgumentError("column $name missing from $path"))
        return i
    end
    return [split(l, ',') for l in lines[2:end]], col
end

"""
    load_rumschlag2022(; time = :exact) -> (survival, infection)

The two assays of Rumschlag et al. (2022): a `Vector{RumschlagCount}` (199 units) and a
`Vector{RumschlagTadpole}` (197 tadpoles). Pond codes are mapped as in the deposit's
metadata (A = Trout Creek, B = Flatwoods, C = Green Swamp). Survival times are the recorded
time of each count (`time.hrs`, 0 to 50 h) with `time = :exact`, or the nominal sampling
times (0, 6, 16, 24, 48 h) with `time = :nominal`. All times are returned in days. The
paper's Table 1 is reproduced with the recorded times, not the nominal ones.
"""
function load_rumschlag2022(; time::Symbol = :exact)
    time in (:nominal, :exact) || throw(ArgumentError("time must be :nominal or :exact"))
    rows, col = _read_rows(rumschlag2022_path(:survival))
    iw, it = col("water.trt"), col(time === :nominal ? "time.point" : "time.hrs")
    il, id = col("live.cells.ul"), col("dead.cells.ul")
    survival = [RumschlagCount(RUMSCHLAG_CODES[r[iw]], parse(Float64, r[it]) / 24,
        parse(Float64, r[il]), parse(Float64, r[id])) for r in rows]
    rows, col = _read_rows(rumschlag2022_path(:infection))
    iw, it, ii = col("water.trt"), col("time.trt"), col("infect.TM")
    infection = [RumschlagTadpole(RUMSCHLAG_CODES[r[iw]], parse(Float64, r[it]) / 24,
        r[ii] == "1") for r in rows]
    return (survival = survival, infection = infection)
end

"""
    RumschlagModel(; m = 4, d = 1, β = 1, data = load_rumschlag2022(), τ = 1.0)

Likelihood of the Rumschlag et al. (2022) assays (paper eqs 7 and 10) for one of the 15
models of their Table 1. `m`, `d` and `β` give the number of estimates of each:

- `1`: one estimate shared by all four waters;
- `2`: one for artificial spring water and one for the three ponds together;
- `4`: one per water.

In Table 1, every parameter that is not shared has the same number of estimates (2 or 4). The
default is the paper's selected model. `τ` is the exposure time (days).

The free parameters, all on the log scale, are listed by [`parameter_names`](@ref):

- the rates `m` and `d` (day⁻¹);
- `β` (µl zoospore⁻¹ day⁻¹);
- `ZL0`, the initial live zoospores (µl⁻¹), shared by both assays;
- `SD_L` and `SD_D`, the standard deviations of the log counts.

The model implements the LogDensityProblems interface. The log density is the
log-likelihood, with a flat prior on the log parameters. Fit it with
[`maximum_likelihood`](@ref).
"""
struct RumschlagModel
    survival::Vector{RumschlagCount}
    infection::Vector{RumschlagTadpole}
    k::NamedTuple{(:m, :d, :β),NTuple{3,Int}}
    τ::Float64
end
function RumschlagModel(; m::Integer = 4, d::Integer = 1, β::Integer = 1,
        data = load_rumschlag2022(), τ::Real = 1.0)
    all(in((1, 2, 4)), (m, d, β)) || throw(ArgumentError("m, d and β must each be 1, 2 or 4"))
    return RumschlagModel(collect(data.survival), collect(data.infection),
        (m = Int(m), d = Int(d), β = Int(β)), float(τ))
end

_rumschlag_labels(k) = k == 1 ? ("",) : k == 2 ? ("_ASW", "_pond") :
    Tuple("_" * string(w) for w in RUMSCHLAG_WATERS)
_rumschlag_group(k, water::Symbol) =
    k == 1 ? 1 : k == 2 ? (water === :ASW ? 1 : 2) : findfirst(==(water), RUMSCHLAG_WATERS)

parameter_names(model::RumschlagModel) = Tuple(vcat(
    [Symbol(string(p) * l) for p in (:m, :d, :β) for l in _rumschlag_labels(model.k[p])],
    [:ZL0, :SD_L, :SD_D]))

LogDensityProblems.dimension(model::RumschlagModel) = sum(model.k) + 3
LogDensityProblems.capabilities(::Type{RumschlagModel}) = LogDensityProblems.LogDensityOrder{0}()

"""
    nobs(model::RumschlagModel)

Number of observations used in AICc: two log counts per survival unit plus one per tadpole
(595 for the deposited data).
"""
nobs(model::RumschlagModel) = 2 * length(model.survival) + length(model.infection)

# natural-scale parameters from log parameters x
function _rumschlag_unpack(model::RumschlagModel, x)
    km, kd, kβ = model.k
    v = exp.(x)
    return (m = v[1:km], d = v[(km + 1):(km + kd)], β = v[(km + kd + 1):(km + kd + kβ)],
        ZL0 = v[end - 2], SD_L = v[end - 1], SD_D = v[end])
end

"""
    rumschlag_parameters(model::RumschlagModel, x, water::Symbol) -> RumschlagParameters

The [`RumschlagParameters`](@ref) of water source `water` at log parameters `x` of `model`
(for example `fit.x` from [`maximum_likelihood`](@ref)).
"""
function rumschlag_parameters(model::RumschlagModel, x::AbstractVector, water::Symbol)
    water in RUMSCHLAG_WATERS || throw(ArgumentError("unknown water :$water"))
    q = _rumschlag_unpack(model, x)
    g(p) = _rumschlag_group(model.k[p], water)
    return RumschlagParameters(; m = q.m[g(:m)], d = q.d[g(:d)], β = q.β[g(:β)])
end

function LogDensityProblems.logdensity(model::RumschlagModel, x::AbstractVector)
    length(x) == LogDensityProblems.dimension(model) ||
        throw(DimensionMismatch("expected $(LogDensityProblems.dimension(model)) parameters"))
    all(isfinite, x) || return -Inf
    q = _rumschlag_unpack(model, x)
    ps = Dict(w => rumschlag_parameters(model, x, w) for w in RUMSCHLAG_WATERS)
    ll = 0.0
    for c in model.survival
        ZL, ZD = pool_state(ZoosporePool(ps[c.water]), q.ZL0, c.t)
        ll += logpdf(Normal(log(ZL + 0.1), q.SD_L), log(c.live + 0.1))
        ll += logpdf(Normal(log(ZD + 0.1), q.SD_D), log(c.dead + 0.1))
    end
    for h in model.infection
        pr = infection_probability(ps[h.water], q.ZL0, h.lag; τ = model.τ)
        ll += h.infected ? log(pr) : log1p(-pr)
    end
    return isfinite(ll) ? ll : -Inf
end

const RUMSCHLAG_START = (m = 1.8, d = 1.1, β = 5e-3, ZL0 = 145.0, SD_L = 0.5, SD_D = 0.5)

function _rumschlag_start(model::RumschlagModel)
    s = RUMSCHLAG_START
    return log.(vcat(fill(s.m, model.k.m), fill(s.d, model.k.d), fill(s.β, model.k.β),
        s.ZL0, s.SD_L, s.SD_D))
end

"""
    maximum_likelihood(model::RumschlagModel; starts = 4, iters = 3000, spread = 0.3,
        rng = Random.Xoshiro(1))

Maximum likelihood fit of a [`RumschlagModel`](@ref), by Nelder–Mead on the log
parameters. It starts from a default point and from `starts − 1` points drawn around it
(normal, standard deviation `spread`), and polishes each start with a smaller simplex.

Returns a `NamedTuple`:

- `model` and `x` (the log parameters);
- `estimates`, the natural-scale parameters by name;
- `loglik`, `df` (number of parameters) and `aicc`, with `n = nobs(model)`;
- `logliks`, the best value from each start, to check that the starts agree.

The lifespans of the paper are `24/m` hours ([`zoospore_lifespan`](@ref)). Use
[`rumschlag_parameters`](@ref) to get each water's parameters.
"""
function maximum_likelihood(model::RumschlagModel; starts::Integer = 4, iters::Integer = 3000,
        spread::Real = 0.3, rng::AbstractRNG = Random.Xoshiro(1))
    f(x) = -LogDensityProblems.logdensity(model, x)
    x0 = _rumschlag_start(model)
    best, bestf = x0, Inf
    logliks = Float64[]
    for s in 1:starts
        x = s == 1 ? x0 : x0 .+ spread .* randn(rng, length(x0))
        x, fx = nelder_mead(f, x; iters = iters, step = 0.5)
        x, fx = nelder_mead(f, x; iters = iters ÷ 2, step = 0.05)
        push!(logliks, -fx)
        fx < bestf && ((best, bestf) = (x, fx))
    end
    k, n = length(best), nobs(model)
    aicc = 2bestf + 2k + 2k * (k + 1) / (n - k - 1)
    est = NamedTuple{parameter_names(model)}(Tuple(exp.(best)))
    return (model = model, x = best, estimates = est, loglik = -bestf, df = k, aicc = aicc,
        logliks = logliks)
end

"""
    rumschlag_model_selection(; data = load_rumschlag2022(), kwargs...) -> Vector{NamedTuple}

Fit the 15 models of Rumschlag et al. (2022) Table 1 and rank them by AICc. The models are
all the ways of giving each of `m`, `d` and `β` one shared estimate, except that the
parameters that vary all use the same split (2 or 4). Keywords are passed to
[`maximum_likelihood`](@ref).

Each entry has the fields of the fit, plus `k` (estimates of `m`, `d` and `β`), `Δaicc` and
`weight` (the Akaike weight). Entries are sorted by AICc.
"""
function rumschlag_model_selection(; data = load_rumschlag2022(), kwargs...)
    designs = NamedTuple{(:m, :d, :β),NTuple{3,Int}}[]
    for vary in Iterators.product((false, true), (false, true), (false, true))
        ks = any(vary) ? (2, 4) : (1,)
        for k in ks
            push!(designs, NamedTuple{(:m, :d, :β)}(Tuple(v ? k : 1 for v in vary)))
        end
    end
    fits = [maximum_likelihood(RumschlagModel(; m = k.m, d = k.d, β = k.β, data); kwargs...)
            for k in designs]
    sort!(fits; by = f -> f.aicc)
    a0 = first(fits).aicc
    w = [exp(-(f.aicc - a0) / 2) for f in fits]
    w ./= sum(w)
    return [merge(f, (k = f.model.k, Δaicc = f.aicc - a0, weight = w[i]))
            for (i, f) in enumerate(fits)]
end


"""
    rumschlag_covariance(fit) -> Matrix{Float64}

Asymptotic covariance matrix of the **log** parameters of a [`maximum_likelihood`](@ref)
fit of a [`RumschlagModel`](@ref), in the order of [`parameter_names`](@ref). It is the
inverse of minus a central finite-difference Hessian of the log-likelihood, with step `h` on
the log scale.
"""
function rumschlag_covariance(fit::NamedTuple; h::Real = 1e-4)
    f(x) = LogDensityProblems.logdensity(fit.model, x)
    x, n = fit.x, length(fit.x)
    H = zeros(n, n)
    for i in 1:n, j in i:n
        ei, ej = zeros(n), zeros(n)
        ei[i] = h; ej[j] = h
        H[i, j] = H[j, i] = (f(x + ei + ej) - f(x + ei - ej) - f(x - ei + ej) + f(x - ei - ej)) /
                            (4h^2)
    end
    return Matrix(inv(Symmetric(-H)))
end

# 95% Wald interval for a rate r with standard error se, carried to the duration 24/r (hours)
# and reported as the paper does: the estimate ± half the width of the transformed interval
function _paper_halfwidth(r, se)
    lo, hi = r - 1.96se, r + 1.96se
    lo > 0 || return Inf
    return (24 / lo - 24 / hi) / 2
end

"""
    rumschlag_summary(fit; ZL0_infection = 1e5) -> NamedTuple

The quantities that Rumschlag et al. (2022) report, computed from a
[`maximum_likelihood`](@ref) fit and its [`rumschlag_covariance`](@ref):

- `lifespan`: the expected lifespan `24/m` (hours) of each water, by water;
- `persistence`: the expected time `24/d` (hours) from death to decomposition, by water;
- `β`: the transmission rate on the paper's scale, `β ZL0 / ZL0_infection`.

Each value is a `(estimate, halfwidth)` pair. The paper writes its 95% intervals as
"estimate ± halfwidth". These are reproduced when the halfwidth is half the width of the Wald
interval of the rate (`m ± 1.96 se`), carried to `24/m` hours.

The paper's β does not match `β` in this package's parametrization. The infection data
identify only the product `β ZL0`. The paper's value (7.1 × 10⁻⁶) equals that product divided
by 10⁵ to two figures. So the paper appears to have used a fixed initial concentration of 10⁵
(units unstated) in eq. 10 for the infection assay. `ZL0_infection = 1e5` is that inferred
convention, not a value printed in the paper.
"""
function rumschlag_summary(fit::NamedTuple; ZL0_infection::Real = 1e5)
    model = fit.model
    names = parameter_names(model)
    V = rumschlag_covariance(fit)
    idx(s) = findfirst(==(s), names)
    function rate(p, water)   # (rate, se of the rate) for parameter p in water
        k = model.k[p]
        s = Symbol(string(p) * _rumschlag_labels(k)[_rumschlag_group(k, water)])
        r = fit.estimates[s]
        return r, r * sqrt(V[idx(s), idx(s)])
    end
    duration(p) = NamedTuple{RUMSCHLAG_WATERS}(Tuple(begin
        r, se = rate(p, w)
        (24 / r, _paper_halfwidth(r, se))
    end for w in RUMSCHLAG_WATERS))
    βs = NamedTuple{RUMSCHLAG_WATERS}(Tuple(begin
        k = model.k.β
        s = Symbol("β" * _rumschlag_labels(k)[_rumschlag_group(k, w)])
        i, j = idx(s), idx(:ZL0)
        b = fit.estimates[s] * fit.estimates.ZL0 / ZL0_infection
        (b, 1.96 * b * sqrt(V[i, i] + V[j, j] + 2V[i, j]))
    end for w in RUMSCHLAG_WATERS))
    return (lifespan = duration(:m), persistence = duration(:d), β = βs)
end

# One fit per model design per session, for the validation registry
# (data/validation/targets.toml), whose calls would otherwise refit for every target.
const _RUMSCHLAG_REFERENCE = Dict{NamedTuple,Any}()
_rumschlag_reference_fit(; m = 4, d = 1, β = 1) =
    get!(() -> maximum_likelihood(RumschlagModel(; m, d, β); starts = 2),
        _RUMSCHLAG_REFERENCE, (m = m, d = d, β = β))
