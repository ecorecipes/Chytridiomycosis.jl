# Refits for Wilber et al. (2017) Ecology Letters 20:1169 (doi:10.1111/ele.12814): the
# within-season vital rates (Table 1) from the Wilber 2016 temperature experiment
# (data/wilber2016), and the transmission functions with a constant zoospore pool (Table 2,
# rows 1–3) from the mesocosm experiment (Dryad doi:10.5061/dryad.p614h, CC0;
# data/wilber2017). This reproduces published analyses; it is not a new model.
# Wiki: lab/wiki/models/bd-wilber2017-ipm.md
#
# The authors fitted everything in Stan and report posterior means. Their priors
# (stat_analysis/stan_files/*.stan):
# - survival: logit s = b0 + b1 z, z the z-scored ln load at 20 °C; b ~ Normal(0, 5);
# - growth: Normal(b0 + b1 x + b2 T, sd = sqrt(σ exp(2 δ x))); b, δ ~ Normal(0, 5),
#   σ (a variance) ~ half-Cauchy(0, 1);
# - initial load: Normal(b0 + b1 T, sd = sqrt(σ exp(2 δ T))); same priors;
# - loss: logistic on x and T, flat priors (the Normal(0, 5) lines are commented out);
# - transmission: Bernoulli(1 − exp(−Λ)), Λ = b0 Δt (+ b1 I Δt or b1 (I/A) Δt);
#   b ~ half-Cauchy(0, 1).
# Data preparation follows ipm_vital_rates_bayesian.R: the vital-rate subsets are those of
# Wilber 2016 (see src/fitting/wilber2016.jl); the transmission data keep Day < 35, swab
# intervals of 4–8 days, complete rows (R's complete.cases: no missing numeric field) and
# frogs uninfected at the first swab: 333 transitions, 124 of them infections.
#
# Posterior means are computed here without MCMC: by tensor-grid quadrature (two passes, the
# second on the mode ± 8 posterior sd) for up to three parameters, and by importance sampling
# from a multivariate t centred on the mode for more. The dynamic-pool models (Table 2, rows
# 4–6) need a latent pool fitted to tank-level series from the raw file DensExpt.csv, which is
# only in the authors' GitHub repository, without a licence; they are not refitted.

"""
    wilber2017_path()

Path of the bundled data of the Wilber et al. (2017) mesocosm transmission experiment
(`data/wilber2017/formatted_dd_data.csv`; Dryad doi:10.5061/dryad.p614h, CC0). See
`data/wilber2017/README.md`.
"""
wilber2017_path() = joinpath(pkgdir(@__MODULE__), "data", "wilber2017", "formatted_dd_data.csv")

"""
    load_wilber2017() -> Vector{NamedTuple}

The 1 378 swab-to-swab transitions of adult frogs in the 16 mesocosms of Wilber et al.
(2017) ([`wilber2017_path`](@ref)), one named tuple per transition:

- `load`, `load_next`: qPCR loads (ZE per swab, natural scale; 0 uninfected; `NaN` missing,
  e.g. the frog died);
- `interval`: days between the swabs; `day`: day of the second swab;
- `swab_id`, `individual`, `tank`, `density` (adults initially in the tank: 1, 4, 8, 16);
- `weight` (g), `svl` (mm, kept as text: one entry is `"31.9?"`), `sex` (`""` if missing);
- `ze_tadpoles`, `ze_frogs`: total ZE on the tank's tadpoles and frogs at the first swab;
  `num_infected`: infected frogs in the tank then; `cum_ze_frogs`, `cum_ze_tadpoles`:
  cumulative totals.

Missing numbers are `NaN`.
"""
function load_wilber2017()
    lines = filter(!isempty ∘ strip, readlines(wilber2017_path()))
    header = split(lines[1], ',')
    header == ["size", "sizeNext", "swab_id", "time", "Day", "individual", "density",
        "weight", "svl", "sex", "ze_tadpoles", "tank_num", "ze_frogs", "num_infected",
        "cum_ze_frogs", "cum_ze_tadpoles"] ||
        throw(ArgumentError("unexpected columns in $(wilber2017_path())"))
    num(s) = isempty(strip(s)) ? NaN : parse(Float64, s)
    return map(lines[2:end]) do l
        v = split(l, ',')
        (load = num(v[1]), load_next = num(v[2]), swab_id = String(v[3]),
         interval = num(v[4]), day = num(v[5]), individual = String(v[6]),
         density = num(v[7]), weight = num(v[8]), svl = String(strip(v[9])),
         sex = String(strip(v[10])), ze_tadpoles = num(v[11]), tank = num(v[12]),
         ze_frogs = num(v[13]), num_infected = num(v[14]), cum_ze_frogs = num(v[15]),
         cum_ze_tadpoles = num(v[16]))
    end
end

# The 333 transitions of the transmission fits (ipm_vital_rates_bayesian.R). R reads `svl`
# and `sex` as text (because of "31.9?"), so complete.cases ignores them.
function _wilber2017_transmission_rows(rows)
    complete(r) = all(v -> !(v isa Float64 && isnan(v)), values(r))
    return filter(r -> r.day < 35 && 3 < r.interval < 9 && complete(r) && r.load == 0, rows)
end

const WILBER2017_TRANSMISSION_FITS = (:const, :const_dd, :const_fd)

"""
    Wilber2017Transmission(model; rows = load_wilber2017())

The log posterior of a Table 2 transmission function with a constant zoospore pool, fitted to
the 333 transitions of frogs uninfected at the first swab (a LogDensityProblems model):
infection by the next swab is Bernoulli with probability `1 − exp(−Λ)`, `Δt` the days
between swabs (4–8), and

- `:const`: `Λ = β0 Δt`;
- `:const_dd`: `Λ = (β0 + β1 I) Δt`, `I` the infected frogs in the tank at the first swab;
- `:const_fd`: `Λ = (β0 + β1 I/A) Δt`, `A` the initial number of frogs in the tank;

with half-Cauchy(0, 1) priors on `β ≥ 0`, as in the authors' Stan models.
`parameter_names(m)`, `pointwise_loglikelihood(m, θ)` and [`posterior_summary`](@ref)
(posterior means and WAIC, Table 2) work on it.
"""
struct Wilber2017Transmission
    model::Symbol
    infected::Vector{Bool}
    Δt::Vector{Float64}
    I::Vector{Float64}
    A::Vector{Float64}
end

function Wilber2017Transmission(model::Symbol; rows = load_wilber2017())
    model in WILBER2017_TRANSMISSION_FITS ||
        throw(ArgumentError("model must be one of $WILBER2017_TRANSMISSION_FITS"))
    r = _wilber2017_transmission_rows(rows)
    return Wilber2017Transmission(model, [x.load_next > 0 for x in r], [x.interval for x in r],
        [x.num_infected for x in r], [x.density for x in r])
end

parameter_names(m::Wilber2017Transmission) = m.model === :const ? (:β0,) : (:β0, :β1)
LogDensityProblems.dimension(m::Wilber2017Transmission) = length(parameter_names(m))
LogDensityProblems.capabilities(::Type{Wilber2017Transmission}) =
    LogDensityProblems.LogDensityOrder{0}()

function _transmission_force(m::Wilber2017Transmission, θ, i)
    m.model === :const && return θ[1] * m.Δt[i]
    x = m.model === :const_dd ? m.I[i] : m.I[i] / m.A[i]
    return (θ[1] + θ[2] * x) * m.Δt[i]
end

# log P(y | Λ) for y ~ Bernoulli(1 − exp(−Λ)), written stably.
_bernoulli_cloglog(y::Bool, Λ) = y ? log(-expm1(-Λ)) : -Λ

"""
    pointwise_loglikelihood(m::Wilber2017Transmission, θ) -> Vector{Float64}

The log-likelihood of each of the transitions in `m` at parameters `θ` (for WAIC).
"""
function pointwise_loglikelihood(m::Wilber2017Transmission, θ::AbstractVector)
    length(θ) == LogDensityProblems.dimension(m) ||
        throw(DimensionMismatch("expected $(LogDensityProblems.dimension(m)) parameters"))
    return [_bernoulli_cloglog(m.infected[i], _transmission_force(m, θ, i))
            for i in eachindex(m.infected)]
end

function LogDensityProblems.logdensity(m::Wilber2017Transmission, θ::AbstractVector)
    any(<(0), θ) && return -Inf
    lp = sum(t -> -log1p(t^2), θ)              # half-Cauchy(0, 1), up to a constant
    ll = 0.0
    for i in eachindex(m.infected)
        Λ = _transmission_force(m, θ, i)
        (m.infected[i] && Λ == 0) && return -Inf
        ll += _bernoulli_cloglog(m.infected[i], Λ)
    end
    return lp + ll
end

const WILBER2017_COMPONENTS = (:survival, :growth, :loss, :initial_load)

"""
    Wilber2017Regression(component; rows = load_wilber2016())

The log posterior of one within-season vital-rate function of Wilber et al. (2017, Table 1),
refitted to the Wilber 2016 temperature experiment with the authors' Stan priors (a
LogDensityProblems model). `component` and parameters:

- `:survival`: `(b00, b10)`, logistic regression on the z-scored ln load of the 20 °C frogs,
  Normal(0, 5) priors;
- `:growth`: `(b01, b11, b21, v01, c01)`, normal with mean `b01 + b11 x + b21 T` and variance
  `v01 exp(2 c01 x)`, Normal(0, 5) priors and half-Cauchy(0, 1) on `v01`;
- `:loss`: `(b02, b12, b22)`, logistic regression on `x` and `T`, flat priors;
- `:initial_load`: `(b03, b13, v03, c03)`, normal with mean `b03 + b13 T` and variance
  `v03 exp(2 c03 T)`, priors as for growth.

The data subsets are those of [`fit_wilber2016`](@ref). `zscore` holds the mean and standard
deviation used to z-score the survival loads. [`posterior_summary`](@ref) gives posterior
means.
"""
struct Wilber2017Regression
    component::Symbol
    X::Matrix{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    names::Tuple
    zscore::NamedTuple
end

function Wilber2017Regression(component::Symbol; rows = load_wilber2016())
    component in WILBER2017_COMPONENTS ||
        throw(ArgumentError("component must be one of $WILBER2017_COMPONENTS"))
    D = _wilber_datasets(rows)
    if component === :survival
        x = D.survival.X[:, 2]
        μ = sum(x) / length(x)
        s = sqrt(sum(abs2, x .- μ) / (length(x) - 1))
        X = hcat(ones(length(x)), (x .- μ) ./ s)
        return Wilber2017Regression(component, X, D.survival.y, Float64[], (:b00, :b10),
            (mean = μ, sd = s))
    elseif component === :loss
        return Wilber2017Regression(component, D.loss.X, D.loss.y, Float64[],
            (:b02, :b12, :b22), (mean = NaN, sd = NaN))
    elseif component === :growth
        g = D.growth                       # columns 1, T, x
        X = g.X[:, [1, 3, 2]]
        return Wilber2017Regression(component, X, g.y, g.z, (:b01, :b11, :b21, :v01, :c01),
            (mean = NaN, sd = NaN))
    end
    c = D.initial_load
    return Wilber2017Regression(component, c.X, c.y, c.z, (:b03, :b13, :v03, :c03),
        (mean = NaN, sd = NaN))
end

parameter_names(m::Wilber2017Regression) = m.names
LogDensityProblems.dimension(m::Wilber2017Regression) = length(m.names)
LogDensityProblems.capabilities(::Type{Wilber2017Regression}) =
    LogDensityProblems.LogDensityOrder{0}()

function LogDensityProblems.logdensity(m::Wilber2017Regression, θ::AbstractVector)
    length(θ) == length(m.names) ||
        throw(DimensionMismatch("expected $(length(m.names)) parameters"))
    normal5(b) = -b^2 / 50
    if m.component === :survival
        return _logistic_loglik(m.X, m.y, θ) + normal5(θ[1]) + normal5(θ[2])
    elseif m.component === :loss
        return _logistic_loglik(m.X, m.y, θ)
    end
    k = size(m.X, 2)
    β = θ[1:k]
    v, c = θ[k+1], θ[k+2]
    v > 0 || return -Inf
    lp = sum(normal5, β) + normal5(c) - log1p(v^2)
    μ = m.X * β
    lv = log(v) .+ 2c .* m.z
    return lp - 0.5 * sum(@. log(2π) + lv + (m.y - μ)^2 / exp(lv))
end

# Damped Newton on the log density with finite-difference derivatives; returns the mode and
# the inverse negative Hessian there.
function _fd_grad_hess(f, x)
    d = length(x)
    h = [1e-4 * max(abs(v), 1e-2) for v in x]
    fx = f(x)
    g = zeros(d)
    H = zeros(d, d)
    e(i) = (v = zeros(d); v[i] = h[i]; v)
    for i in 1:d
        fp, fm = f(x .+ e(i)), f(x .- e(i))
        g[i] = (fp - fm) / (2h[i])
        H[i, i] = (fp - 2 * fx + fm) / h[i]^2
        for j in 1:i-1
            H[i, j] = H[j, i] = (f(x .+ e(i) .+ e(j)) - f(x .+ e(i) .- e(j)) -
                                 f(x .- e(i) .+ e(j)) + f(x .- e(i) .- e(j))) / (4h[i] * h[j])
        end
    end
    return fx, g, H
end

function _mode(f, x0; iters = 100)
    x = collect(Float64, x0)
    for _ in 1:iters
        fx, g, H = _fd_grad_hess(f, x)
        step = -(Symmetric(H) \ g)
        all(isfinite, step) || break
        t = 1.0
        while t > 1e-8 && !(f(x .+ t .* step) >= fx)
            t /= 2
        end
        x .+= t .* step
        maximum(abs, t .* step) < 1e-10 && break
    end
    _, _, H = _fd_grad_hess(f, x)
    return x, Matrix(inv(Symmetric(-H)))
end

function _start(m::Wilber2017Transmission)
    return m.model === :const ? [0.05] : [0.03, m.model === :const_dd ? 0.03 : 0.3]
end

function _start(m::Wilber2017Regression)
    if m.component in (:survival, :loss)
        return _logistic_fit(m.X, m.y).coef
    end
    r = _gls_fit(m.X, m.y, m.z; reml = false)
    return vcat(r.coef, r.σ2, r.c)
end

"""
    posterior_summary(m; n = nothing, width = 8, draws = 20_000, rng = Random.Xoshiro(2017))
        -> NamedTuple

Posterior summaries of a [`Wilber2017Transmission`](@ref) or [`Wilber2017Regression`](@ref)
model without MCMC. The posterior is located by Newton's method on the log density; for up to
three parameters its moments come from tensor-grid quadrature with `n` points per dimension
(default 400, 150 and 50 for one, two and three parameters) on the mode ± `width` sd
(Laplace), and for more from importance sampling with `draws` draws from a multivariate t
(5 d.f.) centred on the mode. Returns `names`, `mean`, `sd`, `mode` and, for a transmission
model, `waic` (`−2 Σᵢ (log E[pᵢ] − Var[log pᵢ])`, deviance scale, as the authors' code
computes it) and `p_waic`; `ess` is the effective sample size of the quadrature or
importance weights. [`posterior_draws`](@ref) resamples the same nodes.
"""
function posterior_summary(m::Union{Wilber2017Transmission,Wilber2017Regression};
        n::Union{Nothing,Integer} = nothing, width::Real = 8, draws::Integer = 20_000,
        rng::AbstractRNG = Random.Xoshiro(2017))
    nodes, w, mode = _posterior_nodes(m; n, width, draws, rng)
    μ = nodes' * w
    v = ((nodes .- μ') .^ 2)' * w
    out = (names = parameter_names(m), mean = μ, sd = sqrt.(v), mode = mode,
        ess = 1 / sum(abs2, w))
    m isa Wilber2017Transmission || return out
    keep = findall(>(1e-14), w)
    wk = w[keep] ./ sum(w[keep])
    nobs = length(m.infected)
    E = zeros(nobs)                    # posterior E[p_i]
    M1 = zeros(nobs)                   # E[log p_i]
    M2 = zeros(nobs)                   # E[(log p_i)^2]
    θ = zeros(size(nodes, 2))
    for (k, g) in enumerate(keep)
        θ .= view(nodes, g, :)
        for i in 1:nobs
            l = _bernoulli_cloglog(m.infected[i], _transmission_force(m, θ, i))
            E[i] += wk[k] * exp(l)
            M1[i] += wk[k] * l
            M2[i] += wk[k] * l^2
        end
    end
    lpd = log.(E)
    pw = M2 .- M1 .^ 2
    return merge(out, (waic = -2 * sum(lpd .- pw), p_waic = sum(pw)))
end

# Quadrature or importance nodes and normalised weights of the posterior, and its mode.
function _posterior_nodes(m; n = nothing, width = 8, draws = 20_000,
        rng = Random.Xoshiro(2017))
    f(θ) = LogDensityProblems.logdensity(m, θ)
    mode, Σ = _mode(f, _start(m))
    d = length(mode)
    sd = sqrt.(max.(diag(Σ), 0.0))
    if d <= 3
        nodes, logw = _grid_nodes(f, m, mode, sd, n === nothing ? (400, 150, 50)[d] : n, width)
    else
        nodes, logw = _importance_nodes(f, mode, Σ, draws, rng)
    end
    w = exp.(logw .- maximum(logw))
    w ./= sum(w)
    return nodes, w, mode
end

"""
    posterior_draws(m, ndraws; rng = Random.default_rng(), kwargs...) -> Matrix{Float64}

`ndraws` approximate posterior draws (`ndraws × dimension`) of a
[`Wilber2017Transmission`](@ref) or [`Wilber2017Regression`](@ref) model, resampled with
replacement from the weighted quadrature or importance nodes of
[`posterior_summary`](@ref) (`kwargs` are passed to it; its `rng` is the one for the
importance proposal). Quadrature draws are jittered uniformly within their grid cell.
"""
function posterior_draws(m::Union{Wilber2017Transmission,Wilber2017Regression},
        ndraws::Integer; rng::AbstractRNG = Random.default_rng(), kwargs...)
    nodes, w, _ = _posterior_nodes(m; kwargs...)
    c = cumsum(w)
    d = size(nodes, 2)
    step = d <= 3 ? [_cell_width(view(nodes, :, j)) for j in 1:d] : zeros(d)
    out = Matrix{Float64}(undef, ndraws, d)
    for s in 1:ndraws
        g = min(searchsortedfirst(c, rand(rng) * c[end]), length(c))
        for j in 1:d
            out[s, j] = nodes[g, j] + step[j] * (rand(rng) - 0.5)
        end
    end
    return out
end

function _cell_width(x)
    u = sort!(unique(x))
    return length(u) > 1 ? u[2] - u[1] : 0.0
end

function _grid_nodes(f, m, mode, sd, n, width)
    d = length(mode)
    lo = mode .- width .* sd
    hi = mode .+ width .* sd
    m isa Wilber2017Transmission && (lo = max.(lo, 0.0))
    axes_ = [range(lo[i], hi[i]; length = n + 1) for i in 1:d]
    mids = [(a[1:end-1] .+ a[2:end]) ./ 2 for a in axes_]
    pts = vec(collect(Iterators.product(mids...)))
    nodes = [p[i] for p in pts, i in 1:d]
    logw = [f(nodes[g, :]) for g in axes(nodes, 1)]
    return nodes, logw
end

function _importance_nodes(f, mode, Σ, draws, rng)
    d = length(mode)
    ν = 5
    L = LinearAlgebra.cholesky(Symmetric(1.5^2 .* Σ)).L
    nodes = Matrix{Float64}(undef, draws, d)
    logw = Vector{Float64}(undef, draws)
    Linv = inv(L)
    for s in 1:draws
        zt = randn(rng, d)
        u = sum(abs2, randn(rng, ν))
        x = mode .+ L * (zt .* sqrt(ν / u))
        nodes[s, :] .= x
        q = sum(abs2, Linv * (x .- mode))
        logq = -(ν + d) / 2 * log1p(q / ν)          # t density up to a constant
        logw[s] = f(x) - logq
    end
    logw[.!isfinite.(logw)] .= -Inf
    return nodes, logw
end

# --- Fig. 2: R0 as the median over posterior draws (validation targets wilber2017_fig2*) -------
#
# R0_environmental.R plots, at each temperature and density, the median of R0 over 500 draws:
# one random posterior row per vital-rate function and per transmission function. Here the
# vital-rate rows come from our refitted posteriors (posterior_draws), and β0, β1 of the
# density-dependent dynamic-pool model from lognormals with the Table 2 means and the posterior
# CVs of the preprint appendix (cv_β0, cv_β1; the latent-pool posterior is not refitted). The
# draws are cached per (ndraws, seed).
const _WILBER2017_DRAWS = Dict{Tuple{Int,Int},Vector{Wilber2017Parameters}}()

function _wilber2017_posterior_parameters(ndraws::Integer, seed::Integer)
    get!(_WILBER2017_DRAWS, (ndraws, seed)) do
        rng = Random.Xoshiro(seed)
        p = Wilber2017Parameters(:wilber2017; transmission = :zpool_dd)
        tr = parameter_set(:wilber2017; set = "zpool_dd")
        D = [posterior_draws(Wilber2017Regression(c), ndraws; rng, draws = 100_000)
             for c in WILBER2017_COMPONENTS]
        names = reduce(vcat, [collect(parameter_names(Wilber2017Regression(c)))
                              for c in WILBER2017_COMPONENTS])
        lnorm(m, cv) = (σ = sqrt(log1p(cv^2)); LogNormal(log(m) - σ^2 / 2, σ))
        β0 = rand(rng, lnorm(tr.β0, tr.cv_β0), ndraws)
        β1 = rand(rng, lnorm(tr.β1, tr.cv_β1), ndraws)
        map(1:ndraws) do k
            v = reduce(vcat, [d[k, :] for d in D])
            _with(p; (s => x for (s, x) in zip(names, v))..., β0 = β0[k], β1 = β1[k])
        end
    end
end

# Median over posterior draws of ln R0 (`kind = :total` or `:host`) at temperature T and
# susceptible density S.
function _wilber2017_lnR0_median(T::Real, S::Real; kind::Symbol = :total,
        ndraws::Integer = 1000, seed::Integer = 2)
    qs = _wilber2017_posterior_parameters(ndraws, seed)
    return median(log(getfield(basic_reproduction_number(q; temperature = T, S), kind))
                  for q in qs)
end
