# Refit of the vital-rate regressions of Wilber et al. (2016) Methods Ecol Evol 7:1182
# (doi:10.1111/2041-210X.12561) to the deposited data (Dryad doi:10.5061/dryad.06sd9;
# data/wilber2016). This reproduces a published analysis; it is not a new model.
# Wiki: lab/wiki/models/bd-wilber2016-ipm.md
#
# Data preparation follows ipm_vital_rate_parameter_analysis.R exactly:
# - loads are natural-log transformed except zeros (a load of 1 ZE would become 0 and be
#   treated as uninfected, as in R; there is none);
# - only 3-day intervals are used, and 26 °C is excluded;
# - survival s(x): infected frogs (x ≠ 0) that did not lose infection, at 20 °C only (no
#   deaths at 4 or 12 °C); logistic regression on x;
# - growth G: the same frogs alive at t + 1; generalized least squares, mean b01 + b21 T +
#   b11 x, variance σ² exp(2 c x) (nlme::gls, varExp(form = ~size), REML);
# - infection φ: uninfected survivors; logistic regression of becoming infected on T;
# - initial load G0: survivors that became infected, dropping ln ZE > 5 ("likely PCR
#   error"); gls, mean b03 + b13 T, variance σ² exp(2 c T) (varExp(form = ~temp), REML);
# - loss l: infected survivors; loss = (load_next == 0), dropping losses from ln ZE > 5;
#   logistic regression on x and T.
# Logistic regressions are fitted by Newton–Raphson (IRLS), as glm does; gls by profiling the
# (restricted) likelihood over the variance exponent, as nlme does. The tests of the Results
# section are likelihood-ratio tests between nested ML fits, as in the R script.

const WILBER_COMPONENTS = (:survival, :growth, :loss, :initial_load, :infection)

# Logistic regression by Newton–Raphson. Returns coefficients, covariance and log-likelihood.
function _logistic_fit(X::AbstractMatrix, y::AbstractVector; iters::Integer = 100)
    β = zeros(size(X, 2))
    for _ in 1:iters
        η = X * β
        μ = _logistic.(η)
        w = μ .* (1 .- μ)
        H = X' * (w .* X)
        Δ = H \ (X' * (y .- μ))
        β .+= Δ
        maximum(abs, Δ) < 1e-12 && break
    end
    μ = _logistic.(X * β)
    w = μ .* (1 .- μ)
    V = inv(Symmetric(X' * (w .* X)))
    return (coef = β, vcov = Matrix(V), loglik = _logistic_loglik(X, y, β))
end

function _logistic_loglik(X, y, β)
    η = X * β
    # y log μ + (1 − y) log(1 − μ), written stably
    return sum(@. y * η - log1p(exp(-abs(η))) - max(η, 0))
end

# Weighted least squares at variance exponent c: Var(ε_i) = σ² exp(2 c z_i).
function _gls_at(X, y, z, c, reml::Bool)
    n, k = size(X)
    w = exp.(-2c .* z)
    XtW = X' .* w'
    H = Symmetric(XtW * X)
    β = H \ (XtW * y)
    r = y .- X * β
    rss = sum(w .* r .^ 2)
    df = reml ? n - k : n
    σ2 = rss / df
    # profiled (restricted) log-likelihood, constants included as in nlme
    ll = -0.5 * (df * (log(2π) + 1 + log(σ2)) - sum(log.(w)))
    reml && (ll -= 0.5 * logdet(H))
    return (coef = β, σ2 = σ2, loglik = ll, H = H)
end

# GLS with an exponential variance function, maximising the (restricted) likelihood over c.
function _gls_fit(X::AbstractMatrix, y::AbstractVector, z::AbstractVector; reml::Bool = true,
        varexp::Bool = true)
    varexp || return merge(_gls_fit_c(X, y, z, 0.0, reml), (c = 0.0,))
    zmax = maximum(abs, z)
    a, b = -5 / zmax, 5 / zmax
    g = (sqrt(5) - 1) / 2
    f(c) = -_gls_at(X, y, z, c, reml).loglik
    x1, x2 = b - g * (b - a), a + g * (b - a)
    f1, f2 = f(x1), f(x2)
    while b - a > 1e-11
        if f1 < f2
            b, x2, f2 = x2, x1, f1
            x1 = b - g * (b - a); f1 = f(x1)
        else
            a, x1, f1 = x1, x2, f2
            x2 = a + g * (b - a); f2 = f(x2)
        end
    end
    c = (a + b) / 2
    return merge(_gls_fit_c(X, y, z, c, reml), (c = c,))
end

function _gls_fit_c(X, y, z, c, reml)
    r = _gls_at(X, y, z, c, reml)
    return (coef = r.coef, σ2 = r.σ2, vcov = r.σ2 .* Matrix(inv(r.H)), loglik = r.loglik)
end

# Regression data sets of the R script, from the rows of load_wilber2016().
function _wilber_datasets(rows; survival_temperatures = (20,))
    lt(v) = (isnan(v) || v == 0) ? v : log(v)
    d = [(x = lt(r.load), xn = lt(r.load_next), T = Float64(r.temperature),
          surv = r.survived) for r in rows if r.interval == 3 && r.temperature != 26]
    nozeros = filter(r -> r.x != 0 && !(r.xn == 0), d)
    surv = filter(r -> r.T in survival_temperatures, nozeros)
    grow = filter(r -> !isnan(r.xn), nozeros)
    alive = filter(r -> r.surv, d)
    inf = filter(r -> r.x == 0, alive)
    init = filter(r -> r.x == 0 && r.xn != 0 && r.xn <= 5, alive)
    loss = filter(r -> r.x != 0 && !(r.x > 5 && r.xn == 0), alive)
    col(v) = collect(Float64, v)
    return (
        survival = (X = hcat(ones(length(surv)), col(r.x for r in surv)),
                    y = col(r.surv for r in surv), names = (:b00, :b10)),
        growth = (X = hcat(ones(length(grow)), col(r.T for r in grow), col(r.x for r in grow)),
                  y = col(r.xn for r in grow), z = col(r.x for r in grow),
                  names = (:b01, :b21, :b11)),
        loss = (X = hcat(ones(length(loss)), col(r.x for r in loss), col(r.T for r in loss)),
                y = col(r.xn == 0 for r in loss), names = (:b02, :b12, :b22)),
        initial_load = (X = hcat(ones(length(init)), col(r.T for r in init)),
                        y = col(r.xn for r in init), z = col(r.T for r in init),
                        names = (:b03, :b13)),
        infection = (X = hcat(ones(length(inf)), col(r.T for r in inf)),
                     y = col(r.xn != 0 for r in inf), names = (:b04, :b14)))
end

"""
    WilberFit

The vital-rate regressions of Wilber et al. (2016) refitted by [`fit_wilber2016`](@ref).
Fields `survival`, `growth`, `loss`, `initial_load` and `infection` are named tuples with
the coefficients `coef` (keyed by the Table 1 symbols), their covariance `vcov` (in the
order of `coef`), the number of observations `n` and the log-likelihood `loglik` (restricted
for the two gls fits); `growth` and `initial_load` also hold the variance scale `σ2`
(`v01`, `v03`) and exponent `c` (`c01`, `c03`). `tests` holds the statistics reported in the
Results section (likelihood-ratio χ² with 1 d.f., and the t value of temperature in `G0`).
"""
struct WilberFit
    survival::NamedTuple
    growth::NamedTuple
    loss::NamedTuple
    initial_load::NamedTuple
    infection::NamedTuple
    tests::NamedTuple
end

_keyed(names, v) = NamedTuple{names}(Tuple(v))

"""
    fit_wilber2016(rows = load_wilber2016(); survival_temperatures = (20,)) -> WilberFit

Refit the five vital-rate regressions of Wilber et al. (2016, Table 1) to the deposited data,
following the authors' R script (see the comments in `src/fitting/wilber2016.jl`): logistic
regressions for survival, loss of infection and infection, and generalized least squares
with an exponential variance function (REML) for growth and initial load.

`survival_temperatures = (4, 12, 20)` fits survival to all temperatures (Appendix S2's
comparison) instead of 20 °C only.

The `tests` are those of the Results section: survival on load (`survival_load`), growth on
load, temperature and non-constant variance (`growth_load`, `growth_temperature`,
`growth_variance`), loss on temperature and load (`loss_temperature`, `loss_load`), the t
value of temperature in the initial-load mean (`initial_temperature_t`) and the variance LRT
(`initial_variance`), and infection on temperature (`infection_temperature`).
"""
function fit_wilber2016(rows = load_wilber2016(); survival_temperatures = (20,))
    D = _wilber_datasets(rows; survival_temperatures)
    lr(a, b) = 2 * (a.loglik - b.loglik)

    s = _logistic_fit(D.survival.X, D.survival.y)
    s0 = _logistic_fit(D.survival.X[:, 1:1], D.survival.y)

    g = D.growth
    gR = _gls_fit(g.X, g.y, g.z; reml = true)
    gML = _gls_fit(g.X, g.y, g.z; reml = false)
    g_noT = _gls_fit(g.X[:, [1, 3]], g.y, g.z; reml = false)
    g_nox = _gls_fit(g.X[:, [1, 2]], g.y, g.z; reml = false)
    g_cv = _gls_fit(g.X, g.y, g.z; reml = false, varexp = false)

    l = _logistic_fit(D.loss.X, D.loss.y)
    l_noT = _logistic_fit(D.loss.X[:, 1:2], D.loss.y)
    l_nox = _logistic_fit(D.loss.X[:, [1, 3]], D.loss.y)

    c = D.initial_load
    cR = _gls_fit(c.X, c.y, c.z; reml = true)
    cML = _gls_fit(c.X, c.y, c.z; reml = false)
    c_cv = _gls_fit(c.X, c.y, c.z; reml = false, varexp = false)

    f = _logistic_fit(D.infection.X, D.infection.y)
    f0 = _logistic_fit(D.infection.X[:, 1:1], D.infection.y)

    tests = (survival_load = lr(s, s0), growth_load = lr(gML, g_nox),
        growth_temperature = lr(gML, g_noT), growth_variance = lr(gML, g_cv),
        loss_temperature = lr(l, l_noT), loss_load = lr(l, l_nox),
        initial_temperature_t = cR.coef[2] / sqrt(cR.vcov[2, 2]),
        initial_temperature_df = length(c.y) - 2, initial_variance = lr(cML, c_cv),
        infection_temperature = lr(f, f0))
    pack(D, r) = (coef = _keyed(D.names, r.coef), vcov = r.vcov, n = length(D.y),
                  loglik = r.loglik)
    return WilberFit(pack(D.survival, s), merge(pack(g, gR), (σ2 = gR.σ2, c = gR.c)),
        pack(D.loss, l), merge(pack(c, cR), (σ2 = cR.σ2, c = cR.c)), pack(D.infection, f),
        tests)
end

function WilberParameters(fit::WilberFit; s0::Real = 1.0, Δt::Real = 3.0)
    return WilberParameters(; fit.survival.coef..., fit.growth.coef..., v01 = fit.growth.σ2,
        c01 = fit.growth.c, fit.loss.coef..., fit.initial_load.coef...,
        v03 = fit.initial_load.σ2, c03 = fit.initial_load.c, fit.infection.coef..., s0, Δt)
end

"""
    WilberRegression(component; rows = load_wilber2016(), survival_temperatures = (20,))

The log-likelihood of one vital-rate regression of Wilber et al. (2016) on the deposited
data, as a LogDensityProblems model (for Bayesian refits, e.g. with Turing through
`LogDensityProblems`). `component` is `:survival`, `:loss` or `:infection` (logistic
regressions; parameters the coefficients) or `:growth` or `:initial_load` (normal regressions
with variance `σ² exp(2 c z)`; parameters the coefficients, then `log σ²` and `c`).
`parameter_names(m)` lists them. The log density is the full (not restricted) likelihood; the
maximum-likelihood estimates of [`fit_wilber2016`](@ref)'s tests maximise it.
"""
struct WilberRegression
    component::Symbol
    X::Matrix{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    names::Tuple
end

function WilberRegression(component::Symbol; rows = load_wilber2016(),
        survival_temperatures = (20,))
    component in WILBER_COMPONENTS ||
        throw(ArgumentError("component must be one of $WILBER_COMPONENTS"))
    D = getfield(_wilber_datasets(rows; survival_temperatures), component)
    if haskey(D, :z)
        v = component === :growth ? (:log_v01, :c01) : (:log_v03, :c03)
        return WilberRegression(component, D.X, D.y, D.z, (D.names..., v...))
    end
    return WilberRegression(component, D.X, D.y, Float64[], D.names)
end

parameter_names(m::WilberRegression) = m.names
LogDensityProblems.dimension(m::WilberRegression) = length(m.names)
LogDensityProblems.capabilities(::Type{WilberRegression}) = LogDensityProblems.LogDensityOrder{0}()

function LogDensityProblems.logdensity(m::WilberRegression, θ::AbstractVector)
    length(θ) == length(m.names) ||
        throw(DimensionMismatch("expected $(length(m.names)) parameters"))
    k = size(m.X, 2)
    isempty(m.z) && return _logistic_loglik(m.X, m.y, θ)
    β = θ[1:k]
    logσ2, c = θ[k+1], θ[k+2]
    μ = m.X * β
    lv = logσ2 .+ 2c .* m.z
    return -0.5 * sum(@. log(2π) + lv + (m.y - μ)^2 / exp(lv))
end

"""
    growth_rate_bootstrap(fit::WilberFit; temperatures = 4:20, n = 1000,
                          mesh = IPMMesh(:wilber2016; n = 50), rng = Random.default_rng(),
                          as_coded = false) -> Matrix{Float64}

The parametric bootstrap of the population growth rate of Wilber et al. (2016, Fig. 4b):
`n` draws of the regression coefficients of each vital-rate function from a multivariate
normal with the fitted covariance (variance parameters `v`, `c` and `s0` held fixed), and
`λ` of the IPM at each temperature for each draw. Returns an `n × length(temperatures)`
matrix. Each draw is shared across temperatures (the authors redraw for every temperature;
the marginal distribution at each temperature is the same).

`as_coded = true` reproduces a slip in the authors' `monte_carlo_params`: the survival slope
draw is stored under a name that is never read, so only the survival intercept varies (with
its marginal variance). The authors' Fig. 4b was made this way; the wide interquartile
ranges at high temperature come from it.
"""
function growth_rate_bootstrap(fit::WilberFit; temperatures = 4:20,
        n::Union{Nothing,Integer} = nothing, mesh::Union{Nothing,IPMMesh} = nothing,
        rng::AbstractRNG = Random.default_rng(), as_coded::Bool = false)
    q = parameter_set(:wilber2016)
    n = n === nothing ? round(Int, q.n_bootstrap) : n
    mesh = mesh === nothing ? IPMMesh(:wilber2016; n = round(Int, q.mesh_n_bootstrap)) : mesh
    p0 = WilberParameters(fit)
    comps = (fit.survival, fit.growth, fit.loss, fit.initial_load, fit.infection)
    dists = [MvNormal(collect(values(c.coef)), Symmetric(c.vcov)) for c in comps]
    out = Matrix{Float64}(undef, n, length(temperatures))
    for i in 1:n
        kw = Pair{Symbol,Float64}[]
        for (c, d) in zip(comps, dists)
            draw = rand(rng, d)
            for (s, v) in zip(keys(c.coef), draw)
                (as_coded && s === :b10) && continue
                push!(kw, s => v)
            end
        end
        p = _with(p0; kw...)
        for (j, T) in enumerate(temperatures)
            out[i, j] = population_growth_rate(p; temperature = T, mesh)
        end
    end
    return out
end
