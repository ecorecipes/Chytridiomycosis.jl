# Bayesian fit of the frog load equation to the Knapp et al. 2022 swabs (new): the prior and
# posterior, pointwise log-likelihoods, PSIS-LOO and WAIC, and the Turing extension. The long
# NUTS run that the wiki and docs report is not repeated here; one very short chain checks that
# the extension samples.
using Turing: Turing
using LogDensityProblems: LogDensityProblems
using Distributions: LogNormal, Normal, logpdf
using Random: Xoshiro

@testset "Bayesian Knapp fit (new)" begin
    LDP = LogDensityProblems
    m = KnappLoadModel(:decay)
    x = knapp_start(m)

    # prior: one lognormal per parameter, on the natural scale
    pr = knapp_prior(m)
    @test keys(pr) == parameter_names(m)
    @test all(d -> d isa LogNormal, values(pr))
    @test exp(pr.ι.μ) ≈ 1e4 && pr.τ.σ == 1.5
    @test keys(knapp_prior(KnappLoadModel(:applied)))[end] === :K

    # posterior = normal prior on the log parameters + log-likelihood
    post = KnappPosterior(m)
    @test LDP.dimension(post) == length(x)
    @test LDP.capabilities(KnappPosterior) isa LDP.LogDensityOrder{0}
    @test parameter_names(post) == parameter_names(m)
    lp = sum(logpdf(Normal(d.μ, d.σ), xi) for (d, xi) in zip(values(pr), x))
    @test LDP.logdensity(post, x) ≈ lp + LDP.logdensity(m, x)
    # the normal density of log θ is the lognormal density of θ times its Jacobian θ
    @test lp ≈ sum(logpdf(d, exp(xi)) + xi for (d, xi) in zip(values(pr), x))
    narrow = merge(pr, (τ = LogNormal(log(300), 0.1),))
    @test LDP.logdensity(KnappPosterior(m; prior = narrow), x) - LDP.logdensity(post, x) ≈
          logpdf(Normal(log(300), 0.1), x[end]) - logpdf(Normal(pr.τ.μ, pr.τ.σ), x[end])
    @test_throws ArgumentError KnappPosterior(m; prior = Base.structdiff(pr, NamedTuple{(:τ,)}))
    @test_throws ArgumentError KnappPosterior(m; prior = merge(pr, (τ = Normal(),)))
    @test_throws DimensionMismatch LDP.logdensity(post, x[1:3])

    # pointwise log-likelihood: one entry per informative swab, summing to the log-likelihood
    obs = observed_swabs(m)
    @test length(obs) == 371 && issorted(obs)
    ll = pointwise_loglikelihood(m, x)
    @test length(ll) == length(obs) && all(isfinite, ll)
    @test sum(ll) ≈ LDP.logdensity(m, x)
    @test pointwise_loglikelihood(post, x) == ll
    X = permutedims(hcat(x, x .+ 0.01))
    LL = pointwise_loglikelihood(m, X)
    @test size(LL) == (2, length(obs)) && LL[1, :] == ll
    @test all(==(-Inf), pointwise_loglikelihood(m, fill(NaN, length(x))))
    @test LDP.logdensity(m, fill(NaN, length(x))) == -Inf
    @test_throws DimensionMismatch pointwise_loglikelihood(m, X[:, 1:3])

    # PSIS-LOO against PSIS.jl 0.9.10 (an independent implementation) on a fixed matrix
    S = 200
    L = [-(1 + i) * abs(sin(0.37s + i))^(1 + i) - 0.1i for s in 1:S, i in 1:3]
    L[:, 3] .= [-8 * (s / S)^6 for s in 1:S]
    loo = psis_loo(L)
    @test loo.k ≈ [-1.505769449462016, -1.0702811157670455, 1.3715272111300607] rtol = 1e-8
    @test loo.pointwise ≈ [-1.3441972020050397, -2.033095918466973, -4.184696388726367] rtol = 1e-8
    @test loo.elpd ≈ sum(loo.pointwise)
    lpd = [log(sum(exp.(L[:, i])) / S) for i in 1:3]
    @test loo.p ≈ sum(lpd) - loo.elpd
    @test_throws ArgumentError psis_loo(L[1:5, :])
    # identical draws: leaving one out changes nothing
    flat = repeat([-1.0 -2.0], 50)
    @test psis_loo(flat).pointwise ≈ [-1.0, -2.0]
    @test psis_loo(flat).p ≈ 0 atol = 1e-12

    # WAIC: log pointwise predictive density minus the posterior variance
    w = waic(L)
    v = [sum(abs2, L[:, i] .- sum(L[:, i]) / S) / (S - 1) for i in 1:3]
    @test w.pointwise ≈ lpd .- v && w.p ≈ sum(v)
    @test waic(flat).elpd ≈ -3.0
    d = elpd_difference(loo, loo)
    @test d.difference == 0 && d.se == 0
    @test elpd_difference(w, loo).difference ≈ w.elpd - loo.elpd
    @test_throws DimensionMismatch elpd_difference(loo, waic(L[:, 1:2]))

    # Turing extension: the log joint is the posterior log density, and NUTS runs
    tm = turing_model(post)
    @test Turing.DynamicPPL.logjoint(tm, (; x = x)) ≈ LDP.logdensity(post, x)
    mn = KnappLoadModel(:none)
    xn = Chytridiomycosis.maximum_likelihood(mn; starts = 1, iters = 1500).x
    chn = Turing.sample(Xoshiro(1), turing_model(KnappPosterior(mn)), Turing.NUTS(10, 0.8), 10;
        initial_params = Turing.InitFromParams((; x = xn)), progress = false)
    A = Array(chn)
    @test size(A, 1) == 10 && length(A) == 10 * length(xn)
    @test all(isfinite, A)
end
