# Bayesian fit of the frog load equation to the Knapp et al. 2022 swabs (NEW analysis, not a
# published result): NUTS for each mechanism, convergence diagnostics, posterior summaries and
# PSIS-LOO / WAIC comparison. This is the run reported in docs/src/api/knapp_bayes.md and the
# wiki page models/new-knapp-load-fit.md ("Bayesian fit"). About 10 minutes with 4 threads.
#
# Needs an environment with this package (dev'ed), Turing and MCMCDiagnosticTools:
#   julia --project=<env> -e 'using Pkg; Pkg.develop(path = "."); Pkg.add(["Turing", "MCMCDiagnosticTools"])'
#   julia --threads=4 --project=<env> scripts/knapp_bayes.jl [outdir]
using Chytridiomycosis, Random, Printf, DelimitedFiles
using Turing: Turing, NUTS, MCMCThreads, sample, InitFromParams
using MCMCDiagnosticTools: ess, rhat

const C = Chytridiomycosis
outdir = isempty(ARGS) ? mktempdir() : mkpath(ARGS[1])
mechs = (:none, :decay, :applied, :bacteria)
# target acceptance: 0.99 removes the divergences of `applied` (16 at 0.9). `bacteria` keeps
# divergences at either value (84 at 0.9, 154 at 0.99): with its effect switched off (Δc → 0),
# h and K are free and the posterior is a funnel. Its elpd is the same in both runs.
δ = (none = 0.9, decay = 0.9, applied = 0.99, bacteria = 0.9)
nchains, ndraws, nadapt = 4, 1000, 1000
quant(v, p) = (s = sort(v); s[clamp(round(Int, p * length(s)), 1, length(s))])

res = Dict{Symbol,Any}()
for mech in mechs
    m = KnappLoadModel(mech)
    post = KnappPosterior(m)
    names = parameter_names(m)
    fit = C.maximum_likelihood(m)                       # starting points for the chains
    rng = Xoshiro(2026)
    inits = [InitFromParams((; x = fit.x .+ 0.3 .* randn(rng, length(fit.x)))) for _ in 1:nchains]
    t = @elapsed chn = sample(Xoshiro(1), turing_model(post), NUTS(nadapt, δ[mech]),
        MCMCThreads(), ndraws, nchains; initial_params = inits, progress = false)
    A = Array(chn)                                      # draw × chain × parameter, log scale
    R, E, Et = rhat(A), ess(A; kind = :bulk), ess(A; kind = :tail)
    ndiv = sum(Array(chn[:numerical_error]))
    X = reshape(A, :, length(names))
    writedlm(joinpath(outdir, "draws_$(mech).tsv"), vcat(permutedims(string.(collect(names))), X))
    ll = pointwise_loglikelihood(m, X)
    loo, w = psis_loo(ll), waic(ll)
    res[mech] = (; m, fit, X, loo, waic = w)
    @printf("\n=== %s: NUTS %.0f s, %d divergences, max R-hat %.3f, min bulk ESS %.0f, min tail ESS %.0f\n",
        mech, t, ndiv, maximum(R), minimum(E), minimum(Et))
    @printf("ML log L %.1f | elpd_loo %.1f (se %.1f), p_loo %.1f, max k %.2f | elpd_waic %.1f\n",
        fit.loglik, loo.elpd, loo.se, loo.p, maximum(loo.k), w.elpd)
    prior = knapp_prior(m)
    println("param        ML      median        2.5%       97.5%  prior 2.5% prior 97.5%  sd/prior sd")
    for (j, n) in enumerate(names)
        v = exp.(X[:, j]); d = prior[n]
        sd = sqrt(sum(abs2, X[:, j] .- sum(X[:, j]) / size(X, 1)) / (size(X, 1) - 1))
        @printf("%-4s %11.3g %11.3g %11.3g %11.3g %11.3g %11.3g  %.2f\n", n, fit.estimates[n],
            quant(v, 0.5), quant(v, 0.025), quant(v, 0.975), exp(d.μ - 1.96d.σ),
            exp(d.μ + 1.96d.σ), sd / d.σ)
    end
end

println("\n=== Mechanisms compared by PSIS-LOO (elpd relative to the best; ΔAIC from ML)")
best = argmax(k -> res[k].loo.elpd, collect(mechs))
bestaic = minimum(res[k].fit.aic for k in mechs)
for k in mechs
    d = elpd_difference(res[k].loo, res[best].loo)
    @printf("%-9s ΔAIC %6.1f | elpd_loo %8.1f (se %.1f) | Δelpd %6.1f (se %.1f)\n", k,
        res[k].fit.aic - bestaic, res[k].loo.elpd, res[k].loo.se, d.difference, d.se)
end
