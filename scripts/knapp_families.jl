# Compare the distribution of Bd swab loads around the median in the N16 fit of the frog load
# equation to the Knapp et al. 2022 Dusy Basin swabs (NEW MODEL, not a published analysis).
# Every mechanism is fitted by maximum likelihood with a lognormal, gamma and Weibull Bd
# observation (same parameters, so log-likelihoods compare directly), from the same starts as
# the original fit. Writes data/knapp2022/family_comparison.toml.
#
#   julia --project -t 1 scripts/knapp_families.jl
using Chytridiomycosis, Random, TOML

const STARTS = 8
const ITERS = 2000
out = Dict{String,Any}("comment" => "OUR fits (N16, a new model), not published results; " *
    "maximum likelihood from $STARTS starts, $ITERS Nelder-Mead iterations each; " *
    "written by scripts/knapp_families.jl")
for mech in (:none, :decay, :applied, :bacteria)
    res = Dict{String,Any}()
    for fam in (:lognormal, :gamma, :weibull)
        Random.seed!(20120719)
        t = @elapsed fit = maximum_likelihood(KnappLoadModel(mech; family = fam);
                                              starts = STARTS, iters = ITERS)
        k = length(fit.x)
        res[string(fam)] = Dict("loglik" => fit.loglik, "aic" => fit.aic, "parameters" => k,
            "sigma" => fit.estimates.σ, "iota" => fit.estimates.ι,
            "starts_within_0.5" => count(l -> l > fit.loglik - 0.5, fit.logliks),
            "seconds" => round(t; digits = 1))
        println(rpad(string(mech), 9), rpad(string(fam), 10), " loglik ", round(fit.loglik; digits = 1),
                "  AIC ", round(fit.aic; digits = 1), "  σ ", round(fit.estimates.σ; sigdigits = 3),
                "  (", count(l -> l > fit.loglik - 0.5, fit.logliks), "/", STARTS, " starts, ",
                round(t; digits = 0), " s)")
        flush(stdout)
    end
    out[string(mech)] = res
end
path = joinpath(pkgdir(Chytridiomycosis), "data", "knapp2022", "family_comparison.toml")
open(io -> TOML.print(io, out), path, "w")
println("wrote ", path)
