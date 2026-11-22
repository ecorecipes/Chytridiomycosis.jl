# Fit the two-trait Bd–bacteria IPM (NEW MODEL) to the joint Bd and J. lividum loads of the
# treated Knapp et al. (2022) cohort, with and without the bacteria acting on Bd growth, and
# profile κ_growth. A new fit, not a published analysis; see src/fitting/bacteria_ipm.jl.
#
#   julia --project scripts/bacteria_ipm_knapp_fit.jl            # about a minute
#
# Prints the estimates of both versions, their log-likelihoods and AIC, a likelihood-ratio
# statistic, and the profile log-likelihood of κ_growth.
using Chytridiomycosis, Printf, Random

function main()
    fits = Dict{Bool,Any}()
    for coupled in (false, true)
        m = KnappJointModel(; coupled)
        f = maximum_likelihood(m; starts = 6, iters = 3000, spread = 0.7)
        fits[coupled] = f
        @printf("\ncoupled = %s: %d swabs, log-likelihood %.2f, AIC %.2f (starts: %s)\n",
            coupled, f.n, f.loglik, f.aic, join((@sprintf("%.2f", l) for l in f.logliks), ", "))
        for (k, v) in pairs(f.estimates)
            @printf("  %-14s %10.4f\n", k, v)
        end
        @printf("  J. lividum half-life %.2f days\n", log(2) / f.estimates.decline)
    end
    lr = 2 * (fits[true].loglik - fits[false].loglik)
    @printf("\nlikelihood ratio (coupled vs not): %.2f on 2 parameters; ΔAIC = %.2f\n", lr,
        fits[true].aic - fits[false].aic)
    # profile κ_growth: re-fit the other parameters at each value
    m = KnappJointModel(; coupled = true)
    best = fits[true]
    println("\nprofile of κ_growth (log-likelihood, difference from the maximum):")
    for κ in (-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 4.0)
        f(z) = -Chytridiomycosis.LogDensityProblems.logdensity(m, [z[1:6]; κ; z[7]])
        z0 = [best.x[1:6]; best.x[8]]
        z, fz = Chytridiomycosis.nelder_mead(f, z0; iters = 3000)
        z, fz = Chytridiomycosis.nelder_mead(f, z; iters = 1500, step = 0.2)
        @printf("  κ_growth = %5.1f: %9.2f  (%6.2f)\n", κ, -fz, -fz - best.loglik)
    end
end

main()
