# Posterior samples of the five hierarchical TPC fits of Gajewski et al. (2021, Ecol Evol
# 11:17920, doi:10.1002/ece3.8379), and their summaries against the paper's Tables 1, 2 and
# Appendix S1 Tables B.3-B.7.
#
# Needs the authors' data (Zenodo doi:10.5281/zenodo.4525954, extracted; see
# data/gajewski2021/README.md and gajewski2021_path()).
#
# Invocation (about 10-30 min per fit):
#   julia --threads=4 --project scripts/gajewski2021_mcmc.jl [form ...]
# with forms among briere2 ratkowsky ikemoto logan10 stinner (default: all five).
#
# Sampler: adaptive random-walk Metropolis (the package's Wilber 2022 sampler) started at the
# posterior mode, 4 chains. The authors used JAGS (5 chains); the posterior is the same, the
# sampler is not. Prints, for each fit, the posterior medians and 95% HDIs of the parameters,
# and the medians of the curve summaries (Topt, maximum r, adjusted Tmin and Tmax) over 1000
# draws, as the paper computes them.
using Chytridiomycosis, LogDensityProblems, Random, Statistics
const C = Chytridiomycosis

forms = isempty(ARGS) ? collect(C.GAJEWSKI_FORMS) : Symbol.(ARGS)
data = load_gajewski2021().constant
Threads.@threads for form in forms
    m = GajewskiModel(; form, data)
    fit = posterior_mode(m)
    Σ0 = Matrix(C.Diagonal((1e-3 .* max.(abs.(fit.x), 1e-3)) .^ 2))
    res = C._w22_metropolis(x -> LogDensityProblems.logdensity(m, x), fit.x, Σ0,
        Random.Xoshiro(2021); chains = 4, warmup = 20_000, iterations = 40_000, thin = 20)
    D = res.draws
    io = IOBuffer()
    println(io, "== $form: acceptance ", round.(res.acceptance; digits = 2),
        ", max R-hat ", round(maximum(res.rhat); digits = 3))
    for (j, n) in enumerate(parameter_names(m))
        lo, hi = C._hdi(D[:, j])
        println(io, rpad(n, 6), " median ", round(median(D[:, j]); sigdigits = 5),
            "  HDI (", round(lo; sigdigits = 5), ", ", round(hi; sigdigits = 5), ")")
    end
    idx = round.(Int, range(1, size(D, 1); length = 1000))
    S = map(idx) do i
        θ = merge(C._gj_named(m, D[i, :]), (R = m.R,))
        c = gajewski_curve(form, θ)
        To = thermal_optimum(c; range = (0.0, 45.0))
        lo, hi = adjusted_limits(c)
        (To, c(To), lo, hi)
    end
    for (k, n) in enumerate(("Topt", "r_max", "Tmin_adj", "Tmax_adj"))
        println(io, rpad(n, 9), " median ", round(median(getindex.(S, k)); digits = 3))
    end
    print(String(take!(io)))
end
