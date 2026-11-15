# Reduced maximum-likelihood refit of Hollanders et al. (2023) to the bundled Dryad data
# (OUR fit of a published model; not published results). Writes
# data/hollanders2023/refit.toml, which hollanders2023_fit() reads so the validation targets
# do not re-optimise inside Pkg.test(). Deterministic: BFGS from hollanders2023_start().
#
#   JULIA_NUM_THREADS=2 julia --project scripts/hollanders2023_refit.jl
using Chytridiomycosis, Dates
const C = Chytridiomycosis

t = @elapsed fit = hollanders2023_fit(; refit = true)
path = C._hollanders2023_write_fit(hollanders2023_refit_path(), fit;
    comment = "OUR reduced maximum-likelihood refit of Hollanders et al. (2023), not published results.\n" *
              "Written by scripts/hollanders2023_refit.jl ($(Dates.today()), Julia $(VERSION), " *
              "$(Threads.nthreads()) thread(s), $(round(t; digits = 0)) s): CTMultistateProblem() defaults,\n" *
              "BFGS from hollanders2023_start(); vcov is the inverse finite-difference observed information.")
println("loglik ", fit.loglik, ", iterations ", fit.iterations, ", max |gradient| ",
        maximum(abs, fit.gradient), "; wrote ", path, " in ", round(t; digits = 0), " s")
