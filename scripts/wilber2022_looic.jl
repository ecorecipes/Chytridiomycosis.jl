# LOOIC of the four Table 1 models of Wilber et al. (2022) (abundance covariate) at the three
# prevalence cutoffs, each the mean over four sampler seeds, written to
# data/wilber2022/looic.toml for the validation targets wilber2022_looic_*.
#
# PSIS-LOO on 79 lakes has a Monte Carlo error of a few tenths to about one LOOIC unit per run
# (a lake or two with Pareto k near or above 0.7), so a single seed can land 2–3 units away
# (seed 2022 of the log-tadpole model at ρ = 0.5: k = 1.32, 1043.3 against a 4-seed mean of
# 1040.2). Averaging four seeds halves our part of the error; the authors' single Stan run
# keeps its own.
#
# Needs the authors' deposit (data/wilber2022/README.md). Run with
#     julia --project=. scripts/wilber2022_looic.jl
# Single-threaded and deterministic: Xoshiro(seed), seeds 1–4, the default solve() settings.
# Takes about 15 minutes on an idle machine (an hour at load 200).

using Chytridiomycosis, Random, TOML

const OUT = joinpath(pkgdir(Chytridiomycosis), "data", "wilber2022", "looic.toml")
const MODELS = (:full, :log_tadpole, :linear_temp, :quadratic_temp)
const CUTOFFS = (0.25, 0.3333, 0.5)
const SEEDS = 1:4

results = Dict{String,Any}()
for model in MODELS
    tab = Dict{String,Any}()
    for ρ in CUTOFFS
        prob = OutbreakHMMProblem(; cutoff = ρ, model)
        vals, ks = Float64[], Float64[]
        for seed in SEEDS
            fit = solve(prob; rng = Xoshiro(seed))
            l = psis_loo(fit)
            push!(vals, round(-2l.elpd; digits = 2))
            push!(ks, round(maximum(l.k); digits = 2))
        end
        m = sum(vals) / length(vals)                     # of the rounded values
        s = sqrt(sum(abs2, vals .- m) / (length(vals) - 1))
        tab["rho" * Chytridiomycosis._w22_rho_label(ρ)] =
            Dict("cutoff" => ρ, "seeds" => collect(SEEDS), "looic" => vals,
                 "max_pareto_k" => ks, "mean" => round(m; digits = 2), "sd" => round(s; digits = 2))
        println("$model ρ=$ρ ", vals, " mean ", round(m; digits = 2), " sd ", round(s; digits = 2))
        flush(stdout)
    end
    results[string(model)] = tab
end

open(OUT, "w") do io
    println(io, "# Wilber et al. (2022) Table 1 models: LOOIC (−2 elpd_loo, PSIS-LOO over 79 lakes) per")
    println(io, "# sampler seed (Xoshiro(1)–Xoshiro(4), default solve settings), written by")
    println(io, "# scripts/wilber2022_looic.jl. Each LOOIC is rounded to 0.01 before the mean and sd.")
    println(io, "# Our refits, not published data.")
    TOML.print(io, results; sorted = true)
end
