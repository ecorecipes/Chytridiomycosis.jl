# N5 (new model, not a published result): identifiability of the Briggs et al. 2010 load
# model from swabs, eDNA and load aggregation. Reproduces the tables on the wiki page
# lab/wiki/models/new-swab-edna-identifiability.md.
#
#   julia --project scripts/n5_identifiability.jl            # from the package root
#   N5_MC=0 julia --project scripts/n5_identifiability.jl    # skip the Monte Carlo (~4 min)
using Chytridiomycosis, Random, LinearAlgebra, Statistics
const LDP = Chytridiomycosis.LogDensityProblems

p = BriggsParameters(:briggs2010)                     # Fig. 2C values from the compendium
# Assay scales: ASSUMED, illustrative only (no published values for the Briggs model).
obs = SurveyObservation(ρ = 10.0, κ = 1.0, σe = 0.5)
Z0 = 1e4                                              # zoospores seeded on day 0 (assumed)
days = [5.0, 10, 15, 20, 30, 40]
B = (:γ, :ν, :f, :μ, :η, :σ)
Zs(n) = Tuple(Symbol("Z0_", i) for i in 1:n)
design(; Ns = [10, 40, 60], swabs = 10, pool = 1, edna = 3, days = days) =
    [SurveyPond(; N, days, swabs = min(swabs, N ÷ pool), pool, edna) for N in Ns]
model(ps, free) = SurveyModel(SurveySample[], ps; p, obs, Z0 = fill(Z0, length(ps)), free)

println("N* = ", round(host_density_threshold(p), digits = 2), "; combinations: ",
    identifiable_combinations(p))
for N in (10, 40, 60)
    s = solve(BriggsMomentProblem(p; N, Z0); saveat = days)
    println("N = $N: mean sporangia ", round.(s.m, sigdigits = 3), ", k ", round.(s.k, sigdigits = 3))
end

println("\n== Structural: rank of the output Jacobian ==")
for (label, ps, free, agg) in [
        ("A swab means, one density (N = 40)", design(; Ns = [40], edna = 0), (B..., :ρ, Zs(1)...), false),
        ("B swab means, three densities", design(; edna = 0), (B..., :ρ, Zs(3)...), false),
        ("C eDNA only, three densities", design(; swabs = 0), (B..., :κ, Zs(3)...), true),
        ("D swab means + eDNA, both scales unknown", design(), (B..., :ρ, :κ, Zs(3)...), false),
        ("E swab means + eDNA, κ known", design(), (B..., :ρ, Zs(3)...), false),
        ("F individual swabs + eDNA, κ known", design(), (B..., :ρ, Zs(3)...), true),
        ("G individual swabs + eDNA, κ and ρ known", design(), (B..., Zs(3)...), true),
        ("H individual swabs only, ρ unknown", design(; edna = 0), (B..., :ρ, Zs(3)...), true),
        ("I individual swabs only, ρ known", design(; edna = 0), (B..., Zs(3)...), true),
        ("J everything at one density, ρ κ known", design(; Ns = [40]), (B..., Zs(1)...), true),
        ("K individual swabs + eDNA, ρ known, κ unknown", design(), (B..., :κ, Zs(3)...), true),
        ("L as G with two densities (10, 40)", design(; Ns = [10, 40]), (B..., Zs(2)...), true)]
    r = sensitivity_rank(model(ps, free); aggregation = agg)
    println(rpad(label, 50), " rank ", r.rank, " of ", r.nparams)
    for v in r.null
        println("    null: ", map(x -> round(x, digits = 2), v))
    end
end

println("\n== Practical: expected-information standard errors (log scale) ==")
function ses(label, ps, free; aggregation = true)
    I = expected_information(model(ps, free); aggregation)
    ev = eigvals(Symmetric(I))
    se = minimum(ev) / maximum(ev) > 1e-13 ? sqrt.(diag(inv(I))) : fill(Inf, length(free))
    println(rpad(label, 46), join(["$n=$(round(s, sigdigits = 2))" for (n, s) in zip(free, se)
                                   if n in B], "  "))
end
F = (B..., :σe)
ses("swab means + eDNA, ρ κ known", design(), (F..., Zs(3)...); aggregation = false)
ses("individual swabs, ρ known, no eDNA", design(; edna = 0), (B..., Zs(3)...))
ses("individual swabs + eDNA, ρ unknown", design(), (F..., :ρ, Zs(3)...))
ses("FULL: individual swabs + eDNA, ρ κ known", design(), (F..., Zs(3)...))
ses("FULL with two densities (10, 40)", design(; Ns = [10, 40]), (F..., Zs(2)...))
ses("FULL with two densities (40, 60)", design(; Ns = [40, 60]), (F..., Zs(2)...))
ses("FULL with μ known", design(), (filter(!=(:μ), F)..., Zs(3)...))
ses("FULL with 20 swabs per visit", design(; swabs = 20), (F..., Zs(3)...))
ses("FULL with 10 eDNA samples per visit", design(; edna = 10), (F..., Zs(3)...))
ses("FULL with 10 visits (days 4–40)", design(; days = 4.0:4:40), (F..., Zs(3)...))
ses("FULL with six ponds (N = 5…60)", design(; Ns = [5, 10, 20, 30, 40, 60]), (F..., Zs(6)...))
ses("six ponds, 10 visits, 10 eDNA", design(; Ns = [5, 10, 20, 30, 40, 60], days = 4.0:4:40,
    edna = 10), (F..., Zs(6)...))
for pool in (1, 2, 5, 10)
    ses("FULL, 10 frogs per visit in tubes of $pool", design(; swabs = 10 ÷ pool, pool),
        (F..., Zs(3)...))
end

if get(ENV, "N5_MC", "1") != "0"
    println("\n== Monte Carlo: maximum likelihood on 40 synthetic data sets (FULL design) ==")
    ps = design()
    free = (F..., Zs(3)...)
    se = sqrt.(diag(inv(expected_information(model(ps, free)))))
    est, gain = Vector{Vector{Float64}}(), Float64[]
    for rep in 1:40
        data = simulate_survey(Xoshiro(100 + rep), p, obs, ps; Z0 = fill(Z0, 3))
        m = SurveyModel(data, ps; p, obs, Z0 = fill(Z0, 3), free)
        fit = maximum_likelihood(m; starts = 2)
        x0 = log.(m.values[m.free])
        push!(gain, fit.loglik - LDP.logdensity(m, x0))
        push!(est, fit.x .- x0)
    end
    E = reduce(hcat, est)
    println("median 2Δloglik over the truth ", round(2median(gain), digits = 1),
        " (χ² with 10 df: median 9.3)")
    println("name   info SE   robust sd   sd   median error (all log scale)")
    for (j, n) in enumerate(free)
        e = E[j, :]
        println(rpad(n, 7), rpad(round(se[j], sigdigits = 2), 10),
            rpad(round(1.4826 * median(abs.(e .- median(e))), sigdigits = 2), 12),
            rpad(round(std(e), sigdigits = 2), 6), round(median(e), sigdigits = 2))
    end
end
