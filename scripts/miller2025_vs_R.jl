# Cross-check of the Julia refugia CTMC (src/models/miller2025.jl) against the authors' R code
# (Miller et al. 2025, Zenodo 10.5281/zenodo.15179344 v1.0, GillespieSSA 0.6.2).
#
# 1. Run the R driver (sources the authors' files read-only):
#      Rscript scripts/miller2025_vs_R.R <lab/papers/code/miller2025/code> <outdir> 10000 4
# 2. Then: julia --threads=4 --project scripts/miller2025_vs_R.jl <outdir>
#
# For each parameter set ("paper" = Table 3 means, "synthetic" = Table S2, "high"), the two
# implementations' runs are compared at each environment, state and survey week (not week
# 0, which is fixed):
# - means: two-sample z statistic;
# - distributions: 10%, 50% and 90% quantiles, and a two-sample Kolmogorov–Smirnov test;
# - the authors' ABC summary statistic (Eq. 2 with week 0, as coded) per environment and in
#   total.
# The runs are independent, so under agreement the z statistics are about N(0, 1) and the
# KS p-values about uniform.
using Chytridiomycosis, Random, Statistics, Printf

const C = Chytridiomycosis
dir = isempty(ARGS) ? error("usage: julia scripts/miller2025_vs_R.jl <outdir of the R driver>") : ARGS[1]

function readcsv(path)
    lines = readlines(path)
    h = split(lines[1], ',')
    rows = [split(l, ',') for l in lines[2:end] if !isempty(l)]
    return Dict(String(h[j]) => [r[j] for r in rows] for j in eachindex(h))
end

# Asymptotic two-sample KS p-value (Kolmogorov distribution), adequate for n, m ≥ 1000.
function ks2(x, y)
    xs, ys = sort(x), sort(y)
    vals = sort!(unique(vcat(xs, ys)))
    D = maximum(abs(searchsortedlast(xs, v) / length(xs) - searchsortedlast(ys, v) / length(ys))
                for v in vals)
    ne = length(xs) * length(ys) / (length(xs) + length(ys))
    λ = (sqrt(ne) + 0.12 + 0.11 / sqrt(ne)) * D
    p = λ < 1e-3 ? 1.0 : clamp(2 * sum((-1)^(k - 1) * exp(-2k^2 * λ^2) for k in 1:100), 0.0, 1.0)
    return D, p
end

names_R = ("S", "I1", "I2", "I3", "SV", "IV1", "IV2", "IV3")
sets = Dict("paper" => (β_sh = 0.675, β_un = 0.401, α = 0.072, ω = 0.569, μ = 0.021),
            "synthetic" => (β_sh = 0.05, β_un = 0.02, α = 0.177, ω = 0.258, μ = 0.029),
            "high" => (β_sh = 1.5, β_un = 1.0, α = 0.3, ω = 2.0, μ = 0.1))
# The "paper" and "synthetic" values are the compendium's; check so the driver and package agree.
let q = parameter_set(:miller2025), s = parameter_set(:miller2025; set = "synthetic")
    @assert (q.β_sh, q.β_un, q.α, q.ω, q.μ) == Tuple(sets["paper"])
    @assert (s.β_sh, s.β_un, s.α, s.ω, s.μ) == Tuple(sets["synthetic"])
end

data = load_miller2025(:exp)
weeks = data.weeks
allz = Float64[]
allp = Float64[]
println("R vs Julia: refugia CTMC (Miller et al. 2025)\n")
for set in ("paper", "synthetic", "high")
    st = readcsv(joinpath(dir, "R_states_$set.csv"))
    er = readcsv(joinpath(dir, "R_errors_$set.csv"))
    nR = length(er["iter_id"])
    θ = sets[set]
    p = RefugiaParameters(:miller2025; θ...)
    rng = Xoshiro(2025)
    sims = Dict(env => [solve(RefugiaCTMCProblem(p; environment = env, saveat = weeks); rng).u
                        for _ in 1:nR] for env in (:shaded, :unshaded))
    @printf("== %s: %s, %d R runs and %d Julia runs per environment\n", set, θ, nR, nR)
    t = parse.(Float64, st["t"])
    sh = st["shaded"] .== "TRUE"
    zs = Float64[]
    worst = (0.0, "")
    for (env, flag) in ((:shaded, true), (:unshaded, false)), (i, nm) in enumerate(names_R),
            (k, w) in enumerate(weeks)
        k == 1 && continue
        xr = parse.(Float64, st[nm][(sh .== flag) .& (t .== w)])
        xj = Float64[u[i, k] for u in sims[env]]
        se = sqrt(var(xr) / length(xr) + var(xj) / length(xj))
        se == 0 && continue
        z = (mean(xj) - mean(xr)) / se
        push!(zs, z)
        abs(z) > abs(worst[1]) && (worst = (z, "$env $nm week $(Int(w)): R mean $(round(mean(xr); digits=3)), Julia $(round(mean(xj); digits=3))"))
        _, pks = ks2(xr, xj)
        push!(allp, pks)
    end
    append!(allz, zs)
    @printf("  state means: %d comparisons, mean z = %.3f, sd z = %.3f, max |z| = %.2f (%s)\n",
            length(zs), mean(zs), std(zs), abs(worst[1]), worst[2])
    @printf("  |z| > 3: %d (expected %.1f under agreement)\n", count(>(3), abs.(zs)),
            length(zs) * 0.0027)
    # infected totals at week 15, quantiles
    for (env, flag) in ((:shaded, true), (:unshaded, false))
        r15 = (sh .== flag) .& (t .== 15)
        IR = parse.(Float64, st["I1"][r15]) .+ parse.(Float64, st["I2"][r15]) .+
             parse.(Float64, st["I3"][r15]) .+ parse.(Float64, st["IV1"][r15]) .+
             parse.(Float64, st["IV2"][r15]) .+ parse.(Float64, st["IV3"][r15])
        IJ = Float64[sum(u[[2, 3, 4, 6, 7, 8], end]) for u in sims[env]]
        @printf("  infected at week 15, %-8s R mean %.3f q(10,50,90) %s | Julia mean %.3f q %s\n",
                env, mean(IR), quantile(IR, [0.1, 0.5, 0.9]), mean(IJ),
                quantile(IJ, [0.1, 0.5, 0.9]))
    end
    # the ABC summary statistic
    eR_sh, eR_un = parse.(Float64, er["err_shaded"]), parse.(Float64, er["err_unshaded"])
    dJ = [refugia_distance(data, a, b) for (a, b) in zip(sims[:shaded], sims[:unshaded])]
    eJ_sh, eJ_un = Float64[x.shaded for x in dJ], Float64[x.unshaded for x in dJ]
    for (lab, a, b) in (("shaded", eR_sh, eJ_sh), ("unshaded", eR_un, eJ_un),
                        ("total", eR_sh .+ eR_un, eJ_sh .+ eJ_un))
        D, pk = ks2(a, b)
        @printf("  Δ %-8s R mean %8.1f sd %6.1f q(2.5,50) %s | Julia mean %8.1f sd %6.1f q %s | KS D = %.4f p = %.3f\n",
                lab, mean(a), std(a), quantile(a, [0.025, 0.5]), mean(b), std(b),
                quantile(b, [0.025, 0.5]), D, pk)
    end
    println()
end
@printf("All sets: %d mean comparisons, mean z %.3f, sd z %.3f, max |z| %.2f; ",
        length(allz), mean(allz), std(allz), maximum(abs.(allz)))
@printf("%d KS tests of state counts, %d with p < 0.01 (expected %.1f)\n", length(allp),
        count(<(0.01), allp), 0.01 * length(allp))
