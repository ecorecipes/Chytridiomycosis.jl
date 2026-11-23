# Compare the package's port of the Barnett et al. (2026, bioRxiv preprint) agent-based model
# with every BehaviorSpace run the authors released (12 files, 29,325 runs). For each author
# run the package runs the same parameters `reps` times (default 8); each output is compared by
# file: means, standard error of the difference, medians and the Kolmogorov–Smirnov distance.
# Needs the authors' repository (barnett2026_path()); writes data/barnett2026/reproduction.toml.
#
#   julia --threads=8 --project scripts/barnett2026_abm.jl [reps]
#
# About 25 minutes of CPU time per replicate set / 8 threads = 4–8 minutes at reps = 8.
# Conventions: −100 is missing in the authors' newer files. In the older ones (column metas_prev)
# intensity and prevalence are 0 when fewer than two metamorphs carry zoosporangia; the
# package's prevalence is set to 0 in the same runs, and intensity is compared where it is
# defined. The newer files round the Bd death fraction to 3 decimals, so the package's is too.
using Chytridiomycosis
using Random: Xoshiro
using Statistics: mean, median, var
using TOML

const REPS = isempty(ARGS) ? 8 : parse(Int, ARGS[1])
const C = Chytridiomycosis

function read_table(path)
    lines = readlines(path)
    k = findfirst(l -> occursin("[run number]", l), lines)
    cells(l) = strip.(split(replace(l, '"' => "", '﻿' => ""), ','))
    h = cells(lines[k])
    rows = [cells(l) for l in lines[(k + 1):end] if !isempty(strip(l))]
    return Dict(name => [something(tryparse(Float64, r[i]), NaN) for r in rows]
                for (i, name) in reverse(collect(enumerate(h))))
end

function ks(x, y)
    sx, sy = sort(x), sort(y)
    return maximum(abs(searchsortedlast(sx, z) / length(sx) - searchsortedlast(sy, z) / length(sy))
                   for z in unique([sx; sy]))
end

const FILES = [
    "control" => "Chytrid_Vaccine_ABM_November no_vaccination-table.csv",
    "infection" => "one-way/Chytrid_Vaccine_ABM_November c_est_vs_coverage-table.csv",
    "transmission" => "one-way/Chytrid_Vaccine_ABM_November c_shedding_vs_coverage-table.csv",
    "growth" => "one-way/Chytrid_Vaccine_ABM_November c_clear_vs_coverage-table.csv",
    "tolerance" => "one-way/Chytrid_Vaccine_ABM_November c_smax_vs_coverage-table.csv",
    "tolerance_infection" => "two-way/anti_infection_vs_tolerance.csv",
    "tolerance_transmission" => "two-way/anti_transmission_vs_tolerance.csv",
    "tolerance_growth" => "two-way/anti_growth_vs_tolerance.csv",
    "adverse_infection" => "hypothetical adverse vaccination/increasing_establishment.csv",
    "adverse_transmission" => "hypothetical adverse vaccination/increased_transmission.csv",
    "adverse_growth" => "hypothetical adverse vaccination/decreasing_clearance.csv",
    "adverse_tolerance" => "hypothetical adverse vaccination/decreasing_tolerance.csv"]

p0 = BarnettParameters(:barnett2026)
result = Dict{String,Any}("reps" => REPS, "commit" => "29f4af6")
for (fi, (name, file)) in enumerate(FILES)
    d = read_table(joinpath(barnett2026_path(), "model output files", split(file, '/')...))
    n = length(d["count metamorphs"])
    old = haskey(d, "metas_prev")
    runs = Matrix{Any}(undef, n, REPS)
    Threads.@threads for q in 1:(n * REPS)
        i, k = (q - 1) % n + 1, (q - 1) ÷ n + 1
        p = C.with(p0; relative_variation = d["relative_variation"][i])
        v = BarnettVaccine(; coverage = d["v-coverage"][i], c_est = d["c_est"][i],
            c_shedding = d["c_shedding"][i], c_clear = d["c_clear"][i], c_smax = d["c_smax"][i])
        runs[i, k] = solve(BarnettABMProblem(p, v); rng = Xoshiro(hash((fi, q))))
    end
    columns = ["survivors" => "count metamorphs", "bd_death_fraction" => "prop_deaths_due_to_bd",
        "median_intensity" => old ? "median_spn_inten_metas" : "median_spn_inten_all_metas",
        "mean_intensity" => old ? "avg_spn_inten_metas" : "avg_spn_inten_all_metas",
        "unvaccinated_median_intensity" => "median_spn_inten_unvax_metas",
        "vaccinated_median_intensity" => "median_spn_inten_vax_metas",
        "prevalence" => old ? "metas_prev" : "all_metas_prev",
        "unvaccinated_fraction" => "prop_surviving_unvax",
        "zoospore_density" => "total-zsp-density"]
    res = Dict{String,Any}("author_runs" => n)
    for (out, colname) in columns
        haskey(d, colname) || continue
        x = [v == -100 ? NaN : v for v in d[colname]]
        old && occursin("intensity", out) && (x = [v == 0 ? NaN : v for v in x])
        old && out == "prevalence" &&
            (x = [d["count metamorphs"][i] == 0 ? NaN : x[i] for i in 1:n])
        y = map(vec(runs)) do r
            v = Float64(getfield(r, Symbol(out)))
            old && out == "prevalence" && r.survivors > 0 && isnan(r.median_intensity) && (v = 0.0)
            !old && out == "bd_death_fraction" ? round(v; digits = 3) : v
        end
        x, y = filter(!isnan, x), filter(!isnan, y)
        (length(x) < 2 || length(y) < 2) && continue
        se = sqrt(var(x) / length(x) + var(y) / length(y))
        se > 0 || continue
        D = ks(x, y)
        m = length(x) * length(y) / (length(x) + length(y))
        res[out] = Dict("author_mean" => mean(x), "package_mean" => mean(y), "se" => se,
            "z" => round((mean(y) - mean(x)) / se; digits = 2), "author_median" => median(x),
            "package_median" => median(y), "ks" => round(D; digits = 4),
            "ks_p" => round(min(1.0, 2exp(-2 * D^2 * m)); sigdigits = 2))
        println(rpad(name, 24), rpad(out, 30), " ", round(mean(x); digits = 3), " / ",
                round(mean(y); digits = 3), "  z ", res[out]["z"], "  medians ",
                round(median(x); digits = 2), " / ", round(median(y); digits = 2), "  KS ",
                res[out]["ks"], " p ", res[out]["ks_p"])
        flush(stdout)
    end
    result[name] = res
end
zs = [r["z"] for s in values(result) if s isa Dict for r in values(s) if r isa Dict]
result["comparisons"] = length(zs)
result["z_above_3"] = count(>(3), abs.(zs))
dir = joinpath(pkgdir(Chytridiomycosis), "data", "barnett2026")
mkpath(dir)
open(joinpath(dir, "reproduction.toml"), "w") do io
    println(io, "# Written by scripts/barnett2026_abm.jl: authors' BehaviorSpace runs against the port.")
    TOML.print(io, result; sorted = true)
end
println(length(zs), " comparisons, ", result["z_above_3"], " with |z| > 3")
