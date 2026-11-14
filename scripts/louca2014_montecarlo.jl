# Monte Carlo reproduction of Louca, Lampo & Doebeli (2014) Proc R Soc B 281:20132783:
# outcome fractions, risk effects (Table 2) and outcome probabilities against Q_o (Fig. 5) and
# against (β_A, μ_Z) (Fig. 4), for the four scenarios (tadpoles suffer/tolerate, Z_o = 0 or
# sampled from 10 to 10^15), with the deterministic DOI-distribution model.
#
#   julia --project -t 2 scripts/louca2014_montecarlo.jl [n_per_scenario] [dt] [exposure|start]
#
# Writes data/louca2014/reproduction.toml (summary) and data/louca2014/montecarlo.csv (one row
# per run: scenario, scaled draws u, Q_o, outcome, t_end). The paper used 10^4 runs per
# scenario and explicit RK2 with dt = 10^-3–10^-2 d; we default to 2000 runs and the
# semi-implicit scheme with dt = 0.05 d (outcome agreement with dt = 0.01 d is recorded).
# Parameter sets are drawn sequentially from one seeded RNG, so the result does not depend on
# the number of threads.
using Chytridiomycosis
using Random
using TOML
using Printf

const SEED = 20132783
n = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 2000
dt = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 0.05
n_check = min(200, n)          # runs re-solved with dt = 0.01 to check the step size
# When the reservoir arrives: at exposure (default; the population is healthy until the second
# year, §2b) or from t = 0 (`start`). The paper does not say.
onset = length(ARGS) >= 3 ? Symbol(ARGS[3]) : :exposure

scenarios = [("suffer_reservoir", :suffer, true), ("tolerate_reservoir", :tolerate, true),
             ("suffer_noreservoir", :suffer, false), ("tolerate_noreservoir", :tolerate, false)]

outdir = joinpath(dirname(louca2014_path()))
mkpath(outdir)
summary = Dict{String,Any}("seed" => SEED, "n_per_scenario" => n, "dt" => dt,
    "method" => "semi_implicit", "years" => 10, "threads" => Threads.nthreads(),
    "reservoir_onset" => string(onset))
rows = String[]
qedges = collect(-8.0:1.0:6.0)                       # log10 Q_o bins (Fig. 5 spans 1e-6–1e4)
βedges = collect(-10.0:1.0:-3.0)                     # log10 β_A bins (Fig. 4)
μedges = collect(-3.0:1.0:0.0)                       # log10 μ_Z bins (Fig. 4)

t_start = time()
for (k, (name, tad, res)) in enumerate(scenarios)
    rng = Xoshiro(SEED + k)
    samples = [louca_sample(rng; tadpoles = tad, reservoir = res) for _ in 1:n]
    out = Vector{Any}(undef, n)
    t0 = time()
    Threads.@threads for i in 1:n
        s = samples[i]
        sol = solve(LoucaDOIProblem(s.p; saveat = 3650.0, reservoir_onset = onset); dt)
        out[i] = (; p = s.p, u = s.u, Qo = reproduction_quotient(s.p), outcome = sol.outcome,
                  t_end = sol.t_end)
    end
    secs = time() - t0
    chk = Vector{Symbol}(undef, n_check)
    Threads.@threads for i in 1:n_check
        chk[i] = solve(LoucaDOIProblem(samples[i].p; saveat = 3650.0, reservoir_onset = onset);
                       dt = 0.01).outcome
    end
    agree = count(i -> chk[i] === out[i].outcome, 1:n_check)
    frac(o) = count(r -> r.outcome === o, out) / n
    re = louca_risk_effects(out)
    pq = louca_outcome_probabilities(out; edges = qedges)
    # Fig. 4: outcome fractions on a (log10 β_A, log10 μ_Z) grid
    grid = Dict{String,Any}()
    for o in (:extinction, :clearance, :persistence)
        M = fill(NaN, length(βedges) - 1, length(μedges) - 1)
        C = zeros(Int, size(M))
        A = zeros(size(M))
        for r in out
            i = searchsortedlast(βedges, log10(r.p.βA))
            j = searchsortedlast(μedges, log10(r.p.μZ))
            (1 <= i <= size(M, 1) && 1 <= j <= size(M, 2)) || continue
            C[i, j] += 1
            A[i, j] += r.outcome === o
        end
        grid[string(o)] = [[C[i, j] > 0 ? round(A[i, j] / C[i, j]; digits = 4) : -1.0
                            for j in axes(M, 2)] for i in axes(M, 1)]
        grid["n"] = [[C[i, j] for j in axes(M, 2)] for i in axes(M, 1)]
    end
    fin = count(r -> isfinite(r.Qo), out)
    nanfix(v) = [isnan(x) ? -1.0 : round(x; digits = 4) for x in v]
    summary[name] = Dict{String,Any}(
        "tadpoles" => string(tad), "reservoir" => res, "n" => n, "cpu_seconds" => round(secs; digits = 1),
        "extinction" => frac(:extinction), "clearance" => frac(:clearance),
        "persistence" => frac(:persistence), "Qo_finite" => fin,
        "dt_check_runs" => n_check, "dt_check_agree" => agree,
        "risk_effects" => Dict(string(s) => round(getfield(re, s); digits = 4) for s in keys(re)),
        "Qo_bins" => Dict("log10_edges" => qedges, "extinction" => nanfix(pq.extinction),
                          "clearance" => nanfix(pq.clearance),
                          "persistence" => nanfix(pq.persistence), "n" => pq.n),
        "fig4_grid" => merge(grid, Dict("log10_βA_edges" => βedges, "log10_μZ_edges" => μedges)))
    for r in out
        push!(rows, join(vcat([name], [@sprintf("%.5f", x) for x in values(r.u)],
                              [@sprintf("%.6g", r.Qo), string(r.outcome), @sprintf("%.2f", r.t_end)],
                              ), ","))
    end
    @info "scenario done" name secs extinction = frac(:extinction) agree
end
summary["wall_seconds"] = round(time() - t_start; digits = 1)

open(joinpath(outdir, "reproduction.toml"), "w") do io
    println(io, "# Louca et al. (2014) Monte Carlo reproduction summary, written by")
    println(io, "# scripts/louca2014_montecarlo.jl. Our simulations, not published data. Do not edit by hand.")
    println(io, "# Fractions are of n runs; risk effects as in louca_risk_effects; bins of -1 are empty.")
    TOML.print(io, summary; sorted = true)
end
# CSV: the columns of u differ between scenarios (Zo, TT may be absent), so write all names.
allnames = (:Zo, :μZ, :μSA, :μST, :λA, :λT, :βA, :βT, :φ, :b, :ηA, :ηT, :TA, :TT, :Tm, :Tbr, :r,
            :τm, :τbr, :σA, :σT, :KA, :KT)
open(joinpath(outdir, "montecarlo.csv"), "w") do io
    # rows hold only the sampled u of their scenario; rewrite them against allnames
    println(io, join(vcat(["scenario"], ["u_" * string(s) for s in allnames], ["Qo", "outcome", "t_end"]), ","))
    for line in rows
        f = split(line, ",")
        name = f[1]
        sc = only(filter(s -> s[1] == name, scenarios))
        keysu = keys(louca_ranges(; tadpoles = sc[2], reservoir = sc[3]))
        vals = Dict(zip(keysu, f[2:(1 + length(keysu))]))
        rest = f[(2 + length(keysu)):end]
        println(io, join(vcat([name], [get(vals, s, "") for s in allnames], rest), ","))
    end
end
println("wrote ", outdir)
