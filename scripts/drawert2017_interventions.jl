# Drawert et al. (2017): one-time interventions in the stochastic Baseline model (§3.3,
# §4.1, Fig. 4, ESM Fig. S10, S14).
#   julia --threads=4 --project scripts/drawert2017_interventions.jl
# Appends to data/drawert2017/interventions_sample.csv (resumable). Run it after
# scripts/drawert2017_baseline.jl: the no-intervention arm of a sampled point is that
# script's ensemble (same seed).
#
# As in the paper, each arm is an independent 30-run ensemble, and an intervention's effect at
# a point is Δ = (mean surviving fraction with it) − (without it). A full analysis (27
# strategies at each of 10 660 points) is out of reach here, so:
# - "grid": all 27 strategies, plus cleaning the environment as coded (all zoospores
#   removed), at the cheap representative points of Groups C (id 182) and D (id 982);
# - "sample": at Baseline-sample points of Groups B (6 points), C (20) and D (40), and at the
#   Group A representative point (id 120), one strategy of each kind with day and fraction
#   drawn at random from the design, cleaning the environment as coded with the same day and
#   fraction, and a placebo: a second no-intervention ensemble with another seed, which
#   measures how often chance alone gives |Δ| ≥ 0.1.
include("drawert2017_common.jl")

pts = load_drawert2017()
q = parameter_set(:drawert2017)
strategies = drawert_interventions()
sample_ids, _, group_of = stratified_sample(pts, BASELINE_SIZES)
pick(g, n) = filter(id -> group_of[id] === g, sample_ids)[1:n]
sample_points = vcat(pick(:B, 6), pick(:C, 20), pick(:D, 40),
                     [round(Int, parameter_set(:drawert2017; set = "groupA").point_id)])
grid_points = [round(Int, parameter_set(:drawert2017; set = s).point_id) for s in ("groupC", "groupD")]

# (id, set, strategy index (0 = none, -1 = placebo), as_coded)
jobs = Tuple{Int,String,Int,Bool}[]
for id in grid_points
    push!(jobs, (id, "grid", 0, false))
    for k in eachindex(strategies)
        push!(jobs, (id, "grid", k, false))
        strategies[k].kind === :clean_environment && push!(jobs, (id, "grid", k, true))
    end
end
for id in sample_points
    rng = Xoshiro(SEED + 2 + id)
    id in sample_ids || push!(jobs, (id, "sample", 0, false))
    for (ik, kind) in enumerate(C.DRAWERT_INTERVENTION_KINDS)
        k = 9 * (ik - 1) + rand(rng, 1:9)
        push!(jobs, (id, "sample", k, false))
        kind === :clean_environment && push!(jobs, (id, "sample", k, true))
    end
    push!(jobs, (id, "sample", -1, false))
end
path = joinpath(OUTDIR, "interventions_sample.csv")
done = Set((parse(Int, r["id"]), r["set"], parse(Int, r["strategy"]), parsebool(r["as_coded"]))
           for r in read_rows(path))
filter!(j -> !(j in done), jobs)
sort!(jobs; by = j -> -pts[j[1]].simulation_time)
println(length(jobs), " ensembles to run on ", Threads.nthreads(), " threads")

header = ["id", "authors_group", "set", "strategy", "kind", "day", "fraction", "as_coded",
          "mean_surviving", "std_surviving", "frac_frogs_extinct", "frac_fungus_extinct",
          "frac_persistence", "seconds"]
lk = ReentrantLock()
Threads.@threads :dynamic for (id, set, k, coded) in jobs
    p = DrawertParameters(:drawert2017; point = id)
    iv = k > 0 ? strategies[k] : nothing
    arm = k == -1 ? 99 : k == 0 ? 0 : 100 + k + (coded ? 50 : 0)
    prob = DrawertIBMProblem(p; intervention = iv, as_coded = coded)
    t = @elapsed s = drawert_summary(prob; rng = point_rng(id, arm))
    kind = k == -1 ? "placebo" : k == 0 ? "none" : string(iv.kind)
    append_row(path, header, fmt.([id, group_of[id], set, k, kind, k > 0 ? iv.day : 0,
        k > 0 ? iv.fraction : 0, coded, s.mean_surviving, s.std_surviving,
        s.frac_frogs_extinct, s.frac_fungus_extinct, s.frac_persistence, t]), lk)
    @printf("point %4d %-6s %-17s %s: %.1f s\n", id, set, kind, coded ? "coded" : "", t)
    flush(stdout)
end
