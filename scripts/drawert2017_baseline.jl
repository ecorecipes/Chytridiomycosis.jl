# Drawert et al. (2017): the stochastic Baseline model on a stratified random sample of the
# authors' 2132 points, 30 runs per point, compared point by point with the authors' own
# 30-run summaries (data/drawert2017/baseline_points.csv, from their C code).
#   julia --threads=4 --project scripts/drawert2017_baseline.jl
# Appends to data/drawert2017/baseline_sample.csv (resumable: finished points are skipped).
# Summaries are written by scripts/drawert2017_summarise.jl.
#
# The sample is stratified by the authors' group (Fig. 3; their 30-run mean and sd of the
# surviving fraction), with the sizes BASELINE_SIZES (drawert2017_common.jl), so that each
# group is compared. Each sampled point stands for N_g/n_g points of its group in weighted
# estimates of population fractions.
# Cost: the authors' table records 1.35e6 CPU seconds for the 2132 points; this simulator is
# about 7 times faster per run, and the sample below costs about 3 CPU-hours.
include("drawert2017_common.jl")

pts = load_drawert2017()
ids, weights, group_of = stratified_sample(pts, BASELINE_SIZES)
path = joinpath(OUTDIR, "baseline_sample.csv")
done = Set(parse(Int, r["id"]) for r in read_rows(path))
todo = filter(!in(done), ids)
println("sample of ", length(ids), " points; ", length(todo), " to run on ",
        Threads.nthreads(), " threads")
# expensive points first, so the threads finish together
sort!(todo; by = id -> -pts[id].simulation_time)

header = ["id", "authors_group", "weight", "N", "mean_surviving", "std_surviving",
          "frac_frogs_extinct", "frac_fungus_extinct", "frac_persistence",
          "mean_days_to_first_death", "std_days_to_first_death", "outcome_class", "group",
          "seconds"]
lk = ReentrantLock()
Threads.@threads :dynamic for id in todo
    p = DrawertParameters(:drawert2017; point = id)
    t = @elapsed s = drawert_summary(DrawertIBMProblem(p); rng = point_rng(id))
    g = group_of[id]
    append_row(path, header, fmt.([id, g, weights[g], p.N, s.mean_surviving, s.std_surviving,
        s.frac_frogs_extinct, s.frac_fungus_extinct, s.frac_persistence,
        s.mean_days_to_first_death, s.std_days_to_first_death, s.outcome_class, s.group, t]), lk)
    @printf("point %4d (group %s): %.1f s\n", id, g, t)
    flush(stdout)
end
