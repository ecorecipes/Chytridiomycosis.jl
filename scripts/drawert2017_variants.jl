# Drawert et al. (2017): the stochastic Nu and Sigma models (ESM §1.1.2–1.1.3), whose load
# dependence is set so that loads level off at S* × Smax (S* = 1.2 or 2; ESM §1.4).
#   julia --threads=4 --project scripts/drawert2017_variants.jl
# Appends to data/drawert2017/variants_sample.csv (resumable).
#
# The Nu and Sigma runs cost about five times a Baseline run (loads linger near Smax), so
# this is a small study with 10 runs per point and variant:
# - a simple random sample of 12 of the authors' points (all groups), S* = 1.2 and 2;
# - the Group C and D points of the Baseline sample (scripts/drawert2017_baseline.jl), the
#   cheap ones, S* = 1.2.
# Outcomes are compared with the Baseline at the same points (the authors' 30 runs and ours).
include("drawert2017_common.jl")

const NRUNS = 10
pts = load_drawert2017()
q = parameter_set(:drawert2017)
rng = Xoshiro(SEED + 1)
random_ids = sort!(shuffle(rng, [r.id for r in pts])[1:12])
sample_ids, _, group_of = stratified_sample(pts, BASELINE_SIZES)
cd_ids = filter(id -> group_of[id] in (:C, :D), sample_ids)

jobs = Tuple{Int,Symbol,Float64,String}[]
for id in random_ids, v in (:nu, :sigma), S in (q.S_star_low, q.S_star_high)
    push!(jobs, (id, v, S, "random"))
end
for id in cd_ids, v in (:nu, :sigma)
    id in random_ids && continue
    push!(jobs, (id, v, q.S_star_low, "groupCD"))
end
path = joinpath(OUTDIR, "variants_sample.csv")
done = Set((parse(Int, r["id"]), Symbol(r["variant"]), parse(Float64, r["S_star"]))
           for r in read_rows(path))
filter!(j -> !((j[1], j[2], j[3]) in done), jobs)
sort!(jobs; by = j -> -pts[j[1]].simulation_time)
println(length(jobs), " variant ensembles to run on ", Threads.nthreads(), " threads")

header = ["id", "authors_group", "sample", "variant", "S_star", "nruns", "mean_surviving",
          "std_surviving", "frac_frogs_extinct", "frac_fungus_extinct", "frac_persistence",
          "outcome_class", "seconds"]
lk = ReentrantLock()
arm(v, S) = (v === :nu ? 1 : 3) + (S == q.S_star_low ? 0 : 1)
Threads.@threads :dynamic for (id, v, S, smp) in jobs
    p = drawert_variant(DrawertParameters(:drawert2017; point = id), v; S_star = S)
    t = @elapsed s = drawert_summary(DrawertIBMProblem(p); nruns = NRUNS,
                                     rng = point_rng(id, arm(v, S)))
    append_row(path, header, fmt.([id, group_of[id], smp, v, S, NRUNS, s.mean_surviving,
        s.std_surviving, s.frac_frogs_extinct, s.frac_fungus_extinct, s.frac_persistence,
        s.outcome_class, t]), lk)
    @printf("point %4d %s S*=%.1f: %.1f s\n", id, v, S, t)
    flush(stdout)
end
