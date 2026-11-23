# Outcome table for decision analysis: the Drawert et al. (2017) individual-based model over
# draws × actions × replicates (× model variants), one CSV row per run.
#
#   julia --threads=8 --project scripts/drawert2017_decision_table.jl [options]
#
#   --draws=100 | all        draws (default 100), stratified by the authors' survival group
#   --actions=ten | paper    none + 3 kinds × days 14, 28, 56 at 75% (default), or none + the
#                            paper's 27 strategies
#   --variants=baseline      comma-separated: baseline, nu, sigma (S* = 1.2)
#   --threshold=0            a frog is infected if its load (sporangia) is above this
#   --nrep=30                replicates per cell
#   --out=DIR                default data/drawert2017
#   --work=DIR               finished (variant, draw) blocks, for resuming (default: a folder
#                            under tempdir() named after the options)
#   --time=50                only time: this many random points, every action, 1 replicate each
#
# Seeds: drawert_run_seed(draw, replicate), the same for every action and variant (common random
# numbers; one stream per run, no per-frog streams). With the same options, Julia version and
# package versions the output is identical, whatever the number of threads. Wall times are in a
# separate file, decision_runs_seconds.csv (same row order), because they are not reproducible.
#
# Writes decision_runs.csv, decision_runs_seconds.csv, decision_draws.csv and
# decision_design.toml. Every run is kept; a failed run is a row with status "error:…".
using Chytridiomycosis, Random, Statistics, TOML, Printf

const C = Chytridiomycosis
opts = Dict("draws" => "100", "actions" => "ten", "variants" => "baseline", "nrep" => "30", "threshold" => "0",
            "out" => dirname(drawert2017_path()), "work" => "", "time" => "")
for a in ARGS
    m = match(r"^--(\w+)=(.*)$", a)
    (m !== nothing && haskey(opts, m[1])) || error("unknown option $a")
    opts[m[1]] = m[2]
end
ndraws = opts["draws"] == "all" ? length(drawert2017_draws()) : parse(Int, opts["draws"])
nrep = parse(Int, opts["nrep"])
threshold = parse(Int, opts["threshold"])
variants = Symbol.(split(opts["variants"], ','))
all(in(C.DRAWERT_VARIANTS), variants) || error("variants must be among $(C.DRAWERT_VARIANTS)")
S_star = parameter_set(:drawert2017).S_star_low
ten = [DrawertIntervention(k, d, 0.75) for k in C.DRAWERT_INTERVENTION_KINDS for d in (14, 28, 56)]
opts["actions"] in ("ten", "paper") || error("actions must be ten or paper")
actions = Any[nothing; opts["actions"] == "ten" ? ten : drawert_interventions()]
draws = C._drawert_draw_sample(ndraws)

if !isempty(opts["time"])
    n = parse(Int, opts["time"])
    pts = shuffle(Xoshiro(1), drawert2017_draws())[1:n]
    for v in variants
        secs = zeros(n)
        C._drawert_decision_run(pts[1], v, nothing, 1; S_star)        # compile
        Threads.@threads :dynamic for i in 1:n
            secs[i] = sum(C._drawert_decision_run(pts[i], v, a, 1; S_star)[2] for a in actions)
        end
        @printf("%s: %d random points, %d actions, 1 replicate: %.2f s per point (median %.2f, max %.2f), %.3f s per run\n",
                v, n, length(actions), mean(secs), median(secs), maximum(secs), mean(secs) / length(actions))
        for g in (:A, :B, :C, :D, :E)
            k = findall(p -> p.group === g, pts)
            isempty(k) || @printf("  group %s: %d points, %.3f s per run\n", g, length(k),
                                  mean(secs[k]) / length(actions))
        end
        flush(stdout)
    end
    exit()
end

work = isempty(opts["work"]) ? joinpath(tempdir(), "drawert2017_decision_$(VERSION)_$(ndraws)_$(opts["actions"])_$(nrep)_t$(threshold)_v2") : opts["work"]
mkpath(work); mkpath(opts["out"])
header = join(C._DRAWERT_DECISION_COLUMNS, ',')
block(v, d) = joinpath(work, "$(v)_$(d.draw_id).csv")
jobs = [(v, d) for v in variants for d in draws if !isfile(block(v, d))]
sort!(jobs; by = j -> -j[2].simulation_time)                          # slowest first
println(length(jobs), " blocks of ", length(actions) * nrep, " runs on ", Threads.nthreads(), " threads")
done = Threads.Atomic{Int}(0)
wall = @elapsed Threads.@threads :greedy for (v, d) in jobs
    io = IOBuffer()
    for a in actions, r in 1:nrep
        row, s = C._drawert_decision_run(d, v, a, r; S_star, infection_threshold = threshold)
        println(io, C._drawert_csv_line(row), ';', @sprintf("%.5f", s))
    end
    tmp = block(v, d) * ".tmp"
    write(tmp, take!(io)); mv(tmp, block(v, d); force = true)
    k = Threads.atomic_add!(done, 1) + 1
    k % 10 == 0 && (println(k, "/", length(jobs), " blocks"); flush(stdout))
end

total = Dict{String,Float64}()
nrows = 0
open(joinpath(opts["out"], "decision_runs.csv"), "w") do f
    open(joinpath(opts["out"], "decision_runs_seconds.csv"), "w") do g
        println(f, header); println(g, "seconds")
        for v in variants, d in draws, l in eachline(block(v, d))
            a, s = split(l, ';')
            println(f, a); println(g, s)
            global nrows += 1
            total[String(v)] = get(total, String(v), 0.0) + parse(Float64, s)
        end
    end
end
nrows == length(variants) * length(draws) * length(actions) * nrep || error("missing cells")
open(joinpath(opts["out"], "decision_draws.csv"), "w") do f
    cols = (:draw_id, :group, :N, :γ, :η, :ν0, :f, :σ0, :μ, :load_growth_rate, :weight)
    println(f, join(cols, ','))
    for d in draws
        println(f, join((getfield(d, c) for c in cols), ','))
    end
end
open(joinpath(opts["out"], "decision_design.toml"), "w") do f
    println(f, "# Written by scripts/drawert2017_decision_table.jl")
    TOML.print(f, Dict("julia" => string(VERSION), "draws" => length(draws),
        "actions" => C._drawert_action_label.(actions), "variants" => String.(variants),
        "S_star" => S_star, "infection_threshold" => threshold,
        "evaluation_days" => collect(C._DRAWERT_EVALUATION_DAYS), "replicates" => nrep, "rows" => nrows); sorted = true)
end
@printf("%d rows; CPU seconds by variant %s; wall %.0f s\n", nrows, total, wall)
