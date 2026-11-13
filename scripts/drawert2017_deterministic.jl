# Drawert et al. (2017): checks on all 2132 sampled points that need no stochastic runs.
#   julia --threads=4 --project scripts/drawert2017_deterministic.jl
# Writes data/drawert2017/deterministic.csv (one row per point) and prints a summary.
#
# Per point:
# - λ: our linearised growth rate vs the authors' column (ESM §1.4);
# - the deterministic Baseline (DrawertODEProblem): outcome, time to first and last death
#   (§3.2: "death of all frogs within 14–275 days");
# - the Nu and Sigma variants (ESM §1.4, S* = 1.2 and 2): deterministic outcome, and the S*
#   reached with ESM eq. (11) as printed;
# - the non-dimensional groups c1–c4 and the nearest-neighbour distance of ESM eq. (10).
include("drawert2017_common.jl")

pts = load_drawert2017()
q = parameter_set(:drawert2017)
cs = [nondimensional_parameters(DrawertParameters(:drawert2017; point = r.id)) for r in pts]
rg = nondimensional_ranges(cs)
rgp = nondimensional_ranges(cs; as_printed = true)
nn = fill(Inf, length(cs))
nnp = fill(Inf, length(cs))
Threads.@threads for i in eachindex(cs)
    for j in eachindex(cs)
        i == j && continue
        nn[i] = min(nn[i], nondimensional_distance(cs[i], cs[j], rg))
        nnp[i] = min(nnp[i], nondimensional_distance(cs[i], cs[j], rgp; as_printed = true))
    end
end

header = ["id", "lambda", "lambda_authors", "det_outcome", "det_first_death", "det_last_death",
          "nu12_outcome", "nu2_outcome", "sigma12_outcome", "sigma2_outcome",
          "nu12_printed_Sstar", "nu2_printed_Sstar", "c1", "c2", "c3", "c4", "nn", "nn_printed"]
rows = Vector{Vector{String}}(undef, length(pts))
Threads.@threads for i in eachindex(pts)
    r = pts[i]
    p = DrawertParameters(:drawert2017; point = r.id)
    s = solve(DrawertODEProblem(p))
    dt = filter(!isnan, s.death_time)
    outs = String[]
    for v in (:nu, :sigma), S in (q.S_star_low, q.S_star_high)
        pv = drawert_variant(p, v; S_star = S)
        push!(outs, string(solve(DrawertODEProblem(pv)).outcome))
    end
    pr = [equilibrium_load_ratio(drawert_variant(p, :nu; S_star = S, as_printed = true))
          for S in (q.S_star_low, q.S_star_high)]
    c = cs[i]
    rows[i] = fmt.([r.id, load_growth_rate(p), r.λ, string(s.outcome),
                    isempty(dt) ? NaN : minimum(dt), isempty(dt) ? NaN : maximum(dt),
                    outs..., pr..., c.c1, c.c2, c.c3, c.c4, nn[i], nnp[i]])
end
path = joinpath(OUTDIR, "deterministic.csv")
open(path, "w") do io
    println(io, join(header, ','))
    foreach(r -> println(io, join(r, ',')), rows)
end

d = read_rows(path)
num(k) = [parse(Float64, x[k]) for x in d]
λ, λa = num("lambda"), num("lambda_authors")
@printf("max |λ − λ_authors| = %.2e (max relative %.2e)\n", maximum(abs.(λ .- λa)),
        maximum(abs.(λ .- λa) ./ λa))
@printf("λ range %.4f–%.4f\n", minimum(λ), maximum(λ))
for k in ("det_outcome", "nu12_outcome", "nu2_outcome", "sigma12_outcome", "sigma2_outcome")
    println(k, ": ", Dict(o => count(x -> x[k] == o, d) for o in unique(x[k] for x in d)))
end
ld = num("det_last_death")
@printf("deterministic last death: %.1f–%.1f days; outside 14–275: %d\n", minimum(ld),
        maximum(ld), count(t -> !(q.t_death_min <= t <= q.t_death_max), ld))
for k in ("nu12_printed_Sstar", "nu2_printed_Sstar")
    v = num(k)
    @printf("%s: S* < 1 for %d of %d points; median S* %.3g\n", k, count(<(1), v), length(v),
            median(v))
end
@printf("nearest-neighbour distance: range of log c %.4f (median %.4f); as printed %.4f (median %.4f)\n",
        minimum(num("nn")), median(num("nn")), minimum(num("nn_printed")), median(num("nn_printed")))
