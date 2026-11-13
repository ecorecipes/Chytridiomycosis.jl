# Drawert et al. (2017): summarise the reproduction runs into
# data/drawert2017/reproduction.toml, which test/drawert2017.jl and the validation targets read.
#   julia --project scripts/drawert2017_summarise.jl
# Inputs (data/drawert2017/): deterministic.csv (scripts/drawert2017_deterministic.jl; not
# committed, 380 kB, ~3 minutes to regenerate), baseline_sample.csv, variants_sample.csv and
# interventions_sample.csv (the other drawert2017_*.jl scripts).
include("drawert2017_common.jl")
using TOML

pts = load_drawert2017()
byid = Dict(r.id => r for r in pts)
q = parameter_set(:drawert2017)
out = Dict{String,Any}()
num(rows, k) = [parse(Float64, r[k]) for r in rows]
r4(x) = round(x; digits = 4)

# --- deterministic checks on all points --------------------------------------------------------
det = read_rows(joinpath(OUTDIR, "deterministic.csv"))
if !isempty(det)
    λ, λa = num(det, "lambda"), num(det, "lambda_authors")
    fd, ld = num(det, "det_first_death"), num(det, "det_last_death")
    out["deterministic"] = Dict(
        "n_points" => length(det),
        "max_abs_lambda_difference" => maximum(abs.(λ .- λa)),
        "frogs_extinct_baseline" => count(r -> r["det_outcome"] == "frogs_extinct", det),
        "frogs_extinct_nu_1.2" => count(r -> r["nu12_outcome"] == "frogs_extinct", det),
        "frogs_extinct_nu_2" => count(r -> r["nu2_outcome"] == "frogs_extinct", det),
        "frogs_extinct_sigma_1.2" => count(r -> r["sigma12_outcome"] == "frogs_extinct", det),
        "frogs_extinct_sigma_2" => count(r -> r["sigma2_outcome"] == "frogs_extinct", det),
        "first_death_day_min" => minimum(ceil.(fd)), "first_death_day_max" => maximum(ceil.(fd)),
        "last_death_min" => r4(minimum(ld)), "last_death_max" => r4(maximum(ld)),
        "last_death_outside_window" => count(t -> !(q.t_death_min <= t <= q.t_death_max), ld),
        "printed_nu_Sstar_below_1_at_1.2" => count(<(1), num(det, "nu12_printed_Sstar")),
        "printed_nu_Sstar_below_1_at_2" => count(<(1), num(det, "nu2_printed_Sstar")),
        "nearest_neighbour_min" => r4(minimum(num(det, "nn"))),
        "nearest_neighbour_median" => r4(median(num(det, "nn"))),
        "nearest_neighbour_printed_min" => r4(minimum(num(det, "nn_printed"))),
        "nearest_neighbour_printed_median" => r4(median(num(det, "nn_printed"))))
end

# --- Baseline sample vs the authors' runs ------------------------------------------------------
bs = read_rows(joinpath(OUTDIR, "baseline_sample.csv"))
classes = ("frog_extinction", "fungal_extinction", "coexistence", "variable")
if !isempty(bs)
    w = num(bs, "weight")
    W = sum(w)
    theirs = [byid[parse(Int, r["id"])] for r in bs]
    their_class = [string(outcome_class(t.frac_frogs_extinct, t.frac_fungus_extinct,
                                        t.frac_persistence)) for t in theirs]
    their_group = [string(survival_group(t.mean_surviving, t.std_surviving)) for t in theirs]
    ours_class = [r["outcome_class"] for r in bs]
    ours_group = [r["group"] for r in bs]
    wfrac(v, x) = r4(sum(w[i] for i in eachindex(v) if v[i] == x; init = 0.0) / W)
    # standard error of a stratified proportion, with the finite-population correction
    strata = [r["authors_group"] for r in bs]
    function wse(v, x)
        s2 = 0.0
        for g in unique(strata)
            ix = findall(==(g), strata)
            n = length(ix)
            Ng = n * w[ix[1]]
            ph = count(i -> v[i] == x, ix) / n
            n > 1 && (s2 += (Ng / W)^2 * (1 - n / Ng) * ph * (1 - ph) / (n - 1))
        end
        return r4(sqrt(s2))
    end
    pf = Float64[]
    pm = Float64[]
    for (r, t) in zip(bs, theirs)
        a = round(Int, 30 * parse(Float64, r["frac_frogs_extinct"]))
        c = round(Int, 30 * t.frac_frogs_extinct)
        push!(pf, fisher_p(a, 30 - a, c, 30 - c))
        push!(pm, welch_p(parse(Float64, r["mean_surviving"]), parse(Float64, r["std_surviving"]),
                          30, t.mean_surviving, t.std_surviving, 30))
    end
    dm = num(bs, "mean_surviving") .- [t.mean_surviving for t in theirs]
    dfd = num(bs, "mean_days_to_first_death") .- [t.mean_days_to_first_death for t in theirs]
    groups = Dict{String,Any}()
    for g in ("A", "B", "C", "D", "E")
        ix = findall(==(g), [r["authors_group"] for r in bs])
        isempty(ix) && continue
        groups[g] = Dict("n" => length(ix),
            "same_group" => count(i -> ours_group[i] == g, ix),
            "mean_surviving_ours" => r4(mean(num(bs[ix], "mean_surviving"))),
            "mean_surviving_authors" => r4(mean(theirs[i].mean_surviving for i in ix)),
            "fisher_p_below_0.05" => count(<(0.05), pf[ix]),
            "welch_p_below_0.05" => count(<(0.05), pm[ix]))
    end
    out["baseline_sample"] = Dict(
        "n_points" => length(bs), "runs_per_point" => 30, "seed" => SEED,
        "class_weighted_ours" => Dict(c => wfrac(ours_class, c) for c in classes),
        "class_weighted_ours_se" => Dict(c => wse(ours_class, c) for c in classes),
        "group_weighted_ours_se" => Dict(g => wse(ours_group, g) for g in ("A", "B", "C", "D", "E")),
        "class_weighted_authors_same_points" => Dict(c => wfrac(their_class, c) for c in classes),
        "group_weighted_ours" => Dict(g => wfrac(ours_group, g) for g in ("A", "B", "C", "D", "E")),
        "same_class" => count(ours_class .== their_class),
        "same_group" => count(ours_group .== their_group),
        "fisher_frogs_extinct_p_below_0.05" => count(<(0.05), pf),
        "welch_mean_surviving_p_below_0.05" => count(<(0.05), pm),
        "mean_surviving_difference_mean" => r4(mean(dm)),
        "mean_surviving_difference_mean_abs" => r4(mean(abs.(dm))),
        "mean_surviving_correlation" => r4(cor(num(bs, "mean_surviving"),
                                               [t.mean_surviving for t in theirs])),
        "days_to_first_death_difference_mean" => r4(mean(dfd)),
        "days_to_first_death_difference_median" => r4(median(dfd)),
        "by_authors_group" => groups,
        "cpu_seconds" => round(sum(num(bs, "seconds")); digits = 1),
        "authors_cpu_seconds_same_points" => round(sum(t.simulation_time for t in theirs); digits = 1))
end

# --- Nu and Sigma models -----------------------------------------------------------------------
vs = read_rows(joinpath(OUTDIR, "variants_sample.csv"))
if !isempty(vs)
    v = Dict{String,Any}()
    for smp in ("random", "groupCD"), var in ("nu", "sigma"), S in ("1.2", "2")
        rows = filter(r -> r["sample"] == smp && r["variant"] == var && r["S_star"] == S, vs)
        isempty(rows) && continue
        ids = [parse(Int, r["id"]) for r in rows]
        v["$(smp)_$(var)_$(S)"] = Dict("n_points" => length(rows),
            "runs_frogs_extinct" => r4(mean(num(rows, "frac_frogs_extinct"))),
            "runs_fungus_extinct" => r4(mean(num(rows, "frac_fungus_extinct"))),
            "runs_persistence" => r4(mean(num(rows, "frac_persistence"))),
            "points_with_persistence" => count(>(0), num(rows, "frac_persistence")),
            "mean_surviving" => r4(mean(num(rows, "mean_surviving"))),
            "baseline_authors_runs_frogs_extinct" => r4(mean(byid[i].frac_frogs_extinct for i in ids)),
            "baseline_authors_runs_persistence" => r4(mean(byid[i].frac_persistence for i in ids)),
            "baseline_authors_mean_surviving" => r4(mean(byid[i].mean_surviving for i in ids)))
    end
    v["runs_per_point"] = parse(Int, vs[1]["nruns"])
    out["variants"] = v
end

# --- interventions -----------------------------------------------------------------------------
iv = read_rows(joinpath(OUTDIR, "interventions_sample.csv"))
if !isempty(iv)
    base = Dict(parse(Int, r["id"]) => parse(Float64, r["mean_surviving"]) for r in bs)
    for r in iv
        r["strategy"] == "0" && (base[parse(Int, r["id"])] = parse(Float64, r["mean_surviving"]))
    end
    Δ(r) = parse(Float64, r["mean_surviving"]) - base[parse(Int, r["id"])]
    filter!(r -> haskey(base, parse(Int, r["id"])), iv)      # skip points still running
    res = Dict{String,Any}()
    # grid at the representative points
    for id in unique(parse(Int, r["id"]) for r in iv if r["set"] == "grid")
        rows = filter(r -> r["set"] == "grid" && parse(Int, r["id"]) == id && r["strategy"] != "0", iv)
        g = Dict{String,Any}("baseline_mean_surviving" => r4(base[id]))
        for kind in ("clean_frogs", "cull_frogs", "clean_environment"), coded in (false, true)
            rk = filter(r -> r["kind"] == kind && parsebool(r["as_coded"]) == coded, rows)
            isempty(rk) && continue
            key = coded ? kind * "_as_coded" : kind
            g[key] = Dict("delta" => r4.(Δ.(rk)),
                          "day" => num(rk, "day"), "fraction" => num(rk, "fraction"),
                          "positive" => count(r -> Δ(r) >= 0.1 - 1e-12, rk),
                          "negative" => count(r -> Δ(r) <= -0.1 + 1e-12, rk))
        end
        res["grid_$id"] = g
    end
    # random strategies at sampled points, by authors' group
    smp = filter(r -> r["set"] == "sample" && r["strategy"] != "0", iv)
    for grp in ("A", "B", "C", "D"), kind in ("clean_frogs", "cull_frogs", "clean_environment",
                                              "clean_environment_as_coded", "placebo")
        coded = endswith(kind, "_as_coded")
        k = coded ? "clean_environment" : kind
        rk = filter(r -> r["authors_group"] == grp && r["kind"] == k &&
                         parsebool(r["as_coded"]) == coded, smp)
        isempty(rk) && continue
        d = Δ.(rk)
        res["sample_$(grp)_$(kind)"] = Dict("n" => length(rk),
            "positive" => count(>=(0.1 - 1e-12), d), "negative" => count(<=(-0.1 + 1e-12), d),
            "mean_delta" => r4(mean(d)))
    end
    out["interventions"] = res
end

path = joinpath(OUTDIR, "reproduction.toml")
open(path, "w") do io
    println(io, "# Drawert et al. (2017) reproduction summary, written by scripts/drawert2017_summarise.jl.")
    println(io, "# Do not edit by hand. Seeds and designs are in scripts/drawert2017_*.jl.")
    TOML.print(io, out; sorted = true)
end
println(read(path, String))
