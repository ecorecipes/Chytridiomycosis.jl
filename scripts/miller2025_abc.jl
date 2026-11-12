# Reproduction of the ABC of Miller et al. (2025) (§2.3, Table 3; Supplement S2, Table S3)
# and the identifiability checks of §3.1. It writes data/miller2025/abc_reproduction.toml,
# which test/miller2025.jl compares with the published numbers. About 10 minutes with 4
# threads.
#
#   julia --threads=4 --project scripts/miller2025_abc.jl [--write]
#
# Our runs are reproductions, not published results. Each ABC run uses 10^6 prior draws and
# accepts 0.2% (2,000), as the paper does.
#
# 1. Experimental data, paper priors (Table 3): 20 replicate runs (seeds 1-20). The spread
#    of the posterior means across seeds is the Monte Carlo error of a single run.
# 2. The same with the priors in the Zenodo code (β ~ U(0,1), ω ~ U(0,10)).
# 3. Synthetic data (Table S3), with μ = 0.029 (code) and μ = 0.033 (paper §S2.2).
# 4. Identifiability:
#    (a) posterior / prior SD ratios and the correlation of log α with log β;
#    (b) simulation-estimation: 20 data sets simulated at known values, refitted, and the
#        bias and 95% coverage of β_sh, β_un, α, αβ_sh, αβ_un and ω recorded;
#    (c) μ free (μ ~ U(0, 1), as in the coded generate_priors): the posterior μ–ω
#        correlation and how much the posterior of μ narrows.
using Chytridiomycosis, Random, Statistics, Printf, TOML

const C = Chytridiomycosis
write_toml = "--write" in ARGS
nrep = 20
exp_data = load_miller2025(:exp)
syn_data = load_miller2025(:syn)
q = parameter_set(:miller2025)
keys6 = (:β_sh, :β_un, :α, :αβ_sh, :αβ_un, :ω)
effects = (:refugia_effect, :vaccination_effect_sh, :vaccination_effect_un)

function replicate_runs(data; prior, μ, seeds, nsamples = 10^6)
    return [refugia_abc_summary(refugia_abc(data; prior, μ, nsamples, seed = s)) for s in seeds]
end

function describe(label, runs, published)
    println("\n== ", label, " (", length(runs), " runs)")
    out = Dict{String,Any}()
    for k in (keys6..., effects...)
        vals = k in effects ? [r[k] for r in runs] : [r[k].mean for r in runs]
        pub = get(published, k, NaN)
        @printf("  %-22s mean over runs %.4f  sd %.4f  range [%.4f, %.4f]  published %.3f\n",
                k, mean(vals), std(vals), minimum(vals), maximum(vals), pub)
        out[string(k)] = Dict("mean" => mean(vals), "sd" => std(vals), "min" => minimum(vals),
                              "max" => maximum(vals))
    end
    for k in keys6
        r1 = runs[1][k]
        @printf("  run 1 %-6s median %.4f  95%% interval (%.4f, %.4f)\n", k, r1.median, r1.lo, r1.hi)
        out[string(k)]["median_run1"] = r1.median
        out[string(k)]["lo_run1"] = r1.lo
        out[string(k)]["hi_run1"] = r1.hi
    end
    return out
end

results = Dict{String,Any}()
t0 = time()

# 1-2. Experimental data
pub = (β_sh = q.β_sh, β_un = q.β_un, α = q.α, αβ_sh = q.αβ_sh, αβ_un = q.αβ_un, ω = q.ω,
       refugia_effect = q.refugia_effect, vaccination_effect_sh = q.vaccination_effect_sh,
       vaccination_effect_un = q.vaccination_effect_un)
results["exp_paper_priors"] = describe("experimental data, paper priors (Table 3)",
    replicate_runs(exp_data; prior = RefugiaPrior(:miller2025), μ = q.μ, seeds = 1:nrep), pub)
results["exp_code_priors"] = describe("experimental data, Zenodo code priors",
    replicate_runs(exp_data; prior = RefugiaPrior(:miller2025; set = "code_priors"), μ = q.μ,
                   seeds = 1:nrep), pub)
# the same at 10^5 draws (the size used in the validation targets): Monte Carlo spread
small = replicate_runs(exp_data; prior = RefugiaPrior(:miller2025), μ = q.μ, seeds = 1:nrep,
                       nsamples = 10^5)
results["exp_paper_priors_1e5"] = describe("experimental data, paper priors, 10^5 draws",
                                           small, pub)

# 3. Synthetic data
qs = parameter_set(:miller2025; set = "synthetic_posterior")
pubs = (β_sh = qs.β_sh, β_un = qs.β_un, α = qs.α, αβ_sh = qs.αβ_sh, αβ_un = qs.αβ_un, ω = qs.ω)
μ_code = parameter_set(:miller2025; set = "synthetic").μ
results["syn_mu_code"] = describe("synthetic data, μ = $(μ_code) (code)",
    replicate_runs(syn_data; prior = RefugiaPrior(:miller2025), μ = μ_code, seeds = 1:5), pubs)
results["syn_mu_paper"] = describe("synthetic data, μ = $(qs.μ) (paper §S2.2)",
    replicate_runs(syn_data; prior = RefugiaPrior(:miller2025), μ = qs.μ, seeds = 1:5), pubs)
@printf("\nloss rate by Eq. (3): experimental %.5f, synthetic %.5f\n", loss_rate(exp_data),
        loss_rate(syn_data))

# 4a. Posterior / prior spread and the α–β ridge (experimental data, one run)
function ridge(res, prior)
    A = res.accepted
    lβ, lα = log.(A[:, 1]), log.(A[:, 3])
    pr = refugia_abc(exp_data; prior, nsamples = 10^5, acceptance = 1.0, seed = 999).accepted
    ratio(post, pri) = std(post) / std(pri)
    return (; sd_ratio_β_sh = ratio(A[:, 1], pr[:, 1]), sd_ratio_α = ratio(A[:, 3], pr[:, 3]),
            sd_ratio_αβ_sh = ratio(A[:, 1] .* A[:, 3], pr[:, 1] .* pr[:, 3]),
            sd_ratio_ω = ratio(A[:, 4], pr[:, 4]),
            cor_logα_logβ_sh = cor(lα, lβ), cor_logα_logβ_un = cor(lα, log.(A[:, 2])),
            cor_logβ_sh_logβ_un = cor(lβ, log.(A[:, 2])))
end
r1 = refugia_abc(exp_data; prior = RefugiaPrior(:miller2025), μ = q.μ, seed = 1)
rg = merge(ridge(r1, RefugiaPrior(:miller2025)),
           # how strongly the posterior itself separates the environments
           (; prob_β_un_lt_β_sh = mean(r1.accepted[:, 2] .< r1.accepted[:, 1]),
            median_ratio_β_un_β_sh = median(r1.accepted[:, 2] ./ r1.accepted[:, 1])))
println("\n== α–β ridge, experimental data (run 1)")
foreach(k -> @printf("  %-22s %.3f\n", k, rg[k]), keys(rg))
results["ridge_exp"] = Dict(string(k) => v for (k, v) in pairs(rg))

# 4b. Simulation-estimation at known values
function simulate_dataset(p, rng; nper = 4)
    w = exp_data.weeks
    counts = zeros(Int, 8, length(w), 2nper)
    for m in 1:2nper
        env = m <= nper ? :shaded : :unshaded
        counts[:, :, m] = solve(RefugiaCTMCProblem(p; environment = env, saveat = w); rng).u
    end
    return RefugiaMesocosmData(copy(w), collect(1:2nper), [m <= nper for m in 1:2nper], counts)
end
prior_width = let pr = refugia_abc(exp_data; nsamples = 10^6, acceptance = 1.0, seed = 998).accepted
    w(v) = quantile(v, 0.975) - quantile(v, 0.025)
    (β_sh = w(pr[:, 1]), β_un = w(pr[:, 2]), α = w(pr[:, 3]), αβ_sh = w(pr[:, 1] .* pr[:, 3]),
     αβ_un = w(pr[:, 2] .* pr[:, 3]), ω = w(pr[:, 4]))
end
function simest(truth; nsets = 20, label)
    p = RefugiaParameters(:miller2025; truth...)
    tv = (β_sh = p.β_sh, β_un = p.β_un, α = p.α, αβ_sh = p.α * p.β_sh, αβ_un = p.α * p.β_un,
          ω = p.ω)
    rng = Xoshiro(4242)
    sums = [refugia_abc_summary(refugia_abc(simulate_dataset(p, rng); prior = RefugiaPrior(:miller2025),
                                            μ = p.μ, seed = 100 + i)) for i in 1:nsets]
    println("\n== simulation-estimation, ", label, ", truth ", tv, ", ", nsets, " data sets")
    out = Dict{String,Any}()
    for k in keys6
        m = [s[k].mean for s in sums]
        cover = count(s -> s[k].lo <= tv[k] <= s[k].hi, sums) / nsets
        width = mean(s[k].hi - s[k].lo for s in sums)
        relerr = mean(abs.(m .- tv[k])) / tv[k]
        @printf("  %-6s truth %.4f  mean of posterior means %.4f  mean |rel. error| %.2f  95%% coverage %.2f  mean CrI width %.4f\n",
                k, tv[k], mean(m), relerr, cover, width)
        # width of the central 95% prior interval of the same quantity
        pw = prior_width[k]
        @printf("         CrI width / prior 95%% width %.3f\n", width / pw)
        out[string(k)] = Dict("truth" => tv[k], "mean_posterior_mean" => mean(m),
                              "cri_over_prior_width" => width / pw,
                              "mean_abs_rel_error" => relerr, "coverage95" => cover,
                              "mean_cri_width" => width)
    end
    return out
end
results["simest_synthetic_truth"] = simest((β_sh = 0.05, β_un = 0.02, α = 0.177, ω = 0.258,
                                             μ = 0.029); label = "Table S2 values")
results["simest_table3_truth"] = simest((β_sh = q.β_sh, β_un = q.β_un, α = q.α, ω = q.ω, μ = q.μ);
                                        label = "Table 3 posterior means")

# 4c. μ free: μ–ω confounding
function mu_free(data, label)
    r = refugia_abc(data; prior = RefugiaPrior(:miller2025; μ_max = 1.0), seed = 7)
    A = r.accepted
    c = cor(A[:, 4], A[:, 5])
    s = refugia_abc_summary(r)
    @printf("\n== μ free (U(0,1)), %s: μ mean %.4f 95%% (%.4f, %.4f); ω mean %.3f 95%% (%.3f, %.3f); cor(ω, μ) %.3f; sd(μ)/sd(prior) %.3f\n",
            label, s.μ.mean, s.μ.lo, s.μ.hi, s.ω.mean, s.ω.lo, s.ω.hi, c, std(A[:, 5]) / (1 / sqrt(12)))
    return Dict("mu_mean" => s.μ.mean, "mu_lo" => s.μ.lo, "mu_hi" => s.μ.hi,
                "omega_mean" => s.ω.mean, "omega_lo" => s.ω.lo, "omega_hi" => s.ω.hi,
                "cor_omega_mu" => c, "sd_ratio_mu" => std(A[:, 5]) / (1 / sqrt(12)))
end
results["mu_free_exp"] = mu_free(exp_data, "experimental data")
results["mu_free_syn"] = mu_free(syn_data, "synthetic data (truth μ = 0.029, ω = 0.258)")

@printf("\nTotal time %.1f min on %d threads\n", (time() - t0) / 60, Threads.nthreads())

if write_toml
    path = joinpath(pkgdir(C), "data", "miller2025", "abc_reproduction.toml")
    open(path, "w") do io
        println(io, "# OUR reproduction of the Miller et al. (2025) ABC, written by scripts/miller2025_abc.jl.")
        println(io, "# Not published results. Each run: 10^6 prior draws (10^5 where named), 0.2% accepted,")
        println(io, "# summary statistic of the authors' code; mean/sd/min/max are over replicate seeds.")
        println(io, "# Generated $(round(Int, time())) (Unix time) with $(Threads.nthreads()) threads, Julia $(VERSION).")
        TOML.print(io, results; sorted = true)
    end
    println("wrote ", path)
end
