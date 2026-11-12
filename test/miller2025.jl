# Miller et al. 2025 (arXiv:2503.06846 v3): CTMC of Bd in Litoria aurea mesocosms with
# shaded or sunlight-heated refugia. The oracles are:
# - the paper (Table 1, Table 3, Table S1, eqs 3-5);
# - the authors' R code (Zenodo 10.5281/zenodo.15179344): their state-change matrix and
#   propensities;
# - the stored output of our ABC reproduction (data/miller2025/abc_reproduction.toml,
#   written by scripts/miller2025_abc.jl), compared with the published posterior.
# The distributional comparison with the authors' R simulator is scripts/miller2025_vs_R.jl.
using Chytridiomycosis
using Test
using Random: Xoshiro
using Distributions: mean, std
using TOML

@testset "Miller 2025 refugia CTMC" begin
    C = Chytridiomycosis
    q = parameter_set(:miller2025)
    p = RefugiaParameters(:miller2025)
    data = load_miller2025(:exp)

    @testset "parameters (Table S1, Table 3)" begin
        @test all(r -> r.verified, parameter_records(:miller2025))
        @test (p.γ1, p.γ2) == (1 / 2.5, q.γ2) && q.γ2 ≈ 1 / 4.5
        @test p.m == (10.0, 100.0, 1.0, 1.0, 10.0, 0.1)
        @test (p.β_sh, p.β_un, p.α, p.ω, p.μ) == (0.675, 0.401, 0.072, 0.569, 0.021)
        @test time_unit(p) === :week
        @test transmissibility(p, :shaded) == p.β_sh && transmissibility(p, :unshaded) == p.β_un
        @test_throws ArgumentError transmissibility(p, :sunny)
        @test RefugiaParameters(:miller2025; β_sh = 0.5).β_sh == 0.5
        s = RefugiaParameters(:miller2025; set = "synthetic")
        @test (s.β_sh, s.β_un, s.α, s.ω, s.μ) == (0.05, 0.02, 0.177, 0.258, 0.029)
        @test refugia_initial_state() == (5, 5, 0, 0, 5, 5, 0, 0)
        # derived records, recomputed (Eqs. 3-5)
        @test round(loss_rate(data); digits = 3) == q.μ
        @test sum(data.counts[:, end, :]) / 8 == q.N_Tf
        @test round(1 - q.β_un / q.β_sh; digits = 3) == q.refugia_effect
        @test round(1 - q.αβ_sh / q.β_sh; digits = 3) == q.vaccination_effect_sh
        @test round(1 - q.αβ_un / q.β_un; digits = 3) == q.vaccination_effect_un
        syn = load_miller2025(:syn)
        qs = parameter_set(:miller2025; set = "synthetic_posterior")
        @test round(sum(syn.counts[:, end, :]) / 8; digits = 1) == qs.N_Tf
        @test round(loss_rate(syn); digits = 3) ≈ qs.μ atol = 0.001   # 0.0335 → "0.033"
        @test RefugiaPrior(:miller2025) == RefugiaPrior(; β_max = 2.0, α_max = 0.5, ω_max = 2.0)
        @test RefugiaPrior(:miller2025; set = "code_priors") ==
              RefugiaPrior(; β_max = 1.0, α_max = 0.5, ω_max = 10.0)
    end

    @testset "data" begin
        @test data.weeks == [0, 1, 2, 4, 6, 8, 10, 15] == refugia_observation_weeks()
        @test data.mesocosm == [2, 3, 5, 6, 8, 9, 11, 12]
        @test data.mesocosm[data.shaded] == [3, 5, 9, 11]
        @test size(data.counts) == (8, 8, 8)
        @test all(==(20), sum(data.counts[:, 1, :]; dims = 1))
        @test all(data.counts[[1, 2, 5], 1, :] .== 5)     # S_U, I_U1, S_V at week 0
        # week 0 has inoculated vaccinated frogs already in I_V2 in three mesocosms
        @test data.mesocosm[vec(data.counts[7, 1, :] .> 0)] == [2, 3, 6]
        syn = load_miller2025(:syn)
        @test syn.mesocosm == 1:8 && syn.shaded == [trues(4); falses(4)]
        @test_throws ArgumentError miller2025_path(:other)
    end

    @testset "transitions (Table 1; authors' nu and propensities)" begin
        ν = refugia_stoichiometry()
        # nu of build_stochastic_system() in fns_stochastic_model.R, row by row
        ν_R = zeros(Int, 8, 16)
        for (j, from, to) in ((1, 1, 2), (2, 1, 0), (3, 2, 3), (4, 2, 0), (5, 3, 4), (6, 3, 0),
                              (7, 4, 5), (8, 4, 0), (9, 5, 6), (10, 5, 0), (11, 6, 7),
                              (12, 6, 0), (13, 7, 8), (14, 7, 0), (15, 8, 5), (16, 8, 0))
            ν_R[from, j] = -1
            to == 0 || (ν_R[to, j] = 1)
        end
        @test ν == ν_R
        @test all(sum(ν; dims = 1) .∈ Ref((0, -1)))          # moves conserve, losses remove one
        @test count(==(-1), sum(ν; dims = 1)) == 8
        x = (3, 2, 1, 1, 4, 1, 2, 1)                            # N = 15
        F = 10 * 2 + 100 * 1 + 1 * 1 + 1 * 1 + 10 * 2 + 0.1 * 1
        a = refugia_rates(p, :shaded, x)
        @test a[1] ≈ 0.675 * 3 * F / 14
        @test a[9] ≈ 0.072 * 0.675 * 4 * F / 14
        @test refugia_rates(p, :unshaded, x)[1] ≈ 0.401 * 3 * F / 14
        @test a[3] ≈ 0.4 * 2 && a[5] ≈ 1 / 4.5 && a[7] ≈ 0.569 && a[15] ≈ 0.569
        @test a[11] ≈ 0.4 && a[13] ≈ 2 / 4.5
        @test a[[2, 4, 6, 8, 10, 12, 14, 16]] ≈ 0.021 .* collect(x)
        # one frog left: the 0/0 propensity is zero, as GillespieSSA's NaN -> 0
        @test refugia_rates(p, :shaded, (1, 0, 0, 0, 0, 0, 0, 0))[1] == 0
        @test all(isfinite, refugia_rates(p, :shaded, (0, 0, 0, 0, 0, 0, 0, 0)))
        @test_throws ArgumentError refugia_rates(p, :shaded, (1, 2, 3))
    end

    @testset "simulation" begin
        prob = RefugiaCTMCProblem(p; environment = :unshaded)
        r1 = solve(prob; rng = Xoshiro(3))
        r2 = solve(prob; rng = Xoshiro(3))
        @test r1.u == r2.u && r1.events == r2.events               # seeded
        @test r1.t == data.weeks && size(r1.u) == (8, 8)
        @test r1.u[:, 1] == collect(refugia_initial_state())
        tr = solve(prob; rng = Xoshiro(3), trajectory = true)
        @test tr.u == r1.u && length(tr.jump_times) == tr.events
        @test issorted(tr.jump_times) && all(<=(15), tr.jump_times)
        # every jump is one column of the state-change matrix
        ν = refugia_stoichiometry()
        xs = vcat([collect(refugia_initial_state())], [collect(s) for s in tr.jump_states])
        @test all(any(xs[i+1] .- xs[i] == ν[:, j] for j in 1:16) for i in 1:tr.events)
        # survey columns equal the jump-path state at each week
        for (k, w) in enumerate(tr.t)
            i = searchsortedlast(tr.jump_times, w)
            @test tr.u[:, k] == xs[i+1]
        end
        totals = vec(sum(tr.u; dims = 1))
        @test all(tr.u .>= 0) && issorted(totals; rev = true) && totals[1] == 20
        # no loss: 20 frogs throughout; vaccinated frogs never return to the U states
        p0 = RefugiaParameters(:miller2025; μ = 0.0, β_sh = 1.5, ω = 2.0)
        for seed in 1:50
            u = solve(RefugiaCTMCProblem(p0); rng = Xoshiro(seed)).u
            @test all(==(20), sum(u; dims = 1))
            @test all(diff(vec(sum(u[1:4, :]; dims = 1))) .<= 0)
        end
        # all frogs lost: the run stops and later surveys see zeros
        pk = RefugiaParameters(:miller2025; μ = 50.0)
        @test all(solve(RefugiaCTMCProblem(pk); rng = Xoshiro(1)).u[:, end] .== 0)
        @test_throws ArgumentError RefugiaCTMCProblem(p; saveat = [0.0, 20.0])
        @test_throws ArgumentError RefugiaCTMCProblem(p; u0 = (1, 2))
    end

    @testset "mean-field ODE (Supplement S3)" begin
        w = data.weeks
        # with β = 0 the chain is linear, so the ODE is its exact mean
        p0 = RefugiaParameters(:miller2025; β_sh = 0.0, β_un = 0.0)
        sims = [solve(RefugiaCTMCProblem(p0); rng = Xoshiro(i)).u for i in 1:4000]
        M, S = mean(sims), std(sims) ./ sqrt(4000)
        O = reduce(hcat, solve(RefugiaODEProblem(p0); saveat = w, abstol = 1e-10, reltol = 1e-10).u)
        @test maximum(abs.(O .- M)[S .> 0] ./ S[S .> 0]) < 4
        @test O[2, end] ≈ 5 * exp(-(0.4 + 0.021) * 15) rtol = 1e-6         # I_U1 decays
        @test sum(O[:, end]) ≈ 20 * exp(-0.021 * 15) rtol = 1e-6
        # with transmission it is close to the simulation mean (not exact: small N)
        sims = [solve(RefugiaCTMCProblem(p); rng = Xoshiro(i)).u for i in 1:2000]
        O = reduce(hcat, solve(RefugiaODEProblem(p); saveat = w).u)
        @test maximum(abs.(O .- mean(sims))) < 1.0
        @test C.to_ode_problem(RefugiaODEProblem(p)).u0 == [5.0, 5, 0, 0, 5, 5, 0, 0]
    end

    @testset "ABC summary statistic (Eq. 2, as coded)" begin
        x0 = repeat(collect(refugia_initial_state()), 1, 8)
        d = refugia_distance(data, x0, x0)
        manual = sum((x0[i, w] - data.counts[i, w, m])^2 for m in 1:8, w in 1:8, i in 1:8)
        @test d.total == manual == d.shaded + d.unshaded
        week0 = sum((x0[i, 1] - data.counts[i, 1, m])^2 for m in 1:8, i in 1:8)
        @test week0 > 0                                      # data week 0 ≠ x0 (IV2 frogs)
        @test refugia_distance(data, x0, x0; include_t0 = false).total == manual - week0
        @test_throws DimensionMismatch refugia_distance(data, x0[:, 1:7], x0)
        # a small run: shape, prior bounds, reproducibility, independence of threads
        prior = RefugiaPrior(:miller2025)
        r = refugia_abc(data; nsamples = 4000, acceptance = 0.05, seed = 11, chunk = 500)
        @test size(r.accepted) == (200, 5) && r.names == (:β_sh, :β_un, :α, :ω, :μ)
        @test all(0 .<= r.accepted[:, 1:2] .<= prior.β_max) && all(r.accepted[:, 3] .<= 0.5)
        @test all(r.accepted[:, 5] .== q.μ) && issorted(r.distance)
        @test r.threshold == r.distance[end]
        @test r.accepted == refugia_abc(data; nsamples = 4000, acceptance = 0.05, seed = 11,
                                        chunk = 500).accepted
        s = refugia_abc_summary(r)
        @test s.refugia_effect ≈ 1 - s.β_un.mean / s.β_sh.mean
        @test s.αβ_sh.lo <= s.αβ_sh.median <= s.αβ_sh.hi
        rμ = refugia_abc(data; prior = RefugiaPrior(:miller2025; μ_max = 1.0), nsamples = 1000,
                         acceptance = 1.0)
        @test 0 < minimum(rμ.accepted[:, 5]) && maximum(rμ.accepted[:, 5]) < 1
    end

    @testset "stored ABC reproduction vs Tables 3 and S3" begin
        path = joinpath(pkgdir(C), "data", "miller2025", "abc_reproduction.toml")
        res = TOML.parsefile(path)
        # A published posterior mean is reproduced when it lies within 3 sd (over our
        # replicate ABC runs at the paper's settings) of our mean, allowing for the paper's
        # rounding to three decimals.
        consistent(block, k, pub) =
            abs(pub - block[k]["mean"]) <= 3 * block[k]["sd"] + 0.0005
        ex = res["exp_paper_priors"]                  # 20 runs, 10^6 draws, Table 3 priors
        for (k, pub) in (("β_sh", q.β_sh), ("β_un", q.β_un), ("α", q.α), ("αβ_sh", q.αβ_sh),
                         ("αβ_un", q.αβ_un), ("ω", q.ω), ("refugia_effect", q.refugia_effect),
                         ("vaccination_effect_sh", q.vaccination_effect_sh),
                         ("vaccination_effect_un", q.vaccination_effect_un))
            @test consistent(ex, k, pub)
        end
        # Shape of the posterior (run 1 medians and 95% intervals) against Table 3
        @test abs(ex["β_sh"]["median_run1"] - 0.505) < 0.1
        @test abs(ex["β_un"]["median_run1"] - 0.192) < 0.05
        @test abs(ex["ω"]["median_run1"] - 0.433) < 0.05
        @test abs(ex["α"]["hi_run1"] - 0.373) < 0.05
        @test abs(ex["β_sh"]["hi_run1"] - 1.876) < 0.1
        # The priors in the Zenodo code (β ~ U(0,1), ω ~ U(0,10)) do not give Table 3: the
        # published run used the priors the paper states.
        cp = res["exp_code_priors"]
        @test !consistent(cp, "β_sh", q.β_sh) && !consistent(cp, "ω", q.ω)
        @test cp["β_sh"]["max"] < 0.4 && cp["ω"]["min"] > 1.5
        # Table S3 (synthetic data). Consistent when μ = 0.033 (the paper's §S2.2 text)…
        qs = parameter_set(:miller2025; set = "synthetic_posterior")
        sp = res["syn_mu_paper"]
        for (k, pub) in (("β_sh", qs.β_sh), ("β_un", qs.β_un), ("α", qs.α), ("αβ_sh", qs.αβ_sh),
                         ("αβ_un", qs.αβ_un), ("ω", qs.ω))
            @test consistent(sp, k, pub)
        end
        # …but not with μ = 0.029, which the code's simulate_samples.R fixes for "syn".
        sc = res["syn_mu_code"]
        @test !consistent(sc, "β_un", qs.β_un) || !consistent(sc, "ω", qs.ω)

        # Identifiability (§3.1). The posterior of αβ is far narrower than its prior, while α
        # and β individually are barely updated, and log α and log β are negatively correlated.
        rg = res["ridge_exp"]
        @test rg["sd_ratio_αβ_sh"] < 0.2 && rg["sd_ratio_β_sh"] > 0.9 && rg["sd_ratio_α"] > 0.5
        @test rg["cor_logα_logβ_sh"] < -0.5
        # Simulation-estimation at the Table S2 values: the 95% intervals of β cover nearly the
        # whole prior; that of αβ is a small fraction of its prior. All cover the truth.
        se = res["simest_synthetic_truth"]
        @test se["β_sh"]["cri_over_prior_width"] > 0.9 && se["β_un"]["cri_over_prior_width"] > 0.8
        @test se["αβ_sh"]["cri_over_prior_width"] < 0.15
        @test all(se[k]["coverage95"] >= 0.9 for k in ("β_sh", "β_un", "α", "αβ_sh", "αβ_un", "ω"))
        # "αβ readily identifiable" (§3.1) holds for the width of its posterior, but its
        # posterior mean is not close to the truth: biased upwards by about 60% (0.0141 vs
        # 0.0089 week⁻¹) in our 20 data sets.
        @test_broken se["αβ_sh"]["mean_abs_rel_error"] < 0.25
        # Loss and recovery (§3.1): with μ free (U(0,1)), μ is weakly identified, and its
        # posterior mean is above the Eq. (3) value and the synthetic truth. ω stays close to
        # its prior. We see little of a μ–ω ridge: the posterior correlation is weak.
        for (k, μtrue) in (("mu_free_exp", q.μ), ("mu_free_syn", 0.029))
            m = res[k]
            @test m["mu_lo"] <= μtrue <= m["mu_hi"] && m["mu_mean"] > 2μtrue
            @test m["omega_hi"] - m["omega_lo"] > 1.5
            @test abs(m["cor_omega_mu"]) < 0.4
        end
    end
end
