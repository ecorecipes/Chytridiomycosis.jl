# Wilber et al. 2016 (Methods Ecol Evol 7:1182): host–parasite IPM of Bd on Rana muscosa, and
# the generic IPM machinery. The oracles are:
# - the paper (eqns 7–13, Table 1, Results, Figs 4–5);
# - the authors' deposited fit (results/IPM_parameters_no_26_linear.rds, in the compendium as
#   set "deposited") and their results: curves decoded from the vector graphics of
#   results/{lambda_plot,var_mean_plot,mean_plot,r0_plot,mean_time_to_death}.pdf (to ~1e-5)
#   and cells of results/transmission_results/*.csv (github.com/mqwilber/ipm_for_parasites);
# - exact results for simple matrices, and the Briggs et al. (2010) load model.
using Chytridiomycosis
using Test
using Random: Xoshiro
using Distributions: median, quantile
using Distributions: Normal
using LinearAlgebra: I, eigvals
using LogDensityProblems: LogDensityProblems

@testset "Wilber 2016 host–parasite IPM" begin
    C = Chytridiomycosis
    q = parameter_set(:wilber2016)
    pT = WilberParameters(:wilber2016)                       # Table 1 as printed
    p = WilberParameters(:wilber2016; set = "deposited")     # what the figures used
    fit = fit_wilber2016()

    @testset "parameters (Table 1, compendium)" begin
        @test all(r -> r.verified || r.symbol in (:β_R0, :S_R0), parameter_records(:wilber2016))
        @test (pT.b00, pT.b10, pT.b11, pT.v01, pT.c01) == (11.824, -0.8605, 0.799, 5.92, -0.049)
        @test (pT.b02, pT.b12, pT.b22, pT.b03, pT.b13, pT.v03, pT.c03) ==
              (1.213, -0.472, -0.151, 0.642, 0.137, 0.59, 0.063)
        @test (pT.b04, pT.b14, pT.s0, pT.Δt) == (-1.66, 0.102, 1.0, 3.0)
        @test time_unit(WilberParameters) === :day
        m = IPMMesh(:wilber2016)
        @test (m.lower, m.upper, m.n) == (-5.0, 18.0, 100)
        @test IPMMesh(:wilber2016; n = 50).n == 50
        # the deposit agrees with Table 1 to the printed digits, except survival
        for s in (:b01, :b11, :b21, :c01, :b02, :b12, :b22, :b03, :b13, :c03, :b14)
            @test round(getfield(p, s); digits = 3) == getfield(pT, s)
        end
        @test round(p.v01; digits = 2) == pT.v01 && round(p.v03; digits = 2) == pT.v03
        @test round(p.b04; digits = 2) == pT.b04
        # Table 1's survival row is not the deposited (20 °C) fit
        @test_broken round(p.b00; digits = 3) == pT.b00
        @test WilberParameters(:wilber2016; b00 = 10.0).b00 == 10.0
        @test_throws ArgumentError WilberParameters(:wilber2016; v01 = -1.0)
        @test_throws ArgumentError WilberParameters(:wilber2016; s0 = 1.5)
        # φ(T): "between 0.22 and 0.6"
        @test round(infection_probability(p, 4); digits = 2) == 0.22
        @test round(infection_probability(p, 20); digits = 1) == 0.6
    end

    @testset "data and refit (Table 1, Results)" begin
        rows = load_wilber2016()
        @test length(rows) == 589
        @test length(unique(r.individual for r in rows)) == 20
        @test sort(unique(r.temperature for r in rows)) == [4, 12, 20, 26]
        @test all(length(unique(r.individual for r in rows if r.temperature == T)) == 5
                  for T in (4, 12, 20, 26))
        @test count(r -> !r.survived, rows) == 10 && all(r -> r.survived == !isnan(r.load_next), rows)
        @test round.(extrema(log(r.load) for r in rows if r.load > 0); digits = 2) == (-1.14, 13.15)
        # the refit reproduces the deposited fit (coefficients and covariances)
        for (s, v) in pairs(merge(fit.survival.coef, fit.growth.coef, fit.loss.coef,
                                  fit.initial_load.coef, fit.infection.coef))
            @test getfield(p, s) ≈ v rtol = 1e-7
        end
        @test fit.growth.σ2 ≈ p.v01 rtol = 1e-7
        @test fit.growth.c ≈ p.c01 rtol = 1e-6
        @test fit.initial_load.σ2 ≈ p.v03 rtol = 1e-7
        @test fit.initial_load.c ≈ p.c03 rtol = 1e-6
        @test fit.survival.vcov ≈ [22.451140902925417 -2.1036090585789786;
                                   -2.1036090585789786 0.19973592310035942] rtol = 1e-4
        @test fit.growth.vcov[3, 3] ≈ 0.0014484360034494586 rtol = 1e-6
        @test fit.loss.vcov[1, 1] ≈ 0.1595741931571224 rtol = 1e-6
        @test (fit.survival.n, fit.growth.n, fit.loss.n, fit.initial_load.n, fit.infection.n) ==
              (105, 258, 304, 43, 166)
        t = fit.tests
        @test round(t.growth_load; digits = 2) == 196.36
        @test round(t.growth_temperature; digits = 2) == 13.56
        @test round(t.growth_variance; digits = 1) == 9.8
        @test round(t.loss_temperature; digits = 3) == 14.555
        @test round(t.loss_load; digits = 3) == 23.701
        @test round(t.initial_temperature_t; digits = 2) == 2.53 && t.initial_temperature_df == 41
        @test round(t.initial_variance; digits = 2) == 6.00
        @test round(t.infection_temperature; digits = 4) == 6.0361
        # survival: the printed χ² = 12.197 and Table 1 come from the all-temperature fit
        @test_broken round(t.survival_load; digits = 3) == 12.197     # 20 °C fit: 8.579
        @test t.survival_load ≈ 8.5789 atol = 1e-3
        fa = fit_wilber2016(; survival_temperatures = (4, 12, 20))
        @test round(fa.survival.coef.b00; digits = 3) == 11.824
        @test round(fa.survival.coef.b10; digits = 4) == -0.8605
        @test fa.tests.survival_load ≈ 12.197 atol = 1e-3
        @test WilberParameters(fit).b00 == fit.survival.coef.b00
    end

    @testset "LogDensityProblems regressions" begin
        rows = load_wilber2016()
        D = C._wilber_datasets(rows)
        for comp in (:loss, :growth, :initial_load)
            m = WilberRegression(comp; rows)
            k = LogDensityProblems.dimension(m)
            @test length(parameter_names(m)) == k
            θ = if comp === :loss
                fit.loss.coef |> collect
            else
                d = getfield(D, comp)
                r = C._gls_fit(d.X, d.y, d.z; reml = false)
                vcat(r.coef, log(r.σ2), r.c)
            end
            f0 = LogDensityProblems.logdensity(m, θ)
            for i in 1:k        # the ML estimate is a stationary point
                e = zeros(k); e[i] = 1e-5
                g = (LogDensityProblems.logdensity(m, θ .+ e) - LogDensityProblems.logdensity(m, θ .- e)) / 2e-5
                @test abs(g) < 1e-3
                e[i] = 0.01
                @test LogDensityProblems.logdensity(m, θ .+ e) < f0
            end
        end
        @test_throws ArgumentError WilberRegression(:fecundity)
        @test_throws DimensionMismatch LogDensityProblems.logdensity(WilberRegression(:survival), [1.0])
    end

    @testset "generic IPM machinery" begin
        m = IPMMesh(0, 1, 4)
        @test meshwidth(m) == 0.25 && meshedges(m) == [0, 0.25, 0.5, 0.75, 1.0]
        @test meshpoints(m) == [0.125, 0.375, 0.625, 0.875] && length(m) == 4
        @test_throws ArgumentError IPMMesh(1, 0, 4)
        @test_throws ArgumentError IPMMesh(0, 1, 0)
        # a kernel with no death or loss keeps every host with :cdf and :boundary
        k = IPMKernel(; survival = x -> 1.0, growth = x -> Normal(x + 0.5, 1.0),
                      initial_load = Normal(0.0, 1.0))
        mesh = IPMMesh(-5, 18, 60)
        d = discretize(k, mesh; integration = :cdf, eviction = :boundary)
        @test all(sum(d.P; dims = 1) .≈ 1) && sum(d.recruitment) ≈ 1
        A = projection_matrix(d, 0.3)
        @test all(sum(A; dims = 1) .≈ 1)
        @test population_growth_rate(A) ≈ 1
        # midpoint rule loses exactly what eviction() measures, to quadrature error
        dm = discretize(k, mesh)
        ev = eviction(k, mesh)
        @test maximum(abs.(1 .- vec(sum(dm.P; dims = 1)) .- ev.growth)) < 3e-3
        @test_throws ArgumentError discretize(k, mesh; integration = :simpson)
        @test_throws ArgumentError discretize(k, mesh; eviction = :reflect)
        @test_throws DomainError projection_matrix(d, 1.5)
        # absorption times: a host surviving each step with probability 0.9 lives 10 steps
        @test absorption_times(fill(0.9, 1, 1)).mean ≈ [10.0]
        @test absorption_times(fill(0.9, 1, 1)).variance ≈ [0.9 / 0.01]    # geometric: q/p²
        # eigen-analysis: elasticities sum to 1; sensitivity matches a finite difference
        B = [0.5 0.2 0.1; 0.3 0.6 0.2; 0.1 0.1 0.6]
        @test sum(elasticity_matrix(B)) ≈ 1
        λ = population_growth_rate(B)
        @test λ ≈ maximum(real, eigvals(B))
        Bp = copy(B); Bp[2, 1] += 1e-7
        @test (population_growth_rate(Bp) - λ) / 1e-7 ≈ sensitivity_matrix(B)[2, 1] rtol = 1e-4
        @test sum(stable_distribution(B)) ≈ 1
        # R0 of a rank-one next-generation matrix equals its spectral radius
        dk = discretize(IPMKernel(p; temperature = 15), IPMMesh(:wilber2016; n = 40))
        F = (dk.uninfected_survival .* dk.recruitment) * fill(0.1, 40)'
        @test basic_reproduction_number(dk, 0.1) ≈ maximum(real, eigvals(F * inv(I - dk.P)))
        @test_throws DimensionMismatch basic_reproduction_number(dk, ones(3))
        mc = mesh_convergence(n -> population_growth_rate(p; temperature = 20,
                                                         mesh = IPMMesh(-5, 18, n)), [25, 50, 100])
        @test mc.converged && mc.relative_change[end] < 1e-6
        @test_throws ArgumentError mesh_convergence(identity, [10])
        s = load_summary(IPMMesh(0, 2, 2), [2.0, 1.0, 1.0])
        @test (s.total, s.prevalence, s.mean, s.variance) == (4.0, 0.5, 1.0, 0.25)
    end

    @testset "density-independent IPM against the deposit (Figs 4–5)" begin
        # λ from the authors' lambda_plot.pdf (4, 12, 16, 20 °C)
        for (T, λ) in ((4, 0.99896), (12, 0.99053), (16, 0.97275), (20, 0.93690))
            @test population_growth_rate(p; temperature = T) ≈ λ atol = 1e-5
        end
        # the 100- and 50-bin meshes agree (the code's comment), and cdf/boundary is close
        @test population_growth_rate(p; temperature = 20, mesh = IPMMesh(:wilber2016; n = 50)) ≈
              population_growth_rate(p; temperature = 20) atol = 1e-5
        @test population_growth_rate(p; temperature = 20, integration = :cdf,
                                     eviction = :boundary) ≈ 0.93690 atol = 2e-4
        # stable distribution conditional on infection (Fig. 5 and its inset)
        for (T, m, vm, xpk, ppk) in ((4, 2.03262, 2.61491, 1.325, 0.05662),
                                     (12, 5.17999, 1.65359, 5.695, 0.02938),
                                     (20, 8.24159, 0.85104, 8.915, 0.03661))
            s = stable_load_distribution(p; temperature = T)
            @test s.mean ≈ m atol = 2e-4
            @test s.vm_ratio ≈ vm atol = 1e-4
            i = argmax(s.probability)
            @test s.x[i] ≈ xpk atol = 1e-3
            @test s.probability[i] ≈ ppk atol = 2e-5
            @test sum(s.probability) ≈ 1
        end
        # mean time to death (Fig. 4a; days), from ln ZE −0.745 and 13.055
        for (T, a, b) in ((4, 2853.9, 1093.8), (12, 335.40, 92.23), (20, 69.52, 13.81))
            mt = mean_time_to_death(p; temperature = T)
            @test mt.infected[19] ≈ a rtol = 1e-3
            @test mt.infected[79] ≈ b rtol = 1e-3
        end
        # R0 (Appendix S1 code, 50 bins; results/r0_plot.pdf)
        m50 = IPMMesh(:wilber2016; n = 50)
        for (T, full, coll) in ((4, 0.22311, 0.22439), (12, 0.96225, 0.97733), (17, 1.68435, 1.71647))
            @test basic_reproduction_number(p; temperature = T, mesh = m50) ≈ full atol = 5e-5
            @test basic_reproduction_number(p; temperature = T, mesh = m50, method = :collapsed) ≈ coll atol = 5e-5
        end
        @test_throws ArgumentError basic_reproduction_number(p; temperature = 10, method = :exact)
        # elasticities (Fig. S5): growth and survival dominate at 20 °C
        e = vital_rate_elasticities(p; temperature = 20, mesh = m50)
        top = first.(sort(collect(pairs(e)); by = x -> -abs(x[2]))[1:4])
        @test Set(top) == Set([:b11, :b21, :b00, :b10])
        @test_throws ArgumentError vital_rate_elasticities(p; temperature = 20, of = :R0)
        # load growth is mean-reverting (b11 < 1): zero at x* = (b01 + b21 T)/(1 − b11)
        xs = (p.b01 + p.b21 * 20) / (1 - p.b11)
        @test abs(load_growth_rate(p; temperature = 20, x = xs)) < 1e-12
        @test xs ≈ log(1e4) atol = 0.05
    end

    @testset "projection, transmission and lowering" begin
        # eqn 13 as coded against the deposited CSVs (pop_size_post_summer_T, prop_uninfected_T)
        for (T, ω, β, N, u) in ((20, 0.01, 0.0, 84.4345088202785, 0.859658216033309),
                                (20, 0.3, 0.00117, 14.2994922417937, 0.0332037763285336),
                                (12, 0.01, 0.00117, 95.2194382243267, 0.818915627798886),
                                (12, 0.3, 0.0, 75.1505888934389, 0.231943846999764))
            s = solve(BdIPMProblem(p; temperature = T, transmission = LoadTransmission(ω, β; as_coded = true)))
            Nend = s.S[end] + sum(s.I[:, end])
            @test Nend ≈ N atol = 1e-8
            @test s.S[end] / Nend ≈ u atol = 1e-10
        end
        # eqn 13 as printed makes β matter where the code's log-of-total does not
        dec(tr) = (s = solve(BdIPMProblem(p; temperature = 20, transmission = tr));
                   1 - (s.S[end] + sum(s.I[:, end])) / 100)
        @test dec(LoadTransmission(0.01, q.β_max)) - dec(LoadTransmission(0.01, 0.0)) > 0.5
        @test dec(LoadTransmission(0.01, q.β_max; as_coded = true)) -
              dec(LoadTransmission(0.01, 0.0; as_coded = true)) < 0.2
        @test_throws ArgumentError LoadTransmission(1.0, 0.0)
        # with β = 0 and ω = φ(T), the density-dependent run is the density-independent one
        s1 = solve(BdIPMProblem(p; temperature = 16))
        s2 = solve(BdIPMProblem(p; temperature = 16,
                                transmission = LoadTransmission(infection_probability(p, 16), 0.0)))
        @test s1.S ≈ s2.S && s1.I ≈ s2.I && all(s1.φ .≈ infection_probability(p, 16))
        @test s1.t == 0:3:120 && length(s1.φ) == 40
        # a constant temperature function gives the constant-temperature run
        s3 = solve(BdIPMProblem(p; temperature = t -> 16.0))
        @test s3.S ≈ s1.S && s3.I ≈ s1.I
        # the long-run decline follows λ
        N = s1.S .+ vec(sum(s1.I; dims = 1))
        @test N[end] / N[end-1] ≈ population_growth_rate(p; temperature = 16) rtol = 2e-3
        # the DiscreteProblem map is the same step
        dp = to_discrete_problem(BdIPMProblem(p; temperature = 16))
        du = similar(dp.u0)
        dp.f(du, dp.u0, dp.p, 0.0)
        @test du ≈ vcat(s1.S[2], s1.I[:, 2])
        @test dp.tspan == (0.0, 40.0)
        @test_throws ArgumentError BdIPMProblem(p; temperature = 16, days = 10)
        @test_throws DimensionMismatch BdIPMProblem(p; temperature = 16, I0 = zeros(3))
    end

    @testset "parametric bootstrap (Fig. 4b)" begin
        B1 = growth_rate_bootstrap(fit; temperatures = [20], n = 60, rng = Xoshiro(1), as_coded = true)
        B2 = growth_rate_bootstrap(fit; temperatures = [20], n = 60, rng = Xoshiro(1))
        iqr(x) = quantile(vec(x), 0.75) - quantile(vec(x), 0.25)
        @test size(B1) == (60, 1)
        @test iqr(B1) > 2 * iqr(B2)          # the authors' slip widens the intervals
        @test 0.90 < median(B1) < 0.99 && 0.92 < median(B2) < 0.95   # Fig. 4b: ≈ 0.94
    end

    @testset "relation to the Briggs et al. (2010) load model" begin
        # In the Briggs model (here a Drawert et al. 2017 Group A point, whose deterministic
        # model is Briggs's), sporangial load grows at the linearised rate a = λ until the
        # lethal load Smax. The same vital rates in the generic IPM (growth by a Δt on ln load
        # with a tiny variance, survival 0 above ln Smax) give that trajectory and the time
        # to reach Smax.
        dp = DrawertParameters(:drawert2017; set = "groupA")
        a = load_growth_rate(dp)
        @test a > 0
        Δt, S0, Smax = 3.0, Float64(dp.S0), Float64(dp.Smax)
        k = IPMKernel(; survival = x -> x < log(Smax) ? 1.0 : 0.0,
                      growth = x -> Normal(x + a * Δt, 1e-3), initial_load = Normal(log(S0), 1e-3))
        h = a * Δt / 8
        mesh = IPMMesh(log(S0) - 4h, log(Smax) + 4h, round(Int, (log(Smax) - log(S0) + 8h) / h))
        d = discretize(k, mesh; integration = :cdf)
        n = d.recruitment ./ sum(d.recruitment)
        ys = meshpoints(mesh)
        for step in 1:3
            n = d.P * n
            @test sum(n .* ys) / sum(n) ≈ log(S0) + a * Δt * step atol = h
        end
        steps = absorption_times(d.P).mean
        j = argmin(abs.(ys .- log(S0)))
        @test steps[j] * Δt ≈ log(Smax / S0) / a atol = 2Δt
        # the Wilber growth function is not exponential: its daily rate falls with load
        @test load_growth_rate(p; temperature = 20, x = 0.0) > load_growth_rate(p; temperature = 20, x = 8.0)
    end
end
