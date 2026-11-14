# Wilber et al. 2017 (Ecol Lett 20:1169): hybrid IPM of Bd on Rana muscosa with tadpoles, a
# zoospore pool, seasonality, density-dependent recruitment and demographic stochasticity,
# and the generic IPM additions (reservoir R0, stochastic projection). The oracles are:
# - the paper (Tables 1–2, Methods, Results, Fig. 2–3; Fig. 2 decoded from its vector graphics);
# - the authors' code (github.com/mqwilber/host_extinction) and the March 2017 preprint in
#   that repository (Fig. 3A, decoded), standing in for Appendices S1–S4 (HTTP 403);
# - R's smooth.spline on the Woodhams et al. (2008) zoospore death rates, run once in R:
#   predict(smooth.spline(c(4, 14.5, 23, 28), c(0.0027, 0.0035, 0.0160, 0.0410)), T);
# - closed forms (Appendix S2 eqn 16 for R0 of an SISZ model; multinomial and Poisson means).
# The validation targets wilber2017_* cover the remaining published numbers.
using Chytridiomycosis
using Test
using Random: Xoshiro
using LogDensityProblems: LogDensityProblems

@testset "Wilber 2017 hybrid IPM" begin
    C = Chytridiomycosis
    q = parameter_set(:wilber2017)
    p = Wilber2017Parameters(:wilber2017)                    # Table 1 as printed
    pc = Wilber2017Parameters(:wilber2017; set = "coded")    # the code's s0 and tadpole survival

    @testset "parameters (Tables 1–2, compendium)" begin
        @test (p.b00, p.b10, p.b01, p.b11, p.b21, p.v01, p.c01) ==
              (5.295, -2.595, 0.012, 0.799, 0.092, 5.92, -0.049)
        @test (p.b02, p.b12, p.b22, p.b03, p.b13, p.v03, p.c03) ==
              (1.213, -0.472, -0.151, 0.642, 0.137, 0.59, 0.063)
        @test (p.s0, p.μA, p.μT, p.sT1, p.sT2, p.sT3) == (0.999, 1.0, 1487.036, 0.987, 0.997, 0.997)
        @test (p.pT1, p.pT2, p.pT3, p.mT1, p.mT2, p.mT3, p.pA, p.λ, p.K) ==
              (1.0, 0.5, 0.0, 0.9, 0.9, 1.0, 0.25, 50.0, 5.0)
        @test isnan(p.ν) && p.Δt == 3.0 && time_unit(Wilber2017Parameters) === :day
        # the coded values round to Table 1
        @test pc.s0 ≈ 0.9^(1 / 122) && round(pc.s0; digits = 3) == p.s0
        @test pc.sT1^122 ≈ 0.2 atol = 1e-5
        @test pc.sT2^122 ≈ 0.7 atol = 1e-5
        @test round(pc.sT1; digits = 3) == p.sT1
        # Table 2
        @test (p.transmission, p.β0, p.β1) == (:zpool_dd, 5.29e-3, 7.52e-2)
        for (tr, β0, β1) in ((:zpool, 1.09e-2, 0.0), (:zpool_fd, 5.77e-3, 0.627),
                (:const, 8.07e-2, 0.0), (:const_dd, 4.18e-2, 5.25e-2), (:const_fd, 4.28e-2, 0.551))
            r = Wilber2017Parameters(:wilber2017; transmission = tr)
            @test (r.β0, r.β1) == (β0, β1)
        end
        @test [parameter_set(:wilber2017; set = s).WAIC for s in
               ("const", "const_dd", "const_fd", "zpool", "zpool_dd", "zpool_fd")] ==
              [401.1, 336.7, 336.6, 345.32, 301.18, 311.06]
        @test all(r -> r.verified || r.symbol in (:surv_mean, :surv_sd, :cv_β0, :cv_β1),
            parameter_records(:wilber2017))
        m = IPMMesh(:wilber2017)
        @test (m.lower, m.upper, m.n) == (-5.0, 18.0, 30)
        @test IPMMesh(:wilber2017; n = round(Int, q.mesh_n_R0)).n == 50
        # the survival z-score is not printed: it is the mean and sd of the 20 °C loads
        zs = Wilber2017Regression(:survival).zscore
        @test zs.mean ≈ q.surv_mean rtol = 1e-12
        @test zs.sd ≈ q.surv_sd rtol = 1e-12
        # on the raw ln load, Table 1's survival is close to the Wilber 2016 deposited fit
        @test p.b00 - p.b10 * p.surv_mean / p.surv_sd ≈ 11.592 atol = 1e-3
        @test p.b10 / p.surv_sd ≈ -0.8763 atol = 1e-4
        k = IPMKernel(p; temperature = 20)
        @test k.survival(p.surv_mean) ≈ 1 / (1 + exp(-p.b00))
        @test k.uninfected_survival == p.s0
        @test Wilber2017Parameters(:wilber2017; b11 = 0.7).b11 == 0.7
        @test_throws ArgumentError Wilber2017Parameters(:wilber2017; transmission = :nope)
        @test_throws ArgumentError Wilber2017Parameters(:wilber2017; v01 = -1.0)
        @test_throws ArgumentError Wilber2017Parameters(:wilber2017; ν = 1.0)
        @test_throws ArgumentError Wilber2017Parameters(:wilber2017; sT1 = 1.5)
    end

    @testset "zoospore survival (Woodhams et al. 2008 spline)" begin
        # R's smooth.spline (df = 4, λ = 2.5e-10) on the four death rates (per hour)
        R = [0.00269999998379, 0.00271576046225, 0.00302663749366, 0.00396859239067,
            0.00792451966134]
        ν = [zoospore_survival(p, T) for T in (4, 8, 12, 16, 20)]
        d = -log.(ν) ./ (24 * p.Δt)
        @test all(abs.(d ./ R .- 1) .< 5e-4)       # R's spline is the interpolant up to its conditioning
        @test ν[1] ≈ exp(-0.0027 * 72) rtol = 1e-12
        @test round.([zoospore_survival(p, T) for T in (4, 12, 20)]; digits = 3) == [0.823, 0.804, 0.565]
        @test issorted(ν; rev = true)
        @test zoospore_survival(Wilber2017Parameters(:wilber2017; ν = 0.05), 20) == 0.05
        @test zoospore_survival(p, 30) < zoospore_survival(p, 28)     # linear extrapolation
    end

    @testset "transmission functions (Table 2)" begin
        @test infection_probability(p; Z = 0, I = 0, A = 10) == 0
        @test infection_probability(p; Z = 1e4, I = 2, A = 10) ≈
              1 - exp(-(p.β0 * log(1e4 + 1) + p.β1 * 2) * 3)
        pf = Wilber2017Parameters(:wilber2017; transmission = :zpool_fd)
        @test infection_probability(pf; Z = 0, I = 2, A = 4) ≈ 1 - exp(-pf.β1 * 0.5 * 3)
        @test infection_probability(pf; Z = 0, I = 0, A = 0) == 0
        pk = Wilber2017Parameters(:wilber2017; transmission = :const)
        @test infection_probability(pk; Z = 0) == 0                      # pool switched on by Z > 0
        @test infection_probability(pk; Z = 1) ≈ 1 - exp(-pk.β0 * 3)
    end

    @testset "demography (Table 1, stable stages)" begin
        A = annual_demography_matrix(pc)
        @test round(population_growth_rate(A); digits = 2) == 1.46
        @test population_growth_rate(A) ≈ 1.46328 atol = 1e-5
        w = stable_distribution(A)
        @test round.(Int, 10 .* w ./ w[4]) == [85, 12, 3, 10]
        prob = HybridIPMProblem(pc; years = 3)
        @test findall(prob.reproduce) == [30, 152, 274]                 # spring, once a year
        @test wilber2017_temperature()(0) ≈ 4 && wilber2017_temperature()(183) ≈ 20
        @test prob.u0 == vcat(85, 12, 3, 10, zeros(30), 0)
        # disease-free (no shedding): Ricker equilibrium of 10–15 adults (Table 1, K)
        s = hybrid_summary(solve(HybridIPMProblem(C._with(pc; μT = 0.0); stochastic = false,
            years = 25)))
        a = s.adults[end-121:end]
        @test 10 <= minimum(a) && maximum(a) <= 15
        @test all(==(0), s.zoospores) && all(==(0), s.infected)
    end

    @testset "generic additions to the IPM machinery" begin
        # Appendix S2 eqn 16: SIS with a pool, one load class
        s0, sI, lI, ν, f, β0, β1, S = 0.99, 0.9, 0.2, 0.7, 5.0, 0.01, 0.05, 6.0
        d = DiscretizedIPM(IPMMesh(0, 1, 1), fill(sI * (1 - lI), 1, 1), [sI * lI], [1.0], s0)
        r = reservoir_reproduction_number(d; host_transmission = S * β1,
            pool_transmission = S * β0, shedding = f, pool_survival = ν)
        @test r.total ≈ s0 * S * β1 / (1 - sI * (1 - lI)) +
                        s0 * S * β0 * f / ((1 - ν) * (1 - sI * (1 - lI)))
        @test r.host ≈ basic_reproduction_number(d, S * β1)
        @test r.environment ≈ r.total - r.host
        @test_throws DomainError reservoir_reproduction_number(d; host_transmission = 1,
            pool_transmission = 1, shedding = 1, pool_survival = 1.0)
        @test_throws DimensionMismatch reservoir_reproduction_number(d; host_transmission = [1, 2],
            pool_transmission = 1, shedding = 1, pool_survival = 0.5)
        # stochastic projection: means T n + F n; conservation without births or deaths
        T = [0.5 0.0 0.0; 0.3 0.9 0.0; 0.0 0.05 1.0]
        F = [0.0 0.0 2.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
        n = [100.0, 50.0, 20.0]
        rng = Xoshiro(1)
        M = reduce(hcat, [demographic_projection(rng, T, F, n) for _ in 1:4000])
        mn = vec(sum(M; dims = 2)) ./ 4000
        @test all(abs.(mn .- (T + F) * n) .< 4 .* sqrt.((T + F) * n ./ 4000) .+ 0.5)
        @test all(isinteger, M)
        P = [0.2 0.0; 0.8 1.0]
        @test sum(demographic_projection(Xoshiro(2), P, zeros(2, 2), [37.0, 5.0])) == 42
        @test_throws DomainError demographic_projection(rng, T, F, [1.5, 0, 0])
        out = zeros(3)
        @test demographic_projection!(Xoshiro(3), out, T, F, n) === out
        @test out == demographic_projection(Xoshiro(3), T, F, n)
    end

    @testset "R0 with an environmental pool (Fig. 2)" begin
        r1 = basic_reproduction_number(p; temperature = 15, S = 1)
        r16 = basic_reproduction_number(p; temperature = 15, S = 16)
        @test r16.total ≈ 16 * r1.total && r16.host ≈ 16 * r1.host     # R0 ∝ S (Fig. 2b)
        @test r1.total > r1.host > 0
        @test r1.invasion ≈ 1 - 1 / r1.total
        # plug-in values at Table 1/2 (the figure plots posterior medians; see the targets)
        @test log(r1.total) ≈ 1.2311 atol = 1e-3
        @test log(r1.host) ≈ 0.1668 atol = 1e-3
        r4 = basic_reproduction_number(p; temperature = 4, S = 4)
        @test log(r4.host) < 0 < log(r4.total)        # the pool lets Bd invade at 4 °C
        # the authors' R0 code omits Δt and sheds ln(eˣ + 1): the model as simulated gives far more
        rm = basic_reproduction_number(p; temperature = 15, S = 1, as_coded = false)
        @test rm.host ≈ r1.host * p.Δt * p.s0 rtol = 0.02
        @test rm.total > 100 * r1.total
        # Fig. 2c shape: host R0 peaks between 16 and 19 °C at S = 4
        Ts = 4:0.5:20
        h = [basic_reproduction_number(p; temperature = T, S = 4).host for T in Ts]
        @test 16 <= Ts[argmax(h)] <= 19
        @test_broken abs(log(basic_reproduction_number(p; temperature = 20, S = 4).host) -
                         validation_target("wilber2017_fig2c_host_20C").published) < 0.03
        # ^ the figure is a posterior median over 500 draws, not the plug-in value; the median of
        #   draws from our refitted posteriors (scripts/wilber2017_fig2.jl) comes within 0.04
    end

    @testset "hybrid model: deterministic and stochastic" begin
        prob = HybridIPMProblem(pc; years = 2, stochastic = false)
        sol = solve(prob)
        @test size(sol.u) == (35, 245) && length(sol.φ) == 244
        @test sol.t[end] == 732
        @test all(>=(-1e-9), sol.u)
        dp = to_discrete_problem(prob)
        du = similar(prob.u0)
        dp.f(du, prob.u0, nothing, 0.0)
        @test du ≈ sol.u[:, 2]
        @test_throws ArgumentError to_discrete_problem(HybridIPMProblem(pc; years = 1))
        # tadpoles shed from the first step; adults get infected; the pool persists
        s = hybrid_summary(sol)
        @test s.zoospores[2] ≈ 100 * pc.μT
        @test maximum(s.prevalence[2:end]) > 0.5
        # without the removal term ψ the pool is (slightly) larger
        s2 = hybrid_summary(solve(HybridIPMProblem(pc; years = 2, stochastic = false,
            zoospore_removal = false)))
        @test all(s2.zoospores .>= s.zoospores .- 1e-6) && s2.zoospores[end] > s.zoospores[end]
        # stochastic: integer states, seeded, mean of one step ≈ the deterministic step
        ps = HybridIPMProblem(pc; years = 1)
        a, b = solve(ps; rng = Xoshiro(5)), solve(ps; rng = Xoshiro(5))
        @test a.u == b.u && all(isinteger, a.u)
        c = C._HybridCache(ps)
        buf = C._hybrid_buffers(35)
        out = zeros(35)
        u1 = sol.u[:, 30]                                     # a state with infected adults
        u1r = round.(u1)
        det = similar(u1r)
        C._hybrid_step!(det, buf, C._HybridCache(prob), u1r, 30, nothing)
        rng = Xoshiro(6)
        acc = zeros(35)
        for _ in 1:2000
            C._hybrid_step!(out, buf, c, u1r, 30, rng)
            acc .+= out
        end
        acc ./= 2000
        @test acc[1] ≈ det[1] rtol = 0.02                     # year-1 tadpoles after reproduction
        @test acc[4] + sum(acc[5:34]) ≈ det[4] + sum(det[5:34]) rtol = 0.02   # adults
        @test acc[end] ≈ det[end] rtol = 0.02                 # zoospores
        # extinction curve (100 runs, 10 years); the paper: c. 50% by 8 years
        ec = extinction_curve(HybridIPMProblem(pc; years = 10); nsims = 100, rng = Xoshiro(7))
        @test issorted(ec.probability) && ec.probability[1] == 0
        @test ec.probability[findfirst(>=(3), ec.years)] == 0
        @test 0.35 <= ec.probability[findfirst(>=(8), ec.years)] <= 0.7
        @test_throws ArgumentError extinction_curve(prob)
    end

    @testset "refits (Tables 1–2)" begin
        rows = load_wilber2017()
        @test length(rows) == 1378 && rows[1].swab_id == "MDF10_1_Jul25"
        @test count(r -> r.svl == "31.9?", rows) == 1
        m = Wilber2017Transmission(:const)
        @test length(m.infected) == 333 && sum(m.infected) == 124
        @test parameter_names(Wilber2017Transmission(:const_dd)) == (:β0, :β1)
        @test LogDensityProblems.logdensity(m, [-0.1]) == -Inf
        @test LogDensityProblems.logdensity(m, [0.08]) ≈
              sum(pointwise_loglikelihood(m, [0.08])) - log1p(0.08^2)
        s = posterior_summary(m)
        @test s.mean[1] ≈ 0.0807 atol = 5e-4                  # Table 2
        @test s.waic ≈ 401.1 atol = 0.15
        sd = posterior_summary(Wilber2017Transmission(:const_dd); n = 50)
        @test s.mean[1] > sd.mean[1]
        @test sd.mean ≈ [0.0418, 0.0525] atol = 6e-4
        @test sd.waic < s.waic - 60                          # density dependence fits better
        sv = posterior_summary(Wilber2017Regression(:survival); n = 80)
        @test sv.mean ≈ [5.295, -2.595] atol = 0.05          # Table 1, within Stan's MC error
        @test sv.mode[1] < sv.mean[1]                         # 4 deaths: a skewed posterior
        g = posterior_summary(Wilber2017Regression(:growth); draws = 5_000)
        @test all(abs.(g.mean .- [0.012, 0.799, 0.092, 5.92, -0.049]) .< [0.05, 0.01, 0.005, 0.3, 0.01])
        @test_broken round(g.mean[2]; digits = 3) == p.b11
        # ^ Table 1 prints the Wilber 2016 REML fit (0.799); the posterior mean that the authors'
        #   simulations used is 0.796 (and v01 5.84, not 5.92)
        @test round(fit_wilber2016().growth.coef.b11; digits = 3) == p.b11
        dr = posterior_draws(m, 2000; rng = Xoshiro(8))
        @test size(dr) == (2000, 1)
        @test sum(dr) / 2000 ≈ s.mean[1] atol = 1e-3
        @test_throws ArgumentError Wilber2017Transmission(:zpool)
        @test_throws ArgumentError Wilber2017Regression(:infection)
    end

    @testset "resistance and tolerance at the stasis load (Fig. 4)" begin
        s = stasis_load(p; temperature = 20)
        @test s.x ≈ (0.012 + 0.092 * 20) / (1 - 0.799)
        @test round(s.survival; digits = 2) == 0.97
        # less resistance (steeper growth) pushes the stasis load past the survival threshold
        @test stasis_load(C._with(p; b11 = 0.9); temperature = 20).survival < 0.05
        # more tolerance (higher survival intercept) raises survival at the same load
        @test stasis_load(C._with(p; b00 = 7.0); temperature = 20).survival > s.survival
        @test stasis_load(C._with(p; b11 = 1.05); temperature = 20).x == Inf
    end

    @testset "sensitivity analysis machinery (Fig. 3)" begin
        # the penalised logistic regression recovers the direction of known effects
        rng = Xoshiro(9)
        X = randn(rng, 2000, 3)
        y = [rand(rng) < 1 / (1 + exp(-(2x[1] - x[2]))) for x in eachrow(X)]
        w = sensitivity_weights(X, y).weights
        @test w[1] > 1.5 && -1.3 < w[2] < -0.7 && abs(w[3]) < 0.2
        @test all(abs.(sensitivity_weights(X, y; C = 1e-4).weights) .< abs.(w))
        r = extinction_sensitivity(pc; nsims = 6, years = 1, rng = Xoshiro(10))
        @test r.names == (:b00, :b10, :b01, :b21, :b11, :v01, :c01, :b03, :b13, :v03, :c03,
            :b02, :b22, :b12, :β0, :β1)
        @test size(r.X) == (6, 16) && length(r.extinct) == 6
        @test all(sign.(r.X) .== sign.([getfield(pc, s) for s in r.names])')
        rz = extinction_sensitivity(Wilber2017Parameters(:wilber2017; transmission = :zpool);
            nsims = 2, years = 1, rng = Xoshiro(11))
        @test length(rz.names) == 15
    end
end
