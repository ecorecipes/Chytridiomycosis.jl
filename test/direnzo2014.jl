# DiRenzo et al. 2014 PLOS ONE 9:e93356 (doi:10.1371/journal.pone.0093356): log-linear mixed
# models of Bd intensity and zoospore output in Atelopus zeteki. Oracles: the printed rates
# (35.4%, SE 0.05; 13.1%, SE 0.04) against the printed t statistics (t50 = 7.018, t46 = 3.164);
# Table 1 against the text (mean intensity at death 4,334,422, SE 1,156,576 in the Results);
# closed-form multivariate normal log densities; recovery of simulated parameters by the ML
# fit. No time series were published, so the fitted trajectories cannot be checked on data.
using Chytridiomycosis
using Test
using Random: Xoshiro
using Distributions: MvNormal, logpdf
using LogDensityProblems

@testset "DiRenzo 2014 Atelopus load and output" begin
    p = DiRenzoParameters(:direnzo2014)

    @testset "compendium and printed rates" begin
        @test all(r -> r.verified, parameter_records(:direnzo2014))
        q = parameter_set(:direnzo2014)
        @test p.r_intensity == q.r_intensity && p.r_output == q.r_output
        @test daily_percent_increase(p) ≈ 100q.r_intensity
        @test daily_percent_increase(p; as_printed = false) ≈ 100expm1(q.r_intensity)
        # the "percent" is 100 r on the ln scale: r/SE agrees with the printed t; ln(1.354)
        # would not
        @test 6.44 < direnzo_t_statistic(p) < 7.87
        @test !(6.44 < log1p(q.r_intensity) / p.se_intensity < 7.87)
        @test 2.91 < direnzo_t_statistic(p; response = :output) < 3.74
        @test direnzo_trajectory(p, 10.0; y0 = 2.0) ≈ 2exp(10q.r_intensity)
        @test_throws ArgumentError direnzo_trajectory(p, 1.0; y0 = 1.0, response = :swab)
    end

    @testset "Table 1 against the text" begin
        frogs = load_direnzo2014()
        @test length(frogs) == 8 && count(f -> f.prior, frogs) == 3
        s = direnzo_summary(frogs)
        @test s.intensity.mean == 4_334_422                # Results and abstract
        @test round(s.intensity.se) == 1_156_576           # Results (abstract: 1,236,431)
        @test (s.intensity.min, s.intensity.max) == (520_436, 9_584_158)
        # the inconsistencies recorded as discrepant targets
        @test s.days.mean ≈ 25.625                         # abstract: 28.88
        @test s.days.min == 18                             # abstract: "between 22 and 33"
        @test s.output.mean ≈ 23.2875                      # text: 23.55
        lr = direnzo_logrank(frogs)
        @test lr.chisq ≈ 1.2974 atol = 1e-4                # Fig. 1: χ² = 0.7
        @test lr.p > 0.05                                  # same conclusion
    end

    @testset "log density" begin
        t = [7.0, 10.0, 14.0, 17.0, 21.0]
        Σ = direnzo_covariance(t; σ_id = 0.8, σ = 0.5, ϕ = 0.6)
        @test Σ[1, 3] ≈ 0.64 + 0.25 * 0.36
        Σt = direnzo_covariance(t; σ_id = 0.8, σ = 0.5, ϕ = 0.6, correlation = :time)
        @test Σt[1, 2] ≈ 0.64 + 0.25 * 0.6^3
        y = [9.0, 10.5, 11.0, 12.5, 13.0]
        m = DiRenzoLoadModel([DiRenzoSeries("a", false, t, y)])
        x = [6.0, 0.35, log(0.8), log(0.5), atanh(0.6)]
        @test LogDensityProblems.dimension(m) == 5
        @test LogDensityProblems.logdensity(m, x) ≈ logpdf(MvNormal(6.0 .+ 0.35 .* t, Σ), y)
        mp = DiRenzoLoadModel([DiRenzoSeries("a", true, t, y)]; prior_effect = true)
        @test LogDensityProblems.logdensity(mp, [6.0, 0.35, 0.5, x[3:end]...]) ≈
              logpdf(MvNormal(6.5 .+ 0.35 .* t, Σ), y)
        @test_throws ArgumentError DiRenzoSeries("a", false, [2.0, 1.0], [1.0, 1.0])

        # recovery: simulate 40 frogs from the printed slope and illustrative variances
        rng = Xoshiro(2014)
        r = p.r_intensity
        series = map(1:40) do i
            ti = collect(7.0:3.5:28.0)
            Σi = direnzo_covariance(ti; σ_id = 0.7, σ = 0.6, ϕ = 0.5)
            DiRenzoSeries(string(i), false, ti, rand(rng, MvNormal(5.0 .+ r .* ti, Σi)))
        end
        fit = maximum_likelihood(DiRenzoLoadModel(series))
        @test fit.estimates.r ≈ r atol = 0.03
        @test fit.estimates.ϕ ≈ 0.5 atol = 0.25
        @test fit.loglik ≈ LogDensityProblems.logdensity(fit.model, fit.x)
    end
end
