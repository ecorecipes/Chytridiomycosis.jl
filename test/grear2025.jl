# Grear et al. 2025 Ecol Monogr 95:e70001 (doi:10.1002/ecm.70001): Bayesian hurdle-lognormal
# model of Bd prevalence and intensity across US amphibians. Oracles: the abstract (11,183
# swabs, 54 species, 609 sites; peak infection at intermediate temperatures; caudates more
# often infected than anurans) and the authors' code (USGS software release
# doi:10.5066/P1MXUJZM: model formula, priors, reference levels, row filters giving 11,183 and
# 11,172 rows). Closed forms: the intercept-only hurdle likelihood equals LoadHurdleModel's;
# finite-difference gradients; a three-tip phylogeny's correlation matrix. The full text was
# not available, so no printed coefficient is checked.
using Chytridiomycosis
using Test
using LogDensityProblems
using Random: Xoshiro, randn

@testset "Grear 2025 hurdle model" begin
    C = Chytridiomycosis

    @testset "data" begin
        s = grear2025_summary()
        @test s.swabs == 11183                  # abstract; 01–04 scripts' row count
        @test s.species == 54                   # abstract (53 named + "Taricha sp.")
        @test s.sites == 594                    # abstract prints 609 (discrepant target)
        @test s.prevalence_caudata > s.prevalence_anura
        d = load_grear2025()
        @test length(d) == 11172                # after the coordinate filter (scripts: 11170/11172)
        @test length(load_grear2025(; complete = false)) == 11189
    end

    @testset "hurdle likelihood against LoadHurdleModel" begin
        rng = Xoshiro(25)
        n = 400
        y = [rand(rng) < 0.4 ? exp(3 + 1.5randn(rng)) : 0.0 for _ in 1:n]
        X = randn(rng, n, 1)
        prob = HurdleRegressionProblem(y, X)
        @test LogDensityProblems.dimension(prob) == 5
        θ = [2.5, 0.0, 0.3, 0.0, log(1.4)]     # a, b, a_hu, b_hu, log σ
        # LoadHurdleModel's π is the probability of infection: logit π = −a_hu
        ll = LogDensityProblems.logdensity(LoadHurdleModel(y), [-0.3, 2.5, log(1.4)])
        hst(s, A) = -2 * log1p((exp(s) / A)^2 / 3) + s
        prior = -0.5 * (2.5 / 20)^2 - 0.5 * 0.3^2 + hst(log(1.4), 2.5)
        @test LogDensityProblems.logdensity(prob, θ) ≈ ll + prior
        c = hurdle_coefficients(prob, θ)
        @test c.prevalence_intercept ≈ -c.hu_intercept
        @test c.sigma ≈ 1.4
    end

    @testset "gradient" begin
        p = Grear2025Problem(load_grear2025(); random = (:species, :region))
        rng = Xoshiro(7)
        θ = 0.1 .* randn(rng, LogDensityProblems.dimension(p))
        lp, g = LogDensityProblems.logdensity_and_gradient(p, θ)
        @test lp ≈ LogDensityProblems.logdensity(p, θ)
        for j in [1, 3, 6, 9, 12, 17, 18, 19, 20, 40, LogDensityProblems.dimension(p)]
            e = zeros(length(θ)); e[j] = 1e-6
            fd = (LogDensityProblems.logdensity(p, θ .+ e) -
                  LogDensityProblems.logdensity(p, θ .- e)) / 2e-6
            @test g[j] ≈ fd rtol = 1e-4 atol = 1e-3
        end
    end

    @testset "phylogenetic correlation" begin
        ph = phylogenetic_correlation("((A_x:1,B_y:1.0):2,C_z:3);")
        @test ph.species == ["A x", "B y", "C z"]
        @test ph.correlation ≈ [1 2/3 0; 2/3 1 0; 0 0 1]
        f = grear2025_phylogeny_path()
        if f !== nothing
            tr = grear2025_phylogeny(f)
            @test length(tr.species) == 53     # the authors' check: 53 tips
            @test all(≈(1), [tr.correlation[i, i] for i in 1:53])
            @test issubset(tr.species, unique(load_grear2025().species))
            # correlated species effects: a positive-definite Cholesky factor exists
            p = Grear2025Problem(; random = (:species,), phylogeny = tr)
            @test p.levels[:species] == tr.species
        else
            @info "Grear 2025 phylogeny not found; skipping its checks"
        end
    end

    @testset "posterior mode: temperature peaks at intermediate values" begin
        fit = C._grear2025_reference_fit()
        opt = grear2025_temperature_optimum(fit)
        @test opt.curvature_prevalence < 0
        @test opt.curvature_intensity < 0
        lo, hi = extrema(fit.problem.data.tmax30)
        @test lo < opt.prevalence < hi
        @test lo < opt.intensity < hi
    end
end
