# Cook et al. 2022 Conserv Sci Pract 4:e12762 (doi:10.1111/csp2.12762): stochastic
# stage-structured projection of Cascades frogs with itraconazole treatment of metamorphs.
#
# Oracles:
#   - SI "Population simulations", Figure S3 and Table S2 (structure, the 28 survival and
#     fecundity values, SY = 0.8 SA, GS = 0.725, PF = 0.52, initial counts);
#   - Results 3.4 and Fig. 4: mean λ 0.9 without treatment and 1.2 with all metamorphs
#     treated every year; mean λ > 1 when 40% of metamorphs are treated in 40% of years,
#     with decreases in either needing compensatory increases in the other;
#   - closed forms: with no stochasticity the simulated means follow the projection matrix of
#     mean rates (one step from a given state).
using Chytridiomycosis
using Test
using Random: Xoshiro

@testset "Cook 2022 treatment projection" begin
    c = parameter_set(:cook2022)
    p = CookParameters(:cook2022)

    @testset "parameters (Table S2)" begin
        @test length(p.adult_survival) == 28 == length(p.fecundity)
        @test p.SU == c.su_mean && p.ST == c.st_mean
        @test p.yopy_ratio == c.yopy_survival_ratio
        @test p.growth == c.subadult_growth && p.female_fraction == c.female_fraction
        @test all(0 .< p.adult_survival .< 1) && all(p.fecundity .> 0)
        prob = CookProjectionProblem(p)
        @test prob.u0 == (1500, 50, 50, 100)
        @test_throws ArgumentError CookProjectionProblem(p; treated_fraction = 1.5)
        @test_throws ArgumentError CookProjectionProblem(p; breeding = :never)
    end

    @testset "one step against the mean matrix" begin
        # with fixed survival and fecundity pools of one value, the expected state after one
        # step equals the projection matrix times the initial state
        q = CookParameters(; SU = 0.2, ST = 0.6, adult_survival = [0.6], fecundity = [30.0],
                           yopy_ratio = 0.8, growth = 0.7, female_fraction = 0.5)
        u0 = (1000, 200, 300, 400)
        for (tr, breeding) in ((0.0, :after_survival), (1.0, :after_survival),
                               (0.0, :before_survival))
            prob = CookProjectionProblem(q; treated_fraction = tr, treated_years = tr,
                                         years = 1, u0, breeding)
            sol = solve(prob; nsims = 40_000, rng = Xoshiro(3))
            m = vec(sum(sol.states[2, :, :]; dims = 2)) ./ 40_000
            A = cook_projection_matrix(q; treated = tr, breeding)
            @test m ≈ A * collect(u0) rtol = 0.01
        end
    end

    @testset "growth rates (Results 3.4, Fig. 4)" begin
        g(q, r; kw...) = cook_growth_rate(solve(CookProjectionProblem(p; treated_fraction = q,
                                                                      treated_years = r, kw...);
                                                nsims = 1000, rng = Xoshiro(1)))
        u, t = g(0, 0), g(1, 1)
        @test round(u.mean; digits = 1) == c.lambda_untreated       # 0.939
        @test t.mean > 1.25                                          # 1.30; paper 1.2 (discrepant)
        th = g(c.threshold_treated_fraction, c.threshold_treated_years)
        @test th.mean > 1
        @test g(0.2, 0.4).mean < 1 && g(0.4, 0.2).mean < 1
        @test t.positive > th.positive > u.positive
        # schedule and breeding switches change the numbers only slightly
        @test g(0.4, 0.4; schedule = :random).mean > 1
        @test round(g(0, 0; breeding = :before_survival).mean; digits = 1) == 1.0
        grid = cook_treatment_grid(p; fractions = [0.0, 1.0], frequencies = [0.0, 1.0],
                                   nsims = 200, rng = Xoshiro(4))
        @test grid.λ[1, 1] < 1 < grid.λ[2, 2]
        @test grid.λ[1, 2] ≈ grid.λ[2, 1] rtol = 0.05   # no treatment either way
        # λ of a run is the geometric mean of its annual adult growth rates
        sol = solve(CookProjectionProblem(p); nsims = 5, rng = Xoshiro(5))
        A = sol.states[:, 4, 1]
        @test sol.λ[1] ≈ exp(sum(log.(A[2:end] ./ A[1:end-1])) / 10)
    end

    @testset "mean matrix" begin
        A0 = cook_projection_matrix(p)
        A1 = cook_projection_matrix(p; treated = 1)
        @test population_growth_rate(A1) > population_growth_rate(A0)
        E = elasticity_matrix(A0)
        @test sum(E) ≈ 1 atol = 1e-8
    end
end
