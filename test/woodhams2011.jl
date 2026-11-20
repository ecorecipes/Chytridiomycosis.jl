# Woodhams et al. 2011 Front Zool 8:8 (doi:10.1186/1742-9994-8-8): the Briggs load model with
# adaptive immunity, treatments and density reduction. Oracles: the Appendix equations and the
# extended Figure 6 legend (PMC3098159): the untreated frog of Fig. 6B hits Smax = 10000 "in 33
# days"; removing the frog from the pool on day 20 only delays death; lowering ν to 0.15 on
# day 20 lets the load fall; both together clear it. Closed forms: S(t) = S0 e^{at} +
# b(e^{at} − 1)/a and S* = −b/a. Equivalence with BriggsLoadProblem when θ = 0. Fig. 6C and 6D
# cannot be checked: r, c (and μ for 6D) are not printed.
using Chytridiomycosis
using Test

@testset "Woodhams 2011 adaptive immunity" begin
    C = Chytridiomycosis
    q = parameter_set(:woodhams2011; set = "fig6b")

    @testset "compendium" begin
        recs = parameter_records(:woodhams2011)
        @test all(r -> r.verified, recs)
        @test all(r -> r.provenance == "A", recs)
        @test Set(r.set for r in recs) == Set(["", "fig6b", "fig6c", "fig6d"])
        @test time_unit(ImmuneFrogParameters) === :day
        # r and c are not published: the keyed constructor refuses to invent them
        @test_throws ArgumentError ImmuneFrogParameters(:woodhams2011; set = "fig6c")
        @test_throws UndefKeywordError ImmunePopulationParameters(:woodhams2011; set = "fig6d")
    end

    p = ImmuneFrogParameters(:woodhams2011; set = "fig6b")
    @testset "single frog without immunity (Fig. 6B)" begin
        @test p.θ == 0 && p.p == q.ν
        @test infection_fraction(p, 5.0) == p.p
        # closed form against the ODE
        sol = solve(ConstantPoolFrogProblem(p; S0 = q.S0, tspan = (0.0, 25.0));
            reltol = 1e-10, abstol = 1e-8).sol
        @test sol(17.0)[1] ≈ constant_pool_load(p, 17.0; S0 = q.S0) rtol = 1e-7
        # (a) untreated: Smax "in 33 days" (32.18, during day 33)
        t_a = time_to_lethal_load(p; S0 = q.S0)
        @test 32 < t_a < 33
        @test solve(ConstantPoolFrogProblem(p; S0 = q.S0)).death_time ≈ t_a rtol = 1e-4
        run(tr) = solve(ConstantPoolFrogProblem(p; S0 = q.S0, tspan = (0.0, 60.0),
            treatments = tr); reltol = 1e-10, abstol = 1e-8)
        # (b) zoospore pool eliminated on day 20: death is delayed (self-reinfection, a > 0)
        b = run([FrogTreatment(q.t_treat; γZ = 0)])
        S20 = constant_pool_load(p, q.t_treat; S0 = q.S0)
        a = p.p * p.ηf - p.σ
        @test b.death_time ≈ q.t_treat + log(p.Smax / S20) / a rtol = 1e-5
        @test b.death_time > t_a
        # (c) ν lowered to 0.15: a < 0, the load falls toward b/|a| = 1500 (not to zero)
        c = run([FrogTreatment(q.t_treat; p = q.ν_treated)])
        pc = ImmuneFrogParameters(; p = q.ν_treated, γZ = p.γZ, ηf = p.ηf, σ = p.σ, Smax = p.Smax)
        @test isnan(c.death_time)
        @test only(constant_pool_equilibria(pc)) ≈ q.ν_treated * p.γZ / (p.σ - q.ν_treated * p.ηf)
        @test c.sol.u[end][1] ≈ 1500 + (S20 - 1500) * exp(-0.05 * 40) rtol = 1e-6
        # (d) both: the load decays to zero
        d = run([FrogTreatment(q.t_treat; γZ = 0, p = q.ν_treated)])
        @test d.sol.u[end][1] ≈ S20 * exp(-0.05 * 40) rtol = 1e-6
        @test d.sol.u[end][1] < c.sol.u[end][1]
        @test isempty(constant_pool_equilibria(p))         # a > 0: no equilibrium
        @test_throws ArgumentError constant_pool_load(
            ImmuneFrogParameters(; p = 0.5, γZ = 1.0, ηf = 1.0, σ = 0.2, Smax = 1e4, θ = 0.1,
                r = 1e-3, c = 0.1), 1.0)
    end

    @testset "single frog with immunity (illustrative r, c)" begin
        # r and c are NOT published; these values are for testing only.
        pc = ImmuneFrogParameters(:woodhams2011; set = "fig6c", r = 1e-3, c = 0.1)
        S★ = only(constant_pool_equilibria(pc))
        k = pc.θ * pc.r / pc.c
        @test pc.p * exp(-k * S★) * (pc.γZ + pc.ηf * S★) ≈ pc.σ * S★ rtol = 1e-8
        res = solve(ConstantPoolFrogProblem(pc; tspan = (0.0, 2000.0),
            treatments = [FrogTreatment(30; S = 0)]); reltol = 1e-9, abstol = 1e-9)
        @test res.sol.u[end][1] ≈ S★ rtol = 1e-3
        @test res.sol.u[end][2] ≈ pc.r * S★ / pc.c rtol = 1e-3
        @test res.sol(30.0)[1] == 0.0 || res.sol(30.0 + 1e-9)[1] < 1   # antifungal on day 30
        # antibodies only slow infection
        @test infection_fraction(pc, 10.0) ≈ pc.p * exp(-1.0)
    end

    @testset "population model (Briggs with immunity)" begin
        b = BriggsParameters(:briggs2010)
        N = 6
        S0 = collect(1.0:N) .* 10
        briggs = solve(BriggsLoadProblem(b; N, S0, Z0 = 50.0, tspan = (0.0, 400.0));
            reltol = 1e-9, abstol = 1e-9)
        imm = solve(ImmunePopulationProblem(ImmunePopulationParameters(b); N, S0, Z0 = 50.0,
            tspan = (0.0, 400.0)); reltol = 1e-9, abstol = 1e-9)
        @test isequal(round.(imm.death_time; digits = 3), round.(briggs.death_time; digits = 3))
        @test imm.sol(150.0)[end] ≈ briggs.sol(150.0)[end] rtol = 1e-5
        # Fig. 6D structure with illustrative μ, r, c (not published)
        q6 = parameter_set(:woodhams2011; set = "fig6d")
        pd = ImmunePopulationParameters(:woodhams2011; set = "fig6d", μ = 0.1, r = 1e-3, c = 0.1)
        @test pd.γ == q6.γ && pd.η == q6.η && pd.p == q6.p
        n0 = Int(q6.N0)
        untreated = solve(ImmunePopulationProblem(pd; N = n0, Z0 = q6.Z0, tspan = (0.0, 100.0)))
        treated = solve(ImmunePopulationProblem(pd; N = n0, Z0 = q6.Z0, tspan = (0.0, 100.0),
            reductions = [DensityReduction(q6.t_treat, Int(q6.N_treated))]))
        @test count(!isnan, treated.removal_time) == n0 - Int(q6.N_treated)
        @test all(isnan, untreated.removal_time)
        # removing frogs slows the pool: the survivors' load is lower afterwards
        @test treated.sol(40.0)[1] < untreated.sol(40.0)[1] || all(isnan, treated.death_time)
    end
end
