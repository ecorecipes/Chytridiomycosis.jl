# Doddington et al. 2013 Ecology 94:1795 (doi:10.1890/12-1270.1): Mallorcan midwife toads and
# Bd zoospores in one pond. Oracles: Table B1 and eqs B.1-B.8 (Appendix B), the authors'
# Berkeley Madonna code (Supplement 1, figshare doi:10.6084/m9.figshare.c.3306075, CC0), Table 1
# of the published version, and closed forms (zoospore decay, disease-free host dynamics).
# Stochastic outcomes against Fig. 3 and Appendix B are validation targets
# (data/validation/targets.toml); the full Fig. 3 grid is scripts/doddington2013_fig3.jl.
using Chytridiomycosis
using Test
using Random: Xoshiro

@testset "Doddington 2013 midwife toads and zoospores" begin
    C = Chytridiomycosis

    @testset "compendium" begin
        recs = parameter_records(:doddington2013)
        @test all(r -> r.verified, recs)
        p = DoddingtonParameters(:doddington2013)
        @test C.time_unit(DoddingtonParameters) === :year
        # the code's fractions of the year round to Table B1's days of the year
        q = DoddingtonParameters(:doddington2013; timing = :table)
        for f in (:tB1, :tB2, :tM1, :tM2, :tOM)
            @test round(getfield(p, f) * 365) == round(getfield(q, f) * 365)
        end
        @test p.tAge == q.tAge
        @test DoddingtonParameters(:doddington2013; set = "fig3a").α ==
              DoddingtonParameters(:doddington2013; set = "appB_b").α
        @test DoddingtonParameters(:doddington2013; set = "fig3b").α <
              DoddingtonParameters(:doddington2013).α
        @test DoddingtonParameters(:doddington2013; ν = 1e-5).ν == 1e-5
        @test_throws ArgumentError DoddingtonParameters(:doddington2013; set = "nope")
        @test_throws ArgumentError DoddingtonParameters(:doddington2013; timing = :x)
        # Fig. 3 axes: 10 log-spaced ticks; the printed labels (3 significant digits)
        g = doddington_grid()
        @test length(g.ν) == length(g.ρ) == 10
        @test round.(g.ν[[2, 5, 9]]; sigdigits = 3) == [3.73e-6, 1.14e-5, 5.01e-5]
        @test round.(g.ρ[[2, 5, 8]]; sigdigits = 3) == [2.06e4, 1.79e5, 1.56e6]
    end

    @testset "exposure trial (Table 1)" begin
        e = doddington_exposure_trial()
        @test round(e.control; digits = 2) == 0.40
        @test round(e.CAPE; digits = 2) == 1.01
        @test round(e.CAPE_excess.rate; digits = 2) == 0.61
        @test round.((e.CAPE_excess.lo, e.CAPE_excess.hi); digits = 2) == (-0.16, 1.37)
        @test round.((e.GPL_excess.lo, e.GPL_excess.hi); digits = 2) == (1.41, 4.54)
        # Table B1's α (2.97) is the rounded excess rate; Table 1 prints 2.98 (discrepant target)
        @test round(e.GPL_excess.rate; digits = 2) == DoddingtonParameters(:doddington2013).α
        # the high-virulence α of Fig. 3a is the upper limit of the BdGPL interval
        @test round(e.GPL_excess.hi; digits = 2) ==
              DoddingtonParameters(:doddington2013; set = "fig3a").α
        @test excess_death_rate(2, 1, 1, 1; z = 0) == (; rate = 1.0, lo = 1.0, hi = 1.0)
    end

    p = DoddingtonParameters(:doddington2013)

    @testset "deterministic model (eqs B.1-B.8)" begin
        # Zoospores alone decay at μZ (B.8 with no infected tadpoles)
        u0 = zeros(14); u0[14] = 1000.0
        prob = DoddingtonODEProblem(p; u0, t_intro = 100.0, tspan = (0.0, 0.05))
        sol = solve(prob)
        @test sol.u[end][14] ≈ 1000 * exp(-p.μZ * 0.05) rtol = 1e-6
        # without Bd the host settles on an annual cycle below K
        prob = DoddingtonODEProblem(p; t_intro = 1e6, tspan = (0.0, 60.0))
        sol = solve(prob)
        A = [sol(t)[1] for t in 50.0:1.0:59.0]
        @test all(a -> 0 < a < p.K, A)
        @test maximum(A) - minimum(A) < 1e-3 * maximum(A)
        @test sol(59.5)[1] ≈ sol(58.5)[1] rtol = 1e-3
        @test all(x -> x == 0, sol.u[end][8:14])
        # the jumps: tadpoles become overwintering tadpoles at tM2
        s = solve(DoddingtonODEProblem(p; t_intro = 1e6, tspan = (0.0, 1.5)))
        @test s(p.tM2 + 1e-6)[2] ≈ 0 atol = 1e-6
        @test s(p.tM2 - 1e-6)[2] > 10
        @test s(p.tM2 + 1e-6)[3] ≈ s(p.tM2 - 1e-6)[3] + s(p.tM2 - 1e-6)[2] rtol = 1e-4
        # ageing at tAge of the next year: J → JM
        @test s(1 + p.tAge + 1e-6)[4] ≈ 0 atol = 1e-6
        # after Bd introduction the infection spreads and juveniles die faster
        sB = solve(DoddingtonODEProblem(p; t_intro = 20.0, tspan = (0.0, 30.0)))
        @test sB(20.0 + 1e-9)[14] ≈ DoddingtonODEProblem(p).Z_intro rtol = 1e-6
        @test sB(30.0)[14] > 0
        # printed (constant εJ) versus coded (seasonal) maturation
        sP = solve(DoddingtonODEProblem(p; maturation = :as_printed, t_intro = 1e6,
            tspan = (0.0, 10.0)))
        @test sP.u[end][1] != sol(10.0)[1]
        @test_throws ArgumentError DoddingtonODEProblem(p; maturation = :other)
        # B.1-B.8 on the right-hand side: mass balance of tadpole infection
        du = zeros(14)
        u = collect(1.0:14.0)
        C.doddington_rhs!(du, u, (p, :as_printed), 0.55)
        λ = p.ν * u[14]
        @test du[8] ≈ λ * u[2] - (p.μL + p.εL) * u[8]
        @test du[14] ≈ p.ρ * (u[8] + u[9]) - p.μZ * u[14]
        @test du[1] ≈ (1 - u[1] / p.K) * p.εJ * (u[5] + u[11] + u[7] + u[13]) - p.μA * u[1]
    end

    @testset "stochastic model (Supplement 1)" begin
        prob = DoddingtonStochasticProblem(p)
        @test prob.dt == 0.0005 && prob.tspan == (0.0, 50.0) && prob.u0[1] == 150
        r1 = solve(prob; rng = Xoshiro(1))
        r2 = solve(prob; rng = Xoshiro(1))
        @test r1.u[end] == r2.u[end]                  # reproducible with a seed
        @test r1.outcome in (:host_extinction, :Bd_extinction, :coexistence)
        @test all(u -> all(>=(0), u), r1.u)
        # zoospores are set to Z_intro at the first step at or after t_intro
        i = findfirst(t -> t > prob.t_intro + 1e-9, r1.t)
        @test r1.u[i - 1][14] == 0 && r1.u[i][14] > 0
        # before Bd, adults fluctuate around the levels of Appendix B (about 50-150)
        A = [u[1] for (t, u) in zip(r1.t, r1.u) if 5 < t < 15]
        @test 30 < sum(A) / length(A) < 150
        # no Bd: never infected
        rn = solve(DoddingtonStochasticProblem(p; t_intro = 1e6, tspan = (0.0, 5.0));
            rng = Xoshiro(3))
        @test all(u -> sum(u[8:14]) == 0, rn.u)
        @test rn.outcome === :Bd_extinction
        # tadpoles move to the overwintering class once a year
        k = findfirst(t -> t >= 1 + p.tM2, rn.t)
        @test rn.u[k][2] == 0 && rn.u[k - 1][2] > 0
        # outcomes
        @test doddington_outcome(zeros(Int, 14)) === :host_extinction
        @test doddington_outcome([1; zeros(Int, 13)]) === :Bd_extinction
        @test doddington_outcome([1; zeros(Int, 12); 5]) === :coexistence
        o = doddington_outcomes(p; runs = 5, rng = Xoshiro(2), tspan = (0.0, 20.0))
        @test o.runs == 5 && o.host_extinction + o.Bd_extinction + o.coexistence == 5
        # the mean of the stochastic model follows the deterministic model before Bd
        T = 3.0
        det = solve(DoddingtonODEProblem(p; t_intro = 1e6, tspan = (0.0, T)))
        sto = [solve(DoddingtonStochasticProblem(p; t_intro = 1e6, tspan = (0.0, T),
            saveat = Inf); rng = Xoshiro(s)).u[end][1] for s in 1:200]
        @test sum(sto) / 200 ≈ det.u[end][1] rtol = 0.05
    end
end
