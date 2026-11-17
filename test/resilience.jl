# Glue to DynamicalResilience.jl and ResilienceCore.jl: named equilibria, analytic and AD
# Jacobians agree, R₀ by the general next-generation method equals the closed forms, and
# resilience falls to zero at the Hopf threshold (Schmidt) and at the fold (load tipping).
using Chytridiomycosis
using Test

const DR = Chytridiomycosis.DynamicalResilience
const RC = Chytridiomycosis.ResilienceCore
const with_ = Chytridiomycosis.with

@testset "resilience glue" begin
    @testset "Schmidt 2017: named equilibria, R₀ and the Hopf threshold" begin
        p = SchmidtParameters(:schmidt2017)
        prob = SchmidtSLIProblem(p)
        for state in (:endemic, :disease_free)
            u = equilibrium_state(prob, state)
            du = zeros(3)
            Chytridiomycosis.schmidt_rhs!(du, u, p, 0.0)
            @test maximum(abs, du) < 1e-9
        end
        @test_throws ArgumentError equilibrium_state(prob, :no_such_state)
        J = RC.linearise(prob; state = :endemic)
        @test J isa RC.ContinuousJacobian
        @test J.J ≈ RC.linearise(to_ode_problem(prob), equilibrium_state(prob, :endemic)).J rtol = 1e-9
        @test RC.resilience(J) ≈ -Chytridiomycosis.equilibrium_stability(p)

        @test DR.reproduction_number(prob).R0 ≈ basic_reproduction_number(p)
        below = SchmidtSLIProblem(with_(p; K = 0.8))                       # K_threshold ≈ 1.0086
        @test DR.reproduction_number(below).R0 < 1
        @test RC.isstable(RC.linearise(below; state = :disease_free))
        @test !RC.isstable(RC.linearise(prob; state = :disease_free))     # R₀ ≈ 24.8 > 1
        @test_throws ArgumentError equilibrium_state(below, :endemic)

        # following the endemic state up in K: stability is lost at a Hopf bifurcation
        Ks = collect(2.0:0.25:14.0)
        start = equilibrium_state(SchmidtSLIProblem(with_(p; K = 2.0)), :endemic)
        path = DR.resilience_along(to_ode_problem(prob), (q, K) -> with_(q; K = K), Ks; u0 = start)
        @test all(path.converged)
        @test path.transition.kind === :hopf
        @test isapprox(path.transition.estimate, hopf_threshold(p); atol = 0.05)
        @test path.metrics.resilience[1] > 0 && path.metrics.resilience[end] < 0
    end

    @testset "Bsal with larvae: R₀ and the disease-free state" begin
        p = SchmidtParameters(:schmidt2017)
        stage = BsalLarvalStage(; carry_over = 0.5, adult_to_larva = 0.05, larva_to_larva = 0.01,
                                larva_to_adult = 0.02)
        for K in (0.5, 25.0)
            prob = BsalLarvalProblem(with_(p; K = K), stage)
            u = equilibrium_state(prob, :disease_free)
            du = zeros(5)
            Chytridiomycosis.bsal_larval_rhs!(du, u, (prob.p, prob.stage), 0.0)
            @test maximum(abs, du) < 1e-9
            r = DR.reproduction_number(prob)
            @test r.R0 ≈ basic_reproduction_number(prob.p, stage) rtol = 1e-10
            @test RC.isstable(RC.linearise(prob)) == (r.R0 < 1)
        end
    end

    @testset "load tipping: the low load disappears at a fold" begin
        m = LoadTipping(; ι = 0.5, λ = 0.25, σ = 0.2, c = 20.0, h = 100.0)
        ιc = critical_immigration(m)                                        # (√20 − √5)² = 5
        @test ιc ≈ 5.0
        low = RC.linearise(m)
        @test RC.isstable(low)
        @test !RC.isstable(RC.linearise(m; state = :threshold))
        @test low.J ≈ RC.linearise(to_ode_problem(m), equilibrium_state(m, :low_load)).J rtol = 1e-10
        ιs = collect(0.5:0.35:6.2)                                         # …, 4.70, 5.05, …
        setι(q, ι) = LoadTipping(ι, q.λ, q.σ, q.c, q.h, q.κ)
        path = DR.resilience_along(to_ode_problem(m), setι, ιs; u0 = equilibrium_state(m, :low_load))
        ok = ιs .< ιc
        @test all(path.converged[ok]) && all(path.stable[ok])
        @test !any(path.stable[.!ok])
        @test path.transition.kind === :fold
        @test 4.7 <= path.transition.estimate <= 5.05
        # resilience along the path is that of the analytic linearisation, and falls towards zero
        @test path.metrics.resilience[ok] ≈ [RC.resilience(RC.linearise(setι(m, ι))) for ι in ιs[ok]] rtol = 1e-8
        @test issorted(path.metrics.resilience[ok]; rev = true)
        @test_throws ArgumentError equilibrium_state(setι(m, 6.0), :low_load)
    end
end
