# N5 (new model, not a published result): swab and eDNA observation model for the Briggs et
# al. 2010 load model, and identifiability of its parameters from these data streams.
using LinearAlgebra: diag, eigvals, Symmetric

@testset "N5 swab and eDNA identifiability (new)" begin
    LDP = Chytridiomycosis.LogDensityProblems
    p = BriggsParameters(:briggs2010)
    # assay scales are illustrative assumptions (no published values), as in the wiki page
    obs = SurveyObservation(ρ = 10.0, κ = 1.0, σe = 0.5)
    days = [5.0, 10, 15, 20, 30, 40]
    B = (:γ, :ν, :f, :μ, :η, :σ)
    design(; Ns = [10, 40, 60], swabs = 10, pool = 1, edna = 3) =
        [SurveyPond(; N, days, swabs = min(swabs, N ÷ pool), pool, edna) for N in Ns]
    Zs(n) = Tuple(Symbol("Z0_", i) for i in 1:n)

    @testset "moments" begin
        prob = BriggsMomentProblem(p; N = 30, Z0 = 1e4, tspan = (0.0, 40.0))
        s = solve(prob; saveat = days)
        # the mean is the deterministic Briggs model with every frog starting at zero
        b = solve(BriggsLoadProblem(p; N = 30, Z0 = 1e4, tspan = (0.0, 40.0));
            saveat = days, reltol = 1e-10, abstol = 1e-10)
        @test all(b.alive)
        @test s.m ≈ [u[1] for u in b.sol.u] rtol = 1e-6
        @test s.Z ≈ [u[end] for u in b.sol.u] rtol = 1e-6
        # the exact solution matches the ODE lowering
        ode = solve(to_ode_problem(prob), Chytridiomycosis.Tsit5(); saveat = days,
            reltol = 1e-10, abstol = 1e-10)
        @test s.w ≈ [u[3] for u in ode.u] rtol = 1e-6
        # at the start the aggregation is the BDI shape ι/λ = γZ0/(Vηf) (N1)
        s0 = solve(prob; saveat = [1e-5])
        @test s0.k[1] ≈ p.γ * 1e4 / (p.V * p.η * p.f) rtol = 1e-3
        @test_throws ArgumentError BriggsMomentProblem(p; N = 0, Z0 = 1.0)
        @test_throws ArgumentError SurveyPond(; N = 5, days, swabs = 6)
        @test_throws ArgumentError SurveyPond(; N = 5, days = [0.0, 1.0])
    end

    @testset "hand derivation: densities identify γ/V, μ, r0, νη(1−f)" begin
        c = identifiable_combinations(p)
        J(N) = Chytridiomycosis.briggs_linear_jacobian(p, N)
        tr(N) = J(N)[1, 1] + J(N)[2, 2]
        dt(N) = J(N)[1, 1] * J(N)[2, 2] - J(N)[1, 2] * J(N)[2, 1]
        g = -(tr(40) - tr(10)) / 30                   # slope of the trace
        amμ = tr(10) + g * 10                          # r0 − μ
        aμ = -(dt(10) - (dt(40) - dt(10)) / 30 * 10)   # −(intercept of det) = r0·μ
        roots = sort([(amμ + sqrt(amμ^2 + 4aμ)) / 2, (amμ - sqrt(amμ^2 + 4aμ)) / 2])
        @test g ≈ c.γ_V
        @test roots ≈ sort([c.r0, -c.μ])               # the pair {r0, −μ}
        bc = -(dt(40) - dt(10)) / 30 - c.r0 * g        # slope of det = −(r0 g + νgη(1−f))
        @test bc / g ≈ c.νη_1mf
    end

    @testset "structural ranks" begin
        rank(ps, free; aggregation = true) = sensitivity_rank(
            SurveyModel(SurveySample[], ps; p, obs, Z0 = fill(1e4, length(ps)), free);
            aggregation)
        A = rank(design(; Ns = [40], edna = 0), (B..., :ρ, Zs(1)...); aggregation = false)
        @test (A.rank, A.nparams) == (3, 8)             # trace, determinant, amplitude
        Bd = rank(design(; edna = 0), (B..., :ρ, Zs(3)...); aggregation = false)
        @test (Bd.rank, Bd.nparams) == (7, 10)
        # every unidentified direction leaves the hand-derived combinations unchanged
        for v in Bd.null
            q(ε) = identifiable_combinations(BriggsParameters(; V = p.V,
                γ = p.γ * exp(ε * v.γ), ν = p.ν * exp(ε * v.ν), f = p.f * exp(ε * v.f),
                μ = p.μ * exp(ε * v.μ), η = p.η * exp(ε * v.η), σ = p.σ * exp(ε * v.σ),
                Smax = p.Smax))
            for n in (:γ_V, :μ, :r0, :νη_1mf)
                @test (q(1e-6)[n] - q(-1e-6)[n]) / 2e-6 / abs(q(0.0)[n]) ≈ 0 atol = 1e-6
            end
        end
        C = rank(design(; swabs = 0), (B..., :κ, Zs(3)...))
        @test (C.rank, C.nparams) == (7, 10)             # eDNA alone: as swab totals
        E = rank(design(), (B..., :ρ, Zs(3)...); aggregation = false)
        @test (E.rank, E.nparams) == (8, 10)             # + eDNA (absolute): ρν
        F = rank(design(), (B..., :ρ, Zs(3)...))
        @test (F.rank, F.nparams) == (9, 10)             # + aggregation: ηf
        @test F.null[1].ρ ≈ -F.null[1].ν atol = 1e-3     # the swab-scale symmetry ρ ↔ 1/ν
        G = rank(design(), (B..., Zs(3)...))
        @test G.rank == G.nparams == 9                   # both assays absolute: all identified
        J1 = rank(design(; Ns = [40]), (B..., Zs(1)...))
        @test (J1.rank, J1.nparams) == (6, 7)            # one density never separates γ, μ
        L = rank(design(; Ns = [10, 40]), (B..., Zs(2)...))
        @test L.rank == L.nparams == 8                   # two densities suffice
    end

    @testset "likelihood and information" begin
        ps = design()
        data = simulate_survey(Xoshiro(3), p, obs, ps; Z0 = fill(1e4, 3))
        @test count(s -> s.kind === :swab, data) == 6 * (10 + 10 + 10)
        @test count(s -> s.kind === :edna, data) == 6 * 3 * 3
        free = (B..., :σe, Zs(3)...)
        m = SurveyModel(data, ps; p, obs, Z0 = fill(1e4, 3), free)
        @test parameter_names(m) == free
        @test LDP.dimension(m) == 10
        @test LDP.capabilities(typeof(m)) isa LDP.LogDensityOrder{0}
        x = Chytridiomycosis.survey_start(m)
        @test isfinite(LDP.logdensity(m, x))
        xf = copy(x); xf[3] = 0.1                          # f > 1
        @test LDP.logdensity(m, xf) == -Inf
        @test_throws DimensionMismatch LDP.logdensity(m, x[1:3])
        @test_throws ArgumentError SurveyModel(data, ps; p, obs, Z0 = fill(1e4, 3), free = (:β,))
        # a pooled tube is the sum of the frogs in it
        d = Chytridiomycosis.survey_swab_distribution(4.0, 2.0, 10.0, 5)
        @test mean(d) ≈ 5 * 10 * 4.0
        # negative-binomial shape information agrees with the variance of the score
        r, M = 3.0, 50.0
        nb = Chytridiomycosis.NegativeBinomial(r, r / (r + M))
        score(y) = r * (sum(1 / (r + j) for j in 0:(y - 1); init = 0.0) + log(r / (r + M)) +
                        (M - y) / (r + M))
        ys = rand(Xoshiro(1), nb, 100_000)
        @test Chytridiomycosis.nb_log_shape_information(r, M) ≈
              sum(score.(ys) .^ 2) / length(ys) rtol = 0.03
        # expected and observed information agree on the diagonal
        I = expected_information(m)
        Io = observed_information(m)
        @test diag(Io)[1:6] ≈ diag(I)[1:6] rtol = 0.05
        # pooling loses information about f (aggregation), not about μ (means)
        se(ps) = sqrt.(diag(inv(expected_information(
            SurveyModel(SurveySample[], ps; p, obs, Z0 = fill(1e4, 3), free)))))
        s1, s10 = se(design()), se(design(; swabs = 1, pool = 10))
        @test s10[3] > 1.4 * s1[3]
        @test s10[4] ≈ s1[4] rtol = 0.05
        # swab totals carry no information on ηf: the matrix is singular
        It = expected_information(m; aggregation = false)
        @test minimum(eigvals(Symmetric(It))) < 1e-10 * maximum(eigvals(Symmetric(It)))
    end

    @testset "fit to synthetic data" begin
        ps = design(; Ns = [10, 60])
        data = simulate_survey(Xoshiro(5), p, obs, ps; Z0 = fill(1e4, 2))
        free = (:f, :η, :σ, :σe, Zs(2)...)
        guess = BriggsParameters(; V = p.V, γ = p.γ, ν = p.ν, f = 0.1, μ = p.μ, η = 10.0,
            σ = 0.3, Smax = p.Smax)
        m = SurveyModel(data, ps; p = guess, obs, Z0 = fill(3e3, 2), free)
        fit = maximum_likelihood(m; starts = 1, rounds = 2, iters = 800)
        mt = SurveyModel(data, ps; p, obs, Z0 = fill(1e4, 2), free)
        @test fit.loglik >= LDP.logdensity(mt, Chytridiomycosis.survey_start(mt)) - 1e-6
        @test 0.5 < fit.estimates.f / p.f < 2
        @test 0.5 < fit.estimates.η / p.η < 2
        prof = profile_likelihood(m, :f, fit.estimates.f .* [1.0, 3.0]; x0 = fit.x,
            rounds = 1, iters = 400)
        @test prof[1] ≈ fit.loglik atol = 0.5
        @test prof[2] < fit.loglik - 1.92                # f = 3 f̂ is excluded
    end
end
