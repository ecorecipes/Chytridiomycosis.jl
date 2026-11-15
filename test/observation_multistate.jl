# Generic multistate capture–recapture layer (src/observation/multistate.jl): transition and
# emission builders, the multi-event likelihood against brute-force enumeration of latent
# paths, and recovery of known parameters from simulated histories.
using Random: Xoshiro
using LinearAlgebra: LinearAlgebra
using LogDensityProblems: LogDensityProblems

@testset "Multistate capture-recapture layer" begin
    @testset "continuous-time transitions" begin
        Q = infection_generator(0.3, 2.6, 0.09, 0.4)
        @test all(abs.(sum(Q; dims = 2)) .< 1e-12)
        for τ in (0.0, 0.5, 1.0, 4.6)
            P = ctmc_transition(Q, τ)
            @test P ≈ exp(Matrix(Q) .* τ) atol = 1e-12
            @test all(sum(P; dims = 2) .≈ 1)
        end
        @test ctmc_transition(Q, 0.7) * ctmc_transition(Q, 1.3) ≈ ctmc_transition(Q, 2.0)
        # equal diagonal (q → 0 branch) still exact
        Q0 = infection_generator(0.0, 0.0, 0.2, 0.2)
        @test ctmc_transition(Q0, 2.0) ≈ exp(Matrix(Q0) .* 2.0)
        # a general generator goes through LinearAlgebra.exp
        G = [-0.5 0.3 0.2 0.0; 0.1 -0.4 0.2 0.1; 0.0 0.3 -0.3 0.0; 0.0 0.0 0.0 0.0]
        @test ctmc_transition(G, 1.5) ≈ exp(G .* 1.5)
        @test_throws DomainError ctmc_transition(Q, -1.0)
    end

    @testset "builders" begin
        Φ = survival_transition([0.8, 0.6], [0.7 0.3; 0.4 0.6])
        @test size(Φ) == (3, 3) && all(sum(Φ; dims = 2) .≈ 1)
        @test Φ[1, 2] ≈ 0.8 * 0.3 && Φ[2, 3] ≈ 0.4
        B = multievent_emission([0.5, 0.4], [0.2 0.0 0.8; 0.0 0.9 0.1])
        @test size(B) == (3, 4) && all(sum(B; dims = 2) .≈ 1)
        @test B[1, 2] ≈ 0.5 * 0.2 && B[3, 1] == 1
        codes = Union{Missing,Int}[1 3 missing; missing missing missing; 2 2 1]
        L = robust_design_likelihoods(B, codes)
        @test L[1, :] ≈ B[:, 1] .* B[:, 3]
        @test L[2, :] == ones(3)
        @test L[3, :] ≈ B[:, 2] .^ 2 .* B[:, 1]
        f(t, j) = t == 1 ? B : 2B
        @test robust_design_likelihoods(f, codes)[3, :] ≈ 8 .* L[3, :]
        q = quasistationary_distribution(Φ)
        M = Φ[1:2, 1:2]
        v = [0.5, 0.5]
        for _ in 1:500
            v = vec(v' * M)
            v ./= sum(v)
        end
        @test q ≈ v
    end

    st = MultiEventStructure(; nstates = 2, assignment = Bool[1 0 1; 0 1 1], baseline = [3, 3],
                             survival = :state, transition = :state, detection = :time,
                             state_assignment = :state)
    rng = Xoshiro(11)
    H = [0 1 0 3 2; 3 3 0 1 missing; 0 0 2 0 0; 1 2 3 0 1; 0 0 0 0 1; 2 0 0 0 0]
    prob = MultiEventProblem(H, st)
    @test size(prob.histories, 1) == 5           # never seen after first on the last occasion dropped
    θ = 0.4 .* randn(rng, LogDensityProblems.dimension(prob))
    par = multievent_parameters(prob, θ)
    @test sum(par.initial) ≈ 1
    @test all(sum(par.transition; dims = 2) .≈ 1)
    @test all(sum(par.assignment; dims = 2) .≈ 1)
    @test par.assignment[1, 2] == 0 && par.assignment[2, 1] == 0

    @testset "likelihood = brute-force sum over latent paths" begin
        T = size(H, 2)
        Φ = [survival_transition(par.survival[:, t], par.transition[:, :, t]) for t in 1:T-1]
        B = [multievent_emission(par.detection[:, t-1], par.assignment) for t in 2:T]
        for r in axes(prob.histories, 1)
            h = prob.histories[r, :]
            f = prob.first[r]
            tot = 0.0
            for path in Iterators.product(ntuple(_ -> 1:3, T - f + 1)...)
                z = collect(path)
                z[1] == 3 && continue
                p = par.initial[z[1]] * par.assignment[z[1], h[f]]
                for t in f+1:T
                    p *= Φ[t-1][z[t-f], z[t-f+1]]
                    h[t] === missing || (p *= B[t-1][z[t-f+1], h[t]+1])
                end
                tot += p
            end
            @test pointwise_loglikelihood(prob, θ)[r] ≈ log(tot)
        end
        @test LogDensityProblems.logdensity(prob, θ) ≈
              sum(prob.counts .* pointwise_loglikelihood(prob, θ))
    end

    @testset "simulation and recovery" begin
        truth = (initial = [0.3, 0.7], survival = repeat([0.65, 0.85], 1, 5),
                 transition = repeat([0.6 0.4; 0.2 0.8], 1, 1, 5),
                 detection = fill(0.6, 2, 5), assignment = [0.5 0.0 0.5; 0.0 0.5 0.5])
        Hs = simulate_capture_histories(Xoshiro(3), st, truth, fill(400, 6))
        @test size(Hs) == (2400, 6)
        @test all(i -> Hs[i, findfirst(>(0), Hs[i, :])] in 1:3, axes(Hs, 1))
        cst = MultiEventStructure(; nstates = 2, assignment = Bool[1 0 1; 0 1 1], baseline = [3, 3],
                                  survival = :state, transition = :state, detection = :constant,
                                  state_assignment = :state)
        fit = solve(MultiEventProblem(Hs, cst))
        @test maximum(abs, fit.gradient) < 1e-2
        @test fit.rank_deficiency == 0
        se = sqrt.(LinearAlgebra.diag(fit.vcov))
        @test all(isfinite, se)
        est = fit.parameters
        @test est.survival[:, 1] ≈ [0.65, 0.85] atol = 0.06
        @test est.transition[1, 2, 1] ≈ 0.4 atol = 0.08
        @test est.transition[2, 1, 1] ≈ 0.2 atol = 0.06
        @test est.detection[1, 1] ≈ 0.6 atol = 0.05
        @test est.assignment[:, 3] ≈ [0.5, 0.5] atol = 0.1
        @test fit.aic ≈ fit.deviance + 2fit.k
        @test qaic(fit, 1.5) ≈ fit.deviance / 1.5 + 2fit.k
    end
end
