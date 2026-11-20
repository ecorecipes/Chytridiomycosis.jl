# Occupancy layer, Mosher et al. 2018 Ecol Appl 28:926 (doi:10.1002/eap.1699) and Adams et
# al. 2010 Ecol Appl 20:289 (doi:10.1890/08-2319.1). Oracles:
# - brute-force sums over latent states (single-season and dynamic likelihoods);
# - the closed-form MLE of the constant single-season model (MacKenzie et al. 2002);
# - Mosher (2017) PhD dissertation (the published Mosher et al. 2018 paper was not readable):
#   Ch. 1 Tables 1.1-1.2 (matrices), Table 1.3 (host occupancy declines 0.90 → 0.71 and
#   0.20 → 0.07 over 4 seasons), Ch. 2 Table 2.3 (parameter counts K of the listed models);
# - recovery of known parameters from simulated data;
# - Adams et al. 2010 Appendix A (Fig. A1 caption: the finite-population ratio → 1).
using Chytridiomycosis
using Test
using Random: Xoshiro

@testset "Occupancy models (Mosher 2018, Adams 2010)" begin
    C = Chytridiomycosis

    @testset "single-season likelihood" begin
        y = [1, 0, missing, 0]
        p = [0.3, 0.5, 0.9, 0.2]
        ψ = 0.6
        @test occupancy_loglikelihood(y, ψ, p) ≈ log(ψ * 0.3 * 0.5 * 0.8)
        @test occupancy_loglikelihood([0, 0, missing], ψ, 0.4) ≈ log(ψ * 0.36 + 0.4)
        @test occupancy_loglikelihood(Union{Int,Missing}[missing, missing], ψ, 0.4) ≈ 0 atol = 1e-14
        Y = [1 0 1; 0 0 0; 0 1 missing]
        @test occupancy_loglikelihood(Y, ψ, 0.4) ≈
              sum(occupancy_loglikelihood(Y[i, :], ψ, 0.4) for i in 1:3)
        @test occupancy_loglikelihood(Y, [0.2, 0.5, 0.7], [0.1 0.2 0.3; 0.4 0.5 0.6; 0.7 0.8 0.9]) ≈
              occupancy_loglikelihood(Y[1, :], 0.2, [0.1, 0.2, 0.3]) +
              occupancy_loglikelihood(Y[2, :], 0.5, [0.4, 0.5, 0.6]) +
              occupancy_loglikelihood(Y[3, :], 0.7, [0.7, 0.8, 0.9])
        @test_throws ArgumentError occupancy_loglikelihood([2, 0], 0.5, 0.5)

        # closed-form MLE with K surveys at every site: p̂/(1 − (1 − p̂)^K) = d/(n K),
        # ψ̂ = n/(N (1 − (1 − p̂)^K)), with n sites detected and d detections
        rng = Xoshiro(11)
        N, K = 400, 4
        occ = rand(rng, N) .< 0.6
        Y = [occ[i] && rand(rng) < 0.35 ? 1 : 0 for i in 1:N, j in 1:K]
        n = count(i -> any(==(1), Y[i, :]), 1:N)
        d = sum(Y)
        f = occupancy_mle(Y)
        @test f.p / (1 - (1 - f.p)^K) ≈ d / (n * K) rtol = 1e-5
        @test f.ψ ≈ n / (N * (1 - (1 - f.p)^K)) rtol = 1e-5
        @test f.naive ≈ n / N
        @test abs(f.ψ - 0.6) < 4f.se_ψ && abs(f.p - 0.35) < 4f.se_p
    end

    @testset "two-species matrices" begin
        p = Mosher2018Parameters(:mosher2018)
        Φ = transition_matrix(p)
        D = detection_matrix(p)
        @test sum(initial_state(p)) ≈ 1
        @test all(sum(Φ; dims = 2) .≈ 1)
        @test all(sum(D[:, 1:4]; dims = 2) .≈ 1)
        # Table 1.1, e.g. AB → B is εAB(1 − εBA); U → AB is γA γB
        @test Φ[4, 2] ≈ p.εAB * (1 - p.εBA)
        @test Φ[3, 4] ≈ p.γA * p.γB
        @test Φ[1, 4] ≈ (1 - p.εA) * p.γBA
        # Table 1.2 with pB = rBa = 0: Bd is never seen without a detected toad
        @test all(D[:, 2] .== 0) && D[2, 3] == 1
        @test D[4, 4] ≈ p.rA * p.rBA
        # partial records are sums of full columns
        @test D[:, 5] ≈ D[:, 1] .+ D[:, 4]
        @test D[:, 6] ≈ D[:, 2] .+ D[:, 3]
        for a in (0, 1, missing), b in (0, 1, missing)
            c = two_species_code(a, b)
            r = two_species_records(c)
            @test isequal(r, (a, b))
        end
        @test two_species_code(1, 1) == 4 && two_species_code(0, 0) == 3
        q = C.with(p; εAB = 0.5)
        @test q.εAB == 0.5 && q.γA == p.γA
    end

    @testset "dynamic likelihood against brute force" begin
        p = Mosher2018Parameters(:mosher2018; set = "rare_susceptible")
        p = C.with(p; pB = 0.3, rBa = 0.2)            # make every state observable
        φ0, Φ, D = initial_state(p), transition_matrix(p), detection_matrix(p)
        codes = Union{Missing,Int}[4 5 missing; 3 6 2; 1 8 7]
        brute = 0.0
        for s1 in 1:4, s2 in 1:4, s3 in 1:4
            s = (s1, s2, s3)
            pr = φ0[s1] * Φ[s1, s2] * Φ[s2, s3]
            for t in 1:3, j in 1:3
                codes[t, j] === missing || (pr *= D[s[t], codes[t, j]])
            end
            brute += pr
        end
        @test dynamic_occupancy_loglikelihood(φ0, Φ, D, codes) ≈ log(brute)
        # a survey-dependent emission function gives the same as a constant matrix
        @test dynamic_occupancy_loglikelihood(φ0, Φ, (t, j) -> D, codes) ≈ log(brute)
        # single-season occupancy is the K = 2 special case with one season
        Φ2 = [1.0 0.0; 0.0 1.0]
        E = [0.7 0.3; 0.0 1.0]                        # occupied, unoccupied; columns 1 = detected
        @test dynamic_occupancy_loglikelihood([0.4, 0.6], Φ2, E, Union{Missing,Int}[1 2 2]) ≈
              occupancy_loglikelihood([1, 0, 0], 0.4, 0.7)
    end

    @testset "dissertation Table 1.3 declines" begin
        # host occupancy falls from 0.90 to 0.71 (SRM) and 0.20 to 0.07 (rare, susceptible)
        s = solve(Mosher2018Problem(Mosher2018Parameters(:mosher2018); nseasons = 4))
        @test s.host[1] ≈ 0.90
        @test round(s.host[4]; digits = 2) == 0.71
        r = solve(Mosher2018Problem(Mosher2018Parameters(:mosher2018; set = "rare_susceptible");
                                    nseasons = 4))
        @test r.host[1] ≈ 0.20
        @test round(r.host[4]; digits = 2) == 0.07
        @test all(sum(s.states; dims = 2) .≈ 1)
        @test s.pathogen ≈ s.states[:, 2] .+ s.states[:, 4]
        @test_throws ArgumentError Mosher2018Problem(Mosher2018Parameters(:mosher2018); nseasons = 0)
    end

    @testset "model structures (dissertation Table 2.3 K)" begin
        rng = Xoshiro(3)
        d = simulate_mosher2018(rng, Mosher2018Parameters(:mosher2018); nsites = 20,
                                nseasons = 3, nsurveys = 2, elevation = 2500 .+ 1000rand(rng, 20))
        d = Mosher2018Data(; toad = d.toad, bd = d.bd, elevation = d.elevation,
                           doy = 150 .+ 60rand(rng, size(d.toad)...),
                           other_amphibians = rand(rng, 0:1, 20), permanent = rand(rng, 0:1, 20))
        det = (pA = :doy2, rA = :doy2, rBA = :elev)
        K(; kw...) = length(parameter_names(Mosher2018Model(d; det..., kw...)))
        @test K(γAB = :elev2, γB = :elev2, γBA = :elev2, εAB = :elev2, εBA = :elev2) == 25
        @test K(γAB = :elev, γB = :elev, γBA = :elev, εAB = :elev, εBA = :elev) == 21
        @test K(γB = :elev2, γBA = :elev2, εBA = :elev2) == 21
        @test K(γB = :elev, γBA = :elev, εBA = :elev) == 19
        @test K(γB = :elev, γBA = :elev) == 18
        @test K(εBA = :elev) == 18
        @test K(γB = :perm, γBA = :perm, εBA = :perm) == 19
        @test K(εBA = :perm) == 18
        @test K(γAB = :amph) == 18
        @test K(εBA = :elev2) == 19
        @test K() == 17
        # the printed null model has K = 14; ours has 13 (one constant per pair)
        @test length(parameter_names(Mosher2018Model(d; null = true, det...))) == 13
        @test_throws ArgumentError Mosher2018Model(d; γA = :doy)
        @test_throws ArgumentError Mosher2018Model(d; foo = :elev)

        # zero slopes give the constant model's likelihood
        mc = Mosher2018Model(d)
        mq = Mosher2018Model(d; det..., γB = :elev2, γBA = :elev2)
        θc = randn(Xoshiro(5), length(parameter_names(mc)))
        θq = zeros(length(parameter_names(mq)))
        for (k, nm) in enumerate(parameter_names(mq))
            j = findfirst(==(nm), parameter_names(mc))
            j === nothing || (θq[k] = θc[j])
        end
        @test mosher2018_loglikelihood(mq, θq) ≈ mosher2018_loglikelihood(mc, θc)
        @test C.LogDensityProblems.logdensity(mc, θc) == mosher2018_loglikelihood(mc, θc)
        sp = mosher2018_site_parameters(mq, θq, 1)
        @test sp.γA ≈ 1 / (1 + exp(-θc[4]))
    end

    @testset "recovery from simulated data" begin
        p = Mosher2018Parameters(:mosher2018)
        d = simulate_mosher2018(Xoshiro(2018), p; nsites = 1000, nseasons = 4, nsurveys = 3)
        m = Mosher2018Model(d)
        f = maximum_likelihood(m)
        @test f.converged
        @test f.k == 14 && f.n == 1000
        est = mosher2018_site_parameters(m, f.θ, 1)
        for s in (:ψA, :ψBA, :γBA, :εAB, :εBA, :pA, :rA, :rBA)
            @test abs(getproperty(est, s) - getproperty(p, s)) < 0.08
        end
        @test f.aicc ≈ -2f.loglik + 2 * 14 + 2 * 14 * 15 / (1000 - 15)
    end

    @testset "data location" begin
        if C._mosher2018_find() === nothing
            withenv("CHYTRIDIOMYCOSIS_MOSHER2018" => mktempdir()) do
                @test_throws ArgumentError mosher2018_path()
            end
        else
            @test isdir(mosher2018_path())
        end
    end

    @testset "Adams 2010 finite population" begin
        @test finite_population_occurrence(0.5, 0.2, 0) == 0
        @test finite_population_occurrence(0.5, 0.2, 3) ≈ 0.5 * (1 - 0.8^3)
        r(p′, N) = finite_population_occurrence(1.0, p′, N)
        @test r(0.1, 50) > r(0.1, 10) && r(0.3, 10) > r(0.1, 10)
        @test r(0.5, 50) ≈ 1 atol = 1e-12
        @test_throws DomainError finite_population_occurrence(1.2, 0.1, 3)
        q = parameter_set(:adams2010; set = "bullfrog")
        @test q.npop == 40 && 0 < q.ψ < 1
    end
end
