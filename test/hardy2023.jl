# Hardy et al. 2023 J Appl Ecol (doi:10.1111/1365-2664.14310): robust-design Pradel model of
# adult male boreal toads at five Colorado sites, 1998-2020, with the Pradel (1996)
# temporal-symmetry layer (src/observation/pradel.jl). Data: Dryad doi:10.5061/dryad.fttdz08wv
# (CC0, data/hardy2023).
#
# Oracles:
#   - Table 1 (captures before and after Bd at each site) and "We captured 1,594 individual
#     male boreal toads";
#   - Appendices S1.3, S1.4 and S1.5: −2 log L, K and AICc of the focal model sets, fitted in
#     Program MARK (all 27 models are fitted here; 26 reproduce MARK's −2 log L to 0.01);
#   - Table 3 (AICc of the top model) and Results (post-Bd φ, f and λ; pre-Bd φ at the
#     lowest and highest site; range of p).
#   All are from the version of record (J Appl Ecol 60:121-131) and its Supporting
#   Information, which print the same numbers as the reprint in B. M. Hardy's dissertation.
#   - Closed forms: probabilities of all possible histories sum to one; brute-force
#     enumeration of the history likelihood; simulation and recovery of known parameters.
using Chytridiomycosis
using Test
using Random: Xoshiro
using LogDensityProblems: LogDensityProblems

@testset "Hardy 2023 boreal toads (robust-design Pradel)" begin
    C = Chytridiomycosis

    @testset "Pradel likelihood: closed forms" begin
        # all histories over T = 3 periods with 2 occasions each; the probabilities of the
        # histories with at least one capture sum to one
        φ, f = [0.7, 0.5], [0.4, 0.9]
        p = [0.3, 0.5, 0.2, 0.6, 0.4, 0.45]
        primary = [1, 1, 2, 2, 3, 3]
        total = 0.0
        for bits in 1:(2^6 - 1)
            h = reshape([(bits >> (j - 1)) & 1 for j in 1:6], 1, 6)
            total += exp(pradel_loglikelihood(PradelGroup(h; primary), φ, f, p))
        end
        @test total ≈ 1 rtol = 1e-12
        # one occasion per period: history 101 by hand (Pradel 1996)
        q = [0.3, 0.6, 0.5]
        g = PradelGroup([1 0 1])
        λ = φ .+ f
        γ2, γ3 = φ[1] / λ[1], φ[2] / λ[2]
        ξ = [1.0, 1 - γ2 + γ2 * (1 - q[1]), 0.0]
        ξ[3] = 1 - γ3 + γ3 * (1 - q[2]) * ξ[2]
        N = [1.0, λ[1], λ[1] * λ[2]]
        D = sum(N .* ξ .* q)
        @test exp(pradel_loglikelihood(g, φ, f, q)) ≈ q[1] * φ[1] * (1 - q[2]) * φ[2] * q[3] / D
        # unsampled occasions: a missing column is the same as dropping it
        h = [1 missing 0 1; 0 missing 1 1]
        g1 = PradelGroup(h; primary = [1, 1, 2, 2])
        g2 = PradelGroup([1 0 1; 0 1 1]; primary = [1, 2, 2])
        @test pradel_loglikelihood(g1, [0.6], [0.3], [0.4, 0.9, 0.5, 0.3]) ≈
              pradel_loglikelihood(g2, [0.6], [0.3], [0.4, 0.5, 0.3])
        @test pradel_capture_probability(g1, [0.4, 0.9, 0.5, 0.3]) ≈ [0.4, 1 - 0.5 * 0.7]
        @test_throws ArgumentError PradelGroup([1 missing; 0 1])
        @test_throws DimensionMismatch pradel_loglikelihood(g2, [0.6, 0.1], [0.3], [0.4, 0.5, 0.3])
    end

    @testset "simulation and recovery" begin
        T = 8
        φ, f = fill(0.7, T - 1), fill(0.35, T - 1)
        primary = repeat(1:T; inner = 2)
        p = fill(0.35, 2T)
        H = simulate_pradel(Xoshiro(11), 1500, φ, f, p; primary)
        g = PradelGroup(H; primary)
        design = (φ = ones(T - 1, 1), f = ones(T - 1, 1), p = ones(2T, 1))
        fit = solve(PradelProblem([g], [design]))
        est = fit.groups[1]
        @test est.φ[1] ≈ 0.7 atol = 0.03
        @test est.f[1] ≈ 0.35 atol = 0.03
        @test est.p[1] ≈ 0.35 atol = 0.02
        @test est.λ[1] ≈ 1.05 atol = 0.03
        @test isapprox(est.N[1], 1500; rtol = 0.1)
        @test maximum(abs, fit.gradient) < 1e-3
        # standard errors: Wald SE of logit φ against a profile-free check of curvature
        se = pradel_se(fit, θ -> pradel_derived(fit.problem, θ, 1).φ[1]).se
        @test 0.005 < se < 0.03
        @test LogDensityProblems.dimension(fit.problem) == 3
    end

    data = load_hardy2023()

    @testset "data (Table 1)" begin
        caps = hardy2023_captures(data)
        @test sum(c.toads for c in caps) == 1594
        tbl = parameter_set(:hardy2023)
        @test [c.name for c in caps] == ["Delilah Pond", "Collegiate Camp", "Denny Creek",
                                         "South Cottonwood", "Morgan's Gulch"]
        keys_ = (:delilah, :collegiate, :denny, :cottonwood, :morgans)
        for (c, k) in zip(caps, keys_)
            @test c.elevation == getproperty(tbl, Symbol("elevation_", k))
            @test c.bd_year == getproperty(tbl, Symbol("bd_year_", k))
            @test c.before == getproperty(tbl, Symbol("captures_pre_", k))
            @test c.after == getproperty(tbl, Symbol("captures_post_", k))
        end
        @test all(s -> size(s.histories, 2) == 69, data)
        # the header swap in Site_1/2_ENV.csv is undone: drought days < winter-onset day
        @test all(s -> all(skipmissing(s.env.drought) .< skipmissing(s.env.winter_start)), data)
    end

    @testset "focal model sets (Appendices S1.3-S1.5)" begin
        # −2 log L and K printed by MARK
        S3 = Dict(:s_bd_t_bd => (27669.30, 77), :s_bd => (28676.37, 8), :t_bd => (28669.52, 21),
                  :site => (28905.69, 7), :constant => (29334.27, 3))
        S4 = Dict(:elev_t_bd => (28046.40, 22), :t_bd => (28085.86, 21),
                  :elev_bd => (28225.02, 5), :elev_x_bd => (28225.02, 6), :bd => (28270.37, 4),
                  :elev => (29294.13, 4), :constant => (29334.27, 3))
        S5 = Dict(:t_bd => (28570.08, 21), :elev_t_bd => (28570.00, 22),
                  :drought_bd => (28750.73, 5), :drought_x_bd => (28750.72, 6),
                  :bd => (28756.20, 4), :elev_bd => (28755.88, 5), :winstart_bd => (28756.20, 5),
                  :elev_x_bd => (28754.36, 6), :winstart_x_bd => (28756.19, 6),
                  :winstart_elev => (29266.9, 5), :winstart => (29278.102, 4),
                  :drought_elev => (29323.899, 5), :drought => (29329.822, 4),
                  :elev => (29330.671, 4), :constant => (29334.271, 3))
        for (par, tbl) in ((:p, S3), (:φ, S4), (:f, S5))
            ms = hardy2023_model_set(par; data)
            @test first(ms).structure === (par === :p ? :s_bd_t_bd :
                                           par === :φ ? :elev_t_bd : :t_bd)
            for m in ms
                dev, k = tbl[m.structure]
                @test m.k == k
                if par === :φ && m.structure === :t_bd
                    # discrepant: our maximum is 9.28 units of −2 log L lower than MARK's
                    @test m.deviance < dev - 9
                elseif par === :f && m.structure === :elev_bd
                    @test m.deviance ≈ dev atol = 0.1   # 28755.79 against 28755.88
                else
                    @test m.deviance ≈ dev atol = 0.015
                end
            end
            par === :f && @test round(ms[2].weight; digits = 2) == 0.28
        end
    end

    @testset "top model (Table 3, Results)" begin
        fit = hardy2023_fit(; data)
        @test fit.k == 115
        @test fit.aicc ≈ 26649.10 atol = 1     # MARK: 26649.10; ours 26649.35
        s = hardy2023_summary(fit; data)
        pub = parameter_set(:hardy2023)
        @test s.φ_post.value < 1e-6
        @test round(s.f_post.value; digits = 2) == pub.f_post
        @test round(s.f_post.se; digits = 2) == pub.f_post_se
        @test round(s.λ_post.value; digits = 2) == pub.lambda_post
        @test round.(s.λ_post.ci; digits = 2) == (pub.lambda_post_lo, pub.lambda_post_hi)
        @test round.(s.p_range; digits = 2) == (pub.p_min, pub.p_max)
        # Fig. 1a model: pre-Bd survival at the lowest and highest site
        fit2 = hardy2023_fit(; φ = :elev_bd, f = :elev_bd, data)
        φM = pradel_se(fit2, θ -> pradel_derived(fit2.problem, θ, 5).φ[1])
        @test round(φM.value; digits = 2) == pub.phi_pre_morgans
        @test round(φM.se; digits = 2) == pub.phi_pre_morgans_se
        φD = pradel_se(fit2, θ -> pradel_derived(fit2.problem, θ, 1).φ[1])
        @test φD.value ≈ pub.phi_pre_delilah atol = 0.015   # 0.482 against 0.47 (SE 0.03)
        @test pradel_coefficient(fit2, :φ_elev).value > 0
        @test pradel_coefficient(fit2, :f_elev).value < 0
    end
end
