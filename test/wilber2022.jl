# Wilber et al. 2022 (J Anim Ecol): outbreak-state HMM, decline regression, path analysis.
# The data-dependent tests need the authors' deposit (not bundled; data/wilber2022/README.md).
using Random: Xoshiro
using LogDensityProblems: LogDensityProblems
using Distributions: Beta as DBeta, Binomial as DBinomial, pdf as dpdf, cdf as dcdf

@testset "Wilber et al. 2022" begin
    @testset "transition matrix and emission weights" begin
        P = outbreak_transition_matrix(0.3, 0.1, 0.6, 0.2)
        @test all(sum(P; dims = 2) .≈ 1)
        @test P[1, :] ≈ [0.7, 0.18, 0.12, 0]
        @test P[2, :] ≈ [0.2, 0.48, 0.32, 0]
        @test P[3, :] == [0, 0, 0.9, 0.1] && P[4, :] == [0, 0, 0, 1]

        @test wilber2022_emission(-1, 0, 0, 0.25) == (1, 1, 1, 1)
        @test wilber2022_emission(3, 0, 0, 0.25) == (0, 0, 0, 1)
        n, y, ρ = 12, 4, 0.25
        e0 = wilber2022_emission(0, n, 0, ρ)
        @test e0[1] == 1 && e0[4] == 0
        @test e0[2] ≈ (1 - (1 - ρ)^(n + 1)) / (n + 1)
        @test e0[3] ≈ (1 - ρ)^(n + 1) / (n + 1)
        # S(a, b) by midpoint quadrature of P(K in class | n, p) × Beta(p; y + 1, n − y + 1)
        function S(lo, hi, a, b; m = 20_000)
            h = (b - a) / m
            ps = a .+ h .* ((1:m) .- 0.5)
            return h * sum(dpdf(DBeta(y + 1, n - y + 1), p) *
                           (dcdf(DBinomial(n, p), hi) - dcdf(DBinomial(n, p), lo - 1)) for p in ps)
        end
        eL, eH = wilber2022_emission(1, n, y, ρ), wilber2022_emission(2, n, y, ρ)
        @test eL[2] ≈ S(1, 3, 0, ρ) rtol = 1e-6                 # n ÷ ⌊1/ρ⌋ = 3
        @test eL[3] ≈ S(1, 3, ρ, 1) rtol = 1e-6
        @test eH[2] ≈ S(4, n, 0, ρ) rtol = 1e-6
        @test eH[3] ≈ S(4, n, ρ, 1) rtol = 1e-6
        # with K = 0 the classes exhaust the binomial: total mass one
        k0 = exp(Chytridiomycosis._lgamma(2n - y + 1) + Chytridiomycosis._lgamma(n + 2) -
                 Chytridiomycosis._lgamma(2n + 2) - Chytridiomycosis._lgamma(n - y + 1))
        @test eL[2] + eL[3] + eH[2] + eH[3] + k0 ≈ 1
        # the code's integer split differs from ⌊nρ⌋ at ρ = 0.3333
        @test wilber2022_emission(1, 3, 1, 0.3333) != wilber2022_emission(1, 3, 1, 0.3333; cutoff_rule = :exact)
        @test wilber2022_emission(1, 8, 2, 0.25) == wilber2022_emission(1, 8, 2, 0.25; cutoff_rule = :exact)
    end

    @testset "compendium" begin
        recs = parameter_records(:wilber2022)
        @test length(recs) == 73 && all(r -> r.verified && r.provenance == "F", recs)
        sets = unique(r.set for r in recs)
        @test all(s -> s in sets, ("abundance_rho0.25", "abundance_rho0.33", "abundance_rho0.5",
                                   "density_rho0.25", "density_rho0.33", "density_rho0.5"))
        p = parameter_set(:wilber2022; set = "abundance_rho0.25")
        @test p.decline_at_7_frogs == 0.84 && p.phi_tadpole == -3.03
    end

    if Chytridiomycosis._wilber2022_find() === nothing
        @info "Wilber et al. 2022 deposit not found; data tests skipped (data/wilber2022/README.md)"
        @test_skip false
    else
        @testset "data (79 lakes, 2004–2019)" begin
            d = wilber2022_data(cutoff = 0.25)
            @test size(d.obs) == (79, 16) && d.years == 2004:2019
            @test [count(==(c), d.obs) for c in -1:3] == [721, 307, 88, 140, 8]
            @test count(d.missing_abundance) == 762
            @test all(r -> count(!=(-1), r) >= 1, eachrow(d.obs))
            @test d.scaling.temperature.mean ≈ 14.34 atol = 0.01   # °C, 11–18 °C range in the text
            @test sum(d.abundance) ≈ 0 atol = 1e-8                 # z-scored
            # observed states are consistent with the swab counts
            ok = [d.obs[i] == -1 || d.obs[i] == 3 ||
                  (d.obs[i] == 0) == (d.num_inf[i] == 0) for i in eachindex(d.obs)]
            @test all(ok)
            d5 = wilber2022_data(cutoff = 0.5)
            @test count(==(2), d5.obs) < count(==(2), d.obs)        # fewer high-prevalence years
            @test_throws ArgumentError wilber2022_data(cutoff = 1.5)
        end

        @testset "HMM likelihood" begin
            prob = OutbreakHMMProblem(cutoff = 0.25)
            @test LogDensityProblems.dimension(prob) == 15
            θ = 0.2 .* randn(Xoshiro(3), 15)
            seqs = outbreak_sequences(prob, θ)
            generic = [hmm_loglikelihood(s.initial, s.transitions, s.likelihoods; emit_first = false)
                       for s in seqs]
            @test Chytridiomycosis.pointwise_loglikelihood(prob, θ) ≈ generic
            m = HMMLogDensity(θ -> outbreak_sequences(prob, θ),
                              θ -> Chytridiomycosis._w22_logprior(prob, θ), 15; emit_first = false)
            @test LogDensityProblems.logdensity(m, θ) ≈ LogDensityProblems.logdensity(prob, θ)
            r = outbreak_rates(prob, θ)
            @test size(r.omega) == (79, 15) && all(0 .< r.phi .< 1)
            p = wilber2022_parameters(prob, θ)
            @test p.a ≈ 1 / (1 + exp(-θ[1])) && p.sigma ≈ exp(θ[end])
            @test keys(p.omega) == (:intercept, :abundance, :temperature, :snow, :tadpole, :temp_snow)
            γ = outbreak_state_probabilities(prob, θ)
            @test all(sum(γ; dims = 3) .≈ 1)
            @test all(γ[:, 1, 1] .≈ 1)                               # every lake starts Bd-free
            v = outbreak_viterbi(prob, θ)
            @test size(v) == (79, 16) && all(v[:, 1] .== 0)
            @test all(v[prob.data.obs .== 3] .== 3)                 # observed E is certain
            @test all(v[prob.data.obs .== 1] .∈ Ref((1, 2)))        # Bd seen: LP or HP
            @test_throws ArgumentError OutbreakHMMProblem(prob.data; model = :log_tadpole)
            @test length(OutbreakHMMProblem(cutoff = 0.25, model = :linear_temp).names) == 11
            @test length(OutbreakHMMProblem(cutoff = 0.25, model = :quadratic_temp).names) == 13
        end

        @testset "posterior vs Fig. 3 (ρ = 0.25, abundance)" begin
            fit = wilber2022_fit(cutoff = 0.25)
            @test maximum(fit.rhat) < 1.05
            @test all(0.15 .< fit.acceptance .< 0.4)
            pub = parameter_set(:wilber2022; set = "abundance_rho0.25")
            for name in (:omega_abundance, :omega_temperature, :omega_snow, :omega_temp_snow,
                         :omega_tadpole, :phi_abundance, :phi_temperature, :phi_snow,
                         :phi_temp_snow)
                @test wilber2022_effect(fit, name).median ≈ getfield(pub, name) atol = 0.05
            end
            @test wilber2022_effect(fit, :phi_tadpole).median ≈ pub.phi_tadpole atol = 0.3
            s = posterior_summary(fit)
            @test s.names == fit.names && all(s.lower .< s.median .< s.upper)
            l = psis_loo(fit)
            @test length(l.pointwise) == 79 && looic(fit) ≈ -2 * l.elpd
            # the cached single-seed fit agrees with the stored four-seed mean (scripts/wilber2022_looic.jl)
            @test abs(looic(fit) - Chytridiomycosis._wilber2022_looic(:full, 0.25).mean) < 2.5
        end

        @testset "decline regression and path analysis" begin
            @test length(Wilber2022DeclineProblem().start) == 82
            @test length(Wilber2022DeclineProblem(subadults = false).start) == 52
            fit = wilber2022_decline()
            @test maximum(fit.rhat) < 1.02
            e = decline_effect(fit)
            @test e.lower > 0                                        # Fig. 5C: above zero
            pd = predicted_decline(fit, 7)
            @test pd.median ≈ 0.84 atol = 0.02
            pa = wilber2022_path_analysis()
            @test pa.n == 111 && pa.populations == 40
            @test pa.host_to_load.lower < 0 < pa.host_to_load.upper
        end
    end
end
