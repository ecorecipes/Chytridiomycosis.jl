# Chowdhury et al. 2025, Ecol Model 501:111028, doi:10.1016/j.ecolmodel.2025.111028.
# Checked against the PDF (§§2-6, Table 2, Table 3, Figs 1, 4, 5, 7, 9) and the online
# supplement (lab/papers/chowdhury2025_si/mmc1.pdf, Appendices A-H). Figure values below were
# read from the figures (pixel positions of the curves); see the compendium notes.
using Chytridiomycosis
using LinearAlgebra: eigvals
using Test

@testset "Chowdhury 2025 stage- and temperature-structured Bsal" begin
    C = Chytridiomycosis
    scen = chowdhury_scenarios()
    pars(T; kw...) = [ChowdhuryParameters(:chowdhury2025; temperature = T, K = s.adults, kw...)
                      for s in scen]
    sims(T; λ = nothing) = map(scen) do s
        p = ChowdhuryParameters(:chowdhury2025; temperature = T, K = s.adults)
        λ === nothing || (p = C.with(p; λ))
        prob = BsalStageTemperatureProblem(p; adults = s.adults, efts = s.efts)
        sol = solve(prob; saveat = 1.0)
        (; prob, sol, tot = [stage_totals(u) for u in sol.u])
    end
    peak(r, f) = maximum(getfield.(r.tot, f))
    peakday(r, f) = argmax(getfield.(r.tot, f)) - 1

    @testset "compendium and parameters" begin
        recs = parameter_records(:chowdhury2025)
        @test all(r -> r.provenance in ("M", "F", "A", "D"), recs)
        # only the values not printed in the paper (derived here) are unverified
        @test all(r -> r.verified || r.provenance == "D", recs)
        @test Set(r.symbol for r in recs if !r.verified) == Set([:χ_a, :c_a, :c_j])
        p6, p14 = ChowdhuryParameters(:chowdhury2025; temperature = 6),
                  ChowdhuryParameters(:chowdhury2025; temperature = 14)
        @test time_unit(p6) === :day
        @test p6.η_a == 1 / (3 * 365) && p6.η_j == 1 / 105 && p6.δ_a == 1 / (10.5 * 365)
        @test (p6.κ_a, p14.κ_a) == (5443.0, 2895.0)               # review: 5,443 and 2,895
        @test (p6.ω_a, p14.ω_a) == (67.91, 177.47)                # review: 67.9 and 177.5
        @test (p6.P_a, p14.P_a) == (0.19, 0.34)
        @test (p6.c_a, p14.c_a) == (464.77, 878.48)
        @test (p6.ξ_a, p14.ξ_a) == (1 / 20, 1 / 25)               # Table 2, as the review says
        @test p6.ω_j == p6.ω_a && p14.P_j == p14.P_a              # "assumed same as adult"
        # Table 3 variant: per capita eft contacts and b_j = 0.005 (§6.2)
        t6 = ChowdhuryParameters(:chowdhury2025; temperature = 6, variant = :table3)
        @test t6.c_j ≈ 8 / (14 * 12) && t6.b_j == 0.005 && t6.c_a == 729.2
        @test parameter_set(:chowdhury2025).c_j ≈ 8 / 14 atol = 1e-3
        # text values that contradict Table 2
        a6 = ChowdhuryParameters(:chowdhury2025; temperature = 6, as_printed = true)
        a14 = ChowdhuryParameters(:chowdhury2025; temperature = 14, as_printed = true)
        @test (a6.ξ_j, a14.ξ_j) == (1 / 25, 1 / 20) && a6.η_a == 1 / (3.5 * 365)
        @test ChowdhuryParameters(:chowdhury2025; temperature = 14, b_a = 0.2).b_a == 0.2
        @test_throws ArgumentError ChowdhuryParameters(:chowdhury2025; temperature = 10)
        @test_throws ArgumentError ChowdhuryParameters(:chowdhury2025; variant = :fig9)
        @test_throws ArgumentError C.with(p6; P_a = 1.5)
        @test [s.adults for s in scen] == [2, 4, 8, 16]
        @test [s.efts for s in scen] == [0.1, 0.2, 0.4, 0.8]
    end

    @testset "model structure" begin
        p = ChowdhuryParameters(:chowdhury2025; temperature = 14, K = 4.0)
        d = stage_disease_free_state(p)
        # Appendix F: the disease-free equilibrium is an equilibrium
        prob = BsalStageTemperatureProblem(p; u0 = [d.S_l, d.S_j, 0, 0, 0, d.S_a, 0, 0, 0, 0, 0])
        du = similar(prob.u0)
        C.chowdhury_rhs!(du, prob.u0, prob, 0.0)
        @test maximum(abs, du) < 1e-10
        @test d.S_a ≈ 0.982 * 4 rtol = 1e-3
        @test d.S_j / d.S_a ≈ 16.2 rtol = 1e-2
        # Appendix F: eigenvalues of the life-cycle Jacobian (F.2) are about -0.0009, -0.0096, -0.02
        J = [-p.η_j 0 p.λ; p.α*p.η_j -(p.δ_j + p.η_a) 0;
             0 p.η_a*(1 - d.S_a / p.K) -p.η_a*d.S_j/p.K-p.δ_a]
        ev = sort(real.(eigvals(J)))
        @test ev[2:3] ≈ [-0.0096, -0.0009] atol = 1e-4
        # the third is -0.0147: the trace of (F.2) is -0.0253, against -0.0305 for the
        # printed eigenvalues (the printed Routh-Hurwitz coefficients do not match (F.2) either)
        @test_broken ev[1] ≈ -0.02 atol = 1e-3
        # the default direct terms are eqs (3a), (3b) and the paper's N_A = A + S_l
        u = [50.0, 0.3, 0.1, 0.2, 0.0, 3.0, 0.5, 0.4, 0.0, 200.0, 1500.0]
        tr = stage_transmission_rates(prob, u)
        A, NA, NT = 3.9, 53.9, 0.6
        βa = p.P_a * p.c_a * A / (p.χ_a + A)
        @test tr.β_a ≈ βa && tr.β_j ≈ p.c_j * p.P_j
        @test tr.direct_a ≈ βa * 3.0 * 0.4 / NA
        @test tr.direct_j ≈ p.c_j * p.P_j * 0.3 * 0.2 / NT
        @test tr.f_a ≈ p.ρ_a * p.b_a * 1500 / (p.κ_a + 1500)
        @test tr.f_j ≈ p.ρ_j * p.b_j * 200 / (p.κ_j + 200)
        # pluggable transmission: without larval dilution the default is HollingContact itself
        nd = BsalStageTemperatureProblem(p; larval_dilution = false)
        @test stage_transmission_rates(nd, u).direct_a ≈
              infection_rate(HollingContact(p.P_a, p.c_a, p.χ_a), 3.0, 0.4, A)
        fd = BsalStageTemperatureProblem(p; adult_contact = FrequencyDependent(0.23))
        @test stage_transmission_rates(fd, u).direct_a ≈ 0.23 * 3.0 * 0.4 / NA
        @test to_ode_problem(fd).p === fd
        # starting state of §5.1
        s = BsalStageTemperatureProblem(p; adults = 4, efts = 0.2)
        @test s.u0 ≈ [0, 0.198, 0.002, 0, 0, 3.96, 0.04, 0, 0, 0, 0]
        @test_throws ArgumentError BsalStageTemperatureProblem(p; u0 = zeros(3))
        @test_throws ArgumentError stage_decline_time(s; stage = :larvae)
    end

    @testset "reproduction numbers (eqs 5a, 5b; Table 3)" begin
        # the numerical next-generation matrix equals the printed formulas
        for p in vcat(pars(6), pars(14), pars(6; variant = :table3))
            d = stage_disease_free_state(p)
            Dj = p.δ_j + p.η_a + p.γ_j + p.μ_j
            Da = p.μ_a + p.δ_a + p.γ_a
            Rj = p.c_j * p.P_j * p.σ_j / (Dj * (p.η_a + p.δ_j + p.σ_j)) +
                 p.ρ_j * d.S_j * p.ω_j * p.σ_j * p.b_j / (p.κ_j * p.ξ_j * Dj * (p.η_a + p.δ_j + p.σ_j))
            Ra = d.S_a^2 * p.c_a * p.P_a * p.σ_a / ((p.σ_a + p.δ_a) * Da * (d.S_l + d.S_a) * (p.χ_a + d.S_a)) +
                 p.ρ_a * p.ω_a * p.b_a * p.σ_a * d.S_a / (p.κ_a * p.ξ_a * (p.σ_a + p.δ_a) * Da)
            r = stage_reproduction_numbers(p)
            @test r.R_j ≈ Rj rtol = 1e-5
            @test r.R_a ≈ Ra rtol = 1e-5
            @test stage_next_generation_matrix(p).R0 ≈ max(Rj, Ra) rtol = 1e-5
            @test basic_reproduction_number(p) == r.R0
        end
        tab3 = Dict(6 => (R_j = [2.35, 4.33, 8.29, 16.21], R_a = [4.82, 9.1, 16.5, 28.1]),
                    14 => (R_j = [1.26, 2.22, 4.15, 8.01], R_a = [16.55, 27.85, 43.72, 65.67]))
        R(T, v) = [stage_reproduction_numbers(p) for p in pars(T; variant = v)]
        # 14 °C: R_a with Table 2 (χ_a = 5.9); R_j with the Table 3 eft values
        @test [r.R_a for r in R(14, :table2)] ≈ tab3[14].R_a rtol = 2e-3
        @test [r.R_j for r in R(14, :table3)] ≈ tab3[14].R_j rtol = 4e-3
        # 6 °C: R_a with the Fig. 1 contact curve (variant :table3)
        @test [r.R_a for r in R(6, :table3)] ≈ tab3[6].R_a rtol = 4e-3
        # 6 °C R_a with Table 2's c_a = 464.77 and χ_a = 12 (the values of Figs 4, 5, 7)
        @test_broken [r.R_a for r in R(6, :table2)] ≈ tab3[6].R_a rtol = 0.02
        # 6 °C R_j: the direct part (the intercept, 0.37) reproduces; the environmental part
        # (slope 0.99 per unit K) is ~3x the formula's 0.33 with any printed b_j, ξ_j
        r6 = R(6, :table3)
        @test r6[1].direct_j ≈ 0.37 atol = 0.005
        @test [r.environmental_j / s.adults for (r, s) in zip(r6, scen)] ≈ fill(0.33, 4) atol = 0.005
        @test_broken [r.R_j for r in r6] ≈ tab3[6].R_j rtol = 0.02
        # Table 2's c_j = 0.571 gives a direct part of 4.43 (6 °C), not 0.37
        @test R(6, :table2)[1].direct_j ≈ 4.43 atol = 0.01
        # §6.2 and Fig. 9: R_a > 1 at 14 °C on the whole grid c_a ∈ [0, 600], K ∈ [1, 16]
        grid14 = [stage_reproduction_numbers(ChowdhuryParameters(:chowdhury2025; temperature = 14,
                    K, c_a)).R_a for c_a in 0:50:600, K in 1:16]
        @test minimum(grid14) > 1
        @test maximum(grid14) ≈ 53.7 atol = 1.5                         # Fig. 9: up to ~52
        # Fig. 9 at 6 °C reaches ~16-17 at c_a = 300, K = 16 with the Fig. 1 curve's χ_a
        @test stage_reproduction_numbers(ChowdhuryParameters(:chowdhury2025; temperature = 6,
                  K = 16, variant = :table3, c_a = 300.0)).R_a ≈ 16.8 atol = 0.3
        # R_j < 1 is attainable for efts (§6.2): low contact, fast spore decay
        @test stage_reproduction_numbers(ChowdhuryParameters(:chowdhury2025; temperature = 14,
                  K = 2, variant = :table3, ξ_j = 0.2)).R_j < 1
    end

    @testset "elasticities: signs stated in §6.1" begin
        p = ChowdhuryParameters(:chowdhury2025; temperature = 6, K = 2.0)
        ej = stage_elasticities(p; number = :R_j)
        ea = stage_elasticities(p; number = :R_a)
        @test all(k -> ej[k] > 0, (:α, :K, :ρ_j, :ω_j, :b_j, :λ))
        @test all(k -> ej[k] < 0, (:μ_j, :κ_j, :ξ_j, :η_a))
        @test all(k -> ea[k] > 0, (:K, :η_j, :ρ_a, :c_a, :P_a, :ω_a, :b_a))
        @test all(k -> ea[k] < 0, (:μ_a, :κ_a, :λ, :χ_a, :ξ_a))
        # direct and environmental shares sum to one; recovery is zero
        @test ej[:c_j] + ej[:ω_j] ≈ 1 atol = 1e-6
        @test ea[:P_a] + ea[:b_a] ≈ 1 atol = 1e-6
        @test ej[:γ_j] == 0
        @test_throws ArgumentError stage_elasticities(p; number = :R0)
    end

    s6, s14 = sims(6), sims(14)

    @testset "adult dynamics (Figs 4, 5, 7)" begin
        # Fig. 4 (e-h): infected adults peak at about half the density
        @test [peak(r, :infected_adults) for r in s6] ≈ [1.0, 2.05, 4.15, 8.3] rtol = 0.06
        @test [peak(r, :infected_adults) for r in s14] ≈ [0.95, 1.9, 3.8, 7.5] rtol = 0.06
        # Fig. 5 (c, d): aquatic zoospore peaks, read from the figure
        za6 = [peak(r, :Z_a) for r in s6]
        za14 = [peak(r, :Z_a) for r in s14]
        @test za6 ≈ [1000, 2050, 4200, 8250] rtol = 0.06
        @test za14 ≈ [2150, 4350, 8700, 17450] rtol = 0.06
        # "peak between 40-60 days" (§5.1); earlier at 14 °C (~32 days) than at 6 °C
        @test all(r -> 40 <= peakday(r, :Z_a) <= 60, s6)
        @test all(r -> 25 <= peakday(r, :Z_a) <= 40, s14)
        # ID50 exceeded by adults only at the highest density at 6 °C, and at 3 of 4 at 14 °C
        κ6, κ14 = s6[1].prob.p.κ_a, s14[1].prob.p.κ_a
        @test (za6 .> κ6) == [false, false, false, true]
        @test (za14 .> κ14) == [false, true, true, true]
        # Fig. 7: β_a at the start
        @test [stage_transmission_rates(r.prob, r.prob.u0).β_a for r in s6] ≈ [12.2, 21.7, 35.2, 50.5] rtol = 0.03
        @test [stage_transmission_rates(r.prob, r.prob.u0).β_a for r in s14] ≈ [74.5, 119.9, 170.7, 216.4] rtol = 0.02
        # Fig. 7 (bottom): peaks of the per capita infection rate from aquatic zoospores
        fa(r) = maximum(stage_transmission_rates(r.prob, u).f_a for u in r.sol.u)
        @test fa.(s6) ≈ [0.012, 0.021, 0.033, 0.045] rtol = 0.05
        @test fa.(s14) ≈ [0.033, 0.045, 0.057, 0.065] rtol = 0.05
        # adults near non-recovery levels within 120 days (§5.1): < 10% left
        @test all(r -> r.tot[121].adults < 0.1 * r.tot[1].adults, s14)
        @test all(r -> r.tot[121].adults < 0.1 * r.tot[1].adults, s6)
        # "infection and decline rates were about 2X slower at 6 °C" (§5.1)
        t6 = [stage_decline_time(r.prob) for r in s6]
        t14 = [stage_decline_time(r.prob) for r in s14]
        @test all(1.5 .< t6 ./ t14 .< 3.0)
        # §7: at 14 °C, 16 adults m⁻² reach an 80% decline 33% faster than 2 m⁻²
        @test_broken 1 - t14[4] / t14[1] ≈ 0.33 atol = 0.05
    end

    @testset "eft dynamics (Figs 4, 5, E.1)" begin
        # With Table 2's λ = 0.82 eggs adult⁻¹ day⁻¹, larvae recruit several efts per m²:
        # the eft totals rise several-fold instead of staying flat and falling (Fig. 4 a-b).
        @test all(r -> peak(r, :efts) > 3 * r.tot[1].efts, s6)
        @test_broken all(r -> peak(r, :efts) < 1.1 * r.tot[1].efts, s6)
        # Fig. 4 (c, d): infected efts peak at 0.38 (6 °C, day ~85) and 0.24 (14 °C, day ~70)
        @test_broken peak(s6[4], :infected_efts) ≈ 0.38 rtol = 0.15
        @test_broken peak(s14[4], :infected_efts) ≈ 0.243 rtol = 0.15
        # Without larval recruitment (λ = 0) the 14 °C eft curves match the figures
        n14, n6 = sims(14; λ = 0.0), sims(6; λ = 0.0)
        @test peak(n14[4], :infected_efts) ≈ 0.243 rtol = 0.06
        @test 65 <= peakday(n14[4], :infected_efts) <= 85
        @test peak(n14[4], :Z_j) ≈ 780 rtol = 0.06                 # Fig. 5b
        @test all(r -> peak(r, :efts) ≈ r.tot[1].efts, n14)
        # and at 6 °C they are slower and lower than Fig. 4c (0.38 at day ~85)
        @test_broken peak(n6[4], :infected_efts) ≈ 0.38 rtol = 0.15
        # the adult side does not depend on larval recruitment over one year
        @test [peak(r, :infected_adults) for r in n14] ≈ [peak(r, :infected_adults) for r in s14] rtol = 0.01
        # Fig. 7 (middle): f_j peaks at 0.011-0.029 (6 °C); ten times what either version gives
        fj(r) = maximum(stage_transmission_rates(r.prob, u).f_j for u in r.sol.u)
        @test fj(s6[4]) ≈ 0.0115 rtol = 0.02
        @test_broken fj(s6[4]) ≈ 0.029 rtol = 0.15
        # "adults declined about 3X faster than efts" (§5.1), only without recruitment
        for n in (n6, n14)
            ratio = [stage_decline_time(r.prob; stage = :efts) / stage_decline_time(r.prob)
                     for r in n]
            @test all(2 .< ratio .< 3.5)
        end
        # with Table 2's recruitment the 6 °C efts never fall by 80% within the year
        @test all(isinf, [stage_decline_time(r.prob; stage = :efts) for r in s6])
        # ID50 not exceeded for efts (§5.1)
        @test all(r -> peak(r, :Z_j) < r.prob.p.κ_j, vcat(n6, n14))
    end
end
