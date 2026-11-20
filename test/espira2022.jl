# Espira et al. 2022 Am Nat 199:E43 (doi:10.1086/717413; published version, author manuscript
# PMC9136953 and online supplement): two hosts, environmental transmission, partially
# overlapping habitats. Oracles: eqs 1-11 of the paper and S1-S18 of the supplement (closed
# forms; S12 is the next-generation matrix and S15 the community R0), a numerical
# next-generation matrix of eq. 1, the stability of the disease-free state of eq. 1 (R0,C > 1
# iff it is unstable), the epidemic ODE, and the worked example of Table S1 (compendium). The
# figures use no published values; the parameters below are illustrative. The paper's
# statements, Figs 3-4 and S1-S2 are validation targets.
using Chytridiomycosis
using Test
using LinearAlgebra: eigvals
using Random: Xoshiro

@testset "Espira 2022 overlapping habitats" begin
    C = Chytridiomycosis
    p = EspiraParameters(; N = (100, 50), γ = (0.1, 0.2), κ = (2.0, 1.0), ρ = (0.01, 0.05),
        π = (0.02, 0.005), α = (5.0, 3.0), ε = (2.0, 1.0), θ = (0.25, 0.5), μ = 0.5)
    @test C.time_unit(EspiraParameters) === :day

    @testset "single species (eqs 2-5)" begin
        η = pickup_efficiency(p)
        a1 = 2.0 * 0.01 * 100 / 2.0
        @test η[1] ≈ a1 / (a1 + 0.5)
        R = species_reproduction_numbers(p)
        @test R[1] ≈ 0.02 * 5.0 / 0.1 * η[1]
        # R0,C with species 2 absent is R0,1
        q = EspiraParameters(; N = (100, 0), γ = p.γ, κ = p.κ, ρ = p.ρ, π = p.π, α = p.α,
            ε = p.ε, θ = p.θ, μ = p.μ)
        @test community_reproduction_number(q) ≈ R[1]
        # limits: frequency-dependent-like (μ → 0) and density-dependent-like (μ large)
        fd = EspiraParameters(; N = p.N, γ = p.γ, κ = p.κ, ρ = p.ρ, π = p.π, α = p.α,
            ε = p.ε, θ = p.θ, μ = 1e-9)
        @test species_reproduction_numbers(fd)[1] ≈ 0.02 * 5.0 / 0.1 rtol = 1e-6
        dd = EspiraParameters(; N = p.N, γ = p.γ, κ = p.κ, ρ = p.ρ, π = p.π, α = p.α,
            ε = p.ε, θ = p.θ, μ = 1e6)
        @test species_reproduction_numbers(dd)[1] ≈ 0.02 * 5.0 / 0.1 * a1 / 1e6 rtol = 1e-5
    end

    @testset "community reproduction number (eq. 6)" begin
        Rc = community_reproduction_number(p)
        G = C.next_generation_matrix(p; reduced = false)
        @test Rc ≈ maximum(abs, eigvals(G)) rtol = 1e-10
        K = C.next_generation_matrix(p)
        @test size(K) == (2, 2)
        @test K[1, 1] * K[2, 2] - K[1, 2] * K[2, 1] ≈ K[1, 1] * K[2, 2] * (1 - 0.25 * 0.5)
        @test Rc >= minimum(species_reproduction_numbers(p))
        # R0,C > 1 iff the disease-free state of eq. 1 is unstable (Jacobian by differences)
        for scale in (0.05, 1.0, 20.0)
            q = EspiraParameters(; N = p.N, γ = p.γ, κ = p.κ, ρ = p.ρ,
                π = scale .* p.π, α = p.α, ε = p.ε, θ = p.θ, μ = p.μ)
            u0 = [q.N[1], q.N[2], 0, 0, 0, 0, 0, 0]
            J = zeros(8, 8)
            h = 1e-6
            for k in 1:8
                up = copy(u0); up[k] += h
                um = copy(u0); um[k] -= h
                f1 = zeros(8); f2 = zeros(8)
                C.espira_rhs!(f1, up, q, 0.0); C.espira_rhs!(f2, um, q, 0.0)
                J[:, k] = (f1 - f2) / 2h
            end
            s = maximum(real, eigvals(J[[3, 4, 7, 8], [3, 4, 7, 8]]))
            @test (community_reproduction_number(q) > 1) == (s > 0)
        end
        # eq. 6 as printed: equal when ε = 1 (θ₁ = θ₂), different otherwise
        q = EspiraParameters(; N = p.N, γ = p.γ, κ = p.κ, ρ = p.ρ, π = p.π, α = p.α,
            ε = (1.0, 1.0), θ = (0.4, 0.4), μ = p.μ)
        @test community_reproduction_number(q; as_printed = true) ≈
              community_reproduction_number(q)
        @test !isapprox(community_reproduction_number(p; as_printed = true),
            community_reproduction_number(p); rtol = 1e-3)
    end

    @testset "supplement (S12, S15, S17, Table S1)" begin
        @test C._espira_S12(p) ≈ C.next_generation_matrix(p) rtol = 1e-12
        R = species_reproduction_numbers(p)
        @test dilution_effect(p) == (R[2] / R[1] < C._espira_S17(p))
        recs = parameter_records(:espira2022)
        @test all(r -> r.verified, recs)
        w = C._espira_worked_example()
        @test w.overlap[1] ≈ w.overlap[2]          # θ₁ε₁ = θ₂ε₂ in the worked example
        @test w.removal_species1 ≈ 1.0             # printed: 1 PU/day
        @test w.λ1 ≈ 0.01 * 0.0001 * (0.75 * 100 + 0.25 * 140)
        # Figs S1/S2 configurations: habitat 1 inside 2 (θ₁ = 1) and 2 inside 1 (θ₂ = 1)
        @test C._espira_panel_ratio(:S2, 0.3, 0.6, 0.5, 0.2) < 1      # η₂ > r: dilution (eq. 9)
        @test C._espira_panel_ratio(:S2, 0.3, 0.4, 0.5, 0.2) > 1
        @test C._espira_panel_ratio(:fig4, 0.3, 0.6, 0.5, 1.0) ≈
              C._espira_panel_ratio(:S1, 0.3, 0.6, 0.5, 1.0)
        @test_throws ArgumentError C._espira_panel_ratio(:x, 0.3, 0.6, 0.5, 0.5)
    end

    @testset "dilution (eqs 8-11)" begin
        # eq. 8 agrees with R0,C < R0,1 when θ₁ε₁ = θ₂ε₂ (here 0.25·2 = 0.5·1)
        @test dilution_effect(p) == dilution_effect(p; criterion = :eq8)
        @test_throws ArgumentError dilution_effect(p; criterion = :x)
        rng = Xoshiro(7)
        @test all(1:500) do _
            q = C._espira_random_parameters(rng)
            dilution_effect(q) == dilution_effect(q; criterion = :eq8)
        end
        # Figs 3-4 axes
        for (η1, η2, r) in ((0.2, 0.7, 0.5), (0.6, 0.3, 0.8), (0.9, 0.9, 0.95))
            q = espira_parameters_from_efficiencies(η1, η2, r)
            @test all(pickup_efficiency(q) .≈ (η1, η2))
            R = species_reproduction_numbers(q)
            @test R[2] / R[1] ≈ r
            # complete overlap: R0,C/R0,1 - 1 = (r - η₂)/(1 - η₂) μ/(a₁ + a₂ + μ)
            a = (η1 / (1 - η1), η2 / (1 - η2))
            @test community_reproduction_number(q) / R[1] - 1 ≈
                  (r - η2) / (1 - η2) / (a[1] + a[2] + 1)
            @test dilution_effect(q) == (η2 > r)                     # eq. 9
            @test dilution_threshold(q) ≈ r
        end
        q = espira_parameters_from_efficiencies(0.3, 0.6, 0.5; θ₂ = 0.5, ε = (2.0, 1.0))
        @test q.θ == (0.25, 0.5)
        @test_throws ArgumentError espira_parameters_from_efficiencies(0.3, 0.6, 0.5;
            θ₂ = 1.0, ε = (1.0, 2.0))
        @test_throws ArgumentError espira_parameters_from_efficiencies(1.0, 0.6, 0.5)
        @test_throws ArgumentError EspiraParameters(; N = p.N, γ = p.γ, κ = p.κ, ρ = p.ρ,
            π = p.π, α = p.α, ε = p.ε, θ = (1.5, 0.5), μ = p.μ)
    end

    @testset "epidemic (eq. 1)" begin
        prob = OverlappingHabitatProblem(p; I0 = (1.0, 0.0), tspan = (0.0, 2000.0))
        sol = solve(prob)
        u = sol.u[end]
        @test u[1] + u[3] + u[5] ≈ 100 rtol = 1e-8     # no births or deaths
        @test u[2] + u[4] + u[6] ≈ 50 rtol = 1e-8
        @test u[3] + u[4] < 1e-6                        # the epidemic is over
        big = community_reproduction_number(p) > 1
        @test (u[5] + u[6] > 5) == big                  # an outbreak only if R0,C > 1
        @test to_ode_problem(prob).tspan == (0.0, 2000.0)
        @test_throws ArgumentError OverlappingHabitatProblem(p; I0 = (200.0, 0.0))
    end
end
