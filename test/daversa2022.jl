# Daversa, Bosch, Manica, Garner & Fenton 2022, Am Nat 200(4), doi:10.1086/720638.
# Equations checked against the PDF (p. 5–6) and the authors' R code (Zenodo 5068208, saved as
# lab/papers/daversa2022_si1.Rmd). The data are Dryad 10.5061/dryad.2v6wwpznk, copied to
# data/daversa2022/. The parameters are posterior medians digitised from Fig. 4a, since the
# paper prints no numbers and the supplemental PDF could not be obtained.
using Chytridiomycosis
using Test

@testset "Daversa 2022 multi-host shared pool" begin
    C = Chytridiomycosis
    p = DaversaParameters(:daversa2022)
    ps = parameter_set(:daversa2022)
    tanks = load_daversa2022(; three_species_only = true)
    all36 = load_daversa2022()

    @testset "compendium and data" begin
        @test p.species == [:salamander, :toad, :frog]
        @test p.γ == 1.0                                    # p. 5: spores last 1 day
        recs = parameter_records(:daversa2022)
        @test all(r -> r.verified, filter(r -> r.provenance != "F", recs))
        @test !any(r -> r.verified, filter(r -> r.provenance == "F", recs))  # digitised
        # Text, p. 6: salamander β ~80% below the anurans, frog shedding ~4× salamander's,
        # recovery roughly equal for all three
        @test 0.1 < p.β[1] / ((p.β[2] + p.β[3]) / 2) < 0.3
        @test 3 < p.λ[3] / p.λ[1] < 5
        @test maximum(p.σ) / minimum(p.σ) < 1.1
        # p. 3: doses of 18,000–24,000 zoospores per litre in 80 L. The code's four inputs are
        # 1.45e6 to 1.95e6, which is just outside the upper bound (1.92e6).
        d = daversa_doses()
        @test d.Z0 == ps.Z0 && length(d.doses) == 3
        @test [t for (t, _) in d.doses] == [10.0, 17.0, 25.0]
        @test all(z -> ps.dose_conc_lo * ps.tank_volume <= z <= 1.02 * ps.dose_conc_hi * ps.tank_volume,
            vcat(d.Z0, last.(d.doses)))
        # Data: 36 tanks, of which 21 have three species. The 15 fixed-design tanks start with
        # 24 animals. The random tanks in design.csv hold 22–32, not 24 as the text says.
        @test length(all36) == 36 && length(tanks) == 21
        fixed = filter(t -> t.treatment != "random", all36)
        @test all(t -> sum(t.H0) == ps.hosts_per_tank, fixed)
        @test extrema(sum(t.H0) for t in all36 if t.treatment == "random") == (22.0, 32.0)
        @test all(t -> t.day_final in (28.0, 29.0), all36)
        @test all(t -> all(>=(0), t.μ), all36)
        # p. 4: 166 of 525 salamanders (31%) infected in the initial phase. That uses the
        # 0.01 GE threshold. The R code uses load > 0, and the prevalences here are per tank.
        @test all(t -> isnan(t.prev_initial[2]) && isnan(t.prev_initial[3]), tanks)
    end

    @testset "equations and solver" begin
        t = tanks[1]
        prob = MultiHostPoolProblem(p, t)
        sol = solve(prob)
        @test C.SciMLBase.successful_retcode(sol)
        # H decays exponentially at the tank's loss rate (dH/dt = −μH)
        @test sol.u[end][1:3] ≈ t.H0 .* exp.(-t.μ .* t.day_final) rtol = 1e-6
        @test sol.u[end][1:3] ≈ t.Hfinal rtol = 1e-6
        # The doses arrive as jumps in Z
        z(tt) = sol(tt)[7]
        @test z(10.0 + 1e-9) - z(10.0 - 1e-9) ≈ 1.94e6 rtol = 1e-3
        # Without hosts the pool decays at γ: Z(t) = Z0 e^{−γt} between doses
        p0 = DaversaParameters(; species = [:x], β = [0.0], λ = [0.0], σ = [0.0], γ = 1.0)
        s0 = solve(MultiHostPoolProblem(p0, [0.0]; Z0 = 100.0, tspan = (0.0, 2.0)))
        @test s0.u[end][3] ≈ 100 * exp(-2) rtol = 1e-6
        # Any number of species: one species, five species
        p5 = DaversaParameters(; species = Symbol.("s", 1:5), β = fill(1e-7, 5),
            λ = fill(1e5, 5), σ = fill(0.5, 5), γ = 1.0)
        s5 = solve(MultiHostPoolProblem(p5, fill(5.0, 5); I0 = [1.0, 0, 0, 0, 0],
            tspan = (0.0, 20.0)))
        @test length(s5.u[end]) == 11 && all(>=(-1e-9), s5.u[end])
        @test_throws ArgumentError MultiHostPoolProblem(p, [1.0, 2.0])
        @test_throws ArgumentError MultiHostPoolProblem(p, [1.0, 1, 1]; I0 = [2.0, 0, 0])
        @test to_ode_problem(prob) isa C.SciMLBase.ODEProblem
    end

    @testset "R0: analytic formula vs next-generation matrix" begin
        for t in tanks, μ in (t.μ, zeros(3))
            R = species_reproduction_numbers(p, t.H0; μ)
            ngm = next_generation_matrix(p, t.H0; μ)
            @test ngm.R0 ≈ sum(R) rtol = 1e-8
            @test basic_reproduction_number(p, t.H0; μ) ≈ ngm.R0 rtol = 1e-8
            # Counting host → spore → host as two generations gives √R0
            @test next_generation_matrix(p, t.H0; μ, spores = :generation).R0^2 ≈ sum(R) rtol = 1e-8
            # Host-to-host block: rank one, so the spectral radius equals the trace (sum)
            K = ngm.K[1:3, 1:3]
            @test sum(K[i, i] for i in 1:3) ≈ ngm.R0 rtol = 1e-8
        end
        # Five species: still the sum
        p5 = DaversaParameters(; species = Symbol.("s", 1:5), β = [1, 2, 3, 4, 5] .* 1e-7,
            λ = [5, 4, 3, 2, 1] .* 1e5, σ = [0.3, 0.5, 0.7, 0.2, 0.9], γ = 1.0)
        H5 = [10.0, 4, 0, 7, 3]
        @test next_generation_matrix(p5, H5; μ = fill(0.01, 5)).R0 ≈
              basic_reproduction_number(p5, H5; μ = fill(0.01, 5)) rtol = 1e-8
        # The threshold: with hosts held constant (μ = 0) and no further doses, a tank with
        # R0_TOT > 1 sustains infection and one with R0_TOT < 1 loses it
        R00(t) = basic_reproduction_number(p, t.H0)
        t_hi, t_lo = argmax(R00, tanks), argmin(R00, tanks)
        @test R00(t_lo) < 1 < R00(t_hi)
        frozen(t) = MultiHostPoolProblem(p, t.H0; I0 = [0.0, 1.0, 1.0], tspan = (0.0, 200.0))
        Itot(t) = sum(solve(frozen(t)).u[end][4:6])
        @test Itot(t_hi) > 2 && Itot(t_lo) < 0.1
    end

    @testset "published results (Fig. 4, Fig. S7)" begin
        R = tank_reproduction_numbers(p)
        Rtot = sum.(R)
        # Fig. 4b medians, digitised: 0.143, 0.247, 0.911. These are point estimates at the
        # median parameters, so a 15% tolerance is allowed.
        @test isapprox(C._median(getindex.(R, 1)), 0.143; rtol = 0.15)
        @test isapprox(C._median(getindex.(R, 2)), 0.247; rtol = 0.15)
        @test isapprox(C._median(getindex.(R, 3)), 0.911; rtol = 0.15)
        # Salamanders contribute least to R0: below frogs in every tank and at most 23% of
        # R0_TOT. They exceed toads in 4 salamander-heavy tanks.
        @test all(r -> r[1] < r[3], R)
        @test count(r -> r[1] < r[2], R) == 17
        @test maximum(r[1] / sum(r) for r in R) < 0.25
        # Fig. 4c: R0_TOT falls as the salamander fraction rises
        ss = [t.H0[1] / (24 - t.H0[1]) for t in tanks]
        @test C._pearson(ss, Rtot) < -0.5
        # Fig. 4c–d: per-tank medians, digitised. Our point estimates track them (r > 0.999)
        # and are 5–10% lower.
        pub = C._daversa_fig4_medians()
        pubv = [pub[t.tank] for t in tanks]
        @test C._pearson(Rtot, pubv) > 0.999
        @test all(0.90 .< Rtot ./ pubv .< 0.95)
        # Fig. 4d and text: "varying the initial ratio of frog and toad abundances had little
        # consequence for R0_TOT". Both our values and the digitised published medians
        # correlate with the frog:toad ratio at r ≈ 0.45 (discrepant target
        # daversa2022_frog_toad_ratio): an interpretation, not a numerical mismatch.
        ft = [t.H0[3] / t.H0[2] for t in tanks]
        @test_broken abs(C._pearson(ft, Rtot)) < 0.3
        @test C._pearson(ft, Rtot) ≈ 0.449 atol = 0.002
        @test C._pearson(ft, pubv) ≈ 0.440 atol = 0.01
        # Fig. S7 / p. 6: observed vs median predicted numbers infected, r = 0.70, R² = 0.50
        r = C._daversa_fit_correlation(p)
        @test isapprox(r, 0.70; atol = 0.02)
        @test isapprox(r^2, 0.50; atol = 0.02)
        # The printed R0 formula (denominator μ_i, no σ_i) does not give Fig. 4b. It gives
        # Inf where no animals were lost (μ = 0). The code and the NGM both use μ_i + σ_i.
        Rp = tank_reproduction_numbers(p; μ = :tank, as_printed = true)
        @test_broken isapprox(C._median(getindex.(Rp, 1)), 0.143; rtol = 0.15)
        @test C._median(getindex.(Rp, 1)) ≈ 4.26 atol = 0.01
        @test any(r -> isinf(r[2]), Rp)
    end

    @testset "likelihood" begin
        ll = daversa_loglikelihood(p, tanks)
        @test isfinite(ll)
        @test ll ≈ -219.02 atol = 0.01                     # regression at the digitised medians
        # Setting frog shedding to zero makes the fit worse
        pbad = DaversaParameters(; species = p.species, β = p.β, λ = [p.λ[1], p.λ[2], 0.0],
            σ = p.σ, γ = p.γ)
        @test daversa_loglikelihood(pbad, tanks) < ll
    end
end
