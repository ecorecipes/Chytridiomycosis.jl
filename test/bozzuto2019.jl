# Bozzuto & Canessa 2019, Glob Ecol Conserv 17:e00551, doi:10.1016/j.gecco.2019.e00551.
# Checked against the PDF (Table 1, eq. 1-2, Results p. 4) and the online Appendix A
# (lab/papers/bozzuto2019_si/bozzuto2019_appendixA.docx: eqs. A1-A7, Figs. A1-A6).
using Chytridiomycosis
using Test

@testset "Bozzuto 2019 seasonal Bsal" begin
    C = Chytridiomycosis
    hib = BozzutoParameters(:bozzuto2019)
    aes = BozzutoParameters(:bozzuto2019; population = :aestivation)
    hib_hi = BozzutoParameters(:bozzuto2019; activity = :high)
    aes_hi = BozzutoParameters(:bozzuto2019; population = :aestivation, activity = :high)
    hib_tab = BozzutoParameters(:bozzuto2019; as_printed = true)

    @testset "compendium" begin
        recs = parameter_records(:bozzuto2019)
        @test all(r -> r.verified, recs)
        @test time_unit(hib) === :year
        @test hib.d == 1 / 8 && aes.d == 1 / 8                     # Table 1: 1/d = 8 y
        @test hib.c2 == 0.95 && hib_tab.c2 == 0.90                  # Fig. A4 vs Table 1
        @test (hib_hi.c1, hib_hi.c2, aes_hi.c1, aes_hi.c2) == (0.85, 0.85, 0.84, 0.84)
        @test_throws ArgumentError BozzutoParameters(:bozzuto2019; population = :tropical)
        @test_throws ArgumentError BozzutoParameters(; c0 = 1, c1 = 1.2, c2 = 0, ℓc = 0, k0 = 1,
            v = 1, t1 = 0.1, t2 = 0.2, d = 0.1, δ0 = 1, δ1 = 0)
        # Derived δ0, δ1 (Appendix A1): infective periods of 40 (or 26) and 5 days, 364-day year
        b = parameter_set(:bozzuto2019)
        for (p, s) in ((hib, "hibernation"), (aes, "aestivation"))
            ps = parameter_set(:bozzuto2019; set = s)
            δmax = b.days_per_year / b.infective_period_warm
            δmin = b.days_per_year / ps.infective_period_cold
            @test (δmin + δmax) / 2 ≈ p.δ0 atol = 1e-10
            @test δmax / p.δ0 - 1 ≈ p.δ1 atol = 5e-4                # Table 1 rounds to 4 dp
        end
        # Aestivation: 1/δ = 40 + 4b days with b = (5 - 40)/(15 - 5)
        @test parameter_set(:bozzuto2019; set = "aestivation").infective_period_cold ==
              40 + 4 * (5 - 40) / (15 - 5)
        # v was chosen so that k0 is the midrange of S(t) (Appendix A2)
        for p in (hib, aes)
            S = [disease_free_susceptibles(p, t) for t in 0:1e-4:1]
            @test (minimum(S) + maximum(S)) / 2 ≈ p.k0 rtol = 2e-3
            @test disease_free_susceptibles(p, 0.0) ≈ disease_free_susceptibles(p, 1.0)
            @test sum(S[1:end-1]) / (length(S) - 1) ≈ p.v * (p.t2 - p.t1) / p.d rtol = 1e-4
        end
    end

    @testset "seasonal functions" begin
        m′ = contact_scaling(hib)
        ts = 0:1e-4:1
        c = [contact_rate(hib, t) for t in ts]
        @test (minimum(c) + maximum(c)) / 2 ≈ hib.c0 rtol = 1e-6   # c0 is the midrange
        @test sum(c[1:end-1]) / (length(ts) - 1) ≈ hib.c0 * seasonal_factor(hib) rtol = 1e-6
        @test contact_rate(hib, 0.3; scale = 1) ≈ m′ * contact_rate(hib, 0.3)
        # Fig. A1: scaled contact peaks ~1.95 in autumn, ~1.5 in spring, ~0 in January
        @test 1.9 < maximum(c) / hib.c0 < 2.0
        @test minimum(c) / hib.c0 < 0.05
        @test argmax(c) / length(ts) > 0.6
        # A1 (i)-(iii): θ(t) has its minimum θ0(1-θ1) at t = 0, maximum θ0(1+θ1) at t = 1/2
        @test bsal_mortality(hib, 0.0) ≈ hib.δ0 * (1 - hib.δ1)
        @test bsal_mortality(hib, 0.5) ≈ hib.δ0 * (1 + hib.δ1)
        @test transmission_probability(hib, 0.37) == hib.p0
        @test larval_inflow(hib, 0.6) == hib.v && larval_inflow(hib, 1.6) == hib.v
        @test larval_inflow(hib, 0.3) == 0 && larval_inflow(hib, 0.8) == 0
        @test larval_inflow(hib, 0.75; inflow = :sinusoidal) ≈ 2 * hib.d * hib.k0  # September peak at t = 3/4
        @test_throws ArgumentError larval_inflow(hib, 0.5; inflow = :pulse)
        # Sinusoidal S(t): eq. (A2b) is the d² → 0 limit, amplitude k1 = d/2π ≈ 0.02
        @test disease_free_susceptibles(hib, 0.0; inflow = :sinusoidal) ≈
              hib.k0 * (1 + hib.d / 2π) rtol = 1e-3
        # contact_scaling against a dense grid
        g = [(1 - hib.c1 * cospi(2t)) * (1 - hib.c2 * cospi(2(t - hib.ℓc))) for t in 0:1e-6:1]
        @test m′ ≈ (minimum(g) + maximum(g)) / 2 rtol = 1e-9
    end

    @testset "R0 and Floquet-R0 (Results p. 4)" begin
        # Non-seasonal R0 = c0 p0 k0/(d + δ0): 8.88 and 2.13
        @test round(basic_reproduction_number(hib), digits = 2) == 8.88
        # Aestivation: Table 1 values give 2.095. The ratio Floquet-R0/R0 matches the paper
        # (0.901 vs 1.92/2.13 = 0.901), so one aestivation input differs by ~1.6% from the
        # one used; which one is not identifiable from the paper.
        @test_broken round(basic_reproduction_number(aes), digits = 2) == 2.13
        @test basic_reproduction_number(aes) ≈ 2.0948 atol = 1e-4
        # Eq. (2), analytic: 8.27 and 1.92
        @test round(floquet_reproduction_number(hib), digits = 2) == 8.27
        @test_broken round(floquet_reproduction_number(aes), digits = 2) == 1.92  # see above
        @test floquet_reproduction_number(aes) ≈ 1.8876 atol = 1e-4
        @test floquet_reproduction_number(aes) / basic_reproduction_number(aes) ≈ 1.92 / 2.13 atol = 2e-3
        # Table 1's c2 = 0.90 gives 8.67, not 8.27
        @test round(floquet_reproduction_number(hib_tab), digits = 2) == 8.67
        # "reduces R0 by approximately 7% and 10%"
        @test round(1 - floquet_reproduction_number(hib) / basic_reproduction_number(hib), digits = 2) == 0.07
        @test round(1 - floquet_reproduction_number(aes) / basic_reproduction_number(aes), digits = 2) == 0.10
        # Numerical, eq. (A3) with S(t) of eq. (A1): 8.33 and 1.92
        @test round(floquet_reproduction_number(hib; method = :numerical), digits = 2) == 8.33
        @test_broken round(floquet_reproduction_number(aes; method = :numerical), digits = 2) == 1.92
        @test floquet_reproduction_number(aes; method = :numerical) ≈ 1.8915 atol = 1e-4
        # High activity (Fig. A4): 8.22 and 1.80
        @test round(floquet_reproduction_number(hib_hi), digits = 2) == 8.22
        @test_broken round(floquet_reproduction_number(aes_hi), digits = 2) == 1.80
        @test floquet_reproduction_number(aes_hi) ≈ 1.7762 atol = 1e-4
    end

    @testset "Floquet theory checks" begin
        for p in (hib, aes, hib_hi)
            # Exact Bacaër–Guernaoui R0 (monodromy + bisection) = paper's numerical method
            Fn = floquet_reproduction_number(p; method = :numerical)
            @test floquet_reproduction_number(p; method = :spectral) ≈ Fn rtol = 1e-7
            # ...and with S = k0 it equals eq. (2) and the time-averaged R0
            Fa = floquet_reproduction_number(p)
            @test floquet_reproduction_number(p; method = :numerical, inflow = :constant) ≈ Fa rtol = 1e-8
            @test averaged_reproduction_number(p) ≈ Fa rtol = 1e-8
            # Floquet exponent sign and Floquet-R0 agree
            @test floquet_exponent(p; method = :numerical) ≈ (p.d + p.δ0) * (Fn - 1) rtol = 1e-8
            # δ1 does not affect Floquet-R0 (p. 3)
            @test floquet_reproduction_number(C.with(p; δ1 = 0.0); method = :numerical) ≈ Fn rtol = 1e-8
            # Eq. (A4a) with k(t) of eq. (A2b) vs the exact sinusoidal S(t): agree to O(d²)
            @test floquet_reproduction_number(p; inflow = :sinusoidal) ≈
                  floquet_reproduction_number(p; method = :numerical, inflow = :sinusoidal) rtol = 1e-3
        end
        # Eq. (A4a)/(A4b) with p1 ≠ 0 against direct integration
        q = C.with(hib; p1 = 0.4)
        @test floquet_reproduction_number(q) ≈
              floquet_reproduction_number(q; method = :numerical, inflow = :constant) rtol = 1e-8
        # Seasonality proper (covariance of c with S) is < 1% in both populations
        @test 1.0 < floquet_reproduction_number(hib; method = :numerical) /
                    averaged_reproduction_number(hib; inflow = :step) < 1.01
        # Generic Floquet multipliers: constant 2×2 system, exp of eigenvalues
        A = [-1.0 0.5; 0.2 -0.3]
        μ = sort(abs.(floquet_multipliers(t -> A, 2)))
        @test μ ≈ sort(exp.(real.(C.eigvals(A)))) rtol = 1e-8
        # Periodic 2×2 NGM with constant F, V: R0 = ρ(F V⁻¹)
        F = [0.0 2.0; 1.5 0.0]
        V = [1.0 0.0; 0.0 3.0]
        @test periodic_reproduction_number(t -> F, t -> V, 2) ≈
              maximum(abs, C.eigvals(F / V)) rtol = 1e-8
    end

    @testset "general analysis (±33%, Fig. 1, Fig. A2)" begin
        for (s, f) in (("fig1a_solid", 4 / 3), ("fig1a_dotted", 4 / 3), ("fig1b", 2 / 3))
            x = parameter_set(:bozzuto2019; set = s)
            @test seasonal_factor(x.c1, x.c2, x.ℓc) ≈ f rtol = 1e-8
        end
        # One-peaked contact rate (c1 = 0) gives Floquet-R0 = R0 (p. 3)
        @test seasonal_factor(0.0, 0.7, 0.3) ≈ 1 rtol = 1e-10
        # The extremes over a grid of (c1, c2, ℓc) are 4/3 and 2/3
        fs = [seasonal_factor(c1, c2, ℓ) for c1 in 0:0.1:1, c2 in 0:0.1:1, ℓ in 0:0.02:0.98]
        @test maximum(fs) ≈ 4 / 3 rtol = 1e-6
        @test minimum(fs) ≈ 2 / 3 rtol = 1e-2
        # Target-proportion range 0.75ϱ ≤ ϱF ≤ 1.5ϱ
        @test all(isapprox.((1 / maximum(fs), 1 / minimum(fs)), (0.75, 1.5); rtol = 1e-2))
    end

    @testset "management (§2.3, Fig. 2)" begin
        # ϱF ≈ 12%, ϱ ≈ 11% (hibernating); ≈ 52%, 47% (aestivating)
        @test round(target_proportion(hib), digits = 2) == 0.12
        @test round(target_proportion(hib; basis = :nonseasonal), digits = 2) == 0.11
        @test_broken round(target_proportion(aes), digits = 2) == 0.52      # R0 offset
        @test round(target_proportion(aes), digits = 2) == 0.53
        @test round(target_proportion(aes; basis = :nonseasonal), digits = 2) == 0.48
        # Avoiding temporary outbreaks: ~5% and ~20% (1/max R0(t)); Fig. 2a-b peaks ~21 and ~5
        @test 0.04 < target_proportion(hib; basis = :peak, inflow = :step) < 0.06
        @test 0.18 < target_proportion(aes; basis = :peak, inflow = :step) < 0.22
        Rt = [instantaneous_reproduction_number(hib, t; inflow = :step) for t in 0:1e-3:1]
        @test 20 < maximum(Rt) < 22
        @test minimum(Rt) < 1.5                          # Fig. 2a trough near month 7
        # The cap at one
        @test target_proportion(C.with(hib; k0 = 1.0)) == 1.0
    end

    @testset "elasticities (Fig. A5)" begin
        e = floquet_elasticities(hib)
        @test e.R0 ≈ 1 rtol = 1e-6
        # Digitised from Fig. A5a: hib (c1, c2, ℓc) ≈ (0.72, -0.92, 2.37)
        @test e.c1 ≈ 0.72 atol = 0.05
        @test e.c2 ≈ -0.92 atol = 0.05
        @test e.ℓc ≈ 2.37 atol = 0.05
        ea = floquet_elasticities(aes)
        @test ea.c1 ≈ -0.75 atol = 0.05
        @test ea.c2 ≈ 0.50 atol = 0.05
        @test ea.ℓc ≈ 1.58 atol = 0.05
        # ℓc ranks first in both populations
        @test all(e -> abs(e.ℓc) > max(abs(e.R0), abs(e.c1), abs(e.c2)), (e, ea))
        # With Table 1's c2 = 0.90 the hibernation c2 elasticity is -0.81, not -0.92
        @test floquet_elasticities(hib_tab).c2 ≈ -0.81 atol = 0.01
    end

    @testset "ODE problem" begin
        prob = SeasonalBsalProblem(hib; tspan = (0.0, 3.0), I0 = 0.0)
        ode = to_ode_problem(prob)
        @test ode.u0 ≈ [disease_free_susceptibles(hib, 0.0), 0.0]
        sol = solve(prob)
        # Without Bsal, S(t) follows eq. (A1)
        @test all(t -> isapprox(sol(t)[1], disease_free_susceptibles(hib, t); rtol = 1e-6),
            0:0.01:3)
        # Invasion threshold: scale the host density (k0 and v) around the target proportion
        function final_I(p, ϱ; years = 30)
            q = C.with(p; k0 = ϱ * p.k0, v = ϱ * p.v)
            s = solve(SeasonalBsalProblem(q; tspan = (0.0, years), I0 = 1e-6); abstol = 1e-20)
            return s.u[end][2]
        end
        ϱF = target_proportion(hib; method = :numerical)
        @test final_I(hib, 0.9ϱF) < 1e-12               # dies out
        @test final_I(hib, 1.1ϱF) > 1e-6                # invades
        # In the linear phase I grows by exp(ξ) per year
        ϱ = 1.02ϱF
        q = C.with(hib; k0 = ϱ * hib.k0, v = ϱ * hib.v)
        s = solve(SeasonalBsalProblem(q; tspan = (0.0, 5.0), I0 = 1e-12); abstol = 1e-20,
            reltol = 1e-10)
        @test log(s(5.0)[2] / s(4.0)[2]) ≈ floquet_exponent(q; method = :numerical) rtol = 1e-3
        # Non-seasonal model: threshold at R0
        ns(ϱ) = solve(SeasonalBsalProblem(C.with(hib; k0 = ϱ * hib.k0); seasonal = false,
            tspan = (0.0, 30.0), I0 = 1e-6); abstol = 1e-20).u[end][2]
        ϱ0 = target_proportion(hib; basis = :nonseasonal)
        @test ns(0.9ϱ0) < 1e-12 && ns(1.1ϱ0) > 1e-6
        # Table 1's unscaled contact formula
        @test to_ode_problem(SeasonalBsalProblem(hib; scaled_contact = false)).p.scale == 1.0
        @test_throws ArgumentError SeasonalBsalProblem(hib; inflow = :pulse)
    end
end
