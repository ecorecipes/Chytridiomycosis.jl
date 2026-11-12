# Thermal performance curves (lab/wiki/models/thermal-performance.md). Curve parameters other
# than published fits are illustrative, chosen to test the algebra, not published values.
@testset "Thermal performance curves (Phase 0)" begin
    # Brière from cardinal temperatures: the optimum and peak are where they were put.
    b = Briere(; Tmin = 2.0, Topt = 18.0, Tmax = 27.0, rmax = 1.0)
    @test thermal_optimum(b) ≈ 18.0 atol = 1e-4
    @test b(18.0) ≈ 1.0 && b(2.0) == 0 && b(27.0) == 0 && b(30.0) == 0
    b1 = Briere(; a = 1e-4, Tmin = 2.0, Tmax = 27.0)                     # Brière-1, m = 2
    @test b1.m == 2
    @test thermal_optimum(b1) ≈ (4 * 27 + 3 * 2 + sqrt(16 * 27^2 + 9 * 2^2 - 16 * 2 * 27)) / 10 atol = 1e-4
    @test_throws ArgumentError Briere(; Tmin = 5.0, Topt = 30.0, Tmax = 25.0, rmax = 1.0)
    @test_throws ArgumentError Briere(; Tmin = 5.0, Tmax = 25.0)

    # Cardinal temperature model: exact cardinal values, and Rosso's domain condition.
    ct = CardinalTPC(; Tmin = 4.0, Topt = 21.0, Tmax = 27.0, rmax = 0.5)
    @test ct(21.0) ≈ 0.5 && ct(4.0) == 0 && ct(27.0) == 0
    @test thermal_optimum(ct) ≈ 21.0 atol = 1e-4
    # Bsal's reported optimum (10-15 degC) and growth range (5-22 degC; Martel 2013) fall outside it
    m = parameter_set(:martel2013)
    @test_throws ArgumentError CardinalTPC(; Tmin = m.T_growth_min_observed,
        Topt = (m.Topt_lo + m.Topt_hi) / 2, Tmax = m.T_no_growth, rmax = 1.0)
    bsal = Briere(; Tmin = 0.0, Topt = (m.Topt_lo + m.Topt_hi) / 2, Tmax = m.T_no_growth, rmax = 1.0)
    @test bsal(m.T_lethal) == 0 && bsal(m.T_growth_max_observed) > 0

    # Gaussian–quadratic: continuous at the optimum, zero from CTmax.
    gq = GaussianQuadratic(; rmax = 2.0, Topt = 20.0, σ = 3.0, CTmax = 26.0)
    @test gq(20.0) == 2.0 && gq(20.0 + 1e-9) ≈ 2.0 && gq(26.0) == 0 && gq(30.0) == 0
    @test gq(14.0) ≈ 2 * exp(-1)

    # Sharpe–Schoolfield: Arrhenius below the inactivation temperatures.
    k = 8.617333262e-5
    ss = SharpeSchoolfield(; r_ref = 1.0, E = 0.65, Eh = 3.0, Th = 30.0)
    @test ss(20.0) ≈ 1 / (1 + exp(3.0 / k * (1 / 303.15 - 1 / 293.15)))
    ssl = SharpeSchoolfield(; r_ref = 1.0, E = 0.65, Eh = 3.0, Th = 30.0, El = 2.0, Tl = 5.0)
    @test ssl(5.0) < ss(5.0) && 15 < thermal_optimum(ssl) < 30
    @test_throws ArgumentError SharpeSchoolfield(; r_ref = 1.0, E = 0.65, Eh = 3.0, Th = 30.0, El = 2.0)

    # Ratkowsky.
    rk = Ratkowsky(; b = 0.05, Tmin = 3.0, c = 0.3, Tmax = 28.0)
    @test rk(3.0) == 0 && rk(28.0) == 0 && rk(15.0) ≈ (0.05 * 12 * (1 - exp(0.3 * -13)))^2

    # AD-friendly: generic in the number type (dual-like check with BigFloat).
    @test b(big"18.0") isa BigFloat && ct(big"21.0") ≈ 0.5

    # Kelly et al. 2024: the 13 Bsal isolate polynomials (S1 Text, Table E).
    isolates = unique(r.set for r in parameter_records(:kelly2024) if !isempty(r.set))
    @test length(isolates) == 13
    eifel = IsolateTPC(:kelly2024, "Eifel2019")
    @test eifel.pathogen == "Bsal" && eifel.source === :kelly2024
    @test eifel(10.0) ≈ evalpoly(10.0, [parameter_set(:kelly2024; set = "Eifel2019")[Symbol("c$i")] for i in 0:4])
    @test eifel(-1.0) == 0 && eifel(31.0) == 0                          # no extrapolation
    @test thermal_optimum(eifel.curve; range = (0.0, 30.0)) ≈ 6.05 atol = 0.05
    @test_throws ArgumentError IsolateTPC(:kelly2024, "Nowhere2030")
    @test_throws ArgumentError IsolateTPC(:woodhams2008, "23C")
    opt = [thermal_optimum(IsolateTPC(:kelly2024, i).curve; range = (0.0, 30.0)) for i in isolates]
    @test minimum(opt) < 7 && maximum(opt) > 24                         # optima span 6 to 24+ degC
    @test PolynomialTPC([-1.0, 0.1])(5.0) == 0 && PolynomialTPC([-1.0, 0.1]; nonnegative = false)(5.0) ≈ -0.5

    # Woodhams et al. 2008: growth rate from the characteristic equation, and rates at
    # four temperatures as an interpolated curve.
    w23 = WoodhamsParameters(:woodhams2008; set = "23C")
    λ = intrinsic_growth_rate(w23)
    @test w23.F * w23.mZ * exp(-λ * w23.s) ≈ λ + w23.mZ + w23.dZ
    @test λ > log(w23.F) / w23.s - 0.01 && round(λ; digits = 3) == 0.044
    Ts, mZ = thermal_data(:woodhams2008, :mZ; set_prefix = "table1")
    @test Ts == [4.0, 14.5, 23.0, 28.0]
    enc = InterpolatedTPC(Ts, mZ)
    @test enc(23.0) == 0.18 && enc(18.75) ≈ (0.052 + 0.18) / 2 && enc(30.0) == 0
    @test InterpolatedTPC(Ts, mZ; extrapolation = :flat)(30.0) == 0.21
    @test_throws ArgumentError thermal_data(:woodhams2008, :mZ)          # two tables at 23 degC
    @test_throws ArgumentError InterpolatedTPC([1.0], [2.0])

    # Seasonal forcing composed with a curve.
    T = SeasonalTemperature(; mean = 15.0, amplitude = 8.0, peak = 200.0)
    @test T(200.0) ≈ 23.0 && T(200.0 + 365 / 2) ≈ 7.0 && T(565.0) ≈ T(200.0)
    k_t = temperature_dependent(eifel, T)
    @test k_t(200.0) ≈ eifel(23.0)
    host = HostTPC("Acris crepitans", "illustrative performance", gq, :none)
    @test host(20.0) == 2.0
end
