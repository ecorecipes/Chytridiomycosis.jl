# Canessa et al. 2018, J Appl Ecol 55:1987, doi:10.1111/1365-2664.13089.
# Checked against the review copy of the paper (Table 1, Figs 2-4, SI eqs S1-S12, Table S1)
# and the archived draws (Dryad doi:10.5061/dryad.0jh5p, data/canessa2018/).
# The full 1,000-draw ensembles come from the session cache (_canessa_reference), which the
# validation targets fill when the whole suite runs.
using Chytridiomycosis
using LinearAlgebra: I as Id, dot
using Chytridiomycosis: median, quantile
using Test

@testset "Canessa 2018 Bsal IPM" begin
    C = Chytridiomycosis
    c = parameter_set(:canessa2018)
    draws = load_canessa2018()
    logistic(x) = 1 / (1 + exp(-x))

    @testset "compendium and data" begin
        recs = parameter_records(:canessa2018)
        unverified = Set(r.symbol for r in recs if !r.verified)
        @test unverified == Set([:growth_series_steps, :Lmax_best, :growth_best, :I0_count])
        @test c.step_length == 3 && c.nsteps == 30
        @test (c.z_lower, c.z_upper, c.mesh_points) == (-10, 15, 300)
        @test c.μG0 ≈ log(200) atol = 1e-4                       # Table S1: 5.2983
        @test time_unit(CanessaParameters(:canessa2018)) === :day
        @test time_unit(CanessaSIParameters(:canessa2018)) === :day
        @test length(draws) == 1000 && [d.id for d in draws] == 1:1000
        @test all(d -> d.transmission[1] == 0, draws)            # T_0 is 0 in every row
        # S0 = (count on the 475 m × 10 m transect) / 0.475 ha
        area = c.transect_length * c.search_width / 1e4
        @test all(d -> isapprox(d.S0 * area, round(d.S0 * area); atol = 1e-8), draws)
        @test median([d.β0 for d in draws]) ≈ c.β0 atol = 5e-5
        @test quantile([d.β0 for d in draws], 0.025) ≈ c.β0_lower rtol = 0.01
        @test quantile([d.β0 for d in draws], 0.975) ≈ c.β0_upper rtol = 0.01
        # Shared with Bozzuto & Canessa (2019): their c0 = β0 (364/3 steps per year)/1.1643
        b = parameter_set(:bozzuto2019)
        hib = BozzutoParameters(:bozzuto2019)
        @test hib.c0 * b.contact_ratio ≈ c.β0 * b.days_per_year / c.step_length rtol = 1e-3
        @test hib.k0 ≈ c.S0 rtol = 0.01                          # 500 against 503 ha⁻¹
        @test b.lifespan == 8                                      # the 8 years of SI line 54
    end

    @testset "fits" begin
        # Table S1 b4-b6 from the inferred best estimates (10 steps, Lmax 10⁴, growth 25)
        b4, b5, b6 = canessa_growth_coefficients(c.Lmax_best, c.growth_best)
        @test (round(b4; digits = 4), round(b5; digits = 4), round(b6; digits = 4)) ==
              (c.b4, c.b5, c.b6)
        # A neighbouring guess does not reproduce them
        @test round(canessa_growth_coefficients(c.Lmax_best, 30.0)[1]; digits = 4) != c.b4
        @test_throws ArgumentError canessa_growth_coefficients(1e4, 1.0)
        # The logistic regression recovers exact logistic data, and takes zeros
        z = [-1.0, log(10), log(100), log(1000), log(1e4)]
        y = logistic.(2.0 .- 0.7 .* z)
        @test all(C._canessa_logistic_fit(z, y) .≈ (2.0, -0.7))
        b2, b3 = C._canessa_logistic_fit(z, collect(draws[1].transmission))
        @test b3 > 0 && isfinite(b2)
        @test_throws ArgumentError C._canessa_logistic_fit(z, [0.5, 1.2, 0.1, 0.1, 0.1])
        # Table S1 best estimates are far from the archived draws
        p = CanessaParameters(:canessa2018)
        @test canessa_survival(p, log(100)) < 0.65
        @test median([d.survival[3] for d in draws]) > 0.89
    end

    @testset "parameters and actions" begin
        acts = canessa_actions()
        @test [a.name for a in acts] == collect(C.CANESSA_ACTION_NAMES)
        a(n) = canessa_action(n)
        @test a(:b).β_low_factor == 0.5 && a(:b).low_load == 100 && a(:b).growth_factor == 1
        @test a(:c).β_low_factor == 0.5 && a(:c).growth_factor == 0.2
        @test a(:d).β_factor ≈ 0.02 && a(:e).β_factor ≈ 1 - 0.8 * 0.98
        @test a(:f).density_factor == 0.5 && a(:g).density_factor ≈ 0.1
        @test a(:h).removal == 0.9
        @test_throws ArgumentError canessa_action(:z)
        @test_throws ArgumentError CanessaAction(:x; removal = 1.0)

        p = CanessaParameters(:canessa2018)
        @test (p.b0, p.b1, p.s0, p.b2, p.b3, p.β0, p.S0) ==
              (c.b0, c.b1, c.s0, c.b2, c.b3, c.β0, c.S0)
        @test (p.b4, p.b5, p.b6, p.Lmax) == (c.b4, c.b5, c.b6, c.Lmax_best)
        pc = CanessaParameters(:canessa2018; action = a(:c))
        @test pc.b4 ≈ canessa_growth_coefficients(c.Lmax_best, 0.2 * c.growth_best)[1]
        @test canessa_growth_mean(pc, 5.0) < canessa_growth_mean(p, 5.0)   # slower growth

        d = draws[7]
        q = CanessaParameters(d)
        @test q.s0 ≈ logistic(q.b0) && q.S0 == d.S0 && q.β0 == d.β0 && q.Lmax == d.Lmax
        @test CanessaParameters(d; s0 = :elicited).s0 == d.survival[1]
        @test_throws ArgumentError CanessaParameters(d; s0 = :other)
        qh = CanessaParameters(d; action = a(:h))
        @test qh.s0 ≈ 0.1 * q.s0 && canessa_survival(qh, 3.0) ≈ 0.1 * canessa_survival(q, 3.0)
        @test CanessaParameters(d; action = a(:g)).S0 ≈ 0.1 * d.S0
        qb = CanessaParameters(d; action = a(:b))
        @test canessa_transmission(qb, log(100)) ≈ 0.5 * canessa_transmission(q, log(100))
        @test canessa_transmission(qb, log(101)) == canessa_transmission(q, log(101))
        @test canessa_growth_mean(q, 12.0) == log(q.Lmax)               # the cap
        @test canessa_growth_variance(q, 12.0) == canessa_growth_variance(q, log(q.Lmax))
        @test canessa_growth_variance(q, 0.0) == c.σ2G
        zs = range(-10, 15; length = 20001)
        @test sum(canessa_initial_load.(Ref(q), zs)) * step(zs) ≈ 1 atol = 1e-8
    end

    @testset "kernel and R0" begin
        q = CanessaParameters(draws[3])
        k = C._canessa_kernel(q, C._CANESSA_DEFAULT_MESH, :truncate)
        kp = C._canessa_kernel(q, C._CANESSA_DEFAULT_MESH, :printed)
        top = findlast(<=(log(q.Lmax)), k.z)
        @test all(iszero, k.P[top+1:end, :])
        # Columns hold the survivors s(z) up to eviction, which the paper does not correct:
        # below z ≈ 0 the quadratic of eq. S9 sends loads under -10, and just below the cap
        # it overshoots ln Lmax so that up to 0.25% of the mass passes 15 (draw 3).
        inside = findall(z -> 0 <= z <= 12, k.z)
        @test all(j -> isapprox(sum(k.P[:, j]), canessa_survival(q, k.z[j]); rtol = 5e-3), inside)
        @test all(j -> isapprox(sum(kp.P[:, j]), canessa_survival(q, kp.z[j]); rtol = 5e-3), inside)
        @test sum(k.P[:, 1]) < 0.5canessa_survival(q, k.z[1])
        @test sum(k.g0) ≈ 1 atol = 1e-10
        # The plain-Julia solve and matrix product agree with LinearAlgebra
        @test C._canessa_fundamental(k.P, k.g0) ≈ (Id - k.P) \ k.g0 rtol = 1e-10
        @test C._canessa_mul(k.P, k.g0) ≈ k.P * k.g0 rtol = 1e-12
        # R0 as the sum over a lifetime of infection opportunities
        R0 = basic_reproduction_number(q)
        n, acc = copy(k.g0), 0.0
        for _ in 1:2000
            acc += dot(k.β, n)
            n = k.P * n
        end
        @test R0 ≈ q.S0 * q.s0 * acc rtol = 1e-10
        @test basic_reproduction_number(q; cap = :printed) < R0       # loads above Lmax die faster
        @test basic_reproduction_number(q; mesh = (-10, 15, 150)) ≈ R0 rtol = 0.01
        @test_throws ArgumentError basic_reproduction_number(q; cap = :other)
        # R0 is linear in S0 and in β0 (actions f, g, d's transmission factor)
        @test basic_reproduction_number(CanessaParameters(draws[3];
            action = canessa_action(:f))) ≈ 0.5R0 rtol = 1e-12
        @test basic_reproduction_number(CanessaParameters(draws[3];
            action = CanessaAction(:x; β_factor = 0.3))) ≈ 0.3R0 rtol = 1e-12
        # The Table S1 point model is far below the ensemble (R0 ≈ 5.1)
        @test 4.9 < basic_reproduction_number(CanessaParameters(:canessa2018)) < 5.3
    end

    @testset "IPM dynamics" begin
        q = CanessaParameters(draws[11])
        area = c.transect_length * c.search_width / 1e4                 # 0.475 ha
        # Per hectare (units = :density): S0 ha⁻¹ and one infected adult per km²
        free = solve(CanessaIPMProblem(q; units = :density, I0 = 0))
        @test free.S ≈ q.S0 .* q.s0 .^ (0:30) && all(iszero, free.I)
        @test free.t == 0:3:90
        sol = solve(CanessaIPMProblem(q; units = :density))
        @test sol.I[1] ≈ 0.01 atol = 1e-8
        @test size(sol.n) == (300, 31) && sum(sol.n[:, end]) ≈ sol.I[end]
        @test issorted(sol.S; rev = true) && sol.S[end] / free.S[end] < 0.05
        # New infections are S s0 φ: susceptible losses beyond natural deaths
        t = 10
        @test sol.S[t] * q.s0 - sol.S[t+1] ≈ sol.S[t] * q.s0 * sol.φ[t]
        @test canessa_density_ratio(q, q; units = :density) ≈ sol.S[end] / (q.S0 * q.s0^30)
        # On the 0.475-ha transect (units = :transect, the default): S0 × 0.475 adults and
        # one infected adult; the same as per hectare with β0 × 0.475 and I0 = 1/0.475 ha⁻¹
        st = solve(CanessaIPMProblem(q))
        @test st.S[1] ≈ area * q.S0
        @test st.I[1] ≈ 1 atol = 1e-8
        qs = CanessaParameters(; (f => getfield(q, f) for f in fieldnames(CanessaParameters))...,
            β0 = area * q.β0)
        sd = solve(CanessaIPMProblem(qs; units = :density, I0 = 1 / area))
        @test st.S ≈ area .* sd.S rtol = 1e-10
        @test canessa_density_ratio(q, q) ≈ st.S[end] / (area * q.S0 * q.s0^30)
        @test canessa_density_ratio(q, q) > canessa_density_ratio(q, q; units = :density)
        qd = CanessaParameters(draws[11]; action = canessa_action(:d))
        @test canessa_density_ratio(qd, q) > 0.99
        @test_throws ArgumentError CanessaIPMProblem(q; I0 = -1)
        @test_throws ArgumentError CanessaIPMProblem(q; units = :other)
        # Movement (Fig. 3)
        surv = canessa_cohort_survival(q)
        @test surv[1] == 1 && surv[end] < 0.01 && surv[end-1] >= 0.01 && issorted(surv; rev = true)
        mv = canessa_dispersal(q)
        @test mv.distance[2] ≈ 1000 * 0.5 * 3 / 365
        @test mv.max_distance ≈ mv.distance[end] && mv.lifespan == length(surv) - 1
        @test mv.mean_distance < mv.max_distance
        @test canessa_dispersal(CanessaParameters(draws[11]; action = canessa_action(:h))).lifespan <= 2
    end

    @testset "ensemble" begin
        sub = draws[1:12]
        ens = canessa_ensemble(canessa_action(:e); draws = sub)
        p5 = CanessaParameters(sub[5]; action = canessa_action(:e))
        @test ens[5].R0 ≈ basic_reproduction_number(p5) rtol = 1e-12
        @test ens[5].ratio ≈ canessa_density_ratio(p5, CanessaParameters(sub[5])) rtol = 1e-12
        ensd = canessa_ensemble(canessa_action(:e); draws = sub[1:2], units = :density)
        @test ensd[2].ratio ≈ canessa_density_ratio(CanessaParameters(sub[2];
            action = canessa_action(:e)), CanessaParameters(sub[2]); units = :density) rtol = 1e-12
        @test_throws ArgumentError canessa_ensemble(; draws = sub[1:1], units = :other)
        @test ens[5].max_distance ≈ canessa_dispersal(p5).max_distance
        s = canessa_summary(ens)
        @test s.R0.lower <= s.R0.median <= s.R0.upper
        @test 0 <= s.P_R0_below_1 <= 1
    end

    @testset "Fig. 2 and Fig. 3 (all 1,000 draws)" begin
        S(n; kw...) = C._canessa_reference(n; kw...).summary
        a = S(:a)
        # Results: median R0 9.6 (5.8-14.3). Ours 9.43 (5.60-14.35) with loads capped at Lmax.
        @test_broken round(a.R0.median; digits = 1) == 9.6
        @test a.R0.median ≈ 9.435 atol = 0.01
        @test_broken round(a.R0.lower; digits = 1) == 5.8         # ours 5.60
        @test round(a.R0.upper; digits = 1) == 14.3
        # The cap exactly as printed: 8.25 (4.97-12.23), every quantile ~14% low
        ap = S(:a; cap = :printed)
        @test_broken round(ap.R0.median; digits = 1) == 9.6
        @test ap.R0.median / a.R0.median ≈ 0.874 atol = 0.01
        # b: 'did not reduce R0' (ours -5%); c: text says R0 increases (ours +2.6%) but
        # Fig. 2 plots c at 8.3, below a.
        @test 0.9 < S(:b).R0.median / a.R0.median < 1
        @test S(:c).R0.median > a.R0.median
        @test_broken abs(S(:c).R0.median - 8.3) < 0.3              # Fig. 2 contradicts the text
        # d: R0 < 1 in every draw, population kept; e: < 2.5% of draws with R0 < 1
        @test S(:d).R0.upper < 1 && S(:d).ratio.lower > 0.99
        @test S(:e).P_R0_below_1 < 0.025 && abs(S(:e).R0.median - 2.4) < 0.25
        # Final densities, simulating on the 0.475-ha transect (the units of β0 implied by eq.
        # S6) with one infected adult: e 0.951 (Fig. 2: 1.0, collapse delayed), f 0.081 (0.083)
        @test abs(S(:e).ratio.median - 1.0) < 0.05
        @test S(:e).ratio.median ≈ 0.951 atol = 0.005
        @test abs(S(:f).R0.median - 4.8) < 0.25
        @test abs(S(:f).ratio.median - 0.083) < 0.03
        # Per hectare (units = :density) the outbreaks run 2.1 times faster: e 0.69, f 0.006
        @test_broken abs(S(:f; units = :density).ratio.median - 0.083) < 0.03
        @test S(:f; units = :density).ratio.median < 0.01
        # g: 'approximately 50%' chance of eradication (ours 59%); h: R0 < 1, ratio 0
        @test 0.4 < S(:g).P_R0_below_1 < 0.6
        @test S(:h).R0.upper < 1 && S(:h).ratio.upper < 1e-3
        # a: effective extirpation within three months
        @test a.ratio.median < 0.05
        # Fig. 3: less than 100 m; slowed growth (c-e) increases distances; h minimal
        @test a.mean_distance.median < 100 && 30 < a.max_distance.median < 55
        @test S(:c).max_distance.median > a.max_distance.median
        @test S(:d).max_distance == S(:e).max_distance == S(:c).max_distance
        @test S(:h).max_distance.upper < 10
    end

    @testset "simplified SI model (eqs S11-S12, Fig. 4)" begin
        p = CanessaSIParameters(:canessa2018)
        @test (p.β, p.S0, p.s0, p.s̄) == (c.β0, c.S0, c.s0, c.s_bar)
        R0 = basic_reproduction_number(p)
        @test R0 ≈ c.S0 * c.s0 * c.β0 / (1 - c.s_bar)
        @test_broken round(R0; digits = 1) == 7.9                  # ours 7.62 (7.69 without s0)
        @test round(required_reduction(p, :β); digits = 2) == 0.87   # Fig. 4a intercept
        @test required_reduction(p, :S0) == required_reduction(p, :β)
        @test round(required_reduction(p, :survival); digits = 2) == 0.75
        @test round(required_reduction(p, (:β, :S0, :survival)); digits = 2) == 0.42
        # Printed error: 'transmission and host survival by 87% each (Fig 4a)'; it is 54%
        @test_broken round(required_reduction(p, (:β, :survival)); digits = 2) == 0.87
        @test required_reduction(p, (:β, :survival)) ≈ 0.5366 atol = 1e-4
        m = required_reduction(p, :β; mS0 = 0.5)
        @test basic_reproduction_number(p; mβ = m, mS0 = 0.5) ≈ 1 atol = 1e-10
        @test required_reduction(p, :β; mS0 = 0.99) == 0
        @test_throws ArgumentError required_reduction(p, :other)
        @test_throws ArgumentError basic_reproduction_number(p; mβ = 1.5)
        # Eq. S6 with S0 per ha gives 0.0039, not 0.007; with S0 as the count on the 0.475-ha
        # transect it gives 0.00697, and the archived draws then imply a valid s̄ in every row
        @test_broken round(final_size_transmission(c.S0, c.S_inf, c.s_bar); digits = 3) == 0.007
        @test final_size_transmission(c.S0, c.S_inf, c.s_bar) ≈ 0.003947 atol = 1e-6
        area = c.transect_length * c.search_width / 1e4
        @test round(final_size_transmission(area * c.S0, c.S_inf, c.s_bar); digits = 3) == 0.007
        sbar(d, A) = 1 - d.β0 * (A * d.S0 - c.S_inf) / log(A * d.S0 / c.S_inf)
        @test all(d -> 0 < sbar(d, area) < 1, draws)
        @test count(d -> sbar(d, 1.0) < 0, draws) > 100
        @test_throws ArgumentError final_size_transmission(7, 503, 0.5)
        # Dynamics: no infection gives S0 s0^t; more reduction keeps more hosts (Fig. 4d)
        free = solve(CanessaSIProblem(p; I0 = 0))
        @test free.S ≈ c.S0 .* c.s0 .^ (0:30)
        r = [final_density_ratio(CanessaSIProblem(p; mβ = m)) for m in 0:0.1:1]
        @test issorted(r) && r[end] ≈ 1
        @test final_density_ratio(CanessaSIProblem(p); quasi_extinction = 0.01) == 0
        sol = solve(CanessaSIProblem(p; mS0 = 0.5))
        @test sol.S[1] ≈ 0.5c.S0 && sol.t == 0:3:90
        # Fig. 4d: reducing S0 alone scales the final density roughly linearly once R0 < 1
        @test final_density_ratio(CanessaSIProblem(p; mS0 = 0.9, mβ = 0)) < 0.1 + 1e-9
    end

    @testset "printed values" begin
        p = CanessaParameters(:canessa2018)
        # Table S1 prints s0 = s(z = 0) = 0.9906; the printed b0 gives 0.9875
        @test_broken round(logistic(c.b0); digits = 4) == c.s0
        # SI line 53: s0 'very similar' to an 8-year life span (0.99897 per 3 days)
        @test_broken round(hazard_to_survival(1 / (365 * 8), 3); digits = 3) == round(c.s0; digits = 3)
        @test p.s0 == 0.9906
    end
end
