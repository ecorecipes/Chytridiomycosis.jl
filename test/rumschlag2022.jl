# Rumschlag et al. 2022 J Anim Ecol 91:170 (doi:10.1111/1365-2656.13612): zoospore
# persistence model. Oracles were read from the PUBLISHED version (Results p. 176, Table 1,
# Figs 2-4); the bioRxiv preprint (doi:10.1101/2021.03.17.435818) prints the same numbers.
# The refit uses the authors' figshare deposit (doi:10.6084/m9.figshare.16814701).
# Wiki: lab/wiki/models/bd-rumschlag2022-zoospore-persistence.md
using Chytridiomycosis
using Test
using Distributions: Chisq, ccdf

@testset "Rumschlag 2022 zoospore persistence" begin
    C = Chytridiomycosis
    waters = (:ASW, :TroutCreek, :Flatwoods, :GreenSwamp)
    published_lifespan = (ASW = 14.4, TroutCreek = 13.3, Flatwoods = 11.4, GreenSwamp = 9.4)

    @testset "live/dead zoospore pool" begin
        pool = ZoosporePool(1.6, 1.1)
        ZL, ZD = pool_state(pool, 100.0, 0.7)
        # closed form against the ODE of the pool
        f(u, _, t) = collect(pool_derivatives(pool, u[1], u[2]))
        sol = solve(C.ODEProblem(f, [100.0, 0.0], (0.0, 0.7)), C.Tsit5(); reltol = 1e-12,
            abstol = 1e-12)
        @test sol.u[end] ≈ [ZL, ZD] rtol = 1e-8
        # with initial dead zoospores, and in the limit d = m
        sol2 = solve(C.ODEProblem(f, [100.0, 30.0], (0.0, 0.7)), C.Tsit5(); reltol = 1e-12,
            abstol = 1e-12)
        @test sol2.u[end] ≈ collect(pool_state(pool, 100.0, 0.7; ZD0 = 30.0)) rtol = 1e-8
        eq = ZoosporePool(1.3, 1.3)
        @test pool_state(eq, 10.0, 2.0)[2] ≈ 10 * 1.3 * 2.0 * exp(-2.6)
        @test pool_state(ZoosporePool(1.3, 1.3 + 1e-12), 10.0, 2.0)[2] ≈
              pool_state(eq, 10.0, 2.0)[2] rtol = 1e-9
        # shedding and uptake
        @test pool_derivatives(pool, 2.0, 3.0; shedding = 5.0, uptake = 1.0) ==
              (5.0 - 1.6 * 2.0 - 1.0, 1.6 * 2.0 - 1.1 * 3.0)
        # what qPCR sees versus what infects
        @test qpcr_zoospores(ZL, ZD) == ZL + ZD
        @test infectious_zoospores(ZL, ZD) == ZL
        @test qpcr_zoospores(pool, 100.0, 0.7) ≈ ZL + ZD
        @test infectious_zoospores(pool, 100.0, 0.7) ≈ ZL
        @test viable_fraction(pool, 0.0) == 1
        @test viable_fraction(ZoosporePool(1.0, 3.0), 50.0) ≈ 1 - 1 / 3 rtol = 1e-6  # d > m: 1 − m/d
        @test viable_fraction(ZoosporePool(3.0, 1.0), 50.0) < 1e-6              # d < m: → 0
        @test_throws ArgumentError ZoosporePool(-1.0, 1.0)
        # Eq. 6 as printed, ZL0 m/(d − m)(e^{−dt} − e^{−mt}), is negative for every t > 0 (in
        # both the preprint and the published version): its bracket is reversed.
        m, d, t = 1.6, 1.1, 0.7
        printed = 100.0 * m / (d - m) * (exp(-d * t) - exp(-m * t))
        @test printed < 0 && printed ≈ -ZD
    end

    @testset "compendium and units" begin
        recs = parameter_records(:rumschlag2022)
        @test all(r -> r.verified, recs)
        @test all(r -> occursin("published version", r.notes), recs)   # version checked
        @test all(r -> r.units == "h", filter(r -> r.symbol in (:lifespan, :persistence), recs))
        for w in waters
            p = RumschlagParameters(:rumschlag2022; water = w)
            @test zoospore_lifespan(p) ≈ published_lifespan[w]
            @test p.m ≈ 24 / published_lifespan[w]                  # per day
            @test persistence_time(p) ≈ 21.6
            @test p.β == 7.1e-6
        end
        @test C.time_unit(RumschlagParameters(:rumschlag2022)) === :day
        @test_throws ArgumentError RumschlagParameters(:rumschlag2022; water = :Tap)
        # m in hours⁻¹ through the units layer: 1/14.4 h⁻¹
        r = only(filter(r -> r.symbol === :lifespan && r.set == "ASW", recs))
        @test in_time_unit(r, :day) ≈ 0.6
    end

    @testset "ODE (eqs 1-4) against the closed forms (eqs 5, 6, 8, 10)" begin
        p = RumschlagParameters(:rumschlag2022; water = :Flatwoods)
        sol = solve(ZoosporePersistenceProblem(p; ZL0 = 1e5, tspan = (0.0, 2.0));
            reltol = 1e-10, abstol = 1e-10)
        for t in (0.25, 1.0, 2.0)
            u = sol(t)
            @test [u[1], u[2]] ≈ collect(pool_state(ZoosporePool(p), 1e5, t)) rtol = 1e-6
            @test u[3] ≈ exp(p.β * 1e5 / p.m * (exp(-p.m * t) - 1)) rtol = 1e-6   # eq. 8
            @test u[3] + u[4] ≈ 1
        end
        # the infection assay: host added at lag t1, removed a day later (eq. 10)
        t1 = 0.5
        s = solve(ZoosporePersistenceProblem(p; ZL0 = 1e5 * exp(-p.m * t1), tspan = (0.0, 1.0));
            reltol = 1e-10, abstol = 1e-12)
        @test s.u[end][4] ≈ infection_probability(p, 1e5, t1) rtol = 1e-6
        @test C.to_ode_problem(ZoosporePersistenceProblem(p; ZL0 = 1.0)) isa C.ODEProblem
        @test_throws ArgumentError ZoosporePersistenceProblem(p; ZL0 = -1.0)
    end

    @testset "published curves (Fig. 4a) with the published values" begin
        # Fig. 4a at lag 0, read off the figure: ASW ≈ 0.29, Green Swamp ≈ 0.225; at 1 day
        # ASW ≈ 0.065 and Green Swamp ≈ 0.02. The printed β gives these only with ZL0 = 1e5.
        pA = RumschlagParameters(:rumschlag2022; water = :ASW)
        pG = RumschlagParameters(:rumschlag2022; water = :GreenSwamp)
        @test isapprox(infection_probability(pA, 1e5, 0.0), 0.29; atol = 0.01)
        @test isapprox(infection_probability(pG, 1e5, 0.0), 0.225; atol = 0.01)
        @test isapprox(infection_probability(pA, 1e5, 1.0), 0.065; atol = 0.01)
        @test isapprox(infection_probability(pG, 1e5, 1.0), 0.02; atol = 0.01)
        # With the fitted initial concentration (≈ 143 zoospores µl⁻¹, Fig. 3a) the printed β
        # gives essentially no infection. The paper does not say what scale its β is on.
        @test_broken isapprox(infection_probability(pA, 143.0, 0.0), 0.29; atol = 0.01)
        @test infection_probability(pA, 143.0, 0.0) < 1e-3
    end

    @testset "ratio to Mitchell et al. 2008" begin
        # Review §9.2 / wiki card: a four-week lifespan (672 h) is 47-71 times longer.
        @test round(Int, 672 / published_lifespan.ASW) == 47
        @test round(Int, 672 / published_lifespan.GreenSwamp) == 71
        # Mitchell's own value is μZ = 13 yr⁻¹, a lifespan of 673.8 h (365-day year): 47-72.
        μZ = parameter_set(:mitchell2008).μZ
        L = 1 / convert_rate(μZ, :year => :hour)
        @test round(L, digits = 1) == 673.8
        @test round(Int, L / published_lifespan.ASW) == 47
        @test_broken round(Int, L / published_lifespan.GreenSwamp) == 71   # gives 72
        @test round(Int, L / published_lifespan.GreenSwamp) == 72
    end

    @testset "deposited data" begin
        data = load_rumschlag2022()
        @test length(data.survival) == 199 && length(data.infection) == 197
        @test Set(c.water for c in data.survival) == Set(waters)
        @test count(h -> h.infected, data.infection) == 32
        @test maximum(c.t for c in data.survival) * 24 < 51
        nominal = load_rumschlag2022(; time = :nominal)
        @test Set(round(Int, 24c.t) for c in nominal.survival) == Set([0, 6, 16, 24, 48])
        @test Set(round(Int, 24h.lag) for h in data.infection) == Set([0, 3, 6, 12, 24])
        @test_throws ArgumentError load_rumschlag2022(; time = :wrong)
        @test_throws ArgumentError rumschlag2022_path(:other)
    end

    @testset "refit of the selected model (Results, p. 176)" begin
        model = RumschlagModel()                     # m by water, one d, one β
        @test C.parameter_names(model) == (:m_ASW, :m_TroutCreek, :m_Flatwoods,
            :m_GreenSwamp, :d, :β, :ZL0, :SD_L, :SD_D)
        @test C.nobs(model) == 2 * 199 + 197
        fit = maximum_likelihood(model; starts = 2)
        @test fit.df == 9
        @test maximum(fit.logliks) - minimum(fit.logliks) < 1e-3
        s = rumschlag_summary(fit)
        for w in waters
            @test round(s.lifespan[w][1], digits = 1) == published_lifespan[w]
            @test rumschlag_parameters(model, fit.x, w).m ≈ 24 / s.lifespan[w][1]
        end
        # "± 95% CI": half the width of the Wald interval of m, carried to 24/m hours
        @test round(s.lifespan.ASW[2], digits = 1) == 2.3
        @test round(s.lifespan.TroutCreek[2], digits = 1) == 2.2
        @test round(s.lifespan.Flatwoods[2], digits = 1) == 1.6
        @test round(s.lifespan.GreenSwamp[2], digits = 1) == 1.1
        @test round(s.persistence.ASW[1], digits = 1) == 21.6
        @test round(s.persistence.ASW[2], digits = 1) == 5.2
        # β on the paper's scale: β ZL0 / 1e5 (an inferred convention; see the compendium)
        @test round(s.β.ASW[1], sigdigits = 2) == 7.1e-6
        # Its 95% half-width is 2.6e-6 in the paper; the delta method on β ZL0 gives 2.54e-6.
        # Ignoring the β-ZL0 covariance (−0.28) would give 2.63e-6. Unexplained, small.
        @test_broken round(s.β.ASW[2], sigdigits = 2) == 2.6e-6
        @test round(s.β.ASW[2], sigdigits = 3) == 2.54e-6                       # regression
        @test s.β.ASW == s.β.GreenSwamp                       # one β for all waters
        # the fitted initial concentration matches Fig. 3a (curves start near 140-145 µl⁻¹)
        @test 135 < fit.estimates.ZL0 < 150
        # model AICc, Table 1: 1379.2 (ours 1379.0)
        @test isapprox(fit.aicc, 1379.2; atol = 0.5)
        # LRT against the AICc-best model (m and d by water): χ² = 7.581, df = 3, p = 0.056
        full_d = maximum_likelihood(RumschlagModel(; d = 4); starts = 2)
        χ2 = 2 * (full_d.loglik - fit.loglik)
        @test isapprox(χ2, 7.581; atol = 1e-3)                  # ours 7.5815
        @test round(ccdf(Chisq(3), χ2), digits = 3) == 0.056
        # Figure 4a with the refit, which needs no scale convention for β
        pA = rumschlag_parameters(model, fit.x, :ASW)
        @test isapprox(infection_probability(pA, fit.estimates.ZL0, 0.0), 0.29; atol = 0.01)
        # with the nominal sampling times the published lifespans are not reproduced
        nom = maximum_likelihood(RumschlagModel(; data = load_rumschlag2022(; time = :nominal));
            starts = 1)
        @test round(24 / nom.estimates.m_ASW, digits = 1) != 14.4
    end

    @testset "model selection (Table 1)" begin
        sel = rumschlag_model_selection(; starts = 1)
        @test length(sel) == 15
        @test length(unique(f.k for f in sel)) == 15
        key(m, d, β) = (m = m, d = d, β = β)
        # Table 1 (published): (estimates of m, d, β) => (df, AICc)
        table1 = Dict(key(4, 4, 1) => (12, 1377.9), key(4, 4, 4) => (15, 1379.0),
            key(4, 1, 1) => (9, 1379.2), key(4, 1, 4) => (12, 1380.2),
            key(2, 2, 1) => (8, 1385.0), key(2, 1, 1) => (7, 1389.1),
            key(2, 2, 2) => (9, 1389.7), key(1, 2, 1) => (7, 1391.8),
            key(1, 2, 2) => (8, 1392.0), key(2, 1, 2) => (8, 1393.9),
            key(1, 4, 1) => (9, 1395.0), key(1, 1, 1) => (6, 1396.5),
            key(1, 1, 2) => (7, 1396.7), key(1, 4, 4) => (12, 1397.2),
            key(1, 1, 4) => (9, 1398.7))
        @test Set(keys(table1)) == Set(f.k for f in sel)
        for f in sel
            df, aicc = table1[f.k]
            @test f.df == df
            if f.k in (key(2, 2, 2), key(2, 1, 2))
                # Ours is 3.3 AICc units LOWER (a better optimum) for both models in which m
                # and β have two estimates; same refit from several starts. Unexplained.
                @test_broken isapprox(f.aicc, aicc; atol = 0.5)
                @test f.aicc < aicc - 3
            else
                @test isapprox(f.aicc, aicc; atol = 0.5)      # ours 0.0-0.4 lower
            end
        end
        # ranking and weights of the top four (Table 1: 0.40, 0.24, 0.22, 0.13)
        @test [f.k for f in sel[1:4]] ==
              [key(4, 4, 1), key(4, 4, 4), key(4, 1, 1), key(4, 1, 4)]
        @test [f.weight for f in sel[1:4]] ≈ [0.40, 0.24, 0.22, 0.13] atol = 0.025
        @test all(f -> f.Δaicc > 2, sel[5:end])
        # the selected model is the simplest within 2 AICc of the best
        within = filter(f -> f.Δaicc < 2, sel)
        @test argmin(f -> f.df, within).k == key(4, 1, 1)
    end
end
