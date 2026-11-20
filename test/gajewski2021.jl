# Gajewski et al. 2021 Ecol Evol 11:17920 (doi:10.1002/ece3.8379): delayed logistic growth of
# Bd optical density with five thermal performance curves. Oracles:
# - the curve and temperature-regime formulas of the authors' code (Zenodo
#   doi:10.5281/zenodo.4525954: TPCurves.R, functions.R, HierarchicalModel_*.bug), written out
#   again here independently of the package;
# - closed forms (Eq. 1 against the ODE of Eq. 3; the Stinner adjusted limits);
# - the printed posterior medians (Appendix S1 Tables B.3-B.7, Table 2) and the printed curve
#   summaries (Table 1: adjusted limits and maximum r; Table 2: Topt; Table B.2: median r over
#   each regime);
# - with the authors' data (not bundled; skipped without it): the penalized deviances of
#   Table 1 and the posterior medians, against the posterior mode of our log density.
using Chytridiomycosis
using Test
using LogDensityProblems

@testset "Gajewski 2021 Bd growth and temperature" begin
    C = Chytridiomycosis
    pset(f) = parameter_set(:gajewski2021; set = string(f))

    @testset "curve forms (TPCurves.R)" begin
        θ = pset(:stinner)
        s = gajewski_curve(:stinner, θ)
        @test s isa Stinner
        for T in (5.0, 12.0, 19.0, 22.5)
            @test s(T) ≈ (T < θ.Topt ? θ.C / (1 + exp(θ.k1 + θ.k2 * T)) :
                                      θ.C / (1 + exp(θ.k1 + θ.k2 * (2θ.Topt - T))))
        end
        @test s(θ.Topt - 3) ≈ s(θ.Topt + 3)              # symmetric about Topt
        # Stinner adjusted limits in closed form: C/(1 + exp(k1 + k2 T)) = 0.01
        lo = (log(θ.C / 0.01 - 1) - θ.k1) / θ.k2
        @test adjusted_limits(s)[1] ≈ lo atol = 1e-6
        @test adjusted_limits(s)[2] ≈ 2θ.Topt - lo atol = 1e-6

        θ = pset(:logan10)
        l = gajewski_curve(:logan10, θ)
        T = 20.0
        @test l(T) ≈ θ.α * (1 / (1 + θ.c * exp(-θ.b * T)) - exp(-(θ.Tmax - T) / θ.ΔT))
        @test l(θ.Tmax) == 0 && l(θ.Tmax + 1) == 0
        @test l(θ.Tmax - 1e-4) < 0                        # unclipped, as coded

        θ = pset(:briere2)
        b = gajewski_curve(:briere2, θ)
        @test b isa Briere
        @test b(20.0) ≈ θ.c * 20 * (20 - θ.Tmin) * abs(θ.Tmax - 20)^(1 / θ.b)
        @test b(θ.Tmin - 1) == 0 && b(θ.Tmax + 1) == 0

        θ = pset(:ratkowsky)
        x = 20.0 + 273.15
        @test gajewski_curve(:ratkowsky, θ)(20.0) ≈
              (θ.c * (x - θ.Tmin) * (1 - exp(θ.k * (x - θ.Tmax))))^2
        x = 20.0 + 273.16
        @test gajewski_curve(:ratkowsky, θ; kelvin_offset = 273.16)(20.0) ≈
              (θ.c * (x - θ.Tmin) * (1 - exp(θ.k * (x - θ.Tmax))))^2

        θ = merge(pset(:ikemoto), (R = parameter_set(:gajewski2021; set = "prior_ikemoto").R,))
        ik = gajewski_curve(:ikemoto, θ)
        x = 18.0 + 273.15
        @test ik(18.0) ≈ θ.ρ * (x / θ.TΦ) * exp(θ.ΔHA / θ.R * (1 / θ.TΦ - 1 / x)) /
              (1 + exp(θ.ΔHL / θ.R * (1 / θ.TL - 1 / x)) + exp(θ.ΔHH / θ.R * (1 / θ.TH - 1 / x)))
        @test θ.R == 1.987
        @test_throws ArgumentError gajewski_curve(:cubic, θ)
    end

    @testset "Eq. 1 is the solution of Eq. 3 at constant temperature" begin
        p = GajewskiParameters(:gajewski2021; form = :briere2)
        @test p.Y0 == 1e-3 && p.K == pset(:briere2).K && p.d == pset(:briere2).d
        @test time_unit(typeof(p)) === :day
        sol = solve(GajewskiGrowthProblem(p; temperature = 17.0); reltol = 1e-10, abstol = 1e-12)
        for t in (0.5, 1.0, 3.0, 7.5, 14.0)
            @test sol(t)[1] ≈ optical_density(p, t, 17.0) rtol = 1e-6
        end
        @test optical_density(p, 0.5, 17.0) == p.Y0              # before the delay
        @test delayed_logistic(0.0, 0.01, 1.0, 0.5) == 0.01
        @test delayed_logistic(1e3, 0.01, 1.0, 0.5) ≈ 1.0
        # the delay of the Ikemoto fit: d printed, exp(d) in the JAGS model
        @test GajewskiParameters(:gajewski2021; form = :ikemoto).d == pset(:ikemoto).d
        @test GajewskiParameters(:gajewski2021; form = :ikemoto, ikemoto_delay = :fit).d ≈
              exp(pset(:ikemoto).d)
        @test GajewskiParameters(:gajewski2021; form = :logan10, α = 0.9).tpc.curve.α == 0.9
        @test_throws ArgumentError GajewskiParameters(:gajewski2021; form = :cubic)
    end

    @testset "temperature regimes (functions.R)" begin
        dla = GajewskiRegime(:DLA)
        # temperature.dla() written out
        function dla_R(t)
            t -= trunc(t)
            t < 0.03848604 && return 18.5
            t < 0.2530747 && return 18.858696 - 9.320158t
            t < 0.3321566 && return 16.5
            t < 0.6361821 && return 11.58362 + 14.80139t
            t < 0.8402777 && return 21.0
            t < 0.8844696 && return 68.53571 - 56.57143t
            return 18.5
        end
        for t in 0.013:0.0371:3.0
            @test dla(t) ≈ dla_R(t)
        end
        @test dla(0.3) == dla(5.3)
        wlf = GajewskiRegime(:WLF)
        @test wlf(0.9) == 25.0 && wlf(0.6) == 28.0
        slip = GajewskiRegime(:WLF; as_coded_kelvin = true)
        @test slip(0.9) ≈ 25 - 273.16 && slip(0.6) == 28.0
        @test GajewskiRegime(:DLA; as_coded_kelvin = true)(0.95) == 18.5
        @test extrema(GajewskiRegime(:WHF).(0:0.001:1)) == (20.0, 21.25)
        @test_throws ArgumentError GajewskiRegime(:XYZ)
    end

    @testset "printed curve summaries from the posterior medians" begin
        # The paper's values are posterior medians of each summary; ours are the summaries of
        # the curve at the posterior medians. Same tolerances as the validation targets.
        for f in (:briere2, :ratkowsky, :ikemoto, :stinner)
            s, p = gajewski_summary(f), pset(f)
            @test s.r_max ≈ p.r_max rtol = 0.01
            @test s.Tmin_adj ≈ p.Tmin_adj atol = 0.2
            @test s.Tmax_adj ≈ p.Tmax_adj atol = 0.2
            @test s.Topt ≈ p.Topt atol = 0.2
            @test s.r_DLA ≈ p.r_DLA rtol = 0.01
            @test s.r_WHF ≈ p.r_WHF rtol = 0.01
            @test s.r_WLF ≈ p.r_WLF rtol = 0.01
        end
        s, p = gajewski_summary(:stinner), pset(:stinner)
        @test s.Tmin_adj ≈ p.Tmin_adj atol = 0.005              # 10.95, two decimals
        @test s.Tmax_adj ≈ p.Tmax_adj atol = 0.005              # 27.52
        @test s.Topt ≈ p.Topt atol = 1e-6
        @test round(s.Topt; digits = 2) == 19.23                   # Table 2
        # Brière-2 optimum from its closed form (the derivative of the curve vanishes)
        b = GajewskiParameters(:gajewski2021; form = :briere2).tpc.curve
        To = gajewski_summary(:briere2).Topt
        h = 1e-5
        @test (b(To + h) - b(To - h)) / 2h ≈ 0 atol = 1e-6
        # Logan-10: the printed α (0.134) cannot give the printed maximum r (0.852); the
        # midpoint of its printed HDI can
        l = gajewski_summary(:logan10)
        @test l.r_max < 0.2
        lo, hi = C._gajewski_interval("logan10", :α)
        lm = gajewski_summary(:logan10; α = (lo + hi) / 2)
        @test lm.r_max ≈ pset(:logan10).r_max rtol = 0.01
        @test lm.Tmin_adj ≈ pset(:logan10).Tmin_adj atol = 0.2
        @test lm.Tmax_adj ≈ pset(:logan10).Tmax_adj atol = 0.2
        @test lm.Topt ≈ pset(:logan10).Topt atol = 0.2
        # the WLF slip makes a difference to the kelvin curves only
        @test gajewski_summary(:ratkowsky).r_WLF != gajewski_summary(:ratkowsky; as_coded = false).r_WLF
        @test gajewski_summary(:briere2).r_WLF == gajewski_summary(:briere2; as_coded = false).r_WLF
    end

    @testset "varying-temperature predictions at the medians" begin
        for f in (:briere2, :ratkowsky, :ikemoto, :stinner)
            sol = gajewski_prediction(f, :DLA)
            K = GajewskiParameters(:gajewski2021; form = f).K
            @test sol(0.0)[1] == 1e-3
            @test 0.85K < sol(14.0)[1] <= K * (1 + 1e-6)      # near carrying capacity
            @test issorted([sol(t)[1] for t in 0:0.25:14])
        end
        # the WLF regime is hotter: every curve grows more slowly there than under DLA
        for f in (:briere2, :ratkowsky, :ikemoto, :stinner)
            @test gajewski_prediction(f, :WLF)(4.0)[1] < gajewski_prediction(f, :DLA)(4.0)[1]
        end
        # without the slip, Ratkowsky grows during the last 3.3 h of each WLF day
        @test gajewski_prediction(:ratkowsky, :WLF; as_coded = false)(4.0)[1] >
              gajewski_prediction(:ratkowsky, :WLF)(4.0)[1]
    end

    data_found = C._gajewski2021_find() !== nothing
    data_found || @info "Gajewski et al. 2021 data not found; skipping the data-based tests"
    data_found && @testset "posterior and the authors' data" begin
        d = load_gajewski2021()
        @test length(d.constant.y) == 4693
        @test sort(unique(d.constant.T)) == [13.0, 15, 17, 19, 21, 23, 25, 26, 27, 28]
        @test sum(length(v.y) for v in values(d.varying)) == 1911
        for f in C.GAJEWSKI_FORMS
            m = GajewskiModel(; form = f, data = d.constant)
            x0 = gajewski_start(m)
            @test LogDensityProblems.dimension(m) == length(x0) == length(parameter_names(m))
            @test isfinite(LogDensityProblems.logdensity(m, x0))
            bad = copy(x0); bad[2] = 0.9                    # K outside its prior
            @test LogDensityProblems.logdensity(m, bad) == -Inf
            fit = posterior_mode(m; rounds = 2, iters = 2500)
            k = length(x0)
            # DIC = D(θ̄) + 2 pD ≈ −2 log L(mode) + 2k: Table 1's penalized deviances
            dic = -2 * C._gajewski_loglik(m, fit.estimates) + 2k
            @test dic ≈ pset(f).penalized_deviance atol = 15
        end
        # the Logan-10 α: the mode lies in the printed HDI, not at the printed median
        fit = posterior_mode(GajewskiModel(; form = :logan10, data = d.constant))
        lo, hi = C._gajewski_interval("logan10", :α)
        @test lo < fit.estimates.α < hi
        # Stinner: the mode is at the printed medians
        fit = posterior_mode(GajewskiModel(; form = :stinner, data = d.constant))
        for n in (:C, :k1, :k2, :Topt, :K, :d)
            @test fit.estimates[n] ≈ pset(:stinner)[n] rtol = 0.02
        end
        # the as-printed likelihood is a different model (lognormal for every form)
        m = GajewskiModel(; form = :stinner, likelihood = :lognormal, data = d.constant)
        @test isfinite(LogDensityProblems.logdensity(m, gajewski_start(m)))
        # Table B.1 (logistic fits to each regime, no TPC) is NOT reproduced: the mode of the
        # authors' model on the deposited data is far from the printed medians
        m = GajewskiLogisticModel(:DLA; data = d.varying)
        fit = posterior_mode(m)
        @test LogDensityProblems.logdensity(m, fit.x) >
              LogDensityProblems.logdensity(m, gajewski_start(m)) + 100
        @test !isapprox(fit.estimates.r, parameter_set(:gajewski2021; set = "logistic_DLA").r;
                        rtol = 0.2)
    end
end
