# Bsal transmission functions: Malagon et al. 2020 Sci Rep 10:5584
# (doi:10.1038/s41598-020-62351-x) and Tompros et al. 2022 Transbound Emerg Dis 69:731
# (doi:10.1111/tbed.14043). Oracles were read from the published articles, Malagon's
# Supplementary Information (Table S1) and Tompros's appendix; the data are the authors'
# deposits (Malagon Supplementary File 2; Tompros TRACE doi:10.7290/xw8yy2RBeu).
# Wiki: lab/wiki/models/bsal-transmission-functions.md
using Chytridiomycosis
using Test
using LogDensityProblems: LogDensityProblems

@testset "Bsal transmission functions" begin
    C = Chytridiomycosis

    @testset "transmission functions" begin
        S, I = 6.0, 2.0
        N = S + I
        @test infection_rate(DensityDependent(0.1), S, I) ≈ 0.1 * S * I
        @test infection_rate(FrequencyDependent(0.23), S, I) ≈ 0.23 * S * I / N
        @test infection_rate(FrequencyDependent(0.23), 6, 2) ≈ 0.345
        @test infection_rate(PowerLaw(0.2, 1.0, 1.0), S, I) ≈ infection_rate(DensityDependent(0.2), S, I)
        # special cases of quadratic saturation (Tompros appendix, below eq. 4)
        @test infection_rate(QuadraticSaturation(10.0, 0.0, 0.0), S, I) ≈
              infection_rate(DensityDependent(0.1), S, I)
        @test infection_rate(QuadraticSaturation(0.0, 1 / 0.23, 0.0), S, I) ≈
              infection_rate(FrequencyDependent(0.23), S, I)
        @test infection_rate(StructuralSaturation(1.0, 2.0, 2.0, 0.0), S, I) ≈
              infection_rate(QuadraticSaturation(1.0, 2.0, 0.0), S, I)
        @test infection_rate(PrevalenceSaturation(1 / 3.39, 0.0, 0.0), S, I) ≈
              infection_rate(PrevalenceDependent(3.39), S, I)
        @test infection_rate(PrevalencePowerLaw(3.0, 1.0, 1.0), S, I) ≈
              infection_rate(PrevalenceDependent(3.0), S, I)
        # Malagon: β f(N) S I / N with f(N) = cN/(K + N)
        h = HollingContact(0.95, 2153.29, 79.55)
        @test infection_rate(h, S, I) ≈ 0.95 * contact_rate(h, N) * S * I / N
        @test contact_rate(h, 79.55) ≈ 2153.29 / 2
        # no infection without susceptibles or infecteds
        for tf in (DensityDependent(1.0), FrequencyDependent(1.0), PowerLaw(1.0, 0.4, 0.4),
                   PrevalencePowerLaw(1.0, 0.4, 0.4), StructuralSaturation(1.0, 1.0, 1.0, 1.0))
            @test infection_rate(tf, 0, 3) == 0
            @test infection_rate(tf, 3, 0) == 0
        end
        @test force_of_infection(FrequencyDependent(0.23), S, I) ≈ 0.23 * I / N
        # frequency dependence: contact rate does not change with N; density dependence does
        @test transmission_coefficient(FrequencyDependent(0.23), 2) ==
              transmission_coefficient(FrequencyDependent(0.23), 32)
        @test transmission_coefficient(DensityDependent(0.1), 30) ≈ 3.0
        @test transmission_coefficient(h, 133) ≈ 0.95 * contact_rate(h, 133)
        @test transmission_coefficient(PowerLaw(0.2, 0.4, 0.4), 10) == Inf
        @test transmission_coefficient(QuadraticSaturation(1.0, 2.0, 3.0), 4) ≈ 4 / (1 + 8 + 48)
        # derivative check of the coefficient against φ at small I
        for tf in (QuadraticSaturation(0.5, 2.0, 0.1), StructuralSaturation(0.3, 2.0, 3.0, 0.4),
                   PrevalenceSaturation(0.3, 0.2, 0.1), HollingContact(0.9, 50.0, 20.0))
            ε = 1e-7
            @test transmission_coefficient(tf, 10.0) ≈ infection_rate(tf, 10.0, ε, 10.0) / ε rtol = 1e-5
        end
    end

    @testset "plugs into the Schmidt S–L–I model" begin
        p = SchmidtParameters(:schmidt2017)
        a = solve(SchmidtSLIProblem(p; tspan = (0.0, 0.5)); abstol = 1e-10, reltol = 1e-10)
        b = solve(TransmissionSLIProblem(p, DensityDependent(p.β); tspan = (0.0, 0.5));
            abstol = 1e-10, reltol = 1e-10)
        @test a.u[end] ≈ b.u[end] rtol = 1e-8
        # frequency dependence with the same initial per capita rate (β K) is the same at the start
        fd = TransmissionSLIProblem(p, FrequencyDependent(p.β * p.K); tspan = (0.0, 0.5))
        du1, du2 = zeros(3), zeros(3)
        C.schmidt_rhs!(du1, a.u[1], p, 0.0)
        C.transmission_sli_rhs!(du2, fd.u0, (p, fd.tf), 0.0)
        @test du1 ≈ du2 rtol = 1e-3
        @test C.SciMLBase.successful_retcode(solve(fd))
    end

    @testset "master equation (Tompros appendix eq. 2)" begin
        tf = FrequencyDependent(0.23)
        for (S0, I0, t) in ((28, 4, 2.0), (5, 5, 3.0), (16, 16, 2.0), (1, 1, 5.0))
            u = susceptible_distribution(tf, S0, I0, t; method = :uniformization)
            e = susceptible_distribution(tf, S0, I0, t; method = :exp)
            sol = solve(SusceptibleChainProblem(tf, S0, I0; tspan = (0.0, t)); abstol = 1e-13,
                reltol = 1e-11)
            @test sum(u) ≈ 1 atol = 1e-12
            @test u ≈ e atol = 1e-12
            @test sol.u[end] ≈ u atol = 1e-9
        end
        # one susceptible and one infected: p_0 = 1 − exp(−β t / 2)
        p = susceptible_distribution(tf, 1, 1, 5.0)
        @test p[1] ≈ -expm1(-0.23 * 5 / 2)
        @test expected_susceptibles(tf, 1, 1, 5.0) ≈ exp(-0.23 * 5 / 2)
        @test susceptible_distribution(tf, 3, 0, 5.0) ≈ [0.0, 0.0, 0.0, 1.0]   # no infecteds
        @test susceptible_distribution(tf, 0, 4, 5.0) == [1.0]
        # large rates fall back to the matrix exponential
        big = susceptible_distribution(DensityDependent(5.0), 20, 20, 3.0)
        @test big[1] ≈ 1 atol = 1e-8
    end

    @testset "Tompros 2022: data and preprocessing" begin
        data = load_tompros2022()
        @test length(data) == 20
        @test sum(t.S[1] + t.I[1] for t in data) == 290            # n = 290 newts
        @test sum(t.I[1] for t in data) == 101                     # Table 1 interior sums to 101
        @test all(t -> t.days == [0.0, 2.0, 5.0], data)
        @test all(t -> t.S[1] + t.I[1] == t.S[end] + t.I[end], data)   # no deaths by day 5
        # the false-negative / false-positive rules (section 2.3): NNPNPP → NNPPPP
        @test clean_infection_series(Bool[0, 0, 1, 0, 1, 1], falses(6)) == Bool[0, 0, 1, 1, 1, 1]
        @test clean_infection_series(Bool[0, 1, 0, 0], falses(4)) == Bool[0, 0, 0, 0]
        # rebuilding the counts from the swabs gives the authors' processed data
        st = tompros_status()
        @test length(st) == 290
        rebuilt = tompros_counts(st)
        @test [(t.tank, t.S, t.I) for t in rebuilt] == [(t.tank, t.S, t.I) for t in data]
        val = load_tompros2022(; data = :validation)
        @test [t.tank for t in val] == [1, 2, 4, 7, 8, 10, 11, 12, 19, 20]
        r58 = tompros_counts(st; days = (5, 8))
        nodeath = filter(t -> t.S[1] + t.I[1] == t.S[2] + t.I[2], r58)
        @test [(t.tank, t.S, t.I) for t in nodeath] == [(t.tank, t.S, t.I) for t in val]
        # ">90% of susceptible newts tested positive within 17 days"; "all positive by day 59"
        sus = [s for (_, s) in st if !s.infected_start]
        @test length(sus) == 189
        firstpos(s) = (j = findfirst(s.positive); j === nothing ? Inf : s.days[j])
        @test count(s -> firstpos(s) <= 17, sus) / length(sus) > 0.9
        @test all(s -> firstpos(s) <= 59, values(st))
        # "Total mortality among treatments ... was 45%": the deposit gives 81/189 = 42.9%
        # of the newts exposed by contact and 156/290 = 53.8% of all newts.
        @test_broken round(count(s -> s.died, sus) / length(sus), digits = 2) == 0.45
        @test count(s -> s.died, sus) == 81
        @test count(s -> s.died, values(st)) == 156
    end

    @testset "Tompros 2022: Table 2 refit" begin
        # the printed estimates give the printed negative log-likelihoods
        pub = (frequency = (63.2, :frequency), quadratic = (62.2, :quadratic),
               power = (63.2, :power), prevalence = (68.7, :prevalence),
               prevalence_power = (68.4, :prevalence_power),
               prevalence_quadratic = (68.7, :prevalence_quadratic))
        for (form, (nll, set)) in pairs(pub)
            m = TomprosModel(form)
            θ = [getproperty(parameter_set(:tompros2022; set = string(set)), n)
                 for n in parameter_names(m)]
            @test round(tompros_negloglik(form, θ), digits = 1) == nll
        end
        # structural saturation: printed estimates give 62.43, printed NLL 62.5 (AIC 132.9
        # implies 62.45 ± 0.025)
        m = TomprosModel(:structural)
        θ = [getproperty(parameter_set(:tompros2022; set = "structural"), n)
             for n in parameter_names(m)]
        @test_broken round(tompros_negloglik(:structural, θ), digits = 1) == 62.5
        @test round(tompros_negloglik(:structural, θ), digits = 2) == 62.43

        sel = tompros_model_selection()
        aic = Dict(f.form => f.aic for f in sel)
        @test first(sel).form === :frequency                        # frequency dependence best
        @test round(aic[:frequency], digits = 1) == 128.3
        @test round(aic[:quadratic], digits = 1) == 130.3
        @test round(aic[:power], digits = 1) == 132.4
        @test round(aic[:structural], digits = 1) == 132.9
        @test round(aic[:prevalence], digits = 1) == 139.3
        @test round(aic[:prevalence_power], digits = 1) == 142.8
        @test round(aic[:prevalence_quadratic], digits = 1) == 143.3
        @test [f.form for f in sel][1:7] == [:frequency, :quadratic, :power, :structural,
            :prevalence, :prevalence_power, :prevalence_quadratic]          # Table 2 order
        # density dependence "failed to converge" in the paper; here it has an MLE and is last
        dd = only(filter(f -> f.form === :density, sel))
        @test last(sel).form === :density
        @test dd.Δaic > 20
        est(f) = only(filter(s -> s.form === f, sel)).estimates
        @test round(est(:frequency).β, digits = 3) == 0.230
        @test round(est(:quadratic).κ₁, digits = 2) == 2.79
        @test round(est(:quadratic).κ₂, digits = 3) == 0.085
        # κ₀ lies on a nearly flat profile (Fig. 4: [0, 34.2]); 0.0227 here, 0.026 printed
        @test_broken round(est(:quadratic).κ₀, digits = 3) == 0.026
        @test round(est(:quadratic).κ₀, digits = 3) == 0.023
        @test round(est(:power).a, digits = 3) == 0.362
        @test round(est(:power).b, digits = 3) == 0.382
        @test round(est(:structural).κ_S, digits = 2) == 2.55
        @test est(:structural).κ₀ == 1e-6                         # at the bound, as printed
        @test round(est(:prevalence).β, digits = 2) == 3.39
        @test round(est(:prevalence_power).β, digits = 2) == 4.29
        @test est(:prevalence_quadratic).κ₁ == 1e-6
        # all starts agree
        @test all(f -> maximum(f.nlls) - f.nll < 1e-3, sel)

        # Fig. 3: 95% of the normalised likelihood in β ∈ [0.181, 0.293]
        ci = likelihood_interval(TomprosModel(:frequency), range(0.0001, 0.5; length = 301))
        @test round(ci.lower, digits = 3) == 0.181
        @test round(ci.upper, digits = 3) == 0.293
        @test_throws ArgumentError likelihood_interval(TomprosModel(:power), 0:0.1:1)

        # section 3.3: "AIC" over days 5 to 8 for the 10 tanks with no deaths
        val = load_tompros2022(; data = :validation)
        q = parameter_set(:tompros2022; set = "quadratic")
        nfd = tompros_negloglik(:frequency, [parameter_set(:tompros2022).β]; data = val)
        nq = tompros_negloglik(:quadratic, [q.κ₀, q.κ₁, q.κ₂]; data = val)
        @test round(2nfd + 2, digits = 1) == 65.9
        @test round(2nq + 6, digits = 1) == 70.6

        # LogDensityProblems interface
        m = TomprosModel(:frequency)
        @test LogDensityProblems.dimension(m) == 1
        @test LogDensityProblems.logdensity(m, [log(0.23)]) ≈ -tompros_negloglik(:frequency, [0.23])
        @test LogDensityProblems.logdensity(m, [NaN]) == -Inf
        @test transmission_function(m, [0.3]) == FrequencyDependent(0.3)
        @test_throws ArgumentError TomprosModel(:nonsense)
    end

    @testset "Tompros 2022: R₀" begin
        p = TomprosParameters(:tompros2022)
        @test p.σ ≈ 0.025
        @test p.μ == 0.0223
        @test round(basic_reproduction_number(p), digits = 1) == 4.9
        @test basic_reproduction_number(FrequencyDependent(p), 10; removal = p.μ + p.σ) ≈
              basic_reproduction_number(p)
        @test C.time_unit(p) === :day
    end

    @testset "Malagon 2020: contacts (Experiment 1)" begin
        data = load_malagon2020()
        @test length(data.contacts) == 672
        cm = contact_means(data.contacts)
        mean_of(h, n) = only(filter(r -> r.habitat === h && r.density == n, cm)).mean
        # the inputs of the authors' Holling fit are the back-transformed means
        for r in data.fit_inputs
            r.density == 133 && continue
            digits = r.habitat === :simple ? 1 : 6
            @test round(mean_of(r.habitat, Int(r.density)); digits = digits) == r.total
        end
        # "plants significantly reduced newt contacts by 3X" at 8 newts; 35% at 4 newts
        @test round(mean_of(:simple, 8) / mean_of(:complex, 8)) == 3
        @test abs(100 * (1 - mean_of(:complex, 4) / mean_of(:simple, 4)) - 35) < 1
        # "per capita contact rates at 8 newts per m² were 5-15X greater than 2 newts per m²":
        # per capita ratios are 4.7 (no plants) and 2.1 (plants); totals give 18.8 and 8.3
        pc(h) = (mean_of(h, 8) / 8) / (mean_of(h, 2) / 2)
        @test_broken 5 <= pc(:complex) && pc(:simple) <= 15
        @test round(pc(:simple), digits = 1) == 4.7 && round(pc(:complex), digits = 1) == 2.1
        # ANOVA F statistics (Results)
        a = contact_anova(data.contacts)
        @test a.density_simple.df == (2, 81) && round(a.density_simple.F, digits = 1) == 199.3
        @test a.density_complex.df == (2, 81) && round(a.density_complex.F, digits = 1) == 59.8
        @test a.habitat_8.df == (1, 58) && round(a.habitat_8.F, digits = 1) == 80.3
        @test a.habitat_4.df == (1, 52) && round(a.habitat_4.F, digits = 2) == 4.89
    end

    @testset "Malagon 2020: Table S1" begin
        tab = malagon_table_s1()
        rec(sym, set) = only(filter(r -> r.symbol === sym && r.set == set,
            parameter_records(:malagon2020))).value
        for h in (:simple, :complex)
            row = first(filter(r -> r.habitat === h, tab))
            @test round(row.c, digits = 2) == rec(:c, string(h))
            @test round(row.K, digits = 2) == rec(:K, string(h))
        end
        @test length(tab) == 12
        for r in tab
            set = "$(r.habitat)_day$(r.disease_state)_$(r.cohabitation)min"
            if set == "complex_day18_10min"
                # printed β = 0.93; the sum of squares falls monotonically to the grid edge,
                # 0.98, where it is 0.240 (the printed SSR). At β = 0.93 the best SSR is 0.257.
                @test_broken r.β == rec(:β, set)
                @test r.β == 0.98
                g = fit_cohabitation(r.c, r.K, [0.2, 0.8, 1.0, 1.0])
                @test round(minimum(g.grid[:, findfirst(==(0.93), 0:0.01:0.98)]), digits = 2) == 0.26
            else
                @test r.β == rec(:β, set)
            end
            @test r.γ == rec(:γ, set)
        end
        # SSR column of Table S1
        @test [round(r.ssr, digits = 2) for r in tab] ==
              [0.04, 0.02, 0.11, 0.03, 0.07, 0.07, 0.07, 0.08, 0.24, 0.16, 0.23, 0.23]
        # the Holling fit as the code does it: c ≤ 1.8e5 per hour, K ≤ 100
        inp = filter(r -> r.habitat === :simple, load_malagon2020().fit_inputs)
        f = fit_contact_function([r.density for r in inp], [r.total for r in inp])
        @test f.c ≈ 24 * f.c_hour
        @test round(f.c_hour, digits = 1) == 89.7                  # contacts newt⁻¹ h⁻¹
    end

    @testset "Malagon 2020: cohabitation model and Fig. 6" begin
        p = MalagonParameters(:malagon2020)
        @test (p.β, p.γ, p.c, p.K) == (0.95, 0.11, 2153.29, 79.55)
        @test p.d == 0.0223
        @test MalagonParameters(:malagon2020; mortality = :code).d ≈ 1 / 42.75
        @test MalagonParameters(:malagon2020; habitat = :complex, disease_state = 18,
            cohabitation = 30).γ == 0.21
        @test_throws ArgumentError MalagonParameters(:malagon2020; habitat = :marsh)
        @test C.time_unit(p) === :day
        # the closed form agrees with the ODE
        for h in (:simple, :complex)
            prob = CohabitationProblem(MalagonParameters(:malagon2020; habitat = h))
            sol = solve(prob; abstol = 1e-12, reltol = 1e-10)
            cf = cohabitation_curve(prob, [5.0, 30.0, 96.0])
            S0 = prob.u0[1]
            @test [sol(t)[4] / S0 for t in (5.0, 30.0, 96.0)] ≈ cf.infectious rtol = 1e-6
            @test [sol(t)[5] / S0 for t in (5.0, 30.0, 96.0)] ≈ cf.dead rtol = 1e-6
        end
        simple = CohabitationProblem(MalagonParameters(:malagon2020; habitat = :simple))
        complex = CohabitationProblem(MalagonParameters(:malagon2020; habitat = :complex))
        # ">95% of newts could become infected after one month" (simple habitat)
        @test cohabitation_curve(simple, 30.0).infectious > 0.95
        # "with >80% mortality in three months"
        @test cohabitation_curve(simple, 90.0).dead > 0.8
        # complex habitat: "infection prevalence was <80% and cumulative mortality was <70%"
        @test cohabitation_curve(complex, 96.0).infectious < 0.8
        @test cohabitation_curve(complex, 96.0).dead < 0.7
        # R₀ = 3.2 and 1.9 from I_simple = 0.95 and I_complex = 0.77 (Anderson & May eq. 6.21)
        @test round(final_size_reproduction_number(0.95), digits = 1) == 3.2
        @test round(final_size_reproduction_number(0.77), digits = 1) == 1.9
        # From the model: complex reaches its asymptote (0.761 → R₀ 1.9). The simple habitat's
        # 0.95 is the value at about one month; its asymptote is 0.988, which gives R₀ 4.5.
        zc = cohabitation_curve(complex, 96.0).infectious
        @test round(final_size_reproduction_number(zc), digits = 1) == 1.9
        zs30 = cohabitation_curve(simple, 30.0).infectious
        @test round(final_size_reproduction_number(zs30), digits = 1) == 3.2
        zs = cohabitation_curve(simple, 96.0).infectious
        @test_broken round(final_size_reproduction_number(zs), digits = 1) == 3.2
        @test round(final_size_reproduction_number(zs), digits = 1) == 4.5
        # Experiment 2 survival: 89% mortality (40/45); median survival 32 days
        surv = filter(r -> r.group === :susceptible, load_malagon2020().survival)
        @test length(surv) == 45 && count(r -> r.dead, surv) == 40
        @test sort([r.days for r in surv])[23] == 32
    end
end
