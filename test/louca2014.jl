# Louca, Lampo & Doebeli 2014 (Proc R Soc B 281:20132783): DOI-distribution model of adults
# and tadpoles, the basic reproduction quotient Q_o and a Monte Carlo extinction-risk
# assessment. The oracles are:
# - the paper (§2–3, Tables 1–2, Figs 2–5) and its supplement (ESM §1–5, Table S1, eqs 1–31);
# - closed forms that the ESM states (r_br integrates to r, g_s(T_s) halves survival, I′(0) = λ,
#   I(K_S) = 0, eq. 9 puts the disease-free equilibrium at (K_A, K_T), sign(ρ) = sign(Q_o − 1));
# - the ESM §4–§5.1 arithmetic (μ_S, T, β̃, β_A, η, μ_Z), from compendium inputs;
# - the stored output of scripts/louca2014_montecarlo.jl (data/louca2014/reproduction.toml).
using Chytridiomycosis
using Test
using Random: Xoshiro
using LinearAlgebra: eigvals
using TOML

@testset "Louca 2014 DOI-distribution model" begin
    C = Chytridiomycosis
    q = parameter_set(:louca2014)
    p3a = LoucaParameters(:louca2014; set = "fig3a")
    Y = q.Y

    @testset "compendium (Table 1, Table S1, ESM §4–5)" begin
        recs = parameter_records(:louca2014)
        @test Set(r.symbol for r in recs if !r.verified) == Set([:Y])   # days per year assumed
        @test time_unit(LoucaParameters) === :day
        @test (q.FA, q.FT, q.TA, q.TT, q.STmax, q.Tm, q.Tbr) == (1e4, 1e4, 1e3, 1e4, 9999, 1e6, 1e6)
        @test (q.μSA, q.ηA, q.r, q.τm, q.τbr, q.σA, q.σT, q.KA, q.KT, q.fA) ==
              (0.01, 2, 100, 365, 60, 0.5, 0.05, 100, 200, 0.2)
        @test (q.m, q.years, q.exposure_time, q.exposure_doi) == (100, 10, 1.5, 1000)
        p2 = LoucaParameters(:louca2014; set = "fig2cd")
        @test (p2.μZ, p2.μSA, p2.βA, p2.τbr, p2.TT, p2.STmax) == (1, 0.2, 1e-3, Y / 2, 500, 1e4)
        @test tadpoles_tolerate(p3a) && !tadpoles_tolerate(LoucaParameters(:louca2014; set = "fig3b"))
        # the lethal DOI is a swab load, comparable with Vredenburg's 10^4 ZE rule, not with
        # the Briggs model's sporangial ceiling
        @test parameter_load(:louca2014, :FA) == parameter_load(:kinney2011, :crash_load)
        @test_throws ArgumentError parameter_load(:louca2014, :FA) < parameter_load(:briggs2010, :Smax)
        rg = louca_ranges(; tadpoles = :tolerate, reservoir = false)
        @test !haskey(rg, :Zo) && !haskey(rg, :TT) && length(rg) == 21
        rg = louca_ranges(; tadpoles = :suffer, reservoir = true)
        @test length(rg) == 23 && rg.Zo.log && !rg.λA.log && rg.Tbr.log
        @test (rg.τm.lo, rg.τm.hi, rg.τbr.hi, rg.φ.hi) == (10, 4Y, Y / 2, Y)
        @test_throws ArgumentError LoucaParameters(:louca2014; set = "fig3a", TT = 500.0) |> LoucaDOIProblem
    end

    @testset "functional responses (ESM §1.4)" begin
        # breeding: wrapped Cauchy that integrates to r and has FWHM τbr
        x = louca_breeding_shape(p3a)
        tt = range(-Y / 2, Y / 2; length = 200_001)
        rb = [louca_breeding_rate(p3a, t) for t in tt]
        @test sum(rb[1:end-1]) * step(tt) ≈ p3a.r rtol = 1e-6
        @test louca_breeding_rate(p3a, p3a.τbr / 2) ≈ (maximum(rb) + minimum(rb)) / 2 rtol = 1e-8
        @test argmax(rb) == 100_001                                     # peak at t = 0
        @test louca_breeding_rate(C.with(p3a; τbr = Y / 2), 17.0) == p3a.r / Y
        @test isinf(louca_breeding_shape(C.with(p3a; τbr = Y / 2)))
        # transmission oscillation
        ps = C.with(p3a; b = 0.5, φ = 100.0)
        @test louca_transmission(ps, :adult, 100.0) ≈ 1.5 * ps.βA
        @test louca_transmission(ps, :tadpole, 100.0 + Y / 2) ≈ 0.5 * ps.βT
        # σ_br, σ_m halve at their tolerances
        @test louca_fecundity_factor(C.with(p3a; Tbr = 30.0), 30.0) == 0.5
        @test louca_metamorphosis_success(C.with(p3a; Tm = 30.0), 30.0) == 0.5
        # g_s(T_s) halves yearly survival; death at F_s; tolerant tadpoles have none
        @test exp(-Y * louca_disease_mortality(p3a, :adult, p3a.TA)) ≈ 0.5
        @test isinf(louca_disease_mortality(p3a, :adult, p3a.FA))
        @test louca_disease_mortality(p3a, :tadpole, 5000.0) == 0
        @test exp(-Y * louca_base_mortality(p3a, :adult)) ≈ p3a.σA
        # I_u(0) = I_d(0) = 0, I′(0) = λ, I(K_S) = 0 and the sign of I around K_S
        for (λ, μS) in ((0.3, 0.1), (-0.2, 0.1), (0.0, 0.1), (0.05, 0.5))
            pp = C.with(p3a; λA = λ, μSA = μS)
            @test louca_load_rates(pp, :adult, 0) == (0.0, 0.0)
            h = 1e-7
            @test louca_load_speed(pp, :adult, h) / h ≈ λ atol = 1e-6
            K = louca_within_host_capacity(pp, :adult)
            if λ > 0
                @test louca_load_speed(pp, :adult, K) ≈ 0 atol = 1e-10
                @test louca_load_speed(pp, :adult, K / 2) > 0 > louca_load_speed(pp, :adult, 2K)
                @test louca_within_host_loss(λ, K) ≈ μS
            else
                @test K == 0 && louca_load_speed(pp, :adult, 10.0) < 0
            end
        end
    end

    @testset "disease-free demography (ESM eqs 8–9)" begin
        p = LoucaParameters(:louca2014; set = "fig2ab")          # constant breeding
        c = louca_competition(p)
        @test c.cA > 0 && c.cT > 0
        # no exposure within the horizon: the population stays at (K_A, K_T)
        for method in (:ode, :semi_implicit)
            s = solve(LoucaDOIProblem(p; years = 3, exposure_time = 1e9, saveat = 365.0); method)
            @test s.NA[end] ≈ p.KA rtol = 1e-6
            @test s.NT[end] ≈ p.KT rtol = 1e-6
            @test s.outcome === :persistence && all(iszero, s.SA)
        end
        @test_throws ArgumentError LoucaDOIProblem(C.with(p; τm = 4Y, σA = 0.1, KT = 10.0))
    end

    @testset "Q_o and early growth (§3a; ESM §2.3, eqs 18–23)" begin
        R = louca_turnover(p3a)
        @test R.RA ≈ p3a.KT / (p3a.KA * p3a.τm) && R.RT ≈ p3a.r / Y * p3a.KA / p3a.KT
        Qo = reproduction_quotient(p3a)
        @test Qo ≈ (p3a.ηA * p3a.βA * p3a.KA / R.RA + p3a.ηT * p3a.βT * p3a.KT / R.RT) / p3a.μZ
        @test reproduction_quotient(LoucaParameters(:louca2014; set = "fig2ab")) ≈ Qo  # same Q_o in Figs 2 and 3
        @test isinf(reproduction_quotient(C.with(p3a; λT = 1.0)))
        # the keyword form generalises to several host types (ESM eq. 19)
        @test reproduction_quotient(; μZ = 0.1, η = [2, 2, 3], β = [1e-4, 1e-4, 1e-5],
                                    N = [100, 200, 50], R = [0.1, 0.2, 0.05], λ = [0, 0, 0]) ≈
              (2e-2 / 0.1 + 4e-2 / 0.2 + 1.5e-3 / 0.05) / 0.1
        # sign of the spectral bound of J̄_o is the sign of Q_o − 1, over sampled parameter sets
        rng = Xoshiro(37)
        n_ok = n_fin = 0
        for _ in 1:300
            s = louca_sample(rng; tadpoles = :suffer, reservoir = false)
            Q = reproduction_quotient(s.p)
            isfinite(Q) || continue
            n_fin += 1
            n_ok += sign(louca_growth_rate(s.p)) == sign(Q - 1)
        end
        @test n_fin >= 50 && n_ok == n_fin    # Q_o is finite in about a quarter of draws
        # ρ_o solves (ρ + μ_Z)ρ = Σ η β N, and eq. (23) approaches the eigenvalue when λ ≈ R̄
        ρo = louca_growth_rate(p3a; method = :leading)
        @test (ρo + p3a.μZ) * ρo ≈ p3a.ηA * p3a.βA * p3a.KA + p3a.ηT * p3a.βT * p3a.KT
        pc = C.with(p3a; λA = R.RA - 1e-4, λT = R.RT - 1e-4)
        @test louca_growth_rate(pc; method = :asymptotic) ≈ louca_growth_rate(pc) rtol = 1e-4
        @test louca_growth_rate(p3a; method = :asymptotic) ≈ louca_growth_rate(p3a) rtol = 0.1
        J = louca_jacobian(p3a)
        @test size(J) == (3, 3) && J[3, 1] == p3a.ηA * p3a.KA && J[1, 3] == p3a.βA
    end

    @testset "parameter estimates (ESM §4, §5.1)" begin
        @test round(louca_within_host_loss(q.fit_λ, q.fit_KS); digits = 2) == 0.25
        @test round(louca_tolerance_from_slope(q.fit_gprime; Y); digits = 2) == 1.73
        σ = [louca_survival_curve(T; xo = q.fit_xo, KS = q.fit_KS, λ = q.fit_λ,
                                  g′0 = q.fit_gprime, μo = q.fit_μo) for T in (0, 30, 60)]
        @test σ[1] == 1 && 1 > σ[2] > σ[3] > 0.3                         # Fig. S4b: ≈0.45 at 60 d
        lo = louca_transmission_from_growth(q.outbreak_ρ; μZ = q.est_μZ_min, η = q.est_η_max,
                                            N = q.outbreak_N_max)
        hi = louca_transmission_from_growth(q.outbreak_ρ; μZ = q.est_μZ_max, η = q.est_η_min,
                                            N = q.outbreak_N_min)
        @test round.((lo, hi); sigdigits = 1) == (7e-6, 2e-3)
        R = q.cruciger_recruits / q.cruciger_N / Y
        lo = louca_endemic_transmission(; μZ = q.est_μZ_min, R, η = q.est_η_max, N = q.cruciger_N)
        hi = louca_endemic_transmission(; μZ = q.est_μZ_max, R, η = q.est_η_min, N = q.cruciger_N)
        @test round.((lo, hi); sigdigits = 1) == (2e-7, 3e-5)
        rec = only(filter(r -> r.symbol === :johnson_survival, parameter_records(:louca2014)))
        @test round(1 / in_time_unit(rec, :day); digits = 2) == q.est_μZ_min   # 7 weeks
        @test q.reservoir_η * q.reservoir_hosts * q.reservoir_doi / q.reservoir_μZ == 1e12
    end

    @testset "binned model: numerics and DOI distributions (ESM §5)" begin
        p = LoucaParameters(:louca2014; set = "fig2cd")
        prob = LoucaDOIProblem(p; years = 4, save_distributions = true, saveat = 5.0)
        so = solve(prob; method = :ode)
        ss = solve(prob; dt = 0.02)
        @test so.outcome === ss.outcome === :persistence
        @test ss.NA[end] ≈ so.NA[end] rtol = 0.01
        @test ss.SA[end] ≈ so.SA[end] rtol = 0.02
        @test so.doi.adult == [(k - 1) * p.FA / p.m for k in 1:p.m]
        @test sum(so.PA) ≈ so.NA[end] && sum(so.doi.adult .* so.PA) / sum(so.PA) ≈ so.SA[end]
        @test size(so.PA_t) == (p.m, length(so.t))
        @test all(>=(0), ss.PA) && all(>=(0), ss.PT)
        # endemic adult DOI distribution: overdispersed, with mass near the lethal DOI
        μ = so.SA[end]
        v = sum(so.doi.adult .^ 2 .* so.PA) / sum(so.PA) - μ^2
        @test v > μ
        # the exposure adds one adult with DOI 10^3 (bin 11 of width 100)
        s1 = solve(LoucaDOIProblem(p; years = 1.5 + 0.5 / Y, saveat = 0.25); dt = 0.01)
        i = findfirst(>(1.5Y), s1.t)
        @test s1.NA[i] ≈ p.KA + 1 rtol = 0.01
        @test s1.SA[i] ≈ q.exposure_doi / (p.KA + 1) rtol = 0.2
        @test to_ode_problem(prob) isa C.ODEProblem
    end

    @testset "Figs 2–3 (Table S1)" begin
        amp(s, f) = (i = findall(>(3Y), s.t); x = getfield(s, f)[i]; maximum(x) - minimum(x))
        avg(s, f) = (i = findall(>(3Y), s.t); x = getfield(s, f)[i]; sum(x) / length(x))
        a = solve(LoucaDOIProblem(LoucaParameters(:louca2014; set = "fig2ab")); method = :ode)
        c = solve(LoucaDOIProblem(LoucaParameters(:louca2014; set = "fig2cd")); method = :ode)
        # Fig. 2: both endemic; tolerant long-lived tadpoles give strong outbreaks and host
        # fluctuations, suffering tadpoles damp them; adults fewer when tadpoles tolerate
        @test a.outcome === c.outcome === :persistence
        @test amp(a, :NT) > 5amp(c, :NT) && amp(a, :SA) > 2amp(c, :SA)
        @test maximum(a.SA) > 0.7 * q.FA                               # Fig. 2b: peaks near 10^4
        @test avg(a, :NA) < avg(c, :NA)
        # Fig. 3: tolerant tadpoles rescue the population; suffering ones go extinct before 2 yr
        s3a = solve(LoucaDOIProblem(p3a); method = :ode)
        s3b = solve(LoucaDOIProblem(LoucaParameters(:louca2014; set = "fig3b")); method = :ode)
        @test s3a.outcome === :persistence && s3a.NA[end] < p3a.KA / 5
        @test s3b.outcome === :extinction && 1.5Y < s3b.t_end < 2Y
        # stochastic model (ESM §5.3), seeded
        st = [solve(LoucaDOIProblem(LoucaParameters(:louca2014; set = "fig3b"); stochastic = true);
                    dt = 0.05, rng = Xoshiro(i)) for i in 1:3]
        @test count(s -> s.outcome === :extinction, st) >= 2
        sc = solve(LoucaDOIProblem(LoucaParameters(:louca2014; set = "fig2cd"); stochastic = true);
                   dt = 0.05, rng = Xoshiro(2))
        @test sc.outcome === :persistence && amp(sc, :NT) > 5amp(c, :NT)   # quasi-cycles
        sco = solve(LoucaDOIProblem(LoucaParameters(:louca2014; set = "fig2cd"); years = 3,
                                    stochastic = true); method = :ode, rng = Xoshiro(2))
        @test sco.outcome === :persistence && minimum(sco.NA) < 0.9 * maximum(sco.NA)
    end

    @testset "Monte Carlo machinery (ESM §5.2)" begin
        rng = Xoshiro(5)
        for (tad, res) in ((:suffer, false), (:tolerate, true))
            s = louca_sample(rng; tadpoles = tad, reservoir = res)
            c = louca_competition(s.p)
            @test c.cA > 0 && c.cT > 0 && all(x -> 0 <= x <= 1, s.u)
            @test tadpoles_tolerate(s.p) == (tad === :tolerate)
            @test (s.p.Zo > 0) == res
        end
        # risk effect of a step response in u is −1.5 (slope of 1{u < 1/2} on u ~ U(0, 1))
        fake = [(u = (μZ = x, b = y), Qo = 10.0^(8x - 4), outcome = x < 0.5 ? :extinction : :clearance)
                for (x, y) in zip(range(0, 1; length = 4001), rand(Xoshiro(1), 4001))]
        re = louca_risk_effects(fake)
        @test re.μZ ≈ -1.5 rtol = 0.01
        @test re.Qo ≈ -1.5 rtol = 0.01
        @test abs(re.b) < 0.1
        pr = louca_outcome_probabilities(fake; edges = [-4.0, 0.0, 4.0])
        @test pr.extinction == [1.0, 0.0] && pr.n == [2000, 2000]
        mc = louca_monte_carlo(3; rng = Xoshiro(9), dt = 0.1)
        @test length(mc) == 3 && all(r -> r.outcome in C.LOUCA_OUTCOMES, mc)
    end

    @testset "stored Monte Carlo (scripts/louca2014_montecarlo.jl)" begin
        rep = load_louca2014()
        sc = ("suffer_reservoir", "tolerate_reservoir", "suffer_noreservoir", "tolerate_noreservoir")
        for k in sc
            d = rep[k]
            @test d["extinction"] + d["clearance"] + d["persistence"] ≈ 1
            @test d["dt_check_agree"] >= 0.97 * d["dt_check_runs"]       # dt = 0.05 vs 0.01
        end
        # clearance is impossible with a reservoir (§3b)
        @test rep["suffer_reservoir"]["clearance"] == 0 == rep["tolerate_reservoir"]["clearance"]
        @test rep["reservoir_onset"] == "exposure"
        ext(k) = rep[k]["extinction"]
        re(k, s) = rep[k]["risk_effects"][s]
        cols = ("suffer_reservoir", "tolerate_reservoir", "suffer_noreservoir", "tolerate_noreservoir")
        # §3d: more extinctions when tadpoles suffer (25% vs 21%) — the ordering reproduces ...
        @test ext("suffer_noreservoir") > ext("tolerate_noreservoir")
        # ... but not the levels (0.2145, 0.1835); not dt, bins, sampling scale or reservoir onset
        @test_broken abs(ext("suffer_noreservoir") - 0.25) < 0.02 && abs(ext("tolerate_noreservoir") - 0.21) < 0.02
        # Table 2: β_A and μ_Z are the top parameters without a reservoir; Q_o beats them
        for k in cols[3:4]
            @test abs(re(k, "Qo")) > max(abs(re(k, "βA")), abs(re(k, "μZ"))) > 0.4
            @test re(k, "βA") > 0 > re(k, "μZ") && re(k, "λA") > 0
        end
        # Table 2 discrepancies (wiki card): T_A (−0.38 … −0.17 published) has ~no effect here ...
        @test_broken all(re(k, "TA") < -0.15 for k in cols)
        # ... nor do K_T, K_A without a reservoir (−0.27/−0.35, −0.22/−0.12 published) ...
        @test_broken re(cols[3], "KT") < -0.15 && re(cols[3], "KA") < -0.15
        # ... Z_o is stronger (1.15/1.12 vs 0.66/0.78) and beats Q_o with a reservoir ...
        @test_broken re(cols[1], "Qo") >= re(cols[1], "Zo")
        # ... and Q_o's risk effect depends on an unstated range (0.90, 0.66, 1.13, 1.02 here)
        @test_broken all(abs(re(k, "Qo") - v) <= 0.1 for (k, v) in zip(cols, (0.66, 0.83, 0.81, 0.71)))
        # β_A is about 0.17 stronger than published (0.32 … 0.39)
        @test_broken all(abs(re(k, "βA") - v) <= 0.1 for (k, v) in zip(cols, (0.32, 0.34, 0.38, 0.39)))
    end

    @testset "reservoir onset" begin
        p = C.with(LoucaParameters(:louca2014; set = "fig2cd"); Zo = 1e3)
        a = LoucaDOIProblem(p; years = 2, saveat = 30.0)
        b = LoucaDOIProblem(p; years = 2, saveat = 30.0, reservoir_onset = :start)
        @test a.u0[end] == 0 && b.u0[end] == p.Zo
        sa = solve(a; dt = 0.05)
        @test all(iszero, sa.SA[sa.t .< 1.5Y])                 # healthy until exposure
        sb = solve(b; dt = 0.05)
        @test sb.SA[findfirst(>(30), sb.t)] > 0                # infected during the burn-in
        @test_throws ArgumentError LoucaDOIProblem(p; reservoir_onset = :never)
    end
end
