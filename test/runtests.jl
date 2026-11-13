using Chytridiomycosis
using Test
using Distributions: LogNormal, NegativeBinomial, mean, pdf, quantile
using Random: Xoshiro

@testset "Chytridiomycosis" begin
    @testset "compendium" begin
        @test :schmidt2017 in compendium_sources()
        recs = parameter_records(:schmidt2017)
        @test all(r -> r.provenance in ("M", "F", "E", "A", "D"), recs)
        @test all(r -> r.verified, recs)
        ps = parameter_set(:schmidt2017)
        @test ps.b == 0.39 && ps.d_I == 52.0 && ps.e == 364.0 && ps.K == 25.0
        @test_throws ArgumentError parameter_records(:no_such_paper)
    end

    # Oracles: Schmidt et al. 2017, R Soc Open Sci 4:160801, Results section (PMC5383822).
    @testset "Schmidt 2017 Bsal S-L-I" begin
        p = SchmidtParameters(:schmidt2017)
        @test round(invasion_threshold(p), digits = 4) == 1.0086         # K_threshold
        @test round(100 * removal_fraction(p)) == 96                      # p_remove ≈ 96%
        eq = endemic_equilibrium(p)
        # Paper reports S_eq = 1.0027; eq. (2.3) with Table 1 values gives 1.0030 (not explained
        # by rounding of e, d_I or β). Open discrepancy: lab/wiki/models/bsal-schmidt2017-sli.md
        @test_broken round(eq.S, digits = 4) == 1.0027
        @test round(eq.S, digits = 4) == 1.0030                           # regression
        du = zeros(3)
        Chytridiomycosis.schmidt_rhs!(du, [eq.S, eq.L, eq.I], p, 0.0)
        @test maximum(abs, du) < 1e-10                                    # really an equilibrium
        @test endemic_equilibrium(Chytridiomycosis.with(p; K = 1.0)) === nothing
        KH = hopf_threshold(p)
        # Paper reports K_H ≈ 9.5; the Jacobian at the (verified) equilibrium gives ≈ 10.21.
        @test_broken isapprox(KH, 9.5; atol = 0.1)
        @test isapprox(KH, 10.21; atol = 0.01)                            # regression
        @test Chytridiomycosis.equilibrium_stability(Chytridiomycosis.with(p; K = 10.0)) < 0
        @test Chytridiomycosis.equilibrium_stability(Chytridiomycosis.with(p; K = 10.5)) > 0
        @test round(wave_speed(p), digits = 2) == 11.15                   # c = 11.15 km/yr
        crash = population_crash(p)
        @test crash.Smin < 0.05 * p.K                                     # "almost zero"
        @test 1 / 12 < crash.tmin < 4 / 12                                # "two to three months"
        # Eq. (2.5): the paper's ≈9% (μ = 1) and ≈0.1% (μ = 0.1) of 258 750 km² after 8 yr
        # equal A(t)/100 — a probable km²/ha slip, recorded in the wiki.
        habitat = 258_750.0
        f1 = scattered_colony_area(p, 8.0; μ = 1.0) / 100 / habitat
        f01 = scattered_colony_area(p, 8.0; μ = 0.1) / 100 / habitat
        @test round(100 * f1) == 9
        @test round(100 * f01, digits = 1) == 0.1
        sol = solve(SchmidtSLIProblem(p; tspan = (0.0, 0.5)))
        @test Chytridiomycosis.SciMLBase.successful_retcode(sol)
    end

    # Briggs, Knapp & Vredenburg 2010 PNAS, Fig. 2C parameters, β = νγ. Equations and values
    # checked against the PDF and SI (Eqs. S1–S2, Table S4).
    @testset "Briggs 2010 Bd load model" begin
        p = BriggsParameters(:briggs2010)
        @test transmission_rate(p) ≈ 5e-4
        Nstar = host_density_threshold(p)
        @test round(Nstar, digits = 2) == 23.15        # Fig. 2C legend: "a density of ∼20 frogs"
        @test load_growth_rate(p, Nstar) ≈ 0 atol = 1e-12
        @test load_growth_rate(p, 20) < 0 < load_growth_rate(p, 30)
        @test host_density_threshold(BriggsParameters(; V = 1.0, γ = 0.01, ν = 0.05, f = 0.9,
            μ = 1.0, η = 17.5, σ = 0.2, Smax = 1e4)) == 0   # self-reinfection alone suffices
        # Below N*: all frogs clear; above N*: all reach Smax, sooner at higher density.
        run(N) = solve(BriggsLoadProblem(p; N = N, S0 = fill(10.0, N), tspan = (0.0, 1095.0));
                       abstol = 1e-8, reltol = 1e-8)
        low = run(10)
        @test all(low.alive) && maximum(low.sol.u[end][1:10]) < 1e-6
        t30, t60, t120 = (maximum(run(N).death_time) for N in (30, 60, 120))
        @test all(!isnan, (t30, t60, t120))
        @test t30 > t60 > t120
        r = run(30)
        @test count(!, r.alive) == 30 && all(iszero, r.sol.u[end][1:30])
        # Heterogeneous loads: the heaviest-infected frog dies first.
        het = solve(BriggsLoadProblem(p; N = 40, S0 = collect(2.0:41.0), tspan = (0.0, 1095.0)))
        @test all(!, het.alive)
        @test argmin(het.death_time) == 40
    end

    # Briggs 2010 stochastic seasonal model (SI Methods; Fig. 3 legend values). The paper
    # gives persistence maps, not numbers, so these are the qualitative statements of its
    # Results; Δt and the Poisson cut-off are our choices.
    @testset "Briggs 2010 seasonal stochastic model" begin
        @test parameter_set(:briggs2010).ν == 0.05
        @test parameter_set(:briggs2010; set = "fig3b").ν == 0.1
        @test_throws ArgumentError parameter_set(:briggs2010; set = "no_such_set")
        pb = BriggsSeasonalParameters(:briggs2010; set = "fig3b")
        pd = BriggsSeasonalParameters(:briggs2010; set = "fig3d")
        pf = BriggsSeasonalParameters(:briggs2010; set = "fig3f")
        @test pb.tadpole === nothing && pb.εZ == 0 && pd.εZ == 1000
        @test pf.tadpole.γ ≈ 100 * pf.adult.γ && !pf.tadpole.lethal
        # With small steps and no Poisson draws, one summer matches the ODE model.
        p2 = BriggsParameters(:briggs2010)
        det = BriggsSeasonalParameters(; V = p2.V, μ = p2.μ, η = p2.η, σ = p2.σ,
            adult = BriggsStage(; γ = p2.γ, ν = p2.ν, f = p2.f, Smax = p2.Smax, θ = 1.0),
            R = 0.0, K = 30, θZ = 1.0, summer = 50.0)
        r = solve(BriggsSeasonalProblem(det; years = 1, Δt = 0.01, poisson_below = 0.0);
            rng = Xoshiro(1))
        ode = solve(BriggsLoadProblem(p2; N = 30, S0 = [10.0; zeros(29)], tspan = (0.0, 50.0));
            abstol = 1e-10, reltol = 1e-10)
        @test isapprox(r.mean_adult_load[1], sum(ode.sol.u[end][1:30]) / 30; rtol = 0.01)
        @test isapprox(r.Z[1], ode.sol.u[end][end]; rtol = 0.01)
        # Unstructured host: persistence only in a narrow band (Fig. 3A).
        withγ(γ) = BriggsSeasonalParameters(; V = pb.V, μ = pb.μ, η = pb.η, σ = pb.σ,
            adult = BriggsStage(; γ, ν = pb.adult.ν, f = pb.adult.f, Smax = pb.adult.Smax,
                θ = pb.adult.θ), R = pb.R, K = pb.K, θZ = pb.θZ)
        pers(p) = persistence_probability(BriggsSeasonalProblem(p); nruns = 40, rng = Xoshiro(7))
        @test 0.3 < pers(pb) < 0.95                    # the starred point of Fig. 3A
        low = solve(BriggsSeasonalProblem(withγ(1e-9)); rng = Xoshiro(1))
        @test !low.persisted && low.hosts_final > 0 && low.bd_final == 0   # Bd dies out
        high = solve(BriggsSeasonalProblem(withγ(1e-2)); rng = Xoshiro(1))
        @test !high.persisted && sum(high.deaths) >= 100                    # frogs driven out
        # External source or long-lived tadpoles: persistence with no adult deaths (Fig. 3D, F).
        rd = solve(BriggsSeasonalProblem(pd); rng = Xoshiro(1))
        rf = solve(BriggsSeasonalProblem(pf); rng = Xoshiro(1))
        @test rd.persisted && sum(rd.deaths) == 0 && all(==(100), rd.adults)
        @test rf.persisted && sum(rf.deaths) == 0
        @test all(<(pf.adult.Smax / 10), rf.mean_adult_load)   # "low to intermediate" loads
        @test pers(pd) == 1.0
    end

    # Mitchell et al. 2008 Proc R Soc B, Table 1 and Fig. 2 (μZ = 13 yr⁻¹, γ = 0): the three
    # outcomes 35 years after Bd is introduced, at the points marked b, c and d in Fig. 2a.
    @testset "Mitchell 2008 pond model" begin
        run(set) = solve(MitchellPondProblem(MitchellParameters(:mitchell2008; set)))
        base, b, c, d = run(""), run("fig2b"), run("fig2c"), run("fig2d")
        @test 9_500 < base.A_disease_free < 11_000        # Fig. 2b–d: about 10 × 10³ adults
        @test b.toads_persist && b.bd_persist             # b: both toads and Bd
        @test c.toads_persist && !c.bd_persist            # c: toads only
        @test c.A[end] > 0.95 * c.A_disease_free          #    ... and the toads recover
        @test minimum(c.A) < 0.2 * c.A_disease_free       #    ... after a deep crash
        @test !d.toads_persist && !d.bd_persist           # d: neither toads nor Bd
        @test base.toads_persist && !base.bd_persist      # Table 1 values: "toads only" region
        # Without Bd the yearly cycle repeats: adults are the same one year apart.
        i0 = findfirst(>=(0), base.t)
        # Two events a year before Bd arrives (birth, metamorphosis), so a year is 2 entries.
        @test base.t[i0 - 2] - base.t[i0 - 4] ≈ 1
        @test isapprox(base.A[i0 - 2], base.A[i0 - 4]; rtol = 1e-6)
    end

    # New model: birth–death–immigration load process (lab/wiki/models/new-load-bdi.md).
    @testset "Load BDI (new)" begin
        m = LoadBDI(0.04, 0.14, 0.2)
        d = load_distribution(m)
        @test d isa NegativeBinomial && d.r ≈ aggregation(m) ≈ 0.04 / 0.14
        # Stationarity: detailed balance of the BDI master equation, (ι + λn) π_n = σ(n+1) π_{n+1}.
        @test all(isapprox((m.ι + m.λ * n) * pdf(d, n), m.σ * (n + 1) * pdf(d, n + 1);
            rtol = 1e-10) for n in 0:200)
        @test mean(d) ≈ mean_load(m) ≈ m.ι / (m.σ - m.λ)
        @test prevalence(m) ≈ 1 - pdf(d, 0)
        @test infection_loss_rate(m) ≈ m.σ * pdf(d, 1) / (1 - pdf(d, 0))
        # Transient: mean follows the Briggs equation; shape k is constant; t → ∞ is stationary.
        briggs_mean(t) = m.ι / (m.λ - m.σ) * (exp((m.λ - m.σ) * t) - 1)
        @test all(mean(load_distribution(m, t)) ≈ briggs_mean(t) for t in (1.0, 10.0, 100.0))
        @test load_distribution(m, 5.0).r == d.r
        @test mean_load(m, 1e4) ≈ mean_load(m) rtol = 1e-8
        grow = LoadBDI(0.04, 0.3, 0.2)                 # epizootic: λ > σ
        @test_throws ArgumentError load_distribution(grow)
        @test mean_load(grow, 60.0) > 10 * mean_load(grow, 30.0) &&
              load_distribution(grow, 60.0).r ≈ 0.04 / 0.3
        @test mean_load(LoadBDI(0.04, 0.2, 0.2), 10.0) ≈ 0.4          # critical case: ι t
        # Swabbing thins the count and keeps the shape.
        s = swab(d, 0.1)
        @test s.r == d.r && mean(s) ≈ 0.1 * mean(d)
        # Enzootic state with an external source (Briggs Fig. 3D values).
        q = parameter_set(:briggs2010; set = "fig3d")
        p = BriggsParameters(; V = q.V, γ = q.γ, ν = q.ν, f = q.f, μ = q.μ, η = q.η, σ = q.σ,
            Smax = q.Smax)
        Z = enzootic_pool(p; N = 100, ε = q.εZ)
        e = enzootic_load_distribution(p; N = 100, ε = q.εZ)
        # Z* balances the pool: input + release by hosts at mean load = loss.
        @test q.εZ + p.η * (1 - p.f) * 100 * mean(e) ≈ (p.μ + p.γ * 100 / p.V) * Z
        Nstar = host_density_threshold(p)
        @test_throws ArgumentError enzootic_pool(p; N = Nstar + 1, ε = q.εZ)
        @test mean(enzootic_load_distribution(p; N = 0.99 * Nstar, ε = q.εZ)) > 20 * mean(e)
        # It predicts the Briggs seasonal simulation it was not fitted to (end-of-summer
        # prevalence and mean load over 10 years, 100 frogs).
        sim = solve(BriggsSeasonalProblem(BriggsSeasonalParameters(:briggs2010; set = "fig3d"));
            rng = Xoshiro(1))
        @test isapprox(mean(sim.infected_adults) / 100, 1 - pdf(e, 0); atol = 0.1)
        @test isapprox(mean(sim.mean_adult_load), mean(e); rtol = 0.5)
    end

    # New model: load tipping, the Briggs equation plus saturating clearance
    # (lab/wiki/models/new-load-tipping.md). c and h are illustrative, not published.
    @testset "Load tipping (new)" begin
        rate(m, S) = m.ι + net_growth_rate(m) * S - m.c * S / (m.h + S)
        slope(m, S) = net_growth_rate(m) - m.c * m.h / (m.h + S)^2
        m = LoadTipping(; ι = 0.5, λ = 0.25, σ = 0.2, c = 20.0, h = 100.0)
        e = load_equilibria(m)
        @test tipping_regime(m) === :bistable
        @test abs(rate(m, e.low)) < 1e-10 && abs(rate(m, e.threshold)) < 1e-10
        @test slope(m, e.low) < 0 < slope(m, e.threshold)       # stable low, unstable threshold
        # The two equilibria meet at ι_crit = (√c − √(R h))² and then vanish.
        ιc = critical_immigration(m)
        @test ιc ≈ (sqrt(20.0) - sqrt(0.05 * 100))^2
        at = load_equilibria(LoadTipping(; ι = ιc * (1 - 1e-9), λ = 0.25, σ = 0.2, c = 20.0, h = 100.0))
        @test isapprox(at.low, at.threshold; rtol = 1e-3) && isapprox(at.low, sqrt(ιc * 100 / 0.05); rtol = 1e-3)
        over = LoadTipping(; ι = 1.01ιc, λ = 0.25, σ = 0.2, c = 20.0, h = 100.0)
        @test tipping_regime(over) === :runaway && isnan(load_equilibria(over).low)
        @test tipping_regime(LoadTipping(; ι = 0.1, λ = 0.5, σ = 0.2, c = 20.0, h = 100.0)) === :runaway  # R h ≥ c
        # No clearance and R < 0: the Briggs / BDI mean load, and no threshold.
        lin = LoadTipping(; ι = 0.04, λ = 0.14, σ = 0.2, c = 0.0, h = 1.0)
        @test tipping_regime(lin) === :stable && critical_immigration(lin) == Inf
        @test load_equilibria(lin).low ≈ mean_load(LoadBDI(0.04, 0.14, 0.2))
        @test load_equilibria(lin).threshold == Inf

        # A population on the Briggs Fig. 2 values. The pool feedback κ(N) cancels the
        # within-host decline exactly at the Briggs threshold N*.
        p = BriggsParameters(:briggs2010)
        c, h, ε = 100.0, 500.0, 1000.0
        Nstar = host_density_threshold(p)
        @test abs(net_growth_rate(LoadTipping(p; c, h, N = Nstar))) < 1e-12
        @test tipping_regime(LoadTipping(p; c, h, N = 20, ε)) === :stable
        @test tipping_regime(LoadTipping(p; c, h, N = 50, ε)) === :bistable
        @test tipping_regime(LoadTipping(p; c, h, N = 70, ε)) === :runaway
        # Fold: with no outside source it is where R = c/h; a source brings it down.
        N0 = critical_host_number(p; c, h)
        @test pool_feedback(p, N0) + p.η * p.ν * p.f - p.σ ≈ c / h rtol = 1e-8
        Nε = critical_host_number(p; c, h, ε)
        @test Nstar < Nε < N0
        @test critical_host_number(p; c = 1e4, h) == Inf        # clearance holds at any density
        @test_throws ArgumentError LoadTipping(p; c, h)
        @test_throws ArgumentError LoadTipping(p; c, h, N = 10, Z = 1.0)

        # The full N-host ODE has the same equilibria and the same threshold.
        N = 50
        eq = load_equilibria(LoadTipping(p; c, h, N, ε))
        sim(S0) = solve(LoadTippingProblem(p; c, h, N, ε, S0, tspan = (0.0, 2000.0));
            abstol = 1e-8, reltol = 1e-8)
        below, above = sim(fill(0.9 * eq.threshold, N)), sim(fill(1.1 * eq.threshold, N))
        @test all(below.alive) && all(isapprox.(below.sol.u[end][1:N], eq.low; rtol = 1e-5))
        @test !any(above.alive)
        # Tipping is collective: what matters is the mean load, not any one host's load.
        kick(k, S) = (S0 = fill(eq.low, N); S0[1:k] .= S; sim(S0))
        @test all(kick(1, 9000.0).alive)          # mean 184 < threshold 320: all recover
        @test !any(kick(2, 9000.0).alive)         # mean 364 > threshold: all 50 die
        @test all(kick(5, 2000.0).alive) && !any(kick(10, 2000.0).alive)
        r = kick(2, 9000.0)
        @test maximum(r.death_time) - minimum(r.death_time) < 1    # and they die together

        # Hosts that differ in clearance (c lognormal across hosts, median 100).
        cs(n, s) = quantile.(LogNormal(log(100.0), s), ((1:n) .- 0.5) ./ n)
        same = pooled_load_equilibria(p; c = fill(c, N), h, ε)       # no variation: as above
        @test same.regime === :bistable && same.low[1] ≈ eq.low && same.threshold[1] ≈ eq.threshold
        @test pooled_load_equilibria(p; c = fill(c, N), h).ι_low == 0     # no source: Bd-free
        @test pooled_load_equilibria(p; c = cs(20, 0.5), h, ε).regime === :stable
        @test pooled_load_equilibria(p; c = cs(80, 0.5), h, ε).regime === :runaway
        cv = cs(N, 0.5)
        het = pooled_load_equilibria(p; c = cv, h, ε)
        @test het.regime === :bistable && issorted(het.low; rev = true)   # weak defences, high load
        hsim(S0) = solve(LoadTippingProblem(p; c = cv, h, N, ε, S0, tspan = (0.0, 5000.0));
            abstol = 1e-8, reltol = 1e-8)
        stay = hsim(het.low)
        @test all(stay.alive) && all(isapprox.(stay.sol.u[end][1:N], het.low; rtol = 1e-5))
        @test all(hsim(0.9 .* het.threshold).alive)
        # Above the threshold the die-off is partial: the best-defended hosts survive, fewer
        # than the Briggs N*, and settle at the low state of the survivors alone.
        tip = hsim(1.1 .* het.threshold)
        ns = count(tip.alive)
        @test 0 < ns < Nstar && tip.alive == (1:N .> N - ns)
        after = pooled_load_equilibria(p; c = cv[tip.alive], h, ε)
        @test after.regime === :stable
        @test all(isapprox.(tip.sol.u[end][1:N][tip.alive], after.low; rtol = 1e-4))
        @test_throws ArgumentError LoadTippingProblem(p; c = cv[1:3], h, N, ε)

        # Recruitment and inherited clearance (stochastic between years).
        rec(K, h2; seed = 1) = solve(LoadTippingRecruitmentProblem(p; c_median = 100.0,
            sdlog = 0.5, h, K, ε, heritability = h2, years = 100); rng = Xoshiro(seed))
        calm = rec(50, 0.0)                         # ceiling below the fold: never tips
        @test all(==(0), calm.deaths) && all(==(50), calm.hosts)
        r0, r9 = rec(80, 0.0), rec(80, 0.9)         # ceiling above the fold
        @test r0.deaths[1] > 40 && 0 < r0.hosts[2] < Nstar          # first tip: partial, to below N*
        @test all(<=(80), r0.hosts) && r0.hosts == rec(80, 0.0).hosts   # bounded; reproducible
        @test count(>(0), r0.deaths) >= 8           # no inheritance: die-offs keep recurring
        @test count(>(0), r9.deaths) <= 3           # strong inheritance: few
        @test r9.median_c[30] > 1.2 * r9.median_c[1]    # the survivors' clearance is passed on
        @test_throws ArgumentError LoadTippingRecruitmentProblem(p; c_median = 100.0, sdlog = 0.5,
            h, K = 80, heritability = 1.5)

        # Stochastic single host: mean first-passage time against a direct linear solve.
        M = 40
        Q = zeros(M, M)
        for n in 0:(M - 1)
            b, d = m.ι + m.λ * n, m.σ * n + m.c * n / (m.h + n)
            Q[n + 1, n + 1] = -(b + d)
            n < M - 1 && (Q[n + 1, n + 2] = b)
            n > 0 && (Q[n + 1, n] = d)
        end
        T = Q \ fill(-1.0, M)
        @test mean_time_to_runaway(m; Smax = M) ≈ T[1] rtol = 1e-8
        @test mean_time_to_runaway(m; Smax = M, from = 10) ≈ T[11] rtol = 1e-8
        times = [mean_time_to_runaway(LoadTipping(; ι, λ = 0.25, σ = 0.2, c = 20.0, h = 100.0);
            Smax = 10_000) for ι in (0.5, 2.0, 4.0, 5.0)]
        @test issorted(times; rev = true) && times[1] > 1e15 && times[end] < 1e3
        @test mean_time_to_runaway(LoadTipping(; ι = 0.0, λ = 0.25, σ = 0.2, c = 20.0, h = 100.0);
            Smax = 100) == Inf
        @test_throws ArgumentError mean_time_to_runaway(LoadTipping(p; c, h, N = 50, ε); Smax = 100)
        # Low-state distribution: detailed balance; and the negative binomial when c = 0.
        d = low_state_distribution(m)
        @test all(isapprox((m.ι + m.λ * n) * pdf(d, n),
            (m.σ * (n + 1) + m.c * (n + 1) / (m.h + n + 1)) * pdf(d, n + 1); rtol = 1e-10) for n in 0:100)
        nb = load_distribution(LoadBDI(0.04, 0.14, 0.2))
        @test all(isapprox(pdf(low_state_distribution(lin; nmax = 2000), n), pdf(nb, n); rtol = 1e-8)
            for n in 0:50)
    end

    # New model: a "treated" class carrying a protective microbiome, on three host–pathogen
    # models (lab/wiki/models/new-microbiome-treatment.md). Treatment values are illustrative.
    @testset "Microbiome treatment (new)" begin
        tr = MicrobiomeTreatment(efficacy = 0.99, half_life = 1.0, direct = 50.0)
        @test reversion_rate(tr) ≈ log(2) && reversion_rate(MicrobiomeTreatment()) == 0
        @test microbiome_reproduction_number(tr; turnover = 0.39) ≈ 50 / (log(2) + 0.39)
        @test equilibrium_coverage(tr; turnover = 0.39) ≈ 1 - (log(2) + 0.39) / 50
        @test equilibrium_coverage(MicrobiomeTreatment(half_life = 1.0, direct = 0.5)) == 0
        @test equilibrium_coverage(MicrobiomeTreatment(half_life = 1.0, environmental = 50.0)) ≈
              equilibrium_coverage(MicrobiomeTreatment(half_life = 1.0, direct = 50.0))
        @test treatment_efficacy(0.6, 0.2, 1.0) ≈ 0.5
        @test_throws ArgumentError MicrobiomeTreatment(efficacy = 1.2)

        # Bsal in salamanders (Schmidt 2017).
        ps = SchmidtParameters(:schmidt2017)
        none = solve(ProbioticSLIProblem(ps, MicrobiomeTreatment(); tspan = (0.0, 1.0)))
        ref = solve(SchmidtSLIProblem(ps); abstol = 1e-10, reltol = 1e-8)
        i = findfirst(>=(0.1), none.t)             # compare during the crash, not at ~1e-10
        @test isapprox(none.hosts[i], sum(ref(none.t[i])); rtol = 1e-3)   # no treatment: Schmidt
        @test critical_coverage(ps, tr) ≈ (1 - 1 / basic_reproduction_number(ps)) / 0.99
        @test critical_coverage(ps, MicrobiomeTreatment(efficacy = 0.9)) > 1     # cannot work
        @test effective_reproduction_number(ps, tr, critical_coverage(ps, tr)) ≈ 1
        sli(t; kw...) = solve(ProbioticSLIProblem(ps, t; coverage = 1.0, kw...))
        @test sli(tr).mean_hosts > 0.99                                    # self-spreading, 99%
        @test sli(MicrobiomeTreatment(efficacy = 0.99, half_life = 1.0)).mean_hosts < 0.2
        @test sli(MicrobiomeTreatment(efficacy = 0.9, half_life = 1.0, direct = 50.0)).mean_hosts < 0.4

        # Surviving infection (tolerance) instead of blocking it: Bsal's reproduction number
        # rises, no coverage stops invasion, and survivors are carriers.
        tol = MicrobiomeTreatment(tolerance = 1.0)
        @test_throws ArgumentError MicrobiomeTreatment(tolerance = -0.1)
        @test effective_reproduction_number(ps, tol, 1.0) ≈
              basic_reproduction_number(ps) * (ps.b + ps.d_I) / ps.b
        @test effective_reproduction_number(ps, tol, 0.0) ≈ basic_reproduction_number(ps)
        @test critical_coverage(ps, tol) == Inf
        @test critical_coverage(ps, MicrobiomeTreatment(efficacy = 0.99, tolerance = 0.9)) > 1
        carriers = sli(tol)
        @test carriers.final_prevalence > 0.9 && carriers.final_hosts < 0.5
        @test sli(MicrobiomeTreatment(tolerance = 0.9)).mean_hosts < 0.05

        # Lower load: less shedding and less death together. With death proportional to load
        # the two nearly cancel in R; a threshold-like death response makes treatment harmful.
        lr = load_reducing_treatment(0.9; half_life = 2.0)
        @test lr.shedding_reduction == 0.9 && lr.tolerance ≈ 0.9 && lr.half_life == 2.0
        @test load_reducing_treatment(0.9; death_exponent = 2).tolerance ≈ 0.99
        @test load_reducing_treatment(0.9; death_exponent = 0).tolerance == 0
        @test_throws ArgumentError load_reducing_treatment(1.5)
        R0s = basic_reproduction_number(ps)
        Rlr(k) = effective_reproduction_number(ps, load_reducing_treatment(0.9; death_exponent = k), 1.0)
        @test Rlr(1) ≈ R0s * 0.1 * (ps.b + ps.d_I) / (ps.b + 0.1 * ps.d_I)
        @test Rlr(0) ≈ 0.1 * R0s && Rlr(0) < Rlr(1) < R0s < Rlr(2)
        @test effective_reproduction_number(ps, MicrobiomeTreatment(tolerance = 1.0,
            shedding_reduction = 1.0), 1.0) == 0
        # No recovery in the Schmidt model: treated hosts infected at the start stay infected,
        # so prevalence stays high even when they shed nothing.
        mute = sli(MicrobiomeTreatment(tolerance = 1.0, shedding_reduction = 1.0))
        @test 0.5 < mute.final_prevalence < carriers.final_prevalence
        @test mute.final_hosts > carriers.final_hosts

        # Trial summaries from the compendium (outcome counts, not model parameters).
        hr = parameter_set(:harris2009ismej)
        @test all(r -> r.verified, parameter_records(:harris2009ismej))
        @test !any(r -> r.verified, parameter_records(:becker2009aem))
        red = mortality_reduction(hr.deaths_treated, hr.n_treated, hr.deaths_control, hr.n_control;
            rng = Xoshiro(1))
        @test red.lower < red.median < red.upper <= 1 && red.median > 0.8
        @test_throws ArgumentError mortality_reduction(7, 6, 5, 6)

        # Bd in toads (Mitchell 2008).
        pb = MitchellParameters(:mitchell2008; set = "fig2b")
        m0, t0 = solve(MitchellPondProblem(pb)), solve(ProbioticPondProblem(pb, MicrobiomeTreatment()))
        @test isapprox(t0.A_final, m0.A[end]; rtol = 1e-4) && t0.bd_persist == m0.bd_persist
        pd = MitchellParameters(:mitchell2008; set = "fig2d")
        prot = MicrobiomeTreatment(efficacy = 0.9, half_life = 1.0)
        @test !solve(ProbioticPondProblem(pd, MicrobiomeTreatment())).toads_persist
        tad = solve(ProbioticPondProblem(pd, MicrobiomeTreatment(efficacy = 0.9); coverage = 0.8))
        @test tad.toads_persist && !tad.bd_persist && tad.A_mean > 0.9 * tad.A_disease_free
        @test maximum(tad.tadpoles_treated) ≈ 0.8
        ad = solve(ProbioticPondProblem(pd, prot; stage = :adults, coverage = 0.8))
        @test !ad.toads_persist && maximum(ad.A_treated) > 0     # treated adults, unprotected tadpoles
        @test solve(ProbioticPondProblem(pd, prot; stage = :adults, coverage = 0.8,
            maternal = 0.5)).toads_persist
        stop5 = solve(ProbioticPondProblem(pd, MicrobiomeTreatment(efficacy = 0.9); coverage = 0.8,
            treat_years = 0:4))
        @test !stop5.toads_persist                               # stopping lets Bd back
        surv = solve(ProbioticPondProblem(pd, MicrobiomeTreatment(tolerance = 0.7); coverage = 0.8))
        @test surv.toads_persist && surv.bd_persist              # toads saved, Bd kept
        @test 0.3 * surv.A_disease_free < surv.A_mean < tad.A_mean
        quiet = solve(ProbioticPondProblem(pd, MicrobiomeTreatment(tolerance = 0.7,
            shedding_reduction = 1.0); coverage = 0.8))
        @test quiet.A_mean >= surv.A_mean                        # less shedding never hurts here
        @test_throws ArgumentError ProbioticPondProblem(pd, prot; stage = :juveniles)

        # Bd in frogs (load tipping model), 80 frogs: above the fold.
        p = BriggsParameters(:briggs2010)
        kw = (c = 100.0, c_treated = 300.0, h = 500.0, N = 80, ε = 1000.0)
        @test critical_coverage(p; kw...) == 14 / 80
        @test isnan(critical_coverage(p; kw..., c_treated = 100.0))
        @test critical_coverage(p; kw..., N = 50) == 0
        frogs(t; cov = 1.0) = solve(ProbioticLoadProblem(p, t; kw..., coverage = cov, days = 400);
            rng = Xoshiro(1))
        un = frogs(MicrobiomeTreatment(); cov = 0.0)
        @test un.alive[end] == 0 && all(==(0), un.treated)
        @test frogs(MicrobiomeTreatment()).alive[end] == 80                 # lasting treatment
        @test frogs(MicrobiomeTreatment(); cov = 0.4).alive[end] == 80      # herd protection
        # one clearance per frog: identical values reproduce the scalar model
        vec = solve(ProbioticLoadProblem(p, MicrobiomeTreatment(); kw..., c = fill(100.0, 80),
            c_treated = fill(300.0, 80), coverage = 0.4, days = 400); rng = Xoshiro(1))
        @test vec.alive == frogs(MicrobiomeTreatment(); cov = 0.4).alive
        @test_throws ArgumentError ProbioticLoadProblem(p, MicrobiomeTreatment(); kw..., c = fill(100.0, 79))
        short = frogs(MicrobiomeTreatment(half_life = 30.0))
        @test short.alive[end] == 0 &&
              minimum(short.death_day) > maximum(un.death_day) + 50        # delayed, not prevented
        @test frogs(MicrobiomeTreatment(half_life = 30.0, direct = 0.1)).alive[end] == 80
    end

    @testset "Bsal with a larval stage (new)" begin
        p = SchmidtParameters(:schmidt2017)
        st = BsalLarvalStage(:wang2021)
        @test st.carry_over == 0.5 && all(r -> r.verified, parameter_records(:wang2021))
        @test_throws ArgumentError BsalLarvalStage(carry_over = 1.5)
        @test disease_free_larvae(p, st) ≈ p.b * p.K * st.larval_period
        # no larval transmission: Schmidt's R0; larvae add to it otherwise
        @test basic_reproduction_number(p, st) ≈ basic_reproduction_number(p)
        inf = BsalLarvalStage(carry_over = 0.5, adult_to_larva = 52.0, larva_to_larva = 5.0)
        @test basic_reproduction_number(p, inf) > basic_reproduction_number(p)
        @test larval_reproduction_number(p, inf) ≈ 5.0 * disease_free_larvae(p, inf) /
              (1 / inf.larval_period + inf.larval_death)
        # without Bsal the population stays at Schmidt's equilibrium
        @test isapprox(solve(BsalLarvalProblem(p, st; I0 = 0.0)).final_adults, 1.0; atol = 1e-6)
        # uninfected larvae are a refuge; full carry-over removes it
        @test solve(BsalLarvalProblem(p, st)).min_adults > 1e-3
        @test solve(BsalLarvalProblem(p, BsalLarvalStage(carry_over = 1.0,
            adult_to_larva = 52.0))).min_adults < 1e-6
        # below Schmidt's invasion threshold, larval transmission lets Bsal invade
        low = Chytridiomycosis.with(p; K = 0.9)
        @test basic_reproduction_number(low) < 1
        @test isapprox(solve(BsalLarvalProblem(low, st; tspan = (0.0, 30.0))).final_adults, 1.0;
            atol = 1e-3)
        both = BsalLarvalStage(carry_over = 0.5, adult_to_larva = 52.0, larva_to_larva = 150.0)
        @test basic_reproduction_number(low, both) > 1
        @test solve(BsalLarvalProblem(low, both; tspan = (0.0, 30.0))).final_adults < 0.9
    end

    @testset "Toads with post-metamorphic infection (new)" begin
        @test all(r -> r.verified, parameter_records(:smith2022))
        @test_throws ArgumentError PostMetamorphicInfection(metamorph_survival = 2.0)
        none = PostMetamorphicInfection()
        for set in ("fig2b", "fig2d")                      # defaults reproduce Mitchell
            p = MitchellParameters(:mitchell2008; set)
            m, r = solve(MitchellPondProblem(p)), solve(ToadStageProblem(p, none))
            @test isapprox(r.A_final, m.A[end]; rtol = 1e-6, atol = 1e-6)
            @test r.toads_persist == m.toads_persist && r.bd_persist == m.bd_persist
            @test r.adult_prevalence == 0
        end
        pd = MitchellParameters(:mitchell2008; set = "fig2d")
        few = solve(ToadStageProblem(pd, PostMetamorphicInfection(metamorph_survival = 0.03)))
        @test few.toads_persist && !few.bd_persist         # a few survivors save the toads
        pc = MitchellParameters(:mitchell2008; set = "fig2c")
        @test !solve(ToadStageProblem(pc, none)).bd_persist
        carr = solve(ToadStageProblem(pc, PostMetamorphicInfection(adult_exposure_days = 14.0,
            adult_relative_shedding = 1.0, adult_clearance = 1.0)))
        @test carr.bd_persist && carr.A_mean < 0.5 * carr.A_disease_free   # adult carriers keep Bd
        @test 0.005 < carr.adult_prevalence < 0.1
        quick = solve(ToadStageProblem(pc, PostMetamorphicInfection(adult_exposure_days = 30.0,
            adult_relative_shedding = 1.0, adult_clearance = 12.0)))
        @test !quick.bd_persist                            # unless adults clear within weeks
    end

    @testset "Microbiome treatment on the stage models (new)" begin
        p = SchmidtParameters(:schmidt2017)
        st = BsalLarvalStage(:wang2021; adult_to_larva = 5.2)
        block = MicrobiomeTreatment(efficacy = 0.99)
        lar(tr; kw...) = solve(ProbioticLarvalProblem(p, st, tr; coverage = 1.0, kw...))
        none = solve(ProbioticLarvalProblem(p, st, MicrobiomeTreatment()))
        @test isapprox(none.min_adults, solve(BsalLarvalProblem(p, st)).min_adults; rtol = 1e-3)
        @test none.final_treated == 0
        q = 0:0.25:9.75
        lost = lar(block; treat = :larvae, treat_times = q)                 # lost at metamorphosis
        kept = lar(block; treat = :larvae, treat_times = q, retention = 1.0)
        @test isapprox(lost.final_adults, none.final_adults; rtol = 0.2) && lost.final_treated < 1e-6
        @test kept.final_adults > 3 * none.final_adults && kept.final_treated > 0.8
        both = lar(block; treat = :both, treat_times = q, retention = 1.0)
        @test both.mean_adults > 0.99                    # both stages, retained: no spread needed
        weeks = lar(MicrobiomeTreatment(efficacy = 0.99, half_life = 0.05); treat = :both,
            treat_times = q, retention = 1.0)
        @test weeks.mean_adults < 0.1                    # a half-life of weeks fails
        @test_throws ArgumentError ProbioticLarvalProblem(p, st, block; treat = :eggs)
        @test_throws ArgumentError ProbioticLarvalProblem(p, st, block; retention = 1.5)

        # Toads when infected metamorphs can survive (Smith 2022): blocking infection no longer
        # removes Bd, and tolerance does better than blocking.
        sm = parameter_set(:smith2022)
        stages = PostMetamorphicInfection(metamorph_survival = 1 - sm.deaths_exposed / sm.n_exposed)
        pd = MitchellParameters(:mitchell2008; set = "fig2d")
        pond(tr) = solve(ProbioticPondProblem(pd, tr; coverage = 0.8, stages))
        blocked, tolerant = pond(MicrobiomeTreatment(efficacy = 0.9)), pond(MicrobiomeTreatment(tolerance = 0.7))
        @test blocked.bd_persist && blocked.A_mean < 0.2 * blocked.A_disease_free
        @test tolerant.A_mean > 0.6 * tolerant.A_disease_free
        # infected adults: untreated, so not with adult treatment; with no treatment the pond
        # model is the toad stage model
        adults = PostMetamorphicInfection(metamorph_survival = stages.metamorph_survival,
            metamorph_exposure_days = 7.0, adult_exposure_days = 30.0,
            adult_relative_shedding = 1.0, adult_clearance = 1.0)
        @test_throws ArgumentError ProbioticPondProblem(pd, MicrobiomeTreatment(); stage = :adults,
            stages = adults)
        ref = solve(ToadStageProblem(pd, adults))
        same = solve(ProbioticPondProblem(pd, MicrobiomeTreatment(); stages = adults))
        @test isapprox(same.A_mean, ref.A_mean; rtol = 1e-4) && same.bd_persist == ref.bd_persist
        @test isapprox(same.A_disease_free, ref.A_disease_free; rtol = 1e-6)
        treated = solve(ProbioticPondProblem(pd, MicrobiomeTreatment(tolerance = 0.7);
            coverage = 0.8, stages = adults))
        @test treated.A_mean > same.A_mean
        chance = solve(ProbioticPondProblem(pd, MicrobiomeTreatment(tolerance = 0.7); coverage = 0.8,
            stages = adults); rng = Xoshiro(2))
        @test isfinite(chance.A_mean) && chance.A_mean >= 0
    end

    @testset "Waning protection (new)" begin
        w = MicrobiomeTreatment(efficacy = 0.9, abundance_half_life = 5 / 365, protective_fold = 100.0)
        @test protection(w, 0.0) == 1 && protection(MicrobiomeTreatment(efficacy = 0.9), 3.0) == 1
        @test protection(w, 1.0) < 1e-6
        d = protection_duration(w)
        @test isapprox(protection(w, d), 0.5; atol = 1e-8)
        @test isapprox(d * 365, 5 * log2(100 * (2 * 1.01 - 1)); rtol = 1e-8)   # ≈ log2(F) half-lives
        @test protection_duration(MicrobiomeTreatment()) == Inf
        @test_throws ArgumentError MicrobiomeTreatment(protective_fold = 0.0)
        # toad tadpoles: one dose wanes within the season; repeated doses restore protection
        pd = MitchellParameters(:mitchell2008; set = "fig2d")
        pond(tr; kw...) = solve(ProbioticPondProblem(pd, tr; coverage = 0.8, kw...))
        lasting = pond(MicrobiomeTreatment(efficacy = 0.9))
        once = pond(w)
        thrice = pond(w; redose_days = [30.0, 60.0])
        @test once.A_mean < thrice.A_mean <= lasting.A_mean * 1.001
        @test_throws ArgumentError pond(w; redose_days = [200.0])
        # Bsal: a waning treatment is as good as none after a year
        ps = SchmidtParameters(:schmidt2017)
        @test solve(ProbioticSLIProblem(ps, MicrobiomeTreatment(efficacy = 0.99,
            abundance_half_life = 5 / 365, protective_fold = 100.0); coverage = 1.0)).final_hosts < 0.01
    end

    @testset "Bacterial abundance per frog (new)" begin
        p = BriggsParameters(:briggs2010)
        kw = (c = 100.0, c_treated = 300.0, h = 500.0, N = 80, ε = 1000.0)
        ab(; half = 5.0, F = 100.0, x = 0.2, kws...) = solve(ProbioticAbundanceProblem(p,
            MicrobiomeTreatment(abundance_half_life = half, protective_fold = F, direct = x); kw..., kws...))
        none = ab(coverage = 0.0)
        old = solve(ProbioticLoadProblem(p, MicrobiomeTreatment(); kw..., days = 365); rng = Xoshiro(1))
        @test none.alive[end] == 0 && maximum(none.death_day) == maximum(old.death_day)   # same die-off
        @test all(isfinite, none.mean_load[1:100])
        same = ab(coverage = 0.75, days = 60)
        @test ab(coverage = 0.75, days = 60, c = fill(100.0, 80), c_treated = fill(300.0, 80)).mean_load ≈
              same.mean_load
        # abundance: the mean declines at the half-life whatever the exchange rate; exchange
        # evens it out between bathed and unbathed frogs
        apart, mixed = ab(coverage = 0.75, x = 0.0, days = 30), ab(coverage = 0.75, x = 0.2, days = 30)
        @test apart.abundance_unbathed[12] == 0 && isapprox(apart.abundance_bathed[11], 2.0^(-2); rtol = 1e-3)
        @test mixed.abundance_unbathed[12] / mixed.abundance_bathed[12] > 0.8
        m(r) = 0.75 * r.abundance_bathed[21] + 0.25 * r.abundance_unbathed[21]
        @test isapprox(m(apart), m(mixed); rtol = 1e-3)
        # one bath only delays the die-off; baths every 90 days prevent it
        once = ab(coverage = 0.75, days = 730)
        @test once.alive[end] == 0 && minimum(filter(>(0), once.death_day)) > maximum(none.death_day) + 20
        @test ab(coverage = 0.75, days = 730, dose_days = 0:90:729).alive[end] == 80
        @test ab(coverage = 0.75, days = 730, dose_days = 0:240:729).alive[end] == 0
        # exchange lets a small bathed group protect the rest
        @test ab(coverage = 0.1, x = 0.0, days = 730, dose_days = 0:90:729).alive[end] < 40
        @test ab(coverage = 0.1, x = 0.2, days = 730, dose_days = 0:90:729).alive[end] == 80
    end

    @testset "Nested model (new)" begin
        p = SchmidtParameters(:schmidt2017)
        nest(bio, tr = MicrobiomeTreatment(); kw...) = solve(NestedChytridProblem(p, bio, tr; kw...))
        @test_throws ArgumentError ChytridBiology(larval_infection = true)
        @test_throws ArgumentError NestedChytridProblem(p, ChytridBiology(); treat = :larvae)
        # with only the latent class it is Schmidt's model, treated or not
        sch = ChytridBiology(latent = true)
        tr = MicrobiomeTreatment(efficacy = 0.99, half_life = 1.0, direct = 5.0, tolerance = 0.3,
            shedding_reduction = 0.2)
        for (t, cov) in ((MicrobiomeTreatment(), 0.0), (tr, 0.7))
            a, b = nest(sch, t; coverage = cov), solve(ProbioticSLIProblem(p, t; coverage = cov))
            @test length(a.adults) == length(b.hosts)
            @test maximum(abs.(a.adults .- b.hosts)) < 1e-4 * p.K     # same trajectory
            @test isapprox(a.final_adults, b.final_hosts; rtol = 1e-3)
        end
        # with larvae and larval infection it is the larval model
        st = BsalLarvalStage(:wang2021; adult_to_larva = 5.2)
        full = ChytridBiology(:wang2021; latent = true, larvae = true, larval_infection = true)
        @test full.carry_over == 0.5
        q = 0:0.25:9.75
        for kw in ((coverage = 0.0,), (coverage = 1.0,),
                (coverage = 1.0, treat = :both, treat_times = q, retention = 1.0))
            a, b = nest(full, tr; kw...), solve(ProbioticLarvalProblem(p, st, tr; kw...))
            @test isapprox(a.mean_adults, b.mean_adults; rtol = 1e-3)
        end
        # growth rate: Schmidt's linearised rate; invasion in every configuration here
        λ = (sqrt((p.e - p.d_I)^2 + 4 * p.e * p.β * p.K) - 2p.b - p.d_I - p.e) / 2
        @test isapprox(invasion_growth_rate(p, sch), λ; rtol = 1e-4)
        @test invasion_growth_rate(p, ChytridBiology()) > invasion_growth_rate(p, sch) > 0
        @test invasion_growth_rate(p, ChytridBiology(environment = true)) <
              invasion_growth_rate(p, ChytridBiology())          # the pool slows the epidemic
        # no pathogen: the population stays at carrying capacity with every switch on
        all_on = ChytridBiology(:wang2021; latent = true, environment = true, load = true,
            recovery = true, larvae = true, larval_infection = true)
        @test isapprox(nest(all_on; I0 = 0.0).final_adults, 1.0; atol = 1e-6)
        @test 0 < nest(all_on).mean_adults < 1
        # fade-out rule: off by default; at Schmidt's density the simplest model's population
        # falls below one animal per 100 km² and is removed, and with a larval refuge it is the
        # pathogen that goes and the host that recovers
        @test isnan(nest(sch).pathogen_fadeout) && isnan(nest(sch).host_extinction)
        @test_throws ArgumentError ChytridBiology(fadeout_density = -1.0)
        lost = nest(ChytridBiology(latent = true, fadeout_density = 1e-4); tspan = (0.0, 30.0))
        @test lost.host_extinction < 5 && lost.final_adults == 0
        saved = nest(ChytridBiology(latent = true, larvae = true, fadeout_density = 1e-4); tspan = (0.0, 30.0))
        @test isnan(saved.host_extinction) && saved.pathogen_fadeout < 5
        @test saved.final_adults > 0.9 && saved.final_prevalence == 0
        # no pathogen: the rule never fires
        quiet = nest(ChytridBiology(latent = true, fadeout_density = 1e-4); I0 = 0.0)
        @test isnan(quiet.pathogen_fadeout) && isapprox(quiet.final_adults, 1.0; atol = 1e-6)
    end

    @testset "Stochastic nested model (new)" begin
        @test_throws ArgumentError StochasticNestedProblem(SchmidtParameters(:schmidt2017),
            ChytridBiology(:wang2021), MicrobiomeTreatment(environmental = 1.0))
        p3 = Chytridiomycosis.with(SchmidtParameters(:schmidt2017); K = 3.0)
        sch = ChytridBiology(latent = true)
        @test_throws ArgumentError StochasticNestedProblem(p3, sch; area = 0.0)
        # reproducible with a seeded generator
        a = solve(StochasticNestedProblem(p3, sch); rng = Xoshiro(7))
        b = solve(StochasticNestedProblem(p3, sch); rng = Xoshiro(7))
        @test a.adults == b.adults && a.infected == b.infected
        # no pathogen: the population stays near carrying capacity and nothing goes extinct
        quiet = solve(StochasticNestedProblem(p3, sch; area = 200.0, introduced = 0); rng = Xoshiro(1))
        @test isapprox(quiet.mean_adults, 1.0; atol = 0.05) && maximum(quiet.infected) == 0
        @test isnan(quiet.pathogen_extinction) && isnan(quiet.host_extinction)
        # a large area recovers the ODE, for Schmidt's model and with every switch on
        full = ChytridBiology(:wang2021; latent = true, environment = true, load = true,
            recovery = true, larvae = true, larval_infection = true)
        for bio in (sch, full)
            ode = solve(NestedChytridProblem(p3, bio)).mean_adults
            sto = solve(StochasticNestedProblem(p3, bio; area = 2000.0, introduced = 20); rng = Xoshiro(2))
            @test isapprox(sto.mean_adults, ode; rtol = 0.1)
        end
        # one infected animal: a major outbreak with probability about 1 − 1/R0 (branching
        # process); here R0 = βK/(b + d_I) ≈ 3
        simple = ChytridBiology()
        R0 = p3.β * p3.K / (p3.b + p3.d_I)
        runs = [solve(StochasticNestedProblem(p3, simple; area = 100.0, tspan = (0.0, 2.0));
            rng = Xoshiro(s)) for s in 1:300]
        major = count(r -> r.peak_infected >= 20, runs) / 300
        @test abs(major - (1 - 1 / R0)) < 0.09
        @test all(r -> r.peak_infected >= 20 || !isnan(r.pathogen_extinction), runs)
        # treatment works through the shared transitions: a lasting block of every adult
        tr = solve(StochasticNestedProblem(p3, sch, MicrobiomeTreatment(efficacy = 1.0);
            coverage = 1.0, introduced = 5); rng = Xoshiro(3))
        @test tr.peak_infected <= 5 && !isnan(tr.pathogen_extinction)
    end

    @testset "Stochastic toads (new)" begin
        T = MicrobiomeTreatment
        pb = MitchellParameters(:mitchell2008; set = "fig2b")
        pd = MitchellParameters(:mitchell2008; set = "fig2d")
        # the deterministic path is unchanged by the refactor into a list of transitions
        m0, t0 = solve(MitchellPondProblem(pb)), solve(ProbioticPondProblem(pb, T()))
        @test isapprox(t0.A_final, m0.A[end]; rtol = 1e-4)
        # whole numbers; reproducible; close to the ODE before Bd arrives
        a = solve(ProbioticPondProblem(pb, T()); rng = Xoshiro(4))
        b = solve(ProbioticPondProblem(pb, T()); rng = Xoshiro(4))
        @test a.A == b.A && all(x -> isinteger(x) && x >= 0, a.A) && all(isinteger, a.Z)
        @test isapprox(a.A_disease_free, t0.A_disease_free; rtol = 0.05)
        # the first years of the Fig. 2d epidemic follow the ODE
        d = solve(ProbioticPondProblem(pd, T()))
        r = solve(ProbioticPondProblem(pd, T()); rng = Xoshiro(1))
        yr5 = findfirst(t -> t >= 5 + pd.t1 - 1e-9, d.t)
        @test isapprox(r.A[yr5], d.A[yr5]; rtol = 0.1)
        # Fig. 2d: without the 50-toad extinction rule, Bd dies out at low toad numbers first
        rs = [solve(ProbioticPondProblem(pd, T()); rng = Xoshiro(s)) for s in 1:10]
        @test !d.toads_persist && all(r -> r.toads_persist && !r.bd_persist, rs)
        # treatment and stages run through the stochastic path
        st = PostMetamorphicInfection(metamorph_survival = 0.29)
        tr = solve(ProbioticPondProblem(pd, T(tolerance = 0.7, abundance_half_life = 5 / 365,
            protective_fold = 100.0); coverage = 0.8, stages = st, redose_days = [30.0, 60.0]); rng = Xoshiro(2))
        @test tr.toads_persist && maximum(tr.tadpoles_treated) > 0
    end

    @testset "Stochastic frogs (new)" begin
        p = BriggsParameters(:briggs2010)
        kw = (c = 100.0, h = 500.0, ε = 1000.0, days = 730)
        frog(n; tr = MicrobiomeTreatment(), s = 1, k...) = solve(StochasticFrogProblem(p, tr; kw..., N = n, k...); rng = Xoshiro(s))
        @test_throws ArgumentError StochasticFrogProblem(p; c = 100.0, h = 0.0, N = 10)
        @test frog(80; s = 3).alive == frog(80; s = 3).alive
        # the same tipping point as the ODE (fold at 67.5 frogs): survive below, die above
        @test all(s -> frog(50; s).alive[end] == 50, 1:5)
        @test all(s -> frog(80; s).alive[end] == 0, 1:5)
        # treatment: 40% treated for life protects all; baths every 90 days protect all
        @test all(s -> frog(80; s, c_treated = 300.0, coverage = 0.4).alive[end] == 80, 1:5)
        bath = MicrobiomeTreatment(abundance_half_life = 5.0, protective_fold = 100.0, direct = 0.1)
        @test all(s -> frog(80; s, tr = bath, c_treated = 300.0, coverage = 0.75,
            dose_days = 0:90:729).alive[end] == 80, 1:3)
        @test frog(80; tr = bath, c_treated = 300.0, coverage = 0.75).alive[end] == 0   # one bath
        # frogs that differ: a partial die-off
        cs = quantile.(LogNormal(log(100.0), 0.5), ((1:80) .- 0.5) ./ 80)
        @test 0 < frog(80; c = cs).alive[end] < 80
    end

    @testset "Fit to Knapp et al. 2022 (new)" begin
        LDP = Chytridiomycosis.LogDensityProblems
        swabs = load_knapp2022()
        @test length(swabs) == 371
        @test count(s -> s.group === :treated, swabs) == 246
        @test count(s -> s.group === :control, swabs) == 44
        @test count(s -> s.group === :wild, swabs) == 81
        @test sort(unique(s.day for s in swabs)) == [-7, 0, 12, 37, 56]
        @test_throws ArgumentError KnappLoadModel(:other)

        m = KnappLoadModel(:applied)
        @test parameter_names(m)[end] === :K
        @test LDP.dimension(m) == 12
        @test LDP.capabilities(KnappLoadModel) isa LDP.LogDensityOrder{0}
        x = knapp_start(m)
        @test isfinite(LDP.logdensity(m, x))
        @test_throws DimensionMismatch LDP.logdensity(m, x[1:3])

        # equilibrium: the positive root, and no suppression gives ι/(−r0)
        r0 = m.r0
        L = Chytridiomycosis.knapp_equilibrium(1e5, r0, 2e5, 1e3)
        @test 1e5 + r0 * L - 2e5 * L / (1e3 + L) ≈ 0 atol = 1e-6
        @test Chytridiomycosis.knapp_equilibrium(1e5, r0, 0.0, 1e3) ≈ 1e5 / -r0

        # the exponential integrator agrees with an adaptive ODE solver
        md = KnappLoadModel(:decay)
        xd = knapp_start(md)
        q = Chytridiomycosis.knapp_params(md, xd)
        f(u, _, t) = q.ι + r0 * u - q.Δc * exp(-t / q.τ) * u / (q.h + u)
        sol = solve(Chytridiomycosis.ODEProblem(f, q.L0, (0.0, 56.0)), Chytridiomycosis.Tsit5();
            saveat = [12.0, 37.0, 56.0], reltol = 1e-10, abstol = 1e-10)
        pred = predicted_loads(md, xd)
        @test pred.treated.bd[3:5] ≈ sol.u rtol = 1e-3
        @test all(pred.wild.bd .≈ pred.wild.bd[1])          # untreated: at equilibrium

        # fits: the mechanisms nest :none, so neither does worse
        fn = maximum_likelihood(KnappLoadModel(:none); starts = 1, iters = 1500)
        fd = maximum_likelihood(md; starts = 1, iters = 1500)
        @test fd.loglik > fn.loglik + 50                     # treated frogs were suppressed
        @test fd.aic ≈ 2 * 12 - 2 * fd.loglik
        @test fd.estimates.L0 < 100
        prof = profile_likelihood(fd, :τ, fd.estimates.τ .* [0.5, 1.0, 2.0]; iters = 300)
        @test prof[2] ≈ fd.loglik atol = 0.5
        @test all(prof .<= fd.loglik + 0.5)
    end

    # Phase 0 foundations
    include("units.jl")
    include("thermal.jl")
    include("validation_targets.jl")
    include("n5_identifiability.jl")
    include("knapp_bayes.jl")
end

include("schmidt2017_spatial.jl")
include("briggs2010_stages.jl")
include("daversa2022.jl")
include("rumschlag2022.jl")
include("miller2025.jl")
include("bozzuto2019.jl")
include("bsal_transmission.jl")
