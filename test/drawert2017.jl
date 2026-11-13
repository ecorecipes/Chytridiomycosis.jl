# Drawert et al. 2017 (J R Soc Interface 14:20170480): stochastic individual-based Bd model
# with one-time interventions. The oracles are:
# - the paper (§3, Tables 1–2, §4) and the SI text (ESM §1, Table S3, eqs 10–12, Table S9);
# - the authors' C code and their table of 2132 points with their own 30-run summaries
#   (figshare 10.6084/m9.figshare.5328127, bundled as data/drawert2017/baseline_points.csv);
# - exact results for linear birth–death processes, and the Briggs et al. (2010) ODE, which
#   the Baseline model's means must follow while no frog dies (an internal check);
# - the stored output of the reproduction scripts (data/drawert2017/reproduction.toml,
#   written by scripts/drawert2017_summarise.jl from scripts/drawert2017_*.jl).
using Chytridiomycosis
using Test
using Random: Xoshiro
using Distributions: Binomial, mean, std, var
using TOML

@testset "Drawert 2017 stochastic IBM" begin
    C = Chytridiomycosis
    q = parameter_set(:drawert2017)
    pts = load_drawert2017()
    rep = TOML.parsefile(joinpath(dirname(drawert2017_path()), "reproduction.toml"))

    @testset "parameters (Table 1, §3.3–3.4, ESM §1.4, Table S9)" begin
        @test all(r -> r.verified, parameter_records(:drawert2017))
        @test (q.Smax, q.S0, q.Z0, q.V, q.tfinal, q.n_runs) == (10000, 100, 0, 1, 1000, 30)
        @test (q.N_min, q.N_max, q.γ_min, q.γ_max) == (5, 200, 1e-6, 1)
        @test (q.η_min, q.η_max, q.σ0_min, q.σ0_max, q.μ_min, q.μ_max) == (5, 20, 0.1, 0.5, 0.01, 1.5)
        @test (q.ν0_min, q.ν0_max, q.f_min, q.f_max) == (0, 1, 0, 1)
        @test (q.λ_min, q.λ_max, q.ε, q.n_points) == (0.035, 0.35, 0.09, 2132)
        @test (q.S_star_low, q.S_star_high) == (1.2, 2.0)
        @test time_unit(DrawertParameters) === :day
        ivs = drawert_interventions()
        @test length(ivs) == 27
        @test Set(iv.day for iv in ivs) == Set([14.0, 28.0, 56.0])
        @test Set(iv.fraction for iv in ivs) == Set([0.25, 0.5, 0.75])
        @test [iv.kind for iv in ivs[1:9:27]] == [:clean_frogs, :cull_frogs, :clean_environment]
        @test_throws ArgumentError DrawertIntervention(:burn, 14, 0.5)
        @test_throws ArgumentError DrawertIntervention(:cull_frogs, 14, 1.5)
        # Table S9 prints γ to 3 significant figures and the rest to 2–3 decimal places
        for (g, id) in zip(("groupA", "groupB", "groupC", "groupD", "groupE"),
                           (120, 1094, 182, 982, 192))
            s = parameter_set(:drawert2017; set = g)
            @test s.point_id == id
            r = pts[id]
            @test r.id == id && s.N == r.N
            for (a, b) in ((s.γ, r.γ), (s.η, r.η), (s.ν0, r.ν0), (s.f, r.f), (s.σ0, r.σ0), (s.μ, r.μ))
                @test isapprox(a, b; rtol = 6e-3) || abs(a - b) <= 5e-4
            end
            p = DrawertParameters(:drawert2017; set = g)
            @test p.N == s.N && p.γ == s.γ && p.Smax == 10000 && p.S0 == 100
        end
        p = DrawertParameters(:drawert2017; point = 182)
        @test (p.N, p.γ, p.η) == (pts[182].N, pts[182].γ, pts[182].η)
        @test model_variant(p) === :baseline
        @test DrawertParameters(:drawert2017; point = 182, N = 10).N == 10
        @test_throws ArgumentError DrawertParameters(:drawert2017; point = 99999)
        @test_throws ArgumentError DrawertParameters(:drawert2017; set = "groupA", point = 1)
    end

    @testset "the authors' table reproduces the paper (§4, Fig. 2a, Fig. 3)" begin
        @test length(pts) == q.n_points && [r.id for r in pts] == 1:2132
        n = length(pts)
        cls = [outcome_class(r.frac_frogs_extinct, r.frac_fungus_extinct, r.frac_persistence)
               for r in pts]
        @test round(count(==(:frog_extinction), cls) / n; digits = 2) == 0.65    # "65%"
        @test round(count(==(:fungal_extinction), cls) / n; digits = 2) == 0.24  # "24%"
        @test round(count(==(:variable), cls) / n; digits = 2) == 0.10           # "approximately 10%"
        # "coexisting in about 0.01% of the cases": the table has 3 points, 0.14%
        @test count(==(:coexistence), cls) == 3
        @test_broken round(count(==(:coexistence), cls) / n; digits = 4) == 0.0001
        grp = [survival_group(r.mean_surviving, r.std_surviving) for r in pts]
        @test [count(==(g), grp) for g in (:A, :B, :C, :D, :E)] == [1386, 109, 221, 396, 20]
        @test round(count(==(:A), grp) / n; digits = 2) == 0.65                  # Group A "65%"
        # Table S9's points are in the groups they represent
        @test [grp[id] for id in (120, 1094, 182, 982, 192)] == [:A, :B, :C, :D, :E]
        @test outcome_class(1, 0, 0) === :frog_extinction && outcome_class(0.5, 0.5, 0) === :variable
        @test survival_group(0.0, 0.0) === :A && survival_group(0.05, 0.02) === :B
        @test survival_group(0.4, 0.2) === :C && survival_group(0.9, 0.05) === :D
        @test survival_group(0.3, 0.05) === :E
        @test intervention_effect(0.35) == 0.3 && intervention_effect(-0.1) == -0.1
        @test intervention_effect(0.05) == 0 && intervention_effect(-0.7) == -0.5
    end

    @testset "growth rate, non-dimensional groups and sampling (ESM §1.4)" begin
        λ = [load_growth_rate(DrawertParameters(:drawert2017; point = r.id)) for r in pts]
        @test maximum(abs.(λ .- [r.λ for r in pts])) < 1.5e-6        # their 6 significant figures
        @test all(q.λ_min .<= λ .<= q.λ_max)
        p = DrawertParameters(:drawert2017; point = 1)
        @test load_growth_rate(p) == load_growth_rate(BriggsParameters(p), p.N)
        b = BriggsParameters(p)
        @test (b.ν, b.σ, b.γ, b.Smax, b.V) == (p.ν0, p.σ0, p.γ, 10000.0, 1.0)
        c = nondimensional_parameters(p)
        @test c.c1 ≈ p.η * p.ν0 * p.f / p.σ0 && c.c4 ≈ p.μ / p.σ0
        @test c.c2 ≈ p.γ * p.η * (1 - p.f) * p.ν0 * p.N / p.σ0^2 && c.c3 ≈ p.γ * p.N / p.σ0
        # ESM eq. (10): with the range of log c2, log c3 as scale the points are spaced just
        # above ε = 0.09; the printed log(c_range) leaves them much farther apart.
        cs = [nondimensional_parameters(DrawertParameters(:drawert2017; point = r.id)) for r in pts]
        nn(r; kw...) = minimum(minimum(nondimensional_distance(cs[i], cs[j], r; kw...)
                                       for j in eachindex(cs) if j != i) for i in eachindex(cs))
        dmin = nn(nondimensional_ranges(cs))
        @test q.ε < dmin < q.ε + 0.0005
        @test nn(nondimensional_ranges(cs; as_printed = true); as_printed = true) > q.ε + 0.001
        @test rep["deterministic"]["nearest_neighbour_min"] == round(dmin; digits = 4)
    end

    @testset "Nu and Sigma models (ESM §1.1.2–1.1.3, eqs 11–12)" begin
        p = DrawertParameters(:drawert2017; set = "groupA")
        @test equilibrium_load_ratio(p) == Inf
        K = p.γ * p.η * (1 - p.f) * p.N / (p.γ * p.N + p.μ) + p.η * p.f
        for S in (q.S_star_low, q.S_star_high)
            pn = drawert_variant(p, :nu; S_star = S)
            ps = drawert_variant(p, :sigma; S_star = S)
            @test model_variant(pn) === :nu && model_variant(ps) === :sigma
            @test equilibrium_load_ratio(pn) ≈ S && equilibrium_load_ratio(ps) ≈ S
            # ESM §1.1.2: S* = ln[(ν0/σ)(γη(1−f)N/(γN + μ) + ηf)]/(Smax ν1)
            @test log(p.ν0 / p.σ0 * K) / (p.Smax * pn.ν1) ≈ S
            # ESM eq. (12)
            @test ps.σ1 ≈ (p.ν0 * K - p.σ0) / (S * p.Smax)
            @test encystment_probability(pn, 0) == p.ν0
            @test encystment_probability(pn, 1000) ≈ p.ν0 * exp(-1000pn.ν1)
            @test sporangium_loss_rate(ps, 1000) ≈ p.σ0 + 1000ps.σ1
            # loads level off at S*·Smax: zero net growth there (quasi-steady pool)
            Sx = S * p.Smax
            @test encystment_probability(pn, Sx) * K ≈ p.σ0
            @test p.ν0 * K ≈ sporangium_loss_rate(ps, Sx)
        end
        @test drawert_variant(drawert_variant(p, :nu), :baseline) == p
        # ESM eq. (11) as printed (no logarithm, no ηf term) does not give S* = 1.2
        pp = drawert_variant(p, :nu; S_star = 1.2, as_printed = true)
        @test pp.ν1 ≈ p.ν0 * p.γ * p.η * (1 - p.f) * p.N / (p.σ0 * (p.γ * p.N + p.μ)) / 12000
        @test_broken isapprox(equilibrium_load_ratio(pp), 1.2; rtol = 0.01)
        @test equilibrium_load_ratio(pp) ≈ 1.4283 rtol = 1e-4       # regression (Group A point)
        @test_throws ArgumentError drawert_variant(p, :both)
        @test_throws ArgumentError drawert_variant(DrawertParameters(:drawert2017; set = "groupA",
                                                                     ν0 = 0.001), :nu)
        @test_throws ArgumentError DrawertIBMProblem(DrawertParameters(:drawert2017;
                                                                       set = "groupA", ν1 = 1e-4, σ1 = 1e-4))
    end

    @testset "deterministic model (eqs 3.1–3.2) = Briggs 2010" begin
        p = DrawertParameters(:drawert2017; point = 982)
        d = solve(DrawertODEProblem(p))
        b = solve(BriggsLoadProblem(BriggsParameters(p); N = p.N, S0 = [p.S0; zeros(p.N - 1)],
                                    tspan = (0.0, 1000.0)); abstol = 1e-8, reltol = 1e-8)
        @test d.outcome === :frogs_extinct && d.survivors == 0 && all(!isnan, d.death_time)
        @test maximum(abs.(d.death_time .- b.death_time)) < 1e-6
        @test to_ode_problem(DrawertODEProblem(p)).u0 == [p.S0; zeros(p.N - 1); 0.0]
        # "deterministic … always extinction — by construction" (Fig. 2a), here for the five
        # representative points and all variants; all 2132 points in reproduction.toml
        for g in ("groupA", "groupB", "groupC", "groupD", "groupE")
            pg = DrawertParameters(:drawert2017; set = g)
            @test solve(DrawertODEProblem(pg)).outcome === :frogs_extinct
            @test solve(DrawertODEProblem(drawert_variant(pg, :sigma; S_star = 1.2))).outcome ===
                  :frogs_extinct
        end
        dt = rep["deterministic"]
        @test dt["n_points"] == 2132
        for k in ("frogs_extinct_baseline", "frogs_extinct_nu_1.2", "frogs_extinct_nu_2",
                  "frogs_extinct_sigma_1.2", "frogs_extinct_sigma_2")
            @test dt[k] == 2132
        end
        # §3.2's "death of all frogs within 14–275 days" is the first death (in whole days)
        @test (dt["first_death_day_min"], dt["first_death_day_max"]) == (q.t_death_min, q.t_death_max)
        @test_broken dt["last_death_max"] <= q.t_death_max
        @test dt["last_death_outside_window"] == 71
        @test dt["printed_nu_Sstar_below_1_at_1.2"] == 270
        # interventions follow determ_frogs.c: the first floor(fN) live frogs, index order
        pc = DrawertParameters(; N = 4, γ = 0.0, η = 10.0, ν0 = 0.5, f = 0.5, σ0 = 1.0, μ = 1.0,
                               Smax = 10000, S0 = 100)
        r = solve(DrawertODEProblem(pc; intervention = DrawertIntervention(:clean_frogs, 1, 0.25)))
        @test r.outcome === :fungus_extinct && r.survivors == 4      # frog 1 cleaned on day 1
        r = solve(DrawertODEProblem(pc; intervention = DrawertIntervention(:cull_frogs, 1, 0.25)))
        @test r.survivors == 3 && r.outcome === :fungus_extinct
        r = solve(DrawertODEProblem(pc; intervention = DrawertIntervention(:clean_environment, 1, 0.5)))
        @test !isnan(r.death_time[1]) && r.survivors == 3            # γ = 0: frog 1 alone
    end

    @testset "stochastic model: mechanics" begin
        p = DrawertParameters(:drawert2017; point = 982)
        prob = DrawertIBMProblem(p; saveat = 0:10.0:60, save_loads = true)
        r1 = solve(prob; rng = Xoshiro(1))
        r2 = solve(prob; rng = Xoshiro(1))
        @test r1.events == r2.events && r1.loads == r2.loads && r1.outcome == r2.outcome
        @test r1.outcome in C.DRAWERT_OUTCOMES
        @test r1.frogs[1] == p.N && r1.sporangia[1] == p.S0 && r1.zoospores[1] == 0
        @test r1.loads[:, 1] == [p.S0; zeros(Int, p.N - 1)]
        @test all(r1.loads .<= p.Smax)
        @test vec(sum(max.(r1.loads, 0); dims = 1)) == r1.sporangia
        @test vec(count(>=(0), r1.loads; dims = 1)) == r1.frogs
        @test r1.survivors == count(isnan, r1.death_time) && !any(r1.culled)
        @test r1.surviving_fraction == r1.survivors / p.N
        @test r1.first_loss == minimum(filter(!isnan, r1.death_time); init = Inf)
        @test_throws ArgumentError DrawertIBMProblem(p; saveat = [2000.0])
        @test_throws ArgumentError DrawertIBMProblem(DrawertParameters(:drawert2017;
                                                                       point = 1, f = 1.5))

        # A frog dies when an infection would take it above Smax (code; main text S_i > Smax).
        # One frog, no loss: a pure birth process from S0 to Smax + 1.
        p1 = DrawertParameters(; N = 1, γ = 0.0, η = 1.0, ν0 = 1.0, f = 1.0, σ0 = 0.0, μ = 0.0,
                               Smax = 50, S0 = 5)
        r = solve(DrawertIBMProblem(p1; tfinal = 100.0); rng = Xoshiro(2))
        @test r.outcome === :frogs_extinct && r.events == 46 && r.t_end == r.death_time[1]
        # Yule process: E[time from 5 to 51] = Σ_{k=5}^{50} 1/k
        ts = [solve(DrawertIBMProblem(p1; tfinal = 100.0); rng = Xoshiro(s)).death_time[1]
              for s in 1:2000]
        @test abs(mean(ts) - sum(1 / k for k in 5:50)) < 4 * std(ts) / sqrt(2000)

        # With f = 1 nothing reaches the pool, and self-reinfection vs loss is a linear
        # birth–death process on one frog: birth rate η f ν0, death rate σ0 per sporangium.
        # From S0 = 1, P(fungus lost) = σ0/(η f ν0) (before Smax is reached, almost surely).
        pb = DrawertParameters(; N = 1, γ = 0.0, η = 1.2, ν0 = 0.5, f = 1.0, σ0 = 0.3, μ = 0.0,
                               Smax = 200, S0 = 1)
        outs = [solve(DrawertIBMProblem(pb); rng = Xoshiro(s)).outcome for s in 1:4000]
        pl = 0.3 / (1.2 * 1.0 * 0.5)
        @test abs(count(==(:fungus_extinct), outs) / 4000 - pl) < 4 * sqrt(pl * (1 - pl) / 4000)
        @test all(o -> o in (:fungus_extinct, :frogs_extinct), outs)
        r = solve(DrawertIBMProblem(pb; saveat = 0:1.0:20); rng = Xoshiro(3))
        @test all(==(0), r.zoospores)

        # The pool alone: ν0 = 0, σ0 = 0 keeps every load fixed, and Z is an immigration–death
        # process with immigration η(1−f)S0 and death μ + γN: Z(t) ~ Poisson with mean
        # m(t) = η(1−f)S0 (1 − e^{−(μ+γN)t})/(μ + γN).
        pz = DrawertParameters(; N = 5, γ = 0.2, η = 3.0, ν0 = 0.0, f = 0.4, σ0 = 0.0, μ = 0.5,
                               Smax = 10000, S0 = 20)
        zs = [solve(DrawertIBMProblem(pz; tfinal = 3.0, saveat = [3.0]); rng = Xoshiro(s)).zoospores[1]
              for s in 1:3000]
        m = 3.0 * 0.6 * 20 * (1 - exp(-(0.5 + 1.0) * 3)) / (0.5 + 1.0)
        @test abs(mean(zs) - m) < 4 * sqrt(m / 3000)
        @test isapprox(var(zs), m; rtol = 0.1)                       # Poisson
    end

    @testset "stochastic Baseline means follow the Briggs ODE (internal check)" begin
        # The Baseline propensities are linear in (S_i, Z) while all frogs live, so the means
        # of the total load and of Z obey the Briggs et al. (2010) ODE exactly until a frog
        # dies. Smax is set high so that none dies in 12 days.
        p = DrawertParameters(; N = 20, γ = 0.02, η = 6.0, ν0 = 0.3, f = 0.3, σ0 = 0.35, μ = 0.8,
                              Smax = 10^7, S0 = 100)
        ts = [4.0, 8.0, 12.0]
        n = 600
        S = zeros(n, 3)
        Z = zeros(n, 3)
        for s in 1:n
            r = solve(DrawertIBMProblem(p; tfinal = 12.0, saveat = ts); rng = Xoshiro(100 + s))
            S[s, :] = r.sporangia
            Z[s, :] = r.zoospores
        end
        b = solve(BriggsLoadProblem(BriggsParameters(p); N = p.N, S0 = [p.S0; zeros(p.N - 1)],
                                    tspan = (0.0, 12.0)); saveat = ts, abstol = 1e-10, reltol = 1e-10)
        Sode = [sum(u[1:end-1]) for u in b.sol.u]
        Zode = [u[end] for u in b.sol.u]
        for k in 1:3
            @test abs(mean(S[:, k]) - Sode[k]) < 4 * std(S[:, k]) / sqrt(n)
            @test abs(mean(Z[:, k]) - Zode[k]) < 4 * std(Z[:, k]) / sqrt(n)
        end
        @test load_growth_rate(p) > 0
    end

    @testset "stochastic Nu and Sigma models: loads level off" begin
        # One frog that reinfects only itself (γ = 0, f = 1): growth ν(S)η − σ(S) stops at
        # S = ln(ην0/σ0)/ν1 (Nu) or (ην0 − σ0)/σ1 (Sigma). Smax is far above.
        base = (; N = 1, γ = 0.0, η = 2.0, ν0 = 0.5, f = 1.0, σ0 = 0.5, μ = 0.0, Smax = 10^6, S0 = 50)
        pn = DrawertParameters(; base..., ν1 = log(2.0) / 400)       # equilibrium 400
        ps = DrawertParameters(; base..., σ1 = 0.5 / 400)            # equilibrium 400
        for p in (pn, ps)
            r = [solve(DrawertIBMProblem(p; tfinal = 80.0, saveat = [80.0]); rng = Xoshiro(s)).sporangia[1]
                 for s in 1:200]
            @test abs(mean(r) - 400) < 12
            @test all(>(0), r)
        end
    end

    @testset "interventions (§3.3; FrogAdditions.c)" begin
        p = DrawertParameters(; N = 40, γ = 0.01, η = 8.0, ν0 = 0.2, f = 0.25, σ0 = 0.3, μ = 0.5,
                              Smax = 10000, S0 = 100)
        ts = [9.0, 10.0]
        # cleaning all frogs and the whole pool on day 10 ends the epidemic then
        r = solve(DrawertIBMProblem(p; intervention = DrawertIntervention(:clean_frogs, 10, 1.0),
                                    saveat = ts, save_loads = true); rng = Xoshiro(4))
        @test r.sporangia[1] > 0 && r.sporangia[2] == 0 && all(r.loads[:, 2] .<= 0)
        r = solve(DrawertIBMProblem(p; intervention = DrawertIntervention(:clean_environment, 10, 0.3),
                                    as_coded = true, saveat = ts); rng = Xoshiro(4))
        @test r.zoospores[1] > 0 && r.zoospores[2] == 0               # as coded: all zoospores
        # the paper's version removes Binomial(Z, fraction): compare the same path
        zb = Int[]
        za = Int[]
        for s in 1:400
            a = solve(DrawertIBMProblem(p; intervention = DrawertIntervention(:clean_environment, 10, 0.3),
                                        saveat = [10.0], tfinal = 10.0); rng = Xoshiro(s))
            b = solve(DrawertIBMProblem(p; saveat = [10.0], tfinal = 10.0); rng = Xoshiro(s))
            push!(za, a.zoospores[1]); push!(zb, b.zoospores[1])
        end
        @test mean(za) / mean(zb) ≈ 0.7 rtol = 0.03
        @test all(za .<= zb)
        # culling removes Binomial(N, fraction) frogs; they do not survive
        culled = [count(solve(DrawertIBMProblem(p; intervention = DrawertIntervention(:cull_frogs, 10, 0.25),
                                                tfinal = 10.0); rng = Xoshiro(s)).culled) for s in 1:400]
        @test mean(culled) ≈ 10 atol = 4 * sqrt(40 * 0.25 * 0.75 / 400)
        @test var(culled) ≈ 40 * 0.25 * 0.75 rtol = 0.25
        rc = solve(DrawertIBMProblem(p; intervention = DrawertIntervention(:cull_frogs, 10, 1.0));
                   rng = Xoshiro(5))
        @test rc.outcome === :frogs_extinct && all(rc.culled) && rc.first_loss <= 10
        # an intervention after the outcome is settled is not applied (the C code stops)
        pl = DrawertParameters(; N = 5, γ = 0.0, η = 1.0, ν0 = 0.1, f = 0.5, σ0 = 2.0, μ = 5.0,
                               Smax = 10000, S0 = 3)
        r = solve(DrawertIBMProblem(pl; intervention = DrawertIntervention(:cull_frogs, 100, 1.0));
                  rng = Xoshiro(6))
        @test r.outcome === :fungus_extinct && r.survivors == 5 && !any(r.culled)
    end

    @testset "reproduction against the authors' runs (§4, Fig. 2a, Fig. 3)" begin
        # One cheap point live: Group D representative (ESM Table S9, id 982), 30 runs, vs the
        # authors' 30 runs (mean 0.942, sd 0.088, fungal extinction in all 30).
        a = pts[982]
        s = drawert_summary(DrawertIBMProblem(DrawertParameters(:drawert2017; point = 982));
                            rng = Xoshiro(982))
        @test s.nruns == 30 && length(s.surviving) == 30
        @test s.group === :D && s.outcome_class === :fungal_extinction
        @test abs(s.mean_surviving - a.mean_surviving) < 3 * sqrt((s.std_surviving^2 +
                                                                    a.std_surviving^2) / 29)
        @test abs(s.mean_days_to_first_death - a.mean_days_to_first_death) < 2

        # Stored: 116 points (stratified by the authors' group), 30 runs each.
        b = rep["baseline_sample"]
        n = b["n_points"]
        @test n == 116 && b["runs_per_point"] == 30
        # same outcome class (Fig. 2a) and group (Fig. 3) as the authors' runs at most points
        @test b["same_class"] / n > 0.85
        @test b["same_group"] / n > 0.8
        # per-point two-sample tests reject at about their nominal 5% rate
        @test b["fisher_frogs_extinct_p_below_0.05"] / n < 0.1
        @test b["welch_mean_surviving_p_below_0.05"] / n < 0.12
        @test abs(b["mean_surviving_difference_mean"]) < 0.02
        @test b["mean_surviving_correlation"] > 0.97
        # weighted estimates of the population fractions (paper: 65%, 24%, ~10%)
        cw = b["class_weighted_ours"]
        @test abs(cw["frog_extinction"] - 0.65) < 0.05
        @test abs(cw["fungal_extinction"] - 0.24) < 0.07
        @test abs(cw["variable"] - 0.10) < 0.07
        @test abs(b["group_weighted_ours"]["A"] - 0.65) < 0.05
    end

    @testset "stored Nu/Sigma and intervention results (§4, §4.1, Fig. 2a, Fig. 4)" begin
        v = rep["variants"]
        @test v["runs_per_point"] == 10
        # load regulation makes coexistence possible (Fig. 2a: ~7–8% of Nu and Sigma points,
        # ~0.1% of Baseline points)
        @test v["random_nu_1.2"]["runs_persistence"] + v["random_sigma_1.2"]["runs_persistence"] > 0
        iv = rep["interventions"]
        # Group D (id 982): culling kills frogs that would have survived; cleaning has
        # little effect (§4.1)
        g = iv["grid_982"]
        @test g["cull_frogs"]["negative"] >= 7
        @test g["clean_frogs"]["negative"] == 0
        # cleaning the environment as coded ignores the fraction removed
        @test g["clean_environment_as_coded"]["positive"] == 0
    end
end
