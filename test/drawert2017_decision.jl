# Additions to the Drawert et al. (2017) model for decision analysis (package extensions; not
# published results). The oracles are:
# - the published model: no mortality, one pulse as a vector and no pulse give the same runs as
#   the plain problem (test/drawert2017.jl checks those against the paper);
# - exact identities: Binomial counts, conservation of frogs, common random numbers up to the
#   first intervention, two pulses against one;
# - the observation distributions of swab_edna.jl, and the importance-weighting property that
#   data have the higher likelihood under the parameters that made them;
# - the stored table (data/drawert2017/decision_*.csv): complete design, shared seeds, the rule
#   for `handled`, and a live re-run of some of its cells.
using Chytridiomycosis
using Test
using TOML
using Random: Xoshiro
using Distributions: mean, var

@testset "Drawert 2017 for decision analysis" begin
    C = Chytridiomycosis
    p = DrawertParameters(; N = 40, γ = 0.01, η = 8.0, ν0 = 0.2, f = 0.25, σ0 = 0.3, μ = 0.5,
                          Smax = 10000, S0 = 100)

    @testset "animals handled" begin
        for kind in C.DRAWERT_INTERVENTION_KINDS, s in 1:20
            r = solve(DrawertIBMProblem(p; intervention = DrawertIntervention(kind, 10, 0.5),
                                        tfinal = 10.0); rng = Xoshiro(s))
            q = only(r.pulses)
            @test q.applied && (q.kind, q.day, q.fraction) === (kind, 10.0, 0.5)
            @test 0 <= q.handled <= q.frogs && r.handled == q.handled
            @test r.treatment_deaths == 0 && !any(r.treatment_death)
            kind === :cull_frogs && @test q.handled == count(r.culled)
            kind === :clean_environment && @test q.handled == 0 && 0 <= q.zoospores_removed <= q.zoospores
            kind === :clean_environment || @test q.zoospores_removed == 0
        end
        all_ = solve(DrawertIBMProblem(p; intervention = DrawertIntervention(:clean_frogs, 10, 1)); rng = Xoshiro(1))
        @test all_.handled == all_.pulses[1].frogs
        hs = [solve(DrawertIBMProblem(p; intervention = DrawertIntervention(:clean_frogs, 1, 0.75),
                                      tfinal = 1.0); rng = Xoshiro(s)).handled for s in 1:400]
        @test mean(hs) ≈ 30 atol = 4 * sqrt(40 * 0.75 * 0.25 / 400)
        @test solve(DrawertIBMProblem(p); rng = Xoshiro(1)).pulses == []
        # not applied once the outcome is settled
        pl = DrawertParameters(; N = 5, γ = 0.0, η = 1.0, ν0 = 0.1, f = 0.5, σ0 = 2.0, μ = 5.0,
                               Smax = 10000, S0 = 3)
        r = solve(DrawertIBMProblem(pl; intervention = DrawertIntervention(:cull_frogs, 100, 1)); rng = Xoshiro(6))
        @test !r.pulses[1].applied && r.pulses[1].frogs == 5 && r.handled == 0
    end

    @testset "treatment mortality" begin
        clean(m) = DrawertIntervention(:clean_frogs, 10, 0.5; mortality = m)
        @test clean(0.0) == DrawertIntervention(:clean_frogs, 10, 0.5)
        a = solve(DrawertIBMProblem(p; intervention = clean(0.0)); rng = Xoshiro(3))
        b = solve(DrawertIBMProblem(p; intervention = DrawertIntervention(:clean_frogs, 10, 0.5)); rng = Xoshiro(3))
        @test a.events == b.events && a.survivors == b.survivors
        r = solve(DrawertIBMProblem(p; intervention = clean(1.0), tfinal = 10.0); rng = Xoshiro(3))
        @test r.treatment_deaths == r.handled == count(r.treatment_death) > 0
        @test all(isnan, r.death_time[r.treatment_death])           # apart from disease deaths
        @test r.survivors == p.N - r.handled - count(!isnan, r.death_time)
        @test r.first_loss <= 10
        ds = [(x = solve(DrawertIBMProblem(p; intervention = clean(0.39), tfinal = 10.0); rng = Xoshiro(s));
               (x.treatment_deaths, x.handled)) for s in 1:300]
        n = sum(last, ds)
        @test sum(first, ds) / n ≈ 0.39 atol = 4 * sqrt(0.39 * 0.61 / n)
        @test_throws ArgumentError DrawertIntervention(:cull_frogs, 10, 0.5; mortality = 0.1)
        @test_throws ArgumentError DrawertIntervention(:clean_frogs, 10, 0.5; mortality = 1.1)
        @test_throws ArgumentError DrawertODEProblem(p; intervention = clean(0.1))
    end

    @testset "several interventions" begin
        one = DrawertIntervention(:clean_frogs, 10, 0.5)
        a = solve(DrawertIBMProblem(p; intervention = one); rng = Xoshiro(4))
        b = solve(DrawertIBMProblem(p; intervention = [one]); rng = Xoshiro(4))
        @test a.events == b.events && a.pulses == b.pulses
        none = solve(DrawertIBMProblem(p; intervention = DrawertIntervention[]); rng = Xoshiro(4))
        @test none.events == solve(DrawertIBMProblem(p); rng = Xoshiro(4)).events
        # frogs and pool cleaned completely on the same day (in either order): Bd is gone
        both = [DrawertIntervention(:clean_environment, 10, 1), DrawertIntervention(:clean_frogs, 10, 1)]
        for iv in (both, reverse(both))
            r = solve(DrawertIBMProblem(p; intervention = iv); rng = Xoshiro(4))
            @test r.outcome === :fungus_extinct && r.t_end == 10
            @test [q.kind for q in r.pulses] == [x.kind for x in iv]
            @test all(q.applied for q in r.pulses)
        end
        # sorted by day; repeated culls add up
        r = solve(DrawertIBMProblem(p; intervention = [DrawertIntervention(:cull_frogs, 8, 0.25),
            DrawertIntervention(:cull_frogs, 4, 0.25)], tfinal = 8.0); rng = Xoshiro(5))
        @test [q.day for q in r.pulses] == [4, 8]
        @test r.handled == count(r.culled) == sum(q.handled for q in r.pulses)
        @test r.pulses[2].frogs <= p.N - r.pulses[1].handled
    end

    @testset "infection state" begin
        pr(; kw...) = DrawertIBMProblem(p; tfinal = 60.0, kw...)
        base = solve(pr(); rng = Xoshiro(8))
        r = solve(pr(; saveat = [5.0, 12.0], save_loads = true, evaluate = [5.0, 12.0]); rng = Xoshiro(8))
        @test r.events == base.events                        # evaluation does not change the run
        for k in 1:2
            @test r.evaluation[k] == (; day = r.t[k], uninfected = count(==(0), r.loads[:, k]),
                                      infected = count(>(0), r.loads[:, k]))
        end
        h = solve(pr(; saveat = [12.0], save_loads = true, evaluate = [12.0], infection_threshold = 200); rng = Xoshiro(8))
        @test h.evaluation[1].infected == count(>(200), h.loads[:, 1])
        for iv in [nothing; DrawertIntervention(:cull_frogs, 10, 0.5); DrawertIntervention(:clean_frogs, 10, 0.5; mortality = 0.3)], s in 1:10
            x = solve(DrawertIBMProblem(p; intervention = iv, evaluate = [1.0, 120.0]); rng = Xoshiro(s))
            @test x.survivors == x.uninfected + x.infected
            @test x.survivors == p.N - count(x.culled) - count(x.treatment_death) - count(!isnan, x.death_time)
            x.outcome === :fungus_extinct && @test x.infected == 0
            x.t_end <= 120 && @test (x.evaluation[2].uninfected, x.evaluation[2].infected) == (x.uninfected, x.infected)
        end
        # a pulse on an evaluation day has been applied
        c = solve(DrawertIBMProblem(p; intervention = DrawertIntervention(:clean_frogs, 10, 1), tfinal = 10.0, evaluate = [10.0]); rng = Xoshiro(1))
        @test c.evaluation[1].infected == 0
        @test_throws ArgumentError DrawertIBMProblem(p; evaluate = [2000.0])
        @test_throws ArgumentError DrawertIBMProblem(p; infection_threshold = -1)
    end

    @testset "draws" begin
        ds = drawert2017_draws()
        @test length(ds) == 2132 && isconcretetype(eltype(ds))
        @test [count(d -> d.group === g, ds) for g in (:A, :B, :C, :D, :E)] == [1386, 109, 221, 396, 20]
        @test DrawertParameters(25) == DrawertParameters(ds[25]) ==
              DrawertParameters(load_drawert2017()[25]) == DrawertParameters(:drawert2017; point = 25)
        @test DrawertParameters(ds[25]; S0 = 5).S0 == 5
        @test load_growth_rate(DrawertParameters(25)) ≈ ds[25].load_growth_rate rtol = 1e-4
        for n in (5, 100, 500, 2132)
            s = C._drawert_draw_sample(n)
            @test length(s) == n && allunique(d.draw_id for d in s)
            @test Set(d.group for d in s) == Set((:A, :B, :C, :D, :E))
            @test sum(d.weight for d in s) ≈ 2132
        end
        @test C._drawert_draw_sample(100) == C._drawert_draw_sample(100)
        # SplitMix64 reference value (first output from state 0), and seeds do not collide
        @test C._splitmix64(UInt64(0)) == 0xe220a8397b1dcdaf
        @test allunique(drawert_run_seed(d, r) for d in 1:2132 for r in 1:30)
    end

    @testset "common random numbers" begin
        prob(iv) = DrawertIBMProblem(p; intervention = iv, saveat = [7.0, 13.9], save_loads = true)
        a = solve(prob(nothing); rng = Xoshiro(drawert_run_seed(1, 1)))
        b = solve(prob(DrawertIntervention(:cull_frogs, 14, 0.75)); rng = Xoshiro(drawert_run_seed(1, 1)))
        @test a.loads == b.loads && a.zoospores == b.zoospores
        # saving the state does not change the run
        @test solve(DrawertIBMProblem(p); rng = Xoshiro(9)).events ==
              solve(DrawertIBMProblem(p; saveat = [5.0]); rng = Xoshiro(9)).events
    end

    @testset "survey data from a run (new)" begin
        obs = SurveyObservation(; ρ = 0.5, κ = 0.01, σe = 0.5)
        days = [7.0, 14.0, 21.0]
        slow(ν0) = DrawertParameters(; N = 40, γ = 0.01, η = 8.0, ν0, f = 0.25, σ0 = 0.3, μ = 0.5,
                                     Smax = 10000, S0 = 100)
        p1, p2 = slow(0.1), slow(0.06)
        run_(q, s) = solve(DrawertIBMProblem(q; tfinal = 21.0, saveat = days, save_loads = true); rng = Xoshiro(s))
        pond = SurveyPond(; N = p.N, days, swabs = 10, pool = 2, edna = 3)
        r = run_(p1, 1)
        x = simulate_survey(Xoshiro(2), p1, r, obs, pond)
        @test length(x) == 3 * 13 && all(s -> s.copies >= 0, x)
        # swab mean ρ × pool × mean load; eDNA median κ Z
        big = [simulate_survey(Xoshiro(s), p1, r, obs, SurveyPond(; N = p.N, days = [21.0], swabs = 20, pool = 2, edna = 20))
               for s in 1:200]
        live = filter(>=(0), r.loads[:, 3])
        sw = [s.copies for y in big for s in y if s.kind === :swab]
        @test mean(sw) ≈ 2 * 0.5 * mean(live) atol = 4 * sqrt(var(sw) / length(sw))
        ed = [log(s.copies) for y in big for s in y if s.kind === :edna]
        r.zoospores[3] > 1000 && @test mean(ed) ≈ log(0.01 * r.zoospores[3]) atol = 4 * 0.5 / sqrt(length(ed))
        # data are more likely, on average, under the parameters that made them
        for (truth, other) in ((p1, p2), (p2, p1))
            Δ = [(y = simulate_survey(Xoshiro(s), truth, run_(truth, 100 + s), obs, pond);
                  survey_loglikelihood(y, truth, obs; nruns = 20) -
                  max(survey_loglikelihood(y, other, obs; nruns = 20), -1e6)) for s in 1:10]
            @test mean(Δ) > 0 && count(>(0), Δ) >= 8
        end
        # log mean of the conditional likelihoods
        ls = [survey_loglikelihood(x, p1, run_(p1, s), obs) for s in 1:5]
        @test survey_loglikelihood(x, p1, obs; seeds = 1:5) ≈ log(mean(exp.(ls .- maximum(ls)))) + maximum(ls)
        @test survey_loglikelihood(SurveySample[], p1, obs) == 0
        # no Bd: no copies, and copies are impossible
        q0 = DrawertParameters(; N = 10, γ = 0.0, η = 1.0, ν0 = 0.1, f = 0.5, σ0 = 5.0, μ = 5.0,
                               Smax = 10000, S0 = 1)
        y0 = simulate_survey(Xoshiro(1), q0, run_(q0, 1), obs, SurveyPond(; N = 10, days, swabs = 2, edna = 1))
        @test all(s -> s.copies == 0, y0) && survey_loglikelihood(y0, q0, run_(q0, 1), obs) == 0
        @test survey_loglikelihood([SurveySample(1, 21.0, :swab, 1, 5.0)], q0, run_(q0, 1), obs) == -Inf
        @test_throws ArgumentError simulate_survey(Xoshiro(1), p1, solve(DrawertIBMProblem(p1; saveat = days, tfinal = 21.0); rng = Xoshiro(1)), obs, pond)
        @test_throws ArgumentError simulate_survey(Xoshiro(1), p1, r, obs, SurveyPond(; N = 7, days, swabs = 1))
    end

    @testset "stored table" begin
        dir = dirname(drawert2017_path())
        des = TOML.parsefile(joinpath(dir, "decision_design.toml"))
        t = load_drawert_decisions()
        @test length(t) == des["rows"] == des["draws"] * length(des["actions"]) *
                           des["replicates"] * length(des["variants"])
        @test allunique((x.variant, x.draw_id, x.action, x.replicate) for x in t)
        @test Set(x.action for x in t) == Set(des["actions"]) && Set(x.group for x in t) == Set(string.('A':'E'))
        @test all(x -> x.status == "ok", t)
        @test all(x -> x.seed == drawert_run_seed(x.draw_id, x.replicate), t)
        @test all(x -> x.treatment_deaths == 0, t)
        @test any(x -> x.action == "none", t)
        @test all(x -> x.survivors == x.uninfected + x.infected <= x.initial_frogs, t)
        @test all(x -> x.surviving_fraction == x.survivors / x.initial_frogs, t)
        @test all(x -> x.survivors <= x.uninfected_d365 + x.infected_d365 <=
                       x.uninfected_d120 + x.infected_d120 <= x.initial_frogs, t)
        @test all(x -> x.outcome != "fungus_extinct" || x.infected == 0, t)
        # culled frogs are not survivors
        @test all(x -> x.kind != "cull_frogs" || x.survivors <= x.initial_frogs - x.handled, t)
        @test all(x -> x.action == "none" ? x.handled == 0 && x.frogs_at_intervention === missing :
                       x.kind == "clean_environment" || !x.applied ? x.handled == 0 :
                       0 <= x.handled <= x.frogs_at_intervention, t)
        frogs = [x for x in t if x.applied === true && x.kind != "clean_environment"]
        n = sum(x.frogs_at_intervention for x in frogs)
        @test sum(x.handled for x in frogs) / n ≈ 0.75 atol = 4 * sqrt(0.75 * 0.25 / n)
        @test all(x -> x.seconds !== missing && x.seconds >= 0, t)
        # re-run the cheapest cells of the table (the same Julia version gives identical rows)
        if des["julia"] == string(VERSION)
            draws = Dict(d.draw_id => d for d in drawert2017_draws())
            acts = Dict(C._drawert_action_label(a) => a for a in [nothing; drawert_interventions()])
            for x in sort(t; by = x -> x.seconds)[1:200]
                y, _ = C._drawert_decision_run(draws[x.draw_id], Symbol(x.variant), acts[x.action],
                                               x.replicate; S_star = des["S_star"])
                @test all(isequal(getfield(x, c), getfield(y, c)) for c in C._DRAWERT_DECISION_COLUMNS)
            end
        end
    end
end
