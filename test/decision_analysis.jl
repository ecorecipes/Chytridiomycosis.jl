# The DecisionAnalysis.jl adapter (ext/ChytridiomycosisDecisionAnalysisExt.jl) for the Drawert
# et al. (2017) problems. The oracles are this package's own solver and stored table:
# - a run through the verbs is the run of `solve` with the same generator;
# - exact identities: counts add up, a pulse that has not happened yet changes nothing (common
#   random numbers), a deterministic model has no spread between replicates;
# - the stored table read through the adapter has the cell means of the CSV.
using Chytridiomycosis
using DecisionAnalysis
using Test
using Random: Xoshiro
using Distributions: mean

@testset "DecisionAnalysis adapter (Drawert 2017)" begin
    DA = DecisionAnalysis
    p = DrawertParameters(; N = 30, γ = 0.01, η = 8.0, ν0 = 0.2, f = 0.25, σ0 = 0.3, μ = 0.5,
                          Smax = 2000, S0 = 100)
    prob = DrawertIBMProblem(p; tfinal = 200.0, evaluate = [60.0, 120.0])
    clean = DrawertIntervention(:clean_frogs, 14, 0.75)

    @testset "parameters" begin
        q = DA.with_parameters(prob, (; draw_id = 7, group = "B", N = 12, γ = 0.02, weight = 3.0))
        @test q.p.N == 12 && q.p.γ == 0.02
        @test q.p.η == p.η && q.p.Smax == p.Smax && q.p.S0 == p.S0
        @test q.tfinal == prob.tfinal && q.evaluate == prob.evaluate
        @test prob.p.N == 30                                        # not mutated
        @test DA.with_parameters(prob, (; N = 12.0)).p.N === 12     # a column of floats
        @test_throws ArgumentError DA.with_parameters(prob, (; draw_id = 7, weight = 1.0))
        @test_throws InexactError DA.with_parameters(prob, (; N = 12.5))
        @test_throws ArgumentError DA.with_parameters(prob, (; ν0 = 1.5))
        row = first(load_drawert_decision_draws())
        @test DA.with_parameters(prob, row).p.N == row.N
        ode = DA.with_parameters(DrawertODEProblem(p; tfinal = 50.0), (; μ = 0.7))
        @test ode.p.μ == 0.7 && ode.tfinal == 50.0
    end

    @testset "actions" begin
        @test isempty(DA.apply_action(prob, nothing).interventions)
        @test isempty(DA.apply_action(prob, DrawertIntervention[]).interventions)
        @test DA.apply_action(prob, clean).interventions == [clean]
        nt = DA.apply_action(prob, (; kind = :clean_frogs, day = 14, fraction = 0.75))
        @test nt.interventions == [clean]
        m = only(DA.apply_action(prob, (; kind = :clean_frogs, day = 14, fraction = 0.75,
                                        mortality = 0.1)).interventions)
        @test m.mortality == 0.1
        two = DA.apply_action(prob, [DrawertIntervention(:cull_frogs, 28, 0.25), clean])
        @test [iv.day for iv in two.interventions] == [14.0, 28.0]
        # an action replaces the problem's own intervention
        @test isempty(DA.apply_action(DA.apply_action(prob, clean), nothing).interventions)
        @test_throws ArgumentError DA.apply_action(prob, :clean)
        @test_throws ArgumentError DA.apply_action(prob, (; kind = :bathe, day = 14, fraction = 0.5))
        ode = DrawertODEProblem(p; tfinal = 50.0)
        @test DA.apply_action(ode, clean).intervention == clean
        @test DA.apply_action(ode, nothing).intervention === nothing
        @test_throws ArgumentError DA.apply_action(ode, [clean, clean])
    end

    @testset "one run is the solver's run" begin
        acted = DA.apply_action(prob, clean)
        direct = solve(DrawertIBMProblem(p; intervention = clean, tfinal = 200.0,
                                         evaluate = [60.0, 120.0]); rng = Xoshiro(11))
        run = DA.simulate(acted; rng = Xoshiro(11))
        @test run.survivors == direct.survivors && run.events == direct.events
        o = DA.outcome(acted, run)
        @test keys(o) == (:initial_frogs, :survivors, :uninfected, :infected, :handled,
                          :treatment_deaths, :frogs_extinct, :fungus_extinct, :persistence,
                          :uninfected_d60, :infected_d60, :uninfected_d120, :infected_d120)
        @test o.initial_frogs == 30 && o.survivors == direct.survivors
        @test o.uninfected + o.infected == o.survivors
        @test o.frogs_extinct + o.fungus_extinct + o.persistence == 1
        @test o.handled == direct.handled
        @test o.uninfected_d60 == direct.evaluation[1].uninfected
        @test o.infected_d120 == direct.evaluation[2].infected
    end

    actions = (; none = nothing, clean = clean, cull = DrawertIntervention(:cull_frogs, 14, 0.5),
               late = DrawertIntervention(:clean_frogs, 199.999, 0.75))
    draws = [(; N = 20, ν0 = 0.2), (; N = 30, ν0 = 0.15), (; N = 40, ν0 = 0.25)]

    @testset "outcomes over draws × actions × replicates" begin
        t = outcomes(prob, actions, draws; replicates = 6, seed = 2017, on_error = :throw)
        @test run_counts(t) == (runs = 72, ok = 72, failed = 0, absent = 0)
        @test action_labels(t) == [:none, :clean, :cull, :late]
        @test parameter_draws(t).N == [20, 30, 40]
        again = outcomes(prob, actions, draws; replicates = 6, seed = 2017)
        @test again.values == t.values                               # reproducible
        @test outcomes(prob, actions, draws; replicates = 6, seed = 2018).values != t.values
        k(name) = findfirst(==(name), outcome_names(t))
        v(action, name) = t.values[findfirst(==(action), action_labels(t)), :, :, k(name)]
        @test all(v(:none, :initial_frogs) .== [20, 30, 40])         # [draw, replicate]
        @test v(:none, :uninfected) .+ v(:none, :infected) == v(:none, :survivors)
        @test all(v(:none, :handled) .== 0)
        @test all(v(:cull, :survivors) .<= v(:cull, :initial_frogs) .- v(:cull, :handled))
        # Common random numbers: on days 60 and 120 the late pulse has not happened, so the run
        # is the baseline's
        @test v(:late, :uninfected_d60) == v(:none, :uninfected_d60)
        @test v(:late, :infected_d120) == v(:none, :infected_d120)
        @test any(v(:clean, :handled) .> 0)
        payoff = payoff_table(t, IndividualsUtility(:survivors))
        @test size(payoff.Q) == (4, 3)
        @test all(payoff.Q[1, :] .<= [20, 30, 40])
    end

    @testset "the deterministic model" begin
        ode = DrawertODEProblem(p; tfinal = 200.0)
        t = outcomes(ode, (; none = nothing, clean = clean), draws; replicates = 2, seed = 1,
                     on_error = :throw)
        @test outcome_names(t) == [:initial_frogs, :survivors, :frogs_extinct, :fungus_extinct,
                                   :persistence]
        @test t.values[:, :, 1, :] == t.values[:, :, 2, :]           # no spread
        q = DrawertParameters(; N = 30, γ = 0.01, η = 8.0, ν0 = 0.15, f = 0.25, σ0 = 0.3, μ = 0.5,
                              Smax = 2000, S0 = 100)
        direct = solve(DrawertODEProblem(q; intervention = clean, tfinal = 200.0))
        @test t.values[2, 2, 1, 2] == direct.survivors
    end

    @testset "failed runs are counted" begin
        bad = [(; N = 20, ν0 = 0.2), (; N = 20, ν0 = 1.5)]
        t = outcomes(prob, (; none = nothing), bad; replicates = 2, seed = 3)
        @test run_counts(t).failed == 2
        @test_throws ArgumentError outcomes(prob, (; none = nothing), bad; replicates = 2, seed = 3,
                                            on_error = :throw)
    end

    @testset "the stored table" begin
        t = drawert_outcome_table()
        runs = load_drawert_decisions()
        stored = load_drawert_decision_draws()
        @test run_counts(t) == (runs = length(runs), ok = length(runs), failed = 0, absent = 0)
        @test length(state_labels(t)) == length(stored) == 100
        @test :none in action_labels(t)
        @test draw_weights(t) ≈ [d.weight for d in stored] ./ sum(d.weight for d in stored)
        @test sum(d.weight for d in stored) ≈ 2132
        @test parameter_draws(t).N == [d.N for d in stored]
        payoff = payoff_table(t, IndividualsUtility(:survivors))
        a = findfirst(==(:none), action_labels(t))
        d = 17
        cell = [r.survivors for r in runs if r.action == "none" && r.draw_id == stored[d].draw_id]
        @test Float64(payoff.Q[a, d]) ≈ mean(cell)
        small = drawert_outcome_table(; outcomes = (:survivors,))
        @test outcome_names(small) == [:survivors]
        @test_throws ArgumentError drawert_outcome_table(; variant = "nu")
        @test t.replicates[1:2] == [:r1, :r2] && length(t.replicates) == 30
    end

    @testset "the stored seeds: stored runs reproduced, new actions joined" begin
        runs = load_drawert_decisions()
        draws = load_drawert_decision_draws()
        @test drawert_run_seed(draws[3], 2) == drawert_run_seed(draws[3].draw_id, 2)
        # The smallest ponds, to keep this quick; two replicates; the stored `none` and one
        # stored clean, run again through the adapter with the stored seeds.
        few = sort(draws; by = d -> d.N)[1:4]
        ids = [d.draw_id for d in few]
        stored_clean = DrawertIntervention(:clean_frogs, 28, 0.75)
        two_cleans = [DrawertIntervention(:clean_frogs, 14, 0.75), stored_clean]
        belief = SampleBelief((; (k => [getfield(d, k) for d in few] for k in keys(first(few)))...);
                              weights = :weight, labels = ids)
        base = DrawertIBMProblem(DrawertParameters(first(few)); evaluate = (120.0, 365.0))
        again = outcomes(base, ["none" => nothing, "clean_frogs_d28_f75" => stored_clean,
                                "two_cleans" => two_cleans], belief;
                         replicates = 2, run_seed = drawert_run_seed, on_error = :throw)
        names = (:survivors, :uninfected, :infected, :uninfected_d120, :infected_d120,
                 :uninfected_d365, :infected_d365, :initial_frogs, :handled, :treatment_deaths)
        for (a, label) in enumerate(("none", "clean_frogs_d28_f75")), (d, id) in enumerate(ids),
            r in 1:2
            row = only(x for x in runs if x.variant == "baseline" && x.action == label &&
                       x.draw_id == id && x.replicate == r)
            for name in names
                k = findfirst(==(name), outcome_names(again))
                @test again.values[a, d, r, k] == getfield(row, name)
            end
        end
        # Joined to the stored rows of the same draws and replicates.
        part = OutcomeTable([(; x..., replicate = string("r", x.replicate)) for x in runs
                             if x.variant == "baseline" && x.draw_id in ids && x.replicate <= 2];
                            draw = :draw_id, action = :action, replicate = :replicate,
                            outcomes = collect(names), status = :status, parameters = few,
                            key = :draw_id, weights = :weight)
        new = outcomes(base, ["two_cleans" => two_cleans], belief; replicates = 2,
                       run_seed = drawert_run_seed, on_error = :throw)
        joined = combine_actions(part, new)
        @test action_labels(joined) == vcat(action_labels(part), :two_cleans)
        @test outcome_names(joined) == collect(names)
        @test joined.values[end, :, :, :] ==
              again.values[3, :, :, [findfirst(==(n), outcome_names(again)) for n in names]]
        # Common random numbers with the stored runs: until day 14 nothing differs, so a pond
        # whose outcome is settled before then has the stored outcome under the new action.
        none = findfirst(==(:none), action_labels(joined))
        settled = [(d, r) for d in 1:4, r in 1:2
                   if only(x.first_loss for x in runs if x.variant == "baseline" &&
                           x.action == "none" && x.draw_id == ids[d] && x.replicate == r) < 14 &&
                      joined.values[none, d, r, 1] == 0]
        for (d, r) in settled
            @test joined.values[end, d, r, 1] == joined.values[none, d, r, 1]
        end
        @test_throws ArgumentError drawert_outcome_table(; variant = "nu",
                                                         actions = ["x" => two_cleans])
    end

    @testset "DrawertSurveyStudy (AbstractStudyDesign)" begin
        obs = SurveyObservation(; ρ = 1.0, κ = 1.0, σe = 0.5)
        pond = SurveyPond(; N = 30, days = [7.0], swabs = 2, pool = 5, edna = 1)
        study = drawert_survey_study(obs, pond; nruns = 3, template = p)
        @test study isa DecisionAnalysis.AbstractStudyDesign
        small_draws = [(; N = 20, ν0 = 0.2), (; N = 30, ν0 = 0.15), (; N = 40, ν0 = 0.25),
                       (; N = 25, ν0 = 0.3)]
        belief = SampleBelief((; N = [d.N for d in small_draws],
                               ν0 = [d.ν0 for d in small_draws]))
        dat = simulate_study(Xoshiro(7), study, small_draws[2])
        @test dat isa Vector{SurveySample}
        @test !isempty(dat)
        ll = study_loglikelihood(study, dat, small_draws[2])
        @test ll isa Float64 && (isfinite(ll) || ll == -Inf)
        ss = SummaryStatistics(belief, study; seed = 9)
        @test ss isa SummaryStatistics
        is_study = ImportanceSampling(belief, study; seed = 9, datasets = 2)
        @test is_study isa ImportanceSampling
        @test size(is_study.loglikelihood) == (2, 4)
        @test_throws ArgumentError study_posterior(Xoshiro(1), study, dat, (;))
    end

    @testset "Canessa 2018 adapters (IPM, SI, canessa_outcome_table)" begin
        d_all = load_canessa2018()[1:4]
        mesh = (-10.0, 15.0, 40)
        # CanessaIPMProblem verbs
        ipm = CanessaIPMProblem(CanessaParameters(d_all[1]); nsteps = 10, mesh)
        ipm_d2 = DA.with_parameters(ipm, d_all[2])
        @test ipm_d2.p.S0 == d_all[2].S0
        ipm_nt = DA.with_parameters(ipm, (; beta0 = 0.01, S0 = 400.0))
        @test ipm_nt.p.β0 == 0.01 && ipm_nt.p.S0 == 400.0
        @test_throws ArgumentError DA.with_parameters(ipm, (; unused = 1.0))
        ipm_f = DA.apply_action(ipm, :f)
        @test ipm_f.p.S0 ≈ 0.5 * ipm.p.S0
        @test_throws ArgumentError DA.apply_action(ipm, :d) # growth_factor != 1
        r_ipm = DA.simulate(ipm_f; rng = Xoshiro(1))
        o_ipm = DA.outcome(ipm_f, r_ipm)
        @test o_ipm.final_susceptibles == r_ipm.S[end]
        @test o_ipm.R0 ≈ basic_reproduction_number(ipm_f.p; mesh)
        t_ipm = outcomes(ipm, [:a => :a, :f => :f, :g => :g],
                         [(; β0 = 0.005, S0 = 400.0), (; β0 = 0.009, S0 = 600.0)];
                         replicates = 1, seed = 1, on_error = :throw)
        @test action_labels(t_ipm) == [:a, :f, :g]
        @test run_counts(t_ipm).ok == 6

        # CanessaSIProblem verbs
        si = CanessaSIProblem(CanessaSIParameters(:canessa2018); nsteps = 20)
        si_draws = [(; β = 0.005, S0 = 400.0), (; β = 0.009, S0 = 600.0)]
        si_acts = [:none => nothing, :thin50 => (; mS0 = 0.5), :trans50 => (; mβ = 0.5)]
        t_si = outcomes(si, si_acts, si_draws; replicates = 1, seed = 1, on_error = :throw)
        @test action_labels(t_si) == [:none, :thin50, :trans50]
        @test run_counts(t_si).ok == 6
        @test_throws ArgumentError DA.apply_action(si, (; bogus = 0.1))

        # canessa_outcome_table includes growth_factor != 1 actions (e.g. :d)
        tbl = canessa_outcome_table(; actions = [:a, :d, :f], draws = d_all, nsteps = 10, mesh)
        @test tbl isa OutcomeTable
        @test action_labels(tbl) == [:a, :d, :f]
        @test length(state_labels(tbl)) == 4
        ens_d = canessa_ensemble(canessa_action(:d); draws = d_all, nsteps = 10, mesh)
        k_R0 = findfirst(==(:R0), outcome_names(tbl))
        @test [tbl.values[2, i, 1, k_R0] for i in 1:4] ≈ [e.R0 for e in ens_d]
    end

    @testset "Cook 2022 adapter (CookProjectionProblem)" begin
        cp = CookParameters(:cook2022)
        cprob = CookProjectionProblem(cp; years = 5)
        c_acts = [:none => nothing, :treat50 => (; treated_fraction = 0.5, treated_years = 1.0)]
        c_draws = [(; SU = 0.10, ST = 0.45), (; SU = 0.15, ST = 0.55)]
        ct = outcomes(cprob, c_acts, c_draws; replicates = 3, seed = 2022, on_error = :throw)
        @test action_labels(ct) == [:none, :treat50]
        @test run_counts(ct) == (runs = 12, ok = 12, failed = 0, absent = 0)
        @test :final_adults in outcome_names(ct) && :lambda in outcome_names(ct)
        # Oracle parity with direct solve
        acted = DA.apply_action(DA.with_parameters(cprob, c_draws[1]), c_acts[2][2])
        direct = solve(acted; nsims = 1, rng = Xoshiro(99))
        o = DA.outcome(acted, DA.simulate(acted; rng = Xoshiro(99)))
        @test o.final_adults == direct.states[end, 4, 1] && o.lambda == direct.λ[1]
        @test_throws ArgumentError DA.apply_action(cprob, (; bad = 1.0))
    end

    @testset "Heard 2024 adapter (HeardMetapopulationProblem)" begin
        hp = HeardParameters(:heard2024; species = :raniformis, epoch = :contemporary,
                             mets = 20, mets_cv = 0.3, mortality_cv = 0.3)
        hp_short = HeardParameters(hp.adult_survival, hp.subadult_survival, hp.carrying_capacity,
                                   hp.prop_dispersing, hp.mets, hp.mets_cv, hp.mortality_cv,
                                   hp.met_survival, hp.sex_ratio, hp.subadult_breeding_fraction,
                                   hp.initial_fraction, 5, hp.cell_size,
                                   hp.distance_decay, hp.max_distance)
        hprob = HeardMetapopulationProblem(hp_short, 1)
        h_acts = [:none => nothing, :boost => (; survival_multiplier = 1.25)]
        h_draws = [(; adult_survival = 0.30, mets = 18.0), (; adult_survival = 0.40, mets = 24.0)]
        ht = outcomes(hprob, h_acts, h_draws; replicates = 2, seed = 2024, on_error = :throw)
        @test action_labels(ht) == [:none, :boost]
        @test run_counts(ht).ok == 8
        @test :final_females in outcome_names(ht) && :persisted in outcome_names(ht)
        # Coupled subadult survival follows adult_survival
        hp_upd = DA.with_parameters(hprob, (; adult_survival = 0.55))
        @test hp_upd.p.adult_survival == 0.55 && hp_upd.p.subadult_survival == 0.55
        @test_throws ArgumentError DA.apply_action(hprob, (; unknown = 1))
    end

    @testset "Turner 2025 adapter (TurnerPVAProblem)" begin
        tprob = TurnerPVAProblem(:turner2025; month = :august, pattern = :random, replicates = 1)
        t_acts = [:aug_rand => (; month = :august, pattern = :random),
                  :oct_max1 => (; month = :october, pattern = :max1dry)]
        t_draws = [(; p_wet = 0.65), (; p_wet = 0.80)]
        tt = outcomes(tprob, t_acts, t_draws; replicates = 2, seed = 2025, on_error = :throw)
        @test action_labels(tt) == [:aug_rand, :oct_max1]
        @test run_counts(tt).ok == 8
        @test :final_females in outcome_names(tt) && :quasi_extinct_20 in outcome_names(tt)
        @test_throws ArgumentError DA.apply_action(tprob, :not_a_pattern)
        @test_throws ArgumentError DA.apply_action(tprob, (; month = :july))
    end
end
