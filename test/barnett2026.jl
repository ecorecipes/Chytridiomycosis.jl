# Barnett et al. 2026 (bioRxiv PREPRINT, doi:10.64898/2026.05.15.725113): agent-based model of
# Bd with an imperfect vaccine (Bd-metabolite prophylaxis), and the field-trial GLMMs. The
# oracles are:
# - the authors' NetLogo code and BehaviorSpace experiments
#   (github.com/kmbarnett100/protection_to_amplification, commit 29f4af6): parameter values,
#   vaccine constants, pond geometry, and the statistics of their 250 control runs;
# - the stored comparison with all 29,325 author runs (data/barnett2026/reproduction.toml,
#   written by scripts/barnett2026_abm.jl);
# - the preprint text (Methods) and Supplementary Table S1;
# - closed forms for the Bd-free model: survival to day 90 under the natural mortalities and
#   the metamorphosis schedule;
# - with the authors' field data (not bundled; skipped without it): the marginal means, standard
#   errors and intervals in their "analysis files", and a glmmTMB 1.1 run of their R script
#   (log-likelihoods, interaction terms), on 2026-11-23.
using Chytridiomycosis
using Test
using TOML
using Random: Xoshiro

@testset "Barnett 2026 prophylaxis ABM (preprint)" begin
    C = Chytridiomycosis
    q = parameter_set(:barnett2026)
    p = BarnettParameters(:barnett2026)

    @testset "parameters (Methods, Table S1, NetLogo)" begin
        @test all(r -> r.verified, parameter_records(:barnett2026))
        @test (p.tad_mort, p.meta_mort, p.t_movement, p.m_land) == (0.06, 0.02, 0.25, 0.1)
        @test (p.meta_start, p.meta_prob, p.meta_end, p.last_day) == (55, 0.11, 74, 90)
        @test (p.tadpoles_per_patch, p.infected_per_patch, p.init_load) == (500, 1, 100)
        @test p.infected_per_patch / p.tadpoles_per_patch == q.init_prevalence == 0.002
        @test (p.exposure, p.est, p.z_mort) == (0.25, 0.25, 2.0)
        @test (p.clearance, p.shedding) == (0.2, 17.8)
        @test (p.self_reinfection, p.shed_own, p.shed_neighbour) == (0.1, 0.4, 0.58)
        @test (p.smax, p.s_k) == (562, 10000)
        # Table S1 derives smax as 10,000 zoospores per day / 17.8 per zoosporangium
        @test round(10000 / p.shedding) == p.smax
        @test (p.world_size, p.pond_mean, p.pond_sd, p.natural_imm) == (21, 135, 18, 0)
        @test (q.n_reps, q.n_control, q.coverage_two_way) == (25, 250, 0.75)
        @test C.time_unit(BarnettParameters) === :day
        @test_throws ArgumentError BarnettABMProblem(C.with(p; shed_own = 1.5))
        @test_throws ArgumentError BarnettABMProblem(C.with(p; est = 1.5))
        @test_throws ArgumentError BarnettABMProblem(C.with(p; infected_per_patch = 501))
    end

    @testset "vaccine constants (BehaviorSpace)" begin
        up = [0.095, 0.182, 0.262, 0.336, 0.405, 0.47, 0.531, 0.588, 0.642, 0.693]
        down = -[0.11, 0.223, 0.357, 0.511, 0.693, 0.916, 1.204, 1.609, 2.303, 16.118]
        es = 0.1:0.1:1.0
        @test [barnett_vaccine(:infection, e; coverage = 0.5).c_est for e in es] == down
        @test [barnett_vaccine(:transmission, e; coverage = 0.5).c_shedding for e in es] == down
        @test [barnett_vaccine(:growth, e; coverage = 0.5).c_clear for e in es] == up
        @test [barnett_vaccine(:tolerance, e; coverage = 0.5).c_smax for e in es] == up
        @test [barnett_vaccine(:infection, e; coverage = 0.5, adverse = true).c_est
               for e in es] == up
        @test [barnett_vaccine(:tolerance, e; coverage = 0.5, adverse = true).c_smax
               for e in es] == down
        w = barnett_vaccine(:tolerance => 0.8, :growth => 0.2; coverage = q.coverage_two_way)
        @test (w.c_smax, w.c_clear, w.c_est, w.coverage) == (0.588, 0.182, 0.0, 0.75)
        @test barnett_vaccine(:infection, 0.0; coverage = 1).c_est == 0
        @test_throws ArgumentError barnett_vaccine(:resistance, 0.3; coverage = 1)
        @test_throws ArgumentError barnett_vaccine(:infection, 1.5; coverage = 1)
        @test_throws ArgumentError BarnettVaccine(; coverage = 1.2)
        @test_throws ArgumentError BarnettVaccine(; c_est = -Inf)
        @test C._nlround(2.5) == 3 && C._nlround(0.49) == 0     # NetLogo rounds halves up
    end

    @testset "pond" begin
        for s in 1:50
            nbrs = C._barnett_pond(p, Xoshiro(s))
            @test 25 <= length(nbrs) <= 60
            @test all(!isempty, nbrs)
            @test all(i in nbrs[j] for (i, l) in enumerate(nbrs) for j in l)   # symmetric
        end
    end

    @testset "Bd-free cohort: demography" begin
        # With no infection, survivors are Binomial(N, s) with s from the schedule: on tick t
        # (age t + 1) tadpoles die, then metamorphs die, then tadpoles older than 55 days
        # metamorphose with probability 0.11, all of them at age 75.
        p0 = C.with(p; infected_per_patch = 0)
        ptad, pmet = 1.0, 0.0
        for age in 1:p0.last_day
            ptad *= 1 - p0.tad_mort
            pmet *= 1 - p0.meta_mort
            if age > p0.meta_start
                m = age > p0.meta_end ? ptad : p0.meta_prob * ptad
                ptad -= m; pmet += m
            end
        end
        s = ptad + pmet
        rs = [solve(BarnettABMProblem(p0; save_daily = true); rng = Xoshiro(k)) for k in 1:10]
        N = sum(r.n_perimeter for r in rs) * p0.tadpoles_per_patch
        @test abs(sum(r.survivors for r in rs) - N * s) < 4 * sqrt(N * s * (1 - s))
        r = rs[1]
        @test r.tadpoles[end] == 0 && r.metamorphs[end] == r.survivors
        @test all(==(0), r.infected) && all(==(0), r.zoospores)
        @test r.prevalence == 0 && isnan(r.median_intensity)
        @test r.bd_death_fraction == 0 && r.zoospore_density == 0
    end

    @testset "one run" begin
        r = solve(BarnettABMProblem(p; save_daily = true); rng = Xoshiro(7))
        @test r.infected[1] >= 1 && r.zoospore_density > 0
        # about 1 run in 10 ends with no metamorph, whatever the seed: undefined outputs are NaN
        @test isequal(r.unvaccinated_fraction, r.survivors == 0 ? NaN : 1.0)
        @test r.realized_coverage == 0
        @test isnan(r.vaccinated_prevalence)
        @test isequal(r.unvaccinated_median_intensity, r.median_intensity)
        r2 = solve(BarnettABMProblem(p; save_daily = true); rng = Xoshiro(7))
        @test r2.survivors == r.survivors && r2.infected == r.infected
        # vaccination skips the infected tadpoles (0.2%); tolerance alone changes no metamorph
        # before day 55
        v = solve(BarnettABMProblem(p, barnett_vaccine(:tolerance, 1.0; coverage = 0.6);
                                    save_daily = true); rng = Xoshiro(8))
        @test 0.59 < v.realized_coverage <= 0.6
        @test all(==(0), v.metamorphs[1:55])
        # no establishment and no infected hosts' reinfection: Bd dies out
        z = solve(BarnettABMProblem(C.with(p; est = 0.0); save_daily = true); rng = Xoshiro(9))
        @test z.infected[60] == 0 && z.prevalence == 0
    end

    @testset "comparison with the authors' runs (stored)" begin
        d = TOML.parsefile(joinpath(pkgdir(Chytridiomycosis), "data", "barnett2026",
                                    "reproduction.toml"))
        cells = [r for (_, s) in d if s isa Dict for (_, r) in s if r isa Dict]
        @test sum(s["author_runs"] for (_, s) in d if s isa Dict) == 29_325
        @test length(cells) == d["comparisons"] >= 80
        # outside the one file that differs (increased shedding), chance-level deviations
        rest = [r for (k, s) in d if s isa Dict && k != "adverse_transmission"
                for (_, r) in s if r isa Dict]
        @test maximum(abs(r["z"]) for r in rest) < 3.5
        @test count(abs(r["z"]) < 2 for r in rest) / length(rest) > 0.85
        @test all(r["ks"] < 0.1 for r in cells)
    end

    @testset "summaries" begin
        c = barnett_summary(p; nreps = 60, rng = Xoshiro(11))
        # the authors' 250 control runs: mean survivors 128.8 (SD 85.3), die-off 9.2%
        @test abs(c.mean.survivors - 128.84) < 4 * 85.3 * sqrt(1 / 60 + 1 / 250)
        @test c.median.prevalence == 1
        @test 0 <= c.dieoff < 0.3
        e = barnett_relative_change(p, barnett_vaccine(:infection, 0.8; coverage = 0.8);
                                    control = c, nreps = 15, rng = Xoshiro(13))
        @test e.intensity < -80 && e.survivors > 0 && e.unvaccinated_intensity < -50
        @test e.scenario.nreps == 15 && e.control === c
    end
end

@testset "Barnett 2026 field trial (preprint)" begin
    C = Chytridiomycosis
    found = C._barnett2026_find() !== nothing
    found || @info "Barnett et al. 2026 data not found; skipping the data-based tests"
    if found
        d = load_barnett2026()
        @test length(d.field.load) == 1672 && length(unique(d.field.pond)) == 12
        @test length(d.challenge.load) == 487 && length(unique(d.challenge.pond)) == 11
        @test_throws ArgumentError BarnettBACIModel(:x; data = d)
        # kind => (log-likelihood, interaction estimate, SE, link-scale means and SEs)
        oracle = (
            field_intensity = (-3800.28526, -1.8271, 0.5701,
                (5.269169, 7.393138, 6.502555, 6.764130), (1.498902, 0.815468, 1.484365, 0.838389)),
            field_prevalence = (-794.0360558, -0.1133, 0.2664,
                (-2.095372, -0.578475, -1.308265, -0.281467), nothing),
            challenge_intensity = (-3241.96101, -0.2035, 0.5941,
                (9.936640, 9.027983, 10.541149, 9.816443), (0.480209, 0.496288, 0.474435, 0.512204)),
            challenge_prevalence = (-241.3945273, -0.2783, 0.6246,
                (2.802612, -0.092241, 2.517240, -0.622401), nothing))
        for (kind, o) in pairs(oracle)
            m = BarnettBACIModel(kind; data = d)
            @test length(parameter_names(m)) == C.LogDensityProblems.dimension(m)
            r = C._barnett_baci(kind)
            @test r.loglik ≈ o[1] atol = 1e-4
            @test r.loglik ≈ C.LogDensityProblems.logdensity(m, r.x) atol = 1e-6
            @test r.interaction.estimate ≈ o[2] atol = 5e-4
            @test r.interaction.se ≈ o[3] atol = 5e-4
            links = map(x -> x.link, values(r.means))
            @test all(isapprox.(links, o[4]; atol = 1e-4))
            o[5] === nothing || @test all(isapprox.(map(x -> x.se, values(r.means)), o[5]; atol = 2e-4))
        end
        @test C._barnett_baci(:field_intensity).interaction.p ≈ 0.001352 atol = 5e-6
        # probabilities and interval of the authors' prev_field_marginal_predictions.csv
        b = C._barnett_baci(:field_prevalence).means.sham_before
        @test b.response ≈ 0.359242 atol = 1e-4
        @test b.response_lower ≈ 0.062642 atol = 1e-4
        @test b.response_upper ≈ 0.824669 atol = 1e-4
    end
end
