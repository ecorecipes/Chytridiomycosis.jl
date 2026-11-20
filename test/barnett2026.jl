# Barnett et al. 2026 (bioRxiv PREPRINT, doi:10.64898/2026.05.15.725113): agent-based model of
# Bd with an imperfect vaccine (Bd-metabolite prophylaxis). The oracles are:
# - the preprint text (Methods) and Supplementary Table S1, for every parameter value and
#   process (no model code was released);
# - closed forms for the Bd-free model: survival to day 90 under the natural mortalities and
#   the metamorphosis schedule, and Binomial/Poisson bookkeeping identities;
# - the preprint's qualitative statements about its scenario figures (validation targets).
using Chytridiomycosis
using Test
using Random: Xoshiro

@testset "Barnett 2026 prophylaxis ABM (preprint)" begin
    C = Chytridiomycosis
    q = parameter_set(:barnett2026)
    p = BarnettParameters(:barnett2026)

    @testset "parameters (Methods, Table S1)" begin
        @test all(r -> r.verified, parameter_records(:barnett2026))
        @test (p.tad_mort, p.meta_mort, p.t_movement, p.m_land) == (0.06, 0.02, 0.25, 0.1)
        @test (p.meta_start, p.meta_prob, p.meta_end, p.last_day) == (55, 0.11, 74, 90)
        @test (p.init_prevalence, p.init_load) == (0.002, 100)
        @test (p.exposure, p.est, p.z_mort) == (0.25, 0.25, 2.0)
        @test (p.maturation_days, p.clearance, p.shedding) == (4, 0.2, 17.8)
        @test (p.shed_own, p.shed_neighbour, p.self_reinfection) == (0.4, 0.5, 0.1)
        @test (p.smax, p.s_k) == (562, 10000)
        # Table S1 derives smax as 10,000 zoospores per day / 17.8 per zoosporangium
        @test round(10000 / p.shedding) == p.smax
        @test (q.n_reps, q.n_control, q.coverage_two_way) == (25, 250, 0.75)
        @test C.time_unit(BarnettParameters) === :day
        @test_throws ArgumentError BarnettABMProblem(C.with(p; shed_own = 0.5))
        @test_throws ArgumentError BarnettABMProblem(C.with(p; est = 1.5))
    end

    @testset "vaccine scenarios" begin
        v = barnett_vaccine(:infection, 0.3; coverage = 0.5)
        @test (v.coverage, v.infection, v.transmission, v.growth, v.tolerance) ==
              (0.5, -0.3, 0.0, 0.0, 0.0)
        @test barnett_vaccine(:growth, 0.3; coverage = 1).growth == 0.3
        @test barnett_vaccine(:tolerance, 0.3; coverage = 1).tolerance == 0.3
        @test barnett_vaccine(:transmission, 0.3; coverage = 1).transmission == -0.3
        @test barnett_vaccine(:infection, 0.3; coverage = 1, adverse = true).infection == 0.3
        @test barnett_vaccine(:tolerance, 0.3; coverage = 1, adverse = true).tolerance == -0.3
        w = barnett_vaccine(:tolerance => 0.8, :growth => 0.2; coverage = q.coverage_two_way)
        @test (w.tolerance, w.growth, w.coverage) == (0.8, 0.2, 0.75)
        @test_throws ArgumentError barnett_vaccine(:resistance, 0.3; coverage = 1)
        @test_throws ArgumentError barnett_vaccine(:infection, 1.5; coverage = 1)
        @test_throws ArgumentError BarnettVaccine(; coverage = 1.2)
        # exponential modulation x0 (1 + e)^imm; a 100% decrease gives 0
        @test C._barnett_mult(0.0, 1.07) == 1
        @test C._barnett_mult(-1.0, 0.93) == 0
        @test C._barnett_mult(1.0, 1.0) ≈ 2
        @test C._barnett_mult(-0.5, 1.0) ≈ 0.5
    end

    @testset "Bd-free cohort: demography" begin
        # No infection: survivors at day 90 are Binomial(N, s) with s from the schedule.
        p0 = C.with(p; init_prevalence = 0.0, n_tadpoles = 20000)
        r = solve(BarnettABMProblem(p0; save_daily = true); rng = Xoshiro(1))
        # survival: tadpole days 1..54 then metamorphosis from day 55 (after that day's
        # deaths) with probability 0.11 per day, all by day 74; metamorphs die at 0.02.
        ptad, pmet = 1.0, 0.0
        for d in 1:p0.last_day
            ptad *= 1 - p0.tad_mort
            pmet *= 1 - p0.meta_mort
            if d >= p0.meta_start
                m = d >= p0.meta_end ? ptad : p0.meta_prob * ptad
                ptad -= m; pmet += m
            end
        end
        s = ptad + pmet
        N = p0.n_tadpoles
        @test abs(r.survivors - N * s) < 4 * sqrt(N * s * (1 - s))
        @test r.tadpoles[end] == 0
        @test r.metamorphs[end] == r.survivors
        @test all(==(0), r.infected) && all(==(0), r.zoospores)
        @test isnan(r.final_intensity) && r.final_prevalence == 0
    end

    @testset "one run" begin
        r = solve(BarnettABMProblem(p; save_daily = true); rng = Xoshiro(7))
        @test r.infected[1] >= 1                 # 0.2% of 1000 tadpoles start infected
        @test maximum(r.zoospores) > 0
        @test length(r.loads) == r.survivors == r.vaccinated_survivors + r.unvaccinated_survivors
        @test r.n_vaccinated == 0
        # metamorph survivors carry loads below the lethal threshold
        @test all(<(p.smax), r.loads)
        # reproducible with a fixed seed
        r2 = solve(BarnettABMProblem(p; save_daily = true); rng = Xoshiro(7))
        @test r2.survivors == r.survivors && r2.zoospores == r.zoospores
        # full coverage with 100% anti-infection resistance stops all new infections:
        # only the initially infected (who are not vaccinated) carry Bd
        v = barnett_vaccine(:infection, 1.0; coverage = 1.0)
        rv = solve(BarnettABMProblem(p, v); rng = Xoshiro(8))
        @test rv.n_vaccinated == p.n_tadpoles - round(Int, p.init_prevalence * p.n_tadpoles)
    end

    @testset "tolerance amplifies intensity (Fig. 1c)" begin
        c = barnett_summary(p; nreps = 100, rng = Xoshiro(11))
        @test 0 <= c.prevalence <= 1 && 0 <= c.metamorph_prevalence <= 1
        d = barnett_relative_change(p, barnett_vaccine(:tolerance, 1.0; coverage = 1.0);
                                    control = c, rng = Xoshiro(12))
        @test d.intensity > 0
        e = barnett_relative_change(p, barnett_vaccine(:infection, 0.8; coverage = 0.8);
                                    control = c, rng = Xoshiro(13))
        @test e.intensity < 0 && e.survivors > 0
        @test_throws ArgumentError barnett_relative_change(p, BarnettVaccine(); metric = :x,
                                                           control = c)
    end
end
