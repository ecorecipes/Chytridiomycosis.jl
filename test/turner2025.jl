# Turner et al. 2025 Conserv Sci Pract 7:e70022 (doi:10.1111/csp2.70022): age-classified PVA
# of Litoria raniformis under environmental-watering scenarios. The oracles are
# - the published paper (read 2026-11-20): Methods 2.2 values (checked in the compendium) and
#   the extinction and quasi-extinction probabilities (validation targets, mostly discrepant);
# - the authors' code and archive (Zenodo doi:10.5281/zenodo.10921748, CC BY 4.0): the first
#   rows of the archived survival arrays; the archived extinction probabilities and minimum
#   population sizes (in data/validation/targets.toml);
# - one-off R runs (2026-11-19) of the authors' nls fits and integrate() calls:
#   integrate(eggfunc, 0, 1.02)/1.02 = 0.917641850605293, integrate(eggfunc, 0, 0.5) =
#   0.499574951491747;
# - closed forms: the stationary wet fraction of the capped watering rules;
# - the abstract (extinction falls from August to October watering).
using Chytridiomycosis
using Test
using Random

@testset "Turner 2025 watering PVA" begin
    C = Chytridiomycosis

    @testset "Turner density dependence and watering" begin
        p = TurnerParameters(:turner2025)
        recs = parameter_records(:turner2025)
        vr = Dict((r.set, r.symbol) => r.verified for r in recs)
        @test vr[("", :years)] && vr[("", :fecundity_sd)] && vr[("", :p_wet)]
        @test !vr[("", :spawn_capacity)] && !vr[("", :egg_a)]
        @test p.winter_daily_survival^p.winter_duration ≈ 0.692 atol = 5e-4  # paper, 2.2.1
        @test TurnerParameters(:turner2025; as_printed = true).fecundity_mean == 3244
        @test p.fecundity_mean == 3245
        @test p.years == 52
        # fitted curves pass near the authors' anchor points (compendium notes)
        @test turner_egg_multiplier(p, 0.01) ≈ 1.0 atol = 0.05
        @test turner_egg_multiplier(p, 1.0) ≈ 0.08 atol = 0.05
        @test turner_tadpole_multiplier(p, 1.0) ≈ 0.10 atol = 0.02
        @test turner_tadpole_multiplier(p, 0.3) ≈ 1.0 atol = 0.05
        @test turner_juvenile_multiplier(p, 0.6) ≈ 0.4 atol = 0.15
        # closed-form integral against R's integrate()
        K = p.spawn_capacity * p.clutch_size_mean * p.sex_ratio
        @test C._turner_logistic_integral(p.egg_a, p.egg_b, p.egg_c, 1.02) / 1.02 ≈
              0.917641850605293 rtol = 1e-6
        @test C._turner_logistic_integral(p.egg_a, p.egg_b, p.egg_c, 0.5) ≈
              0.499574951491747 rtol = 1e-6
        @test turner_egg_ceiling(p, 0.5K) ≈ 0.5K * 0.499574951491747 / 0.5 rtol = 1e-6
        @test turner_egg_ceiling(p, 2K) ≈
              K * 0.917641850605293 + turner_egg_multiplier(p, 1.02) * K rtol = 1e-6
        @test turner_egg_ceiling(p, 0.0) == 0
        # watering rules: year 1 wet, run caps, stationary wet fractions
        rng = Xoshiro(9)
        for (pat, cap, frac) in ((:max1dry, 1, 1 / 1.3), (:max2dry, 2, 1 / 1.39),
                                 (:random, 52, 0.7))
            W = [turner_watering_sequence(rng, pat, 52) for _ in 1:4000]
            @test all(w -> w[1], W)
            longest(w) = maximum(accumulate((r, x) -> x ? 0 : r + 1, w; init = 0))
            @test maximum(longest, W) <= cap
            @test C.mean(C.mean(w[2:end]) for w in W) ≈ frac atol = 0.01
        end
        @test_throws ArgumentError turner_watering_sequence(rng, :drought, 10)
    end

    @testset "Turner PVA" begin
        d = load_turner2025(; month = :august, pattern = :random)
        @test size(d.survival) == (52, 4) && size(d.start) == (6, 1000)
        # the as-coded series is the archived array's last replicate (first year)
        @test d.survival[1, 1] ≈ 0.193913058584897
        @test d.survival[1, 3] ≈ 0.626924448214021
        @test all(sum(d.start; dims = 1) .> 0)
        @test_throws ArgumentError load_turner2025(; month = :june)
        prob = TurnerPVAProblem(:turner2025; month = :october, pattern = :max1dry,
            replicates = 300)
        sol = solve(prob; rng = Xoshiro(1))
        @test size(sol.population) == (53, 300) && size(sol.terminal) == (6, 300)
        @test all(sol.population .>= 0) && all(isinteger, sol.population)
        @test turner_extinction_probability(sol) <= 0.02
        # extinct replicates stay at zero
        e = TurnerPVASolution([5.0 3; 0 2; 0 1], zeros(6, 2))
        @test turner_extinction_probability(e) == 0.5
        @test turner_minimum_population(e) == [0.0, 1.0]
        @test turner_quasi_extinction_probability(e, 1) == 0.5
        @test turner_quasi_extinction_probability(e, 20) == 1
        # watering later in the season lowers extinction risk (abstract)
        ext(m) = turner_extinction_probability(solve(TurnerPVAProblem(:turner2025;
            month = m, replicates = 3000); rng = Xoshiro(4)))
        @test ext(:august) > ext(:october) + 0.05
        # wetness given explicitly: all dry years means no recruitment, so extinction
        p = TurnerParameters(:turner2025)
        dry = falses(200, 52)
        dry[:, 1] .= true
        sd = solve(TurnerPVAProblem(p; survival = d.survival, start = d.start, wetness = dry,
            replicates = 200); rng = Xoshiro(2))
        @test turner_extinction_probability(sd) == 1
        @test_throws ArgumentError TurnerPVAProblem(p; survival = zeros(3, 4), start = d.start)
        # per-replicate survival rows (as intended), if the authors' archive is extracted
        dir = joinpath(pkgdir(Chytridiomycosis), "..", "..", "papers", "turner2025",
            "zenodo", "TurnerBdModel")
        dir = get(ENV, "CHYTRIDIOMYCOSIS_TURNER2025", dir)
        f = joinpath(dir, "referenceData", "AugustWateringArray_RandomDryYears.csv")
        if isfile(f)
            A = load_turner2025_survival_array(:august, :random; path = dir)
            @test size(A) == (10_000, 52, 5) && A[10_000, 1, 2] == d.survival[1, 1]
            s2 = solve(TurnerPVAProblem(:turner2025; replicates = 2000, survival_array = A);
                rng = Xoshiro(3))
            @test 0.1 < turner_extinction_probability(s2) < 0.22
        else
            @info "Turner 2025 survival arrays not extracted; skipping the per-replicate test"
        end
    end
end
