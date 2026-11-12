# Briggs, Knapp & Vredenburg 2010 PNAS: the six-stage R. muscosa / R. sierrae model of Fig. 4
# (SI Methods, "Model with full R. muscosa/R. sierrae stage structure"). The paper shows one
# run per panel and describes it in words, so the oracles are those descriptions and the
# levels read off Fig. 4; the stochastic details the paper leaves open are listed at the top
# of src/models/briggs2010_stages.jl.
using Chytridiomycosis
using Test
using Random: Xoshiro

@testset "Briggs 2010 six-stage model (Fig. 4)" begin
    C = Chytridiomycosis
    pa = BriggsLifeStages(:briggs2010; set = "fig4a")
    pb = BriggsLifeStages(:briggs2010; set = "fig4b")
    pc = BriggsLifeStages(:briggs2010; set = "fig4c")

    @testset "parameters (Fig. 4 legend)" begin
        for set in ("fig4a", "fig4b")
            q = parameter_set(:briggs2010; set)
            @test q.γ_subadult ≈ 10q.γ_adult && q.γ_tadpole ≈ 100q.γ_adult   # derived records
        end
        @test pa.adult.γ == 1e-6 && pb.adult.γ == 1e-4 && pc.adult.γ == 1e-3
        @test pa.sub1.Smax == 1000 && pa.tad1.Smax == 1000 && pc.sub1.Smax == 10_000
        @test !pa.tad1.lethal && !pa.tad3.lethal && pa.sub2.lethal && pa.adult.lethal
        @test pa.adult.ν == 0.1 && pa.tad1.ν == 0.1 && pa.sub1.f == 0.1
        @test pa.adult.f == 0.1 && pa.η == 17.5 && pa.σ == 0.2 && pa.μ == 1.0 && pa.V == 1.0
        @test (pa.θZ, pa.m, pa.ω, pa.pF, pa.R, pa.K) == (0.5, 0.5, 0.9, 0.25, 100.0, 100)
        @test pa.adult.θ == 0.9 && pa.sub1.θ == 0.7 && pa.tad2.θ == 0.7 && pa.summer == 100
    end

    @testset "Bd-free stage structure" begin
        n = disease_free_stages(pa)
        @test n.adult == 100 && n.tad1 ≈ 2250 && n.tad2 ≈ 1575
        @test n.tad3 ≈ 551.25 && n.sub1 ≈ 843.4125 && n.sub2 ≈ 590.38875
        # Fig. 4 before invasion: tad1 2000-2800, tad2 1500-2000, tad3 ≈ 600, sub1 700-1000,
        # sub2 500-700, adults at K = 100
        @test 2000 < n.tad1 < 2800 && 500 < n.tad3 < 700 && 700 < n.sub1 < 1000 && 500 < n.sub2 < 700
        free = solve(BriggsStagedProblem(pa; years = 30, invade_year = 99); rng = Xoshiro(4))
        @test all(==(100), free.counts.adult)
        @test isapprox(sum(free.counts.sub1) / 30, n.sub1; rtol = 0.1)
        @test isapprox(sum(free.counts.tad1) / 30, n.tad1; rtol = 0.1)
        @test all(iszero, free.Z) && free.extinct_year === nothing && !free.persisted
        bad = BriggsLifeStages(; V = 1.0, μ = 1.0, η = 17.5, σ = 0.2, tad1 = pa.tad1,
            tad2 = pa.tad2, tad3 = pa.tad3, sub1 = pa.sub1, sub2 = pa.sub2, adult = pa.adult,
            θZ = 0.5, m = 0.5, ω = 0.9, pF = 0.0, R = 100.0, K = 100)
        @test_throws ArgumentError disease_free_stages(bad)
    end

    @testset "one summer matches the unstructured model" begin
        # deterministic steps (no Poisson draws), 30 adults and no juveniles
        ad = BriggsStage(; γ = 0.01, ν = 0.05, f = 0.05, Smax = 1e4, θ = 1.0)
        none = BriggsStage(; γ = 0.0, ν = 0.05, f = 0.05, Smax = 1e4, θ = 1.0, lethal = false)
        ls = BriggsLifeStages(; V = 1.0, μ = 1.0, η = 17.5, σ = 0.2, tad1 = none, tad2 = none,
            tad3 = none, sub1 = none, sub2 = none, adult = ad, θZ = 0.1, m = 0.5, ω = 0.9,
            pF = 0.0, R = 0.0, K = 30, summer = 50.0)
        r6 = solve(BriggsStagedProblem(ls; years = 1, invade_year = 1, Δt = 0.01,
            poisson_below = 0.0, init = (tad1 = 0, tad2 = 0, tad3 = 0, sub1 = 0, sub2 = 0,
            adult = 30)); rng = Xoshiro(1))
        sp = BriggsSeasonalParameters(; V = 1.0, μ = 1.0, η = 17.5, σ = 0.2, adult = ad,
            R = 0.0, K = 30, θZ = 0.1, summer = 50.0)
        r1 = solve(BriggsSeasonalProblem(sp; years = 1, Δt = 0.01, poisson_below = 0.0);
            rng = Xoshiro(1))
        @test r6.mean_load.adult[1] ≈ r1.mean_adult_load[1]
        @test r6.Z[1] ≈ r1.Z[1]
    end

    @testset "Fig. 4A: occasional subadult recruitment, adults never lethally loaded" begin
        r = solve(BriggsStagedProblem(pa; years = 50, record_year = 50); rng = Xoshiro(1))
        after = 21:50
        @test sum(r.deaths.adult) == 0                          # adult loads stay sublethal
        @test all(r.mean_load.adult[after] .< 1000)
        @test 50 < sum(r.mean_load.adult[after]) / 30 < 500     # Fig. 4A: ≈ 100-200
        @test all(r.mean_load.tad1[after] .== pa.tad1.Smax)     # tadpoles at their cap
        @test sum(r.deaths.sub1) > 1000                         # subadults die of Bd ...
        @test any(diff(r.counts.adult[after]) .> 20)            # ... but some years recruit
        @test minimum(r.counts.adult[after]) < 60 && r.counts.adult[end] > 0
        @test all(1e7 .< r.Z[after] .< 1e8)                     # Fig. 4A: 2-5 × 10⁷
        @test r.persisted && r.extinct_year === nothing
        w = r.within                                            # within-season (cf. year 51)
        @test size(w.mean_load) == (100, 6) && all(w.Z .> 1e6)
        @test maximum(w.max_load[:, 6]) < pa.adult.Smax
    end

    @testset "Fig. 4B: subadults never recruit, slow decline to extinction" begin
        r = solve(BriggsStagedProblem(pb; years = 50); rng = Xoshiro(1))
        @test all(iszero, r.counts.sub2[21:50]) && all(iszero, r.counts.sub1[21:50])
        @test sum(r.deaths.adult) == 0
        @test issorted(r.counts.adult[20:50]; rev = true)       # no recruitment at all
        @test r.counts.adult[30] > 30                           # extant for many years ...
        @test r.counts.adult[50] < 15                           # ... but declining
        @test all(1e6 .< r.Z[21:40] .< 2e6)                     # Fig. 4B: ≈ 1.5 × 10⁶
    end

    @testset "Fig. 4C: rapid extinction within 2-3 years" begin
        for seed in 1:3
            r = solve(BriggsStagedProblem(pc; years = 40); rng = Xoshiro(seed))
            @test r.deaths.adult[20] == 100                     # all adults die in year 20
            @test 21 <= r.extinct_year <= 23
            @test !r.persisted
            @test 5e7 < maximum(r.Z) < 3e8                      # Fig. 4C: peak ≈ 10⁸
        end
    end
end
