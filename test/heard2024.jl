# Heard et al. 2024 J Anim Ecol 93:1670 (doi:10.1111/1365-2656.14177): catch-curve survival
# and spatial metapopulation of Australian Litoria. The full text could not be retrieved; the
# oracles are
# - the authors' code and archive (Zenodo doi:10.5281/zenodo.10457229, CC BY 4.0): sample
#   sizes in Age_structure_analysis.R; survival values entered in Run_STEPS_simulations.R;
#   the archived simulation outputs (persistence fractions, in data/validation/targets.toml);
# - the R package steps 1.3.0 source (dispersal neighbourhood incl. the id2coord shift);
# - closed forms: the geometric likelihood and its MLE, the lognormal moments of the survival
#   SD (against the authors' Monte Carlo recipe);
# - the abstract (halving of adult survival; persistence contracts post-Bd).
using Chytridiomycosis
using Test
using Random
using Distributions: Geometric, logpdf, Normal
using LogDensityProblems

@testset "Heard 2024 age truncation and metapopulation" begin
    C = Chytridiomycosis

    @testset "Heard data and compendium" begin
        d = load_heard2024()
        @test length(d) == 403
        nrec(sp, ep) = count(r -> r.species === sp && r.epoch === ep, d)
        # counts in Age_structure_analysis.R (L. raniformis: 45 historical, 254 contemporary)
        @test nrec(:aurea, :historical) == 55 && nrec(:aurea, :contemporary) == 49
        @test nrec(:raniformis, :historical) == 45 && nrec(:raniformis, :contemporary) == 254
        @test length(unique(r.site for r in d)) == 15
        @test all(r -> r.age >= 1, d)
        L = heard2024_landscapes()
        @test length(L.habitat) == 20 && all(h -> size(h) == (25, 25), L.habitat)
        @test sum(L.habitat[1]) == 1 && sum(L.habitat[20]) == 25
        @test L.connectivity[1] == 0 && L.connectivity[end] == 85
        recs = parameter_records(:heard2024)
        @test !any(r -> r.verified, recs)          # paper not available: code values only
        g = heard_scenario_grid()
        @test g.mets == 3:3:30 && length(g.mets_cv) == 10 && g.mortality_cv[end] ≈ 1.0
        for (sp, ep, S) in ((:raniformis, :contemporary, "Lran_contemporary"),
                            (:aurea, :historical, "Laur_historical"))
            p = HeardParameters(:heard2024; species = sp, epoch = ep)
            v = parameter_set(:heard2024; set = S)
            @test p.adult_survival == v.adult_survival == p.subadult_survival
            @test p.carrying_capacity == v.carrying_capacity
        end
        @test_throws ArgumentError HeardParameters(:heard2024; species = :fletcheri)
    end

    @testset "catch-curve log density" begin
        m = HeardCatchCurve()
        @test LogDensityProblems.dimension(m) == 19
        @test length(C.parameter_names(m)) == 19
        θ = [0.3, -0.2, -1.1, -1.0, zeros(15)...]
        S = heard_survival(θ)
        # with zero site effects the likelihood is a sum of geometric log pmfs
        Sf(i) = 1 / (1 + exp(-(θ[1] + θ[2] * m.species[i] + θ[3] * m.epoch[i])))
        ll = sum(logpdf(Geometric(1 - Sf(i)), m.age[i]) for i in eachindex(m.age))
        @test heard_loglikelihood(m, θ) ≈ ll rtol = 1e-12
        @test S.S_c_Lr ≈ 1 / (1 + exp(-(0.3 - 0.2 - 1.1)))
        # priors: Normal(0, 1/sqrt(0.001)) on α, β; σ uniform on (0, 10) through a logistic
        σ = 10 / (1 + exp(1.0))
        lp = sum(logpdf(Normal(0, 1 / sqrt(0.001)), x) for x in θ[1:3]) +
             log(σ / 10 * (1 - σ / 10)) + 15 * logpdf(Normal(0, σ), 0.0)
        @test LogDensityProblems.logdensity(m, θ) ≈ lp + ll rtol = 1e-10
        @test_throws DimensionMismatch LogDensityProblems.logdensity(m, zeros(3))
        # one-group maximum likelihood is ā/(1 + ā): the geometric MLE
        a = m.age[(m.species .== 1) .& (m.epoch .== 1)]
        Shat = sum(a) / length(a) / (1 + sum(a) / length(a))
        f(x) = sum(logpdf(Geometric(1 - x), k) for k in a)
        @test f(Shat) > f(Shat - 0.01) && f(Shat) > f(Shat + 0.01)
        # a short posterior sample recovers the survivals the authors simulated with
        fit = solve(m; rng = Xoshiro(1), warmup = 2000, iterations = 8000, thin = 4)
        s = heard_survival_summary(fit)
        for (k, set) in ((:S_c_Lr, "Lran_contemporary"), (:S_h_Lr, "Lran_historical"),
                         (:S_c_La, "Laur_contemporary"), (:S_h_La, "Laur_historical"))
            @test abs(s[k].mean - parameter_set(:heard2024; set).adult_survival) < 0.04
        end
        @test 0.4 < s.S_c_Lr.mean / s.S_h_Lr.mean < 0.6        # "a halving" (abstract)
    end

    @testset "Heard metapopulation" begin
        p = HeardParameters(:heard2024; species = :raniformis, epoch = :contemporary,
            mets = 24, mets_cv = 0.5, mortality_cv = 0.5)
        A = heard_projection_matrix(p)
        @test A ≈ [24 * 0.2 * 0.5 * 0.1 24 * 0.5 * 0.1; 0.319 0.319]
        @test heard_lambda(p) ≈ maximum(real, C.eigvals(A))
        # survival SD: exact lognormal moments against the authors' Monte Carlo recipe
        rng = Xoshiro(3)
        h = 0.5 .* (-log(0.319)) .* randn(rng, 400_000) .+ (-log(0.319))
        @test heard_survival_sd(0.319, 0.5) ≈ sqrt(C.var(exp.(-h))) rtol = 0.01
        @test_throws ArgumentError heard_survival_sd(1.0, 0.5)
        # dispersal neighbourhoods: weights sum to 1, the origin is included, and the steps
        # id2coord quirk shifts last-column cells only when as_coded
        H = falses(25, 25)
        H[3, 25] = H[4, 1] = H[3, 24] = true
        patches = [(r, c) for r in 1:25 for c in 1:25 if H[r, c]]
        for coded in (true, false)
            dest = C._heard_destinations(HeardMetapopulationProblem(p, H; as_coded = coded),
                patches)
            @test all(d -> sum(d[2]) ≈ 1, dest)
            i = findfirst(==((3, 25)), patches)
            j = findfirst(==((4, 1)), patches)
            w = Dict(zip(dest[i]...))
            # as coded, (3, 25) sits at (4, 0): (4, 1) is one cell away and the origin
            # itself is out of range (1637 m); otherwise the reverse
            @test coded ? (haskey(w, j) && !haskey(w, i)) : (haskey(w, i) && !haskey(w, j))
        end
        buf = zeros(Int, 5)
        C._heard_multinomial!(buf, Xoshiro(1), 37, [0.1, 0.0, 0.5, 0.2, 0.2])
        @test sum(buf) == 37 && buf[2] == 0
        # one simulation: initial state round(K/2)/2 per stage and patch, conserved shapes
        sol = solve(HeardMetapopulationProblem(p, 20); rng = Xoshiro(2))
        @test size(sol.adults) == (51, 25)
        @test all(==(42), sol.subadults[1, :]) && all(==(42), sol.adults[1, :])
        # the landscape-wide multinomial rounding of steps can push a patch above K; only
        # the total is bounded
        @test all(sum(sol.subadults .+ sol.adults; dims = 2) .<= 25 * p.carrying_capacity)
        solx = solve(HeardMetapopulationProblem(p, 20; as_coded = false); rng = Xoshiro(2))
        @test all(solx.subadults .+ solx.adults .<= p.carrying_capacity)
        @test sol.extinct == (sol.total[end] == 0)
        # a single isolated patch: dispersers can only return to the origin
        sol1 = solve(HeardMetapopulationProblem(p, 1); rng = Xoshiro(2))
        @test size(sol1.adults, 2) == 1
        # pre-Bd survival persists far more often than post-Bd (abstract; archived outputs)
        rng = Xoshiro(5)
        hist = heard_persistence_grid(; species = :raniformis, epoch = :historical,
            mets = [12], mets_cv = [0.3, 0.6], mortality_cv = [0.3, 0.6], landscapes = 15:20,
            rng)
        cont = heard_persistence_grid(; species = :raniformis, epoch = :contemporary,
            mets = [12], mets_cv = [0.3, 0.6], mortality_cv = [0.3, 0.6], landscapes = 15:20,
            rng)
        @test hist > 0.8 && cont < 0.1
        @test 0 <= heard_persistence(p; landscapes = 18:20, rng) <= 1
        @test 0 <= heard_persistence(p; landscapes = 18:20, as_coded = false, rng) <= 1
    end

end
