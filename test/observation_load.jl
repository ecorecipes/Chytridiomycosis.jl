# Load observation layer: units map, load likelihoods, Taylor's power law, observation
# operators (src/observation/load_*.jl).
using Chytridiomycosis, Test
using Distributions: Distributions, LogNormal, Gamma, NegativeBinomial, Normal, logpdf, logcdf,
    cdf, pdf, mean, var
using LogDensityProblems: LogDensityProblems
using Random: Xoshiro

const CL = Chytridiomycosis

@testset "Load observation layer" begin
    @testset "units map: Longo et al. 2013 strains" begin
        S = longo_strains()
        @test S == ["JEL404", "MexMkt", "JEL427", "PAB01", "LFT001_01", "CLFT023", "CLFT024", "LBAbercrom"]
        n = strain_copy_number.(S)
        @test n == [39, 10, 65, 59, 125, 144, 118, 107]
        q = parameter_set(:longo2013)
        @test extrema(n) == (q.copies_per_zoospore_min, q.copies_per_zoospore_max)
        @test all(r -> r.verified, parameter_records(:longo2013))
        @test_throws ArgumentError strain_copy_number("JEL423")
        # standard curves: efficiency and inversion
        c = QPCRStandardCurve(:longo2013, "JEL404")
        @test (c.slope, c.intercept, c.basis) == (-3.13, 31.6, :zoospores)
        @test round(100 * qpcr_efficiency(c)) == 109
        @test qpcr_efficiency(-1 / log10(2)) ≈ 1
        @test standard_quantity(c, cycle_threshold(c, 37.5)) ≈ 37.5
        @test cycle_threshold(c, 1) == c.intercept
        @test QPCRStandardCurve(:longo2013, "LBAbercrom"; standards = :dna).basis === :dna_ng_per_uL
        @test_throws ArgumentError QPCRStandardCurve(3.3, 30.0, :zoospores)
        @test_throws ArgumentError QPCRStandardCurve(:longo2013, "JEL404"; standards = :amplicon)
        @test_throws DomainError cycle_threshold(c, 0)
        @test CL._longo_efficiency_consistent(:A) == CL._longo_efficiency_consistent(:B) == 8
        @test CL._longo_highest_ct("MexMkt") && !CL._longo_highest_ct("JEL404")
        # restandardising ZE between strains goes through ITS1 copies
        x = ZEPerSwab(1000)
        @test load_value(restandardize(x, "JEL404" => "CLFT023")) ≈ 1000 * 39 / 144
        @test load_value(restandardize(x, "JEL404" => "MexMkt")) ≈ 3900
        @test restandardize(restandardize(x, "PAB01" => "CLFT024"), "CLFT024" => "PAB01") ≈ x
        @test restandardize(x, 169 => 169) ≈ x            # Chestnut et al.'s own factor
        @test_throws ArgumentError restandardize(x, 0 => 10)
        # DNA masses (Longo Methods)
        @test dsdna_mass(400) ≈ 4.105e-10 rtol = 1e-3
        @test dsdna_mass(400; bp_mass = 660) ≈ 400 * 660 / 6.02214076e23 * 1e9
        @test genome_equivalents(0.3; genome_bp = 46e6, bp_mass = 660) ≈ 5951 rtol = 1e-3
        @test_throws ArgumentError dsdna_mass(0)
        # viability of shed zoospores (Maguire et al. 2016, Table 2)
        @test shed_viability() ≈ 1 / (1 + exp(-1.29))
        @test shed_viability(18) > shed_viability(0)
        @test shed_viability(35; reference = false) ≈ 1 / (1 + exp(-(1.29 + 0.2 + 0.01 * 35)))
        @test parameter_set(:chestnut2014).copies_per_zoospore > q.copies_per_zoospore_max
    end

    @testset "hurdle likelihood" begin
        rng = Xoshiro(11)
        π0, μ0, σ0 = 0.6, log(200.0), 1.8
        y = [rand(rng) < π0 ? exp(μ0 + σ0 * randn(rng)) : 0.0 for _ in 1:600]
        m = LoadHurdleModel(y)
        @test LogDensityProblems.dimension(m) == 3
        @test parameter_names(m) == (:logit_prevalence, :meanlog, :log_sdlog)
        # the log density is the sum of the hand-written two-part terms
        x = [0.3, 5.0, log(2.0)]
        p = 1 / (1 + exp(-0.3))
        manual = sum(v == 0 ? log(1 - p) : log(p) + logpdf(LogNormal(5.0, 2.0), v) for v in y)
        @test LogDensityProblems.logdensity(m, x) ≈ manual
        @test sum(pointwise_loglikelihood(m, x)) ≈ manual
        @test_throws DimensionMismatch LogDensityProblems.logdensity(m, [0.0, 1.0])
        # maximum likelihood recovers the generating values; the two-part MLE is closed-form
        f = maximum_likelihood(m)
        pos = log.(filter(>(0), y))
        @test f.prevalence ≈ count(>(0), y) / length(y) atol = 1e-4
        @test f.distribution.μ ≈ sum(pos) / length(pos) atol = 1e-3
        @test f.distribution.σ ≈ sqrt(sum(abs2, pos .- sum(pos) / length(pos)) / length(pos)) atol = 1e-3
        @test f.aic ≈ 6 - 2f.loglik
        # lognormal data: the lognormal beats the gamma and the Weibull
        @test f.aic < maximum_likelihood(LoadHurdleModel(y; family = :gamma)).aic
        @test f.aic < maximum_likelihood(LoadHurdleModel(y; family = :weibull)).aic
        # censoring below a limit of detection
        lod = 50.0
        yc = [v < lod ? 0.0 : v for v in y]
        mc = LoadHurdleModel(yc; lod)
        d = LogNormal(5.0, 2.0)
        manualc = sum(v == 0 ? log((1 - p) + p * cdf(d, lod)) : log(p) + logpdf(d, v) for v in yc)
        @test LogDensityProblems.logdensity(mc, x) ≈ manualc
        fc = maximum_likelihood(mc)
        @test fc.prevalence ≈ π0 atol = 0.08          # recovers infected hosts below the LOD
        @test fc.distribution.μ ≈ μ0 atol = 0.25
        # the count scale: every family gives P(K = ⌈y⌉), so the NB is comparable
        mn = LoadHurdleModel(y; family = :negbinomial)
        @test mn.scale === :count
        nb = NegativeBinomial(0.5, 0.5 / (0.5 + 300))
        xn = [0.3, log(300.0), log(0.5)]
        manualn = sum(v == 0 ? log((1 - p) + p * pdf(nb, 0)) : log(p) + logpdf(nb, ceil(v)) for v in y)
        @test LogDensityProblems.logdensity(mn, xn) ≈ manualn
        ml = LoadHurdleModel(y; scale = :count)
        manuall = sum(v == 0 ? log(1 - p) :
                      log(p) + log(cdf(d, ceil(v)) - cdf(d, ceil(v) - 1)) for v in y)
        @test LogDensityProblems.logdensity(ml, x) ≈ manuall rtol = 1e-8
        @test maximum_likelihood(ml).aic < maximum_likelihood(mn).aic
        @test_throws ArgumentError LoadHurdleModel(y; family = :negbinomial, scale = :continuous)
        @test_throws ArgumentError LoadHurdleModel([-1.0, 2.0])
        @test_throws ArgumentError LoadHurdleModel(y; family = :poisson)
        # typed loads
        mt = LoadHurdleModel(ITS1CopiesPerSwab.(y); lod = ITS1CopiesPerSwab(1.0))
        @test mt.y ≈ y && mt.lod ≈ 1.0
        @test hurdle_distribution(m, x).prevalence ≈ p
        # positive-only fits as Schrock et al. (closed-form lognormal MLE)
        cmp = compare_load_families(y)
        @test first(cmp).family === :lognormal && first(cmp).ΔAIC == 0
        @test [c.family for c in cmp] ⊆ [:lognormal, :gamma, :weibull, :exponential]
        g = fit_positive_loads(rand(Xoshiro(2), Gamma(2.0, 50.0), 2000), :gamma)
        @test g.distribution.α ≈ 2 rtol = 0.1
        @test fit_positive_loads(y, :exponential).distribution.θ ≈ sum(filter(>(0), y)) / count(>(0), y)
        @test_throws ArgumentError fit_positive_loads([0.0, 1.0], :lognormal)
    end

    @testset "Taylor's power law" begin
        # an exact power law is recovered exactly
        ms = exp10.(range(-1, 6; length = 12))
        s = solve(TaylorPowerLawProblem(ms, 3.0 .* ms .^ 1.7))
        @test s.slope ≈ 1.7 && s.a ≈ 3.0 && s.r2 ≈ 1 && s.n == 12
        # with noise, the interval covers the generating slope
        rng = Xoshiro(5)
        vs = 2.0 .* ms .^ 1.9 .* exp10.(0.3 .* randn(rng, length(ms)))
        s2 = solve(TaylorPowerLawProblem(ms, vs))
        @test abs(s2.slope - 1.9) < 4 * s2.se_slope
        @test s2.ci_slope[2] - s2.slope ≈ s2.se_slope * Distributions.quantile(Distributions.TDist(10), 0.975)
        # lognormal loads with a common spread: slope 2 (Schrock et al.)
        @test CL._lognormal_taylor_slope() ≈ 2 atol = 1e-12
        # zeros included: Var/mean² = (e^{σ²} − π)/π
        d = LogNormal(3.0, 1.2)
        mv = load_moments(0.4, d)
        @test mv.mean ≈ 0.4 * mean(d)
        @test mv.variance / mv.mean^2 ≈ (exp(1.2^2) - 0.4) / 0.4
        # groups of swabs: Schrock et al.'s filters (≥ 10 swabs, ≥ 3 positives)
        groups = [[rand(rng) < 0.5 ? exp(μ + 1.5randn(rng)) : 0.0 for _ in 1:20] for μ in 1:15]
        push!(groups, [0.0, 0.0, 5.0], zeros(30))
        p = TaylorPowerLawProblem(groups)
        @test length(p.means) == 15
        @test p.variances[1] ≈ var(groups[1])
        @test_throws ArgumentError TaylorPowerLawProblem([1.0, 2.0], [1.0, 2.0])
        @test_throws ArgumentError TaylorPowerLawProblem([1.0, 2.0, 0.0], [1.0, 2.0, 3.0])
        @test_throws ArgumentError solve(p; level = 1.5)
    end

    @testset "Schrock et al. group summaries (synthetic file)" begin
        @test_throws ArgumentError load_schrock2024(joinpath(tempdir(), "no_such_schrock.csv"))
        @test endswith(schrock2024_path(), joinpath("data", "schrock2024", "analysis_aggregation_dataset.csv"))
        path = tempname() * ".csv"
        rng = Xoshiro(9)
        open(path, "w") do io
            println(io, "site,location,region,species,life_stage,N,num_infected,mean_natural_with_zeros,var_natural_with_zeros,dataset,prev,phase,site_id")
            for (ds, b) in (("sierra_nevada", 1.9), ("brazil", 2.06)), k in 1:20
                m = exp10(4 * rand(rng))
                st = isodd(k) ? "adult" : "larva"
                println(io, "s$k,usa,r,\"Rana, sp.\",$st,20,$(k % 7),$m,$(2 * m^b),$ds,NA,",
                    ds == "sierra_nevada" ? (isodd(k) ? "enzootic" : "epizootic") : "Not Available", ",NA")
            end
            println(io, "x,usa,r,sp,unknown,20,5,NA,NA,east_bay,NA,Not Available,NA")
        end
        gs = load_schrock2024(path)
        @test length(gs) == 41 && gs[1].species == "Rana, sp." && isnan(gs[end].mean)
        sl = schrock_taylor_slopes(gs)
        @test sort(collect(keys(sl))) == [:brazil, :sierra_nevada]
        @test sl[:sierra_nevada].slope ≈ 1.9 && sl[:brazil].slope ≈ 2.06
        @test sl[:brazil].n == count(g -> g.dataset === :brazil && g.num_infected >= 3, gs)
        @test haskey(schrock_taylor_slopes(gs; by = :phase), :enzootic)
        @test sort(collect(keys(schrock_taylor_slopes(gs; by = :life_stage)))) == [:adult, :larva]
        rm(path)
    end

    @testset "Schrock et al. group summaries (Dryad file)" begin
        # quoted fields may hold commas, doubled quotes and line breaks (the Dryad file's phase
        # labels do), and CRLF endings are accepted
        rows = Chytridiomycosis._schrock_csv_rows("a,b,c\r\n1,\"x, \"\"y\"\"\",\"two\nlines\"\r\n\n2,z,w\n")
        @test rows == [["a", "b", "c"], ["1", "x, \"y\"", "two\nlines"], ["2", "z", "w"]]
        gs = load_schrock2024()
        @test length(gs) == 1861
        @test count(g -> occursin('\n', g.phase), gs) > 0
        f = schrock_family_fractions()
        @test f.n == 524
        @test sum(values(f.clearly_best)) + f.several ≤ 1 + 1e-12
        @test all(v -> 0 ≤ v ≤ 1, values(f.as_good)) && f.as_good[:lognormal] > f.as_good[:weibull]
        @test_throws ArgumentError schrock_family_fractions(joinpath(tempdir(), "no_such_schrock.csv"))
    end

    @testset "observation operators" begin
        # a hand-built distribution: one uniform and one point component with noise
        d = SwabLoadDistribution(0.5, [1.0, 3.0], [2.0, 3.0], [1.0, 3.0]; σ = 0.4, lod = 10.0)
        @test d.weights == [0.25, 0.75]
        # P(non-detect) + ∫ density over detected loads = 1 (density on the natural scale)
        zs = range(1.0, 7.0; length = 60_001)
        dens = [exp(swab_logpdf(d, exp10(z))) * exp10(z) * log(10) for z in zs]
        @test sum(dens) * step(zs) + exp(swab_logpdf(d, 0.0)) ≈ 1 atol = 1e-5
        @test swab_logpdf(d, 5.0) == swab_logpdf(d, 0.0)        # below the LOD: non-detect
        @test swab_logpdf(d, ZEPerSwab(300)) == swab_logpdf(d, 300.0)
        @test swab_loglikelihood(d, [0.0, 300.0]) ≈ swab_logpdf(d, 0.0) + swab_logpdf(d, 300.0)
        # simulation agrees with detection_probability and load_moments
        r = rand(Xoshiro(1), d, 200_000)
        @test count(>(0), r) / length(r) ≈ detection_probability(d) atol = 0.005
        mv = load_moments(d)
        @test sum(r) / length(r) ≈ mv.mean rtol = 0.03
        @test var(r) ≈ mv.variance rtol = 0.1
        # without noise the uniform component's detection is its share above the LOD
        d0 = SwabLoadDistribution(1.0, [0.0], [2.0], [1.0]; lod = 10.0)
        @test detection_probability(d0) ≈ 0.5
        @test_throws ArgumentError swab_logpdf(SwabLoadDistribution(1.0, [1.0], [1.0], [1.0]), 10.0)
        @test_throws ArgumentError SwabLoadDistribution(1.2, [1.0], [1.0], [1.0])

        # host loads in ZE per swab
        hl = swab_distribution(ZEPerSwab.([0.0, 0.0, 100.0, 1e4]); σ = 0.2)
        @test hl.prevalence == 0.5 && hl.lower == [2.0, 4.0]
        # Briggs sporangia: refused without a stated conversion, warned with one
        sp = SporangiaPerHost.([0.0, 50.0, 500.0])
        @test_throws ArgumentError swab_distribution(sp)
        b = @test_logs (:warn,) swab_distribution(sp; ze_per_sporangium = 2.0, σ = 0.3)
        @test b.lower ≈ log10.([100.0, 1000.0])
        # a mean load with an assumed lognormal spread keeps the arithmetic mean
        ml = swab_distribution(ZEPerSwab(220.0); sdlog10 = 0.95, prevalence = 0.7)
        @test load_moments(ml).mean ≈ 0.7 * 220
        @test load_moments(ml).variance ≈ 0.7 * 220^2 * exp((0.95 * log(10))^2) - (0.7 * 220)^2 rtol = 1e-10

        # Wilber et al. 2016 IPM: the state re-expressed, no extra assumption by default
        p16 = WilberParameters(:wilber2016; set = "deposited")
        sol = solve(BdIPMProblem(p16; temperature = 16))
        k = length(sol.t)
        w = swab_distribution(sol)
        @test w.prevalence ≈ sum(sol.I[:, k]) / (sol.S[k] + sum(sol.I[:, k]))
        @test w.lower[1] ≈ sol.mesh.lower / log(10) && w.upper[end] ≈ sol.mesh.upper / log(10)
        @test detection_probability(w) ≈ w.prevalence
        @test detection_probability(swab_distribution(sol; lod = 1e4)) < w.prevalence
        mlog = sum(sol.I[:, k] .* meshpoints(sol.mesh)) / sum(sol.I[:, k])
        @test sum(w.weights .* (w.lower .+ w.upper) ./ 2) ≈ mlog / log(10)
        @test_throws BoundsError swab_distribution(sol; step = k + 1)

        # Wilber et al. 2017 hybrid IPM: adults
        p17 = Wilber2017Parameters(:wilber2017; set = "coded")
        h = solve(HybridIPMProblem(p17; years = 1, stochastic = false))
        hs = hybrid_summary(h)
        j = findlast(>(0), hs.infected)
        hw = swab_distribution(h; step = j)
        @test hw.prevalence ≈ hs.prevalence[j]

        # Louca et al. 2014: binned DOI, bin 1 uninfected
        pl = LoucaParameters(:louca2014; set = "fig2cd")
        out = solve(LoucaDOIProblem(pl; years = 2, save_distributions = true, saveat = 30.0); dt = 0.05)
        lw = swab_distribution(pl, out; σ = 0.3)
        @test lw.prevalence ≈ out.prevalence_A[end]
        @test length(lw.weights) == pl.m - 1
        @test swab_distribution(pl, out; stage = :tadpole, step = 2).prevalence ≈ out.prevalence_T[2]
        @test_throws ArgumentError swab_distribution(pl, out; stage = :egg)
    end
end
