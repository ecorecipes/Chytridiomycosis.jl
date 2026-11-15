# Hollanders et al. 2023 (src/models/hollanders2023.jl): bundled Dryad data (CC0), the
# continuous-time robust-design likelihood against the generic HMM layer, the load submodel,
# and (when the authors' draws are present) the reporting of their posterior.
using LogDensityProblems: LogDensityProblems
using Random: Xoshiro

@testset "Hollanders et al. 2023" begin
    d = load_hollanders2023()
    @test d.sites == ["Brindle", "Tuntable", "Bat Cave"] && d.nind == [462, 136, 88]
    @test size(d.ch) == (462, 23, 3, 3) && size(d.intervals) == (22, 3)
    @test all(ismissing, d.ch[137:end, :, :, 2])                 # padding rows
    @test all(x -> x isa ITS1CopiesPerSwab, skipmissing(d.loads))
    @test d.intervals[1, 1] == 0.5 && d.intervals[5, 2] == 5.452

    s = hollanders2023_summary_data(d)
    @test sum(s.individuals) == 686
    @test s.primary_detections == 1609
    @test s.detections == 2384
    @test s.infected_detections == 248 && length(s.loads) == 248
    @test round(s.minimum; digits = 3) == 1.425 && round(s.maximum; digits = 2) == 6.18
    @test round(s.mean; digits = 2) == 3.17

    @test load_hazard_ratio(3.174; status = 0.0, load = 2.0, mean = 3.174, sd = 0.9) == 1
    @test load_hazard_ratio(3.174 + 1.8; status = 0.1, load = 2.0, mean = 3.174, sd = 0.9) ≈ exp(2.1)

    # Gauss–Hermite rule integrates normal moments exactly
    z, w = Chytridiomycosis._gauss_hermite(5)
    @test sum(w) ≈ 1 && sum(w .* z .^ 2) ≈ 1 && sum(w .* z .^ 4) ≈ 3

    lm = hollanders2023_load_model(d)
    @test lm.temperature < 0 && lm.rain > 0                       # authors: −0.359, 0.231
    @test 0.3 < lm.sd < 0.5 && lm.sd_individual < 0.35            # authors: 0.41, 0.24

    # the fast recursion is the HMM with states (U, I at each load node, dead)
    prob = CTMultistateProblem(d; nodes = (load = 3, individual = 2, detection = 3))
    @test LogDensityProblems.dimension(prob) == 26
    @test all(i -> i.first < 23, prob.inds)
    θ = hollanders2023_start() .+ 0.2 .* randn(Xoshiro(1), 26)
    gh = Chytridiomycosis._h23_nodes(prob)
    for ind in prob.inds[[1, 50, 300, 500, end]]
        tot = 0.0
        sexes = ind.sex === missing ? ((1, prob.female_fraction), (0, 1 - prob.female_fraction)) :
                ((ind.sex, 1.0),)
        for (a, u) in enumerate(gh.individual[1]), (b, e) in enumerate(gh.detection[1]),
            (sx, ps) in sexes
            init, Ps, L = Chytridiomycosis._h23_hmm(prob, θ, ind, u, exp(θ[26]) * e, sx)
            tot += gh.individual[2][a] * gh.detection[2][b] * ps *
                   exp(hmm_loglikelihood(init, Ps, L))
        end
        @test Chytridiomycosis._h23_loglik(prob, θ, ind, gh) ≈ log(tot)
    end
    @test isfinite(LogDensityProblems.logdensity(prob, θ))
    # swabbing on the last primary changes the likelihood only through Bd+ frogs seen then
    p2 = CTMultistateProblem(d; nodes = (load = 3, individual = 2, detection = 3), swab_final = true)
    @test LogDensityProblems.logdensity(p2, θ) != LogDensityProblems.logdensity(prob, θ)

    if Chytridiomycosis._hollanders2023_draws_find() !== nothing
        p = hollanders2023_posterior()
        @test round(p.phi_status.median; digits = 3) == -0.167
        @test round(p.phi_load.median; digits = 3) == 1.953
        @test round(p.annual_hazard.median; digits = 3) == 0.783
        @test round(p.gain.median; digits = 3) == 0.283
        @test round(p.clearance.median; digits = 3) == 2.623
        @test round.(extrema(p.mortality_sites); digits = 3) == (0.077, 0.104)
    else
        @info "Hollanders et al. 2023 posterior draws not found; skipping their checks"
    end
end
