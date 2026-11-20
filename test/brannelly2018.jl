# Brannelly et al. 2018 BMC Ecol 18:34 (doi:10.1186/s12898-018-0189-5): conditional
# Arnason–Schwarz multistate model and POPAN for Acris crepitans. Oracles: Table 1 (capture
# counts), Table 2 (deviances, parameter counts, ranking), Results (model-averaged monthly
# survival, POPAN super-population size), read from the published full text (PMC6137908).
# Closed forms: a hand-computed Arnason–Schwarz likelihood of one history, and POPAN
# probabilities that sum to one over all histories. Data: bundled Additional files (CC0).
using Chytridiomycosis
using Test
using LogDensityProblems
using Random: Xoshiro

@testset "Brannelly 2018 multistate CMR and POPAN" begin
    C = Chytridiomycosis

    @testset "data (Table 1)" begin
        d = load_brannelly2018()
        s = brannelly2018_summary(d)
        @test size(d.histories) == (315, 10)
        @test s.capture_events == 378       # Table 1: 378 capture instances
        @test s.positive == 63              # Results: 63 Bd-positive captures
        @test s.negative == 315             # Table 1 (labels swapped): 315
        @test s.recaptured == 50            # Table 1: 50 recaptured
        @test s.individuals == 315          # Table 1 prints 266 (discrepant target)
        @test d.intervals == [1, 2, 1, 1, 2, 2, 1, 1, 1]
        @test sum(d.intervals) == 12        # Nov 2010 to Nov 2011
        dc = load_brannelly2018(; occasions = :calendar_month, intervals = :equal)
        @test size(dc.histories, 2) == 10 && all(==(1), dc.intervals)
        @test dc.histories != d.histories   # 28 Sep / 1 Oct grouped differently
    end

    @testset "Arnason–Schwarz likelihood by hand" begin
        st = brannelly2018_structure(; survival = :state, detection = :state,
                                     transition = :state)
        H = [1 0 2; 2 2 0; 1 1 1]
        τ = [1.0, 2.0]
        prob = ArnasonSchwarzProblem(H, st; intervals = τ)
        @test LogDensityProblems.dimension(prob) == 6       # S×2, ψ×2, p×2; no initial state
        θ = [1.0, 0.5, -1.0, -0.5, 0.2, 0.8]
        par = arnason_schwarz_parameters(prob, θ)
        S, p = par.survival[:, 1], par.detection[:, 1]
        ψ = par.transition[:, :, 1]
        @test ψ[1, 2] ≈ 1 / (1 + exp(1.0)) && ψ[2, 1] ≈ 1 / (1 + exp(0.5))
        # history 1 0 2: state 1 at t1; t2 unseen in a or b; seen in state 2 at t3
        Φ(t) = [S[k]^τ[t] * ψ[k, l] for k in 1:2, l in 1:2]
        h1 = sum(Φ(1)[1, a] * (1 - p[a]) * Φ(2)[a, 2] * p[2] for a in 1:2)
        # history 2 2 0: seen in 2 at t2, then dead or alive and missed at t3
        h2 = Φ(1)[2, 2] * p[2] * (1 - sum(Φ(2)[2, b] * p[b] for b in 1:2))
        h3 = Φ(1)[1, 1] * p[1] * Φ(2)[1, 1] * p[1]
        @test LogDensityProblems.logdensity(prob, θ) ≈ log(h1) + log(h2) + log(h3)
        @test pointwise_loglikelihood(prob, θ) ≈ [log(h1), log(h2), log(h3)]
        @test_throws ArgumentError ArnasonSchwarzProblem(H, C.palomar2023_structure())
    end

    @testset "Table 2 model set" begin
        ms = brannelly2018_models()
        @test length(ms) == 10
        # the parameter counts of Table 2 follow from the designs with 10 occasions
        nfree(m) = LogDensityProblems.dimension(ArnasonSchwarzProblem(load_brannelly2018().histories,
                                                                       m.structure))
        @test all(m -> nfree(m) == m.parameters, ms)
        sel = brannelly2018_model_selection()
        # the three best models, in the printed order (Table 2)
        best3 = sort(sel; by = r -> r.aic)[1:3]
        @test [r.name for r in best3] == [m.name for m in ms[1:3]]
        @test all(r -> maximum(abs, r.fit.gradient) < 1e-2, sel)
        # deviances: printed 487.1455 (best) and 581.9903 (null); refit within 0.4%
        @test sel[1].deviance ≈ 487.1455 rtol = 0.004
        @test sel[10].deviance ≈ 581.9903 rtol = 0.001
        @test sel[1].deviance ≈ 488.8656 atol = 0.01      # regression
        # infected frogs more likely to be captured (Results, Fig. 3b)
        p = sel[1].fit.parameters.detection
        @test all(p[2, :] .> p[1, :])
    end

    @testset "model-averaged survival" begin
        s = brannelly2018_survival()
        v = parameter_set(:brannelly2018)
        @test round(s.uninfected; digits = 2) == round(v.survival_uninfected; digits = 2)
        @test round(s.infected; digits = 2) == round(v.survival_infected; digits = 2)
        @test sum(s.weights) ≈ 1
        @test akaike_weights([1.0, 1.0]) ≈ [0.5, 0.5]
        @test akaike_weights([0.0, 2.0])[1] ≈ 1 / (1 + exp(-1))
    end

    @testset "POPAN" begin
        # P(ω) summed over all 2^T − 1 non-empty histories plus P(0) is one
        T = 4
        H = reduce(vcat, [reverse(digits(k; base = 2, pad = T))' for k in 1:2^T-1])
        prob = PopanProblem(H; survival = :time, detection = :time, entry = :time,
                            intervals = [1.0, 2.0, 0.5])
        θ = [0.5, 1.0, -0.3, 0.1, -0.4, 0.3, 0.0, 0.2, -0.5, 0.7, log(10.0)]
        @test LogDensityProblems.dimension(prob) == length(θ)
        par = popan_parameters(prob, θ)
        @test sum(par.entry) ≈ 1
        @test par.N ≈ prob.u + 10
        # with one animal per history the log-likelihood is
        # lgamma(N+1) − lgamma(N−u+1) + Σ log P(ω) + (N − u) log P0; recover P(ω) by differencing
        lg = C.Distributions.SpecialFunctions.loggamma
        ll = LogDensityProblems.logdensity(prob, θ)
        # P0 from the closed form
        φ = par.survival .^ prob.intervals
        p = par.detection
        ψ = zeros(T); ψ[1] = par.entry[1]
        for t in 1:T-1
            ψ[t+1] = ψ[t] * (1 - p[t]) * φ[t] + par.entry[t+1]
        end
        P0 = 1 - sum(ψ .* p)
        ΣlogP = ll - lg(par.N + 1) + lg(par.N - prob.u + 1) - (par.N - prob.u) * log(P0)
        # brute force: enumerate entry time and survival paths
        function Pω(h)
            tot = 0.0
            for b in 1:T, d in b:T            # present from b, alive until d (dies after d)
                pr = par.entry[b] * prod(φ[b:d-1]; init = 1.0) * (d < T ? 1 - φ[d] : 1.0)
                for t in 1:T
                    present = b <= t <= d
                    pr *= present ? (h[t] == 1 ? p[t] : 1 - p[t]) : (h[t] == 1 ? 0.0 : 1.0)
                end
                tot += pr
            end
            return tot
        end
        @test ΣlogP ≈ sum(log(Pω(H[i, :])) for i in axes(H, 1))
        @test P0 ≈ Pω(zeros(Int, T))
        @test P0 + sum(Pω(H[i, :]) for i in axes(H, 1)) ≈ 1
        # simulated data: N recovered
        rng = Xoshiro(2018)
        Ntrue, Tn = 400, 6
        φt, pt = 0.8, 0.4
        pent = fill(1 / Tn, Tn)
        Hs = zeros(Int, Ntrue, Tn)
        for i in 1:Ntrue
            b = findfirst(cumsum(pent) .>= rand(rng))
            alive = true
            for t in b:Tn
                t > b && (alive &= rand(rng) < φt)
                alive || break
                Hs[i, t] = rand(rng) < pt
            end
        end
        fit = solve(PopanProblem(Hs; survival = :constant, detection = :constant,
                                 entry = :time))
        @test fit.N ≈ Ntrue rtol = 0.15
        @test fit.N_ci[1] < Ntrue < fit.N_ci[2]
        @test fit.parameters.survival[1] ≈ φt atol = 0.08
    end

    @testset "POPAN fit to the data" begin
        f = brannelly2018_popan()
        @test f.N ≈ 1394.08 rtol = 0.01                  # regression; printed 1854.33
        @test f.N > brannelly2018_summary().individuals
        @test maximum(abs, f.gradient) < 1e-2
    end
end
