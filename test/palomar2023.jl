# Palomar et al. 2023 multi-event model (src/models/palomar2023.jl). The capture histories were
# not deposited (data/palomar2023/README.md): these tests check the model pattern, the
# compendium arithmetic and fitting on simulated histories.
using Random: Xoshiro
using LogDensityProblems: LogDensityProblems

@testset "Palomar et al. 2023" begin
    st = palomar2023_structure()
    @test st.nstates == 2 && size(st.assignment) == (2, 3)
    @test st.assignment == Bool[1 0 1; 0 1 1] && st.baseline == [3, 3]
    @test (st.survival, st.transition, st.detection, st.state_assignment) ==
          (:state_plus_time, :state_time, :time, :state)
    ms = palomar2023_models()
    @test length(ms) == 13 && ms[1].qaic == 4293.7 && issorted([m.qaic for m in ms])
    @test ms[8].structure.survival === :constant && ms[12].structure.transition === :state

    par = palomar2023_parameters()
    @test size(par.survival) == (2, 12) && par.survival[:, 1] == [0.754, 0.883]
    @test par.transition[:, :, 1] ≈ [0.775 0.225; 0.583 0.417]
    @test par.assignment ≈ [0.119 0 0.881; 0 1 0]
    # quasi-stationary prevalence by hand: left eigenvector of the alive block
    M = [0.754*0.775 0.754*0.225; 0.883*0.583 0.883*0.417]
    λ = (M[1, 1] + M[2, 2] + sqrt((M[1, 1] - M[2, 2])^2 + 4M[1, 2] * M[2, 1])) / 2
    @test palomar2023_implied_prevalence() ≈ 1 / (1 + (λ - M[1, 1]) / M[2, 1])
    @test palomar2023_implied_prevalence() ≈ 0.7136 atol = 1e-4

    # model 1 dimension: π 1 + S (2 + 11) + ψ (2 × 12) + p 12 + δ 2 = 52 for 13 occasions
    H = simulate_capture_histories(Xoshiro(5), st, par, fill(20, 13))
    prob = MultiEventProblem(H, st)
    @test LogDensityProblems.dimension(prob) == 52
    @test all(x -> x in 0:3, H)
    # with δ_U = 1 a seen uninfected toad is never "unknown"
    @test isfinite(LogDensityProblems.logdensity(prob, zeros(52)))

    # unassessed occasions: every seen animal gives the baseline (unknown) event there
    st0 = palomar2023_structure(; survival = :state, transition = :state, detection = :constant,
                                state_assignment = :constant)
    Hs = [0 3 1 0; 2 3 0 2; 3 0 3 1]
    p0 = MultiEventProblem(Hs, st0)
    p1 = MultiEventProblem(Hs, st0; unassessed = [2])
    θ = [0.2, 0.5, -0.3, 0.4, -1.0, 0.3, 0.1]
    @test LogDensityProblems.dimension(p0) == length(θ)
    @test LogDensityProblems.logdensity(p1, θ) > LogDensityProblems.logdensity(p0, θ)
    Hbad = [0 1 1 0]                  # a recorded status on an unassessed occasion is impossible
    @test LogDensityProblems.logdensity(MultiEventProblem(Hbad, st0; unassessed = [3]), θ) == -Inf

    # model selection runs on simulated histories (two cheap structures)
    sel = palomar2023_model_selection(simulate_capture_histories(Xoshiro(9), st0,
        (initial = [0.3, 0.7], survival = repeat([0.7, 0.85], 1, 4),
         transition = repeat([0.6 0.4; 0.3 0.7], 1, 1, 4), detection = fill(0.6, 2, 4),
         assignment = [0.5 0 0.5; 0 0.5 0.5]), fill(250, 5));
        models = [(model = 1, structure = st0, qaic = NaN),
                  (model = 2, structure = palomar2023_structure(; survival = :constant,
                       transition = :state, detection = :constant, state_assignment = :constant),
                   qaic = NaN)])
    @test length(sel) == 2 && sel[1].delta_qaic == 0 && sel[2].delta_qaic >= 0
    @test sel[1].qaic ≈ sel[1].fit.deviance / 1.474 + 2sel[1].fit.k
end
