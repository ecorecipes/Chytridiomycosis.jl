# Generic HMM layer (src/observation/hmm.jl) against brute-force enumeration of all paths.
using Random: Xoshiro
using LogDensityProblems: LogDensityProblems

@testset "Generic HMM layer" begin
    rowstoch(rng, K) = (A = rand(rng, K, K) .+ 0.05; A ./ sum(A; dims = 2))
    function brute(init, Ps, L; emit_first = true)
        nT, nK = size(L)
        paths = Iterators.product(ntuple(_ -> 1:nK, nT)...)
        joint = Dict{NTuple{nT,Int},Float64}()
        for x in paths
            p = init[x[1]] * (emit_first ? L[1, x[1]] : 1.0)
            for t in 2:nT
                p *= Ps[t-1][x[t-1], x[t]] * L[t, x[t]]
            end
            joint[x] = p
        end
        return joint
    end

    rng = Xoshiro(7)
    K, T = 3, 5
    init = rand(rng, K); init ./= sum(init)
    Ps = [rowstoch(rng, K) for _ in 1:T-1]
    L = rand(rng, T, K)
    L[3, :] .= 1.0                                   # a missing observation

    for emit_first in (true, false)
        joint = brute(init, Ps, L; emit_first)
        Z = sum(values(joint))
        @test hmm_loglikelihood(init, Ps, L; emit_first) ≈ log(Z)
        γ = hmm_posterior(init, Ps, L; emit_first)
        γb = zeros(T, K)
        for (x, p) in joint, t in 1:T
            γb[t, x[t]] += p / Z
        end
        @test γ ≈ γb
        @test all(sum(γ; dims = 2) .≈ 1)
        v = hmm_viterbi(init, Ps, L; emit_first)
        best = argmax(joint)
        @test v.path == collect(best)
        @test v.logprob ≈ log(joint[best])
        f = hmm_forward(init, Ps, L; emit_first)
        @test sum(f.logc) ≈ f.loglik
        @test f.filtered[T, :] ≈ γ[T, :]                   # filtering = smoothing at T
        β = hmm_backward(Ps, L, f.logc)
        @test size(β) == (T, K) && all(β[T, :] .== 1)
    end

    # time-homogeneous transitions give the same as a repeated vector
    P = rowstoch(rng, K)
    @test hmm_loglikelihood(init, P, L) ≈ hmm_loglikelihood(init, fill(P, T - 1), L)

    # long sequences do not underflow
    Llong = fill(1e-3, 3000, K)
    @test isfinite(hmm_loglikelihood(init, P, Llong))
    @test hmm_loglikelihood(init, P, Llong) ≈ 3000 * log(1e-3)
    # a zero-probability sequence
    Lz = copy(L); Lz[2, :] .= 0
    @test hmm_loglikelihood(init, Ps, Lz) == -Inf
    @test_throws DomainError hmm_posterior(init, Ps, Lz)
    @test_throws DimensionMismatch hmm_forward(init[1:2], Ps, L)
    @test_throws DimensionMismatch hmm_forward(init, Ps[1:2], L)

    # categorical emissions, with missing observations
    B = [0.9 0.1; 0.2 0.8; 0.5 0.5]
    e = CategoricalEmission(B)
    Le = emission_likelihoods(e, [1, missing, 2, 2, 1])
    @test Le[2, :] == ones(3)
    @test Le[3, :] == B[:, 2]
    joint = brute(init, Ps, Le)
    @test hmm_loglikelihood(init, Ps, Le) ≈ log(sum(values(joint)))

    # covariate-dependent transitions and the LogDensityProblems wrapper
    two(θ, x) = (p = 1 / (1 + exp(-(θ[1] + θ[2] * x[1]))); [1-p p; 0.3 0.7])
    X = reshape(collect(range(-1, 1; length = 4)), 4, 1)
    ct = CovariateTransitions(two, X)
    θ = [0.2, -0.5]
    @test length(ct(θ)) == 4 && all(P -> all(sum(P; dims = 2) .≈ 1), ct(θ))
    obs = [1, 2, 2, 1, 2]
    E2 = CategoricalEmission([0.8 0.2; 0.3 0.7])
    seqs(θ) = [(initial = [0.5, 0.5], transitions = ct(θ), likelihoods = emission_likelihoods(E2, obs)),
               ([0.9, 0.1], ct(θ), emission_likelihoods(E2, reverse(obs)))]
    m = HMMLogDensity(seqs, θ -> -sum(abs2, θ) / 2, 2)
    @test LogDensityProblems.dimension(m) == 2
    ll = Chytridiomycosis.pointwise_loglikelihood(m, θ)
    @test ll[1] ≈ log(sum(values(brute([0.5, 0.5], ct(θ), emission_likelihoods(E2, obs)))))
    @test LogDensityProblems.logdensity(m, θ) ≈ sum(ll) - sum(abs2, θ) / 2

    # sampled paths follow the chain
    paths = [hmm_sample(Xoshiro(k), init, P, 2) for k in 1:20_000]
    freq = [count(p -> p[2] == j, paths) for j in 1:K] ./ length(paths)
    @test freq ≈ vec(init' * P) atol = 0.015
    @test all(p -> length(p) == 2, paths)
end
