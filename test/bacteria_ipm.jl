# Two-trait Bd–bacteria IPM (src/novel/bacteria_ipm.jl, NEW MODEL). There is no published
# version, so the oracles are structural:
#   - with no bacteria (and transmission :wilber2016) the Bd side is exactly the Wilber et al.
#     (2016) IPM, BdIPMProblem on the same mesh;
#   - with no coupling (every κ, ρ and r zero) the infected frogs per Bd load bin are the
#     Wilber IPM's whatever the bacteria do;
#   - with survival 1 and no Bd eviction the number of frogs is conserved;
#   - a bath of coverage f gives a fraction f × (bath mass above the floor) carrying bacteria;
#   - the protective effects act in the expected direction, and the tolerance effect raises R0
#     while the growth effect lowers it (tolerance amplification, cf. Barnett et al. 2026);
#   - R0 is zero without transmission and grows with the number of uninfected frogs for
#     density-dependent transmission.
using Chytridiomycosis
using Test
using LinearAlgebra: dot

@testset "Bd–bacteria IPM (new model)" begin
    C = Chytridiomycosis
    p = BacteriaIPMParameters()
    xmesh = IPMMesh(:wilber2016; n = 50)
    ymesh = IPMMesh(0.0, 14.0, 28)
    nx, ny = xmesh.n, ymesh.n
    T = 12.0

    @testset "mesh" begin
        m = IPMMesh2D(xmesh, ymesh)
        @test length(m) == nx * ny && size(m) == (nx, ny)
        xs, ys = meshpoints(m)
        @test xs == meshpoints(xmesh) && ys == meshpoints(ymesh)
        out = zeros(ny)
        below, above = C._normal_bin_masses!(out, 3.0, 2.0, ymesh)
        @test below + sum(out) + above ≈ 1 atol = 1e-12
    end

    @testset "no bacteria: the Wilber et al. (2016) IPM" begin
        prob = BacteriaIPMProblem(p; temperature = T, xmesh, ymesh, uninfected = 100,
            infected = 0, days = 60, transmission = :wilber2016)
        sol = solve(prob)
        ref = solve(BdIPMProblem(p.bd; temperature = T, mesh = xmesh, S0 = 100, days = 60))
        R = C._bipm_ranges(sol.mesh)
        @test sol.u[R.S0, :] ≈ ref.S rtol = 1e-12
        @test sol.u[R.I0, :] ≈ ref.I rtol = 1e-12
        @test all(iszero, sol.u[R.S, :]) && all(iszero, sol.u[R.I, :])
    end

    @testset "no coupling: Bd marginal is the Wilber IPM" begin
        q = C.with(p; κ_survival = 0.0, κ_growth = 0.0, κ_clearance = 0.0, ρ = 0.0)
        prob = BacteriaIPMProblem(q; temperature = T, xmesh, ymesh, uninfected = 100,
            days = 60, transmission = :wilber2016, baths = [(day = 0, coverage = 0.7),
                (day = 30, coverage = 0.4)])
        sol = solve(prob)
        ref = solve(BdIPMProblem(q.bd; temperature = T, mesh = xmesh, S0 = 100, days = 60))
        R = C._bipm_ranges(sol.mesh)
        Ix = sol.u[R.I0, :] .+ reshape(sum(reshape(sol.u[R.I, :], nx, ny, :); dims = 2), nx, :)
        @test Ix ≈ ref.I rtol = 1e-10
        @test sol.u[R.S0, :] .+ vec(sum(sol.u[R.S, :]; dims = 1)) ≈ ref.S rtol = 1e-10
        s = bacteria_ipm_summary(sol)
        @test s.carrying[1] == 0 && s.carrying[2] > 0.5       # bathed at the start
        @test s.carrying[end] < s.carrying[11]                 # the bacteria decline
    end

    @testset "conservation and baths" begin
        bd = C._with(p.bd; b00 = 60.0)                        # survival 1 at every load
        # Without Bd only the bacteria change, with exact bin integrals: frogs are conserved.
        nobd = C.with(p; bd, pool = C._with(p.pool; β0 = 0.0, β1 = 0.0),
            colonisation = 0.05, exchange = 0.2)
        s = bacteria_ipm_summary(solve(BacteriaIPMProblem(nobd; temperature = T, xmesh,
            ymesh, uninfected = 100, days = 90, baths = [(day = 0, coverage = 0.5)])))
        @test maximum(abs.(s.frogs .- 100)) < 1e-9
        @test s.carrying[2] > s.carrying[1] && all(s.prevalence .== 0)
        # With Bd, the midpoint rule on the Bd axis loses frogs (as in Wilber et al. 2016, who
        # use no eviction correction). Without coupling the loss is the published IPM's.
        q = C.with(p; bd, colonisation = 0.05, exchange = 0.2, κ_survival = 0.0,
            κ_growth = 0.0, ρ = 0.0)
        prob = BacteriaIPMProblem(q; temperature = T, xmesh, ymesh, uninfected = 90,
            infected = 10, days = 90, baths = [(day = 0, coverage = 0.5)],
            transmission = :wilber2016)
        s = bacteria_ipm_summary(solve(prob))
        ref = solve(BdIPMProblem(bd; temperature = T, mesh = xmesh, S0 = 90,
            I0 = prob.u0[C._bipm_ranges(prob.mesh).I0], days = 90))
        @test 100 - s.frogs[end] ≈ 100 - (ref.S[end] + sum(ref.I[:, end])) rtol = 1e-8
        @test all(>=(0), s.zoospores)
        u = copy(prob.u0)
        pc = C._bacteria_pieces(q, prob.mesh, T)
        C._bath!(u, 0.5, pc, prob.mesh)
        R = C._bipm_ranges(prob.mesh)
        carrying = sum(u[R.S]) + sum(u[R.I])
        @test carrying ≈ 0.5 * 100 * sum(pc.bath) rtol = 1e-12
        @test u[R.S0] + sum(u[R.S]) + sum(u[R.I0]) + sum(u[R.I]) ≈ 100 rtol = 1e-12
    end

    @testset "protection acts in the expected direction" begin
        run(q) = bacteria_ipm_summary(solve(BacteriaIPMProblem(q; temperature = T, xmesh,
            ymesh, uninfected = 98, infected = 2, days = 90,
            baths = [(day = 0, coverage = 1.0), (day = 30, coverage = 1.0)])))
        none = run(C.with(p; κ_survival = 0.0, κ_growth = 0.0, ρ = 0.0))
        tol = run(C.with(p; κ_survival = 3.0, κ_growth = 0.0, ρ = 0.0))
        res = run(C.with(p; κ_survival = 0.0, κ_growth = 2.0, ρ = 0.0))
        @test tol.frogs[end] > none.frogs[end]                 # tolerance: more survive
        @test res.mean_bd[11] < none.mean_bd[11]               # resistance: lower loads
        # Bd lowers the bacteria (ρ > 0): carriers decline faster
        hit = run(C.with(p; κ_survival = 0.0, κ_growth = 0.0, ρ = 0.1))
        @test hit.carrying[11] < none.carrying[11]
    end

    @testset "R0" begin
        pc = C._bacteria_pieces(p, IPMMesh2D(xmesh, ymesh), T)
        function prob_with(q; carrying = 0.0, A = 100.0)
            N = 1 + ny + nx + nx * ny + 1
            u0 = zeros(N)
            u0[1] = A * (1 - carrying)
            u0[2:1+ny] .= A * carrying .* pc.bath ./ sum(pc.bath)
            BacteriaIPMProblem(q; temperature = T, xmesh, ymesh, u0)
        end
        R0 = basic_reproduction_number(prob_with(p))
        @test R0 > 1
        nopool = C.with(p; pool = C._with(p.pool; β0 = 0.0, β1 = 0.0))
        @test basic_reproduction_number(prob_with(nopool)) == 0
        @test basic_reproduction_number(prob_with(p; A = 200.0)) > R0
        tol = C.with(p; κ_survival = 3.0, κ_growth = 0.0)
        res = C.with(p; κ_survival = 0.0, κ_growth = 2.0)
        @test basic_reproduction_number(prob_with(tol; carrying = 1.0)) >
              basic_reproduction_number(prob_with(tol))
        @test basic_reproduction_number(prob_with(res; carrying = 1.0)) <
              basic_reproduction_number(prob_with(res))
        @test_throws ArgumentError basic_reproduction_number(BacteriaIPMProblem(p;
            temperature = T, xmesh, ymesh, transmission = :wilber2016))
    end

    @testset "arguments" begin
        @test_throws ArgumentError BacteriaIPMParameters(r = 1.0)
        @test_throws ArgumentError BacteriaIPMParameters(κ_infection = 2.0)
        @test_throws ArgumentError BacteriaIPMProblem(p; temperature = T, transmission = :x)
        @test_throws ArgumentError BacteriaIPMProblem(p; temperature = T, days = 10)
        @test_throws ArgumentError BacteriaIPMProblem(p; temperature = T,
            baths = [(day = 0, coverage = 2.0)])
    end
end
