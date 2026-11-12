# Schmidt et al. 2017, R Soc Open Sci 4:160801: the spatio-temporal (§2.2, eq. 2.4) and
# scattered-colony (eq. 2.5, ESM fig. S3) versions. Checked against the PDF
# (lab/papers/schmidt2017.pdf) and ESM (lab/papers/schmidt2017_si.pdf).
using Chytridiomycosis
using Test
using Random: Xoshiro

@testset "Schmidt 2017 spatial and jump versions" begin
    C = Chytridiomycosis
    p = SchmidtParameters(:schmidt2017)
    ps = parameter_set(:schmidt2017)
    @test ps.I0 == 0.01 && ps.habitat_area == 258_750.0
    @test parameter_set(:schmidt2017; set = "figS3_low").μ == 0.1

    @testset "eq. 2.4 as printed and ESM fig. S1c-f" begin
        # c = √(2D(√((e - d_I)² + 4eβK) - 2b - d_I - e)), exactly as printed in the PDF
        c24 = sqrt(2p.D * (sqrt((p.e - p.d_I)^2 + 4p.e * p.β * p.K) - 2p.b - p.d_I - p.e))
        @test wave_speed(p) ≈ c24
        # Fig. S1e: A(t = 1) ≈ 560 km² at the point estimates: eq. 2.5 is in km²
        @test 555 < scattered_colony_area(p, 1.0) < 565
        # Fig. S1c-f: +20 % in D gives ≈ +9.5 % speed and +20 % area; +20 % in K ≈ +13 % area
        up(; kw...) = C.with(p; kw...)
        @test wave_speed(up(D = 1.2p.D)) / wave_speed(p) ≈ sqrt(1.2)
        @test scattered_colony_area(up(D = 1.2p.D), 1.0) / scattered_colony_area(p, 1.0) ≈ 1.2
        @test 1.12 < scattered_colony_area(up(K = 1.2p.K), 1.0) / scattered_colony_area(p, 1.0) < 1.14
    end

    @testset "reaction-diffusion front (method of lines)" begin
        # planar 1-D domain, Table 1 values: grid 4 m, front measured over 0.25-0.45 yr
        prob = SchmidtPDEProblem(p; L = 6.0, n = 1500, tspan = (0.0, 0.45))
        @test size(prob.u0) == (1500, 3) && maximum(prob.u0[:, 3]) == ps.I0
        ode = to_ode_problem(prob)
        @test ode isa C.SciMLBase.ODEProblem
        sol = solve(prob; saveat = 0.0:0.01:0.45, abstol = 1e-10, reltol = 1e-8)
        @test C.SciMLBase.successful_retcode(sol)
        c = front_speed(prob, sol; window = (0.25, 0.45))
        @test isapprox(c, wave_speed(p); rtol = 0.01)                # 11.15 km/yr
        @test isapprox(c, C.discrete_wave_speed(p, prob.L / prob.n); rtol = 0.01)
        @test C.discrete_wave_speed(p, 1e-4) ≈ wave_speed(p) rtol = 1e-4
        # behind the front the population has crashed ("almost zero", figure 1)
        u = sol.u[end]
        x = C.cell_centres(prob)
        xf = front_position(prob, u)
        @test 4.0 < xf < 5.5
        @test maximum(u[x .< xf - 1.0, 1]) < 0.05p.K
        @test all(u[x .> xf + 0.2, 1] .> 0.99p.K)                   # untouched ahead of it
        # concentric spread in the plane (radial geometry) moves at the same speed
        rad = SchmidtPDEProblem(p; L = 6.0, n = 1500, geometry = :radial, tspan = (0.0, 0.45))
        rsol = solve(rad; saveat = 0.0:0.01:0.45, abstol = 1e-10, reltol = 1e-8)
        @test isapprox(front_speed(rad, rsol; window = (0.25, 0.45)), wave_speed(p); rtol = 0.01)
        # lower carrying capacity, coarser grid: slower front, still eq. 2.4
        p5 = C.with(p; K = 5.0)
        pr5 = SchmidtPDEProblem(p5; L = 6.0, n = 600, tspan = (0.0, 0.8))
        s5 = solve(pr5; saveat = 0.0:0.02:0.8)
        @test isapprox(front_speed(pr5, s5; window = (0.4, 0.8)), wave_speed(p5); rtol = 0.01)
        # below K_threshold Bsal does not spread (eq. 2.2 also holds with diffusion)
        plow = C.with(p; K = 0.9)
        @test plow.K < invasion_threshold(plow) && wave_speed(plow) == 0
        prl = SchmidtPDEProblem(plow; L = 2.0, n = 200, tspan = (0.0, 1.0))
        sl = solve(prl)
        @test front_position(prl, sl.u[end]) == 0
        @test sum(sl.u[end][:, 2:3]) < 0.05 * sum(prl.u0[:, 3])      # infection dies out
        # the finite-volume Laplacian conserves hosts (no births, deaths or infection)
        pc = C.with(p; b = 0.1, d = 0.1, β = 0.0)
        for geometry in (:planar, :radial)
            u0 = zeros(100, 3)
            u0[1:10, 1] .= 1.0
            pc_prob = SchmidtPDEProblem(pc; L = 1.0, n = 100, geometry, u0, tspan = (0.0, 0.5))
            g = C.schmidt_grid(1.0, 100, geometry)
            s = solve(pc_prob; abstol = 1e-12, reltol = 1e-10)
            @test sum(s.u[end][:, 1] .* g.vol) ≈ sum(u0[:, 1] .* g.vol) rtol = 1e-8
        end
        @test_throws ArgumentError SchmidtPDEProblem(p; geometry = :spherical)
        @test_throws ArgumentError SchmidtPDEProblem(p; n = 10, u0 = zeros(5, 3))
    end

    @testset "scattered-colony jump process (eq. 2.5)" begin
        # μ = 0: one colony, a disc of radius c t (rasterised at 0.5 km)
        one = solve(SchmidtColonyProblem(p; μ = 0.0); saveat = [8.0], resolution = 0.5)
        @test one.colonies == [1]
        @test one.summed_area[1] ≈ π * (8wave_speed(p))^2
        @test isapprox(one.infected_area[1], one.summed_area[1]; rtol = 0.01)
        # E[colonies] = e^{μt} and E[summed disc area] = eq. 2.5 exactly
        cp = SchmidtColonyProblem(p; μ = 1.0, tspan = (0.0, 3.0))
        rs = [solve(cp; rng = Xoshiro(i), saveat = [3.0], resolution = nothing) for i in 1:2000]
        @test isapprox(sum(r -> r.colonies[1], rs) / 2000, exp(3); rtol = 0.08)
        @test isapprox(sum(r -> r.summed_area[1], rs) / 2000, scattered_colony_area(p, 3.0);
            rtol = 0.08)
        @test all(r -> isnan(r.infected_area[1]), rs)
        # ESM fig. S3 plots A(t)/habitat: 8.97 at μ = 1 and 0.128 at μ = 0.1 after 8 years,
        # reported as "≈9 %" and "≈0.1 %". As fractions, eq. 2.5 exceeds the habitat 9-fold.
        H = ps.habitat_area
        @test round(scattered_colony_area(p, 8.0) / H, digits = 2) == 8.97
        @test round(scattered_colony_area(p, 8.0; μ = 0.1) / H, digits = 3) == 0.128
        # The explicit spatial process: overlaps and the habitat edge cap the infected area.
        big = SchmidtColonyProblem(p; μ = 1.0)
        fr = [solve(big; rng = Xoshiro(i), saveat = [8.0], resolution = 4.0).fraction[1]
              for i in 1:20]
        @test all(<=(1.0), fr)
        @test sum(fr) / 20 > 0.75                                   # far above "≈9 %"
        r1 = solve(big; rng = Xoshiro(3), saveat = 0:8, resolution = 4.0)
        @test all(r1.infected_area .<= r1.summed_area .* 1.02)
        @test issorted(r1.infected_area) && issorted(r1.colonies)
        @test r1.colonies[end] == length(r1.times) && size(r1.positions) == (length(r1.times), 2)
        @test all(0 .<= r1.positions .<= sqrt(H))
        @test_throws ArgumentError SchmidtColonyProblem(p; origin = :nowhere)
        @test_throws ErrorException solve(SchmidtColonyProblem(p; max_colonies = 5);
            rng = Xoshiro(1), resolution = nothing)
    end
end
