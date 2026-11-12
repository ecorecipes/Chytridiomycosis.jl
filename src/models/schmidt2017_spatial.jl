# Schmidt, Bozzuto, Lötters & Steinfartz (2017) R Soc Open Sci 4:160801, the two spatial
# versions of the Bsal S–L–I model.
#
# 1. Spatio-temporal formulation (§2.2): eq. (2.1) with equal diffusion D ∂²j/∂x² added to
#    j = S, L, I, in one spatial dimension. Discretised here by finite volumes (method of
#    lines) with zero-flux boundaries and lowered to an ODEProblem. Besides the paper's planar
#    1-D domain, `geometry = :radial` solves the same equations for a concentric (radially
#    symmetric) spread in the plane, the "concentric spread" the paper invokes to justify 1-D.
# 2. Scattered-colony model (§2.2, eq. 2.5; Shigesada & Kawasaki 1997): each infected colony
#    spreads as a disc of radius c·age, and founds new colonies elsewhere at rate μ (yr⁻¹).
#    Simulated exactly (pure-birth jump process; no JumpProcesses.jl needed), with colonies
#    placed uniformly at random in a square habitat so that overlaps and the habitat edge can
#    be counted, which eq. (2.5) ignores.
#
# Units: time in years, densities in adult females per hectare, space in km (D in km² yr⁻¹,
# c in km yr⁻¹, areas in km²). Wiki: lab/wiki/models/bsal-schmidt2017-sli.md

const SCHMIDT_GEOMETRIES = (:planar, :radial)

"""
    SchmidtPDEProblem(p; L = 8.0, n = 2000, geometry = :planar, tspan = (0.0, 0.5),
                      I0 = parameter_set(:schmidt2017).I0, x0 = L/n, u0 = nothing)

The spatio-temporal Schmidt et al. (2017) model (§2.2): eq. (2.1) with diffusion `p.D`
(km² yr⁻¹) on all three classes, on `[0, L]` km cut into `n` finite volumes with zero-flux
boundaries. `geometry = :planar` is the paper's one-dimensional domain (the infection starts
at the reflecting boundary `x = 0`, i.e. the middle of a symmetric domain); `geometry =
:radial` is concentric spread in the plane, `x` being the distance from the point of
introduction.

The default start is the paper's: hosts at `K` everywhere and `I0` infectious females ha⁻¹
(one per km², compendium record `I0`) within `x0` km of the origin. Pass `u0` (an `n × 3`
matrix with columns `S`, `L`, `I`) to override.

The front is very sharp at the published values (its leading edge decays as `e^{-qx}` with
`q = c/2D ≈ 89 km⁻¹`), so the grid spacing `L/n` should be a few metres; see
`front_speed` for checking the speed against [`wave_speed`](@ref). `solve` returns the
`ODESolution` (states are `n × 3` matrices); `to_ode_problem` gives the `ODEProblem`.
"""
struct SchmidtPDEProblem{P<:SchmidtParameters,T}
    p::P
    L::Float64
    n::Int
    geometry::Symbol
    u0::Matrix{Float64}
    tspan::T
end

function SchmidtPDEProblem(p::SchmidtParameters; L = 8.0, n::Integer = 2000,
        geometry::Symbol = :planar, tspan = (0.0, 0.5),
        I0 = parameter_set(:schmidt2017).I0, x0 = L / n, u0 = nothing)
    geometry in SCHMIDT_GEOMETRIES ||
        throw(ArgumentError("geometry must be one of $SCHMIDT_GEOMETRIES"))
    L > 0 && n >= 3 || throw(ArgumentError("need L > 0 and n ≥ 3"))
    if u0 === nothing
        x = cell_centres(L, n)
        u = zeros(n, 3)
        u[:, 1] .= p.K
        u[x .<= x0 + 1e-12, 3] .= I0
    else
        size(u0) == (n, 3) || throw(ArgumentError("u0 must be an n × 3 matrix (S, L, I)"))
        u = Matrix{Float64}(u0)
    end
    return SchmidtPDEProblem(p, float(L), Int(n), geometry, u, float.(tspan))
end

"""Cell centres (km) of the finite-volume grid of `n` cells on `[0, L]`."""
cell_centres(L, n) = ((1:n) .- 0.5) .* (L / n)
cell_centres(prob::SchmidtPDEProblem) = cell_centres(prob.L, prob.n)

# Face weights for the finite-volume Laplacian: flux through face i+1/2 is
# D wface[i] (u[i+1] - u[i]) / h, and du_i gets that flux divided by vol[i].
struct SchmidtGrid
    h::Float64
    wface::Vector{Float64}   # length n - 1 (interior faces)
    vol::Vector{Float64}     # length n
end

function schmidt_grid(L, n, geometry)
    h = L / n
    if geometry === :planar
        return SchmidtGrid(h, ones(n - 1), fill(h, n))
    else
        # radial: face i+1/2 at r = i h has weight r; the annulus of cell i has area
        # (r_{i+1/2}² - r_{i-1/2}²)/2 = h r_i (per radian)
        return SchmidtGrid(h, collect((1:n-1) .* h), cell_centres(L, n) .* h)
    end
end

struct SchmidtPDEParams{P<:SchmidtParameters}
    p::P
    grid::SchmidtGrid
end

function schmidt_pde_rhs!(du, u, q::SchmidtPDEParams, t)
    p, g = q.p, q.grid
    n = size(u, 1)
    γ = density_dependence(p)
    @inbounds for i in 1:n
        S, Li, I = u[i, 1], u[i, 2], u[i, 3]
        N = S + Li + I
        du[i, 1] = (p.b - p.d) * S - γ * S * N - p.β * S * I
        du[i, 2] = p.β * S * I - (p.d + p.e) * Li - γ * Li * N
        du[i, 3] = p.e * Li - (p.d + p.d_I) * I - γ * I * N
    end
    Dh = p.D / g.h
    @inbounds for j in 1:3, i in 1:n-1
        flux = Dh * g.wface[i] * (u[i+1, j] - u[i, j])
        du[i, j] += flux / g.vol[i]
        du[i+1, j] -= flux / g.vol[i+1]
    end
    return nothing
end

"""
    to_ode_problem(prob::SchmidtPDEProblem)

The method-of-lines `ODEProblem` (state: `n × 3` matrix of `S`, `L`, `I` per cell).
"""
to_ode_problem(prob::SchmidtPDEProblem) =
    ODEProblem(schmidt_pde_rhs!, copy(prob.u0), prob.tspan,
        SchmidtPDEParams(prob.p, schmidt_grid(prob.L, prob.n, prob.geometry)))

CommonSolve.solve(prob::SchmidtPDEProblem, alg = Tsit5(); kwargs...) =
    solve(to_ode_problem(prob), alg; kwargs...)

"""
    front_position(prob::SchmidtPDEProblem, u; level = 0.5) -> km

Position of the invasion front in state `u` (an `n × 3` matrix): the outermost point where
susceptibles fall below `level · K`, linearly interpolated between cell centres. Returns `0`
if nowhere has crashed yet.
"""
function front_position(prob::SchmidtPDEProblem, u; level = 0.5)
    x = cell_centres(prob)
    thr = level * prob.p.K
    S = view(u, :, 1)
    i = findlast(<(thr), S)
    i === nothing && return 0.0
    i == length(S) && return x[end]
    # S[i] < thr ≤ S[i+1]
    return x[i] + (thr - S[i]) / (S[i+1] - S[i]) * (x[i+1] - x[i])
end

"""
    front_speed(prob::SchmidtPDEProblem, sol; window = (tspan midpoint, tspan end),
                level = 0.5) -> km yr⁻¹

Least-squares slope of [`front_position`](@ref) against time over the saved times of `sol`
inside `window`; compare with the analytic travelling-wave speed [`wave_speed`](@ref)
(eq. 2.4). Solve with dense `saveat` inside the window.
"""
function front_speed(prob::SchmidtPDEProblem, sol;
        window = ((prob.tspan[1] + prob.tspan[2]) / 2, prob.tspan[2]), level = 0.5)
    idx = findall(t -> window[1] - 1e-12 <= t <= window[2] + 1e-12, sol.t)
    length(idx) >= 2 || throw(ArgumentError("need at least two saved times in the window"))
    ts = sol.t[idx]
    xs = [front_position(prob, sol.u[k]; level) for k in idx]
    tm, xm = sum(ts) / length(ts), sum(xs) / length(xs)
    return sum((ts .- tm) .* (xs .- xm)) / sum((ts .- tm) .^ 2)
end

"""
    discrete_wave_speed(p, h)

Linear spreading speed of the semi-discrete (finite-volume, spacing `h` km) planar model:
`min_q (λ(q) + 2D(cosh(qh) - 1)/h²)/q` with `λ` the infected-subsystem growth rate at
`S = K`. Tends to [`wave_speed`](@ref) as `h → 0`; used to separate discretisation error from
model error in tests.
"""
function discrete_wave_speed(p::SchmidtParameters, h)
    λ = (sqrt((p.e - p.d_I)^2 + 4 * p.e * p.β * p.K) - 2p.b - p.d_I - p.e) / 2
    λ > 0 || return zero(λ)
    s(q) = (λ + 2p.D * (cosh(q * h) - 1) / h^2) / q
    # golden-section search around the continuum minimiser q* = √(λ/D)
    a, b = 0.05 * sqrt(λ / p.D), 5 * sqrt(λ / p.D)
    φ = (sqrt(5) - 1) / 2
    for _ in 1:200
        c1, c2 = b - φ * (b - a), a + φ * (b - a)
        s(c1) < s(c2) ? (b = c2) : (a = c1)
    end
    return s((a + b) / 2)
end

# ---------------------------------------------------------------------------------------
# Scattered-colony (long-distance jump) model

"""
    SchmidtColonyProblem(p; μ = p.μ, c = wave_speed(p), tspan = (0.0, 8.0),
                         area = parameter_set(:schmidt2017).habitat_area,
                         origin = :centre, max_colonies = 10^6)

The scattered-colony model behind eq. (2.5) of Schmidt et al. (2017), as an explicit jump
process. At time 0 one infected colony exists; every colony spreads locally as a disc whose
radius grows at `c` km yr⁻¹ (eq. 2.4), and each colony independently sends infected
salamanders that found a new colony at rate `μ` yr⁻¹. New colonies land uniformly at random
in a square habitat of `area` km² (default: the paper's 258 750 km² of fire-salamander
habitat in Germany). `origin` is `:centre`, `:random` or an `(x, y)` position in km.

The expected number of colonies is `e^{μt}` and the expected summed disc area is exactly
eq. (2.5); the simulated union of discs inside the habitat is the area actually infected,
which eq. (2.5) overestimates once colonies overlap or reach the habitat edge.

`solve(prob; rng, saveat, resolution = 1.0)` (`saveat` defaults to nine equally spaced times
over `tspan`, i.e. years 0, 1, …, 8) returns a NamedTuple with the founding
`times` and `positions` (km) of all colonies, and at each time in `saveat`: `colonies`,
`summed_area` (km², the no-overlap area of eq. 2.5), `infected_area` (km², the union of discs
clipped to the habitat, rasterised at `resolution` km; `NaN` when `resolution = nothing`)
and `fraction = infected_area / area`.
"""
struct SchmidtColonyProblem{T}
    μ::Float64
    c::Float64
    tspan::T
    area::Float64
    origin::Union{Symbol,NTuple{2,Float64}}
    max_colonies::Int
end

function SchmidtColonyProblem(p::SchmidtParameters; μ = p.μ, c = wave_speed(p),
        tspan = (0.0, 8.0), area = parameter_set(:schmidt2017).habitat_area,
        origin = :centre, max_colonies::Integer = 10^6)
    μ >= 0 && c >= 0 && area > 0 || throw(ArgumentError("need μ ≥ 0, c ≥ 0, area > 0"))
    o = origin isa Symbol ? origin : (Float64(origin[1]), Float64(origin[2]))
    o isa Symbol && !(o in (:centre, :random)) &&
        throw(ArgumentError("origin must be :centre, :random or an (x, y) position"))
    return SchmidtColonyProblem(float(μ), float(c), float.(tspan), float(area), o,
        Int(max_colonies))
end

function CommonSolve.solve(prob::SchmidtColonyProblem; rng::AbstractRNG = Random.default_rng(),
        saveat = range(prob.tspan[1], prob.tspan[2]; length = 9), resolution = 1.0)
    t0, T = prob.tspan
    side = sqrt(prob.area)
    start = prob.origin === :centre ? (side / 2, side / 2) :
            prob.origin === :random ? (side * rand(rng), side * rand(rng)) : prob.origin
    times = [t0]
    xs = [start[1]]
    ys = [start[2]]
    # Pure-birth (Yule) process: with n colonies the next founding comes after Exp(μn).
    t = t0
    while prob.μ > 0
        t += Random.randexp(rng) / (prob.μ * length(times))
        t > T && break
        length(times) >= prob.max_colonies &&
            error("more than max_colonies = $(prob.max_colonies) colonies by t = $t")
        push!(times, t)
        push!(xs, side * rand(rng))
        push!(ys, side * rand(rng))
    end
    ts = collect(float.(saveat))
    colonies = [count(<=(s), times) for s in ts]
    summed = [sum(π * (prob.c * (s - ti))^2 for ti in times if ti <= s; init = 0.0) for s in ts]
    infected = resolution === nothing ? fill(NaN, length(ts)) :
        colony_union_area(times, xs, ys, prob.c, side, float(resolution), ts)
    return (times = times, positions = hcat(xs, ys), saveat = ts, colonies = colonies,
        summed_area = summed, infected_area = infected, fraction = infected ./ prob.area)
end

# Union of growing discs inside the square [0, side]², on a raster of cell size `res`:
# each cell's infection time is min_i (t_i + |x - x_i|/c); the area at time s is the number
# of cells infected by s times res².
function colony_union_area(times, xs, ys, c, side, res, ts)
    m = max(1, round(Int, side / res))
    h = side / m
    T = maximum(ts)
    arrival = fill(Inf, m, m)
    for k in eachindex(times)
        R = c * (T - times[k])
        R >= 0 || continue
        ilo = max(1, floor(Int, (xs[k] - R) / h) + 1)
        ihi = min(m, floor(Int, (xs[k] + R) / h) + 1)
        jlo = max(1, floor(Int, (ys[k] - R) / h) + 1)
        jhi = min(m, floor(Int, (ys[k] + R) / h) + 1)
        @inbounds for j in jlo:jhi, i in ilo:ihi
            dx = (i - 0.5) * h - xs[k]
            dy = (j - 0.5) * h - ys[k]
            a = times[k] + sqrt(dx^2 + dy^2) / max(c, eps())
            a < arrival[i, j] && (arrival[i, j] = a)
        end
    end
    return [count(<=(s), arrival) * h^2 for s in ts]
end
