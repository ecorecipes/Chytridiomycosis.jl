# Schmidt, Bozzuto, Lötters & Steinfartz (2017) R Soc Open Sci 4:160801.
# Bsal S–L–I model for adult female fire salamanders with logistic density dependence.
# Units: time in years, densities in adult females per hectare, space in km.
# Wiki: lab/wiki/models/bsal-schmidt2017-sli.md

"""
    SchmidtParameters(; b, d, d_I, e, β, K, D = 0.0625, μ = 1.0)
    SchmidtParameters(:schmidt2017)

Parameters of the Schmidt et al. (2017) model (eq. 2.1): recruitment `b`, natural death `d`,
Bsal death `d_I`, latent→infectious rate `e` (all yr⁻¹), mass-action transmission `β`
(ha yr⁻¹), carrying capacity `K` (females ha⁻¹), diffusion `D` (km² yr⁻¹) and scattered-colony
rate `μ` (yr⁻¹). The density-dependence strength is `γ = (b - d)/K`.
"""
Base.@kwdef struct SchmidtParameters{T<:Real}
    b::T
    d::T
    d_I::T
    e::T
    β::T
    K::T
    D::T = 0.0625
    μ::T = 1.0
end

function SchmidtParameters(key::Symbol)
    p = parameter_set(key)
    return SchmidtParameters(; b = p.b, d = p.d, d_I = p.d_I, e = p.e, β = p.β, K = p.K,
        D = p.D, μ = p.μ)
end

"""Return a copy of `p` with fields replaced, e.g. `with(p; K = 10.0)`."""
with(p::SchmidtParameters; kw...) =
    SchmidtParameters(; (f => get(kw, f, getfield(p, f)) for f in fieldnames(SchmidtParameters))...)

density_dependence(p::SchmidtParameters) = (p.b - p.d) / p.K

"""In-place right-hand side of eq. (2.1); state `u = [S, L, I]`."""
function schmidt_rhs!(du, u, p::SchmidtParameters, t)
    S, L, I = u
    N = S + L + I
    γ = density_dependence(p)
    du[1] = (p.b - p.d) * S - γ * S * N - p.β * S * I
    du[2] = p.β * S * I - (p.d + p.e) * L - γ * L * N
    du[3] = p.e * L - (p.d + p.d_I) * I - γ * I * N
    return nothing
end

"""
    invasion_threshold(p) -> K_threshold

Eq. (2.2): Bsal invades a host population at carrying capacity iff `K > K_threshold`,
with `K_threshold = (b + e)(b + d_I)/(β e)`.
"""
invasion_threshold(p::SchmidtParameters) = (p.b + p.e) * (p.b + p.d_I) / (p.β * p.e)

"""Proportion of the host population to remove to prevent an outbreak: `1 - K_threshold/K`."""
removal_fraction(p::SchmidtParameters) = max(zero(p.K), 1 - invasion_threshold(p) / p.K)

"""
    basic_reproduction_number(p) = K / K_threshold

R₀ of the linearised infected subsystem at the disease-free state `S = K` (derived here,
not reported in the paper).
"""
basic_reproduction_number(p::SchmidtParameters) = p.K / invasion_threshold(p)

"""
    endemic_equilibrium(p) -> (S, L, I)

The non-trivial equilibrium of eq. (2.1), eq. (2.3) of the paper:
`I = (b - d)(Kβe - (b + e)(b + d_I)) / (β(Kβe - b(b - d)))`, `L = (b + d_I - βI) I / e`,
`S = (b + e - βI)(b + d_I - βI)/(βe)`. (The `I²` terms cancel in `S + L + I = K_threshold - bI/e`.)
Returns `nothing` if `K ≤ K_threshold`.
"""
function endemic_equilibrium(p::SchmidtParameters)
    Kt = invasion_threshold(p)
    p.K > Kt || return nothing
    I = (p.b - p.d) * (p.K * p.β * p.e - (p.b + p.e) * (p.b + p.d_I)) /
        (p.β * (p.K * p.β * p.e - p.b * (p.b - p.d)))
    L = (p.b + p.d_I - p.β * I) * I / p.e
    S = (p.b + p.e - p.β * I) * (p.b + p.d_I - p.β * I) / (p.β * p.e)
    return (S = S, L = L, I = I)
end

"""Analytic Jacobian of eq. (2.1) at state `u = (S, L, I)`."""
function schmidt_jacobian(p::SchmidtParameters, u)
    S, L, I = u
    N = S + L + I
    γ = density_dependence(p)
    return [(p.b-p.d)-γ*N-γ*S-p.β*I   -γ*S                  -γ*S-p.β*S
            p.β*I-γ*L                 -(p.d+p.e)-γ*N-γ*L    p.β*S-γ*L
            -γ*I                      p.e-γ*I               -(p.d+p.d_I)-γ*N-γ*I]
end

"""Largest real part of the Jacobian's eigenvalues at the endemic equilibrium."""
function equilibrium_stability(p::SchmidtParameters)
    eq = endemic_equilibrium(p)
    eq === nothing && return nothing
    return maximum(real, eigvals(schmidt_jacobian(p, (eq.S, eq.L, eq.I))))
end

"""
    hopf_threshold(p; Kmax = 1000p.K) -> K_H

Carrying capacity at which the endemic equilibrium loses stability (Hopf bifurcation, with
`K` as the bifurcation parameter), found by bisection. The paper reports `K_H ≈ 9.5`.
"""
function hopf_threshold(p::SchmidtParameters; Kmax = 1000 * p.K, tol = 1e-10)
    lo = invasion_threshold(p) * (1 + 1e-6)
    hi = Kmax
    s(K) = equilibrium_stability(with(p; K = K))
    s(lo) < 0 || error("endemic equilibrium unstable just above the invasion threshold")
    s(hi) > 0 || return nothing
    while hi - lo > tol * hi
        mid = (lo + hi) / 2
        s(mid) < 0 ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end

"""
    wave_speed(p) -> c (km yr⁻¹)

Speed of the travelling wave of the reaction–diffusion model with equal diffusion `D` for
all classes (eq. 2.4): `c = 2√(Dλ)`, where
`λ = (√((e - d_I)² + 4eβK) - 2b - d_I - e)/2` is the growth rate of the linearised infected
subsystem at `S = K`.
"""
function wave_speed(p::SchmidtParameters)
    λ = (sqrt((p.e - p.d_I)^2 + 4 * p.e * p.β * p.K) - 2p.b - p.d_I - p.e) / 2
    return λ > 0 ? 2 * sqrt(p.D * λ) : zero(λ)
end

"""
    scattered_colony_area(p, t; c = wave_speed(p), μ = p.μ) -> km²

Eq. (2.5) as printed: total infected area `A(t) = (2πc²/μ)((e^{μt} - 1)/μ - t)`.
Note: the paper's reported proportions of German habitat (≈9% at μ = 1, ≈0.1% at μ = 0.1,
after 8 years, of 258 750 km²) equal this value divided by 100; see the wiki page.
"""
scattered_colony_area(p::SchmidtParameters, t; c = wave_speed(p), μ = p.μ) =
    2π * c^2 / μ * ((exp(μ * t) - 1) / μ - t)

"""
    SchmidtSLIProblem(p; u0 = [K, 0, 0.01], tspan = (0.0, 1.0))

The temporal model (eq. 2.1). The default `u0` is the paper's numerical experiment: a host
population at `K` with one infectious individual per km² (0.01 ha⁻¹).
"""
struct SchmidtSLIProblem{P<:SchmidtParameters,U,T}
    p::P
    u0::U
    tspan::T
end
SchmidtSLIProblem(p::SchmidtParameters; u0 = [p.K, zero(p.K), 0.01], tspan = (0.0, 1.0)) =
    SchmidtSLIProblem(p, collect(float.(u0)), float.(tspan))

to_ode_problem(prob::SchmidtSLIProblem) = ODEProblem(schmidt_rhs!, prob.u0, prob.tspan, prob.p)

CommonSolve.solve(prob::SchmidtSLIProblem, alg = Tsit5(); kwargs...) =
    solve(to_ode_problem(prob), alg; kwargs...)

"""
    population_crash(p; tspan = (0, 2), dt = 1e-4) -> (Smin, tmin)

Minimum susceptible density and when it is reached (the paper's figure 1 summary).
"""
function population_crash(p::SchmidtParameters; tspan = (0.0, 2.0), dt = 1e-4, kw...)
    sol = solve(SchmidtSLIProblem(p; tspan = tspan, kw...); abstol = 1e-10, reltol = 1e-8)
    ts = range(tspan..., step = dt)
    S = [sol(t)[1] for t in ts]
    i = argmin(S)
    return (Smin = S[i], tmin = ts[i])
end
