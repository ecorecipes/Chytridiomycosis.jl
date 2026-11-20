# Espira, Brouwer, Han, Foufopoulos & Eisenberg (2022) Am Nat 199(2):E43-E56,
# doi:10.1086/717413. Two host species with susceptible, infectious and recovered classes;
# each species sheds into its own habitat, and the habitats overlap. Checked against the
# published version (lab/papers/espira2022/espira2022.pdf), the author manuscript
# (PMC9136953; its equations 1-11 are identical to the published ones) and the online
# supplement (espira2022_si.pdf, identical in text to the manuscript supplement NIHMS1800855).
#
# Model (eq. 1 = eq. S5), i = 1, 2 and j the other species:
#   dS_i/dt = -κ_i ρ_i π_i S_i (E_i + θ_i E_j)
#   dI_i/dt =  κ_i ρ_i π_i S_i (E_i + θ_i E_j) - γ_i I_i
#   dR_i/dt =  γ_i I_i
#   dE_i/dt =  α_i I_i/ε_i - E_i (κ_i ρ_i N_i/ε_i + θ_j κ_j ρ_j N_j/ε_i + μ)
# κ contact rate with the environment (day⁻¹), ρ area per contact (pick-up), π infection
# probability per pathogen, α shedding (pathogens day⁻¹), γ recovery, μ pathogen decay, ε_i the
# area of habitat i, θ_i = area(E_i ∩ E_j)/area(E_i). E_i is a concentration (per area). No
# births, deaths, migration or pathogen growth.
#
# Reproduction numbers:
#   η_i    = (κ_i ρ_i N_i/ε_i)/(κ_i ρ_i N_i/ε_i + μ)                       (eq. 4)
#   R0,i   = (π_i α_i/γ_i) η_i                                              (eq. 5)
#   R0,C   = spectral radius of the next-generation matrix of eq. 1          (eq. 6, S15)
#   dilution (R0,C < R0,1) iff η_2 > [(a_2 + b_12 + θ_2 μ)/(a_2 + b_12 + μ)] R0,2/R0,1   (eq. 8)
# with a_2 = κ_2ρ_2N_2/ε_2 and b_12 = θ_1κ_1ρ_1N_1/ε_2.
#
# Next-generation matrix (supplement eqs S6-S12): F and V on [I1, I2, E1, E2] as here; the
# host block of F V⁻¹ (S12) is
#   K_ii = α_i π_i κ_i ρ_i N_i / (γ_i (κ_i ρ_i N_i + θ_j κ_j ρ_j N_j + μ ε_i)),
#   K_ij = α_j θ_i π_i κ_i ρ_i N_i / (γ_j (κ_j ρ_j N_j + θ_i κ_i ρ_i N_i + μ ε_j)),
# so det K = K_11 K_22 (1 - θ_1 θ_2) (S14) and R0,C is S15. This is what
# `community_reproduction_number` computes.
#
# Eq. 6 as printed (in the published version and in the manuscript, and eq. 7 in the same
# way) has denominators γ_1(κ_1ρ_1N_1 + θ_1κ_2ρ_2N_2 + μ): μ without ε_1, and θ_1 where S15
# has θ_2. The supplement's own S15 and S16 have μ ε_i and θ_j, so eqs 6-7 are misprints of
# S15-S16; they coincide when ε_1 = ε_2 = 1 (then θ_1 = θ_2). `as_printed = true` evaluates
# eq. 6 as printed.
#
# Figures S1 and S2 (the supplement's "figs. A1, A2"): the embedded heat-map rasters of S1, S2
# and the published Fig. 4 are the same image (S1 and S2 byte-identical; Fig. 4 differs only by
# its CMYK-to-RGB conversion, mean absolute difference 1.8/255). Under eq. 1, R0,C/R0,1 on the
# figure axes depends on θ_1 and θ_2 separately (θ_1 ε_1 = θ_2 ε_2 gives b_21 = θ_1 a_2 and
# b_12 = θ_2 a_1), so S1 (θ_1 = 1) and S2 (θ_2 = 1) cannot equal Fig. 4 (θ_1 = θ_2); see the
# targets espira2022_figS1 and espira2022_figS2.

"""
    EspiraParameters(; N, γ, κ, ρ, π, α, ε, θ, μ)

Parameters of the two-species environmental transmission model of Espira et al. (2022),
eq. 1. Every argument but `μ` is a pair, one value per species:

- `N`: population size;
- `γ`: recovery rate (day⁻¹);
- `κ`: contact rate with the environment (contacts day⁻¹);
- `ρ`: pick-up, the area covered per contact times the fraction of pathogens picked up there;
- `π`: probability of infection per pathogen picked up;
- `α`: shedding rate (pathogens day⁻¹ per infectious host);
- `ε`: area of the species' habitat;
- `θ`: habitat overlap, `θ_i = area(E_i ∩ E_j)/area(E_i)`, in `[0, 1]`. By this definition
  `θ_1 ε_1 = θ_2 ε_2`; the constructor does not enforce it, since eq. 1 does not need it.
- `μ`: pathogen decay rate (day⁻¹), the same in both habitats.

The paper is analytical and its figures use no printed parameter values; the only numbers are
the worked example of the supplement (Table S1, in the compendium under `:espira2022`), which
has no recovery or shedding rates. There is therefore no keyed constructor; see
[`espira_parameters_from_efficiencies`](@ref) for the axes of Figs 3 and 4.
"""
struct EspiraParameters
    N::NTuple{2,Float64}
    γ::NTuple{2,Float64}
    κ::NTuple{2,Float64}
    ρ::NTuple{2,Float64}
    π::NTuple{2,Float64}
    α::NTuple{2,Float64}
    ε::NTuple{2,Float64}
    θ::NTuple{2,Float64}
    μ::Float64
    function EspiraParameters(N, γ, κ, ρ, π, α, ε, θ, μ)
        pr(x) = (Float64(x[1]), Float64(x[2]))
        v = map(pr, (N, γ, κ, ρ, π, α, ε, θ))
        all(x -> all(>=(0), x), v) && μ >= 0 ||
            throw(ArgumentError("parameters must be non-negative"))
        all(>(0), v[2]) && all(>(0), v[7]) ||
            throw(ArgumentError("γ and ε must be positive"))
        all(<=(1), v[8]) || throw(ArgumentError("θ must lie in [0, 1]"))
        return new(v..., Float64(μ))
    end
end
EspiraParameters(; N, γ, κ, ρ, π, α, ε, θ, μ) = EspiraParameters(N, γ, κ, ρ, π, α, ε, θ, μ)

time_unit(::Type{<:EspiraParameters}) = :day

_other(i) = 3 - i
# pick-up rate of species i's pathogens by species i itself, κ_i ρ_i N_i/ε_i
_self_pickup(p::EspiraParameters, i) = p.κ[i] * p.ρ[i] * p.N[i] / p.ε[i]
# loss rate of pathogens in habitat i: pick-up by both species and decay
_loss(p::EspiraParameters, i) = (j = _other(i);
    _self_pickup(p, i) + p.θ[j] * p.κ[j] * p.ρ[j] * p.N[j] / p.ε[i] + p.μ)

"""
    pickup_efficiency(p::EspiraParameters) -> (η₁, η₂)

The pathogen pick-up efficiency of each species (Espira et al. 2022, eq. 4),
`η_i = (κ_i ρ_i N_i/ε_i)/(κ_i ρ_i N_i/ε_i + μ)`: the fraction of the pathogens a species
sheds that it picks up itself before they decay, when alone. Near 0 transmission is
density-dependent-like, near 1 frequency-dependent-like.
"""
pickup_efficiency(p::EspiraParameters) =
    Tuple(_self_pickup(p, i) / (_self_pickup(p, i) + p.μ) for i in 1:2)

"""
    species_reproduction_numbers(p::EspiraParameters) -> [R0,1, R0,2]

The single-species basic reproduction numbers of Espira et al. (2022), eq. 5:
`R0,i = (π_i α_i/γ_i) η_i`, the reproduction number of species `i` alone (`N_j = 0`).
"""
function species_reproduction_numbers(p::EspiraParameters)
    η = pickup_efficiency(p)
    return [p.π[i] * p.α[i] / p.γ[i] * η[i] for i in 1:2]
end

"""
    next_generation_matrix(p::EspiraParameters; reduced = true)

The next-generation matrix of Espira et al. (2022), eq. 1, at the disease-free state
(`S_i = N_i`). With `reduced = true` it is the 2 × 2 host-to-host matrix `K`, where `K[j, i]`
is the number of species-`j` infections caused, through the environment, by one
species-`i` infection. With `reduced = false` it is `F V⁻¹` on the infected states
`[I₁, I₂, E₁, E₂]`. Both have spectral radius
[`community_reproduction_number`](@ref).
"""
function next_generation_matrix(p::EspiraParameters; reduced::Bool = true)
    # F: new infections; V: transitions, on [I1, I2, E1, E2]
    F = zeros(4, 4)
    V = zeros(4, 4)
    for i in 1:2
        j = _other(i)
        c = p.κ[i] * p.ρ[i] * p.π[i] * p.N[i]
        F[i, 2 + i] = c
        F[i, 2 + j] = c * p.θ[i]
        V[i, i] = p.γ[i]
        V[2 + i, i] = -p.α[i] / p.ε[i]
        V[2 + i, 2 + i] = _loss(p, i)
    end
    G = F / V
    reduced || return G
    return G[1:2, 1:2]   # rows 3:4 of F V⁻¹ are zero, so this block has the same spectral radius
end

"""
    community_reproduction_number(p::EspiraParameters; as_printed = false)

The community basic reproduction number `R0,C` of Espira et al. (2022): the spectral radius
of the next-generation matrix of eq. 1 ([`next_generation_matrix`](@ref)), in closed form

    R0,C = (K₁₁ + K₂₂)/2 + sqrt((K₁₁ + K₂₂)² - 4 K₁₁ K₂₂ (1 - θ₁θ₂))/2,
    K_ii = (π_i α_i/γ_i) (κ_iρ_iN_i/ε_i)/(κ_iρ_iN_i/ε_i + θ_jκ_jρ_jN_j/ε_i + μ).

This is eq. S15 of the supplement. `as_printed = true` evaluates eq. 6 of the main text
exactly as printed, whose `K_ii` denominators are `κ_iρ_iN_i + θ_iκ_jρ_jN_j + μ`: `μ` without
`ε_i`, and the other overlap parameter. S15 shows these are misprints; the two agree when
`ε₁ = ε₂ = 1` and `θ₁ = θ₂`.
"""
function community_reproduction_number(p::EspiraParameters; as_printed::Bool = false)
    K = zeros(2)
    for i in 1:2
        j = _other(i)
        a = p.κ[i] * p.ρ[i] * p.N[i]
        den = as_printed ? a + p.θ[i] * p.κ[j] * p.ρ[j] * p.N[j] + p.μ : _loss(p, i) * p.ε[i]
        K[i] = p.α[i] * p.π[i] * a / (p.γ[i] * den)
    end
    tr = K[1] + K[2]
    return (tr + sqrt(max(tr^2 - 4 * K[1] * K[2] * (1 - p.θ[1] * p.θ[2]), 0.0))) / 2
end

"""
    dilution_threshold(p::EspiraParameters)

The right-hand side of the dilution criterion of Espira et al. (2022), eq. 8:

    [(κ₂ρ₂N₂/ε₂ + θ₁κ₁ρ₁N₁/ε₂ + θ₂μ)/(κ₂ρ₂N₂/ε₂ + θ₁κ₁ρ₁N₁/ε₂ + μ)] R0,2/R0,1.

Adding species 2 dilutes (`R0,C < R0,1`) if and only if `η₂` exceeds it. When habitat 2
lies inside habitat 1 (`θ₂ = 1`) the factor is 1 (eq. 9).
"""
function dilution_threshold(p::EspiraParameters)
    R = species_reproduction_numbers(p)
    a2 = _self_pickup(p, 2)
    b12 = p.θ[1] * p.κ[1] * p.ρ[1] * p.N[1] / p.ε[2]
    return (a2 + b12 + p.θ[2] * p.μ) / (a2 + b12 + p.μ) * R[2] / R[1]
end

"""
    dilution_effect(p::EspiraParameters; criterion = :R0) -> Bool

Whether adding species 2 to species 1 lowers the epidemic potential, `R0,C < R0,1`
(Espira et al. 2022). `criterion = :R0` compares
[`community_reproduction_number`](@ref) with `R0,1` directly; `criterion = :eq8` tests
`η₂ > `[`dilution_threshold`](@ref)`(p)`, the paper's eq. 8.
"""
function dilution_effect(p::EspiraParameters; criterion::Symbol = :R0)
    criterion === :R0 &&
        return community_reproduction_number(p) < species_reproduction_numbers(p)[1]
    criterion === :eq8 && return pickup_efficiency(p)[2] > dilution_threshold(p)
    throw(ArgumentError("criterion must be :R0 or :eq8"))
end

"""
    espira_parameters_from_efficiencies(η₁, η₂, r; θ₂ = 1, ε = (1, 1), μ = 1)
        -> EspiraParameters

Parameters on the axes of Espira et al. (2022), Figs 3 and 4: pick-up efficiencies `η₁`,
`η₂` (in `[0, 1)`), ratio of epidemic potentials `r = R0,2/R0,1`, and overlap `θ₂` (the
fraction of species 2's habitat shared with species 1, the rows of Fig. 4), with the habitat
areas `ε` and the decay rate `μ` fixed as in the figures. `θ₁ = θ₂ ε₂/ε₁` follows from the
definition of θ. Fig. 3 (complete overlap) is `θ₂ = 1` with `ε₁ = ε₂`.

The paper prints no values. This constructor fixes `N = ρ = 1`, `γ = α = 1` and `π₁ = 1`,
sets `κ_i ε_i⁻¹ = μ η_i/(1 - η_i)` so that eq. 4 gives `η_i`, and `π₂` so that eq. 5 gives
`r`; `R0,C/R0,1` depends only on these quantities.
"""
function espira_parameters_from_efficiencies(η1::Real, η2::Real, r::Real; θ₂::Real = 1.0,
        ε = (1.0, 1.0), μ::Real = 1.0)
    0 < η1 < 1 && 0 < η2 < 1 || throw(ArgumentError("η₁ and η₂ must lie in (0, 1)"))
    θ₁ = θ₂ * ε[2] / ε[1]
    θ₁ <= 1 || throw(ArgumentError("θ₁ = θ₂ ε₂/ε₁ = $θ₁ exceeds 1"))
    κ = (μ * η1 / (1 - η1) * ε[1], μ * η2 / (1 - η2) * ε[2])
    return EspiraParameters(; N = (1, 1), γ = (1, 1), κ, ρ = (1, 1),
        π = (1, r * η1 / η2), α = (1, 1), ε, θ = (θ₁, θ₂), μ)
end

function espira_rhs!(du, u, p::EspiraParameters, t)
    S1, S2, I1, I2, R1, R2, E1, E2 = u
    E = (E1, E2)
    S = (S1, S2)
    I = (I1, I2)
    for i in 1:2
        j = _other(i)
        inf = p.κ[i] * p.ρ[i] * p.π[i] * S[i] * (E[i] + p.θ[i] * E[j])
        du[i] = -inf
        du[2 + i] = inf - p.γ[i] * I[i]
        du[4 + i] = p.γ[i] * I[i]
        du[6 + i] = p.α[i] * I[i] / p.ε[i] - E[i] * _loss(p, i)
    end
    return nothing
end

"""
    OverlappingHabitatProblem(p::EspiraParameters; I0 = (1, 0), E0 = (0, 0),
                              tspan = (0.0, 365.0))

The epidemic model of Espira et al. (2022), eq. 1, in days, on the state
`[S₁, S₂, I₁, I₂, R₁, R₂, E₁, E₂]`. Initially `S_i = N_i - I0_i`, `R_i = 0` and the
pathogen concentrations are `E0`. The depletion of the environment uses the constant `N_i`,
as eq. 1 does. `solve(prob)` returns the `ODESolution`; [`to_ode_problem`](@ref) gives the
`ODEProblem`.
"""
struct OverlappingHabitatProblem
    p::EspiraParameters
    u0::Vector{Float64}
    tspan::Tuple{Float64,Float64}
end
function OverlappingHabitatProblem(p::EspiraParameters; I0 = (1.0, 0.0), E0 = (0.0, 0.0),
        tspan = (0.0, 365.0))
    all(i -> 0 <= I0[i] <= p.N[i], 1:2) || throw(ArgumentError("need 0 ≤ I0 ≤ N"))
    all(>=(0), E0) || throw(ArgumentError("E0 must be non-negative"))
    u0 = Float64[p.N[1] - I0[1], p.N[2] - I0[2], I0[1], I0[2], 0, 0, E0[1], E0[2]]
    return OverlappingHabitatProblem(p, u0, (Float64(tspan[1]), Float64(tspan[2])))
end

to_ode_problem(prob::OverlappingHabitatProblem) =
    ODEProblem(espira_rhs!, copy(prob.u0), prob.tspan, prob.p)

CommonSolve.solve(prob::OverlappingHabitatProblem, alg = Tsit5(); abstol = 1e-10,
    reltol = 1e-10, kwargs...) = solve(to_ode_problem(prob), alg; abstol, reltol, kwargs...)

# ---- Helpers for the validation targets (not exported) ----

# Random parameters with log-normal rates. With `geometric = true` the overlaps come from one
# overlap area, θ_i = A/ε_i, as the definition of θ requires; otherwise θ₁, θ₂ are independent.
function _espira_random_parameters(rng::AbstractRNG; geometric::Bool = true)
    ln(s = 1.0) = (exp(s * randn(rng)), exp(s * randn(rng)))
    ε = ln()
    A = rand(rng) * min(ε...)
    θ = geometric ? (A / ε[1], A / ε[2]) : (rand(rng), rand(rng))
    return EspiraParameters(; N = ln(2.0), γ = ln(), κ = ln(), ρ = ln(), π = ln(), α = ln(),
        ε, θ, μ = exp(2 * randn(rng)))
end

# Fraction of `n` random parameter sets for which `f(p)` holds.
function _espira_fraction(f; n::Integer = 2000, seed::Integer = 1, geometric::Bool = true)
    rng = Random.Xoshiro(seed)
    return count(_ -> f(_espira_random_parameters(rng; geometric)), 1:n) / n
end

# R0,C/R0,1 on the axes of Figs 3 and 4.
_espira_ratio(η1, η2, r; kw...) = (q = espira_parameters_from_efficiencies(η1, η2, r; kw...);
    community_reproduction_number(q) / species_reproduction_numbers(q)[1])

# R0,C/R0,1 for the configurations of Fig. 4 (θ₁ = θ₂ = θ, equal habitats), Fig. S1
# (θ₁ = 1, θ₂ = θ: habitat 1 inside habitat 2) and Fig. S2 (θ₂ = 1, θ₁ = θ: habitat 2 inside
# habitat 1), with ε₂ = 1 fixed as in the supplement.
function _espira_panel_ratio(fig::Symbol, η1, η2, r, θ)
    θ1, θ2 = fig === :fig4 ? (θ, θ) : fig === :S1 ? (1.0, θ) : fig === :S2 ? (θ, 1.0) :
             throw(ArgumentError("fig must be :fig4, :S1 or :S2"))
    return _espira_ratio(η1, η2, r; θ₂ = θ2, ε = (θ2 / θ1, 1.0))
end

# Largest |R0,C/R0,1(fig) - R0,C/R0,1(Fig. 4)| over the nine panels (η₁ ∈ {0.1, 0.5, 0.9},
# θ ∈ {0.1, 0.5, 0.9}) and a grid of the axes (η₂, R0,2/R0,1).
function _espira_figure_difference(fig::Symbol)
    d = 0.0
    for θ in (0.1, 0.5, 0.9), η1 in (0.1, 0.5, 0.9), η2 in 0.01:0.02:0.99, r in 0.0:0.02:1.0
        d = max(d, abs(_espira_panel_ratio(fig, η1, η2, r, θ) -
                       _espira_panel_ratio(:fig4, η1, η2, r, θ)))
    end
    return d
end

# The worked example of the supplement (Table S1, eqs S1-S3), from the compendium: the force
# of infection on species 1, the overlap area seen from each habitat, and the three loss rates
# of the pathogens of habitat 1 (PU day⁻¹).
function _espira_worked_example(key::Symbol = :espira2022)
    v = parameter_set(key)
    return (; λ1 = v.π1 * v.κρ1 * ((1 - v.θ1) * v.E1 + v.θ1 * (v.E1 + v.E2)),
        overlap = (v.θ1 * v.ε1, v.θ2 * v.ε2),
        removal_species1 = v.κρ1 * v.N1 * v.E1,
        removal_species2 = v.κρ2 * v.N2 * v.E1,
        decay = v.ε1 * v.μ)
end

# The host block of F V⁻¹ exactly as printed in supplement eq. S12.
function _espira_S12(p::EspiraParameters)
    a(i) = p.κ[i] * p.ρ[i] * p.N[i]
    d1 = p.γ[1] * (a(1) + p.θ[2] * a(2) + p.μ * p.ε[1])
    d2 = p.γ[2] * (a(2) + p.θ[1] * a(1) + p.μ * p.ε[2])
    return [p.α[1] * p.π[1] * a(1)/d1          p.α[2] * p.θ[1] * p.π[1] * a(1)/d2
            p.α[1] * p.θ[2] * p.π[2] * a(2)/d1 p.α[2] * p.π[2] * a(2)/d2]
end

# The right-hand side of supplement eq. S17: dilution iff R0,2/R0,1 is below it.
function _espira_S17(p::EspiraParameters)
    a1, a2 = p.κ[1] * p.ρ[1] * p.N[1], p.κ[2] * p.ρ[2] * p.N[2]
    m = p.μ * p.ε[2]
    return 1 - (p.θ[2] * a2 + p.θ[1] * a1 + p.θ[2] * m) / (a2 + p.θ[1] * a1 + p.θ[2] * m) *
               (m / (a2 + m))
end
