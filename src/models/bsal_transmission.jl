# Bsal transmission functions of two eastern newt (Notophthalmus viridescens) studies:
#
#   Malagon DA, Melara LA, Prosper OF, Lenhart S, Carter ED, Fordyce JA, Peterson AC,
#   Miller DL, Gray MJ (2020) Host density and habitat structure influence host contact rates
#   and Batrachochytrium salamandrivorans transmission. Sci Rep 10:5584.
#   doi:10.1038/s41598-020-62351-x. Supplementary File 1 (model and Table S1) and File 2
#   (data, R and MATLAB code).
#
#   Tompros A, Dean AD, Fenton A, Wilber MQ, Carter ED, Gray MJ (2022) Frequency-dependent
#   transmission of Batrachochytrium salamandrivorans in eastern newts. Transbound Emerg Dis
#   69:731-741. doi:10.1111/tbed.14043. Appendix (eqs 1-10) and the TRACE deposit
#   doi:10.7290/xw8yy2RBeu (data and R code).
#
# Equations were read from the papers, the supplements and the authors' code. Discrepancies
# between them are listed in the docstrings and in
# lab/wiki/models/bsal-transmission-functions.md.
#
# A transmission function gives φ(S, I, N), the total rate of new infections (time⁻¹) with S
# susceptible and I infectious hosts in a population of N. Tompros et al. write
# φ = α_{S,N} S I / N, with α the per capita contact rate times the probability that a
# contact infects. S, I and N are numbers of hosts (Tompros: per 1-m² tank) or densities
# (Schmidt: per ha); the parameters carry the matching units.

"""
    TransmissionFunction

Abstract type of a transmission function. A subtype gives the total rate of new infections
`φ(S, I, N)` (per unit time) through [`infection_rate`](@ref). `S` and `I` are the numbers
(or densities) of susceptible and infectious hosts, and `N` is the population size that
enters the contact rate, which is `S + I` unless the model has other classes.

The subtypes are the eight candidates of Tompros et al. (2022, Table 2) and the Holling
type II contact function of Malagon et al. (2020):

| type | `φ(S, I, N)` |
|---|---|
| [`DensityDependent`](@ref) | `β S I` |
| [`FrequencyDependent`](@ref) | `β S I / N` |
| [`PowerLaw`](@ref) | `β Sᵃ Iᵇ` |
| [`QuadraticSaturation`](@ref) | `S I / (κ₀ + κ₁ N + κ₂ N²)` |
| [`StructuralSaturation`](@ref) | `S I / (κ₀ + κ_S S + κ_I I + κ₂ S I)` |
| [`PrevalenceDependent`](@ref) | `β x y` |
| [`PrevalencePowerLaw`](@ref) | `β xᵃ yᵇ` |
| [`PrevalenceSaturation`](@ref) | `x y / (κ₀ + κ₁ x + κ₂ x²)` |
| [`HollingContact`](@ref) | `β c S I / (K + N)` |

with `x = S/N` and `y = I/N`. Use one as the transmission term of an SLI model with
[`TransmissionSLIProblem`](@ref), of the mesocosm master equation with
[`SusceptibleChainProblem`](@ref), or of the cohabitation experiment with
[`CohabitationProblem`](@ref).
"""
abstract type TransmissionFunction end

"""
    DensityDependent(β)

Density-dependent (mass-action) transmission, `φ = β S I`. The per capita contact rate
`α = β N` grows with population size. `β` has units time⁻¹ host⁻¹ (or time⁻¹ area host⁻¹
when `S` and `I` are densities). This is the transmission term of
[`SchmidtSLIProblem`](@ref).
"""
struct DensityDependent <: TransmissionFunction
    β::Float64
end

"""
    FrequencyDependent(β)

Frequency-dependent transmission, `φ = β S I / N`. The per capita contact rate `α = β`
does not change with population size. `β` has units time⁻¹. Tompros et al. (2022) select
this form, with `β = 0.230` day⁻¹ (see [`TomprosParameters`](@ref)).
"""
struct FrequencyDependent <: TransmissionFunction
    β::Float64
end

"""
    PowerLaw(β, a, b)

Power-law transmission, `φ = β Sᵃ Iᵇ` (Tompros et al. 2022, appendix eq. 7). It is
phenomenological. With `a = b = 1` it is [`DensityDependent`](@ref). The exponents are
dimensionless, and `β` has units time⁻¹.
"""
struct PowerLaw <: TransmissionFunction
    β::Float64
    a::Float64
    b::Float64
end

"""
    QuadraticSaturation(κ₀, κ₁, κ₂)

"Quadratic saturation" of Tompros et al. (2022, appendix eq. 4):
`α = N / (κ₀ + κ₁ N + κ₂ N²)`, so `φ = S I / (κ₀ + κ₁ N + κ₂ N²)`. With `κ₁ = κ₂ = 0` it
is density-dependent (`β = 1/κ₀`); with `κ₀ = κ₂ = 0` it is frequency-dependent
(`β = 1/κ₁`); with only `κ₂ = 0` it is a Holling type II contact rate. The κ's have units
of time (time × host for `κ₀`, time / host for `κ₂`).
"""
struct QuadraticSaturation <: TransmissionFunction
    κ₀::Float64
    κ₁::Float64
    κ₂::Float64
end

"""
    StructuralSaturation(κ₀, κ_S, κ_I, κ₂)

"Structural saturation" of Tompros et al. (2022, appendix eq. 5): saturation depends on
the numbers of susceptible and infected hosts separately,
`φ = S I / (κ₀ + κ_S S + κ_I I + κ₂ S I)`. It equals
[`QuadraticSaturation`](@ref)`(κ₀, κ_S, 0)` when `κ_S = κ_I` and `κ₂ = 0`.
"""
struct StructuralSaturation <: TransmissionFunction
    κ₀::Float64
    κ_S::Float64
    κ_I::Float64
    κ₂::Float64
end

"""
    PrevalenceDependent(β)

"Simple prevalence" transmission of Tompros et al. (2022), `φ = β x y` with `x = S/N` and
`y = I/N`. The total infection rate does not grow with population size. `β` has units
time⁻¹.
"""
struct PrevalenceDependent <: TransmissionFunction
    β::Float64
end

"""
    PrevalencePowerLaw(β, a, b)

Power law in the prevalences, `φ = β xᵃ yᵇ` (Tompros et al. 2022, appendix eq. 8).
"""
struct PrevalencePowerLaw <: TransmissionFunction
    β::Float64
    a::Float64
    b::Float64
end

"""
    PrevalenceSaturation(κ₀, κ₁, κ₂)

Quadratic saturation in the prevalence of susceptibles (Tompros et al. 2022, appendix
eq. 6): `α = (1/N) / (κ₀ + κ₁ x + κ₂ x²)`, so `φ = x y / (κ₀ + κ₁ x + κ₂ x²)`. With
`κ₁ = κ₂ = 0` it is [`PrevalenceDependent`](@ref) with `β = 1/κ₀`.
"""
struct PrevalenceSaturation <: TransmissionFunction
    κ₀::Float64
    κ₁::Float64
    κ₂::Float64
end

"""
    HollingContact(β, c, K)
    HollingContact(p::MalagonParameters)

The transmission term of Malagon et al. (2020): a Holling type II contact rate per host,
`f(N) = c N / (K + N)`, times the probability `β` that a contact with an infectious host
transmits, times the frequency term `S I / N`:

```
φ = β f(N) S I / N = β c S I / (K + N)
```

`c` is the asymptotic number of contacts per host per unit time, `K` the half-saturation
population (hosts, or hosts per m² in the paper) and `β` a probability. In the paper `c` is
in contacts newt⁻¹ day⁻¹ (Supplementary Table S1: 2153.29 without plants, 587.31 with
plants) although the text calls `f` a rate per hour; see [`MalagonParameters`](@ref).
"""
struct HollingContact <: TransmissionFunction
    β::Float64
    c::Float64
    K::Float64
end

"""
    infection_rate(tf::TransmissionFunction, S, I, N = S + I)

Total rate of new infections `φ(S, I, N)` (per unit time) under the transmission function
`tf`. The rate is zero when `S` or `I` is zero. For the prevalence forms `N` must be
positive.

```julia
infection_rate(FrequencyDependent(0.23), 6, 2)   # 0.23 * 6 * 2 / 8 = 0.345 day⁻¹
```
"""
function infection_rate end

infection_rate(tf::TransmissionFunction, S, I) = infection_rate(tf, S, I, S + I)
infection_rate(tf::DensityDependent, S, I, N) = tf.β * S * I
infection_rate(tf::FrequencyDependent, S, I, N) = N > 0 ? tf.β * S * I / N : zero(float(N))
infection_rate(tf::PowerLaw, S, I, N) =
    S > 0 && I > 0 ? tf.β * float(S)^tf.a * float(I)^tf.b : zero(float(N))
infection_rate(tf::QuadraticSaturation, S, I, N) =
    S * I / (tf.κ₀ + tf.κ₁ * N + tf.κ₂ * N^2)
infection_rate(tf::StructuralSaturation, S, I, N) =
    S * I / (tf.κ₀ + tf.κ_S * S + tf.κ_I * I + tf.κ₂ * S * I)
infection_rate(tf::PrevalenceDependent, S, I, N) = tf.β * (S / N) * (I / N)
infection_rate(tf::PrevalencePowerLaw, S, I, N) =
    S > 0 && I > 0 ? tf.β * (S / N)^tf.a * (I / N)^tf.b : zero(float(N))
infection_rate(tf::PrevalenceSaturation, S, I, N) =
    (S / N) * (I / N) / (tf.κ₀ + tf.κ₁ * (S / N) + tf.κ₂ * (S / N)^2)
infection_rate(tf::HollingContact, S, I, N) = tf.β * tf.c * S * I / (tf.K + N)

"""
    force_of_infection(tf::TransmissionFunction, S, I, N = S + I)

Per capita rate at which a susceptible host is infected, `φ(S, I, N) / S`.
"""
force_of_infection(tf::TransmissionFunction, S, I, N = S + I) =
    S > 0 ? infection_rate(tf, S, I, N) / S : zero(float(N))

"""
    transmission_coefficient(tf::TransmissionFunction, N) -> Float64

Rate at which one infectious host infects others in an otherwise susceptible population of
size `N`: `∂φ/∂I` at `S = N`, `I = 0`. Divided by the rate at which infectious hosts are
removed it gives R₀ ([`basic_reproduction_number`](@ref)). It is `β` for
[`FrequencyDependent`](@ref), `β N` for [`DensityDependent`](@ref), `β c N/(K + N)` for
[`HollingContact`](@ref) and `N/(κ₀ + κ₁ N + κ₂ N²)` for [`QuadraticSaturation`](@ref). It
is infinite for the power laws with `b < 1`.
"""
function transmission_coefficient end
transmission_coefficient(tf::DensityDependent, N) = tf.β * N
transmission_coefficient(tf::FrequencyDependent, N) = tf.β
transmission_coefficient(tf::QuadraticSaturation, N) = N / (tf.κ₀ + tf.κ₁ * N + tf.κ₂ * N^2)
transmission_coefficient(tf::StructuralSaturation, N) = N / (tf.κ₀ + tf.κ_S * N)
transmission_coefficient(tf::PrevalenceDependent, N) = tf.β / N
transmission_coefficient(tf::PrevalenceSaturation, N) = (1 / N) / (tf.κ₀ + tf.κ₁ + tf.κ₂)
transmission_coefficient(tf::HollingContact, N) = tf.β * tf.c * N / (tf.K + N)
transmission_coefficient(tf::PowerLaw, N) =
    tf.b < 1 ? Inf : tf.b == 1 ? tf.β * float(N)^tf.a : 0.0
transmission_coefficient(tf::PrevalencePowerLaw, N) =
    tf.b < 1 ? Inf : tf.b == 1 ? tf.β / N : 0.0

"""
    basic_reproduction_number(tf::TransmissionFunction, N; removal)

R₀ of an SIR-type model with transmission `tf` in a susceptible population of size `N`,
when infectious hosts leave at total rate `removal` (death plus recovery, same time unit as
`tf`): `transmission_coefficient(tf, N) / removal`. This is the formula of Tompros et al.
(2022, section 2.7), `R₀ = α_{S,N} / (μ + σ)` with `S = N`.
"""
basic_reproduction_number(tf::TransmissionFunction, N; removal::Real) =
    transmission_coefficient(tf, N) / removal

# ---------------------------------------------------------------------------------------
# A Schmidt-style S–L–I model with a pluggable transmission term

"""
    TransmissionSLIProblem(p::SchmidtParameters, tf::TransmissionFunction;
                           u0 = [p.K, 0, 0.01], tspan = (0.0, 1.0))

The S–L–I model of Schmidt et al. (2017, eq. 2.1) with its mass-action term `β S I`
replaced by `φ(S, I, N)` of the transmission function `tf`, where `N = S + L + I`:

```
dS/dt = (b − d) S − γ S N − φ(S, I, N)
dL/dt = φ(S, I, N) − (d + e) L − γ L N
dI/dt = e L − (d + d_I) I − γ I N
```

`p.β` is not used. Time is in years and densities in adult females per hectare, so the
parameters of `tf` must use the same units (for example `FrequencyDependent(β)` with `β`
in yr⁻¹). With `tf = DensityDependent(p.β)` the problem is the same as
[`SchmidtSLIProblem`](@ref). The transmission functions of the newt studies were fitted to
newts in 1-m² tanks over days; their values are not on the scale of this model and are not
transferred to it here.
"""
struct TransmissionSLIProblem{P<:SchmidtParameters,F<:TransmissionFunction}
    p::P
    tf::F
    u0::Vector{Float64}
    tspan::Tuple{Float64,Float64}
end
function TransmissionSLIProblem(p::SchmidtParameters, tf::TransmissionFunction;
        u0 = [p.K, zero(p.K), 0.01], tspan = (0.0, 1.0))
    length(u0) == 3 || throw(ArgumentError("u0 must be [S, L, I]"))
    return TransmissionSLIProblem(p, tf, collect(float.(u0)),
        (Float64(tspan[1]), Float64(tspan[2])))
end

function transmission_sli_rhs!(du, u, q, t)
    p, tf = q
    S, L, I = u
    N = S + L + I
    γ = density_dependence(p)
    λ = infection_rate(tf, max(S, 0.0), max(I, 0.0), N)
    du[1] = (p.b - p.d) * S - γ * S * N - λ
    du[2] = λ - (p.d + p.e) * L - γ * L * N
    du[3] = p.e * L - (p.d + p.d_I) * I - γ * I * N
    return nothing
end

to_ode_problem(prob::TransmissionSLIProblem) =
    ODEProblem(transmission_sli_rhs!, prob.u0, prob.tspan, (prob.p, prob.tf))

CommonSolve.solve(prob::TransmissionSLIProblem, alg = Tsit5(); kwargs...) =
    solve(to_ode_problem(prob), alg; kwargs...)

# ---------------------------------------------------------------------------------------
# Tompros et al. 2022: master equation for the number of susceptibles in a closed tank

"""
    SusceptibleChainProblem(tf::TransmissionFunction, S0, I0; tspan = (0.0, 5.0))

The probabilistic transmission model of Tompros et al. (2022, appendix eq. 2). A closed
population of `N = S0 + I0` hosts has no deaths or recoveries. The number of susceptibles
`S` falls by one at rate `φ(S, N − S, N)`. The probabilities `p_S(t)`, `S = 0, …, S0`,
follow

```
dp_S/dt = φ(S+1, N−S−1, N) p_{S+1} − φ(S, N−S, N) p_S,     p_{S0}(0) = 1.
```

`solve(prob)` integrates this linear ODE and returns the `ODESolution` of the vector
`[p_0, p_1, …, p_{S0}]`. [`susceptible_distribution`](@ref) gives the same probabilities
at one time from the matrix exponential, which is what the likelihood uses. The time unit
is that of `tf` (days for Tompros et al.).
"""
struct SusceptibleChainProblem{F<:TransmissionFunction}
    tf::F
    S0::Int
    I0::Int
    tspan::Tuple{Float64,Float64}
end
function SusceptibleChainProblem(tf::TransmissionFunction, S0::Integer, I0::Integer;
        tspan = (0.0, 5.0))
    S0 >= 0 && I0 >= 0 || throw(ArgumentError("S0 and I0 must be non-negative"))
    return SusceptibleChainProblem(tf, Int(S0), Int(I0), (Float64(tspan[1]), Float64(tspan[2])))
end

# rates φ_S for S = 0..S0 (index S+1)
_chain_rates(tf, S0, N) = [float(infection_rate(tf, S, N - S, N)) for S in 0:S0]

function _chain_generator(tf, S0, N)
    r = _chain_rates(tf, S0, N)
    n = S0 + 1
    A = zeros(n, n)          # dp/dt = A p
    for i in 1:n
        A[i, i] = -r[i]
        i < n && (A[i, i + 1] = r[i + 1])
    end
    return A
end

function chain_rhs!(dp, p, r, t)
    n = length(p)
    for i in 1:n
        dp[i] = -r[i] * p[i] + (i < n ? r[i + 1] * p[i + 1] : 0.0)
    end
    return nothing
end

function to_ode_problem(prob::SusceptibleChainProblem)
    N = prob.S0 + prob.I0
    p0 = zeros(prob.S0 + 1)
    p0[end] = 1.0
    return ODEProblem(chain_rhs!, p0, prob.tspan, _chain_rates(prob.tf, prob.S0, N))
end

CommonSolve.solve(prob::SusceptibleChainProblem, alg = Tsit5(); kwargs...) =
    solve(to_ode_problem(prob), alg; kwargs...)

"""
    susceptible_distribution(tf::TransmissionFunction, S0, I0, t; method = :auto)
        -> Vector{Float64}

Probabilities `[p_0(t), …, p_{S0}(t)]` of the number of susceptibles after time `t` in a
closed population that starts with `S0` susceptible and `I0` infected hosts
([`SusceptibleChainProblem`](@ref)).

- `method = :uniformization` sums the uniformized chain, `Σₖ Poisson(k; q t) Pᵏ e`, with
  `q` the largest rate and `P = I + A/q`. Every term is non-negative, so small
  probabilities keep their relative accuracy. The series is cut where the Poisson tail is
  below `10⁻¹⁶`.
- `method = :exp` applies the matrix exponential `exp(A t)` of the generator `A`.
- `method = :auto` (default) uses uniformization when `q t ≤ 500` and the matrix
  exponential otherwise.

The two agree to about `10⁻¹²`; the authors integrated the master equation with `lsoda`.
"""
function susceptible_distribution(tf::TransmissionFunction, S0::Integer, I0::Integer, t::Real;
        method::Symbol = :auto)
    method in (:auto, :exp, :uniformization) ||
        throw(ArgumentError("method must be :auto, :exp or :uniformization"))
    S0 == 0 && return [1.0]
    r = _chain_rates(tf, S0, S0 + I0)
    q = maximum(r)
    (q * t == 0) && (p = zeros(S0 + 1); p[end] = 1.0; return p)
    if method === :exp || (method === :auto && q * t > 500)
        P = exp(_chain_generator(tf, S0, S0 + I0) * t)
        return max.(P[:, end], 0.0)
    end
    return _uniformized(r, q, float(t))
end

function _uniformized(r::Vector{Float64}, q::Float64, t::Float64)
    n = length(r)
    v = zeros(n); v[end] = 1.0
    w = similar(v)
    qt = q * t
    weight = exp(-qt)            # Poisson(0; qt); qt ≤ 500 keeps it above floatmin
    p = weight .* v
    cum = weight
    k = 0
    kmax = ceil(Int, qt + 12sqrt(qt) + 40)
    while k < kmax && 1 - cum > 1e-16
        # v ← P v with P = I + A/q (A bidiagonal: diagonal −r, superdiagonal r[i+1])
        for i in 1:n
            w[i] = (1 - r[i] / q) * v[i] + (i < n ? r[i + 1] / q * v[i + 1] : 0.0)
        end
        v, w = w, v
        k += 1
        weight *= qt / k
        cum += weight
        @. p += weight * v
    end
    return max.(p, 0.0)
end

"""
    expected_susceptibles(tf::TransmissionFunction, S0, I0, t) -> Float64

Mean number of susceptibles after time `t` under [`susceptible_distribution`](@ref) (the
solid lines of Tompros et al. 2022, Fig. 1).
"""
function expected_susceptibles(tf::TransmissionFunction, S0::Integer, I0::Integer, t::Real)
    p = susceptible_distribution(tf, S0, I0, t)
    return sum((0:S0) .* p)
end

# ---------------------------------------------------------------------------------------
# Malagon et al. 2020: SEI model of the cohabitation experiment

"""
    MalagonParameters(; β, c, K, γ, d)
    MalagonParameters(:malagon2020; habitat = :simple, disease_state = 12, cohabitation = 10,
                      mortality = :published)

Parameters of the SEI model of Malagon et al. (2020), in **days**:

- `β`: probability that a contact with an infected newt transmits Bsal (dimensionless);
- `c`: asymptotic contact rate, contacts newt⁻¹ day⁻¹;
- `K`: half-saturation density, newts m⁻²;
- `γ`: rate from exposed to infectious (qPCR-positive twice), day⁻¹;
- `d`: disease-induced death rate, day⁻¹.

The keyed constructor reads the compendium. `habitat` is `:simple` (no plants) or
`:complex` (nine artificial plants per m²), which sets `c` and `K`; `disease_state` (12, 18
or 24 days since the infected newt was exposed) and `cohabitation` (10 or 30 minutes) pick
the row of Supplementary Table S1 for `β` and `γ`. The paper's Fig. 6 uses day 12 and
10 minutes.

`mortality = :published` uses the death rate printed in the supplement,
`d = −ln(0.8)/10 = 0.0223` day⁻¹. `mortality = :code` uses `d = 1/42.75` day⁻¹, the
value in the authors' MATLAB code (`SEIParameterStructure4.m`).

**Units of `c`.** The paper defines `f(N) = cN/(K + N)` as contacts per newt per *hour*
and constrains `0 ≤ c ≤ 100`. The values in Table S1 (2153.29 and 587.31) are per *day*:
the authors' code fits the hourly rate and multiplies it by 24. The code's bounds are
`c ≤ 1.8 × 10⁵` (hourly) and `K ≤ 100`, the reverse of the bounds printed in the paper.
"""
Base.@kwdef struct MalagonParameters
    β::Float64
    c::Float64
    K::Float64
    γ::Float64
    d::Float64
end

time_unit(::Type{<:MalagonParameters}) = :day

const MALAGON_HABITATS = (simple = "simple", complex = "complex")

function MalagonParameters(key::Symbol; habitat::Symbol = :simple, disease_state::Integer = 12,
        cohabitation::Integer = 10, mortality::Symbol = :published)
    haskey(MALAGON_HABITATS, habitat) ||
        throw(ArgumentError("habitat must be :simple or :complex; got :$habitat"))
    disease_state in (12, 18, 24) || throw(ArgumentError("disease_state must be 12, 18 or 24"))
    cohabitation in (10, 30) || throw(ArgumentError("cohabitation must be 10 or 30 (minutes)"))
    mortality in (:published, :code) ||
        throw(ArgumentError("mortality must be :published or :code"))
    recs = parameter_records(key)
    function rec(sym, set)
        i = findfirst(r -> r.symbol === sym && r.set == set, recs)
        i === nothing && throw(ArgumentError("$key has no $sym in set \"$set\""))
        return recs[i]
    end
    h = MALAGON_HABITATS[habitat]
    row = "$(h)_day$(disease_state)_$(cohabitation)min"
    d = mortality === :published ? rec(:d, "").value :
        1 / in_time_unit(rec(:death_time_code, ""), :day)
    return MalagonParameters(; β = rec(:β, row).value, c = rec(:c, h).value,
        K = rec(:K, h).value, γ = rec(:γ, row).value, d = d)
end

HollingContact(p::MalagonParameters) = HollingContact(p.β, p.c, p.K)

"""
    contact_rate(p, N)

Contacts per newt per unit time at density `N` under the Holling type II function
`f(N) = c N / (K + N)` of Malagon et al. (2020). `p` is a [`MalagonParameters`](@ref)
(contacts newt⁻¹ day⁻¹) or a [`HollingContact`](@ref).
"""
contact_rate(p::Union{MalagonParameters,HollingContact}, N) = p.c * N / (p.K + N)

"""
    CohabitationProblem(p::MalagonParameters; N = 133, introduced = N/2,
                        duration = 10/1440, tspan = (0.0, 100.0),
                        transmission = HollingContact(p))

The SEI model of Malagon et al. (2020, Methods and Supplementary Information) for their
Experiment 2, in days. Infected newts (`introduced`) are housed with susceptible newts for
`duration` days (10 minutes by default) and then removed:

```
dS/dt = −φ(S, I_c(t), N(t))
dE/dt =  φ(S, I_c(t), N(t)) − γ E
dI/dt =  γ E − d I
```

with `I_c(t) = introduced` for `t < duration` and `0` afterwards, and
`N = S + E + I + I_c`. The default transmission is `φ = β c S I_c / (K + N)`
([`HollingContact`](@ref)). The state is `[S, E, I, C, D]`, where `C` (cumulative newts
that became infectious, `∫γE`) and `D` (cumulative deaths, `∫dI`) are what the paper
plots, as `C/S(0)` and `D/S(0)`.

The defaults follow the authors' code (`SolveSEImodel6.m`, `SEImodel4.m`): the experiment
held 2 newts in 0.015 m², which they scaled to `N = 133` newts m⁻² with half of them
infected. The printed equations multiply by the infectious class `I`, but the code uses the
introduced newts `I_c` and adds them to `N`. The two agree only if `I` includes the
introduced newts, because new infections cannot become infectious within 10 minutes.

`solve(prob)` integrates the model with the removal time `duration` as a `tstop` and
returns the `ODESolution`; [`cohabitation_curve`](@ref) gives the same curves in closed form.

**What Fig. 6 shows.** The paper describes Fig. 6 as a simulated *population* after the
introduction of *one* infected newt. The code instead simulates this cohabitation
experiment, with 66.5 infected newts per 66.5 susceptibles for 10 minutes, and no
transmission afterwards.
"""
struct CohabitationProblem{F<:TransmissionFunction}
    p::MalagonParameters
    tf::F
    introduced::Float64
    duration::Float64
    u0::Vector{Float64}
    tspan::Tuple{Float64,Float64}
end
function CohabitationProblem(p::MalagonParameters; N::Real = 133.0, introduced::Real = N / 2,
        duration::Real = 10 / 1440, tspan = (0.0, 100.0),
        transmission::TransmissionFunction = HollingContact(p))
    0 <= introduced <= N || throw(ArgumentError("need 0 ≤ introduced ≤ N"))
    duration >= 0 || throw(ArgumentError("duration must be non-negative"))
    u0 = [N - introduced, 0.0, 0.0, 0.0, 0.0]
    return CohabitationProblem(p, transmission, float(introduced), float(duration), u0,
        (Float64(tspan[1]), Float64(tspan[2])))
end

function cohabitation_rhs!(du, u, prob, t)
    p = prob.p
    S, E, I = u[1], u[2], u[3]
    Ic = t < prob.duration ? prob.introduced : 0.0
    N = S + E + I + Ic
    λ = Ic > 0 ? infection_rate(prob.tf, max(S, 0.0), Ic, N) : 0.0
    du[1] = -λ
    du[2] = λ - p.γ * E
    du[3] = p.γ * E - p.d * I
    du[4] = p.γ * E
    du[5] = p.d * I
    return nothing
end

to_ode_problem(prob::CohabitationProblem) =
    ODEProblem(cohabitation_rhs!, prob.u0, prob.tspan, prob)

# The removal of the infected newts at t = duration is a discontinuity: pass it as a tstop.
CommonSolve.solve(prob::CohabitationProblem, alg = Tsit5(); kwargs...) =
    solve(to_ode_problem(prob), alg; tstops = [prob.duration], kwargs...)

"""
    cohabitation_curve(prob::CohabitationProblem, t) -> (; infectious, dead)

Closed-form cumulative proportions `C(t)/S(0)` (became infectious) and `D(t)/S(0)` (died)
of the cohabitation model, for a time or a vector of times. During the cohabitation the
population `N` is treated as constant; it changes only through deaths, which are of order
`10⁻⁶` of `N` in 10 minutes. The result agrees with `solve(prob)` to about `10⁻⁶`. The grid
search of [`fit_cohabitation`](@ref) uses this form.
"""
function cohabitation_curve(prob::CohabitationProblem, t)
    p = prob.p
    S0, τ, Ic = prob.u0[1], prob.duration, prob.introduced
    N = S0 + Ic
    λ = S0 > 0 && Ic > 0 && τ > 0 ? infection_rate(prob.tf, S0, Ic, N) / S0 : 0.0
    γ, d = p.γ, p.d
    E_τ = _sei_E(λ, γ, τ)                         # exposed fraction at the end of the window
    I_τ = _sei_I(λ, γ, d, τ)                      # infectious fraction at the end of the window
    S_τ = exp(-λ * τ)
    function at(t)
        if t <= τ
            S, E, I = exp(-λ * t), _sei_E(λ, γ, t), _sei_I(λ, γ, d, t)
        else
            s = t - τ
            S = S_τ
            E = E_τ * exp(-γ * s)
            I = I_τ * exp(-d * s) + E_τ * _decay_transfer(γ, d, s)
        end
        C = 1 - S - E
        return (infectious = C, dead = C - I)
    end
    return t isa Real ? at(t) : (r = at.(t); (infectious = getfield.(r, :infectious),
        dead = getfield.(r, :dead)))
end

# E(t) with S(0) = 1, E(0) = 0, dE = λ e^{-λt} − γE
_sei_E(λ, γ, t) = λ == 0 ? 0.0 : (abs(γ - λ) < 1e-12 ? λ * t * exp(-λ * t) :
    λ * (exp(-λ * t) - exp(-γ * t)) / (γ - λ))
# ∫₀ˢ γ e^{-γu} e^{-d(s-u)} du: infectious from one exposed after time s
_decay_transfer(γ, d, s) = abs(γ - d) < 1e-12 ? γ * s * exp(-γ * s) :
    γ * (exp(-d * s) - exp(-γ * s)) / (γ - d)
# I(t) with S(0) = 1: I(t) = λγ Σᵢ e^{−rᵢt} / Πⱼ≠ᵢ (rⱼ − rᵢ), r = (λ, γ, d). Repeated rates
# are separated by a relative 1e-7, which changes I by far less than the grid resolution.
function _sei_I(λ, γ, d, t)
    (λ == 0 || γ == 0) && return 0.0
    r = [λ, γ, d]
    for i in 2:3, j in 1:(i - 1)
        abs(r[i] - r[j]) < 1e-9 * max(r[j], 1e-12) && (r[i] *= 1 + 1e-7 * i)
    end
    s = 0.0
    for i in 1:3
        s += exp(-r[i] * t) / prod(r[j] - r[i] for j in 1:3 if j != i)
    end
    return λ * γ * s
end

"""
    final_size_reproduction_number(z) -> Float64

R₀ from the final proportion infected `z` through the final-size relation of a closed SIR
epidemic, `1 − z = exp(−R₀ z)` (Anderson & May 1991, eq. 6.21): `R₀ = −ln(1 − z)/z`.
Malagon et al. (2020) apply it to the cumulative proportions infected that they read from
their Fig. 6 (0.95 without plants, 0.77 with plants) to get R₀ = 3.2 and 1.9.
"""
function final_size_reproduction_number(z::Real)
    0 < z < 1 || throw(ArgumentError("z must be in (0, 1)"))
    return -log1p(-z) / z
end

# ---------------------------------------------------------------------------------------
# Tompros et al. 2022: the selected model and its R₀

"""
    TomprosParameters(; β, μ, σ)
    TomprosParameters(:tompros2022)

Parameters of the R₀ calculation of Tompros et al. (2022, section 2.7), in **days**:

- `β`: the frequency-dependent transmission rate (day⁻¹), the best model of Table 2;
- `μ`: the disease-induced mortality rate (day⁻¹), taken from Malagon et al. (2020);
- `σ`: the recovery rate (day⁻¹), `1/40` from a personal communication.

`basic_reproduction_number(p)` gives `β / (μ + σ)`; `FrequencyDependent(p)` is the
transmission function.
"""
Base.@kwdef struct TomprosParameters
    β::Float64
    μ::Float64
    σ::Float64
end

time_unit(::Type{<:TomprosParameters}) = :day

function TomprosParameters(key::Symbol)
    recs = parameter_records(key)
    function rec(sym)
        i = findfirst(r -> r.symbol === sym && isempty(r.set), recs)
        i === nothing && throw(ArgumentError("$key has no $sym"))
        return recs[i]
    end
    return TomprosParameters(; β = rec(:β).value, μ = rec(:μ).value,
        σ = 1 / in_time_unit(rec(:recovery_time), :day))
end

FrequencyDependent(p::TomprosParameters) = FrequencyDependent(p.β)

"""
    basic_reproduction_number(p::TomprosParameters)

R₀ of Tompros et al. (2022): `β / (μ + σ)` (frequency-dependent transmission, so it does
not depend on population size). The paper reports 4.9.
"""
basic_reproduction_number(p::TomprosParameters) =
    basic_reproduction_number(FrequencyDependent(p), 1; removal = p.μ + p.σ)
