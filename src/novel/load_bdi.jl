# NEW MODEL (not from a published chytrid paper): the stochastic counterpart of the Briggs
# et al. (2010) within-host equation is a linear birth–death–immigration (BDI) process, whose
# distribution is negative binomial in closed form (Kendall 1948, Ann Math Stat 19:1).
# Wiki: lab/wiki/models/new-load-bdi.md
#
# For one host, counting sporangia S, with the zoospore pool Z held fixed:
#   immigration  ι = νγZ/V   new sporangia per day from the pool
#   birth        λ = ηνf     per sporangium per day, by self-reinfection
#   death        σ           per sporangium per day
# The mean obeys the Briggs equation dE[S]/dt = ι + (λ − σ) E[S].
# Units: days, sporangia per host.

"""
    LoadBDI(ι, λ, σ)
    LoadBDI(p::BriggsParameters; Z)

Birth–death–immigration process for the sporangial load on one host: immigration `ι`
(sporangia day⁻¹), per-sporangium self-reinfection `λ` and loss `σ` (day⁻¹).
"""
struct LoadBDI{T<:Real}
    ι::T
    λ::T
    σ::T
    function LoadBDI(ι, λ, σ)
        ι, λ, σ = promote(float(ι), float(λ), float(σ))
        (ι >= 0 && λ > 0 && σ > 0) || throw(ArgumentError("need ι ≥ 0, λ > 0, σ > 0"))
        return new{typeof(ι)}(ι, λ, σ)
    end
end
LoadBDI(p::BriggsParameters; Z) = LoadBDI(transmission_rate(p) * Z / p.V, p.η * p.ν * p.f, p.σ)

"""
    aggregation(m::LoadBDI) -> k = ι/λ

Negative-binomial shape of the load distribution: infection from the pool relative to
self-reinfection. It does not change with time, with host density or with swab efficiency.
"""
aggregation(m::LoadBDI) = m.ι / m.λ

# Negative-binomial parameter q at time t for a host with no sporangia at t = 0; t = Inf
# gives the stationary value λ/σ (needs λ < σ).
function bdi_q(m::LoadBDI, t)
    λ, σ = m.λ, m.σ
    if isinf(t)
        λ < σ || throw(ArgumentError("no stationary load distribution when λ ≥ σ"))
        return λ / σ
    end
    λ == σ && return λ * t / (1 + λ * t)
    E = exp((λ - σ) * t)
    return λ * (E - 1) / (λ * E - σ)
end

"""
    load_distribution(m::LoadBDI, t = Inf) -> NegativeBinomial

Distribution of sporangia on a host that was uninfected at time 0, after `t` days at constant
exposure. `t = Inf` is the enzootic stationary distribution (`λ < σ`). For `λ > σ` the mean
grows without bound while the shape `k` stays at `ι/λ`.
"""
function load_distribution(m::LoadBDI, t::Real = Inf)
    m.ι > 0 || throw(ArgumentError("ι = 0: the load stays at zero"))
    return NegativeBinomial(aggregation(m), 1 - bdi_q(m, t))
end

"""Mean sporangia per host at time `t` (the Briggs equation's solution from S(0) = 0)."""
function mean_load(m::LoadBDI, t::Real = Inf)
    q = bdi_q(m, t)
    return aggregation(m) * q / (1 - q)
end

"""Probability that a host carries at least one sporangium at time `t`."""
prevalence(m::LoadBDI, t::Real = Inf) = 1 - (1 - bdi_q(m, t))^aggregation(m)

"""Rate (day⁻¹) at which an uninfected host becomes infected: `ι`."""
infection_gain_rate(m::LoadBDI) = m.ι

"""
    infection_loss_rate(m::LoadBDI)

Stationary rate (day⁻¹) at which an infected host returns to zero load. By detailed balance
it is `ι π₀ / (1 − π₀)`, with `π₀` the uninfected fraction.
"""
function infection_loss_rate(m::LoadBDI)
    π0 = 1 - prevalence(m)
    return m.ι * π0 / (1 - π0)
end

"""
    swab(d::NegativeBinomial, ρ) -> NegativeBinomial

Distribution of the count on a swab that recovers each sporangium independently with
probability `ρ`. The shape `k` is unchanged, so `k` can be estimated without knowing `ρ`.
"""
function swab(d::NegativeBinomial, ρ::Real)
    0 < ρ <= 1 || throw(ArgumentError("ρ must be in (0, 1]"))
    q = 1 - d.p
    return NegativeBinomial(d.r, 1 - ρ * q / (1 - q + ρ * q))
end

"""
    enzootic_pool(p::BriggsParameters; N, ε) -> Z*

Equilibrium zoospore pool when `N` hosts at their mean stationary load share a pool that also
receives `ε` zoospores per day from outside (a reservoir, tadpoles or another host species).
It is finite only below the host-density threshold `host_density_threshold(p)`, and diverges
as `N` approaches it.
"""
function enzootic_pool(p::BriggsParameters; N::Real, ε::Real)
    a = p.σ - p.η * p.ν * p.f
    a > 0 || throw(ArgumentError("self-reinfection alone sustains growth (ηνf ≥ σ)"))
    denom = p.μ + p.γ * N / p.V - p.η * (1 - p.f) * N * transmission_rate(p) / (p.V * a)
    denom > 0 || throw(ArgumentError("N is at or above the host-density threshold"))
    return ε / denom
end

"""
    enzootic_load_distribution(p::BriggsParameters; N, ε) -> NegativeBinomial

Stationary distribution of sporangia per host for `N` hosts below the density threshold, kept
infected by an external zoospore input `ε` (day⁻¹).
"""
enzootic_load_distribution(p::BriggsParameters; N::Real, ε::Real) =
    load_distribution(LoadBDI(p; Z = enzootic_pool(p; N, ε)))
