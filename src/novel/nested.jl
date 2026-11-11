# NEW MODEL (not from a published paper): one nested host–pathogen–microbiome model in which
# each piece of biology is a switch, so that models of different complexity can be compared
# like for like. Wiki: lab/wiki/models/new-nested-model.md
#
# It is anchored on Schmidt et al. (2017): with only `latent` on it is their S–L–I model, and
# with everything off it is an S–I model with the same rates. The switches:
#   latent            infected hosts pass through a latent class (rate e) before shedding
#   environment       a share of transmission goes through a pool of zoospores that outlasts
#                     the hosts that shed them (total transmission is unchanged)
#   load              infections start light (shed little, no disease death) and progress to
#                     heavy (shed fully, die of disease)
#   recovery          infected hosts can clear the infection (from the light class if `load`)
#   larvae            an aquatic larval stage between birth and adulthood
#   larval_infection  larvae are infected by adults, without disease, and carry the infection
#                     through metamorphosis with probability `carry_over`
#   fadeout_density   below this density (per ha) of infected animals the pathogen is removed,
#                     and below it of animals in total the host population is removed; the
#                     deterministic equations otherwise let both fall to densities that mean
#                     far less than one animal and then grow back
# and a treated class as in microbiome_treatment.jl (efficacy, tolerance, shedding reduction,
# loss, spread within a stage, waning with bacterial abundance, retention at metamorphosis).
#
# State: for stage s (1 larva, 2 adult), treatment τ (1 untreated, 2 treated) and class k
# (1 susceptible, 2 latent, 3 light, 4 heavy), x[s, τ, k]; then the pool Z. Larvae use classes
# 1 and 4 only (4 = infected, no disease). Time in years, densities per ha.
# Apart from Schmidt's parameters and the carry-over, every value is ILLUSTRATIVE.

"""
    ChytridBiology(; latent = false, environment = false, load = false, recovery = false,
        larvae = false, larval_infection = false, carry_over = 0, ...)
    ChytridBiology(:wang2021; kwargs...)

Which pieces of biology are switched on in [`NestedChytridProblem`](@ref), and the values they
take when on: `environmental_fraction` (share of transmission through the pool) and
`pool_half_life` (yr); `light_shedding` (shedding of a light infection relative to a heavy
one) and `progression` (light to heavy, yr⁻¹); `recovery_rate` (yr⁻¹); `larval_period` (yr),
`larval_death` (yr⁻¹), `adult_to_larva` (ha yr⁻¹) and `carry_over`. The keyed constructor takes
`carry_over` from the compendium.

`fadeout_density` (per ha; 0 switches it off) is a fade-out rule. Once the pathogen has been
introduced, if infected animals plus the pool's equivalent in shedding animals fall below it,
they are set to zero; if all animals (adults and larvae) fall below it, the host population is
set to zero.
"""
Base.@kwdef struct ChytridBiology
    latent::Bool = false
    environment::Bool = false
    load::Bool = false
    recovery::Bool = false
    larvae::Bool = false
    larval_infection::Bool = false
    environmental_fraction::Float64 = 0.5
    pool_half_life::Float64 = 30 / 365
    light_shedding::Float64 = 0.1
    progression::Float64 = 26.0
    recovery_rate::Float64 = 26.0
    larval_period::Float64 = 0.33
    larval_death::Float64 = 1.0
    adult_to_larva::Float64 = 5.2
    carry_over::Float64 = 0.0
    fadeout_density::Float64 = 0.0
    function ChytridBiology(latent, environment, load, recovery, larvae, larval_infection,
            environmental_fraction, pool_half_life, light_shedding, progression, recovery_rate,
            larval_period, larval_death, adult_to_larva, carry_over, fadeout_density)
        (0 <= environmental_fraction <= 1 && pool_half_life > 0 && 0 <= light_shedding <= 1 &&
         progression > 0 && recovery_rate >= 0 && larval_period > 0 && larval_death >= 0 &&
         adult_to_larva >= 0 && 0 <= carry_over <= 1 && fadeout_density >= 0) ||
            throw(ArgumentError("ChytridBiology: a value is out of range"))
        larval_infection && !larvae &&
            throw(ArgumentError("larval_infection needs larvae = true"))
        return new(latent, environment, load, recovery, larvae, larval_infection,
            environmental_fraction, pool_half_life, light_shedding, progression, recovery_rate,
            larval_period, larval_death, adult_to_larva, carry_over, fadeout_density)
    end
end
function ChytridBiology(key::Symbol; kwargs...)
    v = parameter_set(key)
    return ChytridBiology(; carry_over = v.carried_over / v.n_larvae, kwargs...)
end

nested_index(s, τ, k) = (s - 1) * 8 + (τ - 1) * 4 + k
const NESTED_Z = 17

# Every transition of the model, listed once. `emit(from, to, flow)` is called for each, with
# `flow` in animals per ha per year; `from = 0` is a birth and `to = 0` a death. The ODE sums
# the flows (nested_rhs!) and the stochastic model draws events from them
# (StochasticNestedProblem), so the two cannot disagree. Returns dZ/dt for the pool.
function nested_flows(emit::F, u, q, t) where {F}
    p, bio, tr, retention, clock = q
    ix = nested_index
    x(s, τ, k) = max(u[ix(s, τ, k)], 0.0)
    Z = max(u[NESTED_Z], 0.0)
    g = protection(tr, t - clock[])
    eff, tol, sred = g * tr.efficacy, g * tr.tolerance, g * tr.shedding_reduction
    ω = reversion_rate(tr)
    light = bio.load ? bio.light_shedding : 0.0
    shed(s, τ) = (light * x(s, τ, 3) + x(s, τ, 4)) * (τ == 2 ? 1 - sred : 1.0)
    HA = shed(2, 1) + shed(2, 2)
    φ = bio.environment ? bio.environmental_fraction : 0.0
    μZ = log(2) / bio.pool_half_life
    source = (1 - φ) * HA + φ * μZ * Z                 # infectious pressure from adults
    λA = p.β * source
    λW = bio.larval_infection ? bio.adult_to_larva * source : 0.0
    NA = sum(x(2, τ, k) for τ in 1:2, k in 1:4)
    dd = p.d + density_dependence(p) * NA
    m = 1 / bio.larval_period
    σ = bio.recovery ? bio.recovery_rate : 0.0
    entry = bio.latent ? 2 : (bio.load ? 3 : 4)        # class a new adult infection enters

    # adults
    for τ in 1:2
        emit(ix(2, τ, 1), ix(2, τ, entry), (τ == 2 ? 1 - eff : 1.0) * λA * x(2, τ, 1))
        bio.latent && emit(ix(2, τ, 2), ix(2, τ, bio.load ? 3 : 4), p.e * x(2, τ, 2))
        if bio.load
            emit(ix(2, τ, 3), ix(2, τ, 4), bio.progression * x(2, τ, 3))
            emit(ix(2, τ, 3), ix(2, τ, 1), σ * x(2, τ, 3))
        else
            emit(ix(2, τ, 4), ix(2, τ, 1), σ * x(2, τ, 4))
        end
        emit(ix(2, τ, 4), 0, p.d_I * (τ == 2 ? 1 - tol : 1.0) * x(2, τ, 4))
        for k in 1:4
            emit(ix(2, τ, k), 0, dd * x(2, τ, k))
        end
    end
    # births: susceptible adults, light infections, and treated infected adults in proportion
    # to their tolerance (as in ProbioticSLIProblem)
    scale = bio.larvae ? p.b * (m + bio.larval_death) / m : p.b
    births = scale * (x(2, 1, 1) + x(2, 2, 1) + x(2, 1, 3) + x(2, 2, 3) +
                      tol * (x(2, 2, 2) + x(2, 2, 4)))
    emit(0, ix(bio.larvae ? 1 : 2, 1, 1), births)
    # larvae: infection, death and metamorphosis
    if bio.larvae
        for τ in 1:2
            emit(ix(1, τ, 1), ix(1, τ, 4), (τ == 2 ? 1 - eff : 1.0) * λW * x(1, τ, 1))
            keep = τ == 2 ? retention : 1.0            # share that keeps its treatment label
            for (k, carried) in ((1, 0.0), (4, bio.carry_over))
                n = x(1, τ, k)
                emit(ix(1, τ, k), 0, bio.larval_death * n)
                for (τa, share) in ((τ, keep), (1, 1 - keep))
                    share > 0 || continue
                    emit(ix(1, τ, k), ix(2, τa, entry), carried * share * m * n)
                    emit(ix(1, τ, k), ix(2, τa, 1), (1 - carried) * share * m * n)
                end
            end
        end
    end
    # the microbiome: spread within a stage, and loss
    for s in 1:2
        total = sum(x(s, τ, k) for τ in 1:2, k in 1:4)
        frac = total > 0 ? min(sum(x(s, 2, k) for k in 1:4) / total, 1.0) : 0.0
        for k in 1:4
            emit(ix(s, 1, k), ix(s, 2, k), tr.direct * frac * x(s, 1, k))
            emit(ix(s, 2, k), ix(s, 1, k), ω * x(s, 2, k))
        end
    end
    return bio.environment ? HA - μZ * Z : 0.0
end

function nested_rhs!(du, u, q, t)
    fill!(du, 0.0)
    du[NESTED_Z] = nested_flows(u, q, t) do from, to, f
        from > 0 && (du[from] -= f)
        to > 0 && (du[to] += f)
        return nothing
    end
    return nothing
end

function nested_disease_free(p::SchmidtParameters, bio::ChytridBiology)
    u = zeros(NESTED_Z)
    u[nested_index(2, 1, 1)] = p.K
    bio.larvae && (u[nested_index(1, 1, 1)] = p.b * p.K * bio.larval_period)
    return u
end

"""
    invasion_growth_rate(p::SchmidtParameters, biology::ChytridBiology)

Initial growth rate (yr⁻¹) of the pathogen in an untreated host population at its disease-free
state: the largest real part of the eigenvalues of the infected classes' Jacobian, found by
finite differences. Positive means the pathogen invades.
"""
function invasion_growth_rate(p::SchmidtParameters, bio::ChytridBiology)
    u0 = nested_disease_free(p, bio)
    q = (p, bio, MicrobiomeTreatment(), 0.0, Ref(0.0))
    infected = [nested_index(s, 1, k) for s in 1:2 for k in 2:4]
    push!(infected, NESTED_Z)
    f0, f1 = zeros(NESTED_Z), zeros(NESTED_Z)
    nested_rhs!(f0, u0, q, 0.0)
    h = 1e-6
    J = zeros(length(infected), length(infected))
    for (j, c) in enumerate(infected)
        u = copy(u0); u[c] += h
        nested_rhs!(f1, u, q, 0.0)
        J[:, j] = (f1[infected] .- f0[infected]) ./ h
    end
    return maximum(real, eigvals(J))
end

"""
    NestedChytridProblem(p::SchmidtParameters, biology::ChytridBiology,
        treatment = MicrobiomeTreatment(); treat = :adults, coverage = 0, treat_times = [0.0],
        retention = 0, t_intro = 0.0, I0 = 0.01, tspan = (0.0, 10.0))

One host–pathogen–microbiome model with each piece of biology switched on or off by
`biology`. The population starts at its disease-free state; at each of `treat_times` a
fraction `coverage` of the untreated animals in `treat` (`:adults`, `:larvae`, `:both`) is
treated; `I0` infectious adults per ha arrive at `t_intro`. Time in years.

`solve` returns series `t`, `adults`, `infected` (adults in any infected class) and
`min_adults`, `mean_adults`, `final_adults` (fractions of `K`, after `t_intro`),
`final_prevalence` and `final_treated`; and, for the fade-out rule, `pathogen_fadeout` and
`host_extinction` (the times they happened, or `NaN`).
"""
struct NestedChytridProblem{P<:SchmidtParameters}
    p::P
    biology::ChytridBiology
    treatment::MicrobiomeTreatment
    treat::Symbol
    coverage::Float64
    treat_times::Vector{Float64}
    retention::Float64
    t_intro::Float64
    I0::Float64
    tspan::Tuple{Float64,Float64}
end
function NestedChytridProblem(p::SchmidtParameters, biology::ChytridBiology,
        treatment::MicrobiomeTreatment = MicrobiomeTreatment(); treat::Symbol = :adults,
        coverage = 0.0, treat_times = [0.0], retention = 0.0, t_intro = 0.0, I0 = 0.01,
        tspan = (0.0, 10.0))
    treat in (:adults, :larvae, :both) ||
        throw(ArgumentError("treat must be :adults, :larvae or :both"))
    (treat === :adults || biology.larvae) ||
        throw(ArgumentError("treating larvae needs biology.larvae = true"))
    (0 <= coverage <= 1 && 0 <= retention <= 1) ||
        throw(ArgumentError("coverage and retention must be in [0, 1]"))
    return NestedChytridProblem(p, biology, treatment, treat, float(coverage),
        collect(float.(treat_times)), float(retention), float(t_intro), float(I0),
        (float(tspan[1]), float(tspan[2])))
end

function CommonSolve.solve(prob::NestedChytridProblem, alg = Tsit5(); dt = 0.01,
        abstol = 1e-10, reltol = 1e-8, kwargs...)
    p, (t0, t1) = prob.p, prob.tspan
    clock = Ref(t0)
    q = (p, prob.biology, prob.treatment, prob.retention, clock)
    events = Tuple{Float64,Symbol}[(t, :treat) for t in prob.treat_times if t0 <= t < t1]
    t0 <= prob.t_intro < t1 && push!(events, (prob.t_intro, :introduce))
    push!(events, (t1, :stop))
    sort!(events; by = e -> (e[1], e[2] === :introduce))
    u = nested_disease_free(p, prob.biology)
    stages = prob.treat === :adults ? (2,) : prob.treat === :larvae ? (1,) : (1, 2)
    ts, us = Float64[], Vector{Float64}[]
    t = t0
    # fade-out rule
    limit = prob.biology.fadeout_density
    infected_idx = [nested_index(s, τ, k) for s in 1:2 for τ in 1:2 for k in 2:4]
    host_idx = 1:16
    μZ = log(2) / prob.biology.pool_half_life
    pathogen_gone, host_gone = Ref(NaN), Ref(NaN)
    pathogen_left(v) = sum(max(v[i], 0.0) for i in infected_idx) + μZ * max(v[NESTED_Z], 0.0)
    hosts_left(v) = sum(max(v[i], 0.0) for i in host_idx)
    low_pathogen(v, t) = t > prob.t_intro && isnan(pathogen_gone[]) && 0 < pathogen_left(v) < limit
    low_hosts(v) = isnan(host_gone[]) && hosts_left(v) < limit
    function fade!(integ)
        if low_hosts(integ.u)
            integ.u .= 0.0
            host_gone[] = integ.t
            isnan(pathogen_gone[]) && (pathogen_gone[] = integ.t)
        elseif low_pathogen(integ.u, integ.t)
            integ.u[infected_idx] .= 0.0
            integ.u[NESTED_Z] = 0.0
            pathogen_gone[] = integ.t
        end
        return nothing
    end
    extra = limit > 0 ? (callback = DiscreteCallback((v, t, integ) -> low_hosts(v) || low_pathogen(v, t),
        fade!; save_positions = (false, false)),) : (;)
    for (te, kind) in events
        t = probiotic_advance!(u, t, te, nested_rhs!, q, alg, ts, us; dt, abstol, reltol, extra...,
            kwargs...)
        if kind === :treat
            clock[] = te
            for s in stages, k in 1:4
                i, j = nested_index(s, 1, k), nested_index(s, 2, k)
                moved = prob.coverage * u[i]
                u[i] -= moved
                u[j] += moved
            end
        elseif kind === :introduce
            isnan(host_gone[]) && (u[nested_index(2, 1, 4)] += prob.I0)
        end
    end
    A(v, τs, ks) = sum(v[nested_index(2, τ, k)] for τ in τs, k in ks)
    adults = [A(v, 1:2, 1:4) for v in us]
    infected = [A(v, 1:2, 2:4) for v in us]
    after = findall(>=(prob.t_intro), ts)
    rel = adults[after] ./ p.K
    return (t = ts, adults = adults, infected = infected, min_adults = minimum(rel),
        mean_adults = sum(rel) / length(rel), final_adults = rel[end],
        final_prevalence = adults[end] > 0 ? max(infected[end], 0.0) / adults[end] : 0.0,
        final_treated = adults[end] > 0 ? A(us[end], 2:2, 1:4) / adults[end] : 0.0,
        pathogen_fadeout = pathogen_gone[], host_extinction = host_gone[])
end
