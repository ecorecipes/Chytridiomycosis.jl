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
    function ChytridBiology(latent, environment, load, recovery, larvae, larval_infection,
            environmental_fraction, pool_half_life, light_shedding, progression, recovery_rate,
            larval_period, larval_death, adult_to_larva, carry_over)
        (0 <= environmental_fraction <= 1 && pool_half_life > 0 && 0 <= light_shedding <= 1 &&
         progression > 0 && recovery_rate >= 0 && larval_period > 0 && larval_death >= 0 &&
         adult_to_larva >= 0 && 0 <= carry_over <= 1) ||
            throw(ArgumentError("ChytridBiology: a value is out of range"))
        larval_infection && !larvae &&
            throw(ArgumentError("larval_infection needs larvae = true"))
        return new(latent, environment, load, recovery, larvae, larval_infection,
            environmental_fraction, pool_half_life, light_shedding, progression, recovery_rate,
            larval_period, larval_death, adult_to_larva, carry_over)
    end
end
function ChytridBiology(key::Symbol; kwargs...)
    v = parameter_set(key)
    return ChytridBiology(; carry_over = v.carried_over / v.n_larvae, kwargs...)
end

nested_index(s, τ, k) = (s - 1) * 8 + (τ - 1) * 4 + k
const NESTED_Z = 17

function nested_rhs!(du, u, q, t)
    p, bio, tr, retention, clock = q
    x(s, τ, k) = max(u[nested_index(s, τ, k)], 0.0)
    Z = max(u[NESTED_Z], 0.0)
    fill!(du, 0.0)
    g = protection(tr, t - clock[])
    eff, tol, sred = g * tr.efficacy, g * tr.tolerance, g * tr.shedding_reduction
    ω = reversion_rate(tr)
    light = bio.load ? bio.light_shedding : 0.0
    shed(s, τ) = (light * x(s, τ, 3) + x(s, τ, 4)) * (τ == 2 ? 1 - sred : 1.0)
    HA, HW = shed(2, 1) + shed(2, 2), shed(1, 1) + shed(1, 2)
    φ = bio.environment ? bio.environmental_fraction : 0.0
    μZ = log(2) / bio.pool_half_life
    source = (1 - φ) * HA + φ * μZ * Z                 # infectious pressure from adults
    λ = (p.β * source, bio.larval_infection ? bio.adult_to_larva * source : 0.0)   # adult, larva
    du[NESTED_Z] = bio.environment ? HA - μZ * Z : 0.0
    NA = sum(x(2, τ, k) for τ in 1:2, k in 1:4)
    dd = p.d + density_dependence(p) * NA
    m = 1 / bio.larval_period
    σ = bio.recovery ? bio.recovery_rate : 0.0
    entry = bio.latent ? 2 : (bio.load ? 3 : 4)        # class a new adult infection enters
    add!(s, τ, k, v) = (du[nested_index(s, τ, k)] += v)

    # adults
    for τ in 1:2
        inf = (τ == 2 ? 1 - eff : 1.0) * λ[1] * x(2, τ, 1)
        add!(2, τ, 1, -inf); add!(2, τ, entry, inf)
        if bio.latent
            out = p.e * x(2, τ, 2)
            add!(2, τ, 2, -out); add!(2, τ, bio.load ? 3 : 4, out)
        end
        α = p.d_I * (τ == 2 ? 1 - tol : 1.0)
        if bio.load
            add!(2, τ, 3, -bio.progression * x(2, τ, 3)); add!(2, τ, 4, bio.progression * x(2, τ, 3))
            add!(2, τ, 3, -σ * x(2, τ, 3)); add!(2, τ, 1, σ * x(2, τ, 3))
        else
            add!(2, τ, 4, -σ * x(2, τ, 4)); add!(2, τ, 1, σ * x(2, τ, 4))
        end
        add!(2, τ, 4, -α * x(2, τ, 4))
        for k in 1:4
            add!(2, τ, k, -dd * x(2, τ, k))
        end
    end
    # births: susceptible adults, light infections, and treated infected adults in proportion
    # to their tolerance (as in ProbioticSLIProblem)
    scale = bio.larvae ? p.b * (m + bio.larval_death) / m : p.b
    births = scale * (x(2, 1, 1) + x(2, 2, 1) + x(2, 1, 3) + x(2, 2, 3) +
                      tol * (x(2, 2, 2) + x(2, 2, 4)))
    add!(bio.larvae ? 1 : 2, 1, 1, births)
    # larvae
    if bio.larvae
        for τ in 1:2
            inf = (τ == 2 ? 1 - eff : 1.0) * λ[2] * x(1, τ, 1)
            add!(1, τ, 1, -inf); add!(1, τ, 4, inf)
            keep = τ == 2 ? retention : 1.0            # share that keeps its treatment label
            for (k, carried) in ((1, 0.0), (4, bio.carry_over))
                n = x(1, τ, k)
                add!(1, τ, k, -(m + bio.larval_death) * n)
                for (τa, share) in ((τ, keep), (1, 1 - keep))
                    share > 0 || continue
                    add!(2, τa, entry, carried * share * m * n)
                    add!(2, τa, 1, (1 - carried) * share * m * n)
                end
            end
        end
    end
    # the microbiome: loss, and spread within a stage
    for s in 1:2
        total = sum(x(s, τ, k) for τ in 1:2, k in 1:4)
        frac = total > 0 ? min(sum(x(s, 2, k) for k in 1:4) / total, 1.0) : 0.0
        for k in 1:4
            move = tr.direct * frac * x(s, 1, k) - ω * x(s, 2, k)
            add!(s, 1, k, -move); add!(s, 2, k, move)
        end
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
`final_prevalence` and `final_treated`.
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
    for (te, kind) in events
        t = probiotic_advance!(u, t, te, nested_rhs!, q, alg, ts, us; dt, abstol, reltol, kwargs...)
        if kind === :treat
            clock[] = te
            for s in stages, k in 1:4
                i, j = nested_index(s, 1, k), nested_index(s, 2, k)
                moved = prob.coverage * u[i]
                u[i] -= moved
                u[j] += moved
            end
        elseif kind === :introduce
            u[nested_index(2, 1, 4)] += prob.I0
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
        final_treated = adults[end] > 0 ? A(us[end], 2:2, 1:4) / adults[end] : 0.0)
end
