# NEW MODEL (not from a published paper): the Mitchell et al. (2008) toad–Bd pond model with
# infection after metamorphosis. Wiki: lab/wiki/models/new-toad-postmetamorphic.md
#
# Why: Mitchell's model has Bd infecting tadpoles only, with every infected tadpole removed at
# metamorphosis and adults playing no part. Experiments show that toadlets exposed after
# metamorphosis are infected, that many die and some clear the infection (Smith et al. 2022:
# 72 of 101 died), and surveys find Bd on breeding adults at low prevalence.
#
# Additions to equations (2.1)–(2.5), state (P_X, P_Y, M, A_X, A_Y, Z), years:
#   at metamorphosis (T = t2), with Z the zoospore number then:
#     q_J = 1 − exp(−ν Z × metamorph_exposure_days/365)   uninfected metamorphs infected as
#                                                         they leave the pond
#     new juvenile cohort = P_X (1 − q_J + q_J s_J) + P_Y s_J
#     s_J = metamorph_survival: probability that an infected metamorph clears Bd and lives.
#     (Background juvenile mortality then applies to everyone, as in the paper.)
#   adults are at the pond from T = t1 for adult_exposure_days; while they are there
#     dA_X/dt = ... − ν_A Z A_X + σ_A A_Y           ν_A = adult_relative_susceptibility × ν
#     dA_Y/dt = ν_A Z A_X − (μA + α_A + σ_A) A_Y
#     dZ/dt   = ... + ρ_A A_Y                        ρ_A = adult_relative_shedding × ρ
#   away from the pond infected adults only clear (σ_A) or die (μA + α_A).
#   Infected adults breed like uninfected ones.
# With every addition at its default this is the Mitchell model. Every added value is
# ILLUSTRATIVE unless taken from the compendium.

"""
    PostMetamorphicInfection(; metamorph_survival = 0, metamorph_exposure_days = 0,
        adult_exposure_days = 0, adult_relative_susceptibility = 1,
        adult_relative_shedding = 0, adult_clearance = 0, adult_disease_death = 0)

Infection after metamorphosis for the Mitchell toad model. `metamorph_survival` is the
probability that an infected metamorph clears Bd and survives (0 in Mitchell et al. 2008).
Uninfected metamorphs are exposed to the pond for `metamorph_exposure_days` as they leave.
Adults are at the pond for `adult_exposure_days` from the start of breeding, are infected at
`adult_relative_susceptibility` times the tadpole rate, shed at `adult_relative_shedding`
times the tadpole rate while there, clear infection at `adult_clearance` and die of it at
`adult_disease_death` (both yr⁻¹).
"""
Base.@kwdef struct PostMetamorphicInfection
    metamorph_survival::Float64 = 0.0
    metamorph_exposure_days::Float64 = 0.0
    adult_exposure_days::Float64 = 0.0
    adult_relative_susceptibility::Float64 = 1.0
    adult_relative_shedding::Float64 = 0.0
    adult_clearance::Float64 = 0.0
    adult_disease_death::Float64 = 0.0
    function PostMetamorphicInfection(ms, me, ae, ars, arsh, ac, ad)
        (0 <= ms <= 1 && me >= 0 && ae >= 0 && ars >= 0 && arsh >= 0 && ac >= 0 && ad >= 0) ||
            throw(ArgumentError("need 0 ≤ metamorph_survival ≤ 1 and all other values ≥ 0"))
        return new(ms, me, ae, ars, arsh, ac, ad)
    end
end

function toad_stage_rhs!(du, u, q, t)
    p, s, present = q
    PX, PY, M, AX, AY, Z = (max(x, 0.0) for x in u)
    on = present[] ? 1.0 : 0.0
    inf = on * s.adult_relative_susceptibility * p.ν * Z * AX
    du[1] = -(p.μP + p.ν * Z) * PX
    du[2] = p.ν * Z * PX - (p.μP + p.αP) * PY
    du[3] = -(p.μJ + p.ε) * M
    du[4] = p.ε * M * (1 - (AX + AY) / p.K) - p.μA * AX + s.adult_clearance * AY - inf
    du[5] = inf - (p.μA + s.adult_disease_death + s.adult_clearance) * AY
    du[6] = p.ρ * PY + on * s.adult_relative_shedding * p.ρ * AY +
            p.γ * Z * exp(-p.φ * Z) - p.μZ * Z
    return nothing
end

"""
    ToadStageProblem(p::MitchellParameters, stages::PostMetamorphicInfection; years = 35,
        spinup = 100, Z0 = 1e4, A0 = p.K / 2, toad_limit = 50.0, bd_limit = 100.0)

The Mitchell et al. (2008) pond model with infection after metamorphosis (`stages`). Set up
and extinction limits are as in [`MitchellPondProblem`](@ref), except that zoospores are not
set to zero at breeding while infected adults that shed are present.

`solve` returns series at each event (`t`, `A`, `A_infected`, `Z`), `toads_persist`,
`bd_persist` (zoospores above the limit, or shedding adults still infected), `A_disease_free`,
`A_final`, `A_mean` (adults over the last 10 years) and `adult_prevalence` (the highest
infected fraction of adults in the last 10 years, taken when adults leave the pond).
"""
struct ToadStageProblem
    p::MitchellParameters
    stages::PostMetamorphicInfection
    years::Int
    spinup::Int
    Z0::Float64
    A0::Float64
    toad_limit::Float64
    bd_limit::Float64
end
ToadStageProblem(p::MitchellParameters, stages::PostMetamorphicInfection;
    years::Integer = 35, spinup::Integer = 100, Z0 = 1e4, A0 = p.K / 2, toad_limit = 50.0,
    bd_limit = 100.0) = ToadStageProblem(p, stages, Int(years), Int(spinup), float(Z0),
    float(A0), float(toad_limit), float(bd_limit))

function CommonSolve.solve(prob::ToadStageProblem, alg = Tsit5(); abstol = 1e-6,
        reltol = 1e-8, kwargs...)
    p, s = prob.p, prob.stages
    0 <= p.t1 < p.t2 < 1 || throw(ArgumentError("need 0 ≤ t1 < t2 < 1 (fractions of a year)"))
    stay = s.adult_exposure_days / 365
    stay < 1 || throw(ArgumentError("adult_exposure_days must be less than a year"))
    present = Ref(false)
    q = (p, s, present)
    u = [0.0, 0.0, 0.0, prob.A0, 0.0, 0.0]        # PX, PY, M, AX, AY, Z
    cohorts = Tuple{Float64,Float64}[]
    events = Tuple{Float64,Int,Symbol}[]          # (time, order at ties, kind)
    for y in -prob.spinup:(prob.years - 1)
        push!(events, (y + p.t1, 2, :birth), (y + p.t2, 3, :metamorphosis))
        stay > 0 && push!(events, (y + p.t1 + stay, 1, :adults_leave))
    end
    push!(events, (0.0, 0, :introduce), (float(prob.years), 4, :stop))
    filter!(e -> e[1] <= prob.years, events)
    sort!(events; by = e -> (e[1], e[2]))
    out = (t = Float64[], A = Float64[], A_infected = Float64[], Z = Float64[])
    function step!(t, tn)
        tn > t || return t
        sol = solve(ODEProblem(toad_stage_rhs!, copy(u), (t, tn), q), alg;
            save_everystep = false, abstol, reltol, kwargs...)
        SciMLBase.successful_retcode(sol) || error("toad stage model: solver failed after t = $t")
        u .= max.(sol.u[end], 0.0)
        return tn
    end
    t = float(-prob.spinup)
    A_disease_free = NaN
    prevalence = 0.0
    for (te, _, kind) in events
        while !isempty(cohorts) && first(cohorts)[1] <= te
            tm, size = popfirst!(cohorts)
            t = step!(t, tm)
            u[3] += size
        end
        t = step!(t, te)
        if kind === :birth
            carriers = stay > 0 && s.adult_relative_shedding > 0 && u[5] >= 1
            0 < u[6] < prob.bd_limit && !carriers && (u[6] = 0.0)
            u[1] += p.ι * p.ψ * (u[4] + u[5])
            present[] = stay > 0
        elseif kind === :adults_leave
            present[] = false
            A = u[4] + u[5]
            te >= prob.years - 10 && A > 0 && (prevalence = max(prevalence, u[5] / A))
        elseif kind === :metamorphosis
            qJ = -expm1(-p.ν * u[6] * s.metamorph_exposure_days / 365)
            sJ = s.metamorph_survival
            push!(cohorts, (te + p.ω, (u[1] * (1 - qJ + qJ * sJ) + u[2] * sJ) * exp(-p.μJ * p.ω)))
            u[1] = u[2] = 0.0
        elseif kind === :introduce
            A_disease_free = u[4] + u[5]
            u[6] += prob.Z0
        end
        if 0 < u[4] + u[5] < prob.toad_limit
            u[1] = u[2] = u[3] = u[4] = u[5] = 0.0
            empty!(cohorts)
        end
        push!(out.t, te); push!(out.A, u[4] + u[5]); push!(out.A_infected, u[5])
        push!(out.Z, u[6])
    end
    late = findall(>=(prob.years - 10), out.t)
    carriers = s.adult_relative_shedding > 0 && stay > 0 && u[5] >= 1
    return merge(out, (toads_persist = u[4] + u[5] > 0,
        bd_persist = u[6] >= prob.bd_limit || carriers, A_disease_free = A_disease_free,
        A_final = u[4] + u[5], A_mean = sum(out.A[late]) / length(late),
        adult_prevalence = prevalence))
end
