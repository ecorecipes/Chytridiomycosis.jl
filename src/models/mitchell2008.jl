# Mitchell, Churcher, Garner & Fisher (2008) Proc R Soc B 275:329. Bufo bufo at one pond:
# uninfected and infected tadpoles (P_X, P_Y), uninfected juveniles structured by age (J_X),
# adults (A_X) and free-living zoospores (Z) that may reproduce saprobically.
# Units: years, individuals, zoospores. Wiki: lab/wiki/models/bd-mitchell2008-pond.md
#
#   dP_X/dt = δ(T−t1) ιψ A_X − (μP + νZ) P_X − δ(T−t2) P_X                      (2.1)
#   dP_Y/dt = νZ P_X − (μP + αP) P_Y − δ(T−t2) P_Y                              (2.2)
#   ∂J/∂t + ∂J/∂a = −μJ J (0 < a < ω);  −(μJ + ε) J (a ≥ ω);  J(t,0) = δ(T−t2) P_X   (2.3)
#   dA_X/dt = ε ∫_ω^∞ J da (1 − A_X/K) − μA A_X                                  (2.4)
#   dZ/dt   = ρ P_Y + γ Z exp(−φZ) − μZ Z                                        (2.5)
# T is the time since the start of the year. Infected tadpoles are removed at metamorphosis
# (the paper predicts infected juveniles all die).
#
# How the age structure is solved: juveniles younger than ω only die, at a constant rate, so
# each yearly cohort is tracked in closed form and joins a pool M = ∫_ω^∞ J da when it
# reaches age ω. The ODE state is (P_X, P_Y, M, A_X, Z); births, metamorphosis and cohort
# maturation are applied between integration segments. This is exact for equation (2.3).
#
# Choices the paper leaves open:
#   time of year of the Bd pulse          -> the start of a year (T = 0)
#   when the extinction limits are tested -> at each event (three times a year), and for Bd
#                                            only when no tadpoles are present (before births)
#   what "toad population set to 0" covers -> every toad stage

"""
    MitchellParameters(; t1, t2, ω, ι, ψ, ε, μP, μJ, μA, K, αP, ν, ρ, μZ, γ, φ)
    MitchellParameters(:mitchell2008; set = "")

Parameters of the Mitchell et al. (2008) pond model, in years (`t1`, `t2` as fractions of a
year, `ω` in years). The compendium stores `t1`, `t2` and `ω` in days; the keyed constructor
divides them by 365.
"""
Base.@kwdef struct MitchellParameters
    t1::Float64
    t2::Float64
    ω::Float64
    ι::Float64
    ψ::Float64
    ε::Float64
    μP::Float64
    μJ::Float64
    μA::Float64
    K::Float64
    αP::Float64
    ν::Float64
    ρ::Float64
    μZ::Float64
    γ::Float64
    φ::Float64
end

function MitchellParameters(key::Symbol; set::AbstractString = "")
    p = parameter_set(key; set)
    return MitchellParameters(; t1 = p.t1 / 365, t2 = p.t2 / 365, ω = p.ω / 365, ι = p.ι,
        ψ = p.ψ, ε = p.ε, μP = p.μP, μJ = p.μJ, μA = p.μA, K = p.K, αP = p.αP, ν = p.ν,
        ρ = p.ρ, μZ = p.μZ, γ = p.γ, φ = p.φ)
end

function mitchell_rhs!(du, u, p::MitchellParameters, t)
    PX, PY, M, A, Z = u
    du[1] = -(p.μP + p.ν * Z) * PX
    du[2] = p.ν * Z * PX - (p.μP + p.αP) * PY
    du[3] = -(p.μJ + p.ε) * M
    du[4] = p.ε * M * (1 - A / p.K) - p.μA * A
    du[5] = p.ρ * PY + p.γ * Z * exp(-p.φ * Z) - p.μZ * Z
    return nothing
end

"""
    MitchellPondProblem(p; years = 35, spinup = 100, Z0 = 1e4, A0 = p.K / 2,
                        toad_limit = 50.0, bd_limit = 100.0)

Mitchell et al. (2008) model: run without Bd for `spinup` years from `A0` adults, add `Z0`
zoospores at time 0, and run for `years` more. Toads are set to zero if adults fall below
`toad_limit`, and Bd if zoospores fall below `bd_limit`.

`solve(prob)` returns a NamedTuple of series at every event time (`t`, `PX`, `PY`, `J`, `A`,
`Z`, with values taken just after the event), `toads_persist`, `bd_persist` and
`A_disease_free` (adults just before Bd is introduced).
"""
struct MitchellPondProblem
    p::MitchellParameters
    years::Int
    spinup::Int
    Z0::Float64
    A0::Float64
    toad_limit::Float64
    bd_limit::Float64
end
MitchellPondProblem(p::MitchellParameters; years::Integer = 35, spinup::Integer = 100,
    Z0 = 1e4, A0 = p.K / 2, toad_limit = 50.0, bd_limit = 100.0) =
    MitchellPondProblem(p, years, spinup, Z0, A0, toad_limit, bd_limit)

function mitchell_advance!(u, t, tnext, p, alg; kwargs...)
    tnext > t || return t
    sol = solve(ODEProblem(mitchell_rhs!, u, (t, tnext), p), alg; save_everystep = false,
        kwargs...)
    SciMLBase.successful_retcode(sol) || error("Mitchell model: solver failed after t = $t")
    u .= max.(sol.u[end], 0.0)
    return tnext
end

function CommonSolve.solve(prob::MitchellPondProblem, alg = Tsit5(); abstol = 1e-6,
        reltol = 1e-8, kwargs...)
    p = prob.p
    0 <= p.t1 < p.t2 < 1 || throw(ArgumentError("need 0 ≤ t1 < t2 < 1 (fractions of a year)"))
    u = [0.0, 0.0, 0.0, prob.A0, 0.0]            # PX, PY, M, A, Z
    cohorts = Tuple{Float64,Float64}[]           # (time the cohort reaches age ω, size then)
    events = Tuple{Float64,Symbol}[]
    for y in -prob.spinup:(prob.years - 1)
        push!(events, (y + p.t1, :birth), (y + p.t2, :metamorphosis))
    end
    push!(events, (0.0, :introduce), (float(prob.years), :stop))
    sort!(events; by = first)
    out = (t = Float64[], PX = Float64[], PY = Float64[], J = Float64[], A = Float64[],
        Z = Float64[])
    t = float(-prob.spinup)
    A_disease_free = NaN
    for (te, kind) in events
        # cohorts that reach age ω before this event join the maturing pool on the way
        while !isempty(cohorts) && first(cohorts)[1] <= te
            tm, size = popfirst!(cohorts)
            t = mitchell_advance!(u, t, tm, p, alg; abstol, reltol, kwargs...)
            u[3] += size
        end
        t = mitchell_advance!(u, t, te, p, alg; abstol, reltol, kwargs...)
        if kind === :birth
            0 < u[5] < prob.bd_limit && (u[5] = 0.0)       # no tadpoles now: Bd is at its low
            u[1] += p.ι * p.ψ * u[4]
        elseif kind === :metamorphosis
            push!(cohorts, (te + p.ω, u[1] * exp(-p.μJ * p.ω)))
            u[1] = 0.0
            u[2] = 0.0
        elseif kind === :introduce
            A_disease_free = u[4]
            u[5] += prob.Z0
        end
        if 0 < u[4] < prob.toad_limit                      # toads extinct: every stage
            u[1] = u[2] = u[3] = u[4] = 0.0
            empty!(cohorts)
        end
        immature = sum((s * exp(p.μJ * (tm - te)) for (tm, s) in cohorts); init = 0.0)
        push!(out.t, te); push!(out.PX, u[1]); push!(out.PY, u[2])
        push!(out.J, u[3] + immature); push!(out.A, u[4]); push!(out.Z, u[5])
    end
    return merge(out, (toads_persist = u[4] > 0, bd_persist = u[5] >= prob.bd_limit,
        A_disease_free = A_disease_free))
end
