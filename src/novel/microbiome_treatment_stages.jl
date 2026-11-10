# NEW MODEL (not from a published paper): the microbiome treatment layer on the Bsal model
# with a larval stage (bsal_larval.jl). Wiki: lab/wiki/models/new-microbiome-treatment.md
#
# Every class of bsal_larval.jl is split into untreated (U) and treated (T). State
# (S_U, S_T, L_U, L_T, I_U, I_T, W_SU, W_ST, W_IU, W_IT, E); years; per ha. Treatment acts as
# in microbiome_treatment.jl:
#   efficacy            treated adults and larvae are infected at (1 − efficacy) times the rate
#   tolerance           treated infectious adults die of Bsal at (1 − tolerance) d_I
#   shedding_reduction  treated infected adults and larvae shed (1 − s) times as much
#   reversion ω in every treated class; acquisition a_A = direct T_A/N_A + environmental E for
#   adults and a_W = direct T_W/N_W + environmental E for larvae (the stages do not share skin
#   contact); E follows the treated fraction of adults.
# Two links between the stages, both unknown:
#   retention  probability that a treated larva is still treated after metamorphosis
#   maternal   probability that the larva of a treated mother is born treated
# Larvae of untreated mothers are born untreated. Infected larvae do not die of Bsal.

function probiotic_larval_rhs!(du, u, q, t)
    p, st, tr, retention, maternal = q
    SU, ST, LU, LT, IU, IT, WSU, WST, WIU, WIT = (max(u[i], 0.0) for i in 1:10)
    E = u[11]
    m, out = 1 / st.larval_period, larval_exit(st)
    ω, δ, e, τ, s = reversion_rate(tr), reservoir_rate(tr), tr.efficacy, tr.tolerance,
        tr.shedding_reduction
    NA = SU + ST + LU + LT + IU + IT
    NW = WSU + WST + WIU + WIT
    fA = NA > 0 ? min((ST + LT + IT) / NA, 1.0) : 0.0
    fW = NW > 0 ? min((WST + WIT) / NW, 1.0) : 0.0
    aA = tr.direct * fA + tr.environmental * E
    aW = tr.direct * fW + tr.environmental * E
    shedA, shedW = IU + (1 - s) * IT, WIU + (1 - s) * WIT
    λA = p.β * shedA + st.larva_to_adult * shedW
    λW = st.adult_to_larva * shedA + st.larva_to_larva * shedW
    dd = p.d + density_dependence(p) * NA
    c, r = st.carry_over, retention
    birthsU = p.b * out / m * SU
    birthsT = p.b * out / m * (ST + τ * (LT + IT))
    # adults
    du[1] = m * (WSU + (1 - r) * WST) + (1 - c) * m * (WIU + (1 - r) * WIT) - dd * SU -
            λA * SU + ω * ST - aA * SU
    du[2] = r * m * WST + (1 - c) * r * m * WIT - dd * ST - (1 - e) * λA * ST - ω * ST + aA * SU
    du[3] = c * m * (WIU + (1 - r) * WIT) + λA * SU - (dd + p.e) * LU + ω * LT - aA * LU
    du[4] = c * r * m * WIT + (1 - e) * λA * ST - (dd + p.e) * LT - ω * LT + aA * LU
    du[5] = p.e * LU - (dd + p.d_I) * IU + ω * IT - aA * IU
    du[6] = p.e * LT - (dd + (1 - τ) * p.d_I) * IT - ω * IT + aA * IU
    # larvae
    du[7] = birthsU + (1 - maternal) * birthsT - out * WSU - λW * WSU + ω * WST - aW * WSU
    du[8] = maternal * birthsT - out * WST - (1 - e) * λW * WST - ω * WST + aW * WSU
    du[9] = λW * WSU - out * WIU + ω * WIT - aW * WIU
    du[10] = (1 - e) * λW * WST - out * WIT - ω * WIT + aW * WIU
    du[11] = δ * (fA - E)
    return nothing
end

"""
    ProbioticLarvalProblem(p::SchmidtParameters, stage::BsalLarvalStage, treatment;
        treat = :adults, coverage = 0, treat_times = [0.0], retention = 0, maternal = 0,
        t_intro = 0.0, I0 = 0.01, tspan = (0.0, 10.0))

Bsal with a larval stage ([`BsalLarvalProblem`](@ref)) and a treated class in both stages. At
each of `treat_times` a fraction `coverage` of the untreated animals in `treat` (`:adults`,
`:larvae` or `:both`) is treated. `retention` is the probability that a treated larva is still
treated after metamorphosis; `maternal` the probability that a treated mother's larva is born
treated. Bsal (`I0` infectious adults per ha) arrives at `t_intro`. Time in years.

`solve` returns series `t`, `adults`, `treated_adults`, `larvae`, `treated_larvae`,
`infectious`, and `min_adults`, `mean_adults`, `final_adults` (fractions of `K`, after
`t_intro`), `final_treated` (treated fraction of adults) and `final_prevalence`.
"""
struct ProbioticLarvalProblem{P<:SchmidtParameters}
    p::P
    stage::BsalLarvalStage
    treatment::MicrobiomeTreatment
    treat::Symbol
    coverage::Float64
    treat_times::Vector{Float64}
    retention::Float64
    maternal::Float64
    t_intro::Float64
    I0::Float64
    tspan::Tuple{Float64,Float64}
end
function ProbioticLarvalProblem(p::SchmidtParameters, stage::BsalLarvalStage,
        treatment::MicrobiomeTreatment; treat::Symbol = :adults, coverage = 0.0,
        treat_times = [0.0], retention = 0.0, maternal = 0.0, t_intro = 0.0, I0 = 0.01,
        tspan = (0.0, 10.0))
    treat in (:adults, :larvae, :both) ||
        throw(ArgumentError("treat must be :adults, :larvae or :both"))
    (0 <= coverage <= 1 && 0 <= retention <= 1 && 0 <= maternal <= 1) ||
        throw(ArgumentError("coverage, retention and maternal must be in [0, 1]"))
    return ProbioticLarvalProblem(p, stage, treatment, treat, float(coverage),
        collect(float.(treat_times)), float(retention), float(maternal), float(t_intro),
        float(I0), (float(tspan[1]), float(tspan[2])))
end

function CommonSolve.solve(prob::ProbioticLarvalProblem, alg = Tsit5(); dt = 0.01,
        abstol = 1e-10, reltol = 1e-8, kwargs...)
    p, (t0, t1) = prob.p, prob.tspan
    q = (p, prob.stage, prob.treatment, prob.retention, prob.maternal)
    events = Tuple{Float64,Symbol}[(t, :treat) for t in prob.treat_times if t0 <= t < t1]
    t0 <= prob.t_intro < t1 && push!(events, (prob.t_intro, :introduce))
    push!(events, (t1, :stop))
    sort!(events; by = e -> (e[1], e[2] === :introduce))      # treat before Bsal at a tie
    u = zeros(11)
    u[1], u[7] = p.K, disease_free_larvae(p, prob.stage)
    untreated = prob.treat === :adults ? (1, 3, 5) : prob.treat === :larvae ? (7, 9) :
        (1, 3, 5, 7, 9)
    ts, us = Float64[], Vector{Float64}[]
    t = t0
    for (te, kind) in events
        t = probiotic_advance!(u, t, te, probiotic_larval_rhs!, q, alg, ts, us; dt, abstol,
            reltol, kwargs...)
        if kind === :treat
            for i in untreated
                moved = prob.coverage * u[i]
                u[i] -= moved
                u[i + 1] += moved
            end
        elseif kind === :introduce
            u[5] += prob.I0
        end
    end
    adults = [sum(@view x[1:6]) for x in us]
    treated = [x[2] + x[4] + x[6] for x in us]
    infectious = [x[5] + x[6] for x in us]
    after = findall(>=(prob.t_intro), ts)
    rel = adults[after] ./ p.K
    return (t = ts, adults = adults, treated_adults = treated,
        larvae = [sum(@view x[7:10]) for x in us], treated_larvae = [x[8] + x[10] for x in us],
        infectious = infectious, min_adults = minimum(rel), mean_adults = sum(rel) / length(rel),
        final_adults = rel[end], final_treated = adults[end] > 0 ? treated[end] / adults[end] : 0.0,
        final_prevalence = adults[end] > 0 ? max(infectious[end], 0.0) / adults[end] : 0.0)
end
