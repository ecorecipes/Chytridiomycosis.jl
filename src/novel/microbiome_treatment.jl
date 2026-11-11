# NEW MODEL (not from a published paper): treating amphibians with a skin microbiome selected
# to produce antifungal compounds. A starter model: hosts are either "treated" (carrying the
# microbiome) or "untreated". Wiki: lab/wiki/models/new-microbiome-treatment.md
#
# The same treatment layer is put on three host–pathogen models:
#   ProbioticSLIProblem   Bsal in fire salamanders (Schmidt et al. 2017): treated adults are
#                         infected at a rate reduced by `efficacy`.            time in years
#   ProbioticPondProblem  Bd in common toads (Mitchell et al. 2008): treated tadpoles are
#                         infected at a rate reduced by `efficacy`; adults are not infected
#                         in that model and matter only as carriers.           time in years
#   ProbioticLoadProblem  Bd in frogs (the load tipping model, load_tipping.jl): treated hosts
#                         clear sporangia faster (`c_treated` in place of `c`). time in days
#
# The treatment layer, per host (rates in the host model's time unit):
#   treated -> untreated    at rate ω = log 2 / half_life
#   untreated -> treated    at rate direct × (treated fraction) + environmental × E
#   E (environmental reservoir of the microbiome, scaled so that E = 1 balances a fully
#     treated population) relaxes at rate δ = log 2 / reservoir_half_life
# Protection has two separate endpoints, because trials report better survival without fewer
# infections (lab/wiki/references/manus-microbiome-review-2026.md):
#   efficacy   proportional reduction in the infection rate of a treated host
#   tolerance  proportional reduction in the disease-induced death of a treated, infected host,
#              which stays infectious
#   shedding_reduction  proportional reduction in the pathogen shed by a treated, infected host
# A treatment that lowers pathogen load acts on the last two together
# (`load_reducing_treatment`).
# Protection can also wane with the abundance of the applied bacteria (`abundance_half_life`,
# `protective_fold`): all three effects are multiplied by `protection(treatment, age)`, where
# age is the time since the last application. Hosts that acquire the microbiome from others
# share that clock, so waning is meant for use without transmission.
# NONE of the treatment parameters has a published value here. All numbers are illustrative.

"""
    MicrobiomeTreatment(; efficacy = 0, half_life = Inf, direct = 0, environmental = 0,
                        reservoir_half_life = 1, tolerance = 0, shedding_reduction = 0,
                        abundance_half_life = Inf, protective_fold = Inf)

A protective microbiome carried by "treated" hosts. `efficacy` (0–1) is the proportional
reduction in the infection rate of a treated host. `tolerance` (0–1) is the proportional
reduction in the disease-induced death of a treated host once infected; such a host stays
infectious, shedding `1 − shedding_reduction` times as much as an untreated one. (None of
these is used by [`ProbioticLoadProblem`](@ref), which takes `c_treated`.)
A treated host reverts to untreated with half-life `half_life`. An
untreated host acquires the microbiome at rate `direct ×` (treated fraction of hosts) plus
`environmental × E`, where `E` is an environmental reservoir with half-life
`reservoir_half_life`. Times and rates are in the host model's time unit.

Waning: the applied bacteria decline with `abundance_half_life`, and `protective_fold` is the
starting abundance divided by the abundance that gives half the full effect. Efficacy,
tolerance and shedding reduction are all multiplied by [`protection`](@ref) at the time since
the last application. With the default `abundance_half_life = Inf` there is no waning.
"""
Base.@kwdef struct MicrobiomeTreatment
    efficacy::Float64 = 0.0
    half_life::Float64 = Inf
    direct::Float64 = 0.0
    environmental::Float64 = 0.0
    reservoir_half_life::Float64 = 1.0
    tolerance::Float64 = 0.0
    shedding_reduction::Float64 = 0.0
    abundance_half_life::Float64 = Inf
    protective_fold::Float64 = Inf
    function MicrobiomeTreatment(efficacy, half_life, direct, environmental,
            reservoir_half_life, tolerance, shedding_reduction, abundance_half_life,
            protective_fold)
        (0 <= efficacy <= 1 && 0 <= tolerance <= 1 && 0 <= shedding_reduction <= 1 &&
         half_life > 0 && direct >= 0 && environmental >= 0 && reservoir_half_life > 0 &&
         abundance_half_life > 0 && protective_fold > 0) ||
            throw(ArgumentError("need 0 ≤ efficacy, tolerance, shedding_reduction ≤ 1, " *
                "half-lives and protective_fold > 0 and transmission rates ≥ 0"))
        return new(efficacy, half_life, direct, environmental, reservoir_half_life, tolerance,
            shedding_reduction, abundance_half_life, protective_fold)
    end
end

"""
    protection(t::MicrobiomeTreatment, age)

Fraction of the full treatment effect that remains `age` after an application. Abundance is
`B(age) = B₀ 2^(−age/abundance_half_life)` and the effect is proportional to `B/(B + B₅₀)`,
with `B₀/B₅₀ = protective_fold`, scaled to 1 at age 0:

    protection = (1 + 1/F) / (1 + 2^(age/h)/F)

It falls to about one half after `log2(F)` abundance half-lives.
"""
function protection(t::MicrobiomeTreatment, age::Real)
    (isinf(t.abundance_half_life) || isinf(t.protective_fold)) && return 1.0
    F = t.protective_fold
    return (1 + 1 / F) / (1 + 2.0^(max(age, 0.0) / t.abundance_half_life) / F)
end

"""
    protection_duration(t::MicrobiomeTreatment; level = 0.5)

Time after an application at which [`protection`](@ref) has fallen to `level`.
"""
function protection_duration(t::MicrobiomeTreatment; level::Real = 0.5)
    (isinf(t.abundance_half_life) || isinf(t.protective_fold)) && return Inf
    F = t.protective_fold
    x = F * ((1 + 1 / F) / level - 1)
    return x > 1 ? t.abundance_half_life * log2(x) : 0.0
end

"""
    load_reducing_treatment(load_reduction; death_exponent = 1, kwargs...)

A [`MicrobiomeTreatment`](@ref) that lowers the pathogen load of a treated, infected host to
`1 − load_reduction` times the untreated load. Shedding is taken as proportional to load, and
disease-induced death as proportional to load raised to `death_exponent`:

    shedding_reduction = load_reduction
    tolerance          = 1 − (1 − load_reduction)^death_exponent

`death_exponent = 1` is death proportional to load; above 1, death falls off faster than load
(as with a lethal load threshold); below 1, more slowly. Other keywords (`efficacy`,
`half_life`, `direct`, ...) pass to `MicrobiomeTreatment`.
"""
function load_reducing_treatment(load_reduction::Real; death_exponent::Real = 1.0, kwargs...)
    (0 <= load_reduction <= 1 && death_exponent >= 0) ||
        throw(ArgumentError("need 0 ≤ load_reduction ≤ 1 and death_exponent ≥ 0"))
    return MicrobiomeTreatment(; kwargs..., shedding_reduction = load_reduction,
        tolerance = 1 - (1 - load_reduction)^death_exponent)
end

"""Rate at which a treated host reverts to untreated, `log 2 / half_life`."""
reversion_rate(t::MicrobiomeTreatment) = log(2) / t.half_life
reservoir_rate(t::MicrobiomeTreatment) = log(2) / t.reservoir_half_life

"""
    microbiome_reproduction_number(t; turnover = 0)

Number of hosts that one treated host passes the microbiome to, in an untreated population,
before it loses it or is replaced: `(direct + environmental) / (ω + turnover)`. `turnover`
is the hosts' replacement rate (births of untreated hosts balancing deaths). Above 1 the
microbiome maintains itself without further treatment.
"""
function microbiome_reproduction_number(t::MicrobiomeTreatment; turnover::Real = 0.0)
    loss = reversion_rate(t) + turnover
    gain = t.direct + t.environmental
    return iszero(loss) ? (iszero(gain) ? 1.0 : Inf) : gain / loss
end

"""
    equilibrium_coverage(t; turnover = 0)

Long-run treated fraction in a host population that is always present, with no further
treatment: `1 − 1/R` where `R` is [`microbiome_reproduction_number`](@ref), or 0 if `R ≤ 1`.
Direct and environmental transmission count equally here; they differ only in timing.
"""
function equilibrium_coverage(t::MicrobiomeTreatment; turnover::Real = 0.0)
    R = microbiome_reproduction_number(t; turnover)
    return R > 1 ? 1 - 1 / R : 0.0
end

"""
    treatment_efficacy(treated, untreated, disease_free)

Share of the disease-induced loss that treatment averts, for any outcome measured in the
three situations: `(treated − untreated)/(disease_free − untreated)`.
"""
treatment_efficacy(treated, untreated, disease_free) =
    (treated - untreated) / (disease_free - untreated)

"""
    mortality_reduction(deaths_treated, n_treated, deaths_control, n_control;
                        draws = 20_000, rng)

Proportional reduction in the risk of death in a two-arm challenge trial,
`1 − p_treated/p_control`, as `(median, lower, upper)` with a 95% interval. Each arm's risk has
a Beta(½, ½) prior. A statistical summary of one trial, not a parameter for another host.
"""
function mortality_reduction(deaths_treated::Real, n_treated::Real, deaths_control::Real,
        n_control::Real; draws::Integer = 20_000, rng::AbstractRNG = Random.default_rng())
    (0 <= deaths_treated <= n_treated && 0 <= deaths_control <= n_control) ||
        throw(ArgumentError("deaths must lie between 0 and the group size"))
    pt = rand(rng, Beta(deaths_treated + 0.5, n_treated - deaths_treated + 0.5), draws)
    pc = rand(rng, Beta(deaths_control + 0.5, n_control - deaths_control + 0.5), draws)
    red = sort!(1 .- pt ./ pc)
    q(x) = red[clamp(round(Int, x * draws), 1, draws)]
    return (median = q(0.5), lower = q(0.025), upper = q(0.975))
end

# Integrate `rhs!` from t to tnext, appending the saved states to `ts`, `us`.
function probiotic_advance!(u, t, tnext, rhs!, q, alg, ts, us; dt, kwargs...)
    tnext > t || return t
    sol = solve(ODEProblem(rhs!, copy(u), (t, tnext), q), alg; saveat = dt, kwargs...)
    SciMLBase.successful_retcode(sol) || error("treatment model: solver failed after t = $t")
    if ts !== nothing
        append!(ts, sol.t)
        append!(us, sol.u)
    end
    u .= max.(sol.u[end], 0.0)
    return tnext
end

# ---------------------------------------------------------------------------------------
# Bsal in fire salamanders: Schmidt et al. (2017) S–L–I with every class split into
# untreated (U) and treated (T). State (S_U, S_T, L_U, L_T, I_U, I_T, E); years; adult
# females per ha. With I = I_U + (1 − s) I_T (s = shedding_reduction), N all hosts,
# ε_ = efficacy, τ = tolerance:
#   dS_U/dt = B − d S_U − γ S_U N − β S_U I            + ω S_T − a S_U
#   dS_T/dt =   − d S_T − γ S_T N − (1 − ε_) β S_T I   − ω S_T + a S_U
#   dL_U/dt = β S_U I − (d + e) L_U − γ L_U N          + ω L_T − a L_U
#   dL_T/dt = (1 − ε_) β S_T I − (d + e) L_T − γ L_T N − ω L_T + a L_U
#   dI_U/dt = e L_U − (d + d_I) I_U − γ I_U N          + ω I_T − a I_U
#   dI_T/dt = e L_T − (d + (1 − τ) d_I) I_T − γ I_T N  − ω I_T + a I_U
#   dE/dt   = δ (T/N − E),   T = S_T + L_T + I_T,   a = direct T/N + environmental E
#   B = b (S_U + S_T + τ (L_T + I_T))
# Recruits are untreated. In Schmidt's model only susceptibles breed; here a treated infected
# host also breeds, in proportion to its tolerance (an assumption: with τ = 0 this is Schmidt's
# rule, with τ = 1 a carrier breeds like a healthy animal). Larvae are not in the model.

function probiotic_sli_rhs!(du, u, q, t)
    p, tr, clock = q
    SU, ST, LU, LT, IU, IT = (max(u[i], 0.0) for i in 1:6)
    E = u[7]
    g = protection(tr, t - clock[])
    N = SU + ST + LU + LT + IU + IT
    I = IU + (1 - g * tr.shedding_reduction) * IT
    γ = density_dependence(p)
    ω, δ, e, τ = reversion_rate(tr), reservoir_rate(tr), g * tr.efficacy, g * tr.tolerance
    frac = N > 0 ? min((ST + LT + IT) / N, 1.0) : 0.0
    a = tr.direct * frac + tr.environmental * E
    out = p.d + γ * N
    du[1] = p.b * (SU + ST + τ * (LT + IT)) - out * SU - p.β * SU * I + ω * ST - a * SU
    du[2] = -out * ST - (1 - e) * p.β * ST * I - ω * ST + a * SU
    du[3] = p.β * SU * I - (out + p.e) * LU + ω * LT - a * LU
    du[4] = (1 - e) * p.β * ST * I - (out + p.e) * LT - ω * LT + a * LU
    du[5] = p.e * LU - (out + p.d_I) * IU + ω * IT - a * IU
    du[6] = p.e * LT - (out + (1 - τ) * p.d_I) * IT - ω * IT + a * IU
    du[7] = δ * (frac - E)
    return nothing
end

"""
    ProbioticSLIProblem(p::SchmidtParameters, treatment; coverage = 0, treat_times = [0.0],
                        t_intro = 0.0, I0 = 0.01, tspan = (0.0, 10.0))

Schmidt et al. (2017) Bsal model with a treated class. The population starts at `K`. At each
of `treat_times` a fraction `coverage` of untreated hosts is treated; Bsal (`I0`
infectious per ha) arrives at `t_intro`. Time in years.

`solve` returns series `t`, `hosts`, `treated`, `infectious` and the summaries `min_hosts`,
`mean_hosts` and `final_hosts` after `t_intro`, as fractions of `K`, and `final_prevalence`
(the infectious fraction of hosts at the end).
"""
struct ProbioticSLIProblem{P<:SchmidtParameters}
    p::P
    treatment::MicrobiomeTreatment
    coverage::Float64
    treat_times::Vector{Float64}
    t_intro::Float64
    I0::Float64
    tspan::Tuple{Float64,Float64}
end
function ProbioticSLIProblem(p::SchmidtParameters, treatment::MicrobiomeTreatment;
        coverage = 0.0, treat_times = [0.0], t_intro = 0.0, I0 = 0.01, tspan = (0.0, 10.0))
    0 <= coverage <= 1 || throw(ArgumentError("coverage must be in [0, 1]"))
    return ProbioticSLIProblem(p, treatment, float(coverage), collect(float.(treat_times)),
        float(t_intro), float(I0), (float(tspan[1]), float(tspan[2])))
end

function CommonSolve.solve(prob::ProbioticSLIProblem, alg = Tsit5(); dt = 0.01,
        abstol = 1e-10, reltol = 1e-8, kwargs...)
    p, (t0, t1) = prob.p, prob.tspan
    events = Tuple{Float64,Symbol}[(t, :treat) for t in prob.treat_times if t0 <= t < t1]
    t0 <= prob.t_intro < t1 && push!(events, (prob.t_intro, :introduce))
    push!(events, (t1, :stop))
    sort!(events; by = e -> (e[1], e[2] === :introduce))      # treat before Bsal at a tie
    u = [float(p.K), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]   # SU, ST, LU, LT, IU, IT, E
    ts, us = Float64[], Vector{Float64}[]
    t = t0
    clock = Ref(t0)                                  # time of the last application
    for (te, kind) in events
        t = probiotic_advance!(u, t, te, probiotic_sli_rhs!, (p, prob.treatment, clock), alg,
            ts, us; dt, abstol, reltol, kwargs...)
        if kind === :treat
            clock[] = te
            for i in (1, 3, 5)
                moved = prob.coverage * u[i]
                u[i] -= moved
                u[i + 1] += moved
            end
        elseif kind === :introduce
            u[5] += prob.I0
        end
    end
    hosts = [sum(@view x[1:6]) for x in us]
    infectious = [x[5] + x[6] for x in us]
    after = findall(>=(prob.t_intro), ts)
    rel = hosts[after] ./ p.K
    span = ts[after[end]] - ts[after[1]]
    mean_hosts = span > 0 ? sum((rel[i] + rel[i + 1]) / 2 * (ts[after[i + 1]] - ts[after[i]])
        for i in 1:(length(after) - 1)) / span : rel[end]
    return (t = ts, hosts = hosts, treated = [x[2] + x[4] + x[6] for x in us],
        infectious = infectious, min_hosts = minimum(rel), mean_hosts = mean_hosts,
        final_hosts = rel[end],
        final_prevalence = hosts[end] > 0 ? max(infectious[end], 0.0) / hosts[end] : 0.0)
end

# Relative transmission from a treated host: it is infected (1 − efficacy) times as often,
# sheds (1 − shedding_reduction) times as much, and is infectious for longer, by
# (b + d_I)/(b + (1 − tolerance) d_I) at the disease-free state (where d + γK = b).
function relative_transmission(p::SchmidtParameters, t::MicrobiomeTreatment)
    return (1 - t.efficacy) * (1 - t.shedding_reduction) * (p.b + p.d_I) /
           (p.b + (1 - t.tolerance) * p.d_I)
end

"""
    effective_reproduction_number(p::SchmidtParameters, t::MicrobiomeTreatment, coverage)

Bsal's reproduction number when a fraction `coverage` of hosts is treated:
`R₀ (1 − coverage + coverage × m)`, with
`m = (1 − efficacy)(1 − shedding_reduction)(b + d_I)/(b + (1 − tolerance) d_I)` the
transmission from a treated host relative to an untreated one. It ignores gain and loss of
the microbiome during an infection. Tolerance makes `m` larger and a shedding reduction makes
it smaller; `m > 1` means that treatment raises the reproduction number.
"""
effective_reproduction_number(p::SchmidtParameters, t::MicrobiomeTreatment, coverage::Real) =
    basic_reproduction_number(p) * (1 - coverage + coverage * relative_transmission(p, t))

"""
    critical_coverage(p::SchmidtParameters, t::MicrobiomeTreatment)

Treated fraction needed to stop Bsal invading: `(1 − 1/R₀)/(1 − m)`, with `m` as in
[`effective_reproduction_number`](@ref). A value above 1 means no coverage is enough; `Inf`
means treatment does not lower the reproduction number at all.
"""
function critical_coverage(p::SchmidtParameters, t::MicrobiomeTreatment)
    gain = 1 - relative_transmission(p, t)
    return gain > 0 ? max(0.0, 1 - 1 / basic_reproduction_number(p)) / gain : Inf
end

# ---------------------------------------------------------------------------------------
# Bd in common toads: Mitchell et al. (2008) with treated tadpoles and treated adults.
# State (P_XU, P_XT, P_YU, M, A_U, A_T, Z, E, P_YT); years. Additions to equations (2.1)–(2.5):
#   treated tadpoles are infected at rate (1 − efficacy) ν Z
#   treated infected tadpoles (P_YT) die of disease at (1 − tolerance) α_P, shed zoospores at
#     (1 − shedding_reduction) ρ, and a fraction `tolerance` of them survives metamorphosis
#     (untreated infected metamorphs all die)
#   tadpoles: treated -> untreated at ω; untreated -> treated at a = direct P_T/P + environmental E
#   adults:   treated -> untreated at ω (adults are not infected in this model)
#   dE/dt = δ (P_T/P_ref − E), with P_T = P_XT + P_YT and P_ref = ιψK, the largest yearly cohort
#   at breeding: a fraction `maternal × A_T/A` of tadpoles is born treated, and the reservoir
#                gains `adult_seeding × A_T/A`
# Metamorphs lose the "treated" label (the skin community is assumed to turn over at
# metamorphosis), so new adults are untreated.

function probiotic_pond_rhs!(du, u, q, t)
    p, tr, Pref, clock = q
    PXU, PXT, PYU, M, AU, AT, Z, E, PYT = (max(x, 0.0) for x in u)
    g = protection(tr, t - clock[])
    ω, δ, e, τ = reversion_rate(tr), reservoir_rate(tr), g * tr.efficacy, g * tr.tolerance
    P = PXU + PXT + PYU + PYT
    PT = PXT + PYT
    a = (P > 0 ? tr.direct * min(PT / P, 1.0) : 0.0) + tr.environmental * E
    du[1] = -(p.μP + p.ν * Z) * PXU + ω * PXT - a * PXU
    du[2] = -(p.μP + (1 - e) * p.ν * Z) * PXT - ω * PXT + a * PXU
    du[3] = p.ν * Z * PXU - (p.μP + p.αP) * PYU + ω * PYT - a * PYU
    du[4] = -(p.μJ + p.ε) * M
    du[5] = p.ε * M * (1 - (AU + AT) / p.K) - p.μA * AU + ω * AT
    du[6] = -(p.μA + ω) * AT
    du[7] = p.ρ * (PYU + (1 - g * tr.shedding_reduction) * PYT) + p.γ * Z * exp(-p.φ * Z) - p.μZ * Z
    du[8] = δ * (PT / Pref - E)
    du[9] = (1 - e) * p.ν * Z * PXT - (p.μP + (1 - τ) * p.αP) * PYT - ω * PYT + a * PYU
    return nothing
end

"""
    ProbioticPondProblem(p::MitchellParameters, treatment; stage = :tadpoles, coverage = 0,
        treat_years = 0:(years - 1), maternal = 0, adult_seeding = 0, years = 35,
        spinup = 100, Z0 = 1e4, A0 = p.K / 2, toad_limit = 50, bd_limit = 100,
        stages = PostMetamorphicInfection(), redose_days = [])

Mitchell et al. (2008) toad–Bd model with a treated class. In each year of `treat_years`
(year 0 is the year Bd arrives) a fraction `coverage` is treated: of the newly hatched
tadpoles if `stage = :tadpoles`, or of the adults just before breeding if `stage = :adults`.
Treated adults pass the microbiome on directly (`maternal`: probability that a treated
adult's tadpole is born treated) or through the pond (`adult_seeding`: reservoir added at
breeding when all adults are treated). Time in years.

`stages` adds infection after metamorphosis ([`PostMetamorphicInfection`](@ref)): an infected
metamorph survives with probability `s_J = metamorph_survival`, or `s_J + tolerance (1 − s_J)`
if treated, and metamorphs leaving the pond are exposed for `metamorph_exposure_days`
(treated ones at `1 − efficacy` times the rate). Infected adults are not supported here yet.

`redose_days` (for `stage = :tadpoles`) are days after hatching at which, in each treatment
year, treated tadpoles are dosed again and a fraction `coverage` of the untreated ones is
treated. Each application restarts the waning clock of [`protection`](@ref).

`solve` returns series at each event (`t`, `A`, `A_treated`, `Z`, `tadpoles_treated` as a
fraction of uninfected tadpoles), and `toads_persist`, `bd_persist`, `A_disease_free`,
`A_final` and `A_mean` (adults averaged over the last 10 years).
"""
struct ProbioticPondProblem
    p::MitchellParameters
    treatment::MicrobiomeTreatment
    stage::Symbol
    coverage::Float64
    treat_years::Vector{Int}
    maternal::Float64
    adult_seeding::Float64
    years::Int
    spinup::Int
    Z0::Float64
    A0::Float64
    toad_limit::Float64
    bd_limit::Float64
    stages::PostMetamorphicInfection
    redose_days::Vector{Float64}
end
function ProbioticPondProblem(p::MitchellParameters, treatment::MicrobiomeTreatment;
        stage::Symbol = :tadpoles, coverage = 0.0, years::Integer = 35,
        treat_years = 0:(years - 1), maternal = 0.0, adult_seeding = 0.0,
        spinup::Integer = 100, Z0 = 1e4, A0 = p.K / 2, toad_limit = 50.0, bd_limit = 100.0,
        stages::PostMetamorphicInfection = PostMetamorphicInfection(), redose_days = Float64[])
    stages.adult_exposure_days == 0 || throw(ArgumentError(
        "infected adults (adult_exposure_days > 0) are not supported with treatment yet"))
    stage in (:tadpoles, :adults) || throw(ArgumentError("stage must be :tadpoles or :adults"))
    (0 <= coverage <= 1 && 0 <= maternal <= 1 && adult_seeding >= 0) ||
        throw(ArgumentError("need coverage, maternal in [0, 1] and adult_seeding ≥ 0"))
    return ProbioticPondProblem(p, treatment, stage, float(coverage), collect(Int, treat_years),
        float(maternal), float(adult_seeding), Int(years), Int(spinup), float(Z0), float(A0),
        float(toad_limit), float(bd_limit), stages, collect(Float64, redose_days))
end

function CommonSolve.solve(prob::ProbioticPondProblem, alg = Tsit5(); abstol = 1e-6,
        reltol = 1e-8, kwargs...)
    p, tr = prob.p, prob.treatment
    0 <= p.t1 < p.t2 < 1 || throw(ArgumentError("need 0 ≤ t1 < t2 < 1 (fractions of a year)"))
    clock = Ref(float(-prob.spinup))
    q = (p, tr, p.ι * p.ψ * p.K, clock)
    u = [0.0, 0.0, 0.0, 0.0, prob.A0, 0.0, 0.0, 0.0, 0.0]   # PXU, PXT, PYU, M, AU, AT, Z, E, PYT
    cohorts = Tuple{Float64,Float64}[]
    events = Tuple{Float64,Symbol,Int}[]
    for y in -prob.spinup:(prob.years - 1)
        push!(events, (y + p.t1, :birth, y), (y + p.t2, :metamorphosis, y))
        if prob.stage === :tadpoles && y in prob.treat_years
            for d in prob.redose_days
                0 < d / 365 < p.t2 - p.t1 || throw(ArgumentError(
                    "redose_days must fall inside the tadpole season"))
                push!(events, (y + p.t1 + d / 365, :redose, y))
            end
        end
    end
    push!(events, (0.0, :introduce, 0), (float(prob.years), :stop, prob.years))
    sort!(events; by = first)
    out = (t = Float64[], A = Float64[], A_treated = Float64[], Z = Float64[],
        tadpoles_treated = Float64[])
    step!(t, tn) = probiotic_advance!(u, t, tn, probiotic_pond_rhs!, q, alg, nothing, nothing;
        dt = Float64[], save_everystep = false, abstol, reltol, kwargs...)
    t = float(-prob.spinup)
    A_disease_free = NaN
    for (te, kind, y) in events
        while !isempty(cohorts) && first(cohorts)[1] <= te
            tm, size = popfirst!(cohorts)
            t = step!(t, tm)
            u[4] += size
        end
        t = step!(t, te)
        treating = y in prob.treat_years
        if kind === :birth
            0 < u[7] < prob.bd_limit && (u[7] = 0.0)
            treating && (clock[] = te)
            if treating && prob.stage === :adults
                moved = prob.coverage * u[5]
                u[5] -= moved
                u[6] += moved
            end
            A = u[5] + u[6]
            carried = A > 0 ? u[6] / A : 0.0
            births = p.ι * p.ψ * A
            u[1] += births * (1 - prob.maternal * carried)
            u[2] += births * prob.maternal * carried
            u[8] += prob.adult_seeding * carried
            if treating && prob.stage === :tadpoles
                moved = prob.coverage * u[1]
                u[1] -= moved
                u[2] += moved
            end
        elseif kind === :redose
            clock[] = te
            for (i, j) in ((1, 2), (3, 9))
                moved = prob.coverage * u[i]
                u[i] -= moved
                u[j] += moved
            end
        elseif kind === :metamorphosis
            st = prob.stages
            g = protection(tr, te - clock[])
            exposure = p.ν * u[7] * st.metamorph_exposure_days / 365
            qU, qT = -expm1(-exposure), -expm1(-(1 - g * tr.efficacy) * exposure)
            sU = st.metamorph_survival
            sT = sU + g * tr.tolerance * (1 - sU)
            survivors = u[1] * (1 - qU + qU * sU) + u[2] * (1 - qT + qT * sT) + u[3] * sU + u[9] * sT
            push!(cohorts, (te + p.ω, survivors * exp(-p.μJ * p.ω)))
            u[1] = u[2] = u[3] = u[9] = 0.0
        elseif kind === :introduce
            A_disease_free = u[5] + u[6]
            u[7] += prob.Z0
        end
        if 0 < u[5] + u[6] < prob.toad_limit
            u[1] = u[2] = u[3] = u[4] = u[5] = u[6] = u[9] = 0.0
            empty!(cohorts)
        end
        push!(out.t, te); push!(out.A, u[5] + u[6]); push!(out.A_treated, u[6])
        push!(out.Z, u[7])
        push!(out.tadpoles_treated, u[1] + u[2] > 0 ? u[2] / (u[1] + u[2]) : 0.0)
    end
    late = findall(>=(prob.years - 10), out.t)
    return merge(out, (toads_persist = u[5] + u[6] > 0, bd_persist = u[7] >= prob.bd_limit,
        A_disease_free = A_disease_free, A_final = u[5] + u[6],
        A_mean = sum(out.A[late]) / length(late)))
end

# ---------------------------------------------------------------------------------------
# Bd in frogs: the load tipping model with treated hosts. A treated host has clearance
# `c_treated` in place of `c`. Loads follow the tipping ODE within each day; at the end of
# each day hosts gain or lose the microbiome at random (so runs differ; pass `rng`).

"""
    critical_coverage(p::BriggsParameters; c, c_treated, h, N, ε = 0)

Smallest treated fraction of `N` hosts at which a low-load state exists (the population is
not in the `:runaway` regime of the load tipping model). `NaN` if even full coverage is not
enough.
"""
function critical_coverage(p::BriggsParameters; c, c_treated, h, N::Integer, ε = 0.0)
    safe(k) = pooled_load_equilibria(p; c = vcat(fill(float(c_treated), k),
        fill(float(c), N - k)), h, ε).regime !== :runaway
    safe(0) && return 0.0
    safe(N) || return NaN
    lo, hi = 0, N
    while hi - lo > 1
        mid = (lo + hi) ÷ 2
        safe(mid) ? (hi = mid) : (lo = mid)
    end
    return hi / N
end

"""
    ProbioticLoadProblem(p::BriggsParameters, treatment; c, c_treated, h, N, ε = 0,
                         coverage = 0, treat_days = [0], days = 365)

Load tipping model for `N` frogs sharing a pool, starting uninfected, in which a treated
frog has clearance `c_treated` instead of `c`. On each of `treat_days` every untreated frog
is treated with probability `coverage`. Time in days (so `treatment` rates are per day).

`solve(prob; rng)` returns daily `alive`, `treated` and `mean_load`, and `death_day` for each
frog (`0` for survivors).
"""
struct ProbioticLoadProblem{P<:BriggsParameters}
    p::P
    treatment::MicrobiomeTreatment
    c::Float64
    c_treated::Float64
    h::Float64
    N::Int
    ε::Float64
    coverage::Float64
    treat_days::Vector{Int}
    days::Int
end
function ProbioticLoadProblem(p::BriggsParameters, treatment::MicrobiomeTreatment; c,
        c_treated, h, N::Integer, ε = 0.0, coverage = 0.0, treat_days = [0],
        days::Integer = 365)
    (c >= 0 && c_treated >= 0 && h > 0 && ε >= 0 && 0 <= coverage <= 1 && N > 0) ||
        throw(ArgumentError("need c, c_treated, ε ≥ 0, h > 0, N > 0, coverage in [0, 1]"))
    return ProbioticLoadProblem(p, treatment, float(c), float(c_treated), float(h), Int(N),
        float(ε), float(coverage), collect(Int, treat_days), Int(days))
end

function CommonSolve.solve(prob::ProbioticLoadProblem, alg = Tsit5();
        rng::AbstractRNG = Random.default_rng(), abstol = 1e-6, reltol = 1e-6, kwargs...)
    p, tr = prob.p, prob.treatment
    ids = collect(1:prob.N)                     # frogs still alive
    S = zeros(prob.N)
    treated = falses(prob.N)
    Z = prob.ε / (p.μ + p.γ * prob.N / p.V)
    E = 0.0
    revert = -expm1(-reversion_rate(tr))        # daily probability of losing the microbiome
    decay = exp(-reservoir_rate(tr))
    out = (alive = zeros(Int, prob.days), treated = zeros(Int, prob.days),
        mean_load = fill(NaN, prob.days))
    death_day = zeros(Int, prob.N)
    for day in 1:prob.days
        n = length(ids)
        n == 0 && break
        if (day - 1) in prob.treat_days
            for i in 1:n
                treated[i] || (treated[i] = rand(rng) < prob.coverage)
            end
        end
        c = [treated[i] ? prob.c_treated : prob.c for i in 1:n]
        r = solve(LoadTippingProblem(p; c, h = prob.h, N = n, ε = prob.ε, S0 = S, Z0 = Z,
            tspan = (0.0, 1.0)), alg; abstol, reltol, save_everystep = false, kwargs...)
        SciMLBase.successful_retcode(r.sol) || error("treated load model: solver failed on day $day")
        final = r.sol.u[end]
        Z = max(final[end], 0.0)
        for i in 1:n
            r.alive[i] || (death_day[ids[i]] = day)
        end
        frac = count(treated) / n               # treated fraction during the day
        ids, S, treated = ids[r.alive], max.(final[1:n][r.alive], 0.0), treated[r.alive]
        E = frac + (E - frac) * decay
        acquire = -expm1(-(tr.direct * frac + tr.environmental * E))
        for i in eachindex(ids)
            treated[i] = treated[i] ? rand(rng) >= revert : rand(rng) < acquire
        end
        out.alive[day] = length(ids)
        out.treated[day] = count(treated)
        isempty(ids) || (out.mean_load[day] = sum(S) / length(ids))
    end
    return merge(out, (death_day = death_day,))
end
