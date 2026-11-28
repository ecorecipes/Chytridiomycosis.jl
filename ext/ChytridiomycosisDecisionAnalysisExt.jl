module ChytridiomycosisDecisionAnalysisExt

# The adapter between this package and DecisionAnalysis.jl: methods of its adapter verbs for
# Drawert et al. (2017) (`DrawertIBMProblem`, `DrawertODEProblem`), Canessa et al. (2018)
# (`CanessaIPMProblem`, `CanessaSIProblem`, `canessa_outcome_table`), Cook et al. (2022)
# (`CookProjectionProblem`), Heard et al. (2024) (`HeardMetapopulationProblem`), Turner et al.
# (2025) (`TurnerPVAProblem`), and the Drawert swab/eDNA surveillance design
# (`drawert_survey_study`), so that
#
#     outcomes(prob, actions, draws; replicates, seed)
#
# runs each model over parameter draws × actions × replicates with common random numbers.
# The problems, interventions and solvers are this package's own and need no decision package;
# the methods here only name them in DecisionAnalysis's vocabulary. Outcomes are counts and
# epidemiological indicators: utilities are left to the user.

using Chytridiomycosis: Chytridiomycosis, CANESSA_ACTION_NAMES, CanessaAction, CanessaDraw,
                        CanessaIPMProblem, CanessaParameters, CanessaSIParameters,
                        CanessaSIProblem, CookParameters, CookProjectionProblem,
                        CookProjectionSolution, DrawertIBMProblem, DrawertIntervention,
                        DrawertODEProblem, DrawertParameters, HeardMetapopulationProblem,
                        HeardMetapopulationSolution, HeardParameters, SurveyObservation,
                        SurveyPond, SurveySample, TURNER_MONTHS, TURNER_PATTERNS,
                        TurnerPVAProblem, TurnerPVASolution, TurnerParameters,
                        basic_reproduction_number, canessa_action, canessa_ensemble,
                        drawert_run_seed, heard2024_landscapes, load_canessa2018,
                        load_drawert_decision_draws, load_drawert_decisions, load_turner2025,
                        simulate_survey, solve, survey_loglikelihood
using Random: AbstractRNG, Xoshiro
import DecisionAnalysis

# ---------------------------------------------------------------------------------------------
# Shared parameter-update helper

function _update_struct(::Type{P}, p::P, θ; aliases = NamedTuple()) where {P}
    nt = NamedTuple(θ)
    fields = fieldnames(P)
    resolved = Dict{Symbol,Any}()
    for k in keys(nt)
        target = haskey(aliases, k) ? getfield(aliases, k) : k
        if target in fields
            resolved[target] = getfield(nt, k)
        end
    end
    isempty(resolved) &&
        throw(ArgumentError("the draw names none of the parameters $fields: $(keys(nt))"))
    vals = ntuple(length(fields)) do i
        k = fields[i]
        if haskey(resolved, k)
            ft = fieldtype(P, k)
            v = resolved[k]
            v isa ft ? v : convert(ft, v)
        else
            getfield(p, k)
        end
    end
    return P(vals...), resolved
end

# ---------------------------------------------------------------------------------------------
# Drawert et al. (2017): DrawertIBMProblem & DrawertODEProblem

# `p` with the fields named in the draw `θ` replaced. Columns of `θ` that are not parameters
# (an id, a group, a weight) are ignored, but a draw that names no parameter is a mistake.
function _with(p::DrawertParameters, θ)
    nt = NamedTuple(θ)
    fields = fieldnames(DrawertParameters)
    used = filter(in(fields), keys(nt))
    isempty(used) &&
        throw(ArgumentError("the draw names none of the parameters $fields: $(keys(nt))"))
    value(k) = k in used ? convert(fieldtype(DrawertParameters, k), getfield(nt, k)) : getfield(p, k)
    q = DrawertParameters(; (k => value(k) for k in fields)...)
    Chytridiomycosis._check(q)
    return q
end

"""
    DecisionAnalysis.with_parameters(prob::DrawertIBMProblem, θ)
    DecisionAnalysis.with_parameters(prob::DrawertODEProblem, θ)

The problem with the parameters named in the draw `θ` (a named tuple, one row of the parameter
draws) in place of its own. Any field of `DrawertParameters` can be a column (`N`, `γ`, `η`,
`ν0`, `f`, `σ0`, `μ`, and also `ν1`, `σ1`, `Smax`, `S0`, `Z0`, `V`); the others keep the values
of `prob`. Columns that are not parameters (`draw_id`, `group`, `weight`, …) are ignored, so a
row of `drawert2017_draws()` or `load_drawert_decision_draws()` can be passed as it is; a draw
that names no parameter throws.

`ν1` and `σ1` are **not recomputed**: in the Nu and Sigma models of the paper they depend on the
other parameters (`drawert_variant`), so give them as columns of the draws, computed for each
draw, when the problem is not the Baseline model.
"""
DecisionAnalysis.with_parameters(prob::DrawertIBMProblem, θ) =
    DrawertIBMProblem(_with(prob.p, θ), prob.interventions, prob.tfinal, prob.saveat,
                      prob.save_loads, prob.as_coded, prob.evaluate, prob.infection_threshold)
DecisionAnalysis.with_parameters(prob::DrawertODEProblem, θ) =
    DrawertODEProblem(_with(prob.p, θ); prob.intervention, prob.tfinal)

_pulse(iv::DrawertIntervention) = iv
_pulse(a::NamedTuple) = DrawertIntervention(a.kind, a.day, a.fraction; mortality = get(a, :mortality, 0))
_pulse(a) = throw(ArgumentError("an action is `nothing`, a DrawertIntervention, a named tuple " *
                                "(; kind, day, fraction), or a list of those; got $(typeof(a))"))

_pulses(::Nothing) = DrawertIntervention[]
_pulses(a::Union{AbstractVector,Tuple}) = DrawertIntervention[_pulse(x) for x in a]
_pulses(a) = [_pulse(a)]

"""
    DecisionAnalysis.apply_action(prob::DrawertIBMProblem, action)
    DecisionAnalysis.apply_action(prob::DrawertODEProblem, action)

The problem under `action`, which **replaces** any intervention of `prob`:

- `nothing` (or an empty list): no intervention, the baseline;
- a `DrawertIntervention`, or a named tuple `(; kind, day, fraction)` with an optional
  `mortality`, the arguments of `DrawertIntervention`;
- a vector or tuple of those: several pulses, applied in order of day. Several pulses, and
  `mortality`, are extensions of this package, not in the paper.

The deterministic problem takes at most one pulse, without `mortality`.
"""
function DecisionAnalysis.apply_action(prob::DrawertIBMProblem, action)
    ivs = sort!(_pulses(action); by = iv -> iv.day, alg = MergeSort)
    return DrawertIBMProblem(prob.p, ivs, prob.tfinal, prob.saveat, prob.save_loads, prob.as_coded,
                             prob.evaluate, prob.infection_threshold)
end

function DecisionAnalysis.apply_action(prob::DrawertODEProblem, action)
    ivs = _pulses(action)
    length(ivs) <= 1 ||
        throw(ArgumentError("the deterministic Drawert problem takes at most one pulse"))
    return DrawertODEProblem(prob.p; intervention = isempty(ivs) ? nothing : only(ivs), prob.tfinal)
end

"""
    DecisionAnalysis.simulate(prob::DrawertIBMProblem; rng)
    DecisionAnalysis.simulate(prob::DrawertODEProblem; rng)

One run, `solve(prob; rng)`. The individual-based model draws every random number from `rng`,
in event order: runs of two actions with the same `rng` are identical up to the first pulse on
which they differ and diverge after it, since there are no per-frog streams. The deterministic
model does not touch `rng`.
"""
DecisionAnalysis.simulate(prob::DrawertIBMProblem; rng) = solve(prob; rng)
DecisionAnalysis.simulate(prob::DrawertODEProblem; rng) = solve(prob)

_day_label(d::Real) = isinteger(d) ? string(Int(d)) : replace(string(d), "." => "_")

_classes(outcome::Symbol) = (; frogs_extinct = Int(outcome === :frogs_extinct),
                             fungus_extinct = Int(outcome === :fungus_extinct),
                             persistence = Int(outcome === :persistence))

"""
    DecisionAnalysis.outcome(prob::DrawertIBMProblem, run)
    DecisionAnalysis.outcome(prob::DrawertODEProblem, run)

The named outcomes of one run, all counts:

- `initial_frogs` (`N`), `survivors` (alive at the end; culled frogs do not survive);
- `frogs_extinct`, `fungus_extinct`, `persistence`: 1 for the outcome class of the run, 0 for
  the others, so that their means are probabilities.

For the individual-based model also:

- `uninfected`, `infected`: survivors at or below, and above, the infection threshold of the
  problem (`survivors = uninfected + infected`);
- `handled` (frogs cleaned or culled) and `treatment_deaths`, over all pulses;
- `uninfected_d<day>`, `infected_d<day>` for each day in `prob.evaluate`.

The methods take the problem because a run is a named tuple, not a type of this package.
"""
function DecisionAnalysis.outcome(prob::DrawertIBMProblem, run)
    base = (; initial_frogs = prob.p.N, run.survivors, run.uninfected, run.infected, run.handled,
            run.treatment_deaths, _classes(run.outcome)...)
    days = NamedTuple(Symbol(name, "_d", _day_label(e.day)) => getfield(e, name)
                      for e in run.evaluation for name in (:uninfected, :infected))
    return merge(base, days)
end

DecisionAnalysis.outcome(prob::DrawertODEProblem, run) =
    (; initial_frogs = prob.p.N, run.survivors, _classes(run.outcome)...)

# Stored Drawert table

const _STORED_OUTCOMES = (:survivors, :uninfected, :infected, :uninfected_d120, :infected_d120,
                          :uninfected_d365, :infected_d365, :initial_frogs, :handled,
                          :treatment_deaths)

function Chytridiomycosis.drawert_outcome_table(; variant::AbstractString = "baseline",
        outcomes = _STORED_OUTCOMES, actions = nothing)
    runs = filter(r -> r.variant == variant, load_drawert_decisions())
    isempty(runs) && throw(ArgumentError("the stored table has no runs of the variant $variant"))
    draws = load_drawert_decision_draws()
    # Replicate k is labelled r<k>, as the driver labels it, so that new runs line up.
    stored = DecisionAnalysis.OutcomeTable([(; r..., replicate = string("r", r.replicate)) for r in runs];
                                           draw = :draw_id, action = :action, replicate = :replicate,
                                           outcomes = collect(Symbol, outcomes), status = :status,
                                           parameters = draws, key = :draw_id, weights = :weight)
    (actions === nothing || isempty(actions)) && return stored
    variant == "baseline" ||
        throw(ArgumentError("only the Baseline model can be extended with new actions; the " *
                            "variant $variant has parameters that depend on the others"))
    new = DecisionAnalysis.outcomes(_stored_problem(first(draws)), actions, _stored_belief(draws);
                                    replicates = maximum(r.replicate for r in runs),
                                    run_seed = drawert_run_seed)
    return DecisionAnalysis.combine_actions(stored, new)
end

# The problem of the stored runs (scripts/drawert2017_decision_table.jl): a year, the load
# threshold 0, uninfected and infected frogs on days 120 and 365.
_stored_problem(draw) = DrawertIBMProblem(DrawertParameters(draw);
                                          evaluate = Chytridiomycosis._DRAWERT_EVALUATION_DAYS)

# The stored draws as a belief: labelled by `draw_id` and still carrying it, for the seeds.
_stored_belief(draws) = DecisionAnalysis.SampleBelief(
    Base.structdiff(_columns(draws), NamedTuple{(:weight,)});
    weights = [d.weight for d in draws], labels = [d.draw_id for d in draws])
_columns(rows) = NamedTuple{keys(first(rows))}(Tuple([getfield(r, k) for r in rows]
                                                     for k in keys(first(rows))))

# ---------------------------------------------------------------------------------------------
# Drawert swab / eDNA surveillance study (DecisionAnalysis.AbstractStudyDesign)

struct DrawertSurveyStudy{P<:DrawertParameters,O<:SurveyObservation,S} <: DecisionAnalysis.AbstractStudyDesign
    obs::O
    pond::SurveyPond
    template::P
    seeds::S
end

function Chytridiomycosis.drawert_survey_study(obs::SurveyObservation, pond::SurveyPond;
        nruns::Integer = 30, seeds = 1:nruns,
        template::Union{Nothing,DrawertParameters} = nothing)
    nruns >= 1 || throw(ArgumentError("nruns must be at least 1, got $nruns"))
    tmpl = template === nothing ? DrawertParameters(first(load_drawert_decision_draws())) : template
    return DrawertSurveyStudy(obs, pond, tmpl, seeds)
end

function _survey_params(s::DrawertSurveyStudy, θ)
    p = _with(s.template, θ)
    return p
end

function _survey_pond_for(pond::SurveyPond, N::Integer)
    pond.N == N && return pond
    max_swabs = N ÷ pond.pool
    return SurveyPond(; N, days = pond.days, swabs = min(pond.swabs, max_swabs),
                      pool = pond.pool, edna = pond.edna)
end

function DecisionAnalysis.simulate_study(rng::AbstractRNG, s::DrawertSurveyStudy, θ)
    p = _survey_params(s, θ)
    pond = _survey_pond_for(s.pond, p.N)
    prob = DrawertIBMProblem(p; tfinal = pond.days[end], saveat = pond.days, save_loads = true)
    run = solve(prob; rng)
    return simulate_survey(rng, p, run, s.obs, pond)
end

function DecisionAnalysis.study_loglikelihood(s::DrawertSurveyStudy,
        samples::AbstractVector{SurveySample}, θ)
    p = _survey_params(s, θ)
    return survey_loglikelihood(samples, p, s.obs; seeds = s.seeds)
end

function _drawert_survey_summary(samples::AbstractVector{SurveySample})
    swabs = Float64[x.copies for x in samples if x.kind === :swab]
    ednas = Float64[x.copies for x in samples if x.kind === :edna]
    swab_log = isempty(swabs) ? 0.0 : sum(log1p, swabs) / length(swabs)
    edna_log = isempty(ednas) ? 0.0 : sum(log1p, ednas) / length(ednas)
    swab_pos = isempty(swabs) ? 0.0 : count(>(0), swabs) / length(swabs)
    edna_pos = isempty(ednas) ? 0.0 : count(>(0), ednas) / length(ednas)
    return (; swab_log_mean = swab_log, edna_log_mean = edna_log,
            swab_pos_frac = swab_pos, edna_pos_frac = edna_pos)
end

function DecisionAnalysis.SummaryStatistics(
        belief::Union{DecisionAnalysis.SampleBelief,DecisionAnalysis.DecisionProblem},
        design::DrawertSurveyStudy;
        seed::Union{Nothing,Integer} = nothing,
        regression::DecisionAnalysis.AbstractRegression = DecisionAnalysis.SplineRegression())
    seed === nothing && throw(ArgumentError("give `seed`"))
    rng = Xoshiro(seed)
    b = belief isa DecisionAnalysis.SampleBelief ? belief : belief.belief
    n = length(b.weights)
    cols = b.draws
    ks = keys(cols)
    summaries = map(1:n) do i
        θ = NamedTuple{ks}(Tuple(getfield(cols, k)[i] for k in ks))
        _drawert_survey_summary(DecisionAnalysis.simulate_study(rng, design, θ))
    end
    return DecisionAnalysis.SummaryStatistics(summaries; regression)
end

function DecisionAnalysis.study_posterior(::AbstractRNG, ::DrawertSurveyStudy, data, prior;
        draws::Integer = 1000)
    throw(ArgumentError("DrawertSurveyStudy has no closed-form conjugate posterior; use " *
                        "SummaryStatistics or ImportanceSampling"))
end

# ---------------------------------------------------------------------------------------------
# Canessa et al. (2018): CanessaIPMProblem, CanessaSIProblem, canessa_outcome_table

function DecisionAnalysis.with_parameters(prob::CanessaIPMProblem, θ::CanessaDraw)
    return CanessaIPMProblem(CanessaParameters(θ), prob.I0, prob.nsteps, prob.mesh, prob.cap,
                             prob.units)
end

function DecisionAnalysis.with_parameters(prob::CanessaIPMProblem, θ)
    q, _ = _update_struct(CanessaParameters, prob.p, θ; aliases = (; beta0 = :β0))
    return CanessaIPMProblem(q, prob.I0, prob.nsteps, prob.mesh, prob.cap, prob.units)
end

function _canessa_action_obj(action)
    action === nothing && return canessa_action(:a)
    action isa CanessaAction && return action
    action isa Symbol && return canessa_action(action)
    action isa NamedTuple && return CanessaAction(:custom; action...)
    throw(ArgumentError("an action for CanessaIPMProblem is `nothing`, a Symbol (:a..:h), " *
                        "a CanessaAction, or a NamedTuple of CanessaAction keywords; " *
                        "got $(typeof(action))"))
end

function DecisionAnalysis.apply_action(prob::CanessaIPMProblem, action)
    a = _canessa_action_obj(action)
    a.growth_factor == 1 ||
        throw(ArgumentError("apply_action(::CanessaIPMProblem, action) requires growth_factor == 1 " *
                            "(CanessaParameters stores fitted b4..b6, not raw growth); use " *
                            "canessa_outcome_table for actions that alter within-host growth"))
    p = prob.p
    keep = 1 - a.removal
    q = CanessaParameters(p.b0, p.b1, p.s0 * keep, p.b2, p.b3, p.β0 * a.β_factor,
                          p.μG0, p.σ2G0, p.b4, p.b5, p.b6, p.Lmax, p.σ2G, p.b7,
                          p.S0 * a.density_factor, a.β_low_factor, a.low_load,
                          a.removal, p.step_days)
    return CanessaIPMProblem(q, prob.I0, prob.nsteps, prob.mesh, prob.cap, prob.units)
end

DecisionAnalysis.simulate(prob::CanessaIPMProblem; rng) = solve(prob)

function DecisionAnalysis.outcome(prob::CanessaIPMProblem, run)
    R0 = basic_reproduction_number(prob.p; mesh = prob.mesh, cap = prob.cap)
    return (; initial_susceptibles = run.S[1], final_susceptibles = run.S[end],
            final_infected = run.I[end], final_total = run.S[end] + run.I[end],
            R0, eradicated = Int(R0 < 1))
end

function DecisionAnalysis.with_parameters(prob::CanessaSIProblem, θ)
    q, _ = _update_struct(CanessaSIParameters, prob.p, θ; aliases = (; beta = :β, s_bar = :s̄))
    return CanessaSIProblem(q; mβ = prob.mβ, mS0 = prob.mS0, msurv = prob.msurv,
                            I0 = prob.I0, nsteps = prob.nsteps)
end

function DecisionAnalysis.apply_action(prob::CanessaSIProblem, action)
    if action === nothing
        return CanessaSIProblem(prob.p; mβ = 0.0, mS0 = 0.0, msurv = 0.0,
                                I0 = prob.I0, nsteps = prob.nsteps)
    elseif action isa NamedTuple
        allowed = (:mβ, :mS0, :msurv)
        for k in keys(action)
            k in allowed ||
                throw(ArgumentError("unknown field :$k for CanessaSIProblem action; expected $allowed"))
        end
        return CanessaSIProblem(prob.p; mβ = get(action, :mβ, 0.0),
                                mS0 = get(action, :mS0, 0.0),
                                msurv = get(action, :msurv, 0.0),
                                I0 = prob.I0, nsteps = prob.nsteps)
    end
    throw(ArgumentError("an action for CanessaSIProblem is `nothing` or a NamedTuple " *
                        "(; mβ, mS0, msurv); got $(typeof(action))"))
end

DecisionAnalysis.simulate(prob::CanessaSIProblem; rng) = solve(prob)

function DecisionAnalysis.outcome(prob::CanessaSIProblem, run)
    S_free = prob.p.S0 * prob.p.s0^prob.nsteps
    R0 = basic_reproduction_number(prob.p; mβ = prob.mβ, mS0 = prob.mS0, msurv = prob.msurv)
    return (; initial_susceptibles = run.S[1], final_susceptibles = run.S[end],
            final_infected = run.I[end],
            density_ratio = S_free > 0 ? run.S[end] / S_free : 0.0,
            R0, eradicated = Int(R0 < 1))
end

function _normalize_canessa_action_entry(x)
    x isa Pair && return Symbol(first(x)) => _canessa_action_obj(last(x))
    x isa CanessaAction && return x.name => x
    x isa Symbol && return x => canessa_action(x)
    throw(ArgumentError("each element of `actions` must be a Symbol, CanessaAction, or Pair; got $(typeof(x))"))
end

function Chytridiomycosis.canessa_outcome_table(; actions = CANESSA_ACTION_NAMES,
        draws::AbstractVector{CanessaDraw} = load_canessa2018(),
        units::Symbol = :transect, I0::Union{Nothing,Real} = nothing,
        nsteps::Integer = 30, mesh = Chytridiomycosis._CANESSA_DEFAULT_MESH,
        cap::Symbol = :truncate, v::Real = 0.5, s0::Symbol = :fit)
    pairs_ = [_normalize_canessa_action_entry(a) for a in actions]
    act_labels = Symbol[first(p) for p in pairs_]
    n_act = length(pairs_)
    n_draws = length(draws)
    out_names = Symbol[:final_susceptibles, :density_ratio, :R0, :eradicated,
                       :lifespan, :max_distance, :mean_distance]
    vals = Array{Float64,4}(undef, n_act, n_draws, 1, length(out_names))
    # Baseline disease-free final susceptibles per draw for exact adult count recovery
    S_free = Vector{Float64}(undef, n_draws)
    for i in 1:n_draws
        ref = CanessaParameters(draws[i]; s0)
        S_free[i] = Chytridiomycosis._canessa_units(ref, units, I0).S0 * ref.s0^nsteps
    end
    for (a, (_, act)) in enumerate(pairs_)
        ens = canessa_ensemble(act; draws, units, I0, nsteps, mesh, cap, v, s0)
        for i in 1:n_draws
            r = ens[i]
            vals[a, i, 1, 1] = r.ratio * S_free[i]
            vals[a, i, 1, 2] = r.ratio
            vals[a, i, 1, 3] = r.R0
            vals[a, i, 1, 4] = r.R0 < 1 ? 1.0 : 0.0
            vals[a, i, 1, 5] = r.lifespan
            vals[a, i, 1, 6] = r.max_distance
            vals[a, i, 1, 7] = r.mean_distance
        end
    end
    param_cols = (
        S_0 = [d.survival[1] for d in draws],
        S_10 = [d.survival[2] for d in draws],
        S_100 = [d.survival[3] for d in draws],
        S_1000 = [d.survival[4] for d in draws],
        S_10000 = [d.survival[5] for d in draws],
        T_0 = [d.transmission[1] for d in draws],
        T_10 = [d.transmission[2] for d in draws],
        T_100 = [d.transmission[3] for d in draws],
        T_1000 = [d.transmission[4] for d in draws],
        T_10000 = [d.transmission[5] for d in draws],
        Lmax = [d.Lmax for d in draws],
        growth = [d.growth for d in draws],
        beta0 = [d.β0 for d in draws],
        S0 = [d.S0 for d in draws],
    )
    belief = DecisionAnalysis.SampleBelief(param_cols; labels = [Symbol(:d, d.id) for d in draws])
    return DecisionAnalysis.OutcomeTable(belief.states, act_labels, out_names, [:r1], vals,
                                         fill("ok", n_act, n_draws, 1), belief.draws,
                                         belief.weights)
end

# ---------------------------------------------------------------------------------------------
# Cook et al. (2022): CookProjectionProblem

function DecisionAnalysis.with_parameters(prob::CookProjectionProblem, θ)
    nt = NamedTuple(θ)
    fields = fieldnames(CookParameters)
    used = filter(in(fields), keys(nt))
    isempty(used) &&
        throw(ArgumentError("the draw names none of the CookParameters fields $fields: $(keys(nt))"))
    p = prob.p
    val(k) = k in used ? getfield(nt, k) : getfield(p, k)
    q = CookParameters(; SU = val(:SU), ST = val(:ST),
                       adult_survival = collect(Float64, val(:adult_survival)),
                       fecundity = collect(Float64, val(:fecundity)),
                       yopy_ratio = Float64(val(:yopy_ratio)),
                       growth = Float64(val(:growth)),
                       female_fraction = Float64(val(:female_fraction)),
                       paired = Bool(val(:paired)))
    return CookProjectionProblem(q; treated_fraction = prob.treated_fraction,
                                 treated_years = prob.treated_years, years = prob.years,
                                 u0 = prob.u0, breeding = prob.breeding, schedule = prob.schedule)
end

function DecisionAnalysis.apply_action(prob::CookProjectionProblem, action)
    if action === nothing
        return CookProjectionProblem(prob.p; treated_fraction = 0.0, treated_years = 0.0,
                                     years = prob.years, u0 = prob.u0,
                                     breeding = prob.breeding, schedule = prob.schedule)
    elseif action isa NamedTuple
        allowed = (:treated_fraction, :treated_years, :schedule)
        for k in keys(action)
            k in allowed ||
                throw(ArgumentError("unknown field :$k for CookProjectionProblem action; expected $allowed"))
        end
        return CookProjectionProblem(prob.p;
                                     treated_fraction = get(action, :treated_fraction, 0.0),
                                     treated_years = get(action, :treated_years, 0.0),
                                     years = prob.years, u0 = prob.u0,
                                     breeding = prob.breeding,
                                     schedule = get(action, :schedule, prob.schedule))
    end
    throw(ArgumentError("an action for CookProjectionProblem is `nothing` or a NamedTuple " *
                        "(; treated_fraction, treated_years[, schedule]); got $(typeof(action))"))
end

DecisionAnalysis.simulate(prob::CookProjectionProblem; rng) = solve(prob; nsims = 1, rng)

function DecisionAnalysis.outcome(prob::CookProjectionProblem, run::CookProjectionSolution)
    T_end = size(run.states, 1)
    m_end = run.states[T_end, 1, 1]
    a_end = run.states[T_end, 4, 1]
    tot_end = sum(view(run.states, T_end, :, 1))
    λ = run.λ[1]
    return (; initial_adults = prob.u0[4], final_adults = a_end, final_metamorphs = m_end,
            final_total = tot_end, lambda = λ, positive_growth = Int(λ > 1),
            extinct = Int(a_end == 0))
end

# ---------------------------------------------------------------------------------------------
# Heard et al. (2024): HeardMetapopulationProblem

function DecisionAnalysis.with_parameters(prob::HeardMetapopulationProblem, θ)
    coupled = prob.p.subadult_survival == prob.p.adult_survival
    q, used = _update_struct(HeardParameters, prob.p, θ)
    if coupled && haskey(used, :adult_survival) && !haskey(used, :subadult_survival)
        q = HeardParameters(q.adult_survival, q.adult_survival, q.carrying_capacity,
                            q.prop_dispersing, q.mets, q.mets_cv, q.mortality_cv,
                            q.met_survival, q.sex_ratio, q.subadult_breeding_fraction,
                            q.initial_fraction, q.time_steps, q.cell_size,
                            q.distance_decay, q.max_distance)
    end
    return HeardMetapopulationProblem(q, prob.habitat; as_coded = prob.as_coded)
end

function DecisionAnalysis.apply_action(prob::HeardMetapopulationProblem, action)
    action === nothing && return prob
    action isa NamedTuple ||
        throw(ArgumentError("an action for HeardMetapopulationProblem is `nothing` or a NamedTuple; " *
                            "got $(typeof(action))"))
    allowed = (:adult_survival, :subadult_survival, :carrying_capacity, :prop_dispersing,
               :mets, :survival_multiplier, :capacity_multiplier, :habitat, :landscape)
    for k in keys(action)
        k in allowed ||
            throw(ArgumentError("unknown field :$k for HeardMetapopulationProblem action; expected $allowed"))
    end
    p = prob.p
    coupled = p.subadult_survival == p.adult_survival
    as_ = Float64(get(action, :adult_survival, p.adult_survival))
    ss_ = Float64(get(action, :subadult_survival, coupled ? as_ : p.subadult_survival))
    if haskey(action, :survival_multiplier)
        mult = Float64(action.survival_multiplier)
        as_ = clamp(as_ * mult, 0.0, 1.0)
        ss_ = clamp(ss_ * mult, 0.0, 1.0)
    end
    cap_ = Float64(get(action, :carrying_capacity, p.carrying_capacity))
    if haskey(action, :capacity_multiplier)
        cap_ *= Float64(action.capacity_multiplier)
    end
    disp_ = Float64(get(action, :prop_dispersing, p.prop_dispersing))
    mets_ = Float64(get(action, :mets, p.mets))
    q = HeardParameters(as_, ss_, cap_, disp_, mets_, p.mets_cv, p.mortality_cv,
                        p.met_survival, p.sex_ratio, p.subadult_breeding_fraction,
                        p.initial_fraction, p.time_steps, p.cell_size,
                        p.distance_decay, p.max_distance)
    hab = if haskey(action, :habitat)
        BitMatrix(action.habitat)
    elseif haskey(action, :landscape)
        heard2024_landscapes().habitat[Int(action.landscape)]
    else
        prob.habitat
    end
    return HeardMetapopulationProblem(q, hab; as_coded = prob.as_coded)
end

DecisionAnalysis.simulate(prob::HeardMetapopulationProblem; rng) = solve(prob; rng)

function DecisionAnalysis.outcome(::HeardMetapopulationProblem, run::HeardMetapopulationSolution)
    T_end = size(run.adults, 1)
    ad_end = sum(view(run.adults, T_end, :))
    sub_end = sum(view(run.subadults, T_end, :))
    occ = count(j -> run.adults[T_end, j] + run.subadults[T_end, j] > 0, axes(run.adults, 2))
    return (; initial_females = run.total[1], final_females = run.total[end],
            final_adults = ad_end, final_subadults = sub_end,
            occupied_patches = occ, persisted = Int(!run.extinct), extinct = Int(run.extinct))
end

# ---------------------------------------------------------------------------------------------
# Turner et al. (2025): TurnerPVAProblem

function DecisionAnalysis.with_parameters(prob::TurnerPVAProblem, θ)
    q, _ = _update_struct(TurnerParameters, prob.p, θ)
    return TurnerPVAProblem(q; survival = prob.survival, start = prob.start,
                            pattern = prob.pattern, wetness = prob.wetness,
                            replicates = prob.replicates, as_coded = prob.as_coded)
end

function DecisionAnalysis.apply_action(prob::TurnerPVAProblem, action)
    action === nothing && return prob
    if action isa Symbol
        action in TURNER_PATTERNS ||
            throw(ArgumentError("pattern must be one of $TURNER_PATTERNS; got :$action"))
        return TurnerPVAProblem(prob.p; survival = prob.survival, start = prob.start,
                                pattern = action, wetness = nothing,
                                replicates = prob.replicates, as_coded = prob.as_coded)
    elseif action isa NamedTuple
        allowed = (:month, :pattern, :p_wet)
        for k in keys(action)
            k in allowed ||
                throw(ArgumentError("unknown field :$k for TurnerPVAProblem action; expected $allowed"))
        end
        pat = get(action, :pattern, prob.pattern)
        pat in TURNER_PATTERNS ||
            throw(ArgumentError("pattern must be one of $TURNER_PATTERNS; got :$pat"))
        p = haskey(action, :p_wet) ? first(_update_struct(TurnerParameters, prob.p,
                                                          (; p_wet = action.p_wet))) : prob.p
        if haskey(action, :month)
            m = action.month
            m in TURNER_MONTHS ||
                throw(ArgumentError("month must be one of $TURNER_MONTHS; got :$m"))
            d = load_turner2025(; month = m, pattern = pat)
            return TurnerPVAProblem(p; survival = d.survival, start = d.start,
                                    pattern = pat, wetness = nothing,
                                    replicates = prob.replicates, as_coded = true)
        end
        return TurnerPVAProblem(p; survival = prob.survival, start = prob.start,
                                pattern = pat, wetness = nothing,
                                replicates = prob.replicates, as_coded = prob.as_coded)
    end
    throw(ArgumentError("an action for TurnerPVAProblem is `nothing`, a pattern Symbol " *
                        "($TURNER_PATTERNS), or a NamedTuple (; month, pattern, p_wet); " *
                        "got $(typeof(action))"))
end

function DecisionAnalysis.simulate(prob::TurnerPVAProblem; rng)
    w = prob.wetness === nothing ? nothing : prob.wetness[1:1, :]
    surv = prob.as_coded ? prob.survival : prob.survival[1:1, :, :]
    one = TurnerPVAProblem(prob.p; survival = surv, start = prob.start,
                           pattern = prob.pattern, wetness = w, replicates = 1,
                           as_coded = prob.as_coded)
    return solve(one; rng)
end

function DecisionAnalysis.outcome(::TurnerPVAProblem, run::TurnerPVASolution)
    pop = view(run.population, :, 1)
    min_f = minimum(pop)
    end_f = pop[end]
    return (; initial_females = pop[1], final_females = end_f, min_females = min_f,
            persisted = Int(end_f > 0), extinct = Int(end_f == 0),
            quasi_extinct_20 = Int(min_f < 20))
end

end
