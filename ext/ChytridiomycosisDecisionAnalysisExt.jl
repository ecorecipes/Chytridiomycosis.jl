module ChytridiomycosisDecisionAnalysisExt

# The adapter between this package and DecisionAnalysis.jl: methods of its adapter verbs for the
# two problems of the Drawert et al. (2017) model, so that
#
#     outcomes(prob, actions, draws; replicates, seed)
#
# runs the model over parameter draws × actions × replicates with common random numbers, and the
# stored table of runs as an `OutcomeTable`. The problems, the interventions and the solver are
# this package's own and need no decision package; the methods here only name them in
# DecisionAnalysis's vocabulary. Outcomes are counts of frogs: utilities are left to the user.

using Chytridiomycosis: Chytridiomycosis, DrawertIBMProblem, DrawertIntervention, DrawertODEProblem,
                        DrawertParameters, load_drawert_decision_draws, load_drawert_decisions,
                        solve
import DecisionAnalysis

# ---------------------------------------------------------------------------------------------
# Parameters

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

# ---------------------------------------------------------------------------------------------
# Actions

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

# ---------------------------------------------------------------------------------------------
# Runs and their outcomes

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

# ---------------------------------------------------------------------------------------------
# The stored table

const _STORED_OUTCOMES = (:survivors, :uninfected, :infected, :uninfected_d120, :infected_d120,
                          :uninfected_d365, :infected_d365, :initial_frogs, :handled,
                          :treatment_deaths)

function Chytridiomycosis.drawert_outcome_table(; variant::AbstractString = "baseline",
        outcomes = _STORED_OUTCOMES)
    runs = filter(r -> r.variant == variant, load_drawert_decisions())
    isempty(runs) && throw(ArgumentError("the stored table has no runs of the variant $variant"))
    return DecisionAnalysis.OutcomeTable(runs; draw = :draw_id, action = :action,
                                         replicate = :replicate, outcomes = collect(Symbol, outcomes),
                                         status = :status, parameters = load_drawert_decision_draws(),
                                         key = :draw_id, weights = :weight)
end

end
