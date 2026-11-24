# Table of Drawert et al. (2017) model outcomes by draw × action × replicate, for decision
# analysis (written by scripts/drawert2017_decision_table.jl). These are OUR runs, not published
# results; the model and the actions are the paper's. Utilities, costs and decision rules are
# left to the user.

# SplitMix64 finaliser: a fixed integer hash, the same on every Julia version (Base.hash is not).
function _splitmix64(x::UInt64)
    x += 0x9e3779b97f4a7c15
    x = (x ⊻ (x >> 30)) * 0xbf58476d1ce4e5b9
    x = (x ⊻ (x >> 27)) * 0x94d049bb133111eb
    return x ⊻ (x >> 31)
end

"""
    drawert_run_seed(draw_id, replicate) -> UInt64

The seed of replicate `replicate` at point `draw_id` in the decision tables, for
`Xoshiro(seed)`. It does not depend on the action or the model variant, so a run and its
counterparts under other actions use **common random numbers**: they are identical up to the
first intervention. After that they diverge, because the model draws from one stream in event
order and has no per-frog streams. The seed is a fixed integer hash, the same on every Julia
version; the random numbers made from it can differ between Julia versions.
"""
drawert_run_seed(draw_id::Integer, replicate::Integer) =
    _splitmix64(_splitmix64(UInt64(draw_id)) + UInt64(replicate))

# Stable label of an action.
_drawert_action_label(::Nothing) = "none"
function _drawert_action_label(iv::DrawertIntervention)
    num(x) = isinteger(x) ? string(Int(x)) : string(x)
    s = "$(iv.kind)_d$(num(iv.day))_f$(num(round(100 * iv.fraction; digits = 6)))"
    return iv.mortality > 0 ? "$(s)_m$(num(round(100 * iv.mortality; digits = 6)))" : s
end

# Sample sizes per group: proportional to the square root of the group's size, at least one
# each, adding up to n (largest remainders).
function _drawert_allocation(sizes::AbstractVector{<:Integer}, n::Integer)
    n >= length(sizes) || throw(ArgumentError("need at least one draw per group"))
    n >= sum(sizes) && return collect(Int, sizes)
    x = sqrt.(sizes)
    quota = n .* x ./ sum(x)
    out = min.(max.(floor.(Int, quota), 1), sizes)
    while sum(out) != n
        d = sum(out) < n ? 1 : -1
        ok = [d == 1 ? out[i] < sizes[i] : out[i] > 1 for i in eachindex(out)]
        score = [ok[i] ? d * (quota[i] - out[i]) : -Inf for i in eachindex(out)]
        out[argmax(score)] += d
    end
    return out
end

# n of the 2132 draws, stratified by survival group, with weight = group size / group sample.
# Within a group the draws are those with the smallest seeded integer hashes.
function _drawert_draw_sample(n::Integer; seed::Integer = 20170480)
    draws = drawert2017_draws()
    groups = (:A, :B, :C, :D, :E)
    members = [filter(d -> d.group === g, draws) for g in groups]
    alloc = _drawert_allocation(length.(members), n)
    out = [(; d..., weight = length(m) / k) for (m, k) in zip(members, alloc)
           for d in sort(m; by = d -> _splitmix64(_splitmix64(UInt64(seed)) + UInt64(d.draw_id)))[1:k]]
    return sort!(out; by = d -> d.draw_id)
end

# One row of the table (all columns but the wall time).
const _DRAWERT_DECISION_COLUMNS = (:draw_id, :group, :variant, :action, :kind, :day, :fraction,
    :replicate, :seed, :surviving_fraction, :survivors, :outcome, :first_loss,
    :frogs_at_intervention, :applied, :handled, :treatment_deaths, :initial_frogs, :uninfected,
    :infected, :uninfected_d120, :infected_d120, :uninfected_d365, :infected_d365, :status)

# Fixed evaluation days of the table: within the season, and a year after Bd arrives.
const _DRAWERT_EVALUATION_DAYS = (120.0, 365.0)

# Run one cell: returns the row and the wall time. An error becomes a row with its status.
function _drawert_decision_run(draw, variant::Symbol, iv, replicate::Integer; S_star::Real,
        infection_threshold::Integer = 0)
    seed = drawert_run_seed(draw.draw_id, replicate)
    head = (; draw_id = draw.draw_id, group = String(draw.group), variant = String(variant),
            action = _drawert_action_label(iv), kind = iv === nothing ? "" : String(iv.kind),
            day = iv === nothing ? missing : iv.day,
            fraction = iv === nothing ? missing : iv.fraction, replicate = Int(replicate), seed)
    t0 = time_ns()
    row = try
        p = drawert_variant(DrawertParameters(draw), variant; S_star)
        r = solve(DrawertIBMProblem(p; intervention = iv, infection_threshold,
                                    evaluate = _DRAWERT_EVALUATION_DAYS); rng = Random.Xoshiro(seed))
        (; head..., surviving_fraction = r.surviving_fraction, survivors = r.survivors,
         outcome = String(r.outcome), first_loss = r.first_loss,
         frogs_at_intervention = iv === nothing ? missing : r.pulses[1].frogs,
         applied = iv === nothing ? missing : r.pulses[1].applied, handled = r.handled,
         treatment_deaths = r.treatment_deaths, initial_frogs = p.N, uninfected = r.uninfected,
         infected = r.infected, uninfected_d120 = r.evaluation[1].uninfected,
         infected_d120 = r.evaluation[1].infected,
         uninfected_d365 = r.evaluation[2].uninfected,
         infected_d365 = r.evaluation[2].infected, status = "ok")
    catch err
        err isa InterruptException && rethrow()
        (; head..., surviving_fraction = NaN, survivors = -1, outcome = "", first_loss = NaN,
         frogs_at_intervention = missing, applied = missing, handled = 0, treatment_deaths = 0,
         initial_frogs = draw.N, uninfected = -1, infected = -1, uninfected_d120 = -1,
         infected_d120 = -1, uninfected_d365 = -1, infected_d365 = -1, status = string("error:", nameof(typeof(err))))
    end
    return row, (time_ns() - t0) / 1e9
end

_drawert_cell(x) = x === missing ? "" : string(x)
_drawert_csv_line(row) = join((_drawert_cell(getfield(row, c)) for c in _DRAWERT_DECISION_COLUMNS), ',')

"""
    load_drawert_decisions(; path = data/drawert2017/decision_runs.csv) -> Vector{NamedTuple}

OUR runs of the Drawert et al. (2017) individual-based model over draws × actions ×
replicates (not published results), one named tuple per run, in a Tables.jl row-table form.
Written by `scripts/drawert2017_decision_table.jl`, whose header gives the design. Columns:
`draw_id`, `group`, `variant` (`baseline`, `nu` or `sigma`), `action` (for example `none`,
`clean_frogs_d28_f75`), `kind`, `day`, `fraction`, `replicate`, `seed`
([`drawert_run_seed`](@ref)), `surviving_fraction`, `survivors`, `outcome`
(`frogs_extinct`, `fungus_extinct` or `persistence`), `first_loss`, `frogs_at_intervention`
(alive just before the pulse), `applied`, `handled`, `treatment_deaths`, `initial_frogs`, `uninfected` and `infected` at the
end of the run (live frogs without and with sporangia: `survivors = uninfected + infected`; the
load threshold is 0), the same on days 120 and 365 (`uninfected_d120`, …), `status` (`ok`, or an
error label) and `seconds` (wall time; read from the separate `…_seconds.csv` so that the main
file is reproducible, `missing` if that file is absent).

`kind` is empty, and `day`, `fraction`, `frogs_at_intervention` and `applied` are `missing`, for
`none`. The rule for `handled` is under [`DrawertIBMProblem`](@ref): frogs cleaned or culled,
`Binomial(frogs_at_intervention, fraction)`; 0 when the pulse was not applied (the outcome was
already settled) and for `clean_environment`.

The parameters of the draws, with their sampling `weight` (points of the group per sampled
point), are in `decision_draws.csv`. The draws are a stratified sample of a space-filling
design, not a posterior.
"""
function load_drawert_decisions(; path::AbstractString = joinpath(dirname(drawert2017_path()),
                                                                  "decision_runs.csv"))
    isfile(path) || throw(ArgumentError("no decision table at $path; run " *
                                        "scripts/drawert2017_decision_table.jl"))
    lines = readlines(path)
    Symbol.(split(lines[1], ',')) == collect(_DRAWERT_DECISION_COLUMNS) ||
        throw(ArgumentError("unexpected columns in $path"))
    tpath = replace(path, r"\.csv$" => "_seconds.csv")
    secs = isfile(tpath) ? [parse(Float64, l) for l in readlines(tpath)[2:end]] : nothing
    n = length(lines) - 1
    secs === nothing || length(secs) == n ||
        throw(ArgumentError("$tpath does not have one row per run"))
    optf(x) = isempty(x) ? missing : parse(Float64, x)
    R = @NamedTuple{draw_id::Int, group::String, variant::String, action::String, kind::String,
        day::Union{Missing,Float64}, fraction::Union{Missing,Float64}, replicate::Int,
        seed::UInt64, surviving_fraction::Float64, survivors::Int, outcome::String,
        first_loss::Float64, frogs_at_intervention::Union{Missing,Int},
        applied::Union{Missing,Bool}, handled::Int, treatment_deaths::Int,
        initial_frogs::Int, uninfected::Int, infected::Int, uninfected_d120::Int,
        infected_d120::Int, uninfected_d365::Int, infected_d365::Int, status::String,
        seconds::Union{Missing,Float64}}
    return map(1:n) do i
        v = split(lines[i + 1], ',')
        R((parse(Int, v[1]), String(v[2]), String(v[3]), String(v[4]), String(v[5]), optf(v[6]),
           optf(v[7]), parse(Int, v[8]), parse(UInt64, v[9]), parse(Float64, v[10]),
           parse(Int, v[11]), String(v[12]), parse(Float64, v[13]),
           isempty(v[14]) ? missing : parse(Int, v[14]),
           isempty(v[15]) ? missing : parse(Bool, v[15]), (parse(Int, v[j]) for j in 16:24)...,
           String(v[25]), secs === nothing ? missing : secs[i]))
    end
end

"""
    load_drawert_decision_draws(; path = data/drawert2017/decision_draws.csv) -> Vector{NamedTuple}

The parameter draws of the stored decision table ([`load_drawert_decisions`](@ref)), one named
tuple per draw, in a Tables.jl row-table form: `draw_id`, `group` (the authors' survival group,
a `String`), the parameters `N`, `γ`, `η`, `ν0`, `f`, `σ0`, `μ`, the authors'
`load_growth_rate`, and `weight`, the number of points of the group that the draw stands for.

The draws are a stratified sample of the authors' space-filling design, **not a posterior**:
with these weights an average is an average over the 2132 design points, each counted once.
"""
function load_drawert_decision_draws(; path::AbstractString = joinpath(dirname(drawert2017_path()),
                                                                       "decision_draws.csv"))
    isfile(path) || throw(ArgumentError("no table of draws at $path; run " *
                                        "scripts/drawert2017_decision_table.jl"))
    lines = filter(!isempty ∘ strip, readlines(path))
    names = Tuple(Symbol.(split(lines[1], ',')))
    names == (:draw_id, :group, :N, :γ, :η, :ν0, :f, :σ0, :μ, :load_growth_rate, :weight) ||
        throw(ArgumentError("unexpected columns in $path"))
    return map(lines[2:end]) do line
        v = split(line, ',')
        NamedTuple{names}((parse(Int, v[1]), String(v[2]), parse(Int, v[3]),
                           (parse(Float64, x) for x in v[4:end])...))
    end
end

"""
    drawert_outcome_table(; variant = "baseline", outcomes = (:survivors, :uninfected, …))

The stored decision table ([`load_drawert_decisions`](@ref)) as a `DecisionAnalysis.OutcomeTable`
over draws × actions × replicates, with the parameters and design weights of
[`load_drawert_decision_draws`](@ref). `outcomes` names the columns to keep; the default is
`survivors`, `uninfected`, `infected`, the same on days 120 and 365, `initial_frogs`, `handled`
and `treatment_deaths`. Runs whose `status` is not `ok` are kept and counted as failures.

The method is defined when DecisionAnalysis.jl is loaded (a package extension); without it the
function has no method. The weights are those of a design, not of a posterior.
"""
function drawert_outcome_table end

# --------------------------------------------------------------------------------------------
# Survey data from an individual-based run (observation model: src/novel/swab_edna.jl)

"""
    simulate_survey(rng, p::DrawertParameters, run, obs::SurveyObservation, pond::SurveyPond)
        -> Vector{SurveySample}

Survey data from one run of the stochastic Drawert et al. (2017) model: `run` is the result of
`solve(DrawertIBMProblem(p; saveat, save_loads = true))` with every day of `pond.days` among
`saveat`. **New work** (model N5 applied to this simulator), not a published analysis.

The observation distributions are those of the moment method above, with the run's own state
in place of the linear-phase moments, so they hold from the Drawert model's start (one infected
frog, empty pool) and after frogs die. On each day, with `m` and `v` the mean and (population)
variance of load among live frogs and `w = v − m`:

- a tube of `n` pooled swabs has negative binomial copies, shape `n m²/w`, mean `n ρ m`
  (Poisson if `w ≤ 0`; 0 if no live frog is infected);
- an eDNA sample is lognormal around `κ Z/V` with log-sd `σe`, recorded as 0 below one copy
  (always 0 if `Z = 0`).

With fewer than `swabs × pool` live frogs, as many full tubes are taken as they allow. Use
`saveat` days before any intervention to describe data available for choosing it: `saveat`
does not change the run. `obs` has no default: **the assay scales `ρ`, `κ` and `σe` have no
published values**, so they must be given, and stated as assumptions.
"""
function simulate_survey(rng::AbstractRNG, p::DrawertParameters, run::NamedTuple,
        obs::SurveyObservation, pond::SurveyPond)
    pond.N == p.N || throw(ArgumentError("pond.N must equal p.N"))
    out = SurveySample[]
    for day in pond.days
        m, w, Z, n = _survey_state(run, day)
        for _ in 1:min(pond.swabs, n ÷ pond.pool)
            y = m == 0 ? 0.0 : Float64(rand(rng, _survey_tube(m, w, obs.ρ, pond.pool)))
            push!(out, SurveySample(1, day, :swab, pond.pool, y))
        end
        for _ in 1:pond.edna
            y = Z == 0 ? 0.0 : exp(log(obs.κ * Z / p.V) + obs.σe * randn(rng))
            push!(out, SurveySample(1, day, :edna, 1, y >= 1 ? y : 0.0))
        end
    end
    return out
end

# Mean and excess variance of load among live frogs, pool zoospores and live frogs on a day.
function _survey_state(run::NamedTuple, day::Real)
    hasproperty(run, :loads) ||
        throw(ArgumentError("the run must be solved with save_loads = true"))
    k = findfirst(==(day), run.t)
    k === nothing && throw(ArgumentError("day $day is not among the run's saveat days"))
    n = 0; s = 0.0; s2 = 0.0
    for x in @view run.loads[:, k]
        x >= 0 || continue
        n += 1; s += x; s2 += float(x)^2
    end
    n == 0 && return 0.0, 0.0, run.zoospores[k], 0
    m = s / n
    return m, s2 / n - m^2 - m, run.zoospores[k], n
end

_survey_tube(m, w, ρ, pool) =
    w > 1e-9 * m ? survey_swab_distribution(m, w, ρ, pool) : Poisson(pool * ρ * m)

"""
    survey_loglikelihood(samples, p::DrawertParameters, run, obs::SurveyObservation) -> Float64
    survey_loglikelihood(samples, p::DrawertParameters, obs::SurveyObservation; nruns = 30,
        seeds = 1:nruns) -> Float64

Log-likelihood of survey `samples` from one pond ([`simulate_survey`](@ref)) under the
stochastic Drawert et al. (2017) model. **New work**, not a published analysis.

The first form is conditional on one `run` (solved with `save_loads = true` at the sampled
days). The second integrates over runs by Monte Carlo: it solves `nruns` runs of `p`, with no
intervention and `Xoshiro(seed)` for each of `seeds`, and returns the log of the mean
likelihood. Use it as the log weight of `p` in importance weighting over draws; the estimate of
the likelihood is unbiased, so its noise shrinks as `nruns` grows. Pass the seeds of other runs
([`drawert_run_seed`](@ref)) to condition on the same trajectories. Samples are taken as
independent given the run. The number of tubes is not itself treated as data. `-Inf` means that
no run could have produced the samples: for example, copies where no frog was infected.
"""
function survey_loglikelihood(samples::AbstractVector{SurveySample}, p::DrawertParameters,
        run::NamedTuple, obs::SurveyObservation)
    ll = 0.0
    for s in samples
        m, w, Z, n = _survey_state(run, s.day)
        if s.kind === :swab
            n >= s.pool || return -Inf
            if m == 0
                s.copies == 0 || return -Inf
            else
                ll += logpdf(_survey_tube(m, w, obs.ρ, s.pool), round(Int, s.copies))
            end
        elseif Z == 0
            s.copies == 0 || return -Inf
        else
            ll += knapp_obs(s.copies, log(obs.κ * Z / p.V), obs.σe)
        end
    end
    return ll
end

function survey_loglikelihood(samples::AbstractVector{SurveySample}, p::DrawertParameters,
        obs::SurveyObservation; nruns::Integer = 30, seeds = 1:nruns)
    isempty(samples) && return 0.0
    days = sort!(unique(s.day for s in samples))
    prob = DrawertIBMProblem(p; tfinal = days[end], saveat = days, save_loads = true)
    lls = [survey_loglikelihood(samples, p, solve(prob; rng = Random.Xoshiro(k)), obs)
           for k in seeds]
    top = maximum(lls)
    isfinite(top) || return -Inf
    return top + log(sum(exp(l - top) for l in lls) / length(lls))
end
