# Brannelly LA, Chatfield MWH, Sonn J, Robak M, Richards-Zawacki CL (2018) Fungal infection
# has sublethal effects in a lowland subtropical amphibian population. BMC Ecology 18:34.
# https://doi.org/10.1186/s12898-018-0189-5 (open access, CC BY 4.0; data CC0)
# This implements a published model; it is not a new model.
#
# Read from the full text (Europe PMC, PMC6137908: Methods "Population modeling", Results,
# Tables 1 and 2, Fig. 3) and Additional files 1–3 (bundled in data/brannelly2018). No code
# was published: the authors used M-SURGE (multistate) and Program MARK (POPAN).
#
# Model. Northern cricket frogs (Acris crepitans), Belle Chasse, Louisiana, marked Nov 2010 –
# Nov 2011. Two alive states, 1 = Bd-negative, 2 = Bd-positive (swab qPCR > 1 ZE), and dead.
# Conditional Arnason–Schwarz model: survival Φ, capture p and state transition Ψ, each with
# effects of time in months (t), state at the previous capture (f, "from") and state at the
# current capture (to). Ten models are printed (Table 2); the three within 2 AIC of the best
# (Φ(.)p(to+t)Ψ(f*to+t), Φ(.)p(t)Ψ(f*to+t), Φ(f)p(to+t)Ψ(f*to+t)) were averaged. POPAN
# (Φ, p, pent, each (.) or (t)) ignores disease state and gives the super-population size.
#
# Mapping to the multistate layer (src/observation/multistate.jl, cmr_open.jl). With two
# states Ψ has one free cell per row (Bd− → Bd+, Bd+ → Bd−): "f*to" is one parameter per cell
# (:state), "f*to + t" adds a time effect shared by the two cells on the logit scale
# (:state_plus_time); p(to) is one parameter per current state (:state), "to + t" is
# :state_plus_time. These reproduce every parameter count of Table 2 with 10 occasions
# (e.g. Φ(.)p(to+t)Ψ(f*to+t): 1 + 10 + 10 = 21). "Ψ(i)" (12 parameters, so Ψ has one) is not
# defined in the paper; it is taken as :constant, the same as Ψ(.).
#
# Choices where the paper is silent (documented, not printed):
# - Occasions. "Surveys were conducted twice a month, each month's captures were combined":
#   the 19 survey days of Additional file 3 fall in 10 months. Two groupings are possible for
#   28 Sep and 1 Oct 2011: with the 30 Sep / 1 Oct pair as "September" (the two survey days of
#   a month as stated in Methods; `occasions = :survey_pairs`, the default) or by calendar
#   month (`occasions = :calendar_month`, which puts three days in October). Both give 10
#   occasions. The survey-pair grouping comes closer to the printed deviances (see the
#   validation targets).
# - Intervals. Jan, May and Jul 2011 were not surveyed, so the 9 intervals are 1, 2, 1, 1, 2,
#   2, 1, 1, 1 months. Survival is "monthly" in the paper; we raise it to the interval length
#   (`intervals = :months`, default), as MARK and M-SURGE do when intervals are entered.
#   `intervals = :equal` treats every interval as one step.
# - State within a month. A frog caught on both days of a month is Bd-positive if either
#   swab was positive (the paper averaged the two loads; the average of a positive and a zero
#   load is positive).
# - Additional file 3 codes infected captures as 2 or 3 (3, undocumented, is disease
#   category 2, log10 ZE > 2.455; see data/brannelly2018/README.md); both are Bd-positive.

"""
    brannelly2018_path() -> String

Directory of the bundled Brannelly et al. (2018) Additional files (CC0;
`data/brannelly2018/README.md`).
"""
brannelly2018_path() = normpath(joinpath(@__DIR__, "..", "..", "data", "brannelly2018"))

"""
    Brannelly2018Data

Capture–recapture data of Brannelly et al. (2018) from [`load_brannelly2018`](@ref):

- `ids`: animal IDs (one per row);
- `survey_dates`: the 19 survey days (`Date`s);
- `survey_codes`: `n×19` codes as in Additional file 3 (0 not caught, 1 Bd−, 2 Bd+, 3 Bd+
  with log10 ZE > 2.455);
- `occasion_of_survey`: the occasion (1…10) of each survey day;
- `histories`: `n×10` two-state histories (0 not caught, 1 Bd−, 2 Bd+);
- `intervals`: the 9 interval lengths (months, or ones with `intervals = :equal`);
- `occasions`, `interval_unit`: how occasions and intervals were formed.
"""
struct Brannelly2018Data
    ids::Vector{Int}
    survey_dates::Vector{Date}
    survey_codes::Matrix{Int}
    occasion_of_survey::Vector{Int}
    histories::Matrix{Int}
    intervals::Vector{Float64}
    occasions::Symbol
    interval_unit::Symbol
end

"""
    load_brannelly2018(; path = brannelly2018_path(), occasions = :survey_pairs,
                       intervals = :months) -> Brannelly2018Data

Read Additional file 3 (capture histories over 19 survey days) and collapse it to 10 monthly
occasions, as the paper did. `occasions = :survey_pairs` (default) takes the two survey days
of each month, with 28 Sep and 1 Oct 2011 as September; `:calendar_month` groups by calendar
month. `intervals = :months` (default) gives the interval lengths in months (1, 2, 1, 1, 2,
2, 1, 1, 1); `:equal` gives ones. See the source header for these choices.
"""
function load_brannelly2018(; path::AbstractString = brannelly2018_path(),
                            occasions::Symbol = :survey_pairs, intervals::Symbol = :months)
    occasions in (:survey_pairs, :calendar_month) ||
        throw(ArgumentError("occasions must be :survey_pairs or :calendar_month"))
    intervals in (:months, :equal) || throw(ArgumentError("intervals must be :months or :equal"))
    lines = readlines(joinpath(path, "additional_file3_capture_histories.csv"))
    header = split(lines[2], ',')
    dates = [Date(String(d), dateformat"m/d/y") + Dates.Year(2000) for d in header[8:26]]
    rows = [split(l, ',') for l in lines[3:end]]
    rows = filter(r -> !isempty(strip(r[1])), rows)
    ids = [parse(Int, r[1]) for r in rows]
    codes = [parse(Int, rows[i][7+j]) for i in eachindex(rows), j in 1:19]
    key = occasions === :survey_pairs ? [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 7, 7, 8, 8, 9, 9, 10, 10] :
          (ym = [(Dates.year(d), Dates.month(d)) for d in dates];
           u = unique(ym); [findfirst(==(x), u) for x in ym])
    T = maximum(key)
    H = zeros(Int, length(ids), T)
    for i in eachindex(ids), t in 1:T
        xs = codes[i, key .== t]
        H[i, t] = all(iszero, xs) ? 0 : any(>=(2), xs) ? 2 : 1
    end
    months = [Dates.year(dates[findfirst(==(t), key)]) * 12 +
              Dates.month(dates[findfirst(==(t), key)]) for t in 1:T]
    τ = intervals === :months ? Float64.(diff(months)) : ones(T - 1)
    return Brannelly2018Data(ids, dates, codes, key, H, τ, occasions, intervals)
end

"""
    brannelly2018_summary(data = load_brannelly2018()) -> NamedTuple

Counts of the kind of Table 1 of Brannelly et al. (2018) from the bundled histories:
`individuals` (animals marked), `capture_events` (survey-day captures), `negative`,
`positive` (survey-day captures by swab result), `recaptured` (animals caught on more than
one survey day) and `monthly_events` (captures after collapsing to occasions).
"""
function brannelly2018_summary(d::Brannelly2018Data = load_brannelly2018())
    c = d.survey_codes
    return (individuals = length(d.ids), capture_events = count(>(0), c),
            negative = count(==(1), c), positive = count(>=(2), c),
            recaptured = count(i -> count(>(0), view(c, i, :)) > 1, axes(c, 1)),
            monthly_events = count(>(0), d.histories))
end

"""
    brannelly2018_structure(; survival = :constant, detection = :state_plus_time,
                            transition = :state_plus_time) -> MultiEventStructure

The two-state (1 = *Bd*−, 2 = *Bd*+) multistate structure of Brannelly et al. (2018) for an
[`ArnasonSchwarzProblem`](@ref); states are observed without error. Designs are those of
[`MultiEventStructure`](@ref). The defaults are the best model of Table 2,
Φ(.)p(to + t)Ψ(f * to + t).
"""
brannelly2018_structure(; survival::Symbol = :constant, detection::Symbol = :state_plus_time,
                        transition::Symbol = :state_plus_time) =
    MultiEventStructure(; nstates = 2, assignment = Bool[1 0; 0 1], baseline = [1, 2], survival,
                        detection, transition, state_assignment = :state)

"""
    brannelly2018_models() -> Vector{NamedTuple}

The ten conditional Arnason–Schwarz models of Brannelly et al. (2018, Table 2) in printed
order: `name` (as printed), `structure` (a [`MultiEventStructure`](@ref) from
[`brannelly2018_structure`](@ref)) and `parameters` (the printed number of parameters).
"""
function brannelly2018_models()
    rows = [("Φ(.)p(to+t)Ψ(f*to+t)", :constant, :state_plus_time, :state_plus_time, 21),
            ("Φ(.)p(t)Ψ(f*to+t)", :constant, :time, :state_plus_time, 20),
            ("Φ(f)p(to+t)Ψ(f*to+t)", :state, :state_plus_time, :state_plus_time, 22),
            ("Φ(f)p(t)Ψ(f*to+t)", :state, :time, :state_plus_time, 21),
            ("Φ(t)p(t)Ψ(f*to+t)", :time, :time, :state_plus_time, 28),
            ("Φ(.)p(to+t)Ψ(t)", :constant, :state_plus_time, :time, 20),
            ("Φ(t)p(to+t)Ψ(f*to+t)", :time, :state_plus_time, :state_plus_time, 29),
            ("Φ(.)p(to+t)Ψ(i)", :constant, :state_plus_time, :constant, 12),
            ("Φ(f*t)p(to*t)Ψ(f*to*t)", :state_time, :state_time, :state_time, 54),
            ("Φ(.)p(.)Ψ(.)", :constant, :constant, :constant, 3)]
    return [(name = r[1], structure = brannelly2018_structure(; survival = r[2], detection = r[3],
                                                              transition = r[4]),
             parameters = r[5]) for r in rows]
end

"""
    brannelly2018_popan_models() -> Vector{NamedTuple}

The six POPAN models of Brannelly et al. (2018, Table 2) in printed order: `name` and the
`survival`, `detection` and `entry` designs (`:constant` or `:time`) for a
[`PopanProblem`](@ref).
"""
brannelly2018_popan_models() =
    [(name = "Φ(.)p(t)pent(t)", survival = :constant, detection = :time, entry = :time),
     (name = "Φ(t)p(t)pent(t)", survival = :time, detection = :time, entry = :time),
     (name = "Φ(t)p(.)pent(t)", survival = :time, detection = :constant, entry = :time),
     (name = "Φ(.)p(.)pent(t)", survival = :constant, detection = :constant, entry = :time),
     (name = "Φ(.)p(.)pent(.)", survival = :constant, detection = :constant, entry = :constant),
     (name = "Φ(t)p(t)pent(.)", survival = :time, detection = :time, entry = :constant)]

"""
    brannelly2018_parameters() -> NamedTuple

The published estimates of Brannelly et al. (2018) (compendium `brannelly2018`):
model-averaged monthly survival of *Bd*-negative and *Bd*-positive frogs, the POPAN
super-population size and the annual survival the paper derives from monthly survival.
"""
function brannelly2018_parameters()
    v = parameter_set(:brannelly2018)
    return (survival_uninfected = v.survival_uninfected, survival_infected = v.survival_infected,
            superpopulation = v.superpopulation, annual_survival = v.annual_survival)
end
