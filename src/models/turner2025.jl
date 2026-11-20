# Turner, Heard, Mathwin, Bradshaw, Hall & Wassens (2025) Conserv Sci Pract 7:e70022,
# doi:10.1111/csp2.70022. Code and data: Zenodo doi:10.5281/zenodo.10921748 (CC BY 4.0).
#
# PUBLISHED MODEL, reimplemented from the AUTHORS' CODE (customFunctions/PVA_functionsV2.R,
# pop.func, and Run_scenariosV2.R) and checked against the paper (Methods 2.2, read
# 2026-11-20), which describes the model in outline. The paper's Results (Fig. 2, and the
# abstract) give extinction probabilities 0.13, 0.05, 0.04 (August, September, October; dry
# years independent), 0.025, 0.0004, <0.0001 (at most two dry years) and <0.001 (at most
# one); the archived code and outputs give higher values for August and September and for
# the capped scenarios, so the paper's figures apparently came from an earlier version.
#
# The Bd step (paper section 2.2.1) is NOT implemented: the paper describes it (monthly
# infection probability from a logistic regression on mean midnight water temperature,
# salinity 223 µS/cm and pH 7.9; zoospore load from a linear regression on air
# temperature; daily survival from the Heard et al. 2014 logistic load-survival relation
# with its intercept replaced by a base daily survival ~ (0.996, SD 0.001); 500
# simulated frogs per year; winter survival 0.692) but does not print the regression
# coefficients (they are in Turner, Wassens & Heard 2021 Biol Conserv 264:109391 and Heard
# et al. 2014 Ecol Appl 24:650) or the monthly temperature series, and the archive does not
# contain that code.
#
# Bd's role: the authors turned watering month (inflow temperature ~10, 15, 20 °C) into
# annual adult survival and summer survival of metamorphs through Bd load-survival
# relationships, and stored the result as precomputed arrays (referenceData/
# *WateringArray_*.csv: per replicate and year, wetness, adult survival mean and SD, summer
# survival mean and SD). The code that made the arrays is NOT in the archive, so the Bd
# link itself is not implemented here: the model takes the survival series as input.
#
# Model (female only; age classes 0-1 (row 1 holds EGGS), 1-2, 2-3, 3-4, 4-5, 5+; annual
# step from January; 52 years). For each replicate and year t:
#
#   hatch ~ Normal(0.92575, 0.0645), capped at 1
#   s_meta = Beta(3.645, 10.373) × hatch          egg to metamorph
#   s_0 = s_meta × Normal(m_summer(t), sd_summer(t)) × 0.996^92,  floored at 0
#   s_1..s_4 ~ Normal(m_adult(t), sd_adult(t)), independently, floored at 0
#   F_a = round(0.5 × max(Normal(3245, 897), 0)), independently for ages 1..5
#
#   M = [0 F_1 F_2 F_3 F_4 F_5; s_0 0 ...; 0 s_1 0 ...; ...; 0 0 0 0 s_4 0]   (5+ dies)
#
# Density dependence on the projection matrix, from the numbers at the start of the year:
#   M[2,1] ×= (a_t + b_t x)^(−1/c_t),               x = min(eggs/1417, 1.4)
#   M[3,2] ×= a_j/(exp((b_j − x)/c_j) + 1),          x = min(n_1/500, 0.65)
#   M[k+1,k] ×= the same curve, k = 3..5,            x = min(Σ n_1..n_5/500, 0.65)
# In a dry year (no environmental water) row 1 of M is 0. Then n(t+1) = M n(t); ages 1-5
# are rounded (half to even, as R). Eggs get a ceiling K = 129 × 3244.857 × 0.5: with
# x = eggs/K and f the decreasing logistic egg curve, eggs ×= (1/x)∫₀ˣ f if x ≤ 0.99, else
# eggs = K × mean_{[0,1.02]} f + f(1.02) (eggs − K); then rounded. The integral of the
# logistic is evaluated in closed form (R uses integrate()).
#
# A replicate is extinct when no female aged ≥ 1 is left in year 53 (the authors' measure,
# Run_scenariosV2.R). Start populations are the authors' burnt-in age structures
# (referenceData/Terminal_structure_*.csv), resampled with replacement from the extant ones.
#
# As coded versus as intended (`as_coded = true`, the default):
# - pop.func indexes the survival array with `iter` (= 10,000, the number of replicates)
#   instead of the replicate counter `rpt`, so EVERY replicate uses the year-by-year
#   survival means and SDs of the array's LAST replicate; only wetness uses `rpt`. The
#   bundled data/turner2025/survival_series.csv holds exactly those rows. With
#   `as_coded = false` each replicate uses its own row, which needs the full arrays
#   (load_turner2025_survival_array, from the archive).
# - The 5th adult survival draw is generated but unused (popmat[6,6] = 0), as coded.
# - The 10002-column Terminal_structure files carry a row-number column (1..6) that the R
#   code reads as a start population; we drop it (it is drawn with probability ~1e-4).
#
# Unpublished choices: wetness sequences are not generated in the archived code. When the
# full arrays are not used, turner_watering_sequence draws them by the rule the archived
# arrays follow: year 1 wet, later years wet with probability 0.7, and for `:max2dry` /
# `:max1dry` a year is forced wet after 2 / 1 consecutive dry years.

const TURNER_MONTHS = (:august, :september, :october)
const TURNER_PATTERNS = (:random, :max2dry, :max1dry)

"""
    TurnerParameters(:turner2025; as_printed = false)
    TurnerParameters(; kwargs...)

Life history and density dependence of the *Litoria raniformis* population viability model
of Turner et al. (2025), read from the compendium (the authors' code; the values the paper
prints were checked against it): `years`, `sex_ratio`, fecundity (`fecundity_mean`, `fecundity_sd`, eggs per
female aged 1–5), `hatch_mean`/`hatch_sd`, the beta distribution of survival from hatching
to metamorphosis (`tadpole_alpha`, `tadpole_beta`), winter survival of first-years
(`winter_daily_survival`^`winter_duration`), the egg ceiling (`spawn_capacity` ×
`clutch_size_mean` × `sex_ratio`), the density references `tadpole_capacity`,
`juvenile_capacity`, `adult_capacity` with their caps, the three fitted multiplier curves
(`egg_a/b/c`, `tadpole_a/b/c`, `juvenile_a/b/c`) and `p_wet`, the probability that a year is
watered. Bd-dependent survival is not a parameter: see [`TurnerPVAProblem`](@ref).
`as_printed = true` uses the paper's mean clutch size (3244) instead of the coded 3245.
"""
Base.@kwdef struct TurnerParameters
    years::Int
    sex_ratio::Float64
    fecundity_mean::Float64
    fecundity_sd::Float64
    clutch_size_mean::Float64
    hatch_mean::Float64
    hatch_sd::Float64
    tadpole_alpha::Float64
    tadpole_beta::Float64
    winter_daily_survival::Float64
    winter_duration::Float64
    spawn_capacity::Float64
    tadpole_capacity::Float64
    juvenile_capacity::Float64
    adult_capacity::Float64
    tadpole_ratio_cap::Float64
    juvenile_ratio_cap::Float64
    egg_linear_limit::Float64
    egg_overshoot_point::Float64
    egg_a::Float64
    egg_b::Float64
    egg_c::Float64
    tadpole_a::Float64
    tadpole_b::Float64
    tadpole_c::Float64
    juvenile_a::Float64
    juvenile_b::Float64
    juvenile_c::Float64
    p_wet::Float64
end

function TurnerParameters(key::Symbol; as_printed::Bool = false)
    v = parameter_set(key; set = as_printed ? "as_printed" : "")
    kw = Dict{Symbol,Any}()
    for f in fieldnames(TurnerParameters)
        kw[f] = f === :years ? round(Int, v.years) : getfield(v, f)
    end
    return TurnerParameters(; kw...)
end

time_unit(::Type{<:TurnerParameters}) = :year

with(p::TurnerParameters; kw...) =
    TurnerParameters(; (f => get(kw, f, getfield(p, f)) for f in fieldnames(TurnerParameters))...)

# logistic multiplier a/(exp((b − x)/c) + 1) and its integral from 0 to x
_turner_logistic(a, b, c, x) = a / (exp((b - x) / c) + 1)
_turner_softplus(u) = u > 0 ? u + log1p(exp(-u)) : log1p(exp(u))
_turner_logistic_integral(a, b, c, x) =
    a * c * (_turner_softplus((x - b) / c) - _turner_softplus(-b / c))

"""
    turner_egg_multiplier(p::TurnerParameters, x) -> Float64
    turner_juvenile_multiplier(p::TurnerParameters, x) -> Float64
    turner_tadpole_multiplier(p::TurnerParameters, x) -> Float64

The density-dependence curves of Turner et al. (2025) at density ratio `x`: the logistic
egg-laying curve `a/(exp((b − x)/c) + 1)`, the same form for juvenile and adult survival,
and the Bleasdale curve `(a + b x)^(−1/c)` for first-year survival (no caps applied).
"""
turner_egg_multiplier(p::TurnerParameters, x::Real) =
    _turner_logistic(p.egg_a, p.egg_b, p.egg_c, x)
turner_juvenile_multiplier(p::TurnerParameters, x::Real) =
    _turner_logistic(p.juvenile_a, p.juvenile_b, p.juvenile_c, x)
turner_tadpole_multiplier(p::TurnerParameters, x::Real) =
    (p.tadpole_a + p.tadpole_b * x)^(-1 / p.tadpole_c)
@doc (@doc turner_egg_multiplier) turner_juvenile_multiplier
@doc (@doc turner_egg_multiplier) turner_tadpole_multiplier

"""
    turner_egg_ceiling(p::TurnerParameters, eggs) -> Float64

Number of female eggs after the egg-laying density feedback (before rounding): with
`K = spawn_capacity × clutch_size_mean × sex_ratio` and `x = eggs/K`, eggs are multiplied by
the mean of the egg curve over `[0, x]` (at most 1) if `x ≤ 0.99`; above, the first `K` eggs
get the mean over `[0, 1.02]` and the rest the value at 1.02.
"""
function turner_egg_ceiling(p::TurnerParameters, eggs::Real)
    K = p.spawn_capacity * p.clutch_size_mean * p.sex_ratio
    x = eggs / K
    x > 0 || return float(eggs)
    a, b, c = p.egg_a, p.egg_b, p.egg_c
    if x <= p.egg_linear_limit
        area = min(_turner_logistic_integral(a, b, c, x) / x, 1.0)
        return eggs * area
    end
    x1 = p.egg_overshoot_point
    area1 = _turner_logistic_integral(a, b, c, x1) / x1
    return K * area1 + _turner_logistic(a, b, c, x1) * (eggs - K)
end

"""
    turner_watering_sequence(rng, pattern, years; p_wet = 0.7) -> BitVector

A sequence of watered (`true`) and dry years: year 1 is watered and each later year with
probability `p_wet`; with `pattern = :max2dry` (`:max1dry`) a year is forced wet after two
(one) consecutive dry years, `:random` has no cap. This rule is inferred from the archived
arrays of Turner et al. (2025); the code that drew them is not archived.
"""
function turner_watering_sequence(rng::AbstractRNG, pattern::Symbol, years::Integer;
        p_wet::Real = 0.7)
    pattern in TURNER_PATTERNS ||
        throw(ArgumentError("pattern must be one of $TURNER_PATTERNS; got :$pattern"))
    cap = pattern === :max2dry ? 2 : pattern === :max1dry ? 1 : typemax(Int)
    w = falses(years)
    w[1] = true
    run = 0
    for t in 2:years
        w[t] = run >= cap ? true : rand(rng) < p_wet
        run = w[t] ? 0 : run + 1
    end
    return w
end

# ---- data ------------------------------------------------------------------------------------

"""
    turner2025_path(file::Symbol = :survival_series)

Path of a bundled file of Turner et al. (2025) (Zenodo doi:10.5281/zenodo.10921748,
CC BY 4.0; see `data/turner2025/README.md`): `:survival_series` (the year-by-year survival
means and SDs that the coded model uses, for each watering month and dry-year pattern) or
`:start_populations` (1000 of the authors' burnt-in start populations per scenario).
"""
function turner2025_path(file::Symbol = :survival_series)
    names = (survival_series = "survival_series.csv",
             start_populations = "start_populations.csv")
    haskey(names, file) ||
        throw(ArgumentError("file must be one of $(keys(names)); got :$file"))
    return joinpath(pkgdir(@__MODULE__), "data", "turner2025", names[file])
end

function _turner_check(month, pattern)
    month in TURNER_MONTHS ||
        throw(ArgumentError("month must be one of $TURNER_MONTHS; got :$month"))
    pattern in TURNER_PATTERNS ||
        throw(ArgumentError("pattern must be one of $TURNER_PATTERNS; got :$pattern"))
    return nothing
end

"""
    load_turner2025(; month = :august, pattern = :random)

Bundled inputs of one scenario of Turner et al. (2025): `survival`, a 52 × 4 matrix of the
yearly adult survival mean and SD and metamorph summer survival mean and SD (columns in
that order) used by the coded model, and `start`, a 6 × 1000 matrix of start populations
(eggs, then females aged 1–5). `month` is the month of environmental water delivery
(`:august`, `:september`, `:october`) and `pattern` the dry-year rule (`:random`,
`:max2dry`, `:max1dry`).
"""
function load_turner2025(; month::Symbol = :august, pattern::Symbol = :random)
    _turner_check(month, pattern)
    m, pat = string(month), string(pattern)
    rows = [split(l, ',') for l in readlines(turner2025_path(:survival_series))[2:end]
            if !isempty(strip(l))]
    sel = sort!(filter(r -> r[1] == m && r[2] == pat, rows); by = r -> parse(Int, r[3]))
    survival = [parse(Float64, r[j]) for r in sel, j in 4:7]
    rows = [split(l, ',') for l in readlines(turner2025_path(:start_populations))[2:end]
            if !isempty(strip(l))]
    sel = filter(r -> r[1] == m && r[2] == pat, rows)
    start = [parse(Float64, r[j]) for j in 4:9, r in sel]
    return (survival = survival, start = start)
end

"""
    load_turner2025_survival_array(month, pattern; path = turner2025_archive_path())

The full precomputed array of one scenario from the authors' archive
(`referenceData/<Month>WateringArray_<pattern>.csv`, 40 MB each, not bundled): a
10,000 × 52 × 5 array of wetness (1/0), adult survival mean and SD, and summer survival
mean and SD per replicate and year. Needs the archive extracted under
[`turner2025_archive_path`](@ref); throws an `ArgumentError` otherwise.
"""
function load_turner2025_survival_array(month::Symbol, pattern::Symbol;
        path::AbstractString = turner2025_archive_path())
    _turner_check(month, pattern)
    pat = Dict(:random => "RandomDryYears", :max2dry => "2DryYears", :max1dry => "1DryYears")
    file = joinpath(path, "referenceData",
        uppercasefirst(string(month)) * "WateringArray_" * pat[pattern] * ".csv")
    isfile(file) || throw(ArgumentError("$file not found; extract TurnerBdModel.zip " *
        "(Zenodo doi:10.5281/zenodo.10921748) there or set CHYTRIDIOMYCOSIS_TURNER2025"))
    lines = readlines(file)
    n = length(lines) - 1
    A = Array{Float64,3}(undef, n, 52, 5)
    for (i, l) in enumerate(lines[2:end])
        f = split(l, ',')
        for k in 1:5, t in 1:52
            A[i, t, k] = parse(Float64, f[(k-1)*52+t])
        end
    end
    return A
end

"""
    turner2025_archive_path()

Directory of the extracted authors' archive `TurnerBdModel/` (Zenodo
doi:10.5281/zenodo.10921748): `ENV["CHYTRIDIOMYCOSIS_TURNER2025"]` if set, else
`lab/papers/turner2025/zenodo/TurnerBdModel` beside the package checkout.
"""
turner2025_archive_path() = get(ENV, "CHYTRIDIOMYCOSIS_TURNER2025",
    normpath(joinpath(pkgdir(@__MODULE__), "..", "..", "papers", "turner2025", "zenodo",
        "TurnerBdModel")))

# ---- problem and solution --------------------------------------------------------------------

"""
    TurnerPVAProblem(p::TurnerParameters; survival, start, pattern = :random,
                     wetness = nothing, replicates = 10_000, as_coded = true)
    TurnerPVAProblem(:turner2025; month = :august, pattern = :random, replicates = 10_000,
                     survival_array = nothing)

The stochastic age-classified population viability model of Turner et al. (2025) for
*Litoria raniformis* under one environmental-watering scenario. Inputs:

- `survival`: a `years × 4` matrix (adult survival mean, SD, metamorph summer survival
  mean, SD per year) shared by all replicates (the authors' code as run), or a
  `replicates × years × 5` array with wetness in slice 1 (one row per replicate, as
  intended; `as_coded = false`);
- `start`: a `6 × k` matrix of candidate start populations (eggs, ages 1–5), resampled with
  replacement from the non-empty ones;
- `wetness`: a `replicates × years` `BitMatrix` of watered years; `nothing` draws them with
  [`turner_watering_sequence`](@ref) and `pattern`.

The keyed constructor loads the bundled scenario ([`load_turner2025`](@ref)); pass the full
array from [`load_turner2025_survival_array`](@ref) as `survival_array` to use per-replicate
survival and the authors' own wetness sequences.

`solve(prob; rng = Random.default_rng())` returns a [`TurnerPVASolution`](@ref). Bd enters
only through the precomputed survival series.
"""
struct TurnerPVAProblem
    p::TurnerParameters
    survival::Array{Float64}
    start::Matrix{Float64}
    pattern::Symbol
    wetness::Union{Nothing,BitMatrix}
    replicates::Int
    as_coded::Bool
end

function TurnerPVAProblem(p::TurnerParameters; survival::AbstractArray{<:Real},
        start::AbstractMatrix{<:Real}, pattern::Symbol = :random,
        wetness::Union{Nothing,AbstractMatrix{Bool}} = nothing,
        replicates::Integer = 10_000, as_coded::Bool = true)
    pattern in TURNER_PATTERNS ||
        throw(ArgumentError("pattern must be one of $TURNER_PATTERNS; got :$pattern"))
    size(start, 1) == 6 || throw(ArgumentError("start must have 6 rows"))
    if as_coded
        ndims(survival) == 2 && size(survival) == (p.years, 4) ||
            throw(ArgumentError("as_coded needs a $(p.years) × 4 survival matrix"))
    else
        ndims(survival) == 3 && size(survival, 1) >= replicates &&
            size(survival, 2) == p.years && size(survival, 3) == 5 ||
            throw(ArgumentError("as_coded = false needs a replicates × $(p.years) × 5 array"))
    end
    keep = [j for j in axes(start, 2) if sum(view(start, :, j)) > 0]
    isempty(keep) && throw(ArgumentError("no non-empty start population"))
    w = wetness === nothing ? nothing : BitMatrix(wetness)
    return TurnerPVAProblem(p, Array{Float64}(survival), Matrix{Float64}(start[:, keep]),
        pattern, w, Int(replicates), as_coded)
end

function TurnerPVAProblem(key::Symbol; month::Symbol = :august, pattern::Symbol = :random,
        replicates::Integer = 10_000, survival_array = nothing)
    d = load_turner2025(; month, pattern)
    p = TurnerParameters(key)
    survival_array === nothing &&
        return TurnerPVAProblem(p; survival = d.survival, start = d.start, pattern,
            replicates)
    A = survival_array
    return TurnerPVAProblem(p; survival = A, start = d.start, pattern, replicates,
        wetness = A[1:replicates, :, 1] .== 1, as_coded = false)
end

"""
    TurnerPVASolution

Result of `solve(::TurnerPVAProblem)`: `population`, a `(years + 1) × replicates` matrix of
females aged ≥ 1 (row 1 the start), and `terminal`, the `6 × replicates` final age
structure (eggs first). See [`turner_extinction_probability`](@ref).
"""
struct TurnerPVASolution
    population::Matrix{Float64}
    terminal::Matrix{Float64}
end

"""
    turner_extinction_probability(sol::TurnerPVASolution) -> Float64

Fraction of replicates with no female aged ≥ 1 at the end of the projection (the authors'
`mean(term.pop == 0)`).
"""
turner_extinction_probability(sol::TurnerPVASolution) = mean(==(0), view(sol.population, size(sol.population, 1), :))

"""
    turner_minimum_population(sol::TurnerPVASolution) -> Vector{Float64}

Minimum number of females aged ≥ 1 over the projection, per replicate (the quantity of the
paper's cumulative distribution plots).
"""
turner_minimum_population(sol::TurnerPVASolution) = vec(minimum(sol.population; dims = 1))

"""
    turner_quasi_extinction_probability(sol::TurnerPVASolution, threshold) -> Float64

Fraction of replicates whose minimum number of females aged ≥ 1 falls below `threshold`
(the paper's Fig. 3 uses 20, `parameter_set(:turner2025).quasi_extinction_threshold`).
"""
turner_quasi_extinction_probability(sol::TurnerPVASolution, threshold::Real) =
    mean(<(threshold), turner_minimum_population(sol))

const _TURNER2025_RUNS = Dict{Tuple{Symbol,Symbol},TurnerPVASolution}()

# One seeded run of 10,000 replicates per scenario and session, for the validation targets.
_turner2025_reference(month::Symbol, pattern::Symbol = :random) =
    get!(() -> solve(TurnerPVAProblem(:turner2025; month, pattern);
                     rng = Random.Xoshiro(2025)), _TURNER2025_RUNS, (month, pattern))

function CommonSolve.solve(prob::TurnerPVAProblem; rng::AbstractRNG = Random.default_rng())
    p = prob.p
    T = p.years
    R = prob.replicates
    pop = zeros(T + 1, R)
    term = zeros(6, R)
    betad = Beta(p.tadpole_alpha, p.tadpole_beta)
    winter = p.winter_daily_survival^p.winter_duration
    M = zeros(6, 6)
    n = zeros(6)
    nn = zeros(6)
    ncand = size(prob.start, 2)
    for r in 1:R
        n .= view(prob.start, :, rand(rng, 1:ncand))
        pop[1, r] = sum(view(n, 2:6))
        wet = prob.wetness === nothing ?
              turner_watering_sequence(rng, prob.pattern, T; p_wet = p.p_wet) :
              view(prob.wetness, r, :)
        for t in 1:T
            sum(n) == 0 && break
            if prob.as_coded
                am, asd, sm, ssd = prob.survival[t, 1], prob.survival[t, 2],
                                   prob.survival[t, 3], prob.survival[t, 4]
            else
                am, asd, sm, ssd = prob.survival[r, t, 2], prob.survival[r, t, 3],
                                   prob.survival[r, t, 4], prob.survival[r, t, 5]
            end
            hatch = min(p.hatch_mean + p.hatch_sd * randn(rng), 1.0)
            tomet = rand(rng, betad) * hatch
            summ = sm + ssd * randn(rng)
            s0 = max(tomet * summ * winter, 0.0)
            fill!(M, 0.0)
            M[2, 1] = s0
            for k in 1:5
                s = max(am + asd * randn(rng), 0.0)
                k <= 4 && (M[k+2, k+1] = s)          # the 5th draw is unused, as coded
            end
            for k in 2:6
                f = max(p.fecundity_mean + p.fecundity_sd * randn(rng), 0.0)
                M[1, k] = round(f * p.sex_ratio)
            end
            # density dependence from the numbers at the start of the year
            xt = min(n[1] / p.tadpole_capacity, p.tadpole_ratio_cap)
            M[2, 1] *= turner_tadpole_multiplier(p, xt)
            xj = min(n[2] / p.juvenile_capacity, p.juvenile_ratio_cap)
            M[3, 2] *= turner_juvenile_multiplier(p, xj)
            xa = min(sum(view(n, 2:6)) / p.adult_capacity, p.juvenile_ratio_cap)
            fa = turner_juvenile_multiplier(p, xa)
            for k in 4:6
                M[k, k-1] *= fa
            end
            wet[t] || (M[1, :] .= 0.0)
            mul!(nn, M, n)
            for k in 2:6
                nn[k] = round(nn[k])
            end
            nn[1] = round(turner_egg_ceiling(p, nn[1]))
            n .= nn
            pop[t+1, r] = sum(view(n, 2:6))
        end
        term[:, r] .= n
    end
    return TurnerPVASolution(pop, term)
end
