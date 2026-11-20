# Grear DA, Adams MJ, Backlin AR, et al. (2025) Pan-amphibia distribution of the fungal
# parasite Batrachochytrium dendrobatidis varies with species and temperature. Ecological
# Monographs 95:e70001. https://doi.org/10.1002/ecm.70001 (open access, CC BY-NC 4.0)
# This implements a published model; it is not a new model.
#
# SOURCES. Read from the published full text (lab/papers/grear2025/grear2025.pdf, read
# 2026-11-20: Methods eqs 1–10, "Model fitting and diagnostics", "Imputation for missing
# data", Results, Table 2, Figs 4, 6, 7) and from the authors' code, the USGS software release
# doi:10.5066/P1MXUJZM (CC0; lab/papers/grear2025/code/src, commit b279e8f), whose data file is
# bundled in data/grear2025.
#
# Model (eqs 1–10). A hurdle model of the estimated ITS copy number y per swab:
#   P(y = 0) = 1 − π,  y > 0: π · LogNormal(y | μ, σ)                       (1, 6)
#   logit π = α0 + α1 temp + α2 temp² + α3 precip + α4 sex + α5 age + η_site + ζ_region + ε_spp (2)
#   log-mean μ = β0 + β1 temp + β2 temp² + β3 precip + β4 sex + β5 age + ξ_site + ρ_region + θ_spp (7)
#   site and region effects Normal(0, σ²), species effects Normal(0, σ² Σ_spp) (3–5, 8–10)
# temp: mean daily maximum air temperature over the 30 days before sampling (PRISM, °C); precip:
# total precipitation over the 10 days before sampling (mm); sex: Female (reference), Male,
# undetermined; age: Adult (reference), Juvenile, Larval. Priors: α ~ normal(0, 1),
# β ~ normal(0, 20), half-Student-t(3, 0, 2.5) on "the lognormal variance σ²" and all random
# effects. Fitted in brms (`hurdle_lognormal`, whose hu = 1 − π, so α = −(hu coefficient), as
# 09_Combine_ImputeFits.R transforms them for Figs 6–7). Missing sex (159), life stage (91) and
# genus-only Taricha (42) were multiply imputed: 30 data sets, by-genus multinomial draws of age
# and sex, Taricha sp. assigned T. torosa or T. granulosa with equal probability; each fitted
# with one of 30 phylogenies (Jetz & Pyron 2018), posteriors pooled.
#
# Code formula (04/06/09/11): copyNum ~ Age + Sex + m30max + I(m30max^2) + t10precip
# + (1 | reg/siteID) + (1 | gr(Gesp, cov = phyvcv)), the same in hu; it agrees with eqs 2, 7.
#
# Paper versus code (switches):
# - σ prior. The paper puts half-t(3, 0, 2.5) on the variance σ²; the code leaves brms' default,
#   which is on σ (student_t(3, 0, 2.5); 2.5 because the MAD of log y among positives is 2.42,
#   below 2.5). `sigma_prior = :as_coded` (default) or `:as_printed`.
# - Phylogenies. The paper uses 30 trees, one per imputed data set; the code's file has 100
#   trees and its loop is broken (it calls vcv on an undefined `tree`). `grear2025_phylogeny`
#   reads one tree (`tree = 1` by default).
#
# Choices where neither says (documented):
# - Missing covariates: `missing_covariates = :impute` (default) draws ONE imputed data set as
#   the code does (08_ImputeData.R: for each genus, the observed frequencies of Age×Sex
#   combinations; Spea, with no complete record, gets Adult/Juvenile/Larval with sex Unknown);
#   `:drop` drops the incomplete swabs. The paper pooled 30 imputations.
# - Sample: swabs with a date and coordinates (11,172; the code's own filter). The paper says
#   11,175 after excluding 8 with missing site or date; 11 records here lack a site.
# - A posterior mode is not the posterior the paper summarised. Sites (and species) whose swabs
#   are all negative separate: their hurdle effects and the random-effect standard deviations
#   keep growing, so L-BFGS does not converge to a finite mode (the gradient stays large after
#   thousands of iterations). The fixed effects settle quickly, which is what the targets use;
#   the full NUTS fit is scripts/grear2025_nuts.jl.

"""
    grear2025_path() -> String

Directory of the bundled Grear et al. (2025) swab data (public domain; see
`data/grear2025/README.md`).
"""
grear2025_path() = normpath(joinpath(@__DIR__, "..", "..", "data", "grear2025"))

"""
    Grear2025Data

Swab records of the US survey analysed by Grear et al. (2025), from [`load_grear2025`](@ref):
one entry per swab in each field: `date` (`missing` if not recorded), `latitude`,
`longitude`, `state`, `species` (`Gesp`), `order`, `age`, `sex`, `site` (`siteID`),
`region` (NOAA climate region), `result` (the qPCR result as recorded), `load` (`copyNum`,
Bd quantity per swab, `0` = not detected), `tmax30` (`m30max`, °C) and `precip10`
(`t10precip`, mm). Missing strings are `""` and missing numbers `NaN`.
"""
struct Grear2025Data
    date::Vector{Union{Missing,Date}}
    latitude::Vector{Float64}
    longitude::Vector{Float64}
    state::Vector{String}
    species::Vector{String}
    order::Vector{String}
    age::Vector{String}
    sex::Vector{String}
    site::Vector{String}
    region::Vector{String}
    result::Vector{String}
    load::Vector{Float64}
    tmax30::Vector{Float64}
    precip10::Vector{Float64}
end

Base.length(d::Grear2025Data) = length(d.load)

function Base.getindex(d::Grear2025Data, idx::AbstractVector)
    return Grear2025Data((getfield(d, f)[idx] for f in fieldnames(Grear2025Data))...)
end

# split one CSV line (fields may be quoted)
function _grear_split(line::AbstractString)
    out = String[]
    buf = IOBuffer()
    q = false
    for c in line
        if c == '"'
            q = !q
        elseif c == ',' && !q
            push!(out, String(take!(buf)))
        else
            write(buf, c)
        end
    end
    push!(out, String(take!(buf)))
    return out
end

"""
    load_grear2025(; path = grear2025_path(), complete = true) -> Grear2025Data

Read the bundled swab file. With `complete = true` (default) keep the swabs with a collection
date and coordinates, as the authors' scripts do first (11,172 swabs); `complete = false`
returns all 11,189 rows.
"""
function load_grear2025(; path::AbstractString = grear2025_path(), complete::Bool = true)
    lines = readlines(joinpath(path, "grear2025_swabs.csv"))
    hdr = _grear_split(lines[1])
    col = Dict(h => i for (i, h) in enumerate(hdr))
    rows = [_grear_split(l) for l in lines[2:end] if !isempty(l)]
    str(r, k) = (v = r[col[k]]; v == "NA" ? "" : v)
    num(r, k) = (v = r[col[k]]; v in ("NA", "") ? NaN : parse(Float64, v))
    dt(r) = (v = str(r, "CollectionDate"); isempty(v) ? missing : Date(v, dateformat"m/d/y"))
    d = Grear2025Data([dt(r) for r in rows], [num(r, "Latitude") for r in rows],
                      [num(r, "Longitude") for r in rows], [str(r, "State") for r in rows],
                      [str(r, "Gesp") for r in rows], [str(r, "order") for r in rows],
                      [str(r, "Age") for r in rows], [str(r, "Sex") for r in rows],
                      [str(r, "siteID") for r in rows], [str(r, "reg") for r in rows],
                      [str(r, "Bd.Result") for r in rows], [num(r, "copyNum") for r in rows],
                      [num(r, "m30max") for r in rows], [num(r, "t10precip") for r in rows])
    complete || return d
    keep = [d.date[i] !== missing && !isnan(d.latitude[i]) && !isnan(d.longitude[i])
            for i in 1:length(d)]
    return d[findall(keep)]
end

"""
    grear2025_summary(data = load_grear2025(complete = false)) -> NamedTuple

Counts for comparison with Grear et al. (2025): `swabs` (records with a collection date, the
authors' first filter), `species` (distinct species names among them, counting the
unresolved "Taricha sp." once), `sites` (distinct non-missing `siteID`s, the authors' 200 m
grouping), `model_swabs` (swabs with a date and coordinates, the authors' model data),
`missing_sex`, `missing_age`, `genus_only` (among `model_swabs`), `positive` (dated swabs with
`load > 0`), and `prevalence_caudata`, `prevalence_anura` (raw fraction of dated swabs with
`load > 0` by order).
"""
function grear2025_summary(d::Grear2025Data = load_grear2025(; complete = false))
    dated = findall(i -> d.date[i] !== missing, 1:length(d))
    dd = d[dated]
    sp = filter(!isempty, unique(dd.species))
    prev(o) = (i = findall(==(o), dd.order); count(>(0), dd.load[i]) / length(i))
    keep = findall(i -> d.date[i] !== missing && !isnan(d.latitude[i]) &&
                        !isnan(d.longitude[i]), 1:length(d))
    m = d[keep]
    return (swabs = length(dd), species = length(sp),
            sites = length(unique(filter(!isempty, dd.site))), model_swabs = length(m),
            missing_sex = count(isempty, m.sex), missing_age = count(isempty, m.age),
            genus_only = count(==("Taricha sp."), m.species),
            positive = count(>(0), dd.load), prevalence_caudata = prev("Caudata"),
            prevalence_anura = prev("Anura"))
end

"""
    grear2025_phylogeny_path() -> Union{String,Nothing}

Path of the phylogeny of the authors' software release (`outputSppInData.nex`, not
bundled): `ENV["CHYTRIDIOMYCOSIS_GREAR2025"]` (a directory holding it) or
`lab/papers/grear2025/code/src/` next to the package, else `nothing`.
"""
function grear2025_phylogeny_path()
    cands = String[]
    haskey(ENV, "CHYTRIDIOMYCOSIS_GREAR2025") && push!(cands, ENV["CHYTRIDIOMYCOSIS_GREAR2025"])
    d = pkgdir(@__MODULE__)
    for up in (("..", ".."), ("..", "..", ".."))
        push!(cands, normpath(joinpath(d, up..., "papers", "grear2025", "code", "src")))
    end
    for c in cands
        f = joinpath(c, "outputSppInData.nex")
        isfile(f) && return f
    end
    return nothing
end

"""
    grear2025_phylogeny(path = grear2025_phylogeny_path(); tree = 1) -> (; species, correlation)

The phylogenetic correlation matrix of tree number `tree` (of 100) in the authors' NEXUS file
(`ape::vcv(spTrees[[1]], corr = TRUE)` in their scripts uses the first), via
[`phylogenetic_correlation`](@ref). The paper used 30 trees, one per imputed data set.
Throws if the file is not available.
"""
function grear2025_phylogeny(path = grear2025_phylogeny_path(); tree::Integer = 1)
    path === nothing && throw(ArgumentError("phylogeny not found; set ENV[\"CHYTRIDIOMYCOSIS_GREAR2025\"]"))
    k = 0
    for l in eachline(path)
        m = match(r"^\s*TREE\s+\S+\s*=\s*(.*)$"i, l)
        m === nothing && continue
        k += 1
        k == tree && return phylogenetic_correlation(m.captures[1])
    end
    throw(ArgumentError("tree $tree not found in $path ($k trees)"))
end

const GREAR2025_COVARIATES = ["AgeJuvenile", "AgeLarval", "SexMale", "SexUnknown", "m30max",
                              "m30max^2", "t10precip"]

"""
    grear2025_impute(data; rng = Random.Xoshiro(2025)) -> Grear2025Data

One imputed copy of the swab data, as in the authors' 08_ImputeData.R (Methods,
"Imputation for missing data"): a swab missing Age and/or Sex gets a draw from the observed
frequencies of the Age×Sex combinations of its genus (keeping the recorded one of the two
when only one is missing); in a genus with no complete record (Spea) Age is Adult, Juvenile or
Larval with equal probability and Sex Unknown. Swabs recorded as "Taricha sp." become
*T. torosa* or *T. granulosa* with equal probability. The paper pooled fits to 30 such data
sets.
"""
function grear2025_impute(d::Grear2025Data; rng::AbstractRNG = Random.Xoshiro(2025))
    genus(s) = first(split(s, ' '))
    combos = Dict{String,Vector{Tuple{String,String}}}()
    for i in 1:length(d)
        (isempty(d.age[i]) || isempty(d.sex[i]) || isempty(d.species[i])) && continue
        push!(get!(combos, genus(d.species[i]), Tuple{String,String}[]), (d.age[i], d.sex[i]))
    end
    age, sex, sp = copy(d.age), copy(d.sex), copy(d.species)
    for i in 1:length(d)
        if sp[i] == "Taricha sp."
            sp[i] = rand(rng) < 0.5 ? "Taricha torosa" : "Taricha granulosa"
        end
        (isempty(age[i]) || isempty(sex[i])) || continue
        pool = get(combos, genus(sp[i]), nothing)
        a, s = pool === nothing ? (rand(rng, ("Adult", "Juvenile", "Larval")), "Unknown") :
               pool[rand(rng, 1:length(pool))]
        isempty(age[i]) && (age[i] = a)
        isempty(sex[i]) && (sex[i] = s)
    end
    return Grear2025Data(d.date, d.latitude, d.longitude, d.state, sp, d.order, age, sex,
                         d.site, d.region, d.result, d.load, d.tmax30, d.precip10)
end

"""
    Grear2025Problem(data = load_grear2025(); random = (:species, :site, :region),
                     phylogeny = nothing, missing_covariates = :impute,
                     rng = Random.Xoshiro(2025), sigma_prior = :as_coded)

The hurdle-lognormal model of Grear et al. (2025, eqs 1–10) as a
[`HurdleRegressionProblem`](@ref) (log density with gradient, LogDensityProblems interface).
Both parts have the fixed effects Age (Juvenile, Larval; reference Adult), Sex (Male,
Unknown; reference Female), 30-day mean maximum temperature and its square, and 10-day
precipitation, with the printed priors (normal(0, 20) load part, normal(0, 1) prevalence
part; compendium `grear2025`). `random` picks the random intercepts: `:species` (correlated
by `phylogeny`, a `(; species, correlation)` from [`grear2025_phylogeny`](@ref), or
independent when `nothing`), `:site` and `:region`.

`missing_covariates = :impute` (default) fills missing Age/Sex and "Taricha sp." with one
draw of [`grear2025_impute`](@ref) using `rng`; `:drop` drops those swabs. Swabs of a species
absent from `phylogeny` are dropped. `sigma_prior = :as_coded` puts the half-Student-t(3, 0,
2.5) prior on σ (brms' default, as in the authors' code); `:as_printed` puts it on σ² as the
paper states. The swab data used are in `prob.data`.
"""
struct Grear2025Problem
    model::HurdleRegressionProblem
    data::Grear2025Data
    levels::Dict{Symbol,Vector{String}}
end

function Grear2025Problem(d::Grear2025Data = load_grear2025();
                          random = (:species, :site, :region), phylogeny = nothing,
                          missing_covariates::Symbol = :impute,
                          rng::AbstractRNG = Random.Xoshiro(2025),
                          sigma_prior::Symbol = :as_coded)
    missing_covariates in (:impute, :drop) ||
        throw(ArgumentError("missing_covariates must be :impute or :drop"))
    sigma_prior in (:as_coded, :as_printed) ||
        throw(ArgumentError("sigma_prior must be :as_coded or :as_printed"))
    missing_covariates === :impute && (d = grear2025_impute(d; rng))
    ok = [!isempty(d.age[i]) && !isempty(d.sex[i]) && !isempty(d.species[i]) &&
          d.species[i] != "Taricha sp." && !isnan(d.load[i]) &&
          (phylogeny === nothing || d.species[i] in phylogeny.species) for i in 1:length(d)]
    dd = d[findall(ok)]
    X = hcat(dd.age .== "Juvenile", dd.age .== "Larval", dd.sex .== "Male",
             dd.sex .== "Unknown", dd.tmax30, dd.tmax30 .^ 2, dd.precip10)
    groups = NamedTuple[]
    levels = Dict{Symbol,Vector{String}}()
    for r in random
        if r === :species
            lv = phylogeny === nothing ? sort(unique(dd.species)) : phylogeny.species
            C = phylogeny === nothing ? nothing : phylogeny.correlation
        elseif r === :site
            lv, C = sort(unique(dd.site)), nothing
        elseif r === :region
            lv, C = sort(unique(dd.region)), nothing
        else
            throw(ArgumentError("unknown random effect :$r"))
        end
        field = r === :species ? dd.species : r === :site ? dd.site : dd.region
        pos = Dict(v => k for (k, v) in enumerate(lv))
        push!(groups, (name = r, index = [pos[v] for v in field], levels = length(lv),
                       correlation = C))
        levels[r] = lv
    end
    v = parameter_set(:grear2025)
    m = HurdleRegressionProblem(dd.load, Float64.(X), groups; names = GREAR2025_COVARIATES,
                                prior_b = v.prior_sd_load, prior_a = v.prior_sd_load,
                                prior_b_hu = v.prior_sd_hurdle, prior_a_hu = v.prior_sd_hurdle,
                                prior_scale = v.prior_scale_sd, prior_scale_hu = v.prior_scale_sd,
                                sigma_prior = sigma_prior === :as_coded ? :sd : :variance)
    return Grear2025Problem(m, dd, levels)
end

LogDensityProblems.dimension(p::Grear2025Problem) = LogDensityProblems.dimension(p.model)
LogDensityProblems.capabilities(::Type{<:Grear2025Problem}) = LogDensityProblems.LogDensityOrder{1}()
LogDensityProblems.logdensity(p::Grear2025Problem, θ::AbstractVector) =
    LogDensityProblems.logdensity(p.model, θ)
LogDensityProblems.logdensity_and_gradient(p::Grear2025Problem, θ::AbstractVector) =
    LogDensityProblems.logdensity_and_gradient(p.model, θ)
