# Grear DA, Adams MJ, Backlin AR, et al. (2025) Pan-amphibia distribution of the fungal
# parasite Batrachochytrium dendrobatidis varies with species and temperature. Ecological
# Monographs 95:e70001. https://doi.org/10.1002/ecm.70001 (open access, CC BY-NC 4.0)
# This implements a published model; it is not a new model.
#
# SOURCES. The full text could NOT be read: the publisher's site refuses automated clients,
# and it is not in PubMed Central. What is implemented is read from (1) the abstract (USGS
# Publications Warehouse record 70263881) and (2) the authors' code, the USGS software release
# doi:10.5066/P1MXUJZM (CC0; lab/papers/grear2025/code/src, commit b279e8f), whose data file is
# bundled in data/grear2025. The printed equations, priors and coefficient tables have not
# been seen, so no printed coefficient is entered in the compendium or reproduced; those
# validation targets are pending.
#
# Model (04_tempPrecipModFit.R, 06_randomEffectsFitting.R, 09_Combine_ImputeFits.R,
# 11_OutOfSampleFit.R). brms `hurdle_lognormal` of the Bd quantity per swab `copyNum`, with
# the same formula for the log-load part (mu) and the hurdle part (hu = probability of a
# zero):
#
#   copyNum ~ Age + Sex + m30max + I(m30max^2) + t10precip
#             + (1 | reg/siteID) + (1 | gr(Gesp, cov = phyvcv))
#
# m30max: mean daily maximum air temperature over the 30 days before sampling (PRISM, °C);
# t10precip: total precipitation over the 10 days before sampling (mm); reg: NOAA climate
# region; Gesp: species, with a phylogenetic correlation matrix from one tree of Jetz & Pyron
# (2018) (ape::vcv(tree, corr = TRUE)). Priors: normal(0, 20) on the mu coefficients and
# intercept, normal(0, 1) on the hu coefficients and intercept, brms defaults on σ and the
# random-effect standard deviations. The prevalence coefficients plotted in the paper (α) are
# logit(1 − invlogit(hu coefficient)) = −(hu coefficient) (09_Combine_ImputeFits.R).
# The temperature metric and lag (30-day maximum) and precipitation lag (10 days) were chosen
# by LOO among many fits (04/05); the random-effect structure by LOO/K-fold (06/07).
#
# Reconstructions and choices (the code is not fully consistent, and the paper is unseen):
# - The fixed-effect formula of the final model is not written out in the release (the
#   fitting section of 08_ImputeData.R is empty). It is inferred from the parameter indices
#   and labels of the Figure 6/7 code in 09 (b_Intercept, b_hu_Intercept, then Age Juvenile,
#   Larval, Sex Male, Unknown, m30max, m30max², t10precip for mu, then the same for hu) and
#   from the covariates passed to `predict` in 11. Reference levels Adult and Female (R's
#   alphabetical default).
# - Random effects: 11 predicts with (1|reg/siteNameID) + (1|gr(Gesp, cov=phyvcv)); siteNameID
#   is not a column of the data file, and 04/06 use siteID. We use reg and siteID nested in
#   reg (siteID is unique across regions, so the nesting is (1|reg) + (1|siteID)).
# - Missing Age/Sex (and the unresolved "Taricha sp.") were multiply imputed (30 data sets,
#   08) and the 30 posteriors combined. We drop those swabs instead (`missing_covariates = :drop`).
# - brms' default priors on σ and on the standard deviations are taken to be
#   half-Student-t(3, 0, 2.5) for both parts. brms may widen the scale to the MAD of the
#   response on the link scale; for log copyNum among positives that MAD is below 2.5 here
#   (see test/grear2025.jl), so the scale is 2.5 either way.
# - A MAP (posterior mode) is not the posterior the paper summarised: the full fit (≈1,300
#   parameters, 4 chains of NUTS in Stan) is left to scripts/grear2025_nuts.jl with any
#   LogDensityProblems sampler; the tests use the mode.

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

Counts for comparison with the abstract of Grear et al. (2025): `swabs` (records with a
collection date, the authors' first filter), `species` (distinct species names among them,
counting the unresolved "Taricha sp." once), `sites` (distinct `siteID`s), `positive` (swabs
with `load > 0`), and `prevalence_caudata`, `prevalence_anura` (raw fraction of swabs with
`load > 0` by order).
"""
function grear2025_summary(d::Grear2025Data = load_grear2025(; complete = false))
    dated = findall(i -> d.date[i] !== missing, 1:length(d))
    dd = d[dated]
    sp = filter(!isempty, unique(dd.species))
    prev(o) = (i = findall(==(o), dd.order); count(>(0), dd.load[i]) / length(i))
    return (swabs = length(dd), species = length(sp), sites = length(unique(dd.site)),
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
    grear2025_phylogeny(path = grear2025_phylogeny_path()) -> (; species, correlation)

The phylogenetic correlation matrix of the first tree in the authors' NEXUS file (as
`ape::vcv(spTrees[[1]], corr = TRUE)` in their scripts), via
[`phylogenetic_correlation`](@ref). Throws if the file is not available.
"""
function grear2025_phylogeny(path = grear2025_phylogeny_path())
    path === nothing && throw(ArgumentError("phylogeny not found; set ENV[\"CHYTRIDIOMYCOSIS_GREAR2025\"]"))
    for l in eachline(path)
        m = match(r"^\s*TREE\s+\S+\s*=\s*(.*)$"i, l)
        m === nothing || return phylogenetic_correlation(m.captures[1])
    end
    throw(ArgumentError("no TREE line in $path"))
end

const GREAR2025_COVARIATES = ["AgeJuvenile", "AgeLarval", "SexMale", "SexUnknown", "m30max",
                              "m30max^2", "t10precip"]

"""
    Grear2025Problem(data = load_grear2025(); random = (:species, :site, :region),
                     phylogeny = nothing, missing_covariates = :drop)

The hurdle-lognormal model of Grear et al. (2025) as a [`HurdleRegressionProblem`](@ref)
(log density with gradient, LogDensityProblems interface). Both parts have the fixed effects
Age (Juvenile, Larval; reference Adult), Sex (Male, Unknown; reference Female), 30-day mean
maximum temperature and its square, and 10-day precipitation, with the authors' priors
(normal(0, 20) load part, normal(0, 1) hurdle part; compendium `grear2025`). `random` picks the random intercepts:
`:species` (correlated by `phylogeny`, a `(; species, correlation)` from
[`grear2025_phylogeny`](@ref), or independent when `nothing`), `:site` and `:region`. Swabs
with missing Age or Sex, an unresolved species or a species absent from the phylogeny are
dropped (`missing_covariates = :drop`; the authors imputed them). The swab data are in `prob.data`.
"""
struct Grear2025Problem
    model::HurdleRegressionProblem
    data::Grear2025Data
    levels::Dict{Symbol,Vector{String}}
end

function Grear2025Problem(d::Grear2025Data = load_grear2025();
                          random = (:species, :site, :region), phylogeny = nothing,
                          missing_covariates::Symbol = :drop)
    missing_covariates === :drop || throw(ArgumentError("only missing_covariates = :drop is implemented"))
    ok = [!isempty(d.age[i]) && !isempty(d.sex[i]) && !isempty(d.species[i]) &&
          d.species[i] != "Taricha sp." && !isnan(d.load[i]) &&
          (phylogeny === nothing || d.species[i] in phylogeny.species) for i in 1:length(d)]
    dd = d[findall(ok)]
    n = length(dd)
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
                                prior_b_hu = v.prior_sd_hurdle, prior_a_hu = v.prior_sd_hurdle)
    return Grear2025Problem(m, dd, levels)
end

LogDensityProblems.dimension(p::Grear2025Problem) = LogDensityProblems.dimension(p.model)
LogDensityProblems.capabilities(::Type{<:Grear2025Problem}) = LogDensityProblems.LogDensityOrder{1}()
LogDensityProblems.logdensity(p::Grear2025Problem, θ::AbstractVector) =
    LogDensityProblems.logdensity(p.model, θ)
LogDensityProblems.logdensity_and_gradient(p::Grear2025Problem, θ::AbstractVector) =
    LogDensityProblems.logdensity_and_gradient(p.model, θ)
