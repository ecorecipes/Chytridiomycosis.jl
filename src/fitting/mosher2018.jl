# Maximum likelihood fit of the Mosher et al. (2018) dynamic two-species occupancy model with
# covariates (model and choices: src/models/mosher2018.jl). Program MARK's "dynamic
# conditional-binomial 2-species occupancy" likelihood is reproduced with the generic layer
# in src/observation/occupancy.jl; the maximisation is BFGS (finite-difference gradients)
# on the logit scale, with Wald standard errors from a finite-difference Hessian.
#
# Data: Dryad doi:10.5061/dryad.1c990 (CC0), bundled in data/mosher2018. Layout (checked
# 2026-11-20): boreal_toad_bd_histories.csv has a row of years (2001-2010, one every 20
# columns), a header row ("Site", then "Toad j", "Bd j" for surveys j = 1..10 of each year)
# and one row per site, with "1", "0" or "." (no survey; for Bd, no swab). The covariate file
# has elevation (m), permanence ("permanent"/"semi-permanent"), other amphibians ("yes"/"no"),
# region (NE, NW, SE, SW) and jdate1..jdate100 (day of year of survey j of year y in column
# 10(y−1)+j, 0.5 when there was no survey). 83 sites; Site 33 (never a toad detected) has no
# permanence or other-amphibian value. The paper analysed 82 sites, and dropping Site 33
# reproduces its AICc sample size (below), so load_mosher2018() drops it by default.
#
# Sample size for AICc. Table 3's AICc and deviance (−2 log L) imply n ≈ 494 for every model
# listed: the number of site-years with at least one survey at the 82 sites (494). We use that
# count, which is our inference of what Program MARK used, not something the paper states.

"""
    Mosher2018Data(; toad, bd, elevation = zeros(n), doy = missing,
                     other_amphibians = zeros(n), permanent = zeros(n),
                     region = fill("", n), sites = string.(1:n))

Detection histories and covariates for the Mosher et al. (2018) model. `toad` and `bd` are
`sites × seasons × surveys` arrays of `1`, `0` or `missing` (toad detected; Bd detected on a
swab; `missing` for no survey, or for Bd, no swab). Site covariates: `elevation` (m),
`other_amphibians` (1 if tiger salamanders or boreal chorus frogs were detected at the
site), `permanent` (1 permanent, 0 semi-permanent) and `region` (a label). `doy` is the day
of year of each survey (same shape as `toad`, `missing` allowed).
"""
struct Mosher2018Data
    toad::Array{Union{Missing,Int},3}
    bd::Array{Union{Missing,Int},3}
    elevation::Vector{Float64}
    doy::Array{Union{Missing,Float64},3}
    other_amphibians::Vector{Float64}
    permanent::Vector{Float64}
    region::Vector{String}
    sites::Vector{String}
end

function Mosher2018Data(; toad, bd, elevation = nothing, doy = missing,
                        other_amphibians = nothing, permanent = nothing, region = nothing,
                        sites = nothing)
    size(toad) == size(bd) || throw(DimensionMismatch("toad and bd must have the same size"))
    n = size(toad, 1)
    D = doy === missing ? fill!(Array{Union{Missing,Float64}}(undef, size(toad)), missing) :
        convert(Array{Union{Missing,Float64},3}, doy)
    size(D) == size(toad) || throw(DimensionMismatch("doy must have the size of toad"))
    vec_or(x, d) = x === nothing ? d : collect(x)
    d = Mosher2018Data(convert(Array{Union{Missing,Int},3}, toad),
                       convert(Array{Union{Missing,Int},3}, bd),
                       float.(vec_or(elevation, zeros(n))), D,
                       float.(vec_or(other_amphibians, zeros(n))),
                       float.(vec_or(permanent, zeros(n))),
                       string.(vec_or(region, fill("", n))),
                       string.(vec_or(sites, string.(1:n))))
    all(v -> length(v) == n, (d.elevation, d.other_amphibians, d.permanent, d.region, d.sites)) ||
        throw(DimensionMismatch("site covariates must have one entry per site"))
    return d
end

"""
    mosher2018_path() -> String

Directory holding the Mosher et al. (2018) Dryad files `boreal_toad_bd_histories.csv` and
`boreal_toad_bd_covariates.csv` (doi:10.5061/dryad.1c990, CC0):
`ENV["CHYTRIDIOMYCOSIS_MOSHER2018"]` if set, otherwise the copy bundled in
`data/mosher2018`.
"""
function mosher2018_path()
    p = _mosher2018_find()
    p === nothing && throw(ArgumentError(
        "Mosher et al. 2018 data not found in ENV[\"CHYTRIDIOMYCOSIS_MOSHER2018\"]"))
    return p
end

const MOSHER2018_FILES = ("boreal_toad_bd_histories.csv", "boreal_toad_bd_covariates.csv")

function _mosher2018_find()
    ok(d) = all(f -> isfile(joinpath(d, f)), MOSHER2018_FILES)
    e = get(ENV, "CHYTRIDIOMYCOSIS_MOSHER2018", "")
    isempty(e) || return isdir(e) && ok(e) ? e : nothing
    d = joinpath(pkgdir(@__MODULE__), "data", "mosher2018")
    return ok(d) ? d : nothing
end

function _m18_csv(path)
    rows = [String.(split(l, ',')) for l in eachline(path) if !isempty(strip(l))]
    return [strip.(r) for r in rows]
end

"""
    load_mosher2018(; dir = mosher2018_path(), drop_incomplete = true) -> Mosher2018Data

Read the Mosher et al. (2018) Dryad deposit: 83 sites, 2001-2010, up to 10 surveys a year,
toad detections and Bd swab results (`missing` when no survey or no swab), elevation (m),
permanence (1 permanent), other amphibians (1 present), region and day of year of each
survey. With `drop_incomplete = true` (the default) sites without a permanence or
other-amphibian value are dropped: that is Site 33, which leaves the paper's 82 sites.
"""
function load_mosher2018(; dir::AbstractString = mosher2018_path(), drop_incomplete::Bool = true)
    H = _m18_csv(joinpath(dir, MOSHER2018_FILES[1]))
    Cv = _m18_csv(joinpath(dir, MOSHER2018_FILES[2]))
    H[2][1] == "Site" && H[2][2] == "Toad 1" && H[2][3] == "Bd 1" && length(H[2]) == 201 ||
        throw(ArgumentError("unexpected header in $(MOSHER2018_FILES[1])"))
    [H[1][2 + 20k] for k in 0:9] == string.(2001:2010) ||
        throw(ArgumentError("unexpected year row in $(MOSHER2018_FILES[1])"))
    Cv[1][2:5] == ["Elevation (m)", "Permanence", "Other Amphibians", "Region"] &&
        length(Cv[1]) == 105 || throw(ArgumentError("unexpected header in $(MOSHER2018_FILES[2])"))
    sites = [r[1] for r in H[3:end]]
    sites == [r[1] for r in Cv[2:end]] ||
        throw(ArgumentError("the two files list different sites"))
    n = length(sites)
    rec(x) = x == "." ? missing : x == "1" ? 1 : x == "0" ? 0 :
             throw(ArgumentError("unexpected record $x"))
    toad = Array{Union{Missing,Int}}(undef, n, 10, 10)
    bd = similar(toad)
    doy = Array{Union{Missing,Float64}}(missing, n, 10, 10)
    for i in 1:n, y in 1:10, j in 1:10
        v = H[i + 2]
        toad[i, y, j] = rec(v[1 + 20(y - 1) + 2j - 1])
        bd[i, y, j] = rec(v[1 + 20(y - 1) + 2j])
        d = parse(Float64, Cv[i + 1][5 + 10(y - 1) + j])
        # 0.5 marks no survey; a date without any record (one case) is ignored
        if d != 0.5 && !(toad[i, y, j] === missing && bd[i, y, j] === missing)
            doy[i, y, j] = d
        end
    end
    perm = [r[3] == "permanent" ? 1.0 : r[3] == "semi-permanent" ? 0.0 : NaN for r in Cv[2:end]]
    amph = [r[4] == "yes" ? 1.0 : r[4] == "no" ? 0.0 : NaN for r in Cv[2:end]]
    keep = drop_incomplete ? findall(i -> !isnan(perm[i]) && !isnan(amph[i]), 1:n) : collect(1:n)
    return Mosher2018Data(; toad = toad[keep, :, :], bd = bd[keep, :, :],
        elevation = [parse(Float64, Cv[i + 1][2]) for i in keep], doy = doy[keep, :, :],
        other_amphibians = amph[keep], permanent = perm[keep],
        region = [Cv[i + 1][5] for i in keep], sites = sites[keep])
end

_m18_mean(x) = sum(x) / length(x)
function _m18_sd(x)
    m = _m18_mean(x)
    s = sqrt(sum(abs2, x .- m) / max(length(x) - 1, 1))
    return s > 0 ? s : 1.0
end

# --- model structure -----------------------------------------------------------------------

const MOSHER2018_RATES = (:γA, :γAB, :γB, :γBA, :εA, :εAB, :εB, :εBA, :pA, :rA, :rBA)
const MOSHER2018_PAIRS = ((:γA, :γAB), (:γB, :γBA), (:εA, :εAB), (:εB, :εBA), (:pA, :rA))
const MOSHER2018_FORMULAS = (:const, :elev, :elev2, :amph, :perm, :region, :doy, :doy2)

"""
    Mosher2018Model(data::Mosher2018Data; additive = true, null = false, kw...)

The Mosher et al. (2018) likelihood for `data` with a covariate structure for each of the
rates `γA, γAB, γB, γBA, εA, εAB, εB, εBA, pA, rA, rBA`, given as keywords (default
`:const`). Each is one of `:const`, `:elev`, `:elev2` (linear and quadratic in elevation),
`:amph`, `:perm`, `:region` (site covariates) or, for `pA`, `rA`, `:doy`, `:doy2` (day of
year). Initial occupancy `ψA, ψBA, ψBa` is constant; `pB = rBa = 0`.

With `additive = true`, the members of a pair `(γA, γAB)`, `(γB, γBA)`, `(εA, εAB)`,
`(εB, εBA)`, `(pA, rA)` that have the same non-constant structure share their slopes (the
"additive relationships between parameter pairs" of the paper's Table 3). `null = true` sets each
colonisation or extinction pair equal (one constant per pair), the null model in which rates do not depend on the other
species. Covariates are standardised.

Implements the `LogDensityProblems` interface (the log-likelihood on the logit scale) and
[`parameter_names`](@ref); fit with [`maximum_likelihood`](@ref).
"""
struct Mosher2018Model
    data::Mosher2018Data
    formulas::Dict{Symbol,Symbol}
    intercept::Dict{Symbol,Int}
    slopes::Dict{Symbol,Vector{Tuple{Int,Symbol}}}
    names::Vector{String}
    elev::Vector{Float64}
    amph::Vector{Float64}
    perm::Vector{Float64}
    region::Dict{Symbol,Vector{Float64}}
    doy::Array{Float64,3}
    codes::Array{Union{Missing,Int},3}
end

function Mosher2018Model(data::Mosher2018Data; additive::Bool = true, null::Bool = false, kw...)
    f = Dict{Symbol,Symbol}(r => :const for r in MOSHER2018_RATES)
    for (k, v) in kw
        k in MOSHER2018_RATES || throw(ArgumentError("unknown rate $k"))
        v in MOSHER2018_FORMULAS || throw(ArgumentError("unknown structure $v for $k"))
        v in (:doy, :doy2) && !(k in (:pA, :rA)) &&
            throw(ArgumentError("day of year applies to pA and rA only"))
        f[k] = v
    end
    vital = MOSHER2018_RATES[1:8]
    null && any(r -> f[r] !== :const, vital) &&
        throw(ArgumentError("the null model has constant colonisation and extinction"))
    z(x) = (x .- _m18_mean(x)) ./ _m18_sd(x)
    std_or_one(x) = _m18_sd(x)
    elev, amph, perm = z(data.elevation), z(data.other_amphibians), z(data.permanent)
    regs = sort(unique(data.region))
    region = Dict(Symbol("region_", r) => Float64.(data.region .== r) for r in regs[2:end])
    dv = collect(skipmissing(data.doy))
    doy = if isempty(dv)
        zeros(size(data.doy))
    else
        m, s = sum(dv) / length(dv), std_or_one(dv)
        [ismissing(x) ? 0.0 : (x - m) / s for x in data.doy]
    end
    names = ["ψA", "ψBA", "ψBa"]
    intercept = Dict{Symbol,Int}()
    slopes = Dict{Symbol,Vector{Tuple{Int,Symbol}}}()
    terms(form) = form === :const ? Symbol[] : form === :elev ? [:elev] :
                  form === :elev2 ? [:elev, :elev2] : form === :doy ? [:doy] :
                  form === :doy2 ? [:doy, :doy2] : form === :region ? sort(collect(keys(region))) :
                  [form]
    pairof(r) = findfirst(pr -> r in pr, MOSHER2018_PAIRS)
    for r in MOSHER2018_RATES
        pi = pairof(r)
        partner = pi === nothing ? nothing : (MOSHER2018_PAIRS[pi][1] === r ?
                                               MOSHER2018_PAIRS[pi][2] : MOSHER2018_PAIRS[pi][1])
        if null && r in vital && haskey(intercept, partner)
            intercept[r] = intercept[partner]
            slopes[r] = Tuple{Int,Symbol}[]
            names[intercept[r]] *= "=" * string(r)
            continue
        end
        push!(names, string(r))
        intercept[r] = length(names)
        slopes[r] = Tuple{Int,Symbol}[]
        share = additive && partner !== nothing && haskey(slopes, partner) &&
                f[partner] === f[r] && f[r] !== :const
        for (k, term) in enumerate(terms(f[r]))
            if share
                idx = slopes[partner][k][1]
                names[idx] = replace(names[idx], string(partner, ":") => string(partner, ",", r, ":"))
                push!(slopes[r], (idx, term))
            else
                push!(names, string(r, ":", term))
                push!(slopes[r], (length(names), term))
            end
        end
    end
    codes = [two_species_code(data.toad[i], data.bd[i]) for i in CartesianIndices(data.toad)]
    return Mosher2018Model(data, f, intercept, slopes, names, elev, amph, perm, region, doy,
                           codes)
end

"""
    parameter_names(m::Mosher2018Model)

Names of the logit-scale coefficients of a [`Mosher2018Model`](@ref).
"""
parameter_names(m::Mosher2018Model) = m.names
LogDensityProblems.dimension(m::Mosher2018Model) = length(m.names)
LogDensityProblems.capabilities(::Type{Mosher2018Model}) = LogDensityProblems.LogDensityOrder{0}()
LogDensityProblems.logdensity(m::Mosher2018Model, θ::AbstractVector) = mosher2018_loglikelihood(m, θ)

function _m18_cov(m::Mosher2018Model, term::Symbol, i, t, j)
    term === :elev && return m.elev[i]
    term === :elev2 && return m.elev[i]^2
    term === :amph && return m.amph[i]
    term === :perm && return m.perm[i]
    term === :doy && return m.doy[i, t, j]
    term === :doy2 && return m.doy[i, t, j]^2
    return m.region[term][i]
end

function _m18_rate(m::Mosher2018Model, θ, r::Symbol, i, t = 1, j = 1)
    η = θ[m.intercept[r]]
    for (idx, term) in m.slopes[r]
        η += θ[idx] * _m18_cov(m, term, i, t, j)
    end
    return _expit(η)
end

"""
    mosher2018_site_parameters(m::Mosher2018Model, θ, site; season = 1, survey = 1)

The [`Mosher2018Parameters`](@ref) of one site at coefficients `θ` (detection at the
covariates of the given survey).
"""
function mosher2018_site_parameters(m::Mosher2018Model, θ::AbstractVector, i::Integer;
                                    season::Integer = 1, survey::Integer = 1)
    rate(r) = _m18_rate(m, θ, r, i, season, survey)
    z = zero(rate(:pA))
    return Mosher2018Parameters(; ψA = _expit(θ[1]), ψBA = _expit(θ[2]), ψBa = _expit(θ[3]),
        γA = rate(:γA), γAB = rate(:γAB), γB = rate(:γB), γBA = rate(:γBA), εA = rate(:εA),
        εAB = rate(:εAB), εB = rate(:εB), εBA = rate(:εBA), pA = rate(:pA), rA = rate(:rA),
        rBA = rate(:rBA), pB = z, rBa = z)
end

"""
    mosher2018_parameters_at(m::Mosher2018Model, θ; elevation, doy, other_amphibians = 0,
                             permanent = 1, region = nothing) -> Mosher2018Parameters

The [`Mosher2018Parameters`](@ref) at coefficients `θ` for covariate values in their
original units (elevation in m, day of year), standardised as in the model. `region` is a
region label (the first label alphabetically is the reference). Use it to draw the
elevation curves of the paper's Fig. 5.
"""
function mosher2018_parameters_at(m::Mosher2018Model, θ::AbstractVector; elevation::Real,
                                  doy::Real = _m18_mean(collect(skipmissing(m.data.doy))),
                                  other_amphibians::Real = 0, permanent::Real = 1,
                                  region = nothing)
    d = m.data
    zv(x, v) = (v - _m18_mean(x)) / _m18_sd(x)
    dv = collect(skipmissing(d.doy))
    cov = Dict(:elev => zv(d.elevation, elevation), :amph => zv(d.other_amphibians, other_amphibians),
               :perm => zv(d.permanent, permanent),
               :doy => isempty(dv) ? 0.0 : zv(dv, doy))
    cov[:elev2] = cov[:elev]^2
    cov[:doy2] = cov[:doy]^2
    for k in keys(m.region)
        cov[k] = region !== nothing && k === Symbol("region_", region) ? 1.0 : 0.0
    end
    function rate(r)
        η = θ[m.intercept[r]]
        for (idx, term) in m.slopes[r]
            η += θ[idx] * cov[term]
        end
        return _expit(η)
    end
    z = zero(rate(:pA))
    return Mosher2018Parameters(; ψA = _expit(θ[1]), ψBA = _expit(θ[2]), ψBa = _expit(θ[3]),
        γA = rate(:γA), γAB = rate(:γAB), γB = rate(:γB), γBA = rate(:γBA), εA = rate(:εA),
        εAB = rate(:εAB), εB = rate(:εB), εBA = rate(:εBA), pA = rate(:pA), rA = rate(:rA),
        rBA = rate(:rBA), pB = z, rBa = z)
end

"""
    mosher2018_loglikelihood(m::Mosher2018Model, θ) -> Real

Log-likelihood of all sites at logit-scale coefficients `θ` (ordered as
[`parameter_names`](@ref)).
"""
function mosher2018_loglikelihood(m::Mosher2018Model, θ::AbstractVector)
    φ0 = two_species_initial(_expit(θ[1]), _expit(θ[2]), _expit(θ[3]))
    s = zero(eltype(φ0))
    for i in axes(m.codes, 1)
        rate(r) = _m18_rate(m, θ, r, i)
        Φ = two_species_transition(; γA = rate(:γA), γAB = rate(:γAB), γB = rate(:γB),
                                   γBA = rate(:γBA), εA = rate(:εA), εAB = rate(:εAB),
                                   εB = rate(:εB), εBA = rate(:εBA))
        rBA = rate(:rBA)
        doy_dep = m.formulas[:pA] in (:doy, :doy2) || m.formulas[:rA] in (:doy, :doy2)
        D1 = two_species_detection(; pA = rate(:pA), rA = rate(:rA), rBA)
        emission = doy_dep ?
            ((t, j) -> two_species_detection(; pA = _m18_rate(m, θ, :pA, i, t, j),
                                             rA = _m18_rate(m, θ, :rA, i, t, j), rBA)) : D1
        s += dynamic_occupancy_loglikelihood(φ0, Φ, emission, view(m.codes, i, :, :))
    end
    return s
end

"""
    maximum_likelihood(m::Mosher2018Model; x0 = zeros(dim), iters = 2000) -> NamedTuple

Maximum likelihood fit of a [`Mosher2018Model`](@ref) by BFGS on the logit scale. Returns
`θ` (coefficients), `names`, `se` (Wald, `NaN` where the information is singular, e.g. at a
boundary), `loglik`, `deviance` (`−2 loglik`), `k` (number of coefficients), `n` (the
AICc sample size: site-years with at least one survey, our reading of Program MARK's
default), `aicc` and `converged`.
"""
function maximum_likelihood(m::Mosher2018Model; x0 = nothing, iters::Integer = 2000,
                            gtol::Real = 1e-5)
    k = length(m.names)
    x = x0 === nothing ? zeros(k) : collect(float.(x0))
    nll(θ) = -mosher2018_loglikelihood(m, θ)
    r = _ms_bfgs(nll, x; iters, gtol)
    V, _ = _ms_vcov(-_ms_hessian(nll, r.x))
    n = count(I -> any(!ismissing, view(m.codes, I[1], I[2], :)),
              CartesianIndices((size(m.codes, 1), size(m.codes, 2))))
    ll = -r.f
    aicc = -2ll + 2k + 2k * (k + 1) / (n - k - 1)
    return (θ = r.x, names = m.names, se = [v >= 0 ? sqrt(v) : NaN for v in LinearAlgebra.diag(V)],
            loglik = ll, deviance = -2ll, k = k, n = n, aicc = aicc,
            converged = maximum(abs, r.gradient) < 10gtol)
end

# --- the published model set -----------------------------------------------------------------

const MOSHER2018_DETECTION = (pA = :doy2, rA = :doy2, rBA = :elev)

"""
    mosher2018_models() -> Vector{NamedTuple}

The vital-rate structures of the 10 best models and the two reference models in Table 3 of
Mosher et al. (2018), in the order printed: each entry has `name` and `rates` (keywords for
[`Mosher2018Model`](@ref)) and `null` (the null model, in which each rate does not depend on
the other species). All share the detection structure of the paper: `pA`, `rA` quadratic in
day of year (shared slopes) and `rBA` linear in elevation.
"""
function mosher2018_models()
    e2 = (γAB = :elev2, γB = :elev2, γBA = :elev2, εAB = :elev2, εBA = :elev2)
    e1 = (γAB = :elev, γB = :elev, γBA = :elev, εAB = :elev, εBA = :elev)
    return [(name = "γAB, γB, γBA, εAB, εBA (elev²)", rates = e2, null = false),
            (name = "γAB, γB, γBA, εAB, εBA (elev)", rates = e1, null = false),
            (name = "γB, γBA, εBA (elev²)", rates = (γB = :elev2, γBA = :elev2, εBA = :elev2), null = false),
            (name = "γB, γBA, εBA (elev)", rates = (γB = :elev, γBA = :elev, εBA = :elev), null = false),
            (name = "γB, γBA (elev)", rates = (γB = :elev, γBA = :elev), null = false),
            (name = "εBA (elev)", rates = (εBA = :elev,), null = false),
            (name = "γB, γBA, εBA (perm)", rates = (γB = :perm, γBA = :perm, εBA = :perm), null = false),
            (name = "εBA (perm)", rates = (εBA = :perm,), null = false),
            (name = "γAB (amph)", rates = (γAB = :amph,), null = false),
            (name = "εBA (elev²)", rates = (εBA = :elev2,), null = false),
            (name = "constant, state-dependent", rates = NamedTuple(), null = false),
            (name = "null", rates = NamedTuple(), null = true)]
end

"""
    mosher2018_model(data, i::Integer) -> Mosher2018Model

Model `i` of [`mosher2018_models`](@ref) (the order of the paper's Table 3) for `data`, with
the paper's detection structure.
"""
function mosher2018_model(data::Mosher2018Data, i::Integer)
    s = mosher2018_models()[i]
    return Mosher2018Model(data; null = s.null, MOSHER2018_DETECTION..., s.rates...)
end

"""
    mosher2018_model_selection(data = load_mosher2018(); models = 1:12, kw...) -> Vector

Fit models `models` of [`mosher2018_models`](@ref) by [`maximum_likelihood`](@ref) and return
them sorted by AICc, each a `NamedTuple` with the fields of the fit plus `index`, `name`,
`Δaicc` and `weight` (Akaike weight within the set fitted). Keywords go to
`maximum_likelihood`. Slow (minutes); `scripts/mosher2018_refit.jl` runs it.
"""
function mosher2018_model_selection(data::Mosher2018Data = load_mosher2018();
                                    models = 1:length(mosher2018_models()), kw...)
    fits = [merge(maximum_likelihood(mosher2018_model(data, i); kw...),
                  (index = i, name = mosher2018_models()[i].name)) for i in models]
    sort!(fits; by = f -> f.aicc)
    a0 = first(fits).aicc
    w = [exp(-(f.aicc - a0) / 2) for f in fits]
    w ./= sum(w)
    return [merge(f, (Δaicc = f.aicc - a0, weight = w[k])) for (k, f) in enumerate(fits)]
end

# Reference fits for the validation targets, cached for the session (each takes seconds).
const _M18_CACHE = Dict{Any,Any}()
_mosher2018_reference_data() = get!(load_mosher2018, _M18_CACHE, :data)
_mosher2018_reference_fit(i::Integer) =
    get!(() -> maximum_likelihood(mosher2018_model(_mosher2018_reference_data(), i)), _M18_CACHE, i)
_mosher2018_constant_fit() =
    get!(() -> (m = Mosher2018Model(_mosher2018_reference_data());
                (model = m, fit = maximum_likelihood(m))), _M18_CACHE, :constant)
function _mosher2018_best_at(elevation::Real)
    m = mosher2018_model(_mosher2018_reference_data(), 1)
    return mosher2018_parameters_at(m, _mosher2018_reference_fit(1).θ; elevation)
end
# a rate of the best model over a 10-m grid of the observed elevations
function _mosher2018_best_curve(r::Symbol)
    d = _mosher2018_reference_data()
    es = minimum(d.elevation):10.0:maximum(d.elevation)
    return (elevation = collect(es), value = [getproperty(_mosher2018_best_at(e), r) for e in es])
end
