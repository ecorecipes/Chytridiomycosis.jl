# Maximum likelihood fit of the Mosher et al. (2018) dynamic two-species occupancy model with
# covariates (model and choices: src/models/mosher2018.jl). Program MARK's "dynamic
# conditional-binomial 2-species occupancy" likelihood is reproduced with the generic layer
# in src/observation/occupancy.jl; the maximisation is BFGS (finite-difference gradients)
# on the logit scale, with Wald standard errors from a finite-difference Hessian.
#
# The Dryad deposit (doi:10.5061/dryad.1c990, CC0) could not be downloaded when this was
# written (the Dryad file API needs a login), so there is no CSV reader yet: its column
# layout has not been seen. mosher2018_path() locates the two files once they are saved,
# and Mosher2018Data takes the arrays.

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
`boreal_toad_bd_covariates.csv` (doi:10.5061/dryad.1c990, CC0): `ENV["CHYTRIDIOMYCOSIS_MOSHER2018"]`
if set, otherwise `papers/mosher2018/dryad` in the package directory or one of its parents.
Throws if not found. No reader for the files is provided yet (their layout has not been
checked).
"""
function mosher2018_path()
    p = _mosher2018_find()
    p === nothing && throw(ArgumentError(
        "Mosher et al. 2018 data not found: download doi:10.5061/dryad.1c990 into " *
        "lab/papers/mosher2018/dryad or set ENV[\"CHYTRIDIOMYCOSIS_MOSHER2018\"]"))
    return p
end

const MOSHER2018_FILES = ("boreal_toad_bd_histories.csv", "boreal_toad_bd_covariates.csv")

function _mosher2018_find()
    ok(d) = all(f -> isfile(joinpath(d, f)), MOSHER2018_FILES)
    e = get(ENV, "CHYTRIDIOMYCOSIS_MOSHER2018", "")
    isempty(e) || return isdir(e) && ok(e) ? e : nothing
    d = pkgdir(@__MODULE__)
    while true
        c = joinpath(d, "papers", "mosher2018", "dryad")
        isdir(c) && ok(c) && return c
        parent = dirname(d)
        parent == d && return nothing
        d = parent
    end
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
`(εB, εBA)`, `(pA, rA)` that have the same non-constant structure share their slopes (our
reading of the parameter counts in the source; see the docs). `null = true` sets each
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
    z(x) = (s = std_or_one(x); (x .- sum(x) / length(x)) ./ s)
    std_or_one(x) = (m = sum(x) / length(x); s = sqrt(sum(abs2, x .- m) / max(length(x) - 1, 1));
                     s > 0 ? s : 1.0)
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
boundary), `loglik`, `deviance` (`−2 loglik`), `k` (number of coefficients), `n` (sites),
`aicc` and `converged`.
"""
function maximum_likelihood(m::Mosher2018Model; x0 = nothing, iters::Integer = 2000,
                            gtol::Real = 1e-5)
    k = length(m.names)
    x = x0 === nothing ? zeros(k) : collect(float.(x0))
    nll(θ) = -mosher2018_loglikelihood(m, θ)
    r = _ms_bfgs(nll, x; iters, gtol)
    V, _ = _ms_vcov(-_ms_hessian(nll, r.x))
    n = size(m.codes, 1)
    ll = -r.f
    aicc = -2ll + 2k + 2k * (k + 1) / (n - k - 1)
    return (θ = r.x, names = m.names, se = [v >= 0 ? sqrt(v) : NaN for v in LinearAlgebra.diag(V)],
            loglik = ll, deviance = -2ll, k = k, n = n, aicc = aicc,
            converged = maximum(abs, r.gradient) < 10gtol)
end
