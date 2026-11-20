# Hardy, Muths, Lambert, Schneider, Funk & Bailey (2023) Compensatory recruitment unlikely in
# high-elevation amphibian populations challenged with disease. J Appl Ecol 60:121-131,
# doi:10.1111/1365-2664.14310. Data: Dryad doi:10.5061/dryad.fttdz08wv (CC0), bundled in
# data/hardy2023/. This reproduces a published analysis; it is not a new model.
#
# VERSION READ: the version of record (J Appl Ecol 60:121-131) and its Supporting Information
# (Appendices S1.1-S1.6), checked on 2026-11-20. The implementation was first written from
# Chapter 1 of B. M. Hardy's PhD dissertation (Colorado State University, 2023), which
# reprints the paper; its numbers are identical (only the table, figure and appendix numbers
# differ).
#
# Model: robust-design Pradel model (Pradel 1996) fitted in Program MARK to adult male boreal
# toads (Anaxyrus boreas boreas) at five Colorado sites, 1998-2020, three secondary occasions per
# year (Morgan's Gulch 1-2), p = c (Huggins closed captures). Apparent survival φ (logit), per-
# capita recruitment f (log) and capture probability p (logit) are linear in covariates;
# λ = φ + f and N = n/p* are derived. The likelihood is src/observation/pradel.jl.
#
# Structures (Appendices S1.2-S1.5; the superscript Bd means "until Bd arrival, then one value for
# all sites"):
#   p:  (.), (s) site, (t^Bd) year before Bd, (s^Bd) site before Bd, (s^Bd × t^Bd) site-year
#       before Bd; post-Bd p is a single value shared by all sites.
#   φ, f: (.), (ELEV), (Bd), (t^Bd), (ELEV + t^Bd), (ELEV + Bd), (ELEV × Bd); f also DROUGHT
#       and WINSTART, alone, + Bd, × Bd and + ELEV. ELEV = elevation (m) / 100.
#   With t^Bd the elevation term applies only before Bd (Appendix S1.2: "log(f) = β18 for years
#   after Bd was detected").
#
# Choices where the text is silent (each checked against the printed −2 log L and K where
# possible; see test/hardy2023.jl):
#   - Bd status of the interval t → t + 1 is that of year t (the Bd indicator "is equal to 0
#     for years prior to Bd detection ... and 1 for years including and after Bd detection").
#     The year levels of t^Bd are then the pre-Bd start years 1998-2015 (18 levels), which gives
#     the printed K = 20 for φ(ELEV + t^Bd) and f(ELEV + t^Bd). `bd_interval = :end` uses the
#     Bd status of year t + 1 instead.
#   - Climate covariates for f in the interval t → t + 1 are those of year t − 2, i.e. three
#     years before the recruits are first counted in year t + 1 (the data README: "covariates
#     stop in 2017 because that is the last cohort ... recruited into adult captures in 2020").
#   - In Site_1_ENV.csv and Site_2_ENV.csv the columns "Drought" and "Winter Start" hold the
#     winter-onset day of year (about 270-325) and the drought days (about 60-140)
#     respectively, the reverse of their headers and of the other three files. The loader
#     assigns them by value (winter onset is a day of year after day 250).
#   - Climate covariates are centred and scaled (all sites and years pooled); the paper does
#     not say how they were scaled, and the likelihood does not depend on it.
#   - AICc uses n = total captures, which reproduces MARK's AICc from −2 log L and K.
#   - Site names follow Table 1, matched to the files by elevation (2755, 2951, 3238, 3183,
#     3366 m) and Bd arrival year (2015, 2015, 2016, 2014, 2014).

const HARDY2023_SITES = ("Delilah Pond", "Collegiate Camp", "Denny Creek", "South Cottonwood",
                         "Morgan's Gulch")
const HARDY2023_YEARS = 1998:2020

"""
    hardy2023_path(site::Integer, kind::Symbol = :cmr)

Path of a bundled file of the Hardy et al. (2023) Dryad deposit (doi:10.5061/dryad.fttdz08wv,
CC0; see `data/hardy2023/README.md`): the capture histories (`kind = :cmr`) or the site
covariates (`:env`) of site `1…5`.
"""
function hardy2023_path(site::Integer, kind::Symbol = :cmr)
    1 <= site <= 5 || throw(ArgumentError("site must be 1…5"))
    kind in (:cmr, :env) || throw(ArgumentError("kind must be :cmr or :env"))
    return joinpath(pkgdir(@__MODULE__), "data", "hardy2023",
                    "Site_$(site)_" * uppercase(string(kind)) * ".csv")
end

"""
    Hardy2023Site

One boreal toad site of Hardy et al. (2023), from [`load_hardy2023`](@ref):

- `site` (file number 1…5), `name` (Table 1), `elevation` (m) and `bd_year`, the first year
  with Bd;
- `pit`: the PIT tag of each adult male;
- `histories`: one row per toad, 69 columns (1998–2020, three occasions per year), `1` caught,
  `0` not caught, `missing` not sampled;
- `primary`: the year index (1 = 1998) of each column;
- `env`: a named tuple of vectors `year` (1996–2020), `drought` (days without precipitation in
  the active season), `winter_start` (day of year of the first sustained killing frost; both
  `missing` after 2017) and `bd` (0/1).
"""
struct Hardy2023Site
    site::Int
    name::String
    elevation::Float64
    bd_year::Int
    pit::Vector{String}
    histories::Matrix{Union{Missing,Int}}
    primary::Vector{Int}
    env::NamedTuple
end

function _hardy2023_env(site::Integer)
    lines = filter(!isempty, strip.(readlines(hardy2023_path(site, :env))))
    rows = [split(l, ',') for l in lines[2:end]]
    num(x) = x == "N/A" ? missing : parse(Float64, x)
    year = [parse(Int, r[1]) for r in rows]
    a = [num(r[2]) for r in rows]
    b = [num(r[3]) for r in rows]
    # assign by value: winter onset is a day of year (> 250), drought a count of days (< 250)
    am = collect(skipmissing(a))
    winter_first = sum(am) / length(am) > 250
    drought = winter_first ? b : a
    winter = winter_first ? a : b
    bd = [parse(Int, r[4]) for r in rows]
    elev = parse(Float64, rows[1][5])
    return (year = year, drought = drought, winter_start = winter, bd = bd), elev
end

"""
    load_hardy2023() -> Vector{Hardy2023Site}

The capture histories of adult male boreal toads at the five sites of Hardy et al. (2023)
(1,594 toads; Dryad doi:10.5061/dryad.fttdz08wv, CC0), with elevation, Bd status and the
climate covariates. The CSV columns are checked against the Program MARK history string in
the last column of each file.
"""
function load_hardy2023()
    sites = Hardy2023Site[]
    J = 3 * length(HARDY2023_YEARS)
    for s in 1:5
        lines = filter(!isempty, strip.(readlines(hardy2023_path(s, :cmr))))
        H = Matrix{Union{Missing,Int}}(undef, length(lines) - 1, J)
        pit = String[]
        for (r, l) in enumerate(lines[2:end])
            fld = split(l, ',')
            for j in 1:J
                H[r, j] = fld[j+1] == "." ? missing : parse(Int, fld[j+1])
            end
            m = match(r"/\*\s*(\S+)\s*\*/\s*([01.]+)", fld[end])
            m === nothing && throw(ArgumentError("no MARK history in row $r of site $s"))
            mark = m.captures[2]
            all(j -> (mark[j] == '.' ? missing : mark[j] - '0') === H[r, j], 1:J) ||
                throw(ArgumentError("CSV columns and MARK history differ in row $r of site $s"))
            push!(pit, m.captures[1])
        end
        env, elev = _hardy2023_env(s)
        bd_year = env.year[findfirst(==(1), env.bd)]
        push!(sites, Hardy2023Site(s, HARDY2023_SITES[s], elev, bd_year, pit, H,
                                   repeat(1:length(HARDY2023_YEARS); inner = 3), env))
    end
    return sites
end

"""
    hardy2023_captures(data = load_hardy2023()) -> Vector{NamedTuple}

Captures (capture events, the sum of the histories) at each site before and after Bd arrival,
as in Hardy et al. (2023) Table 1: `name`, `elevation`, `bd_year`, `before`, `after`,
`toads` (individuals).
"""
function hardy2023_captures(data::AbstractVector{Hardy2023Site} = load_hardy2023())
    map(data) do s
        yr = HARDY2023_YEARS[s.primary]
        c(pred) = sum(coalesce(s.histories[i, j], 0) for i in axes(s.histories, 1)
                      for j in axes(s.histories, 2) if pred(yr[j]))
        (name = s.name, elevation = s.elevation, bd_year = s.bd_year,
         before = c(<(s.bd_year)), after = c(>=(s.bd_year)), toads = size(s.histories, 1))
    end
end

const HARDY2023_P = (:constant, :site, :t_bd, :s_bd, :s_bd_t_bd)
const HARDY2023_PHI = (:constant, :elev, :bd, :t_bd, :elev_t_bd, :elev_bd, :elev_x_bd)
const HARDY2023_F = (HARDY2023_PHI..., :drought, :winstart, :drought_bd, :drought_x_bd,
                     :winstart_bd, :winstart_x_bd, :drought_elev, :winstart_elev)

# Linear-predictor design for φ or f over the intervals of every site: returns per-site
# matrices and the column names.
function _hardy2023_interval_design(data, s::Symbol, par::Symbol, bd_interval::Symbol)
    T = length(HARDY2023_YEARS)
    starts = HARDY2023_YEARS[1:T-1]
    isbd(site, y) = (bd_interval === :start ? y : y + 1) >= site.bd_year
    # year levels of t^Bd: pre-Bd intervals of any site
    levels = sort(unique(y for site in data for y in starts if !isbd(site, y)))
    E(site) = site.elevation / 100
    # climate covariates are centred and scaled over all sites and years (an affine change
    # that leaves the likelihood unchanged but helps the optimiser)
    cscale = Dict(map((:drought, :winter_start)) do w
        v = Float64[x for site in data for x in skipmissing(getproperty(site.env, w))]
        m = sum(v) / length(v)
        w => (m, sqrt(sum(abs2, v .- m) / (length(v) - 1)))
    end)
    function clim(site, y, which)
        k = findfirst(==(y - 2), site.env.year)
        v = which === :drought ? site.env.drought[k] : site.env.winter_start[k]
        v === missing && throw(ArgumentError("no $which covariate for $(site.name) in $(y - 2)"))
        μ, σ = cscale[which]
        return (Float64(v) - μ) / σ
    end
    cols = if s === :constant
        [("int", (site, y) -> 1.0)]
    elseif s === :elev
        [("int", (site, y) -> 1.0), ("elev", (site, y) -> E(site))]
    elseif s === :bd
        [("int", (site, y) -> 1.0), ("bd", (site, y) -> float(isbd(site, y)))]
    elseif s === :elev_bd
        [("int", (site, y) -> 1.0), ("elev", (site, y) -> E(site)),
         ("bd", (site, y) -> float(isbd(site, y)))]
    elseif s === :elev_x_bd
        [("int", (site, y) -> 1.0), ("elev", (site, y) -> E(site)),
         ("bd", (site, y) -> float(isbd(site, y))),
         ("elev_bd", (site, y) -> E(site) * isbd(site, y))]
    elseif s in (:t_bd, :elev_t_bd)
        c = Any[(string(l), (site, y) -> float(!isbd(site, y) && y == l)) for l in levels]
        s === :elev_t_bd && push!(c, ("elev", (site, y) -> isbd(site, y) ? 0.0 : E(site)))
        push!(c, ("post_bd", (site, y) -> float(isbd(site, y))))
        c
    elseif s in (:drought, :winstart, :drought_bd, :drought_x_bd, :winstart_bd,
                 :winstart_x_bd, :drought_elev, :winstart_elev)
        par === :f || throw(ArgumentError("climate structures apply to f only"))
        w = startswith(string(s), "drought") ? :drought : :winter_start
        c = Any[("int", (site, y) -> 1.0), (string(w), (site, y) -> clim(site, y, w))]
        if endswith(string(s), "_elev")
            push!(c, ("elev", (site, y) -> E(site)))
        elseif endswith(string(s), "_bd")
            push!(c, ("bd", (site, y) -> float(isbd(site, y))))
            occursin("_x_", string(s)) &&
                push!(c, (string(w, "_bd"), (site, y) -> clim(site, y, w) * isbd(site, y)))
        end
        c
    else
        throw(ArgumentError("unknown structure :$s for $par"))
    end
    X = [[fn(site, y) for y in starts, (_, fn) in cols] for site in data]
    return X, [Symbol(par, "_", nm) for (nm, _) in cols]
end

function _hardy2023_p_design(data, s::Symbol)
    yr(site) = HARDY2023_YEARS[site.primary]
    sampled(site) = [!ismissing(site.histories[1, j]) for j in axes(site.histories, 2)]
    pre(site, y) = y < site.bd_year
    cols = if s === :constant
        [("int", (site, k, y) -> 1.0)]
    elseif s === :site
        [(replace(d.name, r"[^A-Za-z]" => ""), (site, k, y) -> float(site.site == d.site))
         for d in data]
    elseif s === :t_bd
        levels = sort(unique(yr(site)[j] for site in data for j in eachindex(site.primary)
                             if pre(site, yr(site)[j])))
        [[(string(l), (site, k, y) -> float(pre(site, y) && y == l)) for l in levels];
         ("post_bd", (site, k, y) -> float(!pre(site, y)))]
    elseif s === :s_bd
        [[(replace(d.name, r"[^A-Za-z]" => ""), (site, k, y) -> float(pre(site, y) && site.site == d.site))
          for d in data]; ("post_bd", (site, k, y) -> float(!pre(site, y)))]
    elseif s === :s_bd_t_bd
        cells = [(d.site, y) for d in data for y in sort(unique(yr(d)[sampled(d)])) if pre(d, y)]
        [[(string(replace(data[findfirst(x -> x.site == c[1], data)].name, r"[^A-Za-z]" => ""),
                  "_", c[2]),
           (site, k, y) -> float(pre(site, y) && site.site == c[1] && y == c[2])) for c in cells];
         ("post_bd", (site, k, y) -> float(!pre(site, y)))]
    else
        throw(ArgumentError("unknown structure :$s for p"))
    end
    X = [[fn(site, k, yr(site)[k]) for k in eachindex(site.primary), (_, fn) in cols]
         for site in data]
    return X, [Symbol("p_", nm) for (nm, _) in cols]
end

"""
    hardy2023_problem(; φ = :elev_t_bd, f = :elev_t_bd, p = :s_bd_t_bd,
                      data = load_hardy2023(), bd_interval = :start) -> PradelProblem

The robust-design Pradel model of Hardy et al. (2023) for the five boreal toad sites (one
group per site) as a [`PradelProblem`](@ref). The defaults are the paper's top model,
`p(s^Bd × t^Bd) φ(ELEV + t^Bd) f(ELEV + t^Bd)` (Table 3).

Structures (Appendices S1.3–S1.5):

- `p`: `:constant`, `:site`, `:t_bd` (year before Bd), `:s_bd` (site before Bd), `:s_bd_t_bd`
  (site × year before Bd). With `_bd` the post-Bd value is one parameter for all sites.
- `φ`: `:constant`, `:elev`, `:bd`, `:t_bd`, `:elev_t_bd`, `:elev_bd`, `:elev_x_bd`.
- `f`: those of `φ` and `:drought`, `:winstart`, `:drought_bd`, `:drought_x_bd`,
  `:winstart_bd`, `:winstart_x_bd`, `:drought_elev`, `:winstart_elev` (`_bd` additive, `_x_bd`
  interaction; covariates of the year three years before the recruits are counted).

`bd_interval` sets which year's Bd status an interval `t → t + 1` takes (`:start`, year `t`, the
default, or `:end`). Coefficient names are in `prob.names`.
"""
function hardy2023_problem(; φ::Symbol = :elev_t_bd, f::Symbol = :elev_t_bd,
                           p::Symbol = :s_bd_t_bd,
                           data::AbstractVector{Hardy2023Site} = load_hardy2023(),
                           bd_interval::Symbol = :start)
    φ in HARDY2023_PHI || throw(ArgumentError("φ structure must be one of $HARDY2023_PHI"))
    f in HARDY2023_F || throw(ArgumentError("f structure must be one of $HARDY2023_F"))
    p in HARDY2023_P || throw(ArgumentError("p structure must be one of $HARDY2023_P"))
    bd_interval in (:start, :end) || throw(ArgumentError("bd_interval must be :start or :end"))
    Xφ, nφ = _hardy2023_interval_design(data, φ, :φ, bd_interval)
    Xf, nf = _hardy2023_interval_design(data, f, :f, bd_interval)
    Xp, np = _hardy2023_p_design(data, p)
    groups = [PradelGroup(s.histories; primary = s.primary) for s in data]
    designs = [(φ = Xφ[i], f = Xf[i], p = Xp[i]) for i in eachindex(data)]
    return PradelProblem(groups, designs; names = vcat(nφ, nf, np))
end

"""
    hardy2023_fit(; kwargs...) -> PradelFit

Maximum-likelihood fit of [`hardy2023_problem`](@ref)`(; kwargs...)` (BFGS with
finite-difference gradients; a few seconds for the top model).
"""
hardy2023_fit(; kwargs...) = solve(hardy2023_problem(; kwargs...))

"""
    hardy2023_model_set(parameter::Symbol; data = load_hardy2023()) -> Vector{NamedTuple}

Fit one focal model set of Hardy et al. (2023) and rank it by AICc: `parameter` is `:p`
(Appendix S1.3), `:φ` (S1.4) or `:f` (S1.5). Each model varies the focal parameter and keeps the
others constant. Entries (sorted by AICc) have `structure`, `deviance` (−2 log L), `k`,
`aicc`, `Δaicc` and `weight`.
"""
function hardy2023_model_set(parameter::Symbol;
                             data::AbstractVector{Hardy2023Site} = load_hardy2023())
    structs = parameter === :p ? HARDY2023_P : parameter === :φ ? HARDY2023_PHI :
              parameter === :f ? HARDY2023_F :
              throw(ArgumentError("parameter must be :p, :φ or :f"))
    fits = map(structs) do s
        kw = (φ = :constant, f = :constant, p = :constant)
        kw = merge(kw, NamedTuple{(parameter,)}((s,)))
        r = solve(hardy2023_problem(; kw..., data); vcov = false)
        (structure = s, deviance = r.deviance, k = r.k, aicc = r.aicc)
    end
    sort!(collect(fits); by = x -> x.aicc) |> v -> begin
        a0 = first(v).aicc
        w = [exp(-(x.aicc - a0) / 2) for x in v]
        w ./= sum(w)
        [merge(x, (Δaicc = x.aicc - a0, weight = w[i])) for (i, x) in enumerate(v)]
    end
end

"""
    pradel_coefficient(fit::PradelFit, name::Symbol) -> (value, se)

A coefficient of a [`PradelFit`](@ref) on the link scale, with its Wald standard error.
"""
function pradel_coefficient(fit::PradelFit, name::Symbol)
    i = findfirst(==(name), fit.problem.names)
    i === nothing && throw(ArgumentError("no coefficient :$name"))
    return (value = fit.θ[i], se = sqrt(fit.vcov[i, i]))
end

"""
    hardy2023_summary(fit::PradelFit; data = load_hardy2023()) -> NamedTuple

The quantities that Hardy et al. (2023) report from a fit of [`hardy2023_problem`](@ref):

- `φ_post`, `f_post`, `λ_post`: the post-Bd apparent survival, per-capita recruitment and
  realized growth rate (shared by all sites), each `(value, se)`; `λ_post` also has a 95%
  interval `ci` on the log scale;
- `λ_pre`: the realized growth rates `(site, year, value, se)` of the intervals `year → year + 1`
  that end before Bd arrival (Fig. 3), between sampled years of each site;
- `λ_pre_range`: the smallest and largest of those;
- `N`: the estimated number of breeding males of each site by year (Fig. 2), `NaN` in
  unsampled years;
- `p_range`: the range of the estimated capture probabilities at sampled occasions.

`φ_post`, `f_post` and `λ_post` are taken from the post-Bd interval of the site with the
earliest arrival.
"""
function hardy2023_summary(fit::PradelFit; data::AbstractVector{Hardy2023Site} = load_hardy2023())
    prob = fit.problem
    yrs = collect(HARDY2023_YEARS)
    s0 = argmin([s.bd_year for s in data])
    tpost = findfirst(==(data[s0].bd_year), yrs)
    est(q) = pradel_se(fit, θ -> getproperty(pradel_derived(prob, θ, s0), q)[tpost])
    φp, fp, λp = est(:φ), est(:f), est(:λ)
    ci = λp.value .* exp.((-1.96, 1.96) .* (λp.se / λp.value))
    λpre = NamedTuple[]
    for (g, s) in enumerate(data)
        sampled = [any(j -> s.primary[j] == t && !ismissing(s.histories[1, j]),
                       eachindex(s.primary)) for t in eachindex(yrs)]
        for t in 1:length(yrs)-1
            (yrs[t+1] < s.bd_year && sampled[t] && sampled[t+1]) || continue
            r = pradel_se(fit, θ -> pradel_derived(prob, θ, g).λ[t])
            push!(λpre, (site = s.name, year = yrs[t], value = r.value, se = r.se))
        end
    end
    vals = [x.value for x in λpre]
    lo, hi = λpre[argmin(vals)], λpre[argmax(vals)]
    ps = [fit.groups[g].p[j] for g in eachindex(data) for j in eachindex(data[g].primary)
          if prob.groups[g].sampled[j]]
    return (φ_post = φp, f_post = fp, λ_post = merge(λp, (ci = ci,)), λ_pre = λpre,
            λ_pre_range = (lo, hi),
            N = [(site = s.name, year = yrs, N = fit.groups[g].N) for (g, s) in enumerate(data)],
            p_range = extrema(ps))
end

# One fit per structure per session, for the validation registry
const _HARDY2023_REFERENCE = Dict{Any,PradelFit}()
_hardy2023_reference_fit(; φ = :elev_t_bd, f = :elev_t_bd, p = :s_bd_t_bd) =
    get!(() -> hardy2023_fit(; φ, f, p), _HARDY2023_REFERENCE, (φ, f, p))
