# Reproduction of the analyses of Malagon et al. (2020) Sci Rep 10:5584
# (doi:10.1038/s41598-020-62351-x) from the authors' Supplementary File 2 (data, R script and
# MATLAB code; data/malagon2020). This reproduces a published analysis; it is not a new
# model. Wiki: lab/wiki/models/bsal-transmission-functions.md
#
# The pipeline, as the supplement and the code do it:
#   1. Experiment 1: contacts per hour in 1-m² pools at 2, 4 and 8 newts, without and with
#      plants. Each pool-day mean of the four diel observations is log(x + 1) transformed and
#      averaged; the back-transformed means are the "total contacts" of the fit (MasterRun.m).
#   2. A fourth point at 133 newts m⁻² comes from the Experiment 2 cohabitations (2 newts in
#      0.015 m²). Per capita contacts y = total/N are fitted to f(N) = cN/(K + N) by least
#      squares (NonLinFit.m, objfun.m); Table S1 reports c × 24 (per day) and K.
#   3. β and γ are chosen on a grid (β = 0:0.01:0.98, γ = 0:0.01:2) to minimise the squared
#      error between the cumulative proportion infectious C(t)/S(0) of the cohabitation model
#      and the swab data at days 6, 12, 18 and 24 (SolveSEImodel6.m). The model is evaluated
#      on linspace(0, 100, 100) and interpolated linearly, as MATLAB's interp1 does.

const MALAGON_FILES = Dict(:contacts => "contactrates.csv", :survival => "survival.csv",
    :swabs => "cohabitation_swabs.csv", :fit_inputs => "contact_fit_inputs.csv")

"""
    malagon2020_path(name::Symbol)

Path of a bundled data file of Malagon et al. (2020) (Supplementary File 2; see
`data/malagon2020/README.md`). `name` is one of:

- `:contacts`: contacts per hour in Experiment 1 (`contactrates.csv`);
- `:survival`: survival of the Experiment 2 newts (`survival.csv`);
- `:swabs`: cumulative proportions of newts qPCR-positive in Experiment 2, transcribed from
  the authors' `MasterRun.m` (`cohabitation_swabs.csv`);
- `:fit_inputs`: the total contacts per hour that the authors' code fits (`contact_fit_inputs.csv`).
"""
function malagon2020_path(name::Symbol)
    haskey(MALAGON_FILES, name) ||
        throw(ArgumentError("name must be one of $(sort!(collect(keys(MALAGON_FILES))))"))
    return joinpath(pkgdir(@__MODULE__), "data", "malagon2020", MALAGON_FILES[name])
end

function _csv_rows(path)
    lines = filter(!isempty ∘ strip, readlines(path))
    lines[1] = replace(lines[1], '﻿' => "")
    header = strip.(split(lines[1], ','))
    function col(name)
        i = findfirst(==(name), header)
        i === nothing && throw(ArgumentError("column $name missing from $path"))
        return i
    end
    return [strip.(split(l, ',')) for l in lines[2:end]], col
end

const MALAGON_DENSITY_WORDS = Dict("Two" => 2, "Four" => 4, "Eight" => 8)

"""
    load_malagon2020() -> (; contacts, swabs, fit_inputs, survival)

The data of Malagon et al. (2020):

- `contacts`: one `NamedTuple` per video hour of Experiment 1, with `density` (newts per
  1-m² pool), `habitat` (`:simple` without plants, `:complex` with plants), `day`, `diel`
  (1–4), `pool` and `contacts` (contacts in the hour); 672 rows;
- `swabs`: one row per Experiment 2 treatment and swab day, with `disease_state` (days since
  the infected newt was exposed: 12, 18, 24), `cohabitation` (minutes: 10, 30), `day` and
  `proportion` (cumulative proportion of the five susceptible newts positive on two swabs);
- `fit_inputs`: the `(habitat, density, total)` total contacts per hour that the authors'
  code fits to the Holling function;
- `survival`: one row per Experiment 2 newt with `group` (`:susceptible` for newts infected
  by contact, `:infected` for the newts exposed to zoospores), `disease_state` (12, 18, 24,
  or `nothing` for the exposed newts), `contact` (`"1second"`, `"10minutes"`,
  `"30minutes"` or `"NA"`), `days` (days survived, censored at 90) and `dead`.
"""
function load_malagon2020()
    rows, col = _csv_rows(malagon2020_path(:contacts))
    idn, ipl, ic, idy, idi, ipo = col("Density"), col("Plants"), col("Contacts/Hr"),
        col("Day"), col("Diel"), col("Pool")
    contacts = [(density = MALAGON_DENSITY_WORDS[r[idn]],
                 habitat = r[ipl] == "Yes" ? :complex : :simple,
                 day = parse(Int, r[idy]), diel = parse(Int, r[idi]), pool = parse(Int, r[ipo]),
                 contacts = parse(Float64, r[ic])) for r in rows]
    rows, col = _csv_rows(malagon2020_path(:swabs))
    a, b, c, d = col("disease_state"), col("cohabitation_min"), col("day"),
        col("cumulative_proportion")
    swabs = [(disease_state = parse(Int, r[a]), cohabitation = parse(Int, r[b]),
              day = parse(Float64, r[c]), proportion = parse(Float64, r[d])) for r in rows]
    rows, col = _csv_rows(malagon2020_path(:fit_inputs))
    h, n, t = col("habitat"), col("density"), col("total_contacts_per_hour")
    fit_inputs = [(habitat = Symbol(r[h]), density = parse(Float64, r[n]),
                   total = parse(Float64, r[t])) for r in rows]
    rows, col = _csv_rows(malagon2020_path(:survival))
    ig, it, ict, id, idd = col("Group"), col("Treatment"), col("CT"), col("Days Survival"),
        col("Dead")
    survival = [(group = r[ig] == "Susceptible" ? :susceptible : :infected,
                 disease_state = startswith(r[it], "Day") ? parse(Int, r[it][5:end]) : nothing,
                 contact = r[ict], days = parse(Float64, r[id]), dead = r[idd] == "1")
                for r in rows]
    return (contacts = contacts, swabs = swabs, fit_inputs = fit_inputs, survival = survival)
end

_mean(x) = sum(x) / length(x)

# pool-day means of the contacts per hour, by (habitat, density)
function _pool_day_means(contacts)
    groups = Dict{Tuple{Symbol,Int,Int,Int},Vector{Float64}}()
    for c in contacts
        push!(get!(groups, (c.habitat, c.density, c.day, c.pool), Float64[]), c.contacts)
    end
    out = Dict{Tuple{Symbol,Int},Vector{Float64}}()
    for (k, v) in groups
        push!(get!(out, (k[1], k[2]), Float64[]), _mean(v))
    end
    return out
end

"""
    contact_means(contacts = load_malagon2020().contacts)
        -> Vector{NamedTuple}

Mean contacts per hour in each habitat and density of Experiment 1, computed as in the
authors' R script: the four diel observations of a pool on one day are averaged, and the
pool-day means are averaged on the `log(x + 1)` scale and back-transformed
(`exp(mean) − 1`). Each row has `habitat`, `density`, `n` (pool-days), `mean` (the
back-transformed mean, the paper's Fig. 1 bars and the inputs of the Holling fit) and
`arithmetic` (the arithmetic mean).
"""
function contact_means(contacts = load_malagon2020().contacts)
    pd = _pool_day_means(contacts)
    rows = [(habitat = k[1], density = k[2], n = length(v),
             mean = exp(_mean(log1p.(v))) - 1, arithmetic = _mean(v)) for (k, v) in pd]
    return sort!(rows; by = r -> (r.habitat === :complex, r.density))
end

# one-way ANOVA F statistic and degrees of freedom
function _oneway_F(groups::Vector{Vector{Float64}})
    all = reduce(vcat, groups)
    m, k, n = _mean(all), length(groups), length(all)
    ssb = sum(length(g) * (_mean(g) - m)^2 for g in groups)
    ssw = sum(sum((x - _mean(g))^2 for x in g) for g in groups)
    return (F = (ssb / (k - 1)) / (ssw / (n - k)), df = (k - 1, n - k))
end

"""
    contact_anova(contacts = load_malagon2020().contacts) -> NamedTuple

The one-way analyses of variance of Malagon et al. (2020, Results) on `log(x + 1)` of the
pool-day mean contacts per hour, as in their R script:

- `density_simple`, `density_complex`: density effect within each habitat (paper:
  F₂,₈₁ = 199.3 without plants and 59.8 with plants);
- `habitat_8`, `habitat_4`, `habitat_2`: habitat effect at each density (paper: F₁,₅₈ = 80.3
  at 8 newts and F₁,₅₂ = 4.89 at 4 newts).

Each entry is `(F, df)`.
"""
function contact_anova(contacts = load_malagon2020().contacts)
    pd = _pool_day_means(contacts)
    l(h, n) = log1p.(pd[(h, n)])
    dens(h) = _oneway_F([l(h, 2), l(h, 4), l(h, 8)])
    hab(n) = _oneway_F([l(:simple, n), l(:complex, n)])
    return (density_simple = dens(:simple), density_complex = dens(:complex),
        habitat_8 = hab(8), habitat_4 = hab(4), habitat_2 = hab(2))
end

"""
    fit_contact_function(density, total; Kmax = 100) -> (; c_hour, c, K, ssr)

Least-squares fit of the Holling type II contact function `f(N) = c N / (K + N)` to the
per capita contacts `total ./ density` (contacts per newt per hour), as in the authors'
`NonLinFit.m` and `objfun.m`, with `0 ≤ K ≤ Kmax` (the code's bound; the paper prints the
bounds of `c` and `K` the other way round). For each `K` the best `c` is linear, so the fit
is a one-dimensional search over `K`. Returns the hourly `c_hour`, `c = 24 c_hour`
(contacts newt⁻¹ day⁻¹, as in Supplementary Table S1), `K` (newts m⁻²) and the sum of
squared residuals `ssr` (contacts² newt⁻² h⁻²).
"""
function fit_contact_function(density::AbstractVector, total::AbstractVector; Kmax::Real = 100)
    length(density) == length(total) || throw(DimensionMismatch("density and total"))
    N, y = float.(density), float.(total) ./ float.(density)
    function best(K)
        g = N ./ (K .+ N)
        c = sum(g .* y) / sum(g .^ 2)
        return c, sum((c .* g .- y) .^ 2)
    end
    # coarse grid then golden-section refinement on the bracket
    Ks = range(0, Kmax; length = 2001)
    i = argmin([best(K)[2] for K in Ks])
    lo, hi = Ks[max(i - 1, 1)], Ks[min(i + 1, length(Ks))]
    φ = (sqrt(5) - 1) / 2
    for _ in 1:200
        a, b = hi - φ * (hi - lo), lo + φ * (hi - lo)
        best(a)[2] < best(b)[2] ? (hi = b) : (lo = a)
    end
    K = (lo + hi) / 2
    c, ssr = best(K)
    return (c_hour = c, c = 24c, K = K, ssr = ssr)
end

"""
    fit_cohabitation(c, K, observed; days = [6, 12, 18, 24], β = 0:0.01:0.98,
                     γ = 0:0.01:2, N = 133, duration = 10/1440, interpolate = true)
        -> (; β, γ, ssr, grid)

Grid search of Malagon et al. (2020, Supplementary Information; `SolveSEImodel6.m`) for the
transmission probability `β` and the latency rate `γ` (day⁻¹). For each pair it computes
the cumulative proportion infectious `C(t)/S(0)` of the [`CohabitationProblem`](@ref) with
contact parameters `c` (contacts newt⁻¹ day⁻¹) and `K`, and sums the squared differences
from `observed` at `days`. The pair with the smallest sum is returned, with the matrix
`grid` of sums (rows `γ`, columns `β`); ties go to the smallest `β`, then the smallest `γ`.

With `interpolate = true` the model is sampled on `range(0, 100; length = 100)` and
interpolated linearly to `days`, as the authors' code does with `interp1`.
"""
function fit_cohabitation(c::Real, K::Real, observed::AbstractVector;
        days::AbstractVector = [6.0, 12.0, 18.0, 24.0], β = 0:0.01:0.98, γ = 0:0.01:2,
        N::Real = 133.0, duration::Real = 10 / 1440, interpolate::Bool = true)
    length(observed) == length(days) || throw(DimensionMismatch("observed and days"))
    tgrid = collect(range(0.0, 100.0; length = 100))
    grid = zeros(length(γ), length(β))
    for (n, b) in enumerate(β), (m, g) in enumerate(γ)
        prob = CohabitationProblem(MalagonParameters(; β = b, c = c, K = K, γ = g, d = 0.0);
            N = N, duration = duration)
        pred = if interpolate
            C = cohabitation_curve(prob, tgrid).infectious
            [_interp1(tgrid, C, t) for t in days]
        else
            cohabitation_curve(prob, days).infectious
        end
        grid[m, n] = sum(abs2, observed .- pred)
    end
    i = argmin(grid)         # column-major: smallest β, then smallest γ, as MATLAB's find
    return (β = β[i[2]], γ = γ[i[1]], ssr = grid[i], grid = grid)
end

function _interp1(x, y, t)
    j = clamp(searchsortedlast(x, t), 1, length(x) - 1)
    w = (t - x[j]) / (x[j + 1] - x[j])
    return (1 - w) * y[j] + w * y[j + 1]
end

"""
    malagon_table_s1(data = load_malagon2020(); kwargs...) -> Vector{NamedTuple}

Reproduce Supplementary Table S1 of Malagon et al. (2020): for each habitat (`:simple`, no
plants; `:complex`, plants) the Holling fit ([`fit_contact_function`](@ref) of
`data.fit_inputs`), and for each disease state and cohabitation time the grid-search
estimates of `β` and `γ` with their sum of squared residuals
([`fit_cohabitation`](@ref)). Each row has `habitat`, `disease_state`, `cohabitation`,
`β`, `γ`, `ssr`, `c` (per day) and `K`. Keywords are passed to `fit_cohabitation`.
"""
function malagon_table_s1(data = load_malagon2020(); kwargs...)
    rows = NamedTuple[]
    for h in (:simple, :complex)
        inp = filter(r -> r.habitat === h, data.fit_inputs)
        cf = fit_contact_function([r.density for r in inp], [r.total for r in inp])
        for ds in (12, 18, 24), ch in (10, 30)
            obs = sort!(filter(r -> r.disease_state == ds && r.cohabitation == ch, data.swabs);
                by = r -> r.day)
            f = fit_cohabitation(cf.c, cf.K, [r.proportion for r in obs];
                days = [r.day for r in obs], kwargs...)
            push!(rows, (habitat = h, disease_state = ds, cohabitation = ch, β = f.β,
                γ = f.γ, ssr = f.ssr, c = cf.c, K = cf.K))
        end
    end
    return rows
end

# Table S1 once per session, for the validation registry (data/validation/targets.toml).
const _MALAGON_REFERENCE = Ref{Any}(nothing)
function _malagon_reference_table()
    _MALAGON_REFERENCE[] === nothing && (_MALAGON_REFERENCE[] = malagon_table_s1())
    return _MALAGON_REFERENCE[]
end
_malagon_row(h, ds, ch) = only(filter(r -> r.habitat === h && r.disease_state == ds &&
    r.cohabitation == ch, _malagon_reference_table()))
