# Load observation layer, part 3: Taylor's power law, Var = a · mean^b, across groups of
# swabs, as Schrock et al. (2025, Proc R Soc B 292:20242013; preprint 2024) fitted it:
# ordinary least squares of log10(variance) on log10(mean), each group's mean and sample
# variance (n − 1) taken on the natural scale over all swabs INCLUDING the zeros, groups with
# at least 10 swabs and at least 3 positives (their code 01_Aggregation_Groups.R and
# 03_Fig.1_TPL_Plot.R, Zenodo doi:10.5281/zenodo.14291806, MIT). Their group summaries are on
# Dryad (doi:10.5061/dryad.mgqnk998g, CC0) but are not bundled: see data/schrock2024/README.md.

"""
    TaylorPowerLawProblem(means, variances; base = 10)
    TaylorPowerLawProblem(groups; min_swabs = 10, min_infected = 3, base = 10)

Taylor's power law `Var = a · mean^b` across groups, fitted by ordinary least squares of
`log(variance)` on `log(mean)` (logs to `base`).

The second form takes one vector of swab loads per group (zeros for uninfected swabs) and
summarises each group as Schrock et al. (2025) did: mean and sample variance over all swabs
including the zeros, keeping groups with at least `min_swabs` swabs and at least
`min_infected` positive swabs. Their thresholds are in the compendium
(`parameter_set(:schrock2024)`).

`solve(prob; level = 0.95)` returns a `NamedTuple`: `slope` (`b`), `intercept` (`log a`),
`a`, their standard errors `se_slope`, `se_intercept`, the `level` confidence interval of the
slope `ci_slope` (t quantiles, as R's `confint`), `n` (groups), `r2` and `sigma` (residual
standard deviation).

A lognormal load with the same log-scale spread in every group, and every swab infected,
gives `b = 2` exactly. Groups that differ in prevalence do not: with the zeros included,
`Var/mean² = (e^{σ²} − π)/π` grows as the infected fraction `π` falls ([`load_moments`](@ref)).
"""
struct TaylorPowerLawProblem
    means::Vector{Float64}
    variances::Vector{Float64}
    base::Float64
    function TaylorPowerLawProblem(means::AbstractVector{<:Real}, variances::AbstractVector{<:Real};
            base::Real = 10)
        length(means) == length(variances) || throw(DimensionMismatch("means and variances differ in length"))
        length(means) >= 3 || throw(ArgumentError("need at least three groups"))
        all(>(0), means) && all(>(0), variances) ||
            throw(ArgumentError("means and variances must be positive"))
        base > 1 || throw(ArgumentError("base must exceed 1"))
        return new(Float64.(means), Float64.(variances), Float64(base))
    end
end

function TaylorPowerLawProblem(groups::AbstractVector{<:AbstractVector{<:Real}};
        min_swabs::Integer = 10, min_infected::Integer = 3, base::Real = 10)
    ms, vs = Float64[], Float64[]
    for g in groups
        (length(g) >= min_swabs && count(>(0), g) >= min_infected) || continue
        m = sum(g) / length(g)
        push!(ms, m)
        push!(vs, sum(abs2, g .- m) / (length(g) - 1))
    end
    return TaylorPowerLawProblem(ms, vs; base)
end

function CommonSolve.solve(prob::TaylorPowerLawProblem; level::Real = 0.95)
    0 < level < 1 || throw(ArgumentError("level must lie in (0, 1)"))
    x = log.(prob.means) ./ log(prob.base)
    y = log.(prob.variances) ./ log(prob.base)
    n = length(x)
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum(abs2, x .- mx)
    sxx > 0 || throw(ArgumentError("all group means are equal"))
    b = sum((x .- mx) .* (y .- my)) / sxx
    c = my - b * mx
    res = y .- (c .+ b .* x)
    rss = sum(abs2, res)
    s2 = rss / (n - 2)
    se_b = sqrt(s2 / sxx)
    se_c = sqrt(s2 * (1 / n + mx^2 / sxx))
    t = quantile(TDist(n - 2), (1 + level) / 2)
    return (; slope = b, intercept = c, a = prob.base^c, se_slope = se_b, se_intercept = se_c,
        ci_slope = (b - t * se_b, b + t * se_b), n, r2 = 1 - rss / sum(abs2, y .- my),
        sigma = sqrt(s2), level)
end

"""
    load_moments(prevalence, d::UnivariateDistribution) -> (; mean, variance)

Mean and variance of the load of a randomly chosen swab when a fraction `prevalence` of
hosts is infected with loads distributed as `d` and the rest read zero:
`mean = π m`, `variance = π(v + m²) − π² m²`, with `m` and `v` the mean and variance of `d`.
These are the quantities Taylor's power law relates across groups (Schrock et al. 2025
include the zeros). For a lognormal `d` with log-scale sd `σ`, `variance/mean² =
(e^{σ²} − π)/π`.
"""
function load_moments(prevalence::Real, d::UnivariateDistribution)
    0 <= prevalence <= 1 || throw(ArgumentError("prevalence must lie in [0, 1]"))
    m, v = mean(d), var(d)
    return (; mean = prevalence * m, variance = prevalence * (v + m^2) - prevalence^2 * m^2)
end

# Schrock et al. group summaries -------------------------------------------------------------

"""
    SchrockGroup

One group of swabs (species × life stage × site × season × year) from the Dryad summary file
of Schrock et al. (2025): `dataset` (`:brazil`, `:east_bay`, `:sierra_nevada` or `:serdp`,
the eastern US), `species`, `life_stage`, `site`, `N` (swabs), `num_infected`, `mean` and
`variance` of the load over all swabs including zeros (natural scale; `NaN` when missing),
and the epizoological `phase` of Sierra Nevada groups (`""` otherwise).
"""
struct SchrockGroup
    dataset::Symbol
    species::String
    life_stage::String
    site::String
    N::Int
    num_infected::Int
    mean::Float64
    variance::Float64
    phase::String
end

"""
    schrock2024_path()

Where [`load_schrock2024`](@ref) looks for the group summary file of Schrock et al. (2025),
`data/schrock2024/analysis_aggregation_dataset.csv`. The file is not bundled: download it
from Dryad (doi:10.5061/dryad.mgqnk998g, CC0) as `data/schrock2024/README.md` describes.
"""
schrock2024_path() = joinpath(pkgdir(@__MODULE__), "data", "schrock2024", "analysis_aggregation_dataset.csv")

# Split one CSV line, honouring double quotes.
function _csv_fields(line::AbstractString)
    out, buf, q = String[], IOBuffer(), false
    for c in line
        if c == '"'
            q = !q
        elseif c == ',' && !q
            push!(out, String(take!(buf)))
        else
            print(buf, c)
        end
    end
    push!(out, String(take!(buf)))
    return out
end

_schrock_num(s) = (v = tryparse(Float64, strip(s)); v === nothing ? NaN : v)

"""
    load_schrock2024(path = schrock2024_path()) -> Vector{SchrockGroup}

Read the group summaries of Schrock et al. (2025) (`analysis_aggregation_dataset.csv` from
Dryad doi:10.5061/dryad.mgqnk998g). Throws an `ArgumentError` saying where to get the file if
it is missing.
"""
function load_schrock2024(path::AbstractString = schrock2024_path())
    isfile(path) || throw(ArgumentError("no file at $path: download analysis_aggregation_dataset.csv " *
        "from Dryad (doi:10.5061/dryad.mgqnk998g, CC0); see data/schrock2024/README.md"))
    lines = filter(!isempty ∘ strip, readlines(path))
    header = _csv_fields(lines[1])
    function col(name)
        i = findfirst(==(name), header)
        i === nothing && throw(ArgumentError("column $name missing from $path"))
        return i
    end
    ids, isp, ils, isi, iN, ii, im, iv, iph = col("dataset"), col("species"),
        col("life_stage"), col("site"), col("N"), col("num_infected"),
        col("mean_natural_with_zeros"), col("var_natural_with_zeros"), col("phase")
    groups = SchrockGroup[]
    for line in lines[2:end]
        r = _csv_fields(line)
        ph = r[iph] in ("NA", "Not Available", "Not Applicable") ? "" : r[iph]
        push!(groups, SchrockGroup(Symbol(r[ids]), r[isp], r[ils], r[isi],
            round(Int, _schrock_num(r[iN])), round(Int, _schrock_num(r[ii])),
            _schrock_num(r[im]), _schrock_num(r[iv]), ph))
    end
    return groups
end

"""
    schrock_taylor_slopes(groups; by = :dataset, min_infected = 3, level = 0.95)
        -> Dict{Symbol,NamedTuple}

Taylor's power law fitted separately to each `dataset` (or `life_stage` pooled over the
datasets, dropping unknown stages, or `phase` of the Sierra groups) of the Schrock et al.
(2025) group summaries, as their Fig. 1 and life-stage code (`03_Fig.1_TPL_Plot.R`,
`04_TPL_Life_Stage.R`) does: groups
with at least `min_infected` positive swabs, OLS of log10 variance on log10 mean (zeros
included). Each value is the `solve` result of a [`TaylorPowerLawProblem`](@ref). Groups with
a missing or zero mean or variance are dropped.
"""
function schrock_taylor_slopes(groups::AbstractVector{SchrockGroup}; by::Symbol = :dataset,
        min_infected::Integer = 3, level::Real = 0.95)
    by in (:dataset, :life_stage, :phase) || throw(ArgumentError("by must be :dataset, :life_stage or :phase"))
    keep = filter(g -> g.num_infected >= min_infected && g.mean > 0 && g.variance > 0, groups)
    by === :phase && filter!(g -> !isempty(g.phase), keep)
    by === :life_stage && filter!(g -> !(g.life_stage in ("", "unknown")), keep)
    keyof(g) = by === :dataset ? g.dataset : Symbol(getfield(g, by))
    out = Dict{Symbol,Any}()
    for k in unique(keyof.(keep))
        gs = filter(g -> keyof(g) === k, keep)
        length(gs) >= 3 || continue
        out[k] = solve(TaylorPowerLawProblem([g.mean for g in gs], [g.variance for g in gs]); level)
    end
    return out
end

# Validation helper: the TPL slope across lognormal groups with a common log-scale spread
# (σ = 2.0 natural log, every swab infected), from their exact moments. Schrock et al.: "given
# a lognormal distribution, we would expect a TPL slope of 2".
function _lognormal_taylor_slope(; σ = 2.0, meanlogs = range(0, 12; length = 25), prevalence = 1.0)
    mv = [load_moments(prevalence, LogNormal(μ, σ)) for μ in meanlogs]
    return solve(TaylorPowerLawProblem(getfield.(mv, :mean), getfield.(mv, :variance))).slope
end
