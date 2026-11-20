# Logistic fits of Sheets et al. 2021 (Front Vet Sci 8:687084, doi:10.3389/fvets.2021.687084)
# to the authors' deposited viability data (figshare doi:10.6084/m9.figshare.14714841, CC BY
# 4.0; data/sheets2021). This reproduces a published analysis; it is not a new model.
#
# The likelihood is that of the authors' R function obj(): Normal errors whose standard
# deviation b0 + 1e-4 + b1 μ grows with the mean, over the non-missing MTTn readings of one
# isolate at one temperature. The authors optimise abs() of the parameters; we optimise their
# logarithms, which spans the same set. Their start (x0 = 0.005, r = 0.3, K = 1, b0 = 0.01,
# b1 = 0.1) and their 30 random restarts (each parameter × U(0.5, 1.5)) are kept.

"""
    SheetsReading(isolate, T, day, mttn)

One well reading of Sheets et al. (2021): isolate (`"LA"`, `"TN"`, `"VT"`, `"NM"`, `"OH"`),
temperature `T` (°C), `day`, and normalised viability `mttn` (MTT minus heat-killed MTT, OD
at 570 nm; `NaN` when missing).
"""
struct SheetsReading
    isolate::String
    T::Float64
    day::Float64
    mttn::Float64
end

"""
    sheets2021_path()

Path of the bundled viability data of Sheets et al. (2021) (figshare
doi:10.6084/m9.figshare.14714841, CC BY 4.0; see `data/sheets2021/README.md`).
"""
sheets2021_path() = joinpath(pkgdir(@__MODULE__), "data", "sheets2021", "mtt_od.csv")

"""
    load_sheets2021(; path = sheets2021_path()) -> Vector{SheetsReading}

The 5509 readings of the deposited workbook (sheet `Ch1_MTT_OD`), with missing `MTTn` as
`NaN`.
"""
function load_sheets2021(; path::AbstractString = sheets2021_path())
    lines = filter(!isempty ∘ strip, readlines(path))
    hdr = split(lines[1], ',')
    col(n) = findfirst(==(n), hdr)
    id, ii, iT, im = col("Day"), col("Isolate"), col("Temp"), col("MTTn")
    num(s) = s == "NA" ? NaN : parse(Float64, s)
    return [let r = split(l, ',')
                SheetsReading(String(r[ii]), num(r[iT]), num(r[id]), num(r[im]))
            end for l in lines[2:end]]
end

"""
    SheetsGrowthModel(isolate, T; data = load_sheets2021())

The likelihood of Sheets et al.'s (2021) logistic fit for one isolate at one temperature
(°C): `MTTn_i ~ Normal(μ_i, b0 + 1e-4 + b1 μ_i)`, with `μ_i` the logistic curve
[`delayed_logistic`](@ref)`(t_i, x0, K, r)` (no delay). Parameters, all on the log scale:
`x0, r, K, b0, b1` ([`parameter_names`](@ref)). The model implements the LogDensityProblems
interface; the log density is the log-likelihood (a flat prior on the log parameters). Fit
it with [`maximum_likelihood`](@ref).
"""
struct SheetsGrowthModel
    isolate::String
    T::Float64
    t::Vector{Float64}
    y::Vector{Float64}
end

function SheetsGrowthModel(isolate::AbstractString, T::Real; data = load_sheets2021())
    rows = [d for d in data if d.isolate == isolate && d.T == T && !isnan(d.mttn)]
    isempty(rows) && throw(ArgumentError("no readings for isolate \"$isolate\" at $T °C"))
    return SheetsGrowthModel(String(isolate), Float64(T), [d.day for d in rows],
        [d.mttn for d in rows])
end

parameter_names(::SheetsGrowthModel) = (:x0, :r, :K, :b0, :b1)
LogDensityProblems.dimension(::SheetsGrowthModel) = 5
LogDensityProblems.capabilities(::Type{SheetsGrowthModel}) = LogDensityProblems.LogDensityOrder{0}()
nobs(m::SheetsGrowthModel) = length(m.y)

function LogDensityProblems.logdensity(m::SheetsGrowthModel, x::AbstractVector)
    length(x) == 5 || throw(DimensionMismatch("expected 5 parameters"))
    all(isfinite, x) || return -Inf
    x0, r, K, b0, b1 = exp.(x)
    ll = 0.0
    for i in eachindex(m.y)
        μ = delayed_logistic(m.t[i], x0, K, r)
        ll += logpdf(Normal(μ, b0 + 1e-4 + b1 * μ), m.y[i])
    end
    return isfinite(ll) ? ll : -Inf
end

# The start of the authors' fitlogistic() (an algorithm setting, not a published value)
const SHEETS_START = (x0 = 0.005, r = 0.3, K = 1.0, b0 = 0.01, b1 = 0.1)

"""
    maximum_likelihood(m::SheetsGrowthModel; starts = 30, iters = 3000, rng = Random.Xoshiro(2021))

Maximum likelihood fit of a [`SheetsGrowthModel`](@ref) as in the authors' `fitlogistic()`:
Nelder–Mead from their start and from `starts − 1` points with each parameter multiplied by
`U(0.5, 1.5)`, keeping the best. Returns a `NamedTuple` with `model`, `x` (log parameters),
`estimates` (natural scale, by name), `loglik` and `nobs`.
"""
function maximum_likelihood(m::SheetsGrowthModel; starts::Integer = 30, iters::Integer = 3000,
        rng::AbstractRNG = Random.Xoshiro(2021))
    f(x) = -LogDensityProblems.logdensity(m, x)
    p0 = collect(values(SHEETS_START))
    best, bestf = log.(p0), Inf
    for s in 1:starts
        x = log.(s == 1 ? p0 : p0 .* (0.5 .+ rand(rng, 5)))
        x, fx = nelder_mead(f, x; iters, step = 0.5)
        x, fx = nelder_mead(f, x; iters = iters ÷ 2, step = 0.05)
        fx < bestf && ((best, bestf) = (x, fx))
    end
    est = NamedTuple{parameter_names(m)}(Tuple(exp.(best)))
    return (model = m, x = best, estimates = est, loglik = -bestf, nobs = nobs(m))
end

"""
    sheets_fits(; data = load_sheets2021(), kwargs...) -> Vector{NamedTuple}

The logistic fit of every isolate at every temperature (35 fits; Sheets et al. 2021 Table 2
and Fig. 4), each a `NamedTuple` with `isolate`, `T`, the estimates `x0, r, K, b0, b1`,
`loglik` and `nobs`. Keywords go to [`maximum_likelihood`](@ref).
"""
function sheets_fits(; data = load_sheets2021(), kwargs...)
    out = NamedTuple[]
    for iso in SHEETS_ISOLATES, T in SHEETS_TEMPERATURES
        fit = maximum_likelihood(SheetsGrowthModel(iso, T; data); kwargs...)
        push!(out, merge((isolate = iso, T = T), fit.estimates,
            (loglik = fit.loglik, nobs = fit.nobs)))
    end
    return out
end

# One fit per isolate and temperature per session, for the validation registry
const _SHEETS_REFERENCE = Dict{Tuple{String,Float64},Any}()
_sheets_reference_fit(isolate::AbstractString, T::Real) =
    get!(() -> maximum_likelihood(SheetsGrowthModel(isolate, T)), _SHEETS_REFERENCE,
        (String(isolate), Float64(T)))
