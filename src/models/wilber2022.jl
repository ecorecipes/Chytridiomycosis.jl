# Wilber MQ, Knapp RA, Smith TC, Briggs CJ (2022) Host density has limited effects on pathogen
# invasion, disease-induced declines and within-host infection dynamics across a landscape of
# disease. J Anim Ecol 91:2451–2464. https://doi.org/10.1111/1365-2656.13823
# This reproduces a published analysis;
# it is not a new model.
#
# Hypothesis 1 of the paper is a four-state hidden Markov model of the Bd status of Rana
# sierrae/muscosa populations (lakes) surveyed yearly, 2004–2019: U (Bd-free), LP (Bd below
# prevalence ρ), HP (above ρ), E (extirpated). Read from the accepted manuscript (NSF PAR
# 10401712; eqn 1, Fig. 2A, Table 1, Fig. 3) and the authors' deposit (Zenodo
# doi:10.5281/zenodo.7114634, mqwilber/host_density_and_bd v1.0.0: code/analysis_functions.R,
# code/analysisI_failed_invasion.R, code/stan_files/hmm_occupancy_fourstates_covariates_*.stan).
# Appendix S1 (the observation model in full) and S2 (missing covariates) were not obtained:
# Wiley returned HTTP 403 and the PMC copy is behind a browser check. Everything below
# follows the code; where the main text says something else, it is noted.
#
# State process (Fig. 2A; Stan transformed parameters). Rows are the state at t, columns at
# t + 1, order (U, LP, HP, E):
#     U : 1 − a,  a ω,             a (1 − ω),             0
#     LP: φ,      ω (1 − φ),       (1 − ω)(1 − φ),        0
#     HP: 0,      0,               1 − e,                 e
#     E : 0,      0,               0,                     1
# a (arrival) and e (extirpation) are constants (logit scale, Normal(0, 5) priors). φ (LP → U,
# failed invasion) and ω (stay LP rather than go to HP) are probit-linear in covariates of year
# t (eqn 1): Φ(β0 + β1 D_t + β2 Temp_t + β3 SWE_{t+1} + β4 Tad_t + β5 Temp_t SWE_{t+1}).
# SWE is snow water equivalent on 1 April of the NEXT year (code: snow[l, y + 1]), as the text
# says. Priors: intercepts Normal(0, 5), slopes Normal(0, 3), σ_D half-Cauchy(0, 2).
#
# Printed versus coded (all choices below are the code's):
# - D: eqn 1 prints log10(D); the code uses log10(D + 1), D being the larger of the VES count
#   of adults + subadults and the number swabbed, then z-scores it over the 79 lakes × 16 years.
#   Temperature and SWE are z-scored the same way. Fig. 3 effect sizes are on this scale.
# - Missing covariates are filled per lake by the previous/next observed value, or the mean of
#   the two around a gap (impute_vals). When the VES count itself is missing, φ and ω are
#   averaged over D ~ Normal(imputed value, σ_D): E[Φ(η + βD)] = Φ((η + βμ)/√(1 + β²σ_D²)),
#   and for tadpole presence a two-point mixture with weight equal to the filled value
#   (Appendix S2, not obtained). A lake with no tadpole data at all gets 0.5.
# - Every lake starts in U in 2004 with probability 1 and the 2004 observation is not used
#   (stepwise_probs[l][1] = [1, 0, 0, 0], no Pmat[l][1]). The text does not say this.
# - Observation model (Pmat, diagonal): obs U → (1, I(0, ρ), I(ρ, 1), 0) with
#   I(a, b) = ∫_a^b (1 − p)^n dp; obs LP or HP → (0, S(0, ρ), S(ρ, 1), 0) with
#   S(a, b) = ∫_a^b Pr(K ∈ class | n, p) Beta(p; y + 1, n − y + 1) dp, K ~ Binomial(n, p),
#   class LP = 1 ≤ K ≤ n ÷ ⌊1/ρ⌋ and HP = the rest; obs E → (0, 0, 0, 1). These are likelihood
#   weights, not normalised distributions (the uniform prevalence has no 1/(b − a)). The LP/HP
#   bound n ÷ ⌊1/ρ⌋ is integer division: ⌊n/3⌋ for ρ = 0.3333, while the observed state uses
#   y/n ≤ 0.3333, so a population with y/n = 1/3 is observed HP but the classes are split at
#   K ≤ n/3. Keyword `cutoff_rule = :exact` uses ⌊nρ⌋ instead.
# - The adult + subadult count is max(adult) + max(subadult) over the year's surveys, except
#   that when every subadult entry is missing the code returns "missing" even if adults were
#   counted (max_double returns maxx2 = −1). Keyword `subadult_rule = :sum` fixes this.
# - An observed E (no frogs counted, none swabbed) is kept only if it is the last non-missing
#   observation; earlier ones become missing (adjust_threes).
#
# Data: only abundance (not density = count / lake perimeter) is reproducible, because the
# lake areas (data/archival/sites_nov2020.csv) are not in the deposit. The deposit has no
# licence (Zenodo "other-open", no licence file) and is not bundled; see data/wilber2022/.

const WILBER2022_YEARS = 2004:2019
const WILBER2022_DEPOSIT = joinpath("wilber2022_si", "zenodo7114634", "extracted",
                                    "mqwilber-host_density_and_bd-4903d73")

"""
    wilber2022_path() -> String

Root of the authors' deposit for Wilber et al. (2022) (Zenodo
[10.5281/zenodo.7114634](https://doi.org/10.5281/zenodo.7114634), extracted). It is not
bundled with the package (no open licence; see `data/wilber2022/README.md`). The path is
`ENV["CHYTRIDIOMYCOSIS_WILBER2022"]` if set, else `papers/wilber2022_si/zenodo7114634/
extracted/mqwilber-host_density_and_bd-4903d73` under the first ancestor of the package
directory that has it (the lab workspace layout). Throws an `ArgumentError` if neither exists.
"""
function wilber2022_path()
    p = _wilber2022_find()
    p === nothing && throw(ArgumentError(
        "Wilber et al. 2022 deposit not found: download doi:10.5281/zenodo.7114634, extract " *
        "it and set ENV[\"CHYTRIDIOMYCOSIS_WILBER2022\"] to its root (see data/wilber2022/README.md)"))
    return p
end

function _wilber2022_find()
    e = get(ENV, "CHYTRIDIOMYCOSIS_WILBER2022", "")
    isempty(e) || return isdir(e) ? e : nothing
    d = pkgdir(@__MODULE__)
    while true
        c = joinpath(d, "papers", WILBER2022_DEPOSIT)
        isdir(c) && return c
        parent = dirname(d)
        parent == d && return nothing
        d = parent
    end
end

# --- reading the deposit (R's NA = NaN; data.table semantics reproduced) ---------------------

function _w22_csv(path)
    lines = filter(!isempty ∘ strip, readlines(path))
    return split(lines[1], ','), [split(l, ',') for l in lines[2:end]]
end

_w22num(s) = (s = strip(s); isempty(s) || s == "NA" ? NaN : parse(Float64, s))

# One group of rows per (site, year), summarised as in set_up_response_and_covariate_data.
function _w22_surveys(root; include_subadults::Bool, tadpole::Symbol, subadult_rule::Symbol)
    hdr, rows = _w22_csv(joinpath(root, "data", "formatted", "combined_Bd_data_for_analysis_allVES.csv"))
    hdr[[1, 2, 9, 11, 13, 14]] == ["site_id", "year", "bd_load", "adult", "subadult", "tadpole"] ||
        throw(ArgumentError("unexpected columns in combined_Bd_data_for_analysis_allVES.csv"))
    groups = Dict{Tuple{Int,Int},Vector{NTuple{4,Float64}}}()
    for r in rows
        y = parse(Int, r[2])
        y >= first(WILBER2022_YEARS) || continue
        push!(get!(groups, (parse(Int, r[1]), y), NTuple{4,Float64}[]),
              (_w22num(r[9]), _w22num(r[11]), _w22num(r[13]), _w22num(r[14])))
    end
    maxnm(v) = (w = filter(!isnan, v); isempty(w) ? -1.0 : maximum(w))
    out = Dict{Tuple{Int,Int},NamedTuple}()
    for (k, g) in groups
        bd = [x[1] for x in g]
        allna = all(isnan, bd)
        any(isnan, bd) && !allna &&
            throw(ArgumentError("site-year $k mixes swabbed and unswabbed rows (num_inf would be NA)"))
        ninf = allna ? 0 : count(>(0), bd)
        nsamp = allna ? 0 : length(bd)
        tads = [x[4] for x in g]
        tad = if all(isnan, tads)
            -1.0
        elseif tadpole === :presence
            Float64(any(>(0), filter(!isnan, tads)))
        else
            log10(maximum(filter(!isnan, tads)) + 1)
        end
        a1, a2 = maxnm([x[2] for x in g]), maxnm([x[3] for x in g])
        frog = if !include_subadults
            a1
        elseif a1 == -1
            a2
        elseif a2 == -1
            subadult_rule === :code ? a2 : a1     # max_double returns maxx2 (= −1) here
        else
            a1 + a2
        end
        out[k] = (num_inf = ninf, num_samps = nsamp, tadpole = tad, frog = frog,
                  frog_combined = max(frog, Float64(nsamp)))
    end
    return out
end

# impute_vals(): fill leading/trailing gaps with the nearest value, interior gaps with the
# mean of the values on either side; `missingval` marks gaps.
function _w22_impute(row::AbstractVector; missingval = -1.0)
    x = collect(Float64, row)
    obs = [!(v == missingval || isnan(v)) for v in x]
    (any(obs) && !all(obs)) || return x
    idx = findall(obs)
    for t in eachindex(x)
        obs[t] && continue
        prev = findlast(<(t), idx)
        next = findfirst(>(t), idx)
        x[t] = prev === nothing ? row[idx[next]] :
               next === nothing ? row[idx[prev]] : (row[idx[prev]] + row[idx[next]]) / 2
    end
    return x
end

# flat_scale_unflat(): z-score over all cells not equal to −1 (sample sd).
function _w22_scale(M::AbstractMatrix; missingval = -1.0)
    v = filter(!=(missingval), vec(M))
    μ = sum(v) / length(v)
    s = sqrt(sum(abs2, v .- μ) / (length(v) - 1))
    return (M .- μ) ./ s, μ, s
end

function _w22_lake_year_table(path, lakecol, valcol, lakes)
    hdr, rows = _w22_csv(path)
    il, iy, iv = findfirst(==(lakecol), hdr), findfirst(==("year"), hdr), findfirst(==(valcol), hdr)
    want = Set(lakes)
    d = Dict{Tuple{Int,Int},Float64}()
    for r in rows
        l = parse(Int, r[il])
        l in want || continue
        y = parse(Int, r[iy])
        y in WILBER2022_YEARS || continue
        haskey(d, (l, y)) && throw(ArgumentError("duplicate lake-year ($l, $y) in $path"))
        d[(l, y)] = _w22num(r[iv])
    end
    return d
end

# The 79 lakes of Hypothesis 1: rows of analysisI_failedinvasions_excluded.csv with at least
# three non-missing years (prepare_stan_data, num_to_include = 3).
function _w22_analysis_sites(root; num_include::Integer = 3)
    _, rows = _w22_csv(joinpath(root, "data", "candidate_datasets", "analysisI_failedinvasions_excluded.csv"))
    sites = Int[]
    for r in rows
        codes = parse.(Int, r[2:17])
        count(!=(-1), codes) >= num_include && push!(sites, parse(Int, r[1]))
    end
    return sort!(sites)
end

"""
    Wilber2022Data

The data of the outbreak-state HMM of Wilber et al. (2022) for one prevalence cutoff `ρ`,
built by [`wilber2022_data`](@ref) exactly as the authors' `prepare_stan_data`. One row per
lake (`sites`, sorted), one column per year 2004–2019 (`years`). Fields:

- `obs`: observed state, `-1` missing, `0` U (no Bd detected), `1` LP (prevalence `≤ ρ`), `2`
  HP, `3` E (no frogs counted or swabbed; only as the last observation);
- `sample_size`, `num_inf`: frogs swabbed and Bd-positive (`-1` for no survey);
- `abundance`: z-scored `log10(count + 1)` (adults + subadults, or the number swabbed if
  larger; gaps filled); `temperature`: z-scored mean June–September maximum air temperature
  (°C, gridMET); `snow`: z-scored 1 April snow water equivalent (mm, CDEC); `tadpole`:
  presence (0/1, 0.5 if unknown) or z-scored `log10(tadpoles + 1)`;
- `missing_abundance`: `true` where the visual-encounter count is missing;
- `likelihoods`: lake × year × state emission weights (see [`wilber2022_emission`](@ref));
- `scaling`: means and sds used to z-score `abundance` (on `log10(count + 1)`),
  `temperature` (°C) and `snow` (mm);
- `cutoff`, `cutoff_rule`, `tadpole_covariate`.
"""
struct Wilber2022Data
    sites::Vector{Int}
    years::UnitRange{Int}
    obs::Matrix{Int}
    sample_size::Matrix{Int}
    num_inf::Matrix{Int}
    abundance::Matrix{Float64}
    temperature::Matrix{Float64}
    snow::Matrix{Float64}
    tadpole::Matrix{Float64}
    missing_abundance::BitMatrix
    likelihoods::Array{Float64,3}
    scaling::NamedTuple
    cutoff::Float64
    cutoff_rule::Symbol
    tadpole_covariate::Symbol
end

"""
    wilber2022_data(; cutoff = 0.25, tadpole = :presence, include_subadults = true,
                    subadult_rule = :code, cutoff_rule = :code, path = wilber2022_path())
        -> Wilber2022Data

Build the Hypothesis 1 data set of Wilber et al. (2022) from the authors' deposit
([`wilber2022_path`](@ref)): 79 lakes, 2004–2019. `cutoff` is the prevalence ρ between the
low- and high-prevalence states (the paper uses 0.25, 0.3333 and 0.5). `tadpole = :presence`
(the full model) or `:log_abundance` (Table 1, "replace tadpoles p/a with log10 tadpoles").
`subadult_rule = :sum` and `cutoff_rule = :exact` replace two coding quirks of the original
(see [`Wilber2022Data`](@ref) and the source notes); the defaults reproduce the paper.
"""
function wilber2022_data(; cutoff::Real = 0.25, tadpole::Symbol = :presence,
                         include_subadults::Bool = true, subadult_rule::Symbol = :code,
                         cutoff_rule::Symbol = :code, path::AbstractString = wilber2022_path())
    tadpole in (:presence, :log_abundance) ||
        throw(ArgumentError("tadpole must be :presence or :log_abundance"))
    subadult_rule in (:code, :sum) || throw(ArgumentError("subadult_rule must be :code or :sum"))
    cutoff_rule in (:code, :exact) || throw(ArgumentError("cutoff_rule must be :code or :exact"))
    0 < cutoff < 1 || throw(ArgumentError("cutoff must be in (0, 1)"))
    surveys = _w22_surveys(path; include_subadults, tadpole, subadult_rule)
    alllakes = sort!(unique(first.(collect(keys(surveys)))))
    sites = _w22_analysis_sites(path)
    issubset(sites, alllakes) || throw(ArgumentError("analysis lakes missing from the survey data"))
    years = WILBER2022_YEARS
    L, Y = length(sites), length(years)
    get_(l, y, f, fill) = haskey(surveys, (l, y)) ? getfield(surveys[(l, y)], f) : fill
    ninf = [Int(get_(l, y, :num_inf, -1)) for l in sites, y in years]
    nsamp = [Int(get_(l, y, :num_samps, -1)) for l in sites, y in years]
    frog = [get_(l, y, :frog, -1.0) for l in sites, y in years]
    obs = _w22_states(ninf, nsamp, frog, cutoff)

    # covariates, filled per lake before subsetting (as the code does)
    imp(M) = reduce(vcat, permutedims(_w22_impute(r)) for r in eachrow(M))
    abund = imp([get_(l, y, :frog_combined, -1.0) for l in sites, y in years])
    tad = imp([get_(l, y, :tadpole, -1.0) for l in sites, y in years])
    tdat = _w22_lake_year_table(joinpath(path, "data", "environmental_covariates",
        "summer_lake_temperatures.csv"), "lake_id", "max_summer_temp_C", sites)
    sdat = _w22_lake_year_table(joinpath(path, "data", "environmental_covariates",
        "cdec_snowpack_april1.csv"), "lake_id", "SWEmm", sites)
    temp = imp([get(tdat, (l, y), -1.0) for l in sites, y in years])
    snow = imp([get(sdat, (l, y), -1.0) for l in sites, y in years])
    any(==(-1), abund) && throw(ArgumentError("a lake has no abundance data at all"))
    (any(==(-1), temp) || any(==(-1), snow)) &&
        throw(ArgumentError("a lake has no temperature or snow data"))
    abund_s, aμ, as = _w22_scale(log10.(abund .+ 1))
    temp_s, tμ, ts = _w22_scale(temp)
    snow_s, sμ, ss = _w22_scale(snow)
    if tadpole === :presence
        tad_s = copy(tad)
        tad_s[tad .== -1] .= 0.5
    else
        tad_s, _, _ = _w22_scale(tad)
        tad_s[tad .== -1] .= 0.0
    end
    missing_abund = BitMatrix(frog .== -1)
    lik = zeros(L, Y, 4)
    for i in 1:L, t in 1:Y
        lik[i, t, :] .= wilber2022_emission(obs[i, t], nsamp[i, t], ninf[i, t], cutoff;
                                            cutoff_rule)
    end
    scaling = (abundance = (mean = aμ, sd = as), temperature = (mean = tμ, sd = ts),
               snow = (mean = sμ, sd = ss))
    return Wilber2022Data(sites, years, obs, nsamp, ninf, abund_s, temp_s, snow_s, tad_s,
                          missing_abund, lik, scaling, Float64(cutoff), cutoff_rule, tadpole)
end

# convert_to_states() and adjust_threes()
function _w22_states(ninf, nsamp, frog, ρ)
    S = fill(typemin(Int), size(ninf))
    for i in eachindex(S)
        y, n, a = ninf[i], nsamp[i], frog[i]
        S[i] = if y == -1 || (y == 0 && n == 0 && a > 0)
            -1
        elseif y == 0 && n > 0
            0
        elseif y / n <= ρ
            1
        elseif y / n > ρ
            2
        elseif y == 0 && n == 0 && a == 0
            3
        else
            throw(ArgumentError("unclassifiable site-year (num_inf $y, swabbed $n, count $a)"))
        end
    end
    for r in eachrow(S)
        threes = findall(==(3), r)
        isempty(threes) && continue
        last_obs = findlast(!=(-1), r)
        r[threes[1:end-1]] .= -1
        threes[end] < last_obs && (r[threes[end]] = -1)
    end
    return S
end

# ∫_a^b (1 − p)^n dp (undetected_prob in the Stan code)
_w22_undetected(n, a, b) = -((1 - b)^(n + 1) - (1 - a)^(n + 1)) / (n + 1)

# summed_beta_cdfs_sample_dist: Σ_{k=lo}^{hi} ∫_a^b C(n,k) p^k (1−p)^(n−k) Beta(p; y+1, n−y+1) dp
function _w22_summed(lo, hi, n, a, b, y)
    s = 0.0
    lo <= hi || return 0.0
    for k in lo:hi
        lm = (_lgamma(n + 1) + _lgamma(n + 2) + _lgamma(k + y + 1) + _lgamma(2n - k - y + 1)) -
             (_lgamma(k + 1) + _lgamma(n - k + 1) + _lgamma(y + 1) + _lgamma(n - y + 1) +
              _lgamma(2n + 2))
        B = Beta(k + y + 1, 2n - k - y + 1)
        s += exp(lm) * (cdf(B, b) - cdf(B, a))
    end
    return abs(s)
end

_lgamma(x) = Distributions.SpecialFunctions.loggamma(float(x))

"""
    wilber2022_emission(obs, n, y, cutoff; cutoff_rule = :code) -> NTuple{4,Float64}

Emission weights `(U, LP, HP, E)` of one lake-year in the Wilber et al. (2022) code
(transformed data `Pmat`): `obs` is the observed state code (`-1` missing, `0` U, `1` LP, `2`
HP, `3` E), `n` frogs swabbed, `y` positive. Missing gives `(1, 1, 1, 1)`; obs U gives
`(1, I(0, ρ), I(ρ, 1), 0)` with `I(a, b) = ∫_a^b (1 − p)^n dp`; obs LP or HP gives
`(0, S(0, ρ), S(ρ, 1), 0)`, where `S(a, b)` integrates the probability that a binomial count
of `n` lands in the observed class against the `Beta(y + 1, n − y + 1)` density of
prevalence over `[a, b]`; obs E gives `(0, 0, 0, 1)`. The LP class is `1 ≤ K ≤ n ÷ ⌊1/ρ⌋`
(`cutoff_rule = :code`) or `1 ≤ K ≤ ⌊nρ⌋` (`:exact`). The weights are not normalised over
observations, as in the original.
"""
function wilber2022_emission(obs::Integer, n::Integer, y::Integer, ρ::Real;
                             cutoff_rule::Symbol = :code)
    obs == -1 && return (1.0, 1.0, 1.0, 1.0)
    obs == 3 && return (0.0, 0.0, 0.0, 1.0)
    obs == 0 && return (1.0, _w22_undetected(n, 0.0, ρ), _w22_undetected(n, ρ, 1.0), 0.0)
    m = cutoff_rule === :code ? n ÷ floor(Int, 1 / ρ) : floor(Int, n * ρ)
    lo, hi = obs == 1 ? (1, m) : (m + 1, n)
    return (0.0, _w22_summed(lo, hi, n, 0.0, ρ, y), _w22_summed(lo, hi, n, ρ, 1.0, y), 0.0)
end

# --- the state process ------------------------------------------------------------------------

const WILBER2022_MODELS = (:full, :log_tadpole, :linear_temp, :quadratic_temp)

# Covariate terms of the probit predictors of ω and φ, per model (Table 1):
# :full           eqn 1 (hmm_occupancy_fourstates_covariates_imputation.stan)
# :log_tadpole    eqn 1 with z-scored log10(tadpoles + 1) (…_logtadpole.stan)
# :linear_temp    no winter severity (…_linear_temp.stan)
# :quadratic_temp no winter severity, temperature + temperature² (…_nonlinear_temp.stan)
function _w22_terms(model::Symbol)
    model in (:full, :log_tadpole) &&
        return (:intercept, :abundance, :temperature, :snow, :tadpole, :temp_snow)
    model === :linear_temp && return (:intercept, :abundance, :temperature, :tadpole)
    model === :quadratic_temp &&
        return (:intercept, :abundance, :temperature, :temperature2, :tadpole)
    throw(ArgumentError("model must be one of $WILBER2022_MODELS"))
end

"""
    OutbreakHMMProblem(data::Wilber2022Data; model = :full)
    OutbreakHMMProblem(; cutoff = 0.25, model = :full, kwargs...)

The outbreak-state hidden Markov model of Wilber et al. (2022) (Hypothesis 1: failed and
successful Bd invasion), as a LogDensityProblems target on the unconstrained parameter
vector (names in `prob.names`):

- `logit_a`, `logit_ext`: Bd arrival (U → LP or HP) and extirpation (HP → E) per year;
- `omega_*`, `phi_*`: probit coefficients of ω (stay in LP rather than move to HP; also
  the LP share of arrivals) and φ (LP → U, a failed invasion) on the z-scored covariates
  of [`Wilber2022Data`](@ref);
- `log_sigma_abundance`: log sd of the abundance covariate where the count is missing.

`model` is `:full` (eqn 1: abundance, temperature, winter severity, tadpole presence,
temperature × winter severity), `:log_tadpole` (log10 tadpole abundance instead of presence),
`:linear_temp` (no winter severity) or `:quadratic_temp` (no winter severity; temperature and
temperature²), the four models of Table 1. Priors are those of the Stan code: intercepts and
`logit_a`, `logit_ext` Normal(0, 5), slopes Normal(0, 3), σ half-Cauchy(0, 2) (with the
log-Jacobian). The keyword form builds the data with [`wilber2022_data`](@ref) (with
`tadpole = :log_abundance` for `:log_tadpole`).

`solve(prob; chains = 4, warmup = 5000, iterations = 12_000, thin = 4, rng)` samples the
posterior by adaptive random-walk Metropolis (Haario et al. 2001; no MCMC library is a
dependency) and returns an [`OutbreakHMMFit`](@ref) with 12 000 draws by default.
"""
struct OutbreakHMMProblem{M}
    data::Wilber2022Data
    model::Symbol
    names::Vector{Symbol}
end

# the model is a type parameter so that the coefficient names are known at compile time
_w22_terms(::OutbreakHMMProblem{M}) where {M} = _w22_terms(M)

function OutbreakHMMProblem(data::Wilber2022Data; model::Symbol = :full)
    terms = _w22_terms(model)
    want = model === :log_tadpole ? :log_abundance : :presence
    data.tadpole_covariate === want ||
        throw(ArgumentError("model $model needs data built with tadpole = :$want"))
    names = vcat([:logit_a, :logit_ext], [Symbol(:omega_, t) for t in terms],
                 [Symbol(:phi_, t) for t in terms], [:log_sigma_abundance])
    return OutbreakHMMProblem{model}(data, model, names)
end

function OutbreakHMMProblem(; cutoff::Real = 0.25, model::Symbol = :full, kwargs...)
    _w22_terms(model)
    tad = model === :log_tadpole ? :log_abundance : :presence
    return OutbreakHMMProblem(wilber2022_data(; cutoff, tadpole = tad, kwargs...); model)
end

LogDensityProblems.dimension(p::OutbreakHMMProblem) = length(p.names)
LogDensityProblems.capabilities(::Type{<:OutbreakHMMProblem}) = LogDensityProblems.LogDensityOrder{0}()

"""
    Wilber2022Parameters

Natural-scale parameters of the outbreak-state HMM: `a` (yearly Bd arrival probability),
`ext` (yearly extirpation probability from HP), `omega` and `phi` (named tuples of probit
coefficients, keyed by covariate term) and `sigma` (sd of the missing abundance covariate).
Build from an unconstrained vector with [`wilber2022_parameters`](@ref).
"""
struct Wilber2022Parameters{T<:Real,N<:NamedTuple}
    a::T
    ext::T
    omega::N
    phi::N
    sigma::T
end

"""
    wilber2022_parameters(prob::OutbreakHMMProblem, θ) -> Wilber2022Parameters

Map the unconstrained vector `θ` (order `prob.names`) to natural-scale parameters.
"""
function wilber2022_parameters(prob::OutbreakHMMProblem, θ::AbstractVector)
    length(θ) == length(prob.names) ||
        throw(DimensionMismatch("θ has length $(length(θ)), expected $(length(prob.names))"))
    terms = _w22_terms(prob)
    m = length(terms)
    il(x) = 1 / (1 + exp(-x))
    ω = NamedTuple{terms}(Tuple(θ[3:2+m]))
    φ = NamedTuple{terms}(Tuple(θ[3+m:2+2m]))
    a, e, σ = promote(il(θ[1]), il(θ[2]), exp(θ[end]))
    return Wilber2022Parameters(a, e, ω, φ, σ)
end

"""
    outbreak_transition_matrix(a, ext, ω, φ) -> Matrix

The 4×4 yearly transition matrix of Wilber et al. (2022) Fig. 2A, states (U, LP, HP, E), rows
the current state: U → (1 − a, aω, a(1 − ω), 0); LP → (φ, ω(1 − φ), (1 − ω)(1 − φ), 0);
HP → (0, 0, 1 − ext, ext); E absorbing.
"""
function outbreak_transition_matrix(a, e, ω, φ)
    z = zero(promote_type(typeof(a), typeof(e), typeof(ω), typeof(φ)))
    return [1-a      a*ω          a*(1-ω)            z
            φ        ω*(1-φ)      (1-ω)*(1-φ)        z
            z        z            1-e                e
            z        z            z                  one(z)]
end

_Φ(x) = cdf(Normal(), x)

# ω or φ of lake i for the transition from year t to t + 1 (the code's omega[l, y], phi[l, y])
function _w22_rate(d::Wilber2022Data, model::Symbol, β::NamedTuple, σ, i, t)
    T = d.temperature[i, t]
    η = β.intercept + β.temperature * T
    if model in (:full, :log_tadpole)
        S = d.snow[i, t+1]
        η += β.snow * S + β.temp_snow * T * S
    elseif model === :quadratic_temp
        η += β.temperature2 * T^2
    end
    D, tad = d.abundance[i, t], d.tadpole[i, t]
    if !d.missing_abundance[i, t]
        return _Φ(η + β.abundance * D + β.tadpole * tad)
    end
    s = sqrt(1 + (β.abundance * σ)^2)
    model === :log_tadpole && return _Φ((η + β.tadpole * tad + β.abundance * D) / s)
    return tad * _Φ((η + β.tadpole + β.abundance * D) / s) +
           (1 - tad) * _Φ((η + β.abundance * D) / s)
end

"""
    outbreak_rates(prob::OutbreakHMMProblem, θ) -> (; omega, phi)

The lake × year matrices of ω and φ (`Y − 1` columns: column `t` drives the transition from
year `t` to `t + 1`) at parameters `θ`.
"""
function outbreak_rates(prob::OutbreakHMMProblem, θ::AbstractVector)
    p = wilber2022_parameters(prob, θ)
    d = prob.data
    L, Y = size(d.obs)
    ω = [_w22_rate(d, prob.model, p.omega, p.sigma, i, t) for i in 1:L, t in 1:Y-1]
    φ = [_w22_rate(d, prob.model, p.phi, p.sigma, i, t) for i in 1:L, t in 1:Y-1]
    return (omega = ω, phi = φ)
end

"""
    outbreak_sequences(prob::OutbreakHMMProblem, θ) -> Vector

One `(initial, transitions, likelihoods)` named tuple per lake, in the form of the generic HMM
functions ([`hmm_forward`](@ref) with `emit_first = false`, [`HMMLogDensity`](@ref)): every
lake starts in U and the 2004 observation is not used, as in the authors' code.
"""
function outbreak_sequences(prob::OutbreakHMMProblem, θ::AbstractVector)
    p = wilber2022_parameters(prob, θ)
    r = outbreak_rates(prob, θ)
    d = prob.data
    init = [1.0, 0.0, 0.0, 0.0]
    return [(initial = init,
             transitions = [outbreak_transition_matrix(p.a, p.ext, r.omega[i, t], r.phi[i, t])
                            for t in axes(r.omega, 2)],
             likelihoods = d.likelihoods[i, :, :]) for i in axes(d.obs, 1)]
end

function _w22_logprior(prob::OutbreakHMMProblem, θ)
    terms = _w22_terms(prob)
    m = length(terms)
    lp = logpdf(Normal(0, 5), θ[1]) + logpdf(Normal(0, 5), θ[2])
    for (j, t) in enumerate(terms)
        s = t === :intercept ? 5 : 3
        lp += logpdf(Normal(0, s), θ[2+j]) + logpdf(Normal(0, s), θ[2+m+j])
    end
    σ = exp(θ[end])
    return lp + log(2) + logpdf(Cauchy(0, 2), σ) + θ[end]
end

# Scaled forward recursion specialised to the 4-state chain (equal to hmm_loglikelihood on
# outbreak_sequences, without the allocations).
function _w22_loglik(prob::OutbreakHMMProblem, p::Wilber2022Parameters, i)
    d = prob.data
    Y = size(d.obs, 2)
    a, e = p.a, p.ext
    p1, p2, p3, p4 = one(a), zero(a), zero(a), zero(a)
    ll = zero(a)
    for t in 2:Y
        ω = _w22_rate(d, prob.model, p.omega, p.sigma, i, t - 1)
        φ = _w22_rate(d, prob.model, p.phi, p.sigma, i, t - 1)
        q1 = (p1 * (1 - a) + p2 * φ) * d.likelihoods[i, t, 1]
        q2 = (p1 * a * ω + p2 * ω * (1 - φ)) * d.likelihoods[i, t, 2]
        q3 = (p1 * a * (1 - ω) + p2 * (1 - ω) * (1 - φ) + p3 * (1 - e)) * d.likelihoods[i, t, 3]
        q4 = (p3 * e + p4) * d.likelihoods[i, t, 4]
        c = q1 + q2 + q3 + q4
        c > 0 || return oftype(ll, -Inf)
        ll += log(c)
        p1, p2, p3, p4 = q1 / c, q2 / c, q3 / c, q4 / c
    end
    return ll
end

function pointwise_loglikelihood(prob::OutbreakHMMProblem, θ::AbstractVector)
    p = wilber2022_parameters(prob, θ)
    return [_w22_loglik(prob, p, i) for i in axes(prob.data.obs, 1)]
end

function LogDensityProblems.logdensity(prob::OutbreakHMMProblem, θ::AbstractVector)
    all(isfinite, θ) || return -Inf
    lp = _w22_logprior(prob, θ)
    p = wilber2022_parameters(prob, θ)
    for i in axes(prob.data.obs, 1)
        lp += _w22_loglik(prob, p, i)
        isfinite(lp) || return -Inf
    end
    return lp
end

"""
    outbreak_state_probabilities(prob::OutbreakHMMProblem, θ) -> Array{Float64,3}

Smoothed probabilities of the hidden states (U, LP, HP, E) of every lake and year
(lake × year × state) at parameters `θ`, by forward–backward ([`hmm_posterior`](@ref)).
"""
function outbreak_state_probabilities(prob::OutbreakHMMProblem, θ::AbstractVector)
    seqs = outbreak_sequences(prob, θ)
    L, Y = size(prob.data.obs)
    out = zeros(L, Y, 4)
    for (i, s) in enumerate(seqs)
        out[i, :, :] .= hmm_posterior(s.initial, s.transitions, s.likelihoods; emit_first = false)
    end
    return out
end

"""
    outbreak_viterbi(prob::OutbreakHMMProblem, θ) -> Matrix{Int}

The most probable hidden state path of every lake (lake × year; codes `0` U, `1` LP, `2` HP,
`3` E) at parameters `θ` ([`hmm_viterbi`](@ref)).
"""
function outbreak_viterbi(prob::OutbreakHMMProblem, θ::AbstractVector)
    seqs = outbreak_sequences(prob, θ)
    return permutedims(reduce(hcat, [hmm_viterbi(s.initial, s.transitions, s.likelihoods;
                                                 emit_first = false).path .- 1 for s in seqs]))
end
