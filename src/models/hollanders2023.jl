# Hollanders M, Grogan LF, Nock CJ, McCallum HI, Newell DA (2023) Recovered frog populations
# coexist with endemic Batrachochytrium dendrobatidis despite load-dependent mortality.
# Ecological Applications 33:e2724. https://doi.org/10.1002/eap.2724
# This reproduces a published analysis; it is
# not a new model.
#
# Read from the paper (CC BY), Appendix S1 and the authors' code: Zenodo
# doi:10.5281/zenodo.6981761 (m2-cmr.qmd, MIT licence) and the repository
# mhollanders/mfleayi-adults at commit a325a4c (appendix-s1.qmd, mcmc-samples/). Data: Dryad
# doi:10.5061/dryad.g1jwstqtb (CC0), bundled in data/hollanders2023/ (byte-identical to the
# Dryad files by the SHA-256 digests in the Dryad API; see the README there).
#
# Model (m2-cmr.qmd). Robust-design multistate Cormack–Jolly–Seber model of adult Mixophyes
# fleayi at three sites (Brindle, Tuntable, Bat Cave Creeks), 23 primary occasions with up to
# 3 secondary surveys, 2017–2021. Latent states (alive Bd−, alive Bd+, dead); the state
# process between primaries is exp(TRM τ) with the generator of infection_generator and τ the
# interval in units of six weeks (ints.csv). Hazards are log-linear: mortality φ1 (Bd−),
# φ2 = φ1 exp(β_status + β_load x + …) (Bd+), gaining ψ12 and clearing ψ21 (x enters ψ21);
# recapture is logit-linear per secondary survey, p2 = p1 with status and load terms. x is the
# individual's log10 load (ITS1 copies per swab) at the start of the interval, centred and
# divided by two standard deviations (over the 248 observed loads); missing loads are imputed
# from a normal submodel with temperature, rainfall and individual effects. Site and survey
# random effects, body condition (scaled mass index), RJMCMC over all slopes, NIMBLE, 4 × 25 000
# iterations.
#
# Printed versus coded (all the code's):
# - Appendix S1 prints the observation matrix without δ; the code multiplies the Bd+ row by
#   δ_t, which is 1 except at the last primary occasion, when frogs were not swabbed: a Bd+
#   frog seen then is recorded as "seen, Bd−" (keyword `swab_final`).
# - The text says "all pairwise interactions" on mortality; the code has condition × temperature,
#   condition × load and temperature × load only.
# - "Average 6-weekly rainfall over the primary interval" is coded as the centred, scaled
#   interval rainfall divided by the interval length (scl(rain) / τ), not the scaled ratio.
# - The submodels for missing load and condition write dnorm(μ, sigma[1]); NIMBLE's second
#   positional argument is a precision, so the reported "load.sigma[1]" (median 5.85) is the
#   precision of the scaled log10 load (residual sd 0.41), and its half-Cauchy prior is on the
#   precision.
# - sex.mu ~ Beta(17, 3), "centred on the observed sex ratio", has mean 0.85, but sex = 1
#   codes females (14% of known sexes); 15 values are imputed, and sex.mu is also informed by
#   the 671 observed values, so the posterior is near 0.14.
# - Reported slopes (e.g. −0.167 and 1.953) are posterior medians and 95% HDIs *conditional on
#   RJMCMC inclusion*; Fig. 4 and the derived curves use all draws (zeros when excluded).
# - Annual hazard = six-weekly hazard × 52/6 (a 52-week year).
#
# The refit here (CTMultistateProblem) is a reduced version, by maximum likelihood; see its
# docstring for exactly what is kept and what is dropped.

const HOLLANDERS2023_SITES = ["Brindle", "Tuntable", "Bat Cave"]
const HOLLANDERS2023_DRAWS = joinpath("hollanders2023_si", "repo_mfleayi-adults",
                                      "mcmc-samples", "m2-draws.csv")

"""
    hollanders2023_path() -> String

Directory of the bundled Hollanders et al. (2023) mark–recapture data (Dryad
[10.5061/dryad.g1jwstqtb](https://doi.org/10.5061/dryad.g1jwstqtb), CC0): `ch.csv`,
`loads.csv`, `ints.csv`, `dates.csv`, `covs-int.csv`, `covs-surv.csv`, `sex.csv`, `bw.csv`
and `length.csv`, unchanged. See `data/hollanders2023/README.md`.
"""
hollanders2023_path() = normpath(joinpath(@__DIR__, "..", "..", "data", "hollanders2023"))

function _h23_rows(path)
    lines = readlines(path)
    isempty(lines) || (lines[1] = replace(lines[1], '﻿' => ""))
    return [split(strip(l), ',') for l in lines if !isempty(strip(l))]
end

_h23num(s) = (s = strip(s); s == "NA" || isempty(s) ? missing : parse(Float64, s))

"""
    Hollanders2023Data

The *Mixophyes fleayi* robust-design data of Hollanders et al. (2023), as read by
[`load_hollanders2023`](@ref). Dimensions: individual (up to 462), primary occasion (23),
secondary survey (3), site (3: Brindle, Tuntable, Bat Cave Creeks). Fields:

- `sites`, `nind` (462, 136, 88 individuals);
- `ch`: capture codes `[i, t, j, m]`: 0 not seen, 1 seen *Bd*−, 2 seen *Bd*+, `missing`
  where there was no survey (or for padding rows beyond `nind[m]`);
- `loads`: `[i, t, m]` swab load as [`ITS1CopiesPerSwab`](@ref) (`missing` if not
  captured, or captured uninfected or unswabbed); the paper reports log10 ITS copies per swab
  because the copy number of local strains is unknown, so no conversion to ZE is made;
- `intervals`: `[t, m]` length of the interval from primary `t` to `t + 1`, in units of six
  weeks (42 days);
- `dates`: `[t, m]` primary survey dates (the deposit has Tuntable primary 22 as 1/2/20,
  evidently 1 February 2021; kept as given);
- `temperature`, `rainfall`: `[t, m]` mean daily temperature (°C, mean of the SILO daily
  maximum and minimum) and total rainfall (mm) over the interval ending at primary `t`;
- `survey_temperature`, `survey_rainfall`: `[t, j, m]` on the survey day;
- `sex`: `[i, m]` 1 female, 0 male, `missing` unknown;
- `mass`, `length`: `[i, t, m]` body mass (g) and snout–urostyle length (mm).
"""
struct Hollanders2023Data
    sites::Vector{String}
    nind::Vector{Int}
    ch::Array{Union{Missing,Int},4}
    loads::Array{Union{Missing,ITS1CopiesPerSwab},3}
    intervals::Matrix{Float64}
    dates::Matrix{Date}
    temperature::Matrix{Float64}
    rainfall::Matrix{Float64}
    survey_temperature::Array{Float64,3}
    survey_rainfall::Array{Float64,3}
    sex::Matrix{Union{Missing,Int}}
    mass::Array{Union{Missing,Float64},3}
    length::Array{Union{Missing,Float64},3}
end

"""
    load_hollanders2023(; path = hollanders2023_path()) -> Hollanders2023Data

Read the bundled Hollanders et al. (2023) data. The CSV files store R arrays column-major
(`array(unlist(read.csv(...)), dims)` in m2-cmr.qmd); they are reshaped the same way.
"""
function load_hollanders2023(; path::AbstractString = hollanders2023_path())
    nind = [462, 136, 88]
    I, T, J, M = 462, 23, 3, 3
    rd(f) = _h23_rows(joinpath(path, f))
    ch = Array{Union{Missing,Int},4}(missing, I, T, J, M)
    rows = rd("ch.csv")
    length(rows) == I && all(r -> length(r) == T * J * M, rows) ||
        throw(ArgumentError("ch.csv is not $I × $(T * J * M)"))
    for i in 1:I, m in 1:M, j in 1:J, t in 1:T
        v = _h23num(rows[i][t + T * (j - 1) + T * J * (m - 1)])
        ch[i, t, j, m] = v === missing ? missing : Int(v)
    end
    function grid(f)
        r = rd(f)
        A = Array{Union{Missing,Float64},3}(missing, I, T, M)
        for i in 1:I, m in 1:M, t in 1:T
            A[i, t, m] = _h23num(r[i][t + T * (m - 1)])
        end
        return A
    end
    raw = grid("loads.csv")
    loads = Array{Union{Missing,ITS1CopiesPerSwab},3}(missing, I, T, M)
    for k in eachindex(raw)
        raw[k] === missing || (loads[k] = ITS1CopiesPerSwab(raw[k]))
    end
    ints = rd("ints.csv")
    intervals = [parse(Float64, ints[t+1][m]) for t in 1:T-1, m in 1:M]
    dr = rd("dates.csv")
    pd(s) = (p = parse.(Int, split(s, '/')); Date(2000 + p[3], p[2], p[1]))
    dates = [pd(dr[t+1][m]) for t in 1:T, m in 1:M]
    ci = rd("covs-int.csv")
    ci[1] == ["temp.max", "temp.min", "rain"] || throw(ArgumentError("unexpected covs-int.csv"))
    cv(r, k) = parse(Float64, r[k])
    temperature = [(cv(ci[1+t+T*(m-1)], 1) + cv(ci[1+t+T*(m-1)], 2)) / 2 for t in 1:T, m in 1:M]
    rainfall = [cv(ci[1+t+T*(m-1)], 3) for t in 1:T, m in 1:M]
    cs = rd("covs-surv.csv")
    cs[1] == ["temp.max", "temp.min", "vp", "dailyrain"] ||
        throw(ArgumentError("unexpected covs-surv.csv"))
    r3(t, j, m) = cs[1+t+T*(j-1)+T*J*(m-1)]
    stemp = [(cv(r3(t, j, m), 1) + cv(r3(t, j, m), 2)) / 2 for t in 1:T, j in 1:J, m in 1:M]
    srain = [cv(r3(t, j, m), 4) for t in 1:T, j in 1:J, m in 1:M]
    sx = rd("sex.csv")
    sex = Matrix{Union{Missing,Int}}(missing, I, M)
    for i in 1:I, m in 1:M
        v = _h23num(sx[i+1][m])
        sex[i, m] = v === missing ? missing : Int(v)
    end
    return Hollanders2023Data(copy(HOLLANDERS2023_SITES), nind, ch, loads, intervals, dates,
                              temperature, rainfall, stemp, srain, sex, grid("bw.csv"),
                              grid("length.csv"))
end

# the primary-level observed state of m2-cmr.qmd (z.primary): 0 not captured, else the
# largest code over the secondaries (1 Bd−, 2 Bd+)
function _h23_zprimary(d::Hollanders2023Data, i, t, m)
    z = 0
    for j in axes(d.ch, 3)
        c = d.ch[i, t, j, m]
        c === missing || (z = max(z, c))
    end
    return z
end

"""
    hollanders2023_summary_data(d = load_hollanders2023()) -> NamedTuple

The data summaries of the Results ("Sampling summary"), computed as in m2-cmr.qmd:
`individuals` per site, `primary_detections` (individual × primary occasions with a
capture), `detections` (including recaptures in secondary surveys), `infected_detections`
(primary detections of a *Bd*+ frog), `loads` (observed log10 ITS1 copies per swab; their
`mean`, `sd`, `minimum`, `maximum`), and `gains` and `clearances`: the observed infection
changes between consecutive captures of the same frog (the `trans` array: the ordered list of
an individual's observed primary states, gaps ignored).
"""
function hollanders2023_summary_data(d::Hollanders2023Data = load_hollanders2023())
    T = size(d.ch, 2)
    prim = 0
    det = 0
    inf = 0
    gains = 0
    clears = 0
    for m in eachindex(d.nind), i in 1:d.nind[m]
        seq = Int[]
        for t in 1:T
            z = _h23_zprimary(d, i, t, m)
            z > 0 || continue
            prim += 1
            inf += z == 2
            det += count(j -> coalesce(d.ch[i, t, j, m], 0) > 0, axes(d.ch, 3))
            push!(seq, z)
        end
        for k in 1:length(seq)-1
            gains += seq[k] == 1 && seq[k+1] == 2
            clears += seq[k] == 2 && seq[k+1] == 1
        end
    end
    L = [x.log10value for x in skipmissing(d.loads)]
    μ = mean(L)
    return (individuals = copy(d.nind), primary_detections = prim, detections = det,
            infected_detections = inf, loads = L, mean = μ, sd = sqrt(var(L)),
            minimum = minimum(L), maximum = maximum(L), gains = gains, clearances = clears)
end

# centre and divide by two standard deviations (the `scl` of m2-cmr.qmd)
function _h23_scl(x::AbstractArray)
    v = collect(skipmissing(x))
    μ, s = mean(v), sqrt(var(v))
    return map(y -> y === missing ? missing : (y - μ) / (2s), x), μ, s
end

"""
    load_hazard_ratio(log10load; status, load, mean, sd) -> Float64

Ratio of the mortality hazard of an infected frog with swab load `10^log10load` (ITS1
copies per swab) to that of an uninfected frog in Hollanders et al. (2023):
`exp(status + load (log10load − mean) / (2 sd))`, where `status` and `load` are the
log-hazard slopes and `mean`, `sd` the mean and standard deviation of the observed log10
loads used to scale the covariate.
"""
load_hazard_ratio(log10load::Real; status::Real, load::Real, mean::Real, sd::Real) =
    exp(status + load * (log10load - mean) / (2sd))

# --- the authors' posterior draws (not bundled) ---------------------------------------------------

"""
    hollanders2023_draws_path() -> String

The authors' 10 000 posterior draws of the mark–recapture model (`m2-draws.csv`, 54 MB) from
the repository mhollanders/mfleayi-adults (commit a325a4c), which has no licence and is not
bundled. The path is `ENV["CHYTRIDIOMYCOSIS_HOLLANDERS2023_DRAWS"]` if set, else
`papers/hollanders2023_si/repo_mfleayi-adults/mcmc-samples/m2-draws.csv` under the first
ancestor of the package directory that has it. Throws an `ArgumentError` if neither exists.
"""
function hollanders2023_draws_path()
    p = _hollanders2023_draws_find()
    p === nothing && throw(ArgumentError(
        "Hollanders et al. 2023 posterior draws not found: download mcmc-samples/m2-draws.csv " *
        "from mhollanders/mfleayi-adults and set ENV[\"CHYTRIDIOMYCOSIS_HOLLANDERS2023_DRAWS\"]"))
    return p
end

function _hollanders2023_draws_find()
    e = get(ENV, "CHYTRIDIOMYCOSIS_HOLLANDERS2023_DRAWS", "")
    isempty(e) || return isfile(e) ? e : nothing
    d = pkgdir(@__MODULE__)
    while true
        c = joinpath(d, "papers", HOLLANDERS2023_DRAWS)
        isfile(c) && return c
        parent = dirname(d)
        parent == d && return nothing
        d = parent
    end
end

const _H23_DRAWS = Ref{Union{Nothing,Dict{String,Vector{Float64}}}}(nothing)

function _h23_draws()
    _H23_DRAWS[] === nothing || return _H23_DRAWS[]
    path = hollanders2023_draws_path()
    lines = readlines(path)
    # quoted header names contain commas ("eps.site[1, 1]")
    hdr = [m.captures[1] for m in eachmatch(r"\"([^\"]*)\"", lines[1])]
    want = [k for (k, h) in enumerate(hdr) if occursin(r"^(phi|psi12|psi21|p|load)\.(alpha|beta|w)|^eps\.site|^sigma\.survey|^p\.sigma", h)]
    cols = Dict(hdr[k] => Vector{Float64}(undef, length(lines) - 1) for k in want)
    for (r, l) in enumerate(lines[2:end])
        f = split(l, ',')
        for k in want
            cols[hdr[k]][r] = f[k] == "NA" ? NaN : parse(Float64, f[k])
        end
    end
    _H23_DRAWS[] = cols
    return cols
end

# 95% highest-density interval (shortest interval holding ⌊0.95 n⌋ + 1 sorted draws)
function _hdi(v::AbstractVector; mass::Real = 0.95)
    s = sort(v)
    n = length(s)
    k = floor(Int, mass * n)
    w = [s[i+k] - s[i] for i in 1:n-k]
    i = argmin(w)
    return s[i], s[i+k]
end

function _h23_sum(v)
    lo, hi = _hdi(v)
    return (median = median(v), lower = lo, upper = hi)
end

"""
    hollanders2023_posterior(; draws = nothing) -> NamedTuple

Posterior summaries (median and 95% highest-density interval, as the paper reports) computed
from the authors' own MCMC draws ([`hollanders2023_draws_path`](@ref); not bundled). This
reproduces the paper's reporting from its posterior; it is not a refit. Fields (each
`(median, lower, upper)`):

- slopes *conditional on RJMCMC inclusion*, as printed: `phi_condition`, `phi_temperature`,
  `phi_status`, `phi_load`, `gain_temperature`, `gain_rain`, `gain_condition`,
  `clear_temperature`, `clear_rain`, `clear_condition`, `clear_load`, `p_status`;
  `inclusion` holds the inclusion frequencies of the same slopes;
- site-averaged six-weekly rates (geometric mean over the three site intercepts, as
  m2-cmr.qmd's `phi.alpha.mu`): `mortality`, `gain`, `clearance`, and `recapture`
  (probability, mean on the logit scale);
- `annual_hazard` (`mortality × 52/6`), `annual_survival` (`exp(−annual_hazard)`),
  `gain_probability` and `clearance_probability` (`1 − exp(−rate)`, the paper's conversion);
- `mortality_sites`, `gain_sites`, `clearance_sites`: site medians (Brindle, Tuntable, Bat Cave);
- `survey_sd`: the sd of the survey effects on log mortality, gaining and clearing;
- `hazard_ratio`: `(log10load, summary)` pairs of the infected/uninfected mortality ratio at
  `10^3.17`, `10^5` and `10^6` ITS1 copies per swab with all draws (model-averaged, as in
  Fig. 4), scaled by the mean and sd of the bundled loads.
"""
function hollanders2023_posterior(; draws = nothing)
    c = draws === nothing ? _h23_draws() : draws
    cond(b, w) = _h23_sum(c[b][c[w] .== 1])
    sl = (phi_condition = ("phi.beta[1]", "phi.w[1]"), phi_temperature = ("phi.beta[2]", "phi.w[2]"),
          phi_status = ("phi.beta[4]", "phi.w[4]"), phi_load = ("phi.beta[5]", "phi.w[5]"),
          gain_temperature = ("psi12.beta[1]", "psi12.w[1]"), gain_rain = ("psi12.beta[2]", "psi12.w[2]"),
          gain_condition = ("psi12.beta[4]", "psi12.w[4]"),
          clear_temperature = ("psi21.beta[1]", "psi21.w[1]"), clear_rain = ("psi21.beta[2]", "psi21.w[2]"),
          clear_condition = ("psi21.beta[4]", "psi21.w[4]"), clear_load = ("psi21.beta[5]", "psi21.w[5]"),
          p_status = ("p.beta[6]", "p.w[6]"))
    slopes = map(x -> cond(x...), sl)
    inclusion = map(x -> mean(c[x[2]]), sl)
    n = length(c["phi.alpha"])
    sitemean(a, k) = [exp((3log(c[a][r]) + sum(c["eps.site[$k, $m]"][r] for m in 1:3)) / 3) for r in 1:n]
    phi, gain, clear = sitemean("phi.alpha", 1), sitemean("psi12.alpha", 2), sitemean("psi21.alpha", 3)
    lg(p) = log(p / (1 - p))
    rec = [_logistic((3lg(c["p.alpha"][r]) + sum(c["eps.site[4, $m]"][r] for m in 1:3)) / 3) for r in 1:n]
    site(a, k, m) = median(exp.(log.(c[a]) .+ c["eps.site[$k, $m]"]))
    s = hollanders2023_summary_data()
    hr = [(lg10, _h23_sum([load_hazard_ratio(lg10; status = c["phi.beta[4]"][r],
                                              load = c["phi.beta[5]"][r], mean = s.mean, sd = s.sd)
                           for r in 1:n])) for lg10 in (3.17, 5.0, 6.0)]
    annual = phi .* (52 / 6)
    return merge(slopes, (inclusion = inclusion, mortality = _h23_sum(phi), gain = _h23_sum(gain),
        clearance = _h23_sum(clear), recapture = _h23_sum(rec), annual_hazard = _h23_sum(annual),
        annual_survival = _h23_sum(exp.(-annual)), gain_probability = _h23_sum(1 .- exp.(-gain)),
        clearance_probability = _h23_sum(1 .- exp.(-clear)),
        mortality_sites = [site("phi.alpha", 1, m) for m in 1:3],
        gain_sites = [site("psi12.alpha", 2, m) for m in 1:3],
        clearance_sites = [site("psi21.alpha", 3, m) for m in 1:3],
        survey_sd = [median(c["sigma.survey[$k]"]) for k in 1:3], hazard_ratio = hr))
end

# --- the reduced refit -----------------------------------------------------------------------------

const CTMS_NAMES = [:log_mortality_1, :log_mortality_2, :log_mortality_3, :mortality_temperature,
    :mortality_status, :mortality_load, :log_gain_1, :log_gain_2, :log_gain_3, :gain_temperature,
    :gain_rain, :log_clear_1, :log_clear_2, :log_clear_3, :clear_temperature, :clear_rain,
    :clear_load, :logit_p_1, :logit_p_2, :logit_p_3, :p_female, :p_temperature, :p_rain,
    :p_status, :p_load, :log_sd_individual]

struct _H23Ind
    site::Int
    first::Int
    zfirst::Int
    codes::Matrix{Union{Missing,Int}}    # T × J, 1 seen Bd−, 2 seen Bd+, 3 not seen
    x::Vector{Union{Missing,Float64}}    # scaled log10 load per primary
    sex::Union{Missing,Int}
    umean::Float64                       # individual load effect | observed loads
    usd::Float64
end

"""
    CTMultistateProblem(data = load_hollanders2023(); nodes = (load = 5, individual = 3,
                        detection = 7), swab_final = false)

Reduced maximum-likelihood refit of the continuous-time robust-design multistate model of
Hollanders et al. (2023). Each frog's capture history is an HMM over (alive *Bd*−, alive
*Bd*+, dead) between primary occasions, with transition matrices
[`ctmc_transition`](@ref)`(`[`infection_generator`](@ref)`(ψ12, ψ21, φ1, φ2), τ)` over the
unequal intervals `τ` (six-week units) and recapture at each secondary survey. Conditioned on
first capture. Kept as in the authors' code:

- site-specific intercepts of mortality, gaining, clearing (log six-weekly hazard) and
  recapture (logit), covariates centred and divided by two sd as coded;
- mortality: interval temperature, *Bd* status and load (log10 ITS1 copies per swab at the
  start of the interval, scaled); gaining: temperature and rainfall; clearing: temperature,
  rainfall and load; recapture: sex, survey-day temperature, rainfall, *Bd* status and load,
  and a normal individual effect;
- frogs not swabbed on the final primary (`swab_final = false`: a *Bd*+ frog seen then is
  "seen, *Bd*−", the code's δ).

Missing loads (every occasion without a positive swab) are integrated out, not imputed by
MCMC: the load submodel (normal, with temperature, rainfall and a normal individual effect,
fitted first by maximum likelihood to the observed loads; [`hollanders2023_load_model`](@ref))
gives, for each frog, the distribution of its individual effect given its observed loads, and
the latent loads are integrated by Gauss–Hermite quadrature (`nodes.load` per occasion,
`nodes.individual` over the individual load effect, `nodes.detection` over the recapture
individual effect). Unknown sex (15 frogs) is a two-point mixture at the observed female
fraction. Dropped: site and survey random effects (sd 0.33, 0.49, 0.29 on the log hazards),
the secondary-survey effects on recapture, body condition and every interaction, and
RJMCMC. Parameters are maximised over the unconstrained vector `CTMS_NAMES` (26).

`solve(prob)` returns a [`CTMultistateFit`](@ref). It is a LogDensityProblems target whose
log density is the log-likelihood (flat prior).
"""
struct CTMultistateProblem
    data::Hollanders2023Data
    inds::Vector{_H23Ind}
    temperature::Matrix{Float64}          # scaled, [t, m] for the interval ending at t
    rainfall::Matrix{Float64}             # scaled and divided by the interval length
    survey_temperature::Array{Float64,3}  # scaled, [t, j, m]
    load_mean::Float64
    load_sd::Float64
    load_model::NamedTuple
    female_fraction::Float64
    nodes::NamedTuple
    swab_final::Bool
    names::Vector{Symbol}
end

# Gauss–Hermite nodes and weights for E[f(Z)], Z ~ Normal(0, 1) (Golub–Welsch).
function _gauss_hermite(n::Integer)
    n == 1 && return [0.0], [1.0]
    J = LinearAlgebra.SymTridiagonal(zeros(n), [sqrt(k / 1) for k in 1:n-1])
    F = eigen(J)
    return F.values, F.vectors[1, :] .^ 2
end

"""
    hollanders2023_load_model(data = load_hollanders2023()) -> NamedTuple

Maximum-likelihood fit of the authors' imputation submodel to the 248 observed loads: the
scaled log10 load `x = (log10 L − mean)/(2 sd)` of frog `i` at primary `t` is
`Normal(β_T temp_t + β_R rain_t + u_i, σ²)` with `u_i ~ Normal(0, σ_u²)` and no intercept
(the code fixes `load.alpha = 0`). The temperature × rainfall term, included in 55% of the
authors' RJMCMC draws, is left out. Returns `temperature`, `rain`, `sd` (σ), `sd_individual`
(σ_u) and `loglik`, on the scaled-load scale (authors' posterior medians: −0.359, 0.231,
precision 5.85 i.e. σ = 0.41, and σ_u = 0.241).
"""
function hollanders2023_load_model(d::Hollanders2023Data = load_hollanders2023())
    tsc, = _h23_scl(d.temperature)
    rsc = _h23_rain(d)
    lx, = _h23_scl(map(x -> x === missing ? missing : x.log10value, d.loads))
    groups = Vector{Tuple{Vector{Float64},Matrix{Float64}}}()
    for m in eachindex(d.nind), i in 1:d.nind[m]
        ts = [t for t in axes(lx, 2) if lx[i, t, m] !== missing]
        isempty(ts) && continue
        push!(groups, (Float64[lx[i, t, m] for t in ts],
                       hcat(Float64[tsc[t, m] for t in ts], Float64[rsc[t, m] for t in ts])))
    end
    function nll(θ)
        β = θ[1:2]
        s2, u2 = exp(2θ[3]), exp(2θ[4])
        ll = 0.0
        for (x, X) in groups
            r = x .- X * β
            n = length(r)
            a = u2 / (s2 + n * u2)
            ll -= 0.5 * (n * log(2π) + (n - 1) * log(s2) + log(s2 + n * u2) +
                         (sum(abs2, r) - a * sum(r)^2) / s2)
        end
        return -ll
    end
    r = _ms_bfgs(nll, [0.0, 0.0, log(0.4), log(0.2)]; gtol = 1e-7)
    return (temperature = r.x[1], rain = r.x[2], sd = exp(r.x[3]), sd_individual = exp(r.x[4]),
            loglik = -r.f)
end

# rain.int.scl of m2-cmr.qmd: scl(rain) reshaped to [t, m], divided by the length of the
# interval ending at t (1 for the first primary)
function _h23_rain(d::Hollanders2023Data)
    r, = _h23_scl(d.rainfall)
    T = size(r, 1)
    return [r[t, m] / (t == 1 ? 1.0 : d.intervals[t-1, m]) for t in 1:T, m in axes(r, 2)]
end

function CTMultistateProblem(d::Hollanders2023Data = load_hollanders2023();
                             nodes::NamedTuple = (load = 5, individual = 3, detection = 7),
                             swab_final::Bool = false)
    T, J = size(d.ch, 2), size(d.ch, 3)
    tsc, = _h23_scl(d.temperature)
    rsc = _h23_rain(d)
    ssc, = _h23_scl(d.survey_temperature)
    lx, lμ, lsd = _h23_scl(map(x -> x === missing ? missing : x.log10value, d.loads))
    lm = hollanders2023_load_model(d)
    known = collect(skipmissing(d.sex))
    fem = mean(known)
    inds = _H23Ind[]
    for m in eachindex(d.nind), i in 1:d.nind[m]
        zp = [_h23_zprimary(d, i, t, m) for t in 1:T]
        f = findfirst(>(0), zp)
        (f === nothing || f == T) && continue
        codes = Matrix{Union{Missing,Int}}(missing, T, J)
        for t in 1:T, j in 1:J
            c = d.ch[i, t, j, m]
            codes[t, j] = c === missing ? missing : (c == 0 ? 3 : c)
        end
        x = Union{Missing,Float64}[lx[i, t, m] for t in 1:T]
        obs = [t for t in 1:T if x[t] !== missing]
        r = [x[t] - lm.temperature * tsc[t, m] - lm.rain * rsc[t, m] for t in obs]
        v = 1 / (1 / lm.sd_individual^2 + length(obs) / lm.sd^2)
        push!(inds, _H23Ind(m, f, zp[f], codes, x, d.sex[i, m], v * sum(r; init = 0.0) / lm.sd^2,
                            sqrt(v)))
    end
    return CTMultistateProblem(d, inds, Float64.(tsc), rsc, Float64.(ssc), lμ, lsd, lm, fem,
                               nodes, swab_final, copy(CTMS_NAMES))
end

LogDensityProblems.dimension(p::CTMultistateProblem) = length(p.names)
LogDensityProblems.capabilities(::Type{<:CTMultistateProblem}) = LogDensityProblems.LogDensityOrder{0}()

# Emission likelihoods of primary t for state U (scalar) and I at load x (via closure).
@inline function _h23_emit(prob, θ, ind::_H23Ind, t, sex, ε, xI, swabbed)
    m = ind.site
    eU, eI, seen = 1.0, 1.0, false
    for j in axes(ind.codes, 2)
        c = ind.codes[t, j]
        c === missing && continue
        η1 = θ[17+m] + θ[21] * sex + θ[22] * prob.survey_temperature[t, j, m] +
             θ[23] * prob.rainfall[t, m] + ε
        p1 = _logistic(η1)
        p2 = _logistic(η1 + θ[24] + θ[25] * xI)
        if c == 3
            eU *= 1 - p1
            eI *= 1 - p2
        else
            seen = true
            if c == 1
                eU *= p1
                eI *= swabbed ? 0.0 : p2
            else
                eU = 0.0
                eI *= swabbed ? p2 : 0.0
            end
        end
    end
    return eU, eI, seen
end

# Transition rows from primary t (load x at t) to t + 1.
@inline function _h23_trans(prob, θ, ind::_H23Ind, t, x)
    m = ind.site
    tt, rr = prob.temperature[t+1, m], prob.rainfall[t+1, m]
    φ1 = exp(θ[m] + θ[4] * tt)
    φ2 = φ1 * exp(θ[5] + θ[6] * x)
    g = exp(θ[6+m] + θ[10] * tt + θ[11] * rr)
    c = exp(θ[11+m] + θ[15] * tt + θ[16] * rr + θ[17] * x)
    a11, a12, a21, a22 = -(g + φ1), g, c, -(c + φ2)
    P = _ctmc_absorbing2(a11, a12, a21, a22, prob.data.intervals[t, m])
    return P
end

# log-likelihood of one frog (outer quadrature over its load and recapture individual effects
# and, if unknown, its sex)
function _h23_loglik(prob::CTMultistateProblem, θ::AbstractVector, ind::_H23Ind, gh)
    T = size(ind.codes, 1)
    zl, wl = gh.load
    zu, wu = gh.individual
    ze, we = gh.detection
    σε = exp(θ[26])
    lm = prob.load_model
    sexes = ind.sex === missing ? ((1, prob.female_fraction), (0, 1 - prob.female_fraction)) :
            ((ind.sex, 1.0),)
    K = length(zl)
    xs = zeros(T, K)
    total = 0.0
    αI = zeros(K)
    nI = zeros(K)
    for (a, u) in enumerate(zu)
        uval = ind.umean + ind.usd * u
        for t in ind.first:T
            if ind.x[t] === missing
                μ = lm.temperature * prob.temperature[t, ind.site] +
                    lm.rain * prob.rainfall[t, ind.site] + uval
                for k in 1:K
                    xs[t, k] = μ + lm.sd * zl[k]
                end
            else
                xs[t, :] .= ind.x[t]
            end
        end
        # transition rows per occasion and node (independent of the recapture effect)
        PU = zeros(T, 3)          # Σ_k w_k P_k[1, :]
        PI = zeros(T, K, 3)       # P_k[2, :]
        for t in ind.first:T-1
            for k in 1:K
                P = _h23_trans(prob, θ, ind, t, xs[t, k])
                for c in 1:3
                    PU[t, c] += wl[k] * P[1, c]
                    PI[t, k, c] = P[2, c]
                end
            end
        end
        inner = 0.0
        for (b, e) in enumerate(ze), (sx, ps) in sexes
            ε = σε * e
            # initial: known state at first capture, load nodes with weights wl
            aU = ind.zfirst == 1 ? 1.0 : 0.0
            for k in 1:K
                αI[k] = ind.zfirst == 2 ? wl[k] : 0.0
            end
            aD = 0.0
            ll = 0.0
            for t in ind.first+1:T
                mU = aU * PU[t-1, 1]
                mI = aU * PU[t-1, 2]
                mD = aU * PU[t-1, 3] + aD
                for k in 1:K
                    mU += αI[k] * PI[t-1, k, 1]
                    mI += αI[k] * PI[t-1, k, 2]
                    mD += αI[k] * PI[t-1, k, 3]
                end
                swabbed = prob.swab_final || t < T
                eU, _, seen = _h23_emit(prob, θ, ind, t, sx, ε, 0.0, swabbed)
                aU = mU * eU
                c = aU
                for k in 1:K
                    _, eI, _ = _h23_emit(prob, θ, ind, t, sx, ε, xs[t, k], swabbed)
                    αI[k] = mI * wl[k] * eI
                    c += αI[k]
                end
                aD = seen ? 0.0 : mD
                c += aD
                c > 0 || (ll = -Inf; break)
                ll += log(c)
                aU /= c
                aD /= c
                αI ./= c
            end
            inner += we[b] * ps * exp(ll)
        end
        total += wu[a] * inner
    end
    return log(total)
end

function _h23_nodes(prob::CTMultistateProblem)
    return (load = _gauss_hermite(prob.nodes.load), individual = _gauss_hermite(prob.nodes.individual),
            detection = _gauss_hermite(prob.nodes.detection))
end

"""
    pointwise_loglikelihood(prob::CTMultistateProblem, θ) -> Vector

Log-likelihood of each frog's capture history at the unconstrained parameters `θ`.
"""
function pointwise_loglikelihood(prob::CTMultistateProblem, θ::AbstractVector)
    gh = _h23_nodes(prob)
    out = zeros(length(prob.inds))
    nt = Threads.nthreads()
    if nt > 1
        chunks = collect(Iterators.partition(eachindex(prob.inds), cld(length(prob.inds), nt)))
        tasks = [Threads.@spawn begin
                     for i in ch
                         out[i] = _h23_loglik(prob, θ, prob.inds[i], gh)
                     end
                 end for ch in chunks]
        foreach(wait, tasks)
    else
        for i in eachindex(prob.inds)
            out[i] = _h23_loglik(prob, θ, prob.inds[i], gh)
        end
    end
    return out
end

function LogDensityProblems.logdensity(prob::CTMultistateProblem, θ::AbstractVector)
    ll = sum(pointwise_loglikelihood(prob, θ))
    return isnan(ll) ? -Inf : ll
end

# The same likelihood for one frog at one value of each outer variable (load individual
# effect u, recapture effect ε, sex), as an HMM with states (U, I at load node 1…K, dead)
# passed to hmm_forward: the check that the fast recursion above is that HMM.
function _h23_hmm(prob::CTMultistateProblem, θ, ind::_H23Ind, u, ε, sex)
    T = size(ind.codes, 1)
    zl, wl = _gauss_hermite(prob.nodes.load)
    K = length(zl)
    lm = prob.load_model
    xs = [ind.x[t] === missing ?
          lm.temperature * prob.temperature[t, ind.site] + lm.rain * prob.rainfall[t, ind.site] +
          ind.umean + ind.usd * u + lm.sd * zl[k] : Float64(ind.x[t]) for t in 1:T, k in 1:K]
    n = K + 2
    init = zeros(n)
    ind.zfirst == 1 ? (init[1] = 1.0) : (init[2:K+1] .= wl)
    Ps = Matrix{Float64}[]
    L = ones(T - ind.first + 1, n)
    for t in ind.first:T-1
        P = zeros(n, n)
        for k in 1:K
            Pk = _h23_trans(prob, θ, ind, t, xs[t, k])
            P[1, 1] += wl[k] * Pk[1, 1]
            P[1, 2:K+1] .+= wl[k] * Pk[1, 2] .* wl
            P[1, n] += wl[k] * Pk[1, 3]
            P[1+k, 1] = Pk[2, 1]
            P[1+k, 2:K+1] .= Pk[2, 2] .* wl
            P[1+k, n] = Pk[2, 3]
        end
        P[n, n] = 1.0
        push!(Ps, P)
        s = t + 1
        swabbed = prob.swab_final || s < T
        eU, _, seen = _h23_emit(prob, θ, ind, s, sex, ε, 0.0, swabbed)
        L[s-ind.first+1, 1] = eU
        for k in 1:K
            L[s-ind.first+1, 1+k] = _h23_emit(prob, θ, ind, s, sex, ε, xs[s, k], swabbed)[2]
        end
        L[s-ind.first+1, n] = seen ? 0.0 : 1.0
    end
    return init, Ps, L
end

"""
    CTMultistateFit

Maximum-likelihood fit of a [`CTMultistateProblem`](@ref): `problem`, `θ` (unconstrained,
names `problem.names`), `vcov` (inverse observed information; `NaN` if singular), `loglik`,
`gradient` at the optimum and `iterations`. [`hollanders2023_summary`](@ref) gives the
quantities the paper reports, with Wald intervals.
"""
struct CTMultistateFit
    problem::CTMultistateProblem
    θ::Vector{Float64}
    vcov::Matrix{Float64}
    loglik::Float64
    gradient::Vector{Float64}
    iterations::Int
end

"""
    hollanders2023_start() -> Vector{Float64}

A neutral starting point for `solve(::CTMultistateProblem)`: six-weekly hazards of 0.1
(mortality), 0.3 (gaining) and 1 (clearing), recapture probability 0.3, individual sd 1 and
all slopes 0. These are round numbers, not the published estimates.
"""
hollanders2023_start() = vcat(fill(log(0.1), 3), 0.0, 0.0, 0.0, fill(log(0.3), 3), 0.0, 0.0,
                              fill(0.0, 3), 0.0, 0.0, 0.0, fill(log(0.3 / 0.7), 3),
                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

# solve(prob::CTMultistateProblem; x0 = hollanders2023_start(), iters = 500, hessian = true)
# Documented in CTMultistateProblem.
function CommonSolve.solve(prob::CTMultistateProblem; x0::AbstractVector = hollanders2023_start(),
                           iters::Integer = 500, hessian::Bool = true, gtol::Real = 1e-3)
    f(θ) = -LogDensityProblems.logdensity(prob, θ)
    r = _ms_bfgs(f, x0; iters, gtol)
    V = hessian ? first(_ms_vcov(-_ms_hessian(f, r.x))) : fill(NaN, length(r.x), length(r.x))
    return CTMultistateFit(prob, r.x, V, -r.f, r.gradient, r.iterations)
end

"""
    hollanders2023_summary(fit::CTMultistateFit; draws = 4000, rng = Random.Xoshiro(2023))
        -> NamedTuple

The quantities Hollanders et al. (2023) report, from a reduced refit: each a named tuple
`(estimate, lower, upper)` with a 95% interval from `draws` normal draws of the parameters
(mean `fit.θ`, covariance `fit.vcov`; the estimate is at `fit.θ`). Slopes (`phi_temperature`,
`phi_status`, `phi_load`, `gain_temperature`, `gain_rain`, `clear_temperature`, `clear_rain`,
`clear_load`, `p_status`, `p_load`) are on the log-hazard or logit scale per two sd of the
covariate; `mortality`, `gain`, `clearance` are six-weekly hazards at mean covariates,
geometric means over the sites; `annual_hazard = mortality × 52/6`, `annual_survival`,
`gain_probability`, `clearance_probability` (`1 − exp(−rate)`); `recapture` (inverse logit of
the mean site intercept); `mortality_sites`, `gain_sites`, `clearance_sites` (vectors of
estimates); `hazard_ratio` at `10^5` and `10^6` ITS1 copies per swab
([`load_hazard_ratio`](@ref)); and the fitted `load_model`.
"""
function hollanders2023_summary(fit::CTMultistateFit; draws::Integer = 4000,
                                rng::AbstractRNG = Random.Xoshiro(2023))
    θ = fit.θ
    ok = all(isfinite, fit.vcov)
    Θ = if ok
        C = _w22_chol(fit.vcov)
        [θ .+ C * randn(rng, length(θ)) for _ in 1:draws]
    else
        [θ]
    end
    p = fit.problem
    q(g) = (v = [g(t) for t in Θ];
            (estimate = g(θ), lower = ok ? quantile(v, 0.025) : NaN, upper = ok ? quantile(v, 0.975) : NaN))
    sl(j) = q(t -> t[j])
    gm(r) = t -> exp(sum(t[r]) / 3)
    hr(lg) = q(t -> load_hazard_ratio(lg; status = t[5], load = t[6], mean = p.load_mean, sd = p.load_sd))
    return (phi_temperature = sl(4), phi_status = sl(5), phi_load = sl(6), gain_temperature = sl(10),
            gain_rain = sl(11), clear_temperature = sl(15), clear_rain = sl(16), clear_load = sl(17),
            p_female = sl(21), p_status = sl(24), p_load = sl(25),
            mortality = q(gm(1:3)), gain = q(gm(7:9)), clearance = q(gm(12:14)),
            annual_hazard = q(t -> gm(1:3)(t) * 52 / 6),
            annual_survival = q(t -> exp(-gm(1:3)(t) * 52 / 6)),
            gain_probability = q(t -> 1 - exp(-gm(7:9)(t))),
            clearance_probability = q(t -> 1 - exp(-gm(12:14)(t))),
            recapture = q(t -> _logistic(sum(t[18:20]) / 3)),
            mortality_sites = exp.(θ[1:3]), gain_sites = exp.(θ[7:9]), clearance_sites = exp.(θ[12:14]),
            individual_sd = q(t -> exp(t[26])),
            hazard_ratio = [(5.0, hr(5.0)), (6.0, hr(6.0))], load_model = p.load_model)
end
