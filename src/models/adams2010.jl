# Adams MJ, Chelgren ND, Reinitz D, Cole RA, Rachowicz LJ, Galvan S, McCreary B, Pearl CA,
# Bailey LL, Bettaso J, Bull EL, Leu M (2010) Using occupancy models to understand the
# distribution of an amphibian pathogen, Batrachochytrium dendrobatidis. Ecol Appl
# 20:289-302, doi:10.1890/08-2319.1. Supplement: Ecological Archives A020-006 (figshare
# collection doi:10.6084/m9.figshare.c.3294440, CC0): Appendix A (a figure) and Appendix B
# (covariate correlations). No detection histories are deposited.
#
# Checked against the published PDF (lab/papers/adams2010/adams2010.pdf), 2026-11-20.
# Model: single-season occupancy (MacKenzie et al. 2002), fitted in Program MARK 5.1, with a
# population (site) of one anuran species as the site and each animal swabbed once as a
# survey. p = p′ p″ is the probability an animal from a population with Bd tests positive
# (p′ that it carries Bd, p″ that the test finds it; p″ assumed ≈ 1), w (ψ) the probability
# that Bd occurs in the population. Logit links; continuous covariates standardised to mean
# 0, SD 1. Native anurans (72 populations analysed of 78) and bullfrogs (40) were analysed
# separately. Best models: natives p(larv, date, elev, date×elev, hf), w(ralu); bullfrogs
# p(larv, date, stag), w(lat) (Tables 3-6). The likelihood is SingleSeasonOccupancyModel in
# src/observation/occupancy.jl; without the data, only quantities derived from the printed
# coefficients (Tables 4 and 6, in the compendium) can be computed.
#
# Two printed numbers do not follow from the printed coefficients (validation targets):
# - ŵ = 53.4% (95% CI 36.6-69.4%) for natives "based on the top-ranked model", whereas the
#   Table 4 intercept 1.053 (CI −0.035 to 2.141) gives 0.741 (0.491-0.895) for populations
#   other than R. luteiventris. (For bullfrogs, ŵ = 66.4% (42.9-83.8%) is exactly the Table 6
#   intercept at mean latitude.)
# - The native odds ratios for larvae (0.10, CI 0.04-0.20) and per human-footprint point (2.37,
#   CI 1.63-3.43): Table 4 gives exp(−2.432) = 0.088 (0.041-0.188) and, per point,
#   exp(0.872/2.023) = 1.54 (2.39 per standard deviation).
# Appendix A: with N animals in a population, the probability that at least one carries Bd
# is ψ[1 − (1 − p′)^N], near ψ for p′ > 0.05 and N > 50.

"""
    finite_population_occurrence(ψ, p′, N)

Probability `ψ[1 − (1 − p′)^N]` that the pathogen occurs in a population of `N` animals and
at least one of them could test positive, with per-animal probability `p′` (Adams et al.
2010, Appendix A). Divided by `ψ` it is the ratio plotted in their Fig. A1, which approaches
1 as `N` and `p′` increase.
"""
function finite_population_occurrence(ψ::Real, p′::Real, N::Real)
    (0 <= ψ <= 1 && 0 <= p′ <= 1) || throw(DomainError((ψ, p′), "probabilities must be in [0, 1]"))
    N >= 0 || throw(DomainError(N, "N must be non-negative"))
    return ψ * (1 - (1 - p′)^N)
end

const ADAMS2010_SETS = (native = "native", bullfrog = "bullfrog")

"""
    adams2010_coefficients(analysis::Symbol) -> NamedTuple

Logit-scale coefficients of the best a priori model of Adams et al. (2010) for
`analysis = :native` (Table 4: `p` on larvae, date, elevation, date × elevation and human
footprint; `ψ` on *Rana luteiventris*) or `:bullfrog` (Table 6: `p` on larvae, date and
stagnant water; `ψ` on latitude), read from the compendium. Fields `p_<covariate>` and
`ψ_<covariate>` hold `(β, se, scale)`, where `scale` is the standard deviation of a
standardised continuous covariate (`1` for intercepts and indicators).
"""
function adams2010_coefficients(analysis::Symbol)
    haskey(ADAMS2010_SETS, analysis) || throw(ArgumentError("analysis must be :native or :bullfrog"))
    q = parameter_set(:adams2010; set = ADAMS2010_SETS[analysis])
    terms = analysis === :native ?
        (:p_intercept, :p_larv, :p_date, :p_elev, :p_date_elev, :p_hf, :ψ_intercept, :ψ_ralu) :
        (:p_intercept, :p_larv, :p_date, :p_stag, :ψ_intercept, :ψ_lat)
    get1(k) = hasproperty(q, k) ? getproperty(q, k) : 1.0
    return NamedTuple{terms}(map(t -> (β = getproperty(q, Symbol("β_", t)),
                                       se = getproperty(q, Symbol("se_", t)),
                                       scale = get1(Symbol("scale_", t))), terms))
end

"""
    adams2010_odds_ratio(analysis, term; per = 1, level = 0) -> Float64

Odds ratio of a coefficient of [`adams2010_coefficients`](@ref) for a change of `per` units
of the original covariate, `exp(β per / scale)`. With `level = 0.95`, returns the Wald
interval bound pair instead, `exp((β ∓ 1.96 se) per / scale)`.
"""
function adams2010_odds_ratio(analysis::Symbol, term::Symbol; per::Real = 1, level::Real = 0)
    c = getproperty(adams2010_coefficients(analysis), term)
    f = per / c.scale
    level == 0 && return exp(c.β * f)
    z = _normal_quantile((1 + level) / 2)
    return (exp((c.β - z * c.se) * f), exp((c.β + z * c.se) * f))
end

_normal_quantile(q) = Distributions.quantile(Normal(), q)

"""
    adams2010_occupancy(analysis; level = 0) -> Float64

Occupancy probability `ψ` at the intercept of the best model (`:native`: populations other
than *R. luteiventris*; `:bullfrog`: mean latitude), `expit(β)`, or with `level = 0.95` its
Wald interval from the coefficient's standard error.
"""
function adams2010_occupancy(analysis::Symbol; level::Real = 0)
    c = adams2010_coefficients(analysis).ψ_intercept
    level == 0 && return _expit(c.β)
    z = _normal_quantile((1 + level) / 2)
    return (_expit(c.β - z * c.se), _expit(c.β + z * c.se))
end
