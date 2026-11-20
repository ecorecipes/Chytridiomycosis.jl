# Adams MJ, Chelgren ND, Reinitz D, Cole RA, Rachowicz LJ, Galvan S, McCreary B, Pearl CA,
# Bailey LL, Bettaso J, Bull EL, Leu M (2010) Using occupancy models to understand the
# distribution of an amphibian pathogen, Batrachochytrium dendrobatidis. Ecol Appl
# 20:289-302, doi:10.1890/08-2319.1. Supplement: Ecological Archives A020-006 (figshare
# collection doi:10.6084/m9.figshare.c.3294440, CC0): Appendix A (a figure) and Appendix B
# (covariate correlations). No detection histories are deposited.
#
# SOURCE. The article is paywalled and was not read; only its abstract and the two
# appendices were. What is implemented is limited to that:
# - The model is single-season occupancy (MacKenzie et al. 2002) applied to Bd in amphibian
#   populations: a population (native species or the non-native bullfrog, Rana catesbeiana)
#   is "occupied" if Bd occurs in it, and each animal tested is a survey that detects Bd with
#   probability p if it occurs (abstract). Its likelihood is occupancy_loglikelihood
#   (src/observation/occupancy.jl). The covariate models of detection (development stage,
#   day of year, elevation, human activity) and of occurrence were not read and are not
#   implemented.
# - Appendix A: when the number N of animals that could be sampled is finite, the
#   probability that Bd occurs in the population and could be found is ψ[1 − (1 − p′)^N],
#   whose ratio to ψ approaches 1 as N and p′ increase (Fig. A1 caption). The meaning of p′
#   is taken from that caption only (a per-animal probability); the paper's derivation was
#   not read.

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
