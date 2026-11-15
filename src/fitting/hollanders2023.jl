# Reduced refit of Hollanders et al. (2023) Ecol Appl (doi:10.1002/eap.2724) to the bundled
# Dryad data (data/hollanders2023/), and helpers for its validation targets. See
# src/models/hollanders2023.jl for the model, what the reduced refit keeps and drops, and the
# printed-versus-coded notes. This reproduces a published analysis; it is not a new model.
#
# The authors ran NIMBLE for 2–3 days (4 chains × 30 000 iterations, RJMCMC, latent states and
# imputed covariates sampled). Here the latent states and loads are summed and integrated
# exactly (forward algorithm, Gauss–Hermite quadrature), the 26 parameters of the reduced model
# are maximised by BFGS with finite-difference gradients, and intervals are Wald intervals
# from the finite-difference observed information.

const _HOLLANDERS2023_FIT = Ref{Union{Nothing,CTMultistateFit}}(nothing)

"""
    hollanders2023_fit() -> CTMultistateFit

The reduced maximum-likelihood refit ([`CTMultistateProblem`](@ref) with default options,
started from [`hollanders2023_start`](@ref)), computed once per session and cached. It takes
several minutes; set `JULIA_NUM_THREADS` to use more than one thread.
"""
function hollanders2023_fit()
    _HOLLANDERS2023_FIT[] === nothing || return _HOLLANDERS2023_FIT[]
    _HOLLANDERS2023_FIT[] = solve(CTMultistateProblem())
    return _HOLLANDERS2023_FIT[]
end

# largest |ours − printed| over the median and 95% HDI bounds of the status and load effects
function _hollanders2023_slope_deviation()
    s = hollanders2023_posterior()
    printed = ((s.phi_status, (-0.167, -1.714, 1.053)), (s.phi_load, (1.953, -0.011, 3.834)))
    return maximum(maximum(abs.((x.median, x.lower, x.upper) .- p)) for (x, p) in printed)
end

# P(β ≤ 0) under the Wald (normal) approximation of the refit
function _hollanders2023_wald_below_zero(name::Symbol)
    fit = hollanders2023_fit()
    j = findfirst(==(name), fit.problem.names)
    return cdf(Normal(fit.θ[j], sqrt(fit.vcov[j, j])), 0.0)
end
