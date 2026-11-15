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
    hollanders2023_refit_path()

Path of the stored reduced refit, `data/hollanders2023/refit.toml`, written by
`scripts/hollanders2023_refit.jl`: the parameter names, the optimum `theta`, its Wald
covariance `vcov`, the log-likelihood, the gradient at the optimum and the BFGS iterations.
These are our fit, not published results.
"""
hollanders2023_refit_path() = joinpath(pkgdir(@__MODULE__), "data", "hollanders2023", "refit.toml")

"""
    hollanders2023_fit(; refit = false) -> CTMultistateFit

The reduced maximum-likelihood refit ([`CTMultistateProblem`](@ref) with default options,
started from [`hollanders2023_start`](@ref)). By default it is read from the stored result
([`hollanders2023_refit_path`](@ref)), so the validation targets do not re-optimise; the
log-likelihood at the stored optimum is checked against the model in the tests. With
`refit = true`, or if the file is missing, the fit is recomputed (several minutes; set
`JULIA_NUM_THREADS` to use more than one thread). The result is cached for the session.
"""
function hollanders2023_fit(; refit::Bool = false)
    refit && return (_HOLLANDERS2023_FIT[] = solve(CTMultistateProblem()))
    _HOLLANDERS2023_FIT[] === nothing || return _HOLLANDERS2023_FIT[]
    path = hollanders2023_refit_path()
    _HOLLANDERS2023_FIT[] = isfile(path) ? _hollanders2023_stored_fit(path) :
                            solve(CTMultistateProblem())
    return _HOLLANDERS2023_FIT[]
end

function _hollanders2023_stored_fit(path::AbstractString)
    d = TOML.parsefile(path)
    prob = CTMultistateProblem()
    Symbol.(d["names"]) == prob.names ||
        throw(ArgumentError("$path was written for different parameters; rerun scripts/hollanders2023_refit.jl"))
    θ = Float64.(d["theta"])
    V = permutedims(reduce(hcat, [Float64.(r) for r in d["vcov"]]))
    return CTMultistateFit(prob, θ, V, Float64(d["loglik"]), Float64.(d["gradient"]),
                           Int(d["iterations"]))
end

# Write a fit in the format _hollanders2023_stored_fit reads (scripts/hollanders2023_refit.jl).
function _hollanders2023_write_fit(path::AbstractString, fit::CTMultistateFit; comment = "")
    d = Dict{String,Any}("names" => string.(fit.problem.names), "theta" => fit.θ,
        "vcov" => [collect(r) for r in eachrow(fit.vcov)], "loglik" => fit.loglik,
        "gradient" => fit.gradient, "iterations" => fit.iterations)
    open(path, "w") do io
        isempty(comment) || println(io, "# ", replace(comment, "\n" => "\n# "))
        TOML.print(io, d)
    end
    return path
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
