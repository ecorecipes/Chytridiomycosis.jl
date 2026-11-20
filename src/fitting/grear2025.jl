# Posterior-mode fit of the Grear et al. (2025) hurdle model (doi:10.1002/ecm.70001) to the
# bundled survey data, and the temperatures at which its quadratic temperature effects peak.
# The paper summarised a full NUTS posterior (brms/Stan, 30 imputed data sets); a mode is a
# quick check of the direction and shape of the effects, not a reproduction of the printed
# posterior summaries (which we could not read). See src/models/grear2025.jl.

"""
    grear2025_fit(prob = Grear2025Problem(); iters = 5000) -> NamedTuple

Posterior mode ([`posterior_mode`](@ref)) of a [`Grear2025Problem`](@ref), started from the
intercept-only moments. Returns `problem`, `x`, `logdensity`, `gradient`, `iterations` and
`coefficients` ([`hurdle_coefficients`](@ref): load-part `b`, hurdle-part `b_hu`, prevalence
`b_prevalence = −b_hu`, on the uncentred covariates [`Grear2025Problem`](@ref) lists).
"""
function grear2025_fit(prob::Grear2025Problem = Grear2025Problem(); iters::Integer = 5000)
    m = prob.model
    x0 = zeros(LogDensityProblems.dimension(m))
    ly = m.logy[m.positive]
    x0[m.blocks.a] = sum(ly) / length(ly)
    fz = 1 - count(m.positive) / length(m.y)
    x0[m.blocks.a_hu] = log(fz / (1 - fz))
    x0[m.blocks.log_sigma] = log(sqrt(sum(abs2, ly .- x0[m.blocks.a]) / length(ly)))
    for g in m.blocks.groups
        x0[g.log_sd] = x0[g.log_sd_hu] = log(0.5)
    end
    r = posterior_mode(m; x0, iters)
    return (problem = prob, x = r.x, logdensity = r.logdensity, gradient = r.gradient,
            iterations = r.iterations, coefficients = hurdle_coefficients(m, r.x))
end

"""
    grear2025_temperature_optimum(fit) -> NamedTuple

Temperature (30-day mean daily maximum, °C) at which the fitted quadratic temperature effect
peaks, `−b₁/(2b₂)`, for the prevalence (`prevalence`) and for the log load of infected
animals (`intensity`), with the quadratic coefficients (`curvature_prevalence`,
`curvature_intensity`; negative means a peak). The abstract reports "peak infection at
intermediate temperatures" for both.
"""
function grear2025_temperature_optimum(fit::NamedTuple)
    c = fit.coefficients
    i1 = findfirst(==("m30max"), c.covariates)
    i2 = findfirst(==("m30max^2"), c.covariates)
    bp, bl = c.b_prevalence, c.b
    return (prevalence = -bp[i1] / (2bp[i2]), intensity = -bl[i1] / (2bl[i2]),
            curvature_prevalence = bp[i2], curvature_intensity = bl[i2])
end

const _GREAR_CACHE = Dict{Any,Any}()

# The fit used by the validation targets and tests: species random intercepts only
# (independent; the phylogeny is not bundled), 1000 L-BFGS iterations. This is an
# approximate mode: species with no positive swab push their hurdle effects and the species
# standard deviations outwards (no finite mode is reached quickly), but the fixed temperature
# effects settle within a few hundred iterations (optimum 14.9 °C for prevalence and 10.9 °C
# for intensity after 1000 and after 5000 iterations). Cached per session.
function _grear2025_reference_fit()
    get!(_GREAR_CACHE, :reference) do
        grear2025_fit(Grear2025Problem(; random = (:species,)); iters = 1000)
    end
end
