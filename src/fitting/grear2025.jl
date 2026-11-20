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

# The fit used by the validation targets and tests: the printed structure (species, site and
# region random intercepts; independent species, as the phylogeny is not bundled), one
# imputed data set, 1000 L-BFGS iterations. This is NOT a converged mode: sites and species
# whose swabs are all negative (or all positive) separate, and their random effects and the
# standard deviations keep growing (site SD of the hurdle part > 70 on the logit scale after
# 3000 iterations). The prevalence fixed effects are stable from 500 to 3000 iterations
# (temperature optimum 14.6 °C throughout); the intensity ones are not (optimum 11.9, 14.0 and
# 15.9 °C after 500, 1000 and 3000 iterations), so only prevalence is compared with numbers.
# Cached per session.
function _grear2025_reference_fit()
    get!(_GREAR_CACHE, :reference) do
        grear2025_fit(Grear2025Problem(); iters = 1000)
    end
end

const _GREAR2025_PRINTED = Dict("AgeJuvenile" => "5_juvenile", "AgeLarval" => "5_larval",
                                "SexMale" => "4_male", "SexUnknown" => "4_undetermined",
                                "m30max" => "1_temp", "m30max^2" => "2_temp2",
                                "t10precip" => "3_precip")

"""
    grear2025_printed_intervals() -> Dict{Symbol,Tuple{Float64,Float64}}

The 95% credible intervals of the fixed effects printed by Grear et al. (2025, Results and
Figs 6–7), read from the compendium (`data/compendium/grear2025.toml`): keys `:alpha1_temp`,
…, `:beta5_larval` (α prevalence, β log intensity; 1 temperature, 2 temperature², 3
precipitation, 4 sex, 5 life stage). No point estimates are printed.
"""
function grear2025_printed_intervals()
    t = TOML.parsefile(joinpath(COMPENDIUM_DIR, "grear2025.toml"))
    out = Dict{Symbol,Tuple{Float64,Float64}}()
    for r in t["parameter"]
        startswith(r["symbol"], "alpha") || startswith(r["symbol"], "beta") || continue
        u = r["uncertainty"]
        out[Symbol(r["symbol"])] = (Float64(u["lo"]), Float64(u["hi"]))
    end
    return out
end

"""
    grear2025_coefficient_check(fit) -> Vector{NamedTuple}

Compare the fixed effects of a [`grear2025_fit`](@ref) with the printed 95% intervals
([`grear2025_printed_intervals`](@ref)): one row per coefficient with `symbol`, `value`
(ours, prevalence coefficients on the logit-π scale, intensity on the log scale), `lo`, `hi`
and `inside`.
"""
function grear2025_coefficient_check(fit::NamedTuple)
    c = fit.coefficients
    ci = grear2025_printed_intervals()
    rows = NamedTuple[]
    for (part, v) in (("alpha", c.b_prevalence), ("beta", c.b))
        for (j, nm) in enumerate(c.covariates)
            s = Symbol(part, _GREAR2025_PRINTED[nm])
            lo, hi = ci[s]
            push!(rows, (symbol = s, value = v[j], lo = lo, hi = hi, inside = lo <= v[j] <= hi))
        end
    end
    return rows
end
