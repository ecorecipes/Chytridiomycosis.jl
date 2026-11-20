# Refits of Brannelly et al. (2018) BMC Ecol 18:34 (doi:10.1186/s12898-018-0189-5) to the
# bundled capture histories: the ten conditional Arnason–Schwarz models of Table 2, their
# model-averaged monthly survival, and the POPAN models with the super-population size.
# The paper used M-SURGE and Program MARK; these are maximum-likelihood refits with the
# package's layer (src/observation/cmr_open.jl), so differences in data processing (see
# src/models/brannelly2018.jl) show up as small differences in deviance.

"""
    brannelly2018_fit(model = 1; data = load_brannelly2018()) -> ArnasonSchwarzFit

Fit model `model` of [`brannelly2018_models`](@ref) (an index, or a
[`MultiEventStructure`](@ref)) to the capture histories by maximum likelihood with the
interval lengths of `data`.
"""
function brannelly2018_fit(model = 1; data::Brannelly2018Data = load_brannelly2018())
    st = model isa MultiEventStructure ? model : brannelly2018_models()[model].structure
    return solve(ArnasonSchwarzProblem(data.histories, st; intervals = data.intervals))
end

const _BRANNELLY_CACHE = Dict{Any,Any}()

"""
    brannelly2018_model_selection(; data = load_brannelly2018()) -> Vector{NamedTuple}

Fit the ten models of Table 2 and return, in printed order, `name`, `fit` (an
[`ArnasonSchwarzFit`](@ref)), `deviance`, `k` (the printed number of parameters, used for
AIC as in the paper), `aic` (`deviance + 2k`), `delta_aic` and `weight` (Akaike weight).
Results are cached per `data`.
"""
function brannelly2018_model_selection(; data::Brannelly2018Data = load_brannelly2018())
    key = (:cas, data.occasions, data.interval_unit)
    haskey(_BRANNELLY_CACHE, key) && return _BRANNELLY_CACHE[key]
    rows = map(brannelly2018_models()) do m
        fit = brannelly2018_fit(m.structure; data)
        (name = m.name, fit = fit, deviance = fit.deviance, k = m.parameters,
         aic = fit.deviance + 2m.parameters)
    end
    best = minimum(r.aic for r in rows)
    w = akaike_weights([r.aic for r in rows])
    out = [merge(r, (delta_aic = r.aic - best, weight = w[i])) for (i, r) in enumerate(rows)]
    _BRANNELLY_CACHE[key] = out
    return out
end

"""
    brannelly2018_survival(; data = load_brannelly2018(), models = 1:3) -> NamedTuple

Model-averaged monthly survival of *Bd*-negative (`uninfected`) and *Bd*-positive
(`infected`) frogs over `models` (default: the three models within 2 AIC of the best in
Table 2, which the paper averaged), with Akaike weights renormalised over those models
("natural averaging", Burnham & Anderson 2002). Survival is constant over time in these
models. Also returns the `weights`.
"""
function brannelly2018_survival(; data::Brannelly2018Data = load_brannelly2018(),
                                models = 1:3)
    sel = brannelly2018_model_selection(; data)[models]
    w = akaike_weights([r.aic for r in sel])
    S = [r.fit.parameters.survival[:, 1] for r in sel]
    return (uninfected = sum(w[i] * S[i][1] for i in eachindex(w)),
            infected = sum(w[i] * S[i][2] for i in eachindex(w)), weights = w)
end

"""
    brannelly2018_popan(model = 1; data = load_brannelly2018()) -> PopanFit

Fit POPAN model `model` of [`brannelly2018_popan_models`](@ref) (default: the best,
Φ(.)p(t)pent(t)) to the capture histories, ignoring disease state. The fit's `N` and `N_ci`
are the super-population size and its log-based 95% interval. Cached per `data`.
"""
function brannelly2018_popan(model::Integer = 1; data::Brannelly2018Data = load_brannelly2018())
    key = (:popan, model, data.occasions, data.interval_unit)
    haskey(_BRANNELLY_CACHE, key) && return _BRANNELLY_CACHE[key]
    m = brannelly2018_popan_models()[model]
    fit = solve(PopanProblem(data.histories; survival = m.survival, detection = m.detection,
                             entry = m.entry, intervals = data.intervals))
    _BRANNELLY_CACHE[key] = fit
    return fit
end
