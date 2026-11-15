# Model selection for Palomar et al. (2023) Biodivers Conserv (doi:10.1007/s10531-022-02525-3):
# the 13 multi-event models of Table 1 ranked by QAIC with the U-CARE overdispersion
# ĉ = 1.474. The capture histories were not deposited, so this runs on any histories coded as
# in palomar2023_structure (the tests use simulated ones); nothing here reproduces a
# published number until the data are available.

"""
    palomar2023_model_selection(histories; ĉ = 1.474, models = palomar2023_models())
        -> Vector{NamedTuple}

Fit every model of Palomar et al. (2023, Table 1) to `histories` (an `n×T` matrix of event
codes 0–3, coded as in [`palomar2023_structure`](@ref)) by maximum likelihood and rank them by
QAIC (`deviance/ĉ + 2k`, with `k` the number of identifiable parameters). Returns one named
tuple per model, sorted by QAIC: `model`, `fit` (a [`MultiEventFit`](@ref)), `qaic`,
`delta_qaic` and the published `published_qaic`.
"""
function palomar2023_model_selection(histories::AbstractMatrix; ĉ::Real = 1.474,
                                     models = palomar2023_models())
    rows = map(models) do m
        fit = solve(MultiEventProblem(histories, m.structure))
        (model = m.model, fit = fit, qaic = qaic(fit, ĉ), published_qaic = m.qaic)
    end
    best = minimum(r.qaic for r in rows)
    out = [merge(r, (delta_qaic = r.qaic - best,)) for r in rows]
    return sort(out; by = r -> r.qaic)
end

"""
    palomar2023_swabbing_experiment(; rng = Random.Xoshiro(2023), survival = 0.8,
        infection = 0.2, clearance = 0.6, recapture = 0.6, swab_fraction = 0.43,
        unswabbed = [1, 2, 12, 13], initial_infected = 0.28,
        releases = [100, 80, 70, 60, 55, 50, 45, 40, 35, 30, 25, 22, 20],
        structure = palomar2023_structure()) -> NamedTuple

A simulation experiment (our analysis, not a published result) on the design of Palomar et
al. (2023). Capture histories of 632 males over 13 annual occasions are simulated with the
*same* survival in both states, constant infection and clearance probabilities and
recapture, and state assessed by swabbing a random fraction `swab_fraction` of captures,
independent of state, except in the `unswabbed` years (2006, 2007, 2017 and 2018 were not
swabbed; swabs were taken 2008–2016), when every capture is "status unknown". The same
histories are then fitted with `structure` (default: the published model 1) twice: as
published, with a state-dependent δ constant over years (`published`), and with the
unswabbed years declared unassessed (`unassessed`). Each fit is summarised by the mean over
the first `T − 2` intervals (the last is confounded) of `S_I`, `S_U`, `ψ_U→I`, `ψ_I→U`, and
its `δ_I`, `δ_U`; `truth` holds the simulated values and `histories` the data.
"""
function palomar2023_swabbing_experiment(; rng::AbstractRNG = Random.Xoshiro(2023),
        survival::Real = 0.8, infection::Real = 0.2, clearance::Real = 0.6,
        recapture::Real = 0.6, swab_fraction::Real = 0.43,
        unswabbed::AbstractVector{<:Integer} = [1, 2, 12, 13], initial_infected::Real = 0.28,
        releases::AbstractVector{<:Integer} = [100, 80, 70, 60, 55, 50, 45, 40, 35, 30, 25, 22, 20],
        structure::MultiEventStructure = palomar2023_structure())
    T = length(releases)
    Ψ = repeat([1-clearance clearance; infection 1-infection], 1, 1, T - 1)
    truth = (initial = [initial_infected, 1 - initial_infected],
             survival = fill(Float64(survival), 2, T - 1), transition = Ψ,
             detection = fill(Float64(recapture), 2, T - 1),
             assignment = [swab_fraction 0 1-swab_fraction; 0 swab_fraction 1-swab_fraction])
    δ0 = _ms_unassessed(structure)
    δt = [t in unswabbed ? δ0 : truth.assignment for t in 1:T]
    H = simulate_capture_histories(rng, structure, truth, releases; assignment = δt)
    summ(fit) = (S_I = mean(fit.parameters.survival[1, 1:T-2]),
                 S_U = mean(fit.parameters.survival[2, 1:T-2]),
                 infection = mean(fit.parameters.transition[2, 1, 1:T-2]),
                 clearance = mean(fit.parameters.transition[1, 2, 1:T-2]),
                 δ_I = fit.parameters.assignment[1, 1], δ_U = fit.parameters.assignment[2, 2],
                 loglik = fit.loglik)
    published = solve(MultiEventProblem(H, structure))
    corrected = solve(MultiEventProblem(H, structure; unassessed = unswabbed))
    return (truth = (S_I = survival, S_U = survival, infection = infection,
                     clearance = clearance, δ_I = swab_fraction, δ_U = swab_fraction),
            published = summ(published), unassessed = summ(corrected), histories = H)
end
