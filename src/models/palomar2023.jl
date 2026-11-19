# Palomar G, Fernández-Chacón A, Bosch J (2023) Amphibian survival compromised by long-term
# effects of chytrid fungus. Biodiversity and Conservation 32:793–809.
# https://doi.org/10.1007/s10531-022-02525-3
# This implements a published model; it is
# not a new model.
#
# Read from the paper (CC BY) and Online Appendices A (E-SURGE model pattern) and B (U-CARE
# goodness of fit), lab/papers/palomar2023{,_si1,_si2}. The data statement says "Data will be
# available online upon publication"; no deposit was found (Dryad, Zenodo, DataCite, figshare,
# DIGITAL.CSIC searched 2026-11-15), so nothing is refitted and the estimates are pending
# validation targets.
#
# Model. Multi-event Arnason–Schwarz model (Pradel 2005) of male Bufo spinosus at Laguna
# Grande, Peñalara, 13 annual occasions 2006–2018. States I (alive Bd+), U (alive Bd−),
# † (dead). Steps: initial state (π), survival S, state transition ψ (conditional on
# survival), recapture p (B1), state assignment δ (B2). Events: 0 not seen, 1 seen infected,
# 2 seen uninfected, 3 seen with status unknown. Retained model (Table 1, model 1): S state +
# time, ψ state × time, p time, δ state; QAIC with ĉ = 1.474 (U-CARE Test 3G + Test M).
#
# Printed versus appendix. The main text codes events as 1 = "seen alive and Bd-negative",
# 2 = "seen alive and Bd-positive", but its printed B2 sends I to event 1 and U to event 2,
# and its next paragraph says a seen Bd+ toad assessed becomes "event 1"; Appendix A codes
# 1 = infected and 2 = not infected, consistent with B2. The appendix coding is used here.
# Swabbing was a random subsample of captures (588 of 1369, 2008–2016), yet δ is state-
# dependent, and the reported δ_U = 1.000 (no standard error, a boundary) with δ_I = 0.119
# means every capture of unknown status is attributed to state I by the fitted model.

"""
    palomar2023_structure(; survival = :state_plus_time, transition = :state_time,
                          detection = :time, state_assignment = :state) -> MultiEventStructure

The multi-event model pattern of Palomar et al. (2023) (Online Appendix A) for a
[`MultiEventProblem`](@ref): alive states 1 = I (*Bd*+) and 2 = U (*Bd*−); events 0 not
seen, 1 seen infected, 2 seen uninfected, 3 seen with status unknown (the appendix coding);
an encountered I gives event 1 with probability `δ_I`, otherwise 3; an encountered U gives
event 2 with probability `δ_U`, otherwise 3. The defaults are the retained model 1 of
Table 1 ("S: state + time, ψ: state * time, p: time, δ: state"); the other rows of the
table are [`palomar2023_models`](@ref). `state_assignment = :constant` makes δ the same in
both states, as random subsampling of captures for swabbing implies.
"""
palomar2023_structure(; survival::Symbol = :state_plus_time, transition::Symbol = :state_time,
                      detection::Symbol = :time, state_assignment::Symbol = :state) =
    MultiEventStructure(; nstates = 2, assignment = Bool[1 0 1; 0 1 1], baseline = [3, 3],
                        survival, transition, detection, state_assignment)

"""
    palomar2023_models() -> Vector{NamedTuple}

The 13 multi-event models of Palomar et al. (2023, Table 1), in the printed order: `model`
(number), `structure` (a [`MultiEventStructure`](@ref) from [`palomar2023_structure`](@ref))
and the published `qaic`. Fit each to the capture histories with
[`palomar2023_model_selection`](@ref).
"""
function palomar2023_models()
    rows = [(1, :state_plus_time, :state_time, :time, 4293.7),
            (2, :time, :state_time, :time, 4297.0),
            (3, :state_time, :state_time, :time, 4300.2),
            (4, :state_time, :state_time, :state_time, 4305.9),
            (5, :state_time, :state_time, :state_plus_time, 4311.2),
            (6, :state, :state_time, :time, 4325.5),
            (7, :state_time, :state_time, :state, 4331.7),
            (8, :constant, :state_time, :time, 4335.5),
            (9, :state_time, :state_plus_time, :time, 4341.1),
            (10, :state_time, :state_time, :constant, 4344.8),
            (11, :state_time, :time, :time, 4381.0),
            (12, :state_time, :state, :time, 4418.8),
            (13, :state_time, :constant, :time, 4430.2)]
    return [(model = r[1], structure = palomar2023_structure(; survival = r[2], transition = r[3],
                                                              detection = r[4]),
             qaic = r[5]) for r in rows]
end

"""
    palomar2023_parameters(; occasions = 13, initial_infected = 0.5) -> NamedTuple

Natural-scale parameters of the multi-event model at the mean annual estimates of Palomar et
al. (2023) (compendium `palomar2023`), constant over `occasions`, in the layout of
[`multievent_parameters`](@ref) (states I, U): survival `S_I = 0.754`, `S_U = 0.883`;
transitions `ψ_U→I = 0.583`, `ψ_I→U = 0.225`; recapture `p = 0.593`; state assignment
`δ_I = 0.119`, `δ_U = 1`. The initial-state probability at first capture is not reported:
`initial_infected` is an assumption, used only to simulate.
"""
function palomar2023_parameters(; occasions::Integer = 13, initial_infected::Real = 0.5)
    v = parameter_set(:palomar2023)
    T = occasions
    S = repeat([v.survival_infected, v.survival_uninfected], 1, T - 1)
    Ψ = zeros(2, 2, T - 1)
    for t in 1:T-1
        Ψ[:, :, t] = [1-v.clearance v.clearance; v.infection 1-v.infection]
    end
    p = fill(v.recapture_mean, 2, T - 1)
    δ = [v.assignment_infected 0 1-v.assignment_infected
         0 v.assignment_uninfected 1-v.assignment_uninfected]
    return (initial = [initial_infected, 1 - initial_infected], survival = S, transition = Ψ,
            detection = p, assignment = δ)
end

"""
    palomar2023_implied_prevalence() -> Float64

Bd prevalence among living male toads implied by the mean annual estimates of Palomar et al.
(2023): the infected share of the [`quasistationary_distribution`](@ref) of the survival and
transition matrix (`S_I = 0.754`, `S_U = 0.883`, `ψ_U→I = 0.583`, `ψ_I→U = 0.225`). Our
arithmetic, for comparison with the observed swab prevalence of 27.79%.
"""
function palomar2023_implied_prevalence()
    par = palomar2023_parameters(; occasions = 2)
    Φ = survival_transition(par.survival[:, 1], par.transition[:, :, 1])
    return quasistationary_distribution(Φ)[1]
end
