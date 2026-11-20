# Mosher, Bailey, Muths & Huyvaert (2018) Host–pathogen metapopulation dynamics suggest high
# elevation refugia for boreal toads. Ecol Appl 28:926-937, doi:10.1002/eap.1699.
# Data: Dryad doi:10.5061/dryad.1c990 (CC0): detection histories of boreal toads (Anaxyrus
# boreas boreas) and Bd at 83 historic breeding sites in the southern Rocky Mountains,
# 2001-2010, with site and survey covariates.
#
# SOURCE OF THE EQUATIONS. The published article (paywalled) could not be read. The model
# below was read from the first author's PhD dissertation (Mosher BA 2017, "Informing the
# ecology and conservation of amphibians imperiled by chytridiomycosis", Colorado State
# University; Mountain Scholar), whose Chapter 2 ("Complex ecological relationships influence
# extinction risk in an amphibian-pathogen system") is the pre-publication version of the
# paper and whose Chapter 1 (published as Mosher et al. 2018, Ecography 41:32) gives the
# transition and detection matrices (Tables 1.1, 1.2). The dissertation analysed 82 sites;
# the Dryad deposit has 83. Numbers in the published paper may differ from the dissertation.
#
# Model: dynamic two-species occupancy, conditional parameterisation of Richmond et al.
# (2010), species A = boreal toad (host), species B = Bd (pathogen), states (A, B, U, AB);
# the equations are in src/observation/occupancy.jl. Bd is detected only on swabs of
# detected toads, so pB = rBa = 0 and the Bd-only state B is unobservable. Colonisation and
# extinction are constant over years (time-varying rates are not identifiable with the
# unobservable state; dissertation Ch. 1). Initial state probabilities are constant over
# sites. Covariates enter on the logit scale: elevation (linear or quadratic), presence of
# other amphibians, habitat permanence and region for the vital rates; day of year (linear
# or quadratic) for toad detection pA, rA; elevation for Bd detection rBA.
#
# Choices where the source is silent (src/fitting/mosher2018.jl):
# - Covariates are standardised (mean 0, SD 1 over sites, or over surveys for day of year).
#   With an intercept this does not change the maximised likelihood.
# - When both members of a pair (γA, γAB), (γB, γBA), (εA, εAB), (εB, εBA) or (pA, rA)
#   have the same covariate structure, they share slopes and differ in intercept
#   ("additive"). Inferred from the parameter counts K of dissertation Table 2.3, which this
#   rule reproduces for every model listed except the null model (14 printed, 13 here).
# - AICc uses the number of sites as the sample size (the Program MARK default for
#   occupancy data).
# - A performed survey without a recorded day of year is given the mean day of year.

"""
    Mosher2018Parameters(; ψA, ψBA, ψBa, γA, γAB, γB, γBA, εA, εAB, εB, εBA, pA, rA, rBA,
                           pB = 0, rBa = 0)
    Mosher2018Parameters(:mosher2018; set = "")

Parameters of the dynamic two-species (host–pathogen) occupancy model of Mosher et al.
(2018), species A the boreal toad and B Bd, all probabilities per year (vital rates) or per
survey (detection):

- `ψA`, `ψBA`, `ψBa`: initial occupancy of toads, and of Bd with toads present or absent;
- `γA`, `γAB`: toad colonisation with Bd absent or present; `εA`, `εAB`: toad extinction
  with Bd absent or present;
- `γB`, `γBA`: Bd colonisation with toads absent or present; `εB`, `εBA`: Bd extinction with
  toads absent or present;
- `pA`, `rA`: toad detection with Bd absent or present; `rBA`: Bd detection on a detected
  toad; `pB`, `rBa`: Bd detection with toads absent or undetected (0 for swab data).

The keyed constructor reads the scenarios of Mosher (2017, PhD dissertation, Table 1.3) from
the compendium: the base set is the "SRM boreal toad–Bd" scenario chosen to mimic the
southern Rocky Mountain data, and `set = "rare_susceptible"` the rare, susceptible host.
These are simulation scenarios, not estimates; the published paper's estimates were not
available (see the docs page).
"""
Base.@kwdef struct Mosher2018Parameters{T<:Real}
    ψA::T
    ψBA::T
    ψBa::T
    γA::T
    γAB::T
    γB::T
    γBA::T
    εA::T
    εAB::T
    εB::T
    εBA::T
    pA::T
    rA::T
    rBA::T
    pB::T = zero(pA)
    rBa::T = zero(pA)
end

const MOSHER2018_FIELDS = (:ψA, :ψBA, :ψBa, :γA, :γAB, :γB, :γBA, :εA, :εAB, :εB, :εBA, :pA,
                           :rA, :rBA, :pB, :rBa)

function Mosher2018Parameters(key::Symbol; set::AbstractString = "")
    key === :mosher2018 || throw(ArgumentError("unknown source :$key"))
    q = parameter_set(:mosher2018; set)
    vals = map(f -> float(getproperty(q, f)), MOSHER2018_FIELDS)
    return Mosher2018Parameters(vals...)
end

with(p::Mosher2018Parameters; kw...) =
    Mosher2018Parameters(; (f => get(kw, f, getfield(p, f)) for f in MOSHER2018_FIELDS)...)

"""
    initial_state(p::Mosher2018Parameters) -> Vector

Initial state distribution `(A, B, U, AB)`, from [`two_species_initial`](@ref).
"""
initial_state(p::Mosher2018Parameters) = two_species_initial(p.ψA, p.ψBA, p.ψBa)

"""
    transition_matrix(p::Mosher2018Parameters) -> Matrix

Annual 4×4 transition matrix, states `(A, B, U, AB)`, from
[`two_species_transition`](@ref).
"""
transition_matrix(p::Mosher2018Parameters) =
    two_species_transition(; γA = p.γA, γAB = p.γAB, γB = p.γB, γBA = p.γBA, εA = p.εA,
                           εAB = p.εAB, εB = p.εB, εBA = p.εBA)

"""
    detection_matrix(p::Mosher2018Parameters) -> Matrix

Per-survey 4×8 observation matrix, from [`two_species_detection`](@ref).
"""
detection_matrix(p::Mosher2018Parameters) =
    two_species_detection(; pA = p.pA, rA = p.rA, rBA = p.rBA, pB = p.pB, rBa = p.rBa)

"""
    Mosher2018Problem(p::Mosher2018Parameters; nseasons = 10)

Projection of the expected proportions of sites in each state of the Mosher et al. (2018)
model over `nseasons` years from the initial distribution. `solve(prob)` returns a
`NamedTuple` with `season` (1, 2, …), `states` (`nseasons×4`, columns `A`, `B`, `U`, `AB`),
`host` (`A + AB`, toad occupancy), `pathogen` (`B + AB`), `host_only` (`A`) and
`host_without_pathogen` (`A / (A + AB)`, the share of toad sites without Bd).
"""
struct Mosher2018Problem{P<:Mosher2018Parameters}
    p::P
    nseasons::Int
end
Mosher2018Problem(p::Mosher2018Parameters; nseasons::Integer = 10) =
    (nseasons >= 1 || throw(ArgumentError("nseasons must be positive"));
     Mosher2018Problem(p, Int(nseasons)))

function CommonSolve.solve(prob::Mosher2018Problem)
    Φ = transition_matrix(prob.p)
    x = initial_state(prob.p)
    S = zeros(eltype(x), prob.nseasons, 4)
    for t in 1:prob.nseasons
        t > 1 && (x = transpose(Φ) * x)
        S[t, :] .= x
    end
    host = S[:, 1] .+ S[:, 4]
    return (season = collect(1:prob.nseasons), states = S, host = host,
            pathogen = S[:, 2] .+ S[:, 4], host_only = S[:, 1],
            host_without_pathogen = S[:, 1] ./ host)
end

"""
    simulate_mosher2018(rng, p::Mosher2018Parameters; nsites, nseasons, nsurveys,
                        elevation = zeros(nsites)) -> Mosher2018Data

Simulate toad and Bd detection histories of `nsites` sites over `nseasons` years with
`nsurveys` surveys each from constant parameters `p`. Bd is recorded (`0` or `1`) on every
survey; `elevation` is stored as a covariate but does not affect the simulation.
"""
function simulate_mosher2018(rng::AbstractRNG, p::Mosher2018Parameters; nsites::Integer,
                             nseasons::Integer, nsurveys::Integer,
                             elevation = zeros(nsites))
    toad = Array{Union{Missing,Int}}(undef, nsites, nseasons, nsurveys)
    bd = similar(toad)
    φ0, Φ, D = initial_state(p), transition_matrix(p), detection_matrix(p)
    for i in 1:nsites
        s = simulate_dynamic_occupancy(rng, φ0, Φ, D, nseasons, nsurveys; observed = 1:4)
        for t in 1:nseasons, j in 1:nsurveys
            toad[i, t, j], bd[i, t, j] = two_species_records(s.codes[t, j])
        end
    end
    return Mosher2018Data(; toad, bd, elevation = collect(float.(elevation)))
end
