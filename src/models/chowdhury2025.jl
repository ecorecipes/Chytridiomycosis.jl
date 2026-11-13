# Chowdhury MM, Gray MJ, Carter ED, Woodhams DC, McCartney JA, Surles JG, Nelms KM, Gan H,
# Peace A (2025) Fungal pathogen transmission dynamics in North American salamanders:
# mathematical insights for disease management. Ecol Model 501:111028.
# doi:10.1016/j.ecolmodel.2025.111028
# Wiki: lab/wiki/models/bsal-chowdhury2025-stage-temperature.md
#
# Eastern newt (Notophthalmus viridescens) and Bsal at 6 and 14 °C. Larvae (never infected),
# terrestrial efts and aquatic adults, SEIR in efts and adults, and one zoospore pool per
# habitat. Time in days; host densities per m². Checked against the published PDF (§§2–6,
# Table 2, Table 3, Figs 1, 4, 5, 7) and the online supplement (Elsevier mmc1.pdf,
# Appendices A–H; lab/papers/chowdhury2025_si/mmc1.pdf). No code or data were deposited
# ("data will be made available on request").
#
# Model (§3.2.3; state order of Appendix F):
#   S_l' = λA − η_j S_l
#   S_j' = α η_j S_l − β_j S_j I_j/N_T − f_j(Z_j) S_j − (δ_j + η_a) S_j
#   E_j' = f_j(Z_j) S_j + β_j S_j I_j/N_T − (σ_j + δ_j + η_a) E_j
#   I_j' = σ_j E_j − (γ_j + δ_j + μ_j + η_a) I_j
#   R_j' = γ_j I_j − (δ_j + η_a) R_j
#   S_a' = η_a(1 − A/K) S_j − β_a(A) S_a I_a/N_A − f_a(Z_a) S_a − δ_a S_a
#   E_a' = f_a(Z_a) S_a + η_a(1 − A/K) E_j + β_a(A) S_a I_a/N_A − (σ_a + δ_a) E_a
#   I_a' = η_a(1 − A/K) I_j + σ_a E_a − (γ_a + δ_a + μ_a) I_a
#   R_a' = η_a(1 − A/K) R_j + γ_a I_a − δ_a R_a
#   Z_j' = ω_j I_j − ξ_j Z_j,     Z_a' = ω_a I_a − ξ_a Z_a
# with A = S_a + E_a + I_a + R_a, N_T = S_j + E_j + I_j + R_j, N_A = A + S_l (Table 2),
# β_a(A) = P_a c_a A/(χ_a + A) (eq. 3a), β_j = c_j P_j (eq. 3b) and
# f_i(Z) = ρ_i b_i Z/(κ_i + Z) (eq. 4).
#
# What the PDF does not give, and what does not agree (details in the compendium notes):
# - χ_a, the half-saturation density of the adult contact function, is never printed. At
#   14 °C three figures agree on χ_a ≈ 5.9 m⁻² (the Fig. 1 curve, β_a(0) in Fig. 7 and the
#   R_a column of Table 3). At 6 °C, Fig. 7 with Table 2's c_a = 464.77 gives χ_a ≈ 11.7, and
#   this reproduces Figs 4, 5 and 7. The Fig. 1 curve and the 6 °C R_a column of Table 3
#   instead imply c_a ≈ 729 day⁻¹ and χ_a ≈ 23.3 m⁻². `variant = :table3` uses those.
# - Table 3's R_j needs c_j = 8/(14 × 12) = 0.0476 day⁻¹ (8 contacts in 14 days among 12
#   efts, per capita) rather than the 0.571 of §2.2 and Table 2, and b_j = 0.005 (§6.2). With
#   these, 14 °C reproduces; at 6 °C the direct part reproduces but the environmental part
#   is about 3 times larger than the formula gives (unresolved).
# - §2.3 says terrestrial spores persist 25 days at 6 °C and 20 at 14 °C; Table 2 gives
#   ξ_j = 1/20 and 1/25. Appendix C says efts mature after 3.5 years; Table 2 says 3. Table 3
#   reproduces only with the Table 2 values, which are the default. `as_printed = true` uses
#   the text values.
# - The eft trajectories of Figs 4 and 5 (efts flat or declining, never recovering in the
#   10-year runs of Appendix E) are not reproduced with Table 2's λ = 0.82 eggs adult⁻¹ day⁻¹:
#   larvae then recruit several efts per m² within months. They are approximated when
#   larval recruitment is negligible (λ ≈ 0). The adult trajectories are reproduced either way.
#
# Units: Z_a enters f_a against κ_a in zoospores per mL, and Z_j against κ_j in zoospores
# per m², while Z' = ω I − ξ Z with I in hosts per m² and ω in zoospores per host per day.
# The paper equates the two without a volume or area (an implicit 1 mL of water and 1 m² of
# soil per m² of habitat). The package keeps the paper's numbers and reports Z in the units
# of κ.

const CHOWDHURY_TEMPERATURES = (6, 14)
const CHOWDHURY_VARIANTS = (:table2, :table3)
const CHOWDHURY_STATE = (:S_l, :S_j, :E_j, :I_j, :R_j, :S_a, :E_a, :I_a, :R_a, :Z_j, :Z_a)

"""
    ChowdhuryParameters(; λ, α, η_j, η_a, K, δ_j, δ_a, σ_j, σ_a, κ_j, κ_a, γ_j = 0,
                        γ_a = 0, μ_j, μ_a, ρ_j, ρ_a, ξ_j, ξ_a, ω_j, ω_a, P_j, P_a, b_j, b_a,
                        c_j, c_a, χ_a, temperature = NaN)
    ChowdhuryParameters(:chowdhury2025; temperature = 14, K = 2, variant = :table2,
                        as_printed = false, kwargs...)

Parameters of the Chowdhury et al. (2025) stage-structured Bsal model of the eastern newt
(Table 2). Rates are per day and densities per m².

- Life cycle: birth rate `λ` (eggs per adult per day), larval survival to metamorphosis
  `α`, larva-to-eft rate `η_j`, eft-to-adult rate `η_a`, adult carrying capacity `K`, and
  natural death rates `δ_j` (efts) and `δ_a` (adults). Larvae do not die except through `α`.
- Disease, eft (`_j`) and adult (`_a`): rate of leaving the exposed class `σ`, recovery `γ`
  (zero in the paper), Bsal-induced death `μ`, shedding `ω` (zoospores per infected host
  per day), spore decay `ξ`, encystment fraction `ρ`, spore contact rate `b`, ID50 `κ`
  (spores per m² for efts, per mL for adults) and probability of transmission per contact
  `P`.
- Contacts: `c_j`, the per capita eft contact rate (constant), and the adult Holling type II
  contact rate `c_a A/(χ_a + A)` with maximum `c_a` and half-saturation density `χ_a`.

The keyed constructor reads the compendium (`data/compendium/chowdhury2025.toml`) at
`temperature` 6 or 14 °C with carrying capacity `K` (the paper uses 2, 4, 8 and 16 m⁻², see
[`chowdhury_scenarios`](@ref)). Other keywords override single fields.

- `variant = :table2` (default) uses Table 2, with `χ_a` derived from Fig. 7 (`χ_a` is not
  printed). These values reproduce the simulations of Figs 4, 5 and 7.
- `variant = :table3` uses the values that reproduce Table 3's reproduction numbers:
  `c_j = 8/(14 × 12)` per capita, `b_j = 0.005` (§6.2) and, at 6 °C, the contact curve of
  Fig. 1 (`c_a ≈ 729`, `χ_a ≈ 23.3`) instead of Table 2's `c_a = 464.77`.
- `as_printed = true` uses the values printed in the text where they contradict Table 2:
  terrestrial spore persistence of 25 days at 6 °C and 20 days at 14 °C (§2.3; Table 2 has
  the reverse), and an eft stage of 3.5 years (Appendix C; Table 2 has 3).
"""
struct ChowdhuryParameters
    λ::Float64
    α::Float64
    η_j::Float64
    η_a::Float64
    K::Float64
    δ_j::Float64
    δ_a::Float64
    σ_j::Float64
    σ_a::Float64
    κ_j::Float64
    κ_a::Float64
    γ_j::Float64
    γ_a::Float64
    μ_j::Float64
    μ_a::Float64
    ρ_j::Float64
    ρ_a::Float64
    ξ_j::Float64
    ξ_a::Float64
    ω_j::Float64
    ω_a::Float64
    P_j::Float64
    P_a::Float64
    b_j::Float64
    b_a::Float64
    c_j::Float64
    c_a::Float64
    χ_a::Float64
    temperature::Float64
    function ChowdhuryParameters(λ, α, η_j, η_a, K, δ_j, δ_a, σ_j, σ_a, κ_j, κ_a, γ_j, γ_a,
            μ_j, μ_a, ρ_j, ρ_a, ξ_j, ξ_a, ω_j, ω_a, P_j, P_a, b_j, b_a, c_j, c_a, χ_a,
            temperature)
        rates = (λ, η_j, η_a, δ_j, δ_a, σ_j, σ_a, γ_j, γ_a, μ_j, μ_a, ξ_j, ξ_a, ω_j, ω_a,
                 b_j, b_a, c_j, c_a)
        all(>=(0), rates) || throw(ArgumentError("rates must be non-negative"))
        (K > 0 && κ_j > 0 && κ_a > 0 && χ_a >= 0) ||
            throw(ArgumentError("K, κ_j and κ_a must be positive and χ_a non-negative"))
        all(x -> 0 <= x <= 1, (α, ρ_j, ρ_a, P_j, P_a)) ||
            throw(ArgumentError("α, ρ and P must lie in [0, 1]"))
        return new(λ, α, η_j, η_a, K, δ_j, δ_a, σ_j, σ_a, κ_j, κ_a, γ_j, γ_a, μ_j, μ_a, ρ_j,
            ρ_a, ξ_j, ξ_a, ω_j, ω_a, P_j, P_a, b_j, b_a, c_j, c_a, χ_a, temperature)
    end
end

ChowdhuryParameters(; λ, α, η_j, η_a, K, δ_j, δ_a, σ_j, σ_a, κ_j, κ_a, γ_j = 0.0,
        γ_a = 0.0, μ_j, μ_a, ρ_j, ρ_a, ξ_j, ξ_a, ω_j, ω_a, P_j, P_a, b_j, b_a, c_j, c_a, χ_a,
        temperature = NaN) =
    ChowdhuryParameters(λ, α, η_j, η_a, K, δ_j, δ_a, σ_j, σ_a, κ_j, κ_a, γ_j, γ_a, μ_j, μ_a,
        ρ_j, ρ_a, ξ_j, ξ_a, ω_j, ω_a, P_j, P_a, b_j, b_a, c_j, c_a, χ_a, temperature)

function ChowdhuryParameters(key::Symbol; temperature::Real = 14, K::Real = 2.0,
        variant::Symbol = :table2, as_printed::Bool = false, kwargs...)
    T = round(Int, temperature)
    (T == temperature && T in CHOWDHURY_TEMPERATURES) ||
        throw(ArgumentError("temperature must be one of $CHOWDHURY_TEMPERATURES °C"))
    variant in CHOWDHURY_VARIANTS ||
        throw(ArgumentError("variant must be one of $CHOWDHURY_VARIANTS"))
    sets = ["", "$(T)C"]
    variant === :table3 && push!(sets, "$(T)C_table3")
    as_printed && push!(sets, "$(T)C_as_printed")
    recs = parameter_records(key)
    vals = Dict{Symbol,Float64}()
    for s in sets, r in recs
        r.set == s || continue
        vals[r.symbol] = _time_part(r.units) === nothing ? r.value : in_time_unit(r, :day)
    end
    val(k) = haskey(vals, k) ? vals[k] : throw(ArgumentError("$key has no $k at $(T) °C"))
    base = (; λ = val(:λ), α = val(:α), η_j = 1 / val(:larval_period),
        η_a = 1 / val(:eft_period), K = Float64(K), δ_j = val(:δ_j),
        δ_a = 1 / val(:adult_lifespan), σ_j = val(:σ_j), σ_a = val(:σ_a), κ_j = val(:κ_j),
        κ_a = val(:κ_a), γ_j = val(:γ_j), γ_a = val(:γ_a), μ_j = val(:μ_j), μ_a = val(:μ_a),
        ρ_j = val(:ρ_j), ρ_a = val(:ρ_a), ξ_j = val(:ξ_j), ξ_a = val(:ξ_a), ω_j = val(:ω_j),
        ω_a = val(:ω_a), P_j = val(:P_j), P_a = val(:P_a), b_j = val(:b_j), b_a = val(:b_a),
        c_j = val(:c_j), c_a = val(:c_a), χ_a = val(:χ_a), temperature = Float64(T))
    return ChowdhuryParameters(; merge(base, values(kwargs))...)
end

time_unit(::Type{<:ChowdhuryParameters}) = :day

"""Return a copy of `p` with fields replaced, e.g. `with(p; c_j = 0.1)`."""
with(p::ChowdhuryParameters; kw...) =
    ChowdhuryParameters(; (f => get(kw, f, getfield(p, f)) for f in fieldnames(ChowdhuryParameters))...)

"""
    chowdhury_scenarios() -> Vector{NamedTuple{(:adults, :efts)}}

The four initial densities (per m²) simulated by Chowdhury et al. (2025, §5.1, Figs 4, 5
and 7): adults 2, 4, 8 and 16, paired with efts 0.1, 0.2, 0.4 and 0.8. The adult carrying
capacity `K` equals the initial adult density in each run. Read from the compendium.
"""
function chowdhury_scenarios()
    q = parameter_set(:chowdhury2025)
    return [(; adults = getfield(q, Symbol("adults_$i")), efts = getfield(q, Symbol("efts_$i")))
            for i in 1:4]
end

# --- Transmission terms ---------------------------------------------------------------------

_default_adult_contact(p::ChowdhuryParameters) = HollingContact(p.P_a, p.c_a, p.χ_a)
_default_eft_contact(p::ChowdhuryParameters) = FrequencyDependent(p.c_j * p.P_j)

# Total rate of new adult infections by contact. The contact function sees the adults only;
# with `larval_dilution` the frequency term is divided by N_A = A + S_l instead of A, as in
# the paper (Table 2: N_A is adults plus larvae).
function _adult_direct(tf::TransmissionFunction, S, I, A, N_A, larval_dilution::Bool)
    (S > 0 && I > 0 && A > 0) || return 0.0
    r = infection_rate(tf, S, I, A)
    return larval_dilution ? r * A / N_A : r
end
_eft_direct(tf::TransmissionFunction, S, I, N) = (S > 0 && I > 0 && N > 0) ?
    infection_rate(tf, S, I, N) : 0.0

_spore_foi(ρ, b, κ, Z) = Z > 0 ? ρ * b * Z / (κ + Z) : 0.0

"""
    BsalStageTemperatureProblem(p::ChowdhuryParameters; adults = p.K, efts = adults/20,
                                exposed = 0.01, larvae = 0, tspan = (0.0, 360.0),
                                adult_contact = HollingContact(p.P_a, p.c_a, p.χ_a),
                                eft_contact = FrequencyDependent(p.c_j * p.P_j),
                                larval_dilution = true, u0 = nothing)

The stage-structured Bsal model of Chowdhury et al. (2025, §3.2.3), in days, with state
`[S_l, S_j, E_j, I_j, R_j, S_a, E_a, I_a, R_a, Z_j, Z_a]` (larvae; efts; adults; terrestrial
and aquatic zoospores). `solve(prob)` integrates the `ODEProblem` of
[`to_ode_problem`](@ref) with `Tsit5`.

The default start is the paper's (§5.1): `adults` and `efts` per m² with a fraction
`exposed` of each in the exposed class, no larvae and no zoospores. Give `u0` to start
elsewhere.

The direct transmission terms can be replaced by any [`TransmissionFunction`](@ref):

- `adult_contact` gives the adult term as `infection_rate(adult_contact, S_a, I_a, A)`,
  multiplied by `A/N_A` when `larval_dilution` is true (`N_A = A + S_l`). The default,
  `HollingContact(P_a, c_a, χ_a)` with the dilution, is the paper's
  `β_a(A) S_a I_a/N_A` with `β_a(A) = P_a c_a A/(χ_a + A)` (eq. 3a).
- `eft_contact` gives the eft term as `infection_rate(eft_contact, S_j, I_j, N_T)`. The
  default `FrequencyDependent(c_j P_j)` is the paper's `β_j S_j I_j/N_T` (eq. 3b).

For example `adult_contact = FrequencyDependent(0.23)` puts the Tompros et al. (2022)
estimate, frequency-dependent and fitted in 1-m² tanks at 14 °C, in place of the paper's
contact function (a new combination, not a published result).
"""
struct BsalStageTemperatureProblem{F<:TransmissionFunction,G<:TransmissionFunction}
    p::ChowdhuryParameters
    adult_contact::F
    eft_contact::G
    larval_dilution::Bool
    u0::Vector{Float64}
    tspan::Tuple{Float64,Float64}
end

function BsalStageTemperatureProblem(p::ChowdhuryParameters; adults::Real = p.K,
        efts::Real = adults / 20, exposed::Real = 0.01, larvae::Real = 0.0,
        tspan = (0.0, 360.0), adult_contact::TransmissionFunction = _default_adult_contact(p),
        eft_contact::TransmissionFunction = _default_eft_contact(p),
        larval_dilution::Bool = true, u0 = nothing)
    0 <= exposed <= 1 || throw(ArgumentError("exposed must lie in [0, 1]"))
    if u0 === nothing
        u0 = [larvae, (1 - exposed) * efts, exposed * efts, 0, 0,
              (1 - exposed) * adults, exposed * adults, 0, 0, 0, 0]
    end
    length(u0) == 11 || throw(ArgumentError("u0 must have the 11 states $CHOWDHURY_STATE"))
    all(>=(0), u0) || throw(ArgumentError("u0 must be non-negative"))
    return BsalStageTemperatureProblem(p, adult_contact, eft_contact, larval_dilution,
        collect(float.(u0)), (Float64(tspan[1]), Float64(tspan[2])))
end

"""
    stage_transmission_rates(prob::BsalStageTemperatureProblem, u) -> NamedTuple

The transmission quantities of Fig. 7 at state `u`:

- `β_a`: the adult direct transmission rate per susceptible-infectious pair,
  `transmission_coefficient(adult_contact, A)`, which is `P_a c_a A/(χ_a + A)` (eq. 3a)
  for the default contact function;
- `β_j`: `transmission_coefficient(eft_contact, N_T)`, `c_j P_j` by default (eq. 3b);
- `f_j`, `f_a`: the per capita rates of infection from zoospores,
  `ρ b Z/(κ + Z)` (eq. 4);
- `direct_j`, `direct_a`, `environmental_j`, `environmental_a`: the total rates of new
  infections (per m² per day) by each route.
"""
function stage_transmission_rates(prob::BsalStageTemperatureProblem, u)
    p = prob.p
    Sl, Sj, Ej, Ij, Rj, Sa, Ea, Ia, Ra, Zj, Za = u
    A = Sa + Ea + Ia + Ra
    N_T = Sj + Ej + Ij + Rj
    f_j = _spore_foi(p.ρ_j, p.b_j, p.κ_j, Zj)
    f_a = _spore_foi(p.ρ_a, p.b_a, p.κ_a, Za)
    return (; β_a = transmission_coefficient(prob.adult_contact, A),
        β_j = transmission_coefficient(prob.eft_contact, N_T), f_j, f_a,
        direct_j = _eft_direct(prob.eft_contact, Sj, Ij, N_T),
        direct_a = _adult_direct(prob.adult_contact, Sa, Ia, A, A + Sl, prob.larval_dilution),
        environmental_j = f_j * max(Sj, 0.0), environmental_a = f_a * max(Sa, 0.0))
end

function chowdhury_rhs!(du, u, prob, t)
    p = prob.p
    Sl, Sj, Ej, Ij, Rj, Sa, Ea, Ia, Ra, Zj, Za = u
    A = Sa + Ea + Ia + Ra
    N_T = Sj + Ej + Ij + Rj
    dj = _eft_direct(prob.eft_contact, Sj, Ij, N_T)
    da = _adult_direct(prob.adult_contact, Sa, Ia, A, A + Sl, prob.larval_dilution)
    ej = _spore_foi(p.ρ_j, p.b_j, p.κ_j, Zj) * Sj
    ea = _spore_foi(p.ρ_a, p.b_a, p.κ_a, Za) * Sa
    rec = p.η_a * (1 - A / p.K)                   # eq. (2): efts that join the adults
    out_j = p.δ_j + p.η_a
    du[1] = p.λ * A - p.η_j * Sl
    du[2] = p.α * p.η_j * Sl - dj - ej - out_j * Sj
    du[3] = ej + dj - (p.σ_j + out_j) * Ej
    du[4] = p.σ_j * Ej - (p.γ_j + p.μ_j + out_j) * Ij
    du[5] = p.γ_j * Ij - out_j * Rj
    du[6] = rec * Sj - da - ea - p.δ_a * Sa
    du[7] = ea + rec * Ej + da - (p.σ_a + p.δ_a) * Ea
    du[8] = rec * Ij + p.σ_a * Ea - (p.γ_a + p.δ_a + p.μ_a) * Ia
    du[9] = rec * Rj + p.γ_a * Ia - p.δ_a * Ra
    du[10] = p.ω_j * Ij - p.ξ_j * Zj
    du[11] = p.ω_a * Ia - p.ξ_a * Za
    return nothing
end

to_ode_problem(prob::BsalStageTemperatureProblem) =
    ODEProblem(chowdhury_rhs!, prob.u0, prob.tspan, prob)

CommonSolve.solve(prob::BsalStageTemperatureProblem, alg = Tsit5();
        reltol = 1e-8, abstol = 1e-10, kwargs...) =
    solve(to_ode_problem(prob), alg; reltol, abstol, kwargs...)

"""
    stage_totals(u) -> NamedTuple

Summaries of a state vector of [`BsalStageTemperatureProblem`](@ref): `larvae`, `efts`
(`N_T`), `adults` (`A`), `infected_efts` and `infected_adults` (the `I` classes, as in
Fig. 4), `exposed_efts`, `exposed_adults`, and the zoospore pools `Z_j` and `Z_a`.
"""
function stage_totals(u)
    Sl, Sj, Ej, Ij, Rj, Sa, Ea, Ia, Ra, Zj, Za = u
    return (; larvae = Sl, efts = Sj + Ej + Ij + Rj, adults = Sa + Ea + Ia + Ra,
        infected_efts = Ij, infected_adults = Ia, exposed_efts = Ej, exposed_adults = Ea,
        Z_j = Zj, Z_a = Za)
end

"""
    stage_decline_time(prob::BsalStageTemperatureProblem; stage = :adults, fraction = 0.8)

The first time at which the total of `stage` (`:adults` or `:efts`) has fallen by
`fraction` from its initial value, or `Inf` if it does not within `prob.tspan`. Linear
interpolation between daily outputs.
"""
function stage_decline_time(prob::BsalStageTemperatureProblem; stage::Symbol = :adults,
        fraction::Real = 0.8)
    stage in (:adults, :efts) || throw(ArgumentError("stage must be :adults or :efts"))
    0 < fraction < 1 || throw(ArgumentError("fraction must lie in (0, 1)"))
    sol = solve(prob; saveat = range(prob.tspan...; step = 1.0))
    x = [getfield(stage_totals(u), stage) for u in sol.u]
    level = (1 - fraction) * x[1]
    for i in 2:length(x)
        if x[i] <= level
            return sol.t[i-1] + (x[i-1] - level) / (x[i-1] - x[i]) * (sol.t[i] - sol.t[i-1])
        end
    end
    return Inf
end

# --- Disease-free state and reproduction numbers --------------------------------------------

"""
    stage_disease_free_state(p::ChowdhuryParameters) -> (; S_l, S_j, S_a)

The disease-free equilibrium of Appendix F (per m²):

```
S_a* = ((αλ − δ_a) η_a − δ_a δ_j) K / (α η_a λ),
S_l* = λ S_a*/η_j,      S_j* = α λ S_a*/(δ_j + η_a).
```

It exists (is positive) when `α λ η_a > δ_a (η_a + δ_j)`. With Table 2 it has
`S_a* ≈ 0.98 K`, `S_j* ≈ 15.9 K` and `S_l* ≈ 85 K`: far more efts and larvae than the
simulations start with (0.05 K efts, no larvae).
"""
function stage_disease_free_state(p::ChowdhuryParameters)
    X = (p.α * p.λ - p.δ_a) * p.η_a - p.δ_a * p.δ_j
    S_a = X * p.K / (p.α * p.η_a * p.λ)
    return (; S_l = p.λ * S_a / p.η_j, S_j = p.α * p.λ * S_a / (p.δ_j + p.η_a), S_a)
end

"""
    stage_next_generation_matrix(p::ChowdhuryParameters; adult_contact, eft_contact,
                                 larval_dilution = true) -> (; F, V, K, R0)

The next-generation matrix `K = F V⁻¹` at the disease-free state
([`stage_disease_free_state`](@ref)), with infected compartments
`(E_j, I_j, E_a, I_a, Z_j, Z_a)` as in Appendix F. Shedding into the pools and the
maturation of infected efts into adults are transitions (`V`); only infections of hosts are
new (`F`). The direct terms are differentiated numerically, so any
[`TransmissionFunction`](@ref) can be used; the defaults are those of
[`BsalStageTemperatureProblem`](@ref). `R0` is the spectral radius of `K`.

Only the `E_j` and `E_a` rows of `K` are non-zero and adult infections never produce eft
infections, so the eigenvalues are `K[1,1] = R_j` and `K[3,3] = R_a`
([`stage_reproduction_numbers`](@ref)).
"""
function stage_next_generation_matrix(p::ChowdhuryParameters;
        adult_contact::TransmissionFunction = _default_adult_contact(p),
        eft_contact::TransmissionFunction = _default_eft_contact(p),
        larval_dilution::Bool = true)
    d = stage_disease_free_state(p)
    h = 1e-7 * max(d.S_j, d.S_a, 1.0)
    # ∂(new infections)/∂I at the disease-free state; both N_T and A grow with I
    βj = _eft_direct(eft_contact, d.S_j, h, d.S_j + h) / h
    βa = _adult_direct(adult_contact, d.S_a, h, d.S_a + h, d.S_a + h + d.S_l,
        larval_dilution) / h
    rec = p.η_a * (1 - d.S_a / p.K)
    F = zeros(6, 6)
    F[1, 2] = βj
    F[1, 5] = p.ρ_j * p.b_j * d.S_j / p.κ_j
    F[3, 4] = βa
    F[3, 6] = p.ρ_a * p.b_a * d.S_a / p.κ_a
    V = zeros(6, 6)
    V[1, 1] = p.σ_j + p.δ_j + p.η_a
    V[2, 1] = -p.σ_j
    V[2, 2] = p.γ_j + p.δ_j + p.μ_j + p.η_a
    V[3, 1] = -rec
    V[3, 3] = p.σ_a + p.δ_a
    V[4, 2] = -rec
    V[4, 3] = -p.σ_a
    V[4, 4] = p.γ_a + p.δ_a + p.μ_a
    V[5, 2] = -p.ω_j
    V[5, 5] = p.ξ_j
    V[6, 4] = -p.ω_a
    V[6, 6] = p.ξ_a
    K = F / V
    return (; F, V, K, R0 = maximum(abs, eigvals(K)))
end

"""
    stage_reproduction_numbers(p::ChowdhuryParameters; kwargs...) -> NamedTuple

The two eigenvalues of the next-generation matrix (eqs. 5a, 5b; Appendix F), each split
into its direct and environmental routes: `R_j = direct_j + environmental_j` for efts,
`R_a = direct_a + environmental_a` for adults, and `R0 = max(R_j, R_a)`. Keywords as for
[`stage_next_generation_matrix`](@ref). For the default contact functions,

```
R_j = c_j P_j σ_j / (D_j (η_a + δ_j + σ_j))
      + ρ_j S_j* ω_j σ_j b_j / (κ_j ξ_j D_j (η_a + δ_j + σ_j)),
R_a = S_a*² c_a P_a σ_a / ((σ_a + δ_a) D_a (S_l* + S_a*)(χ_a + S_a*))
      + ρ_a ω_a b_a σ_a S_a* / (κ_a ξ_a (σ_a + δ_a) D_a),
```

with `D_j = δ_j + η_a + γ_j + μ_j` and `D_a = μ_a + δ_a + γ_a`. Table 3 is reproduced with
`variant = :table3` at 14 °C (both columns) and for `R_a` at 6 °C.
"""
function stage_reproduction_numbers(p::ChowdhuryParameters; kwargs...)
    m = stage_next_generation_matrix(p; kwargs...)
    Vi = inv(m.V)
    direct_j = m.F[1, 2] * Vi[2, 1]
    environmental_j = m.F[1, 5] * Vi[5, 1]
    direct_a = m.F[3, 4] * Vi[4, 3]
    environmental_a = m.F[3, 6] * Vi[6, 3]
    R_j = direct_j + environmental_j
    R_a = direct_a + environmental_a
    return (; R_j, R_a, R0 = max(R_j, R_a), direct_j, environmental_j, direct_a,
        environmental_a)
end

"""
    basic_reproduction_number(p::ChowdhuryParameters; kwargs...)

`R0 = max(R_j, R_a)`, the spectral radius of the next-generation matrix of Chowdhury et
al. (2025); see [`stage_reproduction_numbers`](@ref).
"""
basic_reproduction_number(p::ChowdhuryParameters; kwargs...) =
    stage_reproduction_numbers(p; kwargs...).R0

"""
    stage_elasticities(p::ChowdhuryParameters; number = :R_a, h = 1e-6) -> Dict{Symbol,Float64}

Elasticities `(θ/R) ∂R/∂θ` of `R_j` or `R_a` ([`stage_reproduction_numbers`](@ref)) with
respect to each parameter θ of `p` (Appendix G, eq. 6), by central differences in `log θ`.
Parameters equal to zero (the recovery rates) have elasticity zero. `ω_j = ω_a` and
`P_j = P_a` in the paper but are varied separately here.
"""
function stage_elasticities(p::ChowdhuryParameters; number::Symbol = :R_a, h::Real = 1e-6)
    number in (:R_j, :R_a) || throw(ArgumentError("number must be :R_j or :R_a"))
    R(q) = getfield(stage_reproduction_numbers(q), number)
    R0 = R(p)
    out = Dict{Symbol,Float64}()
    for f in fieldnames(ChowdhuryParameters)
        f === :temperature && continue
        θ = getfield(p, f)
        if θ == 0
            out[f] = 0.0
            continue
        end
        up = R(with(p; f => θ * exp(h)))
        dn = R(with(p; f => θ * exp(-h)))
        out[f] = (up - dn) / (2h * R0)
    end
    return out
end
