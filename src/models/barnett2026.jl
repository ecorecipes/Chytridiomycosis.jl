# Barnett, McMahon, Shepack, … , Johnson & Civitello (2026), "From protection to
# amplification: imperfect chytridiomycosis prophylaxis increases infections in wild
# amphibians". bioRxiv PREPRINT, not peer reviewed: doi:10.64898/2026.05.15.725113 (v2,
# 2026-08-12). A stochastic, stage-structured, spatial agent-based model (ABM) of Bd in a
# within-season cohort of Pseudacris regilla tadpoles that become metamorphs, with an
# imperfect vaccine (the Bd-metabolite prophylaxis) that changes one or more of: infection
# establishment (anti-infection resistance), zoospore shedding (anti-transmission
# resistance), zoosporangium clearance (anti-growth resistance) and the lethal load
# (tolerance). Units: days (one tick = one day); zoosporangia per host; zoospores per patch.
#
# PORT OF THE AUTHORS' CODE. Chytrid_Vaccine_ABM.nlogo (NetLogo 6.3.0) and the BehaviorSpace
# output are at https://github.com/kmbarnett100/protection_to_amplification (commit 29f4af6,
# 2026-08-12). solve(::BarnettABMProblem) follows `setup` and `go` procedure by procedure, and
# scripts/barnett2026_abm.jl compares the distributions of its outputs with all the authors'
# runs (data/barnett2026/reproduction.toml). Random number streams differ, so the match is in
# distribution, not run by run.
#
# One tick (`go`), in order:
# 1. tadpoles age and die (tad-mort); on tick 0, round(coverage × tadpoles) are vaccinated, and
#    those not infected get imm ~ U(1 − rv/2, 1 + rv/2)·v_efficacy and est·exp(c_est imm);
# 2. metamorphs die (meta-mort); tadpoles older than 55 days metamorphose with probability
#    0.11, all of them once any is older than 74 days; smax = baseline·exp(c_smax imm);
# 3. round(t_movement × tadpoles) tadpoles move to a random perimeter patch; round(m_land ×
#    metamorphs) metamorphs move to a random land patch by the pond, the others to a random
#    perimeter patch (every metamorph, every tick);
# 4. infected tadpoles below s_k, then those at or above s_k, then infected metamorphs:
#    maturation, clearance, shedding and self-reinfection;
# 5. infection from every perimeter patch, with the zoospores left over decaying.
#
# How the code differs from the text of the preprint. The code is followed throughout:
# - Clearance: a zoosporangium survives a day with probability exp(−0.2 exp(c_clear imm)),
#   0.819 at baseline, not 1 − 0.2.
# - Infection from a patch uses ALL hosts in the world: of the zoospores on a patch,
#   Z (1 − e^−h)·Σ(expo est)/h establish (rounded; no binomial draw), with h = z_mort + Σ expo
#   over every host, wherever it is. They are then shared over all hosts in the world, and
#   only those that fall on a host on that patch infect it; the others are lost.
# - The sharing weights (infprob) are fixed at hatching, so a vaccinated host's lower est
#   reduces the number of establishing zoospores for everyone but not its own share.
# - Relative variation ("±10%") is ±5% (a range of 10% of v_efficacy).
# - A host counts as infected (bd = 1) as soon as zoospores establish, but its stages advance
#   only while bd = 1, and bd returns to 0 whenever its zoosporangia number 0, even with
#   prezoosporangia on board. Those wait until the next infection restarts the clock.
# - Self-reinfection is f·est (not rounded) and overwrites the newest stage; the fractional
#   part of a load is dropped at the next clearance.
# - A tadpole at or above s_k loses its establishment probability for good (est = 0). A tadpole
#   that crosses s_k during a tick and stays there is processed twice in that tick, because the
#   second group is formed after the first has been updated.
# - Of a tadpole's released zoospores 10% are self-exposure; of the rest, 40% go to its patch,
#   58% of what is left to one neighbouring perimeter patch (8 neighbours) and the remainder to
#   a patch off the perimeter, where they infect nobody. Metamorphs put all the rest on their
#   own patch.
# - Every perimeter patch starts with 500 tadpoles, 1 of them infected (0.2%). The pond is 1 +
#   round(N(135, 18)) patches: a random non-border patch of a 21 × 21 world and its nearest
#   non-border patches. The perimeter has 28–49 patches.
# - Outputs are the metamorphs at the end of the last tick, except zoospore density, which is
#   the sum over ticks of the mean zoospores per perimeter patch before infection.
#
# Not ported: SimplePond, birth_pulses other than 1, plotting, and zoospores off the perimeter
# (they do not affect any host or output).

const BARNETT_MODES = (:infection, :transmission, :growth, :tolerance)

"""
    BarnettParameters(:barnett2026; kwargs...)

Parameters of the Barnett et al. (2026, **bioRxiv preprint**) agent-based model of Bd in a
cohort of *Pseudacris regilla* tadpoles and metamorphs, read from
`data/compendium/barnett2026.toml` (values of the preprint and of the authors' NetLogo code).
Rates and probabilities are per day.

- `tad_mort`, `meta_mort`: daily probabilities of natural death;
- `t_movement`: fraction of tadpoles moving to a random perimeter patch each day;
- `m_land`: fraction of metamorphs on land each day;
- `meta_start`, `meta_prob`, `meta_end`: a tadpole older than `meta_start` days metamorphoses
  with daily probability `meta_prob`; all do once one is older than `meta_end`;
- `last_day`: number of daily ticks;
- `tadpoles_per_patch`, `infected_per_patch`, `init_load`: tadpoles and infected tadpoles on
  each perimeter patch at the start, and the zoosporangia of each infected one;
- `exposure`, `est`, `z_mort`: exposure per host, establishment probability, zoospore death rate;
- `clearance`: zoosporangium clearance rate (daily survival `exp(−clearance)`);
- `shedding`: zoospores per zoosporangium per day; `self_reinfection`, `shed_own`,
  `shed_neighbour`: the share that re-exposes the host, the share of the rest that goes to the
  tadpole's patch, and the share of what is left that goes to a neighbouring perimeter patch;
- `smax`: lethal load of a metamorph; `s_k`: tadpole load above which no new infections
  establish;
- `v_efficacy`, `relative_variation`: a vaccinated host's response `imm` is uniform on
  `v_efficacy × [1 − relative_variation/2, 1 + relative_variation/2]`;
- `natural_imm`: increase in `imm` on clearing an infection (0 in every published run);
- `world_size`, `pond_mean`, `pond_sd`: side of the square world and the normal distribution of
  the pond size, in patches.
"""
Base.@kwdef struct BarnettParameters
    tad_mort::Float64
    meta_mort::Float64
    t_movement::Float64
    m_land::Float64
    meta_start::Int
    meta_prob::Float64
    meta_end::Int
    last_day::Int
    tadpoles_per_patch::Int
    infected_per_patch::Int
    init_load::Int
    exposure::Float64
    est::Float64
    z_mort::Float64
    clearance::Float64
    shedding::Float64
    self_reinfection::Float64
    shed_own::Float64
    shed_neighbour::Float64
    smax::Float64
    s_k::Float64
    v_efficacy::Float64
    relative_variation::Float64
    natural_imm::Float64
    world_size::Int
    pond_mean::Float64
    pond_sd::Float64
end

function BarnettParameters(key::Symbol; kwargs...)
    q = parameter_set(key)
    int(x) = round(Int, x)
    base = (; tad_mort = q.tad_mort, meta_mort = q.meta_mort, t_movement = q.t_movement,
            m_land = q.m_land, meta_start = int(q.meta_start), meta_prob = q.meta_prob,
            meta_end = int(q.meta_end), last_day = int(q.last_day),
            tadpoles_per_patch = int(q.tadpoles_per_patch),
            infected_per_patch = int(q.infected_per_patch), init_load = int(q.init_load),
            exposure = q.exposure, est = q.est, z_mort = q.z_mort, clearance = q.clearance,
            shedding = q.shedding, self_reinfection = q.self_reinfection,
            shed_own = q.shed_own, shed_neighbour = q.shed_neighbour_code, smax = q.smax,
            s_k = q.s_k, v_efficacy = q.v_efficacy, relative_variation = q.relative_variation,
            natural_imm = q.natural_imm, world_size = int(q.world_size),
            pond_mean = q.pond_mean, pond_sd = q.pond_sd)
    return BarnettParameters(; merge(base, values(kwargs))...)
end

time_unit(::Type{<:BarnettParameters}) = :day

with(p::BarnettParameters; kw...) =
    BarnettParameters(; merge(NamedTuple(f => getfield(p, f) for f in fieldnames(BarnettParameters)),
                              values(kw))...)

function _check(p::BarnettParameters)
    for f in (:tad_mort, :meta_mort, :t_movement, :m_land, :meta_prob, :est, :self_reinfection,
              :shed_own, :shed_neighbour)
        0 <= getfield(p, f) <= 1 || throw(ArgumentError("$f must lie in [0, 1]"))
    end
    0 <= p.infected_per_patch <= p.tadpoles_per_patch ||
        throw(ArgumentError("need 0 ≤ infected_per_patch ≤ tadpoles_per_patch"))
    p.world_size >= 5 || throw(ArgumentError("world_size must be at least 5"))
    p.meta_start <= p.meta_end <= p.last_day ||
        throw(ArgumentError("need meta_start ≤ meta_end ≤ last_day"))
    return p
end

"""
    BarnettVaccine(; coverage = 0, c_est = 0, c_shedding = 0, c_clear = 0, c_smax = 0)

An imperfect vaccine (prophylaxis) in the Barnett et al. (2026, preprint) model. `coverage`
is the fraction of tadpoles selected on the first day, and the constants `c` multiply a
vaccinated host's traits by `exp(c imm)`: establishment probability (`c_est`), shedding
(`c_shedding`), clearance rate (`c_clear`) and lethal load (`c_smax`). See
[`barnett_vaccine`](@ref) for the preprint's scenarios.
"""
Base.@kwdef struct BarnettVaccine
    coverage::Float64 = 0.0
    c_est::Float64 = 0.0
    c_shedding::Float64 = 0.0
    c_clear::Float64 = 0.0
    c_smax::Float64 = 0.0
    function BarnettVaccine(coverage, c_est, c_shedding, c_clear, c_smax)
        0 <= coverage <= 1 || throw(ArgumentError("coverage must lie in [0, 1]"))
        all(isfinite, (c_est, c_shedding, c_clear, c_smax)) ||
            throw(ArgumentError("the vaccine constants must be finite"))
        return new(coverage, c_est, c_shedding, c_clear, c_smax)
    end
end

# The authors' constants for a relative change e (BehaviorSpace): log(1 + e) to 3 decimals,
# except −10% (−0.11) and −100% (log 1e-7 = −16.118).
function _barnett_constant(e::Real)
    -1 <= e <= 1 || throw(ArgumentError("efficacy must lie in [0, 1]"))
    q = parameter_set(:barnett2026)
    isapprox(e, -0.1; atol = 1e-9) && return q.c_decrease_10
    isapprox(e, -1; atol = 1e-9) && return q.c_decrease_100
    return round(log1p(e); digits = 3)
end

"""
    barnett_vaccine(mode, efficacy; coverage, adverse = false) -> BarnettVaccine
    barnett_vaccine(modes::Pair...; coverage, adverse = false) -> BarnettVaccine

The vaccine of one of the preprint's scenarios. `mode` is one of `:infection` (decreases
establishment), `:transmission` (decreases shedding), `:growth` (increases clearance) or
`:tolerance` (increases the lethal load), and `efficacy` the size of the change (0.1–1.0 in
the preprint), converted to the authors' constants (`log(1 ± e)` rounded to 3 decimals; −0.11
for a 10% and −16.118 for a 100% decrease). With `adverse = true` the direction is reversed
(the "hypothetical adverse reactions"). Pairs such as `:tolerance => 0.8, :growth => 0.2`
combine modes; the authors' two-way runs used `coverage_two_way = 0.75` and no variation
between hosts (`relative_variation = 0`).
"""
function barnett_vaccine(modes::Pair{Symbol,<:Real}...; coverage::Real, adverse::Bool = false)
    field = (infection = :c_est, transmission = :c_shedding, growth = :c_clear,
             tolerance = :c_smax)
    kw = Dict{Symbol,Float64}()
    for (mode, e) in modes
        mode in BARNETT_MODES ||
            throw(ArgumentError("mode must be one of $BARNETT_MODES; got :$mode"))
        0 <= e <= 1 || throw(ArgumentError("efficacy must lie in [0, 1]"))
        protective = mode in (:growth, :tolerance) ? e : -e
        kw[field[mode]] = _barnett_constant(adverse ? -protective : protective)
    end
    return BarnettVaccine(; coverage = Float64(coverage), kw...)
end
barnett_vaccine(mode::Symbol, efficacy::Real; kw...) = barnett_vaccine(mode => efficacy; kw...)

"""
    BarnettABMProblem(p::BarnettParameters, vaccine = BarnettVaccine(); save_daily = false)

One run of the Barnett et al. (2026, **bioRxiv preprint**) agent-based model, a port of the
authors' NetLogo code. `solve(prob; rng)` returns the authors' BehaviorSpace outputs at the end
of the last day, as a `NamedTuple`:

- `survivors`: metamorphs alive (`count metamorphs`);
- `prevalence`, `unvaccinated_prevalence`, `vaccinated_prevalence`: fraction infected;
- `median_intensity`, `mean_intensity`, and the same with the prefix `unvaccinated_` or
  `vaccinated_`: load of metamorphs with zoosporangia (at least two of them);
- `unvaccinated_fraction`: fraction of survivors that are unvaccinated;
- `zoospore_density`: sum over days of the zoospores per perimeter patch;
- `bd_death_fraction`: share of metamorph deaths due to Bd; `realized_coverage`: fraction of
  tadpoles immunised; `n_perimeter`: number of perimeter patches;
- with `save_daily = true`, also `day`, `tadpoles`, `metamorphs`, `infected` and `zoospores`
  (on perimeter patches) at the end of each day.

Quantities the authors code as missing (−100; 0 in their older runs) are `NaN`.
"""
struct BarnettABMProblem
    p::BarnettParameters
    vaccine::BarnettVaccine
    save_daily::Bool
end
BarnettABMProblem(p::BarnettParameters, v::BarnettVaccine = BarnettVaccine();
                  save_daily::Bool = false) = BarnettABMProblem(_check(p), v, save_daily)

@inline _nlround(x) = floor(Int, x + 0.5)          # NetLogo's round

# Perimeter patches of a random pond, as lists of the neighbouring perimeter patches
# (8-neighbourhood) of each (procedure create-pond).
function _barnett_pond(p::BarnettParameters, rng::AbstractRNG)
    n = p.world_size
    inner = [(x, y) for x in 1:(n - 2) for y in 1:(n - 2)]
    c = rand(rng, inner)
    others = Random.shuffle!(rng, filter(!=(c), inner))
    sort!(others; by = q -> (q[1] - c[1])^2 + (q[2] - c[2])^2, alg = MergeSort)
    k = clamp(_nlround(p.pond_mean + p.pond_sd * randn(rng)), 1, length(others))
    pond = Set(@view others[1:k])
    push!(pond, c)
    edge(q) = any(d -> (q[1] + d[1], q[2] + d[2]) ∉ pond, ((1, 0), (-1, 0), (0, 1), (0, -1)))
    per = sort!(filter(edge, collect(pond)))
    index = Dict(q => i for (i, q) in enumerate(per))
    return [[index[(q[1] + a, q[2] + b)] for a in -1:1 for b in -1:1
             if (a, b) != (0, 0) && haskey(index, (q[1] + a, q[2] + b))] for q in per]
end

# Move a uniformly random subset of size k to the front of v.
function _partial_shuffle!(rng::AbstractRNG, v::Vector{Int}, k::Int)
    for i in 1:k
        r = rand(rng, i:length(v))
        v[i], v[r] = v[r], v[i]
    end
    return v
end

function CommonSolve.solve(prob::BarnettABMProblem; rng::AbstractRNG = Random.default_rng())
    p, v = prob.p, prob.vaccine
    nbrs = _barnett_pond(p, rng)
    P = length(nbrs)
    N = P * p.tadpoles_per_patch
    patch = repeat(1:P; inner = p.tadpoles_per_patch)      # 0 = on land
    metamorph = falses(N)
    bd = falses(N)
    spn = zeros(N)                                         # zoosporangia (can be fractional)
    pz = zeros(5, N)                                       # pz0 … pz4
    imm = zeros(N)
    est = fill(p.est, N)
    smax = fill(Inf, N)
    vaccinated = falses(N)
    for i in 1:P, k in 1:p.infected_per_patch
        j = (i - 1) * p.tadpoles_per_patch + k
        bd[j] = true
        spn[j] = p.init_load
    end
    tads = collect(1:N)
    mets = Int[]
    Z = zeros(Int, P)
    here = [Int[] for _ in 1:P]
    density = 0
    bd_deaths = 0
    base_deaths = 0
    coverage = 0.0
    nd = prob.save_daily ? p.last_day : 0
    d_tad = zeros(Int, nd); d_met = zeros(Int, nd); d_inf = zeros(Int, nd); d_zoo = zeros(Int, nd)

    function advance!(j)                                   # update-infections
        pz[5, j] = pz[4, j]; pz[4, j] = pz[3, j]; pz[3, j] = pz[2, j]; pz[2, j] = pz[1, j]
        pz[1, j] = 0.0
        spn[j] += pz[5, j]
        return nothing
    end
    function clear!(j)
        n = trunc(Int, spn[j])
        spn[j] = n > 0 ? rand(rng, Binomial(n, exp(-p.clearance * exp(v.c_clear * imm[j])))) : 0
        if spn[j] == 0
            bd[j] = false
            imm[j] += p.natural_imm
        end
        return nothing
    end
    release(j) = _nlround(p.shedding * exp(v.c_shedding * imm[j]) * spn[j])
    function shed_tadpole!(j, reinfect)
        z = release(j)
        z == 0 && (reinfect && (pz[1, j] = 0.0); return nothing)
        f = _nlround(p.self_reinfection * z)
        reinfect && (pz[1, j] = f * est[j])
        own = _nlround(p.shed_own * (z - f))
        near = _nlround(p.shed_neighbour * (z - f - own))
        i = patch[j]
        Z[i] += own
        isempty(nbrs[i]) || (Z[rand(rng, nbrs[i])] += near)
        return nothing
    end

    for tick in 0:(p.last_day - 1)
        age = tick + 1
        filter!(_ -> rand(rng) >= p.tad_mort, tads)
        if tick == 0 && v.coverage > 0
            k = _nlround(v.coverage * length(tads))
            _partial_shuffle!(rng, tads, k)
            rv = p.relative_variation
            for j in @view tads[1:k]
                bd[j] && continue                          # no post-exposure benefit
                imm[j] = p.v_efficacy * (1 - rv / 2 + rv * rand(rng))
                est[j] = p.est * exp(v.c_est * imm[j])
                vaccinated[j] = true
            end
            coverage = isempty(tads) ? 0.0 : count(vaccinated) / length(tads)
        end
        m0 = length(mets)
        filter!(_ -> rand(rng) >= p.meta_mort, mets)
        base_deaths += m0 - length(mets)
        if age > p.meta_start
            all_go = age > p.meta_end
            keep = 0
            for j in tads
                if all_go || rand(rng) < p.meta_prob
                    metamorph[j] = true
                    smax[j] = p.smax * exp(v.c_smax * imm[j])
                    push!(mets, j)
                else
                    keep += 1
                    tads[keep] = j
                end
            end
            resize!(tads, keep)
        end
        k = _nlround(p.t_movement * length(tads))
        _partial_shuffle!(rng, tads, k)
        for j in @view tads[1:k]
            patch[j] = rand(rng, 1:P)
        end
        k = _nlround(p.m_land * length(mets))
        _partial_shuffle!(rng, mets, k)
        for (a, j) in enumerate(mets)
            patch[j] = a <= k ? 0 : rand(rng, 1:P)
        end
        for j in tads
            (bd[j] && spn[j] < p.s_k) || continue
            advance!(j)
            clear!(j)
            shed_tadpole!(j, true)
        end
        for j in tads                                      # formed after the first group
            (bd[j] && spn[j] >= p.s_k) || continue
            est[j] = 0.0
            clear!(j)
            shed_tadpole!(j, false)
        end
        m0 = length(mets)
        filter!(mets) do j
            bd[j] || return true
            advance!(j)
            spn[j] >= smax[j] && return false
            clear!(j)
            z = release(j)
            f = _nlround(p.self_reinfection * z)
            pz[1, j] = f * est[j]
            patch[j] > 0 && (Z[patch[j]] += z - f)
            return true
        end
        bd_deaths += m0 - length(mets)
        density += _nlround(sum(Z) / P)
        # infection-step, over the hosts of the whole world
        n = length(tads) + length(mets)
        sum_est = 0.0
        foreach(empty!, here)
        for list in (tads, mets), j in list
            sum_est += est[j]
            patch[j] > 0 && push!(here[patch[j]], j)
        end
        h = p.z_mort + p.exposure * n
        left = exp(-h)
        for i in 1:P
            Z[i] > 0 || continue
            S = _nlround(Z[i] * (1 - left) * p.exposure * sum_est / h)
            Z[i] = _nlround(Z[i] * left)
            hosts = here[i]
            m = length(hosts)
            (S > 0 && m > 0) || continue
            r = rand(rng, Binomial(S, m / n))
            for (a, j) in enumerate(hosts)
                r == 0 && break
                x = a == m ? r : rand(rng, Binomial(r, 1 / (m - a + 1)))
                pz[1, j] += x
                r -= x
            end
        end
        ni = 0
        for list in (tads, mets), j in list
            pz[1, j] > 0 && (bd[j] = true)
            bd[j] && (ni += 1)
        end
        if prob.save_daily
            d_tad[age] = length(tads); d_met[age] = length(mets)
            d_inf[age] = ni; d_zoo[age] = sum(Z)
        end
    end

    n = length(mets)
    L = spn[mets]
    V = vaccinated[mets]
    inf = bd[mets]
    frac(a, b) = b == 0 ? NaN : a / b
    function intensity(sel)
        x = L[(L .> 0) .& sel]
        return length(x) >= 2 ? (median(x), sum(x) / length(x)) : (NaN, NaN)
    end
    all_, unv, vac = intensity(trues(n)), intensity(.!V), intensity(V)
    res = (; survivors = n, prevalence = frac(count(inf), n),
           unvaccinated_prevalence = frac(count(inf .& .!V), count(!, V)),
           vaccinated_prevalence = frac(count(inf .& V), count(V)),
           median_intensity = all_[1], mean_intensity = all_[2],
           unvaccinated_median_intensity = unv[1], unvaccinated_mean_intensity = unv[2],
           vaccinated_median_intensity = vac[1], vaccinated_mean_intensity = vac[2],
           unvaccinated_fraction = frac(count(!, V), n), zoospore_density = density,
           bd_death_fraction = frac(bd_deaths, bd_deaths + base_deaths),
           realized_coverage = coverage, n_perimeter = P)
    prob.save_daily || return res
    return merge(res, (; day = collect(1:p.last_day), tadpoles = d_tad, metamorphs = d_met,
                       infected = d_inf, zoospores = d_zoo))
end

const BARNETT_OUTPUTS = (:survivors, :prevalence, :unvaccinated_prevalence,
    :vaccinated_prevalence, :median_intensity, :mean_intensity, :unvaccinated_median_intensity,
    :unvaccinated_mean_intensity, :vaccinated_median_intensity, :vaccinated_mean_intensity,
    :unvaccinated_fraction, :zoospore_density, :bd_death_fraction)

"""
    barnett_summary(p::BarnettParameters, vaccine = BarnettVaccine(); nreps, rng) -> NamedTuple

Run `nreps` replicates (default `n_reps = 25` from the compendium) of
[`BarnettABMProblem`](@ref). Returns `runs` (a `NamedTuple` of vectors, one per output),
`mean` and `median` (over the runs where an output is defined, as the authors' figure code
does) and `dieoff`, the fraction of runs with no metamorph left.
"""
function barnett_summary(p::BarnettParameters, vaccine::BarnettVaccine = BarnettVaccine();
        nreps::Integer = round(Int, parameter_set(:barnett2026).n_reps),
        rng::AbstractRNG = Random.default_rng())
    nreps >= 1 || throw(ArgumentError("nreps must be at least 1"))
    prob = BarnettABMProblem(p, vaccine)
    rs = [solve(prob; rng) for _ in 1:nreps]
    runs = NamedTuple{BARNETT_OUTPUTS}(Tuple([Float64(getfield(r, k)) for r in rs]
                                             for k in BARNETT_OUTPUTS))
    stat(f) = map(x -> (y = filter(!isnan, x); isempty(y) ? NaN : f(y)), runs)
    return (; runs, mean = stat(mean), median = stat(median),
            dieoff = count(==(0), runs.survivors) / nreps, nreps = Int(nreps))
end

"""
    barnett_relative_change(p, vaccine; nreps, ncontrol, rng, control = nothing) -> NamedTuple

Percent change of a vaccination scenario relative to the unvaccinated control, as in the
authors' figure code: each run's output minus the **median** of the control runs, divided by
that median, times 100, averaged over the runs where it is defined. (The authors smooth these
run-level values over efficacy and coverage with a GAM; the package reports the mean of the
cell.) Returns `survivors` (final population), `intensity` (median load of infected
metamorphs), `prevalence`, `zoospore_density`, `unvaccinated_intensity` (Fig. 2a–b),
`unvaccinated_survival` (unvaccinated share of the survivors against `1 − coverage`, Fig.
2c–d), and the summaries `scenario` and `control`.

Defaults: `nreps = n_reps = 25` scenario runs and `ncontrol = n_control = 250` control runs
(compendium). Pass a precomputed `control` ([`barnett_summary`](@ref) of the unvaccinated
model) to reuse it across scenarios.
"""
function barnett_relative_change(p::BarnettParameters, vaccine::BarnettVaccine;
        nreps::Integer = round(Int, parameter_set(:barnett2026).n_reps),
        ncontrol::Integer = round(Int, parameter_set(:barnett2026).n_control),
        rng::AbstractRNG = Random.default_rng(), control = nothing)
    c = control === nothing ? barnett_summary(p; nreps = ncontrol, rng) : control
    s = barnett_summary(p, vaccine; nreps, rng)
    rel(a, b) = 100 * (a - b) / b
    u0 = 1 - vaccine.coverage
    return (; survivors = rel(s.mean.survivors, c.median.survivors),
            intensity = rel(s.mean.median_intensity, c.median.median_intensity),
            prevalence = rel(s.mean.prevalence, c.median.prevalence),
            zoospore_density = rel(s.mean.zoospore_density, c.median.zoospore_density),
            unvaccinated_intensity = rel(s.mean.unvaccinated_median_intensity,
                                         c.median.median_intensity),
            unvaccinated_survival = u0 > 0 ? rel(s.mean.unvaccinated_fraction, u0) : NaN,
            scenario = s, control = c)
end
