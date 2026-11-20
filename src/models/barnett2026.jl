# Barnett, McMahon, Shepack, … , Johnson & Civitello (2026), "From protection to
# amplification: imperfect chytridiomycosis prophylaxis increases infections in wild
# amphibians". bioRxiv PREPRINT, not peer reviewed: doi:10.64898/2026.05.15.725113 (v2,
# 2026-08-12). A stochastic, stage-structured, spatial agent-based model (ABM) of Bd in a
# within-season cohort of Pseudacris regilla tadpoles that become metamorphs, with an
# imperfect vaccine (the Bd-metabolite prophylaxis) that changes one or more of: infection
# establishment (anti-infection resistance), zoospore shedding (anti-transmission
# resistance), zoosporangium clearance (anti-growth resistance) and the lethal load
# (tolerance). Checked against the main text (Methods, 'Bd–amphibian–vaccine model' and the
# three scenario sections) and Supplementary Table S1. The authors' model is a NetLogo 6.3.0
# program; NO CODE was released with the preprint and none was found on GitHub or Zenodo, so
# there is no oracle beyond the text, and every choice below is the package's.
# Units: days (one tick = one day); zoosporangia per host; zoospores per patch.
#
# Processes as printed (Table S1):
# - Tadpoles die with probability tad_mort per day, metamorphs with meta_mort. From day
#   meta_start each tadpole metamorphoses with probability meta_prob per day; on day meta_end
#   all remaining tadpoles do. Metamorphs keep their load and immune traits.
# - A fraction t_movement of tadpoles move to a new perimeter pond patch each day. Metamorphs
#   move each day; a fraction m_land are on a terrestrial patch, the rest on perimeter patches.
# - Vaccination on day 0: a fraction `coverage` of tadpoles is selected at random. Selected
#   tadpoles that are already infected keep the baseline traits (metabolites give no
#   post-exposure benefit). A vaccinated host's trait x is x0 · exp(c_x · imm), with c_x a
#   per-scenario constant (0 for unvaccinated hosts) and imm the host's response.
# - Zoospores on a patch are removed by contact with hosts (exposure) and by exponential decay
#   (z_mort); a contacting zoospore establishes with probability est. Transmission is a
#   multinomial draw.
# - An established zoospore passes through 4 prezoosporangium stages, one per day, and then
#   becomes a zoosporangium. Zoosporangia survive each day with probability 1 − clearance
#   (binomial). Each shed `shedding` zoospores per day: 40% to the host's patch, 50% to a
#   neighbouring perimeter or inner patch, 10% re-expose the host itself (establishing with
#   probability est).
# - A tadpole at s_k zoosporangia cannot be further infected. A metamorph dies of Bd when its
#   load reaches smax. Tadpoles do not die of Bd.
# - Zoospores on deep-pond or terrestrial patches never infect anyone.
#
# Choices where the preprint is silent (all can be changed through BarnettParameters):
# - Population and pond size. The number of tadpoles (n_tadpoles = 1000, so that 0.2% is 2
#   tadpoles) and of perimeter patches (n_perimeter = 60, about the ring of Fig. S15) are not
#   given. Perimeter patches form a ring; patch i neighbours i ± 1.
# - Geometry of shedding. From a perimeter patch in a NetLogo grid (8 neighbours), 2 of the
#   8 neighbours are perimeter patches and the rest are deep or terrestrial, so the 50%
#   "neighbouring" share reaches the perimeter ring with probability perimeter_from_perimeter
#   = 2/8 (to i − 1 or i + 1 at random). From a terrestrial patch next to the ring, 3 of the 8
#   neighbours are perimeter patches (perimeter_from_land = 3/8, to i − 1, i, i + 1). The 40%
#   own-patch share of a host on land is lost.
# - Contact as competing hazards. Over one day a zoospore on a perimeter patch with n hosts
#   faces hazards z_mort (death) and exposure per host (contact). It survives with
#   probability exp(−h), h = z_mort + n·exposure, contacts host j with probability
#   exposure·(1 − exp(−h))/h and dies otherwise: a multinomial over n + 2 outcomes.
#   Contacted zoospores are removed whether or not they establish.
# - Shedding is a Poisson draw with mean shedding × modifier × load, split multinomially into
#   the own-patch, neighbour and self shares.
# - Individual response: imm ~ Uniform(v_efficacy (1 − rv), v_efficacy (1 + rv)).
# - Efficacy. A "change of e" to a baseline trait (10–100%) is c = log(1 + e) for an increase
#   and c = log(1 − e) for a decrease, so the vaccinated trait is x0 (1 ± e)^imm; a 100%
#   decrease (c = −Inf) sets the trait to 0. Probabilities are capped at 1.
# - Daily order: natural deaths; metamorphosis; movement; within-host step (clearance of
#   zoosporangia, maturation, shedding and self-reinfection); environmental transmission and
#   zoospore decay; Bd deaths of metamorphs with load ≥ smax.
# - Outputs. The preprint reports "final population size, infection prevalence, average
#   infection intensity, and spillover risk" (zoospore density) without defining them. The
#   package records, at the end of each day: prevalence (load > 0) and mean load of infected
#   hosts, among all live hosts and among live metamorphs, and zoospores per perimeter patch;
#   and returns their season means, the same quantities among survivors at the end, and the
#   number of survivors. barnett_relative_change compares scenario means with control means;
#   its default statistic is the metamorph one (metric = :metamorph), chosen because the field
#   trial sampled metamorphs and only metamorphs die of Bd. Which statistic is used changes
#   the sign of some results (see the barnett2026 targets in data/validation/targets.toml).

const BARNETT_MODES = (:infection, :transmission, :growth, :tolerance)

"""
    BarnettParameters(; tad_mort, meta_mort, t_movement, m_land, meta_start, meta_prob,
        meta_end, last_day, init_prevalence, init_load, exposure, est, z_mort,
        maturation_days, clearance, shedding, shed_own, shed_neighbour, self_reinfection,
        smax, s_k, v_efficacy, relative_variation, n_tadpoles = 1000, n_perimeter = 60,
        perimeter_from_perimeter = 0.25, perimeter_from_land = 0.375)
    BarnettParameters(:barnett2026; kwargs...)

Parameters of the Barnett et al. (2026, **bioRxiv preprint**) agent-based model of Bd in a
cohort of *Pseudacris regilla* tadpoles and metamorphs. Rates and probabilities are per day.

- `tad_mort`, `meta_mort`: daily probabilities of natural death;
- `t_movement`: fraction of tadpoles moving to another perimeter patch each day;
- `m_land`: fraction of metamorphs on land each day;
- `meta_start`, `meta_prob`, `meta_end`: metamorphosis window and daily probability;
- `last_day`: number of daily ticks;
- `init_prevalence`, `init_load`: fraction of tadpoles infected at the start, and their load
  (zoosporangia);
- `exposure`: per-host daily contact hazard with a zoospore on the same patch; `est`:
  establishment probability per contacting zoospore; `z_mort`: zoospore death rate;
- `maturation_days`: prezoosporangium stages before a zoosporangium;
- `clearance`: daily clearance probability of a zoosporangium; `shedding`: zoospores per
  zoosporangium per day; `shed_own`, `shed_neighbour`, `self_reinfection`: shares of shed
  zoospores to the host's patch, a neighbouring patch and the host itself;
- `smax`: lethal load of a metamorph; `s_k`: tadpole load above which no new infections
  establish;
- `v_efficacy`, `relative_variation`: a vaccinated host's response `imm` is uniform on
  `v_efficacy × [1 − relative_variation, 1 + relative_variation]`.

The last four keywords are **not in the preprint** and are package choices: the initial
number of tadpoles, the number of perimeter pond patches (a ring), and the probabilities that
the neighbouring share of shed zoospores reaches the perimeter ring from a perimeter patch
(2 of 8 neighbours) and from a terrestrial patch (3 of 8).

The keyed constructor reads every published value from `data/compendium/barnett2026.toml`.
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
    init_prevalence::Float64
    init_load::Int
    exposure::Float64
    est::Float64
    z_mort::Float64
    maturation_days::Int
    clearance::Float64
    shedding::Float64
    shed_own::Float64
    shed_neighbour::Float64
    self_reinfection::Float64
    smax::Float64
    s_k::Int
    v_efficacy::Float64
    relative_variation::Float64
    n_tadpoles::Int = 1000
    n_perimeter::Int = 60
    perimeter_from_perimeter::Float64 = 0.25
    perimeter_from_land::Float64 = 0.375
end

function BarnettParameters(key::Symbol; kwargs...)
    q = parameter_set(key)
    base = (; tad_mort = q.tad_mort, meta_mort = q.meta_mort, t_movement = q.t_movement,
            m_land = q.m_land, meta_start = round(Int, q.meta_start), meta_prob = q.meta_prob,
            meta_end = round(Int, q.meta_end), last_day = round(Int, q.last_day),
            init_prevalence = q.init_prevalence, init_load = round(Int, q.init_load),
            exposure = q.exposure, est = q.est, z_mort = q.z_mort,
            maturation_days = round(Int, q.maturation_days), clearance = q.clearance,
            shedding = q.shedding, shed_own = q.shed_own, shed_neighbour = q.shed_neighbour,
            self_reinfection = q.self_reinfection, smax = q.smax, s_k = round(Int, q.s_k),
            v_efficacy = q.v_efficacy, relative_variation = q.relative_variation)
    return BarnettParameters(; merge(base, values(kwargs))...)
end

time_unit(::Type{<:BarnettParameters}) = :day

with(p::BarnettParameters; kw...) =
    BarnettParameters(; merge(NamedTuple(f => getfield(p, f) for f in fieldnames(BarnettParameters)),
                              values(kw))...)

function _check(p::BarnettParameters)
    for f in (:tad_mort, :meta_mort, :t_movement, :m_land, :meta_prob, :init_prevalence, :est,
              :clearance, :perimeter_from_perimeter, :perimeter_from_land)
        0 <= getfield(p, f) <= 1 || throw(ArgumentError("$f must lie in [0, 1]"))
    end
    s = p.shed_own + p.shed_neighbour + p.self_reinfection
    isapprox(s, 1; atol = 1e-12) ||
        throw(ArgumentError("shed_own + shed_neighbour + self_reinfection must be 1; got $s"))
    p.n_tadpoles >= 1 && p.n_perimeter >= 3 ||
        throw(ArgumentError("need n_tadpoles ≥ 1 and n_perimeter ≥ 3"))
    p.maturation_days >= 1 || throw(ArgumentError("maturation_days must be at least 1"))
    p.meta_start <= p.meta_end <= p.last_day ||
        throw(ArgumentError("need meta_start ≤ meta_end ≤ last_day"))
    return p
end

"""
    BarnettVaccine(; coverage = 0, infection = 0, transmission = 0, growth = 0,
        tolerance = 0)

An imperfect vaccine (prophylaxis) in the Barnett et al. (2026, preprint) model. `coverage`
is the fraction of tadpoles selected on day 0. The other four are signed relative changes
("efficacy") to a vaccinated host's baseline trait:

- `infection`: change to the establishment probability `est` (anti-infection resistance is
  negative, for example `-0.3` for a 30% decrease);
- `transmission`: change to the shedding rate (anti-transmission resistance is negative);
- `growth`: change to the clearance probability (anti-growth resistance is positive);
- `tolerance`: change to the lethal load `smax` (tolerance is positive).

A change `e` enters as the constant `c = log(1 + e)` of the preprint's exponential modulation
`x0 exp(c imm)`; each must be at least `-1`. Adverse reactions (the preprint's "hypothetical
adverse reactions") use the opposite signs. See [`barnett_vaccine`](@ref).
"""
Base.@kwdef struct BarnettVaccine
    coverage::Float64 = 0.0
    infection::Float64 = 0.0
    transmission::Float64 = 0.0
    growth::Float64 = 0.0
    tolerance::Float64 = 0.0
    function BarnettVaccine(coverage, infection, transmission, growth, tolerance)
        0 <= coverage <= 1 || throw(ArgumentError("coverage must lie in [0, 1]"))
        all(>=(-1), (infection, transmission, growth, tolerance)) ||
            throw(ArgumentError("a relative change must be at least -1"))
        return new(coverage, infection, transmission, growth, tolerance)
    end
end

"""
    barnett_vaccine(mode, efficacy; coverage, adverse = false) -> BarnettVaccine
    barnett_vaccine(modes::Pair...; coverage, adverse = false) -> BarnettVaccine

The vaccine of one of the preprint's scenarios. `mode` is one of `:infection` (decreases
establishment), `:transmission` (decreases shedding), `:growth` (increases clearance) or
`:tolerance` (increases the lethal load), and `efficacy` the size of the change (0.1–1.0 in
the preprint). With `adverse = true` the direction is reversed (increased establishment,
increased shedding, decreased clearance, decreased tolerance). Pairs such as
`:tolerance => 0.8, :growth => 0.2` combine modes (the two-way scenarios, at the coverage
`coverage_two_way = 0.75` of the compendium).
"""
function barnett_vaccine(modes::Pair{Symbol,<:Real}...; coverage::Real, adverse::Bool = false)
    kw = Dict{Symbol,Float64}()
    for (mode, e) in modes
        mode in BARNETT_MODES ||
            throw(ArgumentError("mode must be one of $BARNETT_MODES; got :$mode"))
        0 <= e <= 1 || throw(ArgumentError("efficacy must lie in [0, 1]"))
        protective = mode in (:growth, :tolerance) ? e : -e
        kw[mode] = adverse ? -protective : protective
    end
    return BarnettVaccine(; coverage = Float64(coverage), kw...)
end
barnett_vaccine(mode::Symbol, efficacy::Real; kw...) = barnett_vaccine(mode => efficacy; kw...)

"""
    BarnettABMProblem(p::BarnettParameters, vaccine = BarnettVaccine(); save_daily = false)

One run of the Barnett et al. (2026, **bioRxiv preprint**) agent-based model, from day 0 to
`p.last_day`. `solve(prob; rng)` returns a `NamedTuple`:

- `survivors`: live hosts at the end of the last day;
- season means over the days (end of each day; days where a quantity is undefined are
  skipped): `prevalence` (fraction of live hosts with a zoosporangium load above 0),
  `intensity` (mean load of infected live hosts), `unvaccinated_intensity` (the same for
  unvaccinated hosts), `metamorph_prevalence`, `metamorph_intensity`,
  `metamorph_unvaccinated_intensity` (the same among metamorphs) and `zoospore_density`
  (zoospores per perimeter patch);
- the same at the end of the last day, among survivors: `final_prevalence`,
  `final_intensity`, `final_unvaccinated_intensity`, `final_zoospore_density` (`NaN` when
  undefined);
- `vaccinated_survivors`, `unvaccinated_survivors`, `n_vaccinated` (hosts vaccinated on
  day 0);
- `loads`, `vaccinated`: load and vaccination status of each survivor;
- with `save_daily = true`, also `day`, `tadpoles`, `metamorphs`, `infected` and `zoospores`
  (total on perimeter patches) at the end of each day.
"""
struct BarnettABMProblem
    p::BarnettParameters
    vaccine::BarnettVaccine
    save_daily::Bool
end
BarnettABMProblem(p::BarnettParameters, v::BarnettVaccine = BarnettVaccine();
                  save_daily::Bool = false) = BarnettABMProblem(_check(p), v, save_daily)

# Vaccinated trait multiplier x/x0 = (1 + change)^imm, with (1 + change) = 0 giving 0.
@inline _barnett_mult(change::Float64, imm::Float64) =
    change == 0 ? 1.0 : (1 + change == 0 ? 0.0 : exp(log1p(change) * imm))

@inline _ring(i::Int, n::Int) = mod1(i, n)

# Draw a multinomial with equal probabilities over the k hosts in `hosts`, adding to `out`.
function _equal_split!(out::Vector{Int}, hosts::Vector{Int}, total::Int, rng::AbstractRNG)
    k = length(hosts)
    rem = total
    for (a, j) in enumerate(hosts)
        rem == 0 && break
        x = a == k ? rem : rand(rng, Binomial(rem, 1 / (k - a + 1)))
        out[j] += x
        rem -= x
    end
    return out
end

function CommonSolve.solve(prob::BarnettABMProblem; rng::AbstractRNG = Random.default_rng())
    p, v = prob.p, prob.vaccine
    N, P, M = p.n_tadpoles, p.n_perimeter, p.maturation_days
    alive = trues(N)
    metamorph = falses(N)
    on_land = falses(N)
    patch = rand(rng, 1:P, N)
    load = zeros(Int, N)                  # zoosporangia
    pz = zeros(Int, M, N)                 # prezoosporangia, stage 1 = newest
    vaccinated = falses(N)
    est = fill(p.est, N)
    clear = fill(p.clearance, N)
    shed = fill(p.shedding, N)
    smax = fill(p.smax, N)
    Z = zeros(Int, P)                     # zoospores on perimeter patches

    # initial infections
    n_inf = round(Int, p.init_prevalence * N)
    for j in Random.randperm(rng, N)[1:n_inf]
        load[j] = p.init_load
    end
    # vaccination on day 0 (only uninfected selected tadpoles respond)
    n_sel = round(Int, v.coverage * N)
    rv = p.relative_variation
    for j in Random.randperm(rng, N)[1:n_sel]
        load[j] > 0 && continue
        vaccinated[j] = true
        imm = p.v_efficacy * (1 - rv + 2rv * rand(rng))
        est[j] = min(1.0, p.est * _barnett_mult(v.infection, imm))
        shed[j] = p.shedding * _barnett_mult(v.transmission, imm)
        clear[j] = min(1.0, p.clearance * _barnett_mult(v.growth, imm))
        smax[j] = p.smax * _barnett_mult(v.tolerance, imm)
    end
    n_vaccinated = count(vaccinated)

    nd = prob.save_daily ? p.last_day : 0
    d_tad = zeros(Int, nd); d_met = zeros(Int, nd); d_inf = zeros(Int, nd); d_zoo = zeros(Int, nd)
    newinf = zeros(Int, N)
    hosts_on = [Int[] for _ in 1:P]
    deposit = zeros(Int, P)

    acc_prev = 0.0; acc_int = 0.0; acc_uint = 0.0; acc_zoo = 0.0
    n_prev = 0; n_int = 0; n_uint = 0
    acc_mint = 0.0; acc_muint = 0.0; n_mint = 0; n_muint = 0; acc_mprev = 0.0; n_mprev = 0

    @inline function can_establish(j)
        return metamorph[j] || load[j] < p.s_k
    end

    for day in 1:p.last_day
        # 1. natural mortality
        for j in 1:N
            alive[j] || continue
            rand(rng) < (metamorph[j] ? p.meta_mort : p.tad_mort) && (alive[j] = false)
        end
        # 2. metamorphosis
        if day >= p.meta_start
            for j in 1:N
                (alive[j] && !metamorph[j]) || continue
                if day >= p.meta_end || rand(rng) < p.meta_prob
                    metamorph[j] = true
                end
            end
        end
        # 3. movement
        for j in 1:N
            alive[j] || continue
            if metamorph[j]
                patch[j] = _ring(patch[j] + rand(rng, -1:1), P)
                on_land[j] = rand(rng) < p.m_land
            elseif rand(rng) < p.t_movement
                patch[j] = _ring(patch[j] + rand(rng, (-1, 1)), P)
            end
        end
        # 4. within-host: clearance, maturation, shedding, self-reinfection
        fill!(deposit, 0)
        fill!(newinf, 0)
        for j in 1:N
            alive[j] || continue
            if load[j] > 0 && clear[j] > 0
                load[j] = rand(rng, Binomial(load[j], 1 - clear[j]))
            end
            load[j] += pz[M, j]
            for s in M:-1:2
                pz[s, j] = pz[s - 1, j]
            end
            pz[1, j] = 0
            load[j] > 0 || continue
            μ = shed[j] * load[j]
            μ > 0 || continue
            total = rand(rng, Poisson(μ))
            own = rand(rng, Binomial(total, p.shed_own))
            nb = rand(rng, Binomial(total - own, p.shed_neighbour / (1 - p.shed_own)))
            self = total - own - nb
            if self > 0 && can_establish(j)
                newinf[j] += rand(rng, Binomial(self, est[j]))
            end
            i = patch[j]
            if on_land[j]
                k = rand(rng, Binomial(nb, p.perimeter_from_land))
                deposit[_ring(i + rand(rng, -1:1), P)] += k
            else
                deposit[i] += own
                k = rand(rng, Binomial(nb, p.perimeter_from_perimeter))
                deposit[_ring(i + rand(rng, (-1, 1)), P)] += k
            end
        end
        Z .+= deposit
        # 5. environmental transmission and zoospore decay (competing hazards)
        foreach(empty!, hosts_on)
        for j in 1:N
            (alive[j] && !on_land[j]) && push!(hosts_on[patch[j]], j)
        end
        for i in 1:P
            Zi = Z[i]
            Zi == 0 && continue
            hs = hosts_on[i]
            n = length(hs)
            h = p.z_mort + n * p.exposure
            h > 0 || continue
            gone = rand(rng, Binomial(Zi, -expm1(-h)))
            Z[i] = Zi - gone
            n == 0 && continue
            C = rand(rng, Binomial(gone, n * p.exposure / h))
            C == 0 && continue
            tmp = zeros(Int, n)
            _equal_split!(tmp, collect(1:n), C, rng)
            for (a, j) in enumerate(hs)
                tmp[a] > 0 && can_establish(j) || continue
                newinf[j] += rand(rng, Binomial(tmp[a], est[j]))
            end
        end
        for j in 1:N
            alive[j] && (pz[1, j] += newinf[j])
        end
        # 6. Bd-induced mortality of metamorphs
        for j in 1:N
            if alive[j] && metamorph[j] && load[j] > 0 && load[j] >= smax[j]
                alive[j] = false
            end
        end
        # season accumulators (live hosts at the end of the day)
        na = 0; ni = 0; sl = 0; nui = 0; sul = 0; nmi = 0; sml = 0; numi = 0; suml = 0; nm = 0
        for j in 1:N
            alive[j] || continue
            na += 1
            metamorph[j] && (nm += 1)
            if load[j] > 0
                ni += 1; sl += load[j]
                vaccinated[j] || (nui += 1; sul += load[j])
                if metamorph[j]
                    nmi += 1; sml += load[j]
                    vaccinated[j] || (numi += 1; suml += load[j])
                end
            end
        end
        nm > 0 && (acc_mprev += nmi / nm; n_mprev += 1)
        nmi > 0 && (acc_mint += sml / nmi; n_mint += 1)
        numi > 0 && (acc_muint += suml / numi; n_muint += 1)
        na > 0 && (acc_prev += ni / na; n_prev += 1)
        ni > 0 && (acc_int += sl / ni; n_int += 1)
        nui > 0 && (acc_uint += sul / nui; n_uint += 1)
        acc_zoo += sum(Z) / P
        if prob.save_daily
            d_tad[day] = count(j -> alive[j] && !metamorph[j], 1:N)
            d_met[day] = count(j -> alive[j] && metamorph[j], 1:N)
            d_inf[day] = ni
            d_zoo[day] = sum(Z)
        end
    end

    surv = findall(alive)
    loads = load[surv]
    vacc = vaccinated[surv]
    inf = loads .> 0
    n = length(surv)
    ninf = count(inf)
    unv_inf = inf .& .!vacc
    q(a, k) = k == 0 ? NaN : a / k
    res = (; survivors = n,
           prevalence = q(acc_prev, n_prev), intensity = q(acc_int, n_int),
           zoospore_density = acc_zoo / p.last_day,
           unvaccinated_intensity = q(acc_uint, n_uint),
           metamorph_prevalence = q(acc_mprev, n_mprev),
           metamorph_intensity = q(acc_mint, n_mint),
           metamorph_unvaccinated_intensity = q(acc_muint, n_muint),
           final_prevalence = q(ninf, n), final_intensity = q(sum(loads[inf]; init = 0), ninf),
           final_zoospore_density = sum(Z) / P,
           final_unvaccinated_intensity = q(sum(loads[unv_inf]; init = 0), count(unv_inf)),
           vaccinated_survivors = count(vacc), unvaccinated_survivors = n - count(vacc),
           n_vaccinated, loads, vaccinated = collect(vacc))
    prob.save_daily || return res
    return merge(res, (; day = collect(1:p.last_day), tadpoles = d_tad, metamorphs = d_met,
                       infected = d_inf, zoospores = d_zoo))
end

_barnett_nanmean(x) = (y = filter(!isnan, x); isempty(y) ? NaN : sum(y) / length(y))

"""
    barnett_summary(p::BarnettParameters, vaccine = BarnettVaccine(); nreps, rng) -> NamedTuple

Run `nreps` replicates (default `n_reps = 25` from the compendium) of
[`BarnettABMProblem`](@ref) and return the means over runs of `survivors` and of the season
means `prevalence`, `intensity`, `zoospore_density`, `unvaccinated_intensity`, and of the
end-of-season `final_prevalence`, `final_intensity`, `final_zoospore_density`,
`final_unvaccinated_intensity` (runs where a quantity is undefined are skipped). Also
`median_prevalence` and `median_final_prevalence`, the fraction of survivors that are
unvaccinated (`unvaccinated_fraction`), and the per-run values in `runs`.
"""
function barnett_summary(p::BarnettParameters, vaccine::BarnettVaccine = BarnettVaccine();
        nreps::Integer = round(Int, parameter_set(:barnett2026).n_reps),
        rng::AbstractRNG = Random.default_rng())
    nreps >= 1 || throw(ArgumentError("nreps must be at least 1"))
    prob = BarnettABMProblem(p, vaccine)
    rs = [solve(prob; rng) for _ in 1:nreps]
    keys = (:survivors, :prevalence, :intensity, :zoospore_density, :unvaccinated_intensity,
            :final_prevalence, :final_intensity, :final_zoospore_density,
            :final_unvaccinated_intensity, :metamorph_prevalence, :metamorph_intensity,
            :metamorph_unvaccinated_intensity)
    unvf = [r.survivors == 0 ? NaN : r.unvaccinated_survivors / r.survivors for r in rs]
    runs = merge(NamedTuple{keys}(Tuple([Float64(getfield(r, k)) for r in rs] for k in keys)),
                 (; unvaccinated_fraction = unvf))
    means = NamedTuple{keys}(Tuple(_barnett_nanmean(runs[k]) for k in keys))
    return merge(means, (; median_prevalence = _barnett_median(runs.prevalence),
                         median_final_prevalence = _barnett_median(runs.final_prevalence),
                         unvaccinated_fraction = _barnett_nanmean(unvf), runs,
                         nreps = Int(nreps)))
end

function _barnett_median(x)
    v = sort(filter(!isnan, x))
    isempty(v) && return NaN
    return (v[(end + 1) ÷ 2] + v[end ÷ 2 + 1]) / 2
end

"""
    barnett_relative_change(p, vaccine; metric = :metamorph, nreps, ncontrol, rng,
        control = nothing) -> NamedTuple

Percent change of a vaccination scenario relative to the unvaccinated control, as in the
preprint's figures: `100 (scenario mean − control mean) / control mean` for `survivors`
(final population size), `prevalence`, `intensity`, `zoospore_density` and
`unvaccinated_intensity` (infected unvaccinated hosts against all infected hosts of the
control, Fig. 2a–b). The preprint does not define these statistics, so three readings are
offered: `metric = :metamorph` (default) uses season means over metamorphs (prevalence among
live metamorphs, mean load of infected live metamorphs), the stage the field trial sampled and
the only stage that dies of Bd; `:season` uses season means over all live hosts; `:final`
uses the survivors at the end of the last day. Zoospore density is the season mean for
`:metamorph` and `:season`.
Defaults: `nreps = n_reps = 25` scenario runs and `ncontrol = n_control = 250` control runs
(compendium). Pass a precomputed `control` ([`barnett_summary`](@ref) of the unvaccinated
model) to reuse it across scenarios.

`unvaccinated_survival` is the percent change in the unvaccinated share of survivors relative
to the unvaccinated share at the start, `1 − coverage` (Fig. 2c–d).
"""
function barnett_relative_change(p::BarnettParameters, vaccine::BarnettVaccine;
        metric::Symbol = :metamorph,
        nreps::Integer = round(Int, parameter_set(:barnett2026).n_reps),
        ncontrol::Integer = round(Int, parameter_set(:barnett2026).n_control),
        rng::AbstractRNG = Random.default_rng(), control = nothing)
    metric in (:season, :final, :metamorph) ||
        throw(ArgumentError("metric must be :season, :final or :metamorph"))
    c = control === nothing ? barnett_summary(p; nreps = ncontrol, rng) : control
    s = barnett_summary(p, vaccine; nreps, rng)
    rel(a, b) = 100 * (a - b) / b
    k(x) = metric === :season ? x :
           metric === :final ? Symbol(:final_, x) :
           x === :zoospore_density ? x : Symbol(:metamorph_, x)
    u0 = 1 - vaccine.coverage
    return (; survivors = rel(s.survivors, c.survivors),
            prevalence = rel(s[k(:prevalence)], c[k(:prevalence)]),
            intensity = rel(s[k(:intensity)], c[k(:intensity)]),
            zoospore_density = rel(s[k(:zoospore_density)], c[k(:zoospore_density)]),
            unvaccinated_intensity = rel(s[k(:unvaccinated_intensity)], c[k(:intensity)]),
            unvaccinated_survival = u0 > 0 ? rel(s.unvaccinated_fraction, u0) : NaN,
            scenario = s, control = c)
end
