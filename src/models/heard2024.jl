# Heard, Scroggie, Hollanders & Scheele (2024) J Anim Ecol 93:1670-1683,
# doi:10.1111/1365-2656.14177. Data and code: Zenodo doi:10.5281/zenodo.10457229 (CC BY 4.0).
#
# PUBLISHED MODEL. Checked against the published paper (Methods 2.2; read 2026-11-20) and
# implemented from the authors' code: the paper describes the model in outline and defers
# details (carrying capacity, dispersal proportion and range, initial abundance) to its
# Supplementary Material, which was not available. Everything below follows
# STEPS_model.R and Run_STEPS_simulations.R in the archive and the R package `steps` 1.3.0
# (CRAN archive; the version current when the paper was written), whose source was read for
# the semantics of growth(), kernel_dispersal() and ceiling_density().
#
# Bd's role: Bd acts ONLY through an imposed change in the annual survival of sub-adult and
# adult females, from the pre-Bd (historical, museum specimens) to the post-Bd (contemporary)
# catch-curve estimate (src/fitting/heard2024.jl). There is no infection, load or
# transmission process.
#
# The model (female only, annual step, one habitat patch per landscape cell):
#
#   A = [ F_s  F_a ]    F_s = mets × 0.2 × sex_ratio × met_survival   (sub-adult fecundity)
#       [ S_s  S_a ]    F_a = mets × sex_ratio × met_survival         (adult fecundity)
#                       S_s = S_a = adult survival
#
# Each year, in the order of steps::population_dynamics (change, dispersal, density dependence):
#
# 1. Growth with global (landscape-wide) environmental stochasticity: one Normal draw per
#    non-zero entry of A, SD = (F_s, F_a) × mets_cv for fecundities and survival_sd for both
#    survivals, shared by all patches; entries are truncated at 0, survivals and F_s (a
#    diagonal entry: steps caps every entry on or below the diagonal) at 1. Demographic
#    stochasticity ("full"): new sub-adults = Poisson(F_s n_s) + Poisson(F_a n_a); adults =
#    Binomial(n_s, S_s) + Binomial(n_a, S_a) (sub-adults never stay sub-adults).
# 2. Dispersal: in every occupied patch, round(prop_dispersing × n) individuals of each stage
#    (R's round, half to even) leave and are split multinomially among all habitat cells
#    closer than max_distance, INCLUDING the origin, with weights exp(−d/200 m).
# 3. Ceiling density dependence: patches with more than K sub-adults + adults are scaled down
#    to K, then the whole landscape is "rounded" by steps' round_pop: one multinomial draw of
#    round(total) individuals over all patch × stage cells, with probabilities proportional
#    to their (scaled) numbers. Every year this re-allocates individuals at random among
#    patches and stages while preserving expectations (as coded in steps 1.3.0).
#
# The metapopulation is extinct if no female is left after 50 years.
#
# As-coded details (switch `as_coded = true`, the default, reproduces them):
# - steps' id2coord() puts a cell in the last column of a row at column 0 of the NEXT row, so
#   dispersers leaving such a cell use a neighbourhood shifted by one row and one column.
#   With `as_coded = false` the true cell position is used.
# - the landscape-wide multinomial re-allocation of step 3. With `as_coded = false` only
#   patches above K are rounded (to the nearest integer) and nothing is re-allocated.
#
# Unpublished choices (the code is silent):
# - survival_sd is the SD of exp(−h), h ~ Normal(−log S, −log S × mortality_cv). The R code
#   estimates it from 10,000 random draws; we use the exact moments of that lognormal
#   variable, sqrt(exp(−2μ + 2σ²) − exp(−2μ + σ²)), with μ = −log S, σ = μ × cv.

const HEARD_SPECIES = (:aurea, :raniformis)
const HEARD_EPOCHS = (:historical, :contemporary)

_heard_set(species::Symbol, epoch::Symbol) =
    (species === :aurea ? "Laur_" : "Lran_") * string(epoch)

function _heard_check(species, epoch)
    species in HEARD_SPECIES ||
        throw(ArgumentError("species must be one of $HEARD_SPECIES; got :$species"))
    epoch in HEARD_EPOCHS ||
        throw(ArgumentError("epoch must be one of $HEARD_EPOCHS; got :$epoch"))
    return nothing
end

"""
    HeardParameters(; adult_survival, carrying_capacity, prop_dispersing, mets, mets_cv,
                    mortality_cv, subadult_survival = adult_survival, met_survival = 0.1,
                    sex_ratio = 0.5, subadult_breeding_fraction = 0.2, initial_fraction = 0.5,
                    time_steps = 50, cell_size = 65.0, distance_decay = 200.0,
                    max_distance = 1500.0)
    HeardParameters(:heard2024; species = :raniformis, epoch = :contemporary, mets = 24,
                    mets_cv = 0.5, mortality_cv = 0.5, survival = :as_coded)

Parameters of the female-only, two-stage spatial metapopulation model of Heard et al. (2024)
for *Litoria aurea* (`species = :aurea`) and *L. raniformis* (`:raniformis`):

- `adult_survival`, `subadult_survival`: annual survival of adult and sub-adult females
  (set equal in the paper's simulations);
- `mets`: metamorphs per adult female per year (scenario, 3–30), `met_survival` their
  survival to the sub-adult stage, `sex_ratio` the female fraction, and
  `subadult_breeding_fraction` the fraction of sub-adults that breed;
- `mets_cv`: coefficient of variation of the fecundities (environmental stochasticity),
  `mortality_cv`: coefficient of variation of the mortality hazard (both scenarios, 0.1–1);
- `carrying_capacity`: ceiling on sub-adult plus adult females per patch;
- `prop_dispersing`: proportion of each stage dispersing per year;
- `cell_size`, `distance_decay`, `max_distance`: landscape cell width and dispersal kernel
  `exp(−d/distance_decay)` truncated at `max_distance` (m);
- `initial_fraction`, `time_steps`: each stage starts at `round(K/2) × initial_fraction`
  females per patch; the horizon is `time_steps` years.

The keyed constructor reads the species × epoch set (`epoch = :historical` for pre-Bd,
`:contemporary` for post-Bd survival) from the compendium; `mets`, `mets_cv` and
`mortality_cv` are scenario values the paper varies over a grid (see
[`heard_scenario_grid`](@ref)). `survival = :as_coded` (default) uses the survival the
authors' simulation script enters (0.319, 0.641, 0.352, 0.673 for *L. raniformis* post- and
pre-Bd and *L. aurea* post- and pre-Bd); `:as_printed` uses the posterior means printed in
the paper's Results (0.318, 0.655, 0.352, 0.688).

**Bd acts only through `adult_survival`** (an imposed change between epochs); the model has
no infection process. The carrying capacities and dispersal proportions are given only in
the paper's Supplementary Material and come from the authors' code (unverified records).
"""
Base.@kwdef struct HeardParameters
    adult_survival::Float64
    subadult_survival::Float64 = adult_survival
    carrying_capacity::Float64
    prop_dispersing::Float64
    mets::Float64
    mets_cv::Float64
    mortality_cv::Float64
    met_survival::Float64 = 0.1
    sex_ratio::Float64 = 0.5
    subadult_breeding_fraction::Float64 = 0.2
    initial_fraction::Float64 = 0.5
    time_steps::Int = 50
    cell_size::Float64 = 65.0
    distance_decay::Float64 = 200.0
    max_distance::Float64 = 1500.0
end

function HeardParameters(key::Symbol; species::Symbol = :raniformis,
        epoch::Symbol = :contemporary, mets::Real = 24, mets_cv::Real = 0.5,
        mortality_cv::Real = 0.5, survival::Symbol = :as_coded)
    _heard_check(species, epoch)
    survival in (:as_coded, :as_printed) ||
        throw(ArgumentError("survival must be :as_coded or :as_printed; got :$survival"))
    v = parameter_set(key; set = _heard_set(species, epoch))
    S = survival === :as_coded ? v.adult_survival : v.adult_survival_printed
    return HeardParameters(; adult_survival = S, subadult_survival = S, carrying_capacity = v.carrying_capacity,
        prop_dispersing = v.prop_dispersing, mets = Float64(mets),
        mets_cv = Float64(mets_cv), mortality_cv = Float64(mortality_cv),
        met_survival = v.met_survival, sex_ratio = v.sex_ratio,
        subadult_breeding_fraction = v.subadult_breeding_fraction,
        initial_fraction = v.initial_fraction, time_steps = round(Int, v.time_steps),
        cell_size = v.cell_size, distance_decay = v.distance_decay,
        max_distance = v.max_distance)
end

time_unit(::Type{<:HeardParameters}) = :year

with(p::HeardParameters; kw...) =
    HeardParameters(; (f => get(kw, f, getfield(p, f)) for f in fieldnames(HeardParameters))...)

"""
    heard_scenario_grid(key = :heard2024)

The scenario grid of Heard et al. (2024): a `NamedTuple` of vectors `mets` (3:3:30),
`mets_cv` and `mortality_cv` (0.1:0.1:1), read from the compendium. The paper crosses them
with 20 landscapes ([`heard2024_landscapes`](@ref)), giving 20,000 simulations per species
and epoch.
"""
function heard_scenario_grid(key::Symbol = :heard2024)
    v = parameter_set(key)
    mets = collect(v.mets_min:v.mets_step:v.mets_max)
    cv = [round(x; digits = 10) for x in v.cv_min:v.cv_step:v.cv_max]
    return (mets = mets, mets_cv = cv, mortality_cv = copy(cv))
end

"""
    heard_projection_matrix(p::HeardParameters) -> Matrix{Float64}

The mean 2 × 2 Lefkovitch matrix `[F_s F_a; S_s S_a]` (sub-adults, adults) with
`F_a = mets × sex_ratio × met_survival` and `F_s = subadult_breeding_fraction × F_a`.
"""
function heard_projection_matrix(p::HeardParameters)
    Fa = p.mets * p.sex_ratio * p.met_survival
    Fs = p.mets * p.subadult_breeding_fraction * p.sex_ratio * p.met_survival
    return [Fs Fa; p.subadult_survival p.adult_survival]
end

"""
    heard_survival_sd(survival, mortality_cv) -> Float64

SD of the annual survival `exp(−h)` when the hazard `h ~ Normal(−log S, −log S × cv)`, as
the environmental-stochasticity SD of the survival entries. The authors estimate it from
10,000 random draws; this is the exact value (the moments of a lognormal variable).
"""
function heard_survival_sd(S::Real, cv::Real)
    0 < S < 1 || throw(ArgumentError("survival must be in (0, 1)"))
    μ = -log(S)
    σ = μ * cv
    return sqrt(max(exp(-2μ + 2σ^2) - exp(-2μ + σ^2), 0.0))
end

"""
    heard_lambda(p::HeardParameters) -> Float64

Dominant eigenvalue of [`heard_projection_matrix`](@ref): the deterministic annual growth
rate of a single, closed, uncapped patch.
"""
heard_lambda(p::HeardParameters) = maximum(real, eigvals(heard_projection_matrix(p)))

# ---- landscape and dispersal -----------------------------------------------------------------

"""
    HeardMetapopulationProblem(p::HeardParameters, habitat::AbstractMatrix{Bool};
                               as_coded = true)
    HeardMetapopulationProblem(p, landscape::Integer; as_coded = true)

One simulation of the Heard et al. (2024) spatial metapopulation on a grid of cells, where
`habitat[r, c]` marks a habitat patch (an integer selects one of the 20 bundled landscapes,
[`heard2024_landscapes`](@ref)). `solve(prob; rng = Random.default_rng())` returns a
[`HeardMetapopulationSolution`](@ref).

Each year: growth with global environmental and full demographic stochasticity, kernel
dispersal, then ceiling density dependence, as in the R package `steps` 1.3.0 that the
authors used. With `as_coded = true` (default) two quirks of `steps` are kept: the shifted
dispersal neighbourhood of last-column cells, and the landscape-wide multinomial
re-allocation that `steps` uses to round populations after the ceiling is applied. See the
source header for details.
"""
struct HeardMetapopulationProblem
    p::HeardParameters
    habitat::BitMatrix
    as_coded::Bool
end
HeardMetapopulationProblem(p::HeardParameters, habitat::AbstractMatrix{Bool};
    as_coded::Bool = true) = HeardMetapopulationProblem(p, BitMatrix(habitat), as_coded)
HeardMetapopulationProblem(p::HeardParameters, landscape::Integer; as_coded::Bool = true) =
    HeardMetapopulationProblem(p, heard2024_landscapes().habitat[landscape]; as_coded)

"""
    HeardMetapopulationSolution

Result of `solve(::HeardMetapopulationProblem)`: `patches` (the `(row, col)` of each habitat
cell), `subadults` and `adults` (`(time_steps + 1) × n_patches` matrices of females, row 1
the initial state), `total` (females per year) and `extinct` (no female left at the end).
"""
struct HeardMetapopulationSolution
    patches::Vector{Tuple{Int,Int}}
    subadults::Matrix{Int}
    adults::Matrix{Int}
    total::Vector{Int}
    extinct::Bool
end

# Destinations of dispersers from each patch: (patch indices, kernel weights). Follows
# steps' get_distance_info / get_ids_dists: offsets of whole cells, distance < max_distance,
# the origin included (distance 0).
function _heard_destinations(prob::HeardMetapopulationProblem, patches)
    p = prob.p
    nr, nc = size(prob.habitat)
    index = Dict(rc => i for (i, rc) in enumerate(patches))
    out = Vector{Tuple{Vector{Int},Vector{Float64}}}(undef, length(patches))
    for (i, (r, c)) in enumerate(patches)
        # steps' id2coord: row-major id; a last-column cell maps to column 0 of the next row
        r0, c0 = (prob.as_coded && c == nc) ? (r + 1, 0) : (r, c)
        ids = Int[]
        w = Float64[]
        for (j, (r2, c2)) in enumerate(patches)
            d = p.cell_size * hypot(r2 - r0, c2 - c0)
            d < p.max_distance || continue
            push!(ids, j)
            push!(w, exp(-d / p.distance_decay))
        end
        out[i] = (ids, w ./ sum(w))
    end
    return out
end

# Multinomial draw of n over probabilities ∝ w (sequential binomials).
function _heard_multinomial!(out::AbstractVector{Int}, rng::AbstractRNG, n::Integer,
        w::AbstractVector{<:Real})
    fill!(out, 0)
    rest = float(sum(w))
    left = Int(n)
    for k in eachindex(w)
        left == 0 && break
        if k == lastindex(w) || rest <= w[k]
            out[k] = left
            break
        end
        q = clamp(w[k] / rest, 0.0, 1.0)
        x = q > 0 ? rand(rng, Binomial(left, q)) : 0
        out[k] = x
        left -= x
        rest -= w[k]
    end
    return out
end

function CommonSolve.solve(prob::HeardMetapopulationProblem;
        rng::AbstractRNG = Random.default_rng())
    p = prob.p
    nr, nc = size(prob.habitat)
    patches = [(r, c) for r in 1:nr for c in 1:nc if prob.habitat[r, c]]   # row-major, as R
    np = length(patches)
    T = p.time_steps
    sub = zeros(Int, T + 1, np)
    ad = zeros(Int, T + 1, np)
    total = zeros(Int, T + 1)
    n0 = round(p.carrying_capacity * 0.5) * p.initial_fraction  # R round: half to even
    ns = fill(Float64(n0), np)
    na = fill(Float64(n0), np)
    sub[1, :] .= round.(Int, ns)
    ad[1, :] .= round.(Int, na)
    total[1] = sum(sub[1, :]) + sum(ad[1, :])
    A = heard_projection_matrix(p)
    sdF = (A[1, 1] * p.mets_cv, A[1, 2] * p.mets_cv)
    sdS = heard_survival_sd(p.adult_survival, p.mortality_cv)
    dest = _heard_destinations(prob, patches)
    buf = zeros(Int, np)
    pool = zeros(Int, 2np)
    for t in 1:T
        # 1. growth: global noise on each non-zero entry
        Fs = clamp(A[1, 1] + sdF[1] * randn(rng), 0.0, 1.0)     # diagonal: steps caps at 1
        Ss = clamp(A[2, 1] + sdS * randn(rng), 0.0, 1.0)
        Fa = max(A[1, 2] + sdF[2] * randn(rng), 0.0)
        Sa = clamp(A[2, 2] + sdS * randn(rng), 0.0, 1.0)
        for i in 1:np
            s, a = round(Int, ns[i]), round(Int, na[i])
            s + a == 0 && continue
            off = (Fs * s > 0 ? rand(rng, Poisson(Fs * s)) : 0) +
                  (Fa * a > 0 ? rand(rng, Poisson(Fa * a)) : 0)
            surv = rand(rng, Binomial(s, Ss)) + rand(rng, Binomial(a, Sa))
            ns[i], na[i] = off, surv
        end
        # 2. dispersal (origin numbers from before any arrivals)
        if p.prop_dispersing > 0
            for (stage, x) in ((1, ns), (2, na))
                orig = round.(Int, x)
                new = copy(orig)
                for i in 1:np
                    orig[i] > 0 || continue
                    nd = round(Int, orig[i] * p.prop_dispersing)   # half to even, as R
                    nd == 0 && continue
                    ids, w = dest[i]
                    new[i] -= nd
                    _heard_multinomial!(view(buf, 1:length(ids)), rng, nd, w)
                    for (k, j) in enumerate(ids)
                        new[j] += buf[k]
                    end
                end
                x .= new
            end
        end
        # 3. ceiling density dependence and rounding
        for i in 1:np
            tot = ns[i] + na[i]
            if tot > p.carrying_capacity
                f = p.carrying_capacity / tot
                ns[i] *= f
                na[i] *= f
            end
        end
        if prob.as_coded
            w = vcat(ns, na)
            N = round(Int, sum(w))
            if N > 0
                _heard_multinomial!(pool, rng, N, w)
                ns .= pool[1:np]
                na .= pool[np+1:end]
            end
        else
            ns .= round.(ns)
            na .= round.(na)
        end
        sub[t+1, :] .= round.(Int, ns)
        ad[t+1, :] .= round.(Int, na)
        total[t+1] = sum(sub[t+1, :]) + sum(ad[t+1, :])
    end
    return HeardMetapopulationSolution(patches, sub, ad, total, total[end] == 0)
end

"""
    heard_persistence(p::HeardParameters; landscapes = 1:20, replicates = 1,
                      as_coded = true, rng = Random.default_rng()) -> Float64

Fraction of simulations in which the metapopulation persists for `p.time_steps` years,
over the given landscapes ([`heard2024_landscapes`](@ref)) with `replicates` runs each (the
paper runs one per scenario and landscape).
"""
function heard_persistence(p::HeardParameters; landscapes = 1:20, replicates::Integer = 1,
        as_coded::Bool = true, rng::AbstractRNG = Random.default_rng())
    L = heard2024_landscapes().habitat
    n = 0
    k = 0
    for l in landscapes, _ in 1:replicates
        prob = HeardMetapopulationProblem(p, L[l]; as_coded)
        n += !solve(prob; rng).extinct
        k += 1
    end
    return n / k
end

"""
    heard_persistence_grid(; species = :raniformis, epoch = :contemporary, mets = nothing,
                           mets_cv = nothing, mortality_cv = nothing, landscapes = 1:20,
                           as_coded = true, rng = Random.default_rng())

Fraction of the paper's scenarios (`mets` × `mets_cv` × `mortality_cv` × landscapes, one
simulation each) in which the metapopulation persists. Any of the three grids can be
replaced by a subset; `nothing` uses the full grid of [`heard_scenario_grid`](@ref).
"""
function heard_persistence_grid(; species::Symbol = :raniformis,
        epoch::Symbol = :contemporary, mets = nothing, mets_cv = nothing,
        mortality_cv = nothing, landscapes = 1:20, as_coded::Bool = true,
        rng::AbstractRNG = Random.default_rng())
    g = heard_scenario_grid()
    ms = mets === nothing ? g.mets : mets
    fs = mets_cv === nothing ? g.mets_cv : mets_cv
    vs = mortality_cv === nothing ? g.mortality_cv : mortality_cv
    L = heard2024_landscapes().habitat
    n = 0
    k = 0
    for m in ms, f in fs, v in vs
        p = HeardParameters(:heard2024; species, epoch, mets = m, mets_cv = f,
            mortality_cv = v)
        for l in landscapes
            n += !solve(HeardMetapopulationProblem(p, L[l]; as_coded); rng).extinct
            k += 1
        end
    end
    return n / k
end

# ---- data ------------------------------------------------------------------------------------

"""
    heard2024_path(file::Symbol = :age_structure)

Path of a bundled file of Heard et al. (2024) (Zenodo doi:10.5281/zenodo.10457229,
CC BY 4.0; see `data/heard2024/README.md`): `:age_structure` (skeletochronology ages,
`Age_structure_data.csv`), `:landscapes` (the 20 sample landscapes, converted from
`sample_landscapes.grd`) or `:connectivity` (their connectivity index).
"""
function heard2024_path(file::Symbol = :age_structure)
    names = (age_structure = "age_structure.csv", landscapes = "landscapes.csv",
             connectivity = "landscape_connectivity.csv")
    haskey(names, file) ||
        throw(ArgumentError("file must be one of $(keys(names)); got :$file"))
    return joinpath(pkgdir(@__MODULE__), "data", "heard2024", names[file])
end

"""
    heard2024_landscapes() -> (habitat, connectivity)

The 20 sample landscapes of Heard et al. (2024): `habitat`, a vector of 25 × 25 `BitMatrix`
(row 1 north; `true` = habitat patch, one 65 m cell each), and `connectivity`, the authors'
connectivity index of each landscape (0–85).
"""
function heard2024_landscapes()
    lines = filter(!isempty ∘ strip, readlines(heard2024_path(:landscapes)))
    rows = [parse.(Int, split(l, ',')) for l in lines[2:end]]
    nl = maximum(r[1] for r in rows)
    nr = maximum(r[2] for r in rows)
    nc = length(rows[1]) - 2
    habitat = [falses(nr, nc) for _ in 1:nl]
    for r in rows
        habitat[r[1]][r[2], :] .= r[3:end] .== 1
    end
    conn = Int[]
    for l in readlines(heard2024_path(:connectivity))[2:end]
        isempty(strip(l)) && continue
        push!(conn, parse(Int, split(l, ',')[2]))
    end
    return (habitat = habitat, connectivity = conn)
end

# ---- reference runs for the validation targets -----------------------------------------------

const _HEARD2024_GRID = Dict{Tuple{Symbol,Symbol},Any}()

# The paper's full scenario grid (20,000 simulations) for one species and epoch, once per
# session with a fixed seed. Loops in the order of heard_persistence_grid (mets, mets_cv,
# mortality_cv, landscape) with one generator. Returns every run (`runs`: mets, mets_cv,
# mortality_cv, landscape, persist), the overall persistence fraction and the fraction per
# `mets` value.
function _heard2024_reference_grid(species::Symbol, epoch::Symbol)
    get!(_HEARD2024_GRID, (species, epoch)) do
        g = heard_scenario_grid()
        L = heard2024_landscapes().habitat
        rng = Random.Xoshiro(2024)
        runs = NamedTuple{(:mets, :mets_cv, :mortality_cv, :landscape, :persist),
                          Tuple{Float64,Float64,Float64,Int,Bool}}[]
        for m in g.mets, f in g.mets_cv, v in g.mortality_cv
            p = HeardParameters(:heard2024; species, epoch, mets = m, mets_cv = f,
                mortality_cv = v)
            for l in eachindex(L)
                ok = !solve(HeardMetapopulationProblem(p, L[l]); rng).extinct
                push!(runs, (mets = m, mets_cv = f, mortality_cv = v, landscape = l,
                             persist = ok))
            end
        end
        by = [count(r -> r.mets == m && r.persist, runs) / count(r -> r.mets == m, runs)
              for m in g.mets]
        (runs = runs, overall = count(r -> r.persist, runs) / length(runs), mets = g.mets,
         by_mets = by)
    end
end

# Logistic regression of Heard et al. (2024) Table 1 (Metapop_sims_analysis_GLM.R): post-Bd
# persistence of the scenarios that persisted pre-Bd, on fertility, CV of fertility, CV of
# mortality and landscape connectivity, each scaled to mean 0 and SD 1 (R's scale()), with
# all two-way interactions. Returns the 11 coefficients in the order of Table 1. Pre- and
# post-Bd runs are paired by scenario and landscape (the authors match on connectivity,
# which pairs landscapes 18 and 19, both 77, with landscape 18; a negligible difference).
function _heard2024_glm(species::Symbol)
    pre = _heard2024_reference_grid(species, :historical).runs
    post = _heard2024_reference_grid(species, :contemporary).runs
    conn = heard2024_landscapes().connectivity
    keep = [i for i in eachindex(pre) if pre[i].persist]
    z(v) = (v .- mean(v)) ./ sqrt(var(v))
    F = z([post[i].mets for i in keep])
    Cf = z([post[i].mets_cv for i in keep])
    Cm = z([post[i].mortality_cv for i in keep])
    Co = z(Float64[conn[post[i].landscape] for i in keep])
    y = Float64[post[i].persist for i in keep]
    X = hcat(ones(length(y)), F, Cf, Cm, Co, F .* Cf, F .* Cm, F .* Co, Cf .* Cm, Cf .* Co,
             Cm .* Co)
    β = zeros(size(X, 2))
    for _ in 1:50                                   # IRLS (Newton) for the logit link
        η = X * β
        μ = 1 ./ (1 .+ exp.(-η))
        w = max.(μ .* (1 .- μ), 1e-12)
        step = (X' * (w .* X)) \ (X' * (y .- μ))
        β .+= step
        maximum(abs, step) < 1e-10 && break
    end
    names = (:intercept, :fertility, :cv_fertility, :cv_mortality, :connectivity,
             :fertility_cv_fertility, :fertility_cv_mortality, :fertility_connectivity,
             :cv_fertility_cv_mortality, :cv_fertility_connectivity, :cv_mortality_connectivity)
    return NamedTuple{names}(Tuple(β))
end

# Smallest fecundity (metamorphs per female) at which at least half the scenarios persist.
function _heard2024_mets_half(species::Symbol, epoch::Symbol)
    r = _heard2024_reference_grid(species, epoch)
    i = findfirst(>=(0.5), r.by_mets)
    return i === nothing ? Inf : r.mets[i]
end
