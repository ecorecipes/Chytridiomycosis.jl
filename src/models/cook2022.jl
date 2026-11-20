# Cook, Pope, Cummings & Piovia-Scott (2022) In situ treatment of juvenile frogs for disease
# can reverse population declines. Conserv Sci Pract 4:e12762, doi:10.1111/csp2.12762
# (CC BY 4.0). Stage-structured stochastic projection of Cascades frogs (Rana cascadae) with
# itraconazole treatment of metamorphs. Model read from the Supporting Information
# ("Population simulations", Figure S3, Table S2); values in data/compendium/cook2022.toml.
#
# Stages: metamorph M, young of the previous year Y, subadult S, adult A. One step is a year:
#   Y' ~ Binomial(M_untreated, SU) + Binomial(M_treated, ST)
#   S' ~ Binomial(Y, SY) + (subadults surviving, SS, that do not mature)
#   A' ~ Binomial(A, SA) + (subadults surviving, SS, that mature, GS)
#   M' ~ Poisson(F × females),  females ~ Binomial(A', PF)
# with SA = SS drawn each year from 28 Jolly-Seber estimates, SY = 0.80 SA and F drawn from 28
# fecundities (environmental stochasticity), and binomial demographic stochasticity
# ("selecting the number of individuals surviving/transitioning from a binomial
# distribution"). Start: 1500 metamorphs, 50 YOPY, 50 subadults, 100 adults; 10 years; 1000
# simulations. λ of a simulation is the geometric mean of its annual adult growth rates,
# (A_10 / A_0)^(1/10); the paper reports the mean over simulations.
#
# Choices where the SI is silent (documented, with switches):
#   - Breeding (`breeding = :after_survival`): the metamorphs of year t + 1 are produced by
#     the adults present in year t + 1 (after overwinter survival and maturation), as the
#     adult → metamorph arrow of Fig. S3 suggests for a census at metamorphosis.
#     `:before_survival` uses the adults of year t.
#   - Subadults survive with SS, and survivors mature with probability GS (Fig. S3 draws a
#     survival loop and a transition arrow; GS is called a "transition rate").
#   - SA and F are drawn as a pair (the same one of the 28 population-years; `paired = true`);
#     the SI says only that both are "randomly sampled from a pool".
#   - The number of females is Binomial(A, PF); the SI says "n = the number of adult females".
#   - SU and ST are "sampled from the posterior distribution" of the CJS model, which is not
#     published. The default uses the posterior means of Table S2 (0.12, 0.49) every year;
#     pass a `Distributions` distribution to sample a value each year instead.
#   - Treatment frequency (`schedule = :fixed`): treatment is done in round(f × years) years
#     chosen at random in each simulation; `:random` treats each year with probability f. In a
#     treatment year, each metamorph is treated with probability equal to the treated fraction.
#   - A simulation whose adults die out has λ = 0 (the geometric mean of annual rates with a
#     zero); this happens in fewer than 1 in 1000 simulations without treatment.

"""
    CookParameters(; SU, ST, adult_survival, fecundity, yopy_ratio, growth, female_fraction,
                   paired = true)
    CookParameters(:cook2022; paired = true)

Parameters of the Cascades frog projection of Cook et al. (2022), SI Table S2:

- `SU`, `ST`: overwinter survival of untreated and itraconazole-treated metamorphs. Each is a
  number, or a `Distributions` distribution sampled once a year.
- `adult_survival`: the pool of annual adult (and subadult) survival probabilities, `SA = SS`.
- `fecundity`: the pool of metamorphs produced per adult female, `F`.
- `yopy_ratio`: YOPY survival relative to adult survival (`SY = yopy_ratio × SA`).
- `growth`: annual probability that a surviving subadult becomes an adult (`GS`).
- `female_fraction`: proportion of adults that are female (`PF`).
- `paired`: draw `SA` and `F` from the same population-year (true) or independently.

The keyed constructor reads every value from the compendium (`SU` = 0.12 and `ST` = 0.49, the
posterior means used in the projection).
"""
Base.@kwdef struct CookParameters{TU,TT}
    SU::TU
    ST::TT
    adult_survival::Vector{Float64}
    fecundity::Vector{Float64}
    yopy_ratio::Float64
    growth::Float64
    female_fraction::Float64
    paired::Bool = true
end

function CookParameters(key::Symbol; paired::Bool = true)
    p = parameter_set(key)
    pool(prefix) = [getproperty(p, Symbol(prefix, lpad(i, 2, '0'))) for i in 1:28]
    return CookParameters(; SU = p.su_mean, ST = p.st_mean, adult_survival = pool("sa_"),
        fecundity = pool("f_"), yopy_ratio = p.yopy_survival_ratio, growth = p.subadult_growth,
        female_fraction = p.female_fraction, paired)
end

time_unit(::Type{<:CookParameters}) = :year

"""
    CookProjectionProblem(p::CookParameters; treated_fraction = 0, treated_years = 0,
        years = 10, u0 = nothing, breeding = :after_survival, schedule = :fixed)

A stochastic projection of the Cook et al. (2022) model, SI Fig. S3:
metamorph → YOPY → subadult → adult, with adults producing metamorphs.

- `treated_fraction`: the proportion of metamorphs treated in a treatment year.
- `treated_years`: the proportion of years with treatment (both 0 without treatment, both 1
  for "treatment of all metamorphs in each year").
- `u0`: the initial `(metamorph, yopy, subadult, adult)` counts. The default is read from the
  compendium: 1500, 50, 50 and 100.
- `breeding` and `schedule` select the choices described in the source header:
  `:after_survival` or `:before_survival`, and `:fixed` or `:random`.

`solve(prob; nsims = 1000, rng)` runs the simulations and returns a
[`CookProjectionSolution`](@ref).
"""
struct CookProjectionProblem{P<:CookParameters}
    p::P
    treated_fraction::Float64
    treated_years::Float64
    years::Int
    u0::NTuple{4,Int}
    breeding::Symbol
    schedule::Symbol
end

function CookProjectionProblem(p::CookParameters; treated_fraction::Real = 0.0,
        treated_years::Real = 0.0, years::Integer = 10, u0 = nothing,
        breeding::Symbol = :after_survival, schedule::Symbol = :fixed)
    0 <= treated_fraction <= 1 && 0 <= treated_years <= 1 ||
        throw(ArgumentError("treated_fraction and treated_years must lie in [0, 1]"))
    breeding in (:after_survival, :before_survival) ||
        throw(ArgumentError("breeding must be :after_survival or :before_survival"))
    schedule in (:fixed, :random) || throw(ArgumentError("schedule must be :fixed or :random"))
    if u0 === nothing
        c = parameter_set(:cook2022)
        u0 = (c.initial_metamorphs, c.initial_yopy, c.initial_subadults, c.initial_adults)
    end
    return CookProjectionProblem(p, float(treated_fraction), float(treated_years), Int(years),
        Tuple(round.(Int, u0)), breeding, schedule)
end

"""
    CookProjectionSolution

Result of `solve(::CookProjectionProblem)`:

- `problem`;
- `states`: a `(years + 1) × 4 × nsims` array of stage counts (metamorph, YOPY, subadult,
  adult);
- `λ`: the growth rate of each simulation, `(A_T / A_0)^(1/T)`, which is the geometric mean of
  its annual adult growth rates.

Use [`cook_growth_rate`](@ref) for the mean and the proportion of simulations with `λ > 1`.
"""
struct CookProjectionSolution{P}
    problem::P
    states::Array{Int,3}
    λ::Vector{Float64}
end

_cook_draw(x::Real, rng) = float(x)
_cook_draw(d, rng) = rand(rng, d)

function _cook_simulate!(out::AbstractMatrix{Int}, prob::CookProjectionProblem, rng::AbstractRNG)
    p, T = prob.p, prob.years
    M, Y, S, A = prob.u0
    out[1, :] .= (M, Y, S, A)
    treat = if prob.schedule === :fixed
        k = round(Int, prob.treated_years * T)
        v = falses(T); v[Random.randperm(rng, T)[1:k]] .= true; v
    else
        [rand(rng) < prob.treated_years for _ in 1:T]
    end
    n = length(p.adult_survival)
    for t in 1:T
        i = rand(rng, 1:n)
        sa = p.adult_survival[i]
        f = p.fecundity[p.paired ? i : rand(rng, 1:length(p.fecundity))]
        su, st = _cook_draw(p.SU, rng), _cook_draw(p.ST, rng)
        nT = treat[t] ? rand(rng, Binomial(M, prob.treated_fraction)) : 0
        Yn = rand(rng, Binomial(nT, st)) + rand(rng, Binomial(M - nT, su))
        Ssurv = rand(rng, Binomial(S, sa))
        mature = rand(rng, Binomial(Ssurv, p.growth))
        Sn = rand(rng, Binomial(Y, p.yopy_ratio * sa)) + Ssurv - mature
        An = rand(rng, Binomial(A, sa)) + mature
        breeders = prob.breeding === :after_survival ? An : A
        females = rand(rng, Binomial(breeders, p.female_fraction))
        Mn = rand(rng, Poisson(f * females))
        M, Y, S, A = Mn, Yn, Sn, An
        out[t+1, :] .= (M, Y, S, A)
    end
    return out
end

# solve(prob::CookProjectionProblem; nsims = 1000, rng) -> CookProjectionSolution
function CommonSolve.solve(prob::CookProjectionProblem; nsims::Integer = 1000,
                           rng::AbstractRNG = Random.default_rng())
    T = prob.years
    X = zeros(Int, T + 1, 4, nsims)
    for s in 1:nsims
        _cook_simulate!(view(X, :, :, s), prob, rng)
    end
    λ = [(X[end, 4, s] / X[1, 4, s])^(1 / T) for s in 1:nsims]
    return CookProjectionSolution(prob, X, λ)
end

"""
    cook_growth_rate(sol::CookProjectionSolution) -> (mean, positive)

The mean over simulations of the adult growth rate `λ`, as in Cook et al. (2022), Results
3.4 and Fig. 4a, and the proportion of simulations with `λ > 1` (the numbers in the boxes of
Fig. 4b).
"""
cook_growth_rate(sol::CookProjectionSolution) =
    (mean = sum(sol.λ) / length(sol.λ), positive = count(>(1), sol.λ) / length(sol.λ))

"""
    cook_treatment_grid(p::CookParameters; fractions = 0:0.2:1, frequencies = 0:0.2:1,
        nsims = 1000, rng = Random.default_rng(), kwargs...) -> NamedTuple

Mean `λ` and the proportion of simulations with `λ > 1` over a grid of treated fractions
(rows) and treatment frequencies (columns), as in Cook et al. (2022) Fig. 4b. Other keywords
are passed to [`CookProjectionProblem`](@ref). Returns `(fractions, frequencies, λ, positive)`.
"""
function cook_treatment_grid(p::CookParameters; fractions = 0:0.2:1, frequencies = 0:0.2:1,
                             nsims::Integer = 1000, rng::AbstractRNG = Random.default_rng(),
                             kwargs...)
    L = zeros(length(fractions), length(frequencies))
    P = similar(L)
    for (i, q) in enumerate(fractions), (j, r) in enumerate(frequencies)
        g = cook_growth_rate(solve(CookProjectionProblem(p; treated_fraction = q,
                                                         treated_years = r, kwargs...);
                                   nsims, rng))
        L[i, j], P[i, j] = g.mean, g.positive
    end
    return (fractions = collect(fractions), frequencies = collect(frequencies), λ = L,
            positive = P)
end

"""
    cook_projection_matrix(p::CookParameters; treated = 0.0, breeding = :after_survival)
        -> Matrix{Float64}

The 4 × 4 projection matrix of mean rates (stages metamorph, YOPY, subadult, adult), with the
mean of the survival and fecundity pools and a fraction `treated` of metamorphs treated each
year. With `breeding = :after_survival`, metamorphs come from the adults after survival and
maturation, so the fecundity terms are `PF E[F SA]` (adults) and `PF E[F SA] GS` (subadults).
`E[F SA]` is the mean product when `paired`, and the product of the means otherwise.
[`population_growth_rate`](@ref), [`sensitivity_matrix`](@ref) and
[`elasticity_matrix`](@ref) apply to it.

The paper reports only the stochastic simulation. The matrix ignores environmental and
demographic variance, so its dominant eigenvalue is larger than the mean simulated `λ`. It is
provided for sensitivity analysis and is not a published result.
"""
function cook_projection_matrix(p::CookParameters; treated::Real = 0.0,
                                breeding::Symbol = :after_survival)
    mean_(x) = x isa Real ? float(x) : Distributions.mean(x)
    sa = sum(p.adult_survival) / length(p.adult_survival)
    f = sum(p.fecundity) / length(p.fecundity)
    fsa = p.paired ? sum(p.fecundity .* p.adult_survival) / length(p.fecundity) : f * sa
    sm = (1 - treated) * mean_(p.SU) + treated * mean_(p.ST)
    g, pf = p.growth, p.female_fraction
    A = zeros(4, 4)
    A[2, 1] = sm
    A[3, 2] = p.yopy_ratio * sa
    A[3, 3] = sa * (1 - g)
    A[4, 3] = sa * g
    A[4, 4] = sa
    if breeding === :after_survival
        A[1, 3] = pf * fsa * g
        A[1, 4] = pf * fsa
    else
        A[1, 4] = pf * f
    end
    return A
end
