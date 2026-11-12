# Daversa, Bosch, Manica, Garner & Fenton (2022) Am Nat 200(4). doi:10.1086/720638
# Host species share one zoospore pool in a mesocosm (salamander, toad and tree frog larvae).
# Units: time in days; hosts and zoospores are numbers per tank (80 L), not densities.
# Wiki: lab/wiki/models/bd-daversa2022-multihost.md
#
# Model (p. 5, checked against the PDF and against the authors' R code, Zenodo 5068208):
#   dH_i/dt = -μ_i H_i
#   dI_i/dt = β_i (H_i - I_i) Z - I_i (μ_i + σ_i)
#   dZ/dt   = Σ_k λ_k I_k + Φ(t) - Z (γ + Σ_k H_k β_k)
# H_i are all hosts of species i, I_i are the infected ones, Φ(t) are pulse inoculations, and
# the term Z Σ H_k β_k is the loss of zoospores to contacts with hosts.
#
# R0 (p. 6): R0_i = β_i λ_i H_i / (μ_i (γ + Σ_k H_k β_k)) as printed, and R0_TOT = Σ_i R0_i.
# The authors' code (R0calc) divides by (μ_i + σ_i), the mean time infected, and so does the
# next-generation matrix of the equations above. The printed formula has dropped σ_i, and
# with it the Fig. 4b values cannot be reproduced. `species_reproduction_numbers` uses μ_i + σ_i
# by default, and `as_printed = true` gives the printed form.

"""
    DaversaParameters(; species, β, λ, σ, γ)
    DaversaParameters(:daversa2022; set = "")

Parameters of the Daversa et al. (2022) shared-pool model for any number of host species.

- `species`: the species names.
- `β`: transmission rates, in zoospore⁻¹ day⁻¹. This is the infection rate per susceptible
  host per zoospore in the tank, and also each host's rate of removing zoospores.
- `λ`: shedding rates, in zoospores per infected host per day.
- `σ`: recovery rates, in day⁻¹. Recovered hosts become susceptible again.
- `γ`: zoospore mortality, in day⁻¹.

The keyed constructor reads the three species (salamander, toad, frog) and `γ` from the
compendium. Its β, λ and σ are posterior medians digitised from Fig. 4a: the paper prints no
numbers, so treat them as approximate. Host loss rates `μ` vary from tank to tank and belong
to the problem, not to these parameters.
"""
struct DaversaParameters
    species::Vector{Symbol}
    β::Vector{Float64}
    λ::Vector{Float64}
    σ::Vector{Float64}
    γ::Float64
    function DaversaParameters(species, β, λ, σ, γ)
        n = length(species)
        length(β) == length(λ) == length(σ) == n ||
            throw(ArgumentError("species, β, λ and σ must have the same length"))
        all(>=(0), β) && all(>=(0), λ) && all(>=(0), σ) && γ >= 0 ||
            throw(ArgumentError("rates must be non-negative"))
        return new(collect(Symbol, species), collect(Float64, β), collect(Float64, λ),
            collect(Float64, σ), Float64(γ))
    end
end
DaversaParameters(; species, β, λ, σ, γ) = DaversaParameters(species, β, λ, σ, γ)

const DAVERSA_SPECIES = (:salamander, :toad, :frog)
const DAVERSA_CODES = ("S", "T", "F")

# A rate from the compendium, converted to day⁻¹ from its recorded units (a named set
# overrides the base set).
function _daversa_rate(key::Symbol, set::AbstractString, sym::AbstractString)
    recs = filter(r -> r.symbol === Symbol(sym), parameter_records(key))
    i = something(findfirst(r -> !isempty(set) && r.set == set, recs),
        findfirst(r -> r.set == "", recs))
    return in_time_unit(recs[i], :day)
end

function DaversaParameters(key::Symbol; set::AbstractString = "")
    isempty(set) || any(r -> r.set == set, parameter_records(key)) ||
        throw(ArgumentError("$key has no parameter set \"$set\""))
    r(s) = _daversa_rate(key, set, s)
    return DaversaParameters(; species = collect(DAVERSA_SPECIES),
        β = [r("β_" * c) for c in DAVERSA_CODES], λ = [r("λ_" * c) for c in DAVERSA_CODES],
        σ = [r("σ_" * c) for c in DAVERSA_CODES], γ = r("γ"))
end

nspecies(p::DaversaParameters) = length(p.species)

"""In-place right-hand side; state `u = [H₁…Hₙ, I₁…Iₙ, Z]`, parameters `(p, μ)`."""
function daversa_rhs!(du, u, pμ, t)
    p, μ = pμ
    n = nspecies(p)
    Z = u[2n + 1]
    removal = p.γ
    shed = 0.0
    @inbounds for i in 1:n
        H, I = u[i], u[n + i]
        du[i] = -μ[i] * H
        du[n + i] = p.β[i] * (H - I) * Z - I * (μ[i] + p.σ[i])
        removal += H * p.β[i]
        shed += p.λ[i] * I
    end
    du[2n + 1] = shed - Z * removal
    return nothing
end

"""
    MultiHostPoolProblem(p::DaversaParameters, H0; I0 = zeros, Z0 = 0, μ = zeros,
                         tspan = (0.0, 30.0), doses = Tuple{Float64,Float64}[])

The Daversa et al. (2022) model for one tank: `n = length(p.species)` host species sharing
one zoospore pool.

- `H0`: initial numbers of hosts.
- `I0`: initial numbers infected.
- `Z0`: initial zoospores.
- `μ`: per-species host loss rates, in day⁻¹ (death plus metamorphosis). The paper estimates
  them per tank as `ln(H(0)/H(T))/T`.
- `doses`: `(day, zoospores)` pulses, which make up the inoculation term Φ(t).

The state is `[H₁…Hₙ, I₁…Iₙ, Z]`. `solve(prob)` returns the `ODESolution`, and
`to_ode_problem(prob)` the `ODEProblem` with the dosing callback attached. Use
[`daversa_doses`](@ref) for the inoculation schedule of the published fit.
"""
struct MultiHostPoolProblem
    p::DaversaParameters
    H0::Vector{Float64}
    I0::Vector{Float64}
    Z0::Float64
    μ::Vector{Float64}
    tspan::Tuple{Float64,Float64}
    doses::Vector{Tuple{Float64,Float64}}
end

function MultiHostPoolProblem(p::DaversaParameters, H0::AbstractVector;
        I0::AbstractVector = zeros(length(H0)), Z0::Real = 0.0,
        μ::AbstractVector = zeros(length(H0)), tspan = (0.0, 30.0),
        doses = Tuple{Float64,Float64}[])
    n = nspecies(p)
    length(H0) == length(I0) == length(μ) == n ||
        throw(ArgumentError("H0, I0 and μ need one entry per species ($n)"))
    all(0 .<= I0 .<= H0) || throw(ArgumentError("need 0 ≤ I0 ≤ H0"))
    return MultiHostPoolProblem(p, collect(Float64, H0), collect(Float64, I0), Float64(Z0),
        collect(Float64, μ), (Float64(tspan[1]), Float64(tspan[2])),
        [(Float64(t), Float64(z)) for (t, z) in doses])
end

function to_ode_problem(prob::MultiHostPoolProblem)
    u0 = vcat(prob.H0, prob.I0, prob.Z0)
    times = [t for (t, _) in prob.doses if prob.tspan[1] < t <= prob.tspan[2]]
    iz = length(u0)
    if isempty(times)
        return ODEProblem(daversa_rhs!, u0, prob.tspan, (prob.p, prob.μ))
    end
    function dose!(integ)
        for (t, z) in prob.doses
            t == integ.t && (integ.u[iz] += z)
        end
        return nothing
    end
    cb = DiscreteCallback((u, t, integ) -> t in times, dose!; save_positions = (true, true))
    return ODEProblem(daversa_rhs!, u0, prob.tspan, (prob.p, prob.μ); callback = cb,
        tstops = times)
end

function CommonSolve.solve(prob::MultiHostPoolProblem, alg = Tsit5(); abstol = 1e-8,
        reltol = 1e-8, kwargs...)
    return solve(to_ode_problem(prob), alg; abstol, reltol, kwargs...)
end

"""
    daversa_doses(key = :daversa2022) -> (; Z0, doses)

The inoculations used in the authors' fit, from the compendium: `Z0` zoospores at day 0,
then `doses` as `(day, zoospores)` pulses on days 10, 17 and 25.

These come from the authors' R code. The paper text instead says three inoculations, on
days 1, 11 and 18, of 18,000 to 24,000 zoospores L⁻¹ in 80 L.
"""
function daversa_doses(key::Symbol = :daversa2022)
    ps = parameter_set(key)
    return (; Z0 = ps.Z0, doses = [(ps.dose_t1, ps.dose_n1), (ps.dose_t2, ps.dose_n2),
        (ps.dose_t3, ps.dose_n3)])
end

# --- Reproduction numbers ----------------------------------------------------------------

"""
    species_reproduction_numbers(p::DaversaParameters, H; μ = zeros, as_printed = false)

The contribution of each species to Bd's basic reproduction number in a tank with host
numbers `H` and loss rates `μ` (day⁻¹):

    R0_i = β_i λ_i H_i / ((μ_i + σ_i) (γ + Σ_k H_k β_k))

This is the formula in the authors' code. It is also what the next-generation matrix of the
model gives (see [`next_generation_matrix`](@ref)). Each term multiplies three things: the
infections a zoospore causes in species i during its lifetime, the time an infected host stays
infected, and its shedding.

The formula printed on p. 6 has `μ_i` alone in the denominator, which drops recovery. Pass
`as_printed = true` to use it. It overstates R0 about twentyfold when σ ≈ 0.7 day⁻¹.

The total `R0_TOT` is the sum of the contributions ([`basic_reproduction_number`](@ref)).
"""
function species_reproduction_numbers(p::DaversaParameters, H::AbstractVector;
        μ::AbstractVector = zeros(length(H)), as_printed::Bool = false)
    n = nspecies(p)
    length(H) == length(μ) == n || throw(ArgumentError("H and μ need one entry per species"))
    pool = p.γ + sum(H[k] * p.β[k] for k in 1:n)
    return [H[i] == 0 ? 0.0 :
            p.β[i] * p.λ[i] * H[i] / ((as_printed ? μ[i] : μ[i] + p.σ[i]) * pool) for i in 1:n]
end

"""
    basic_reproduction_number(p::DaversaParameters, H; μ = zeros, as_printed = false)

`R0_TOT = Σ_i R0_i`, from [`species_reproduction_numbers`](@ref). Daversa et al. (2022)
define the total as this sum. Every transmission passes through the one pool, so the
host-to-host next-generation matrix has rank one. Its spectral radius is therefore its trace,
which is the same sum.
"""
basic_reproduction_number(p::DaversaParameters, H::AbstractVector; kwargs...) =
    sum(species_reproduction_numbers(p, H; kwargs...))

"""
    next_generation_matrix(p::DaversaParameters, H; μ = zeros, spores = :transition)
        -> (; K, R0)

The next-generation matrix `K = F V⁻¹` of the model at the disease-free state, computed
numerically from the right-hand side. It does not use the R0 formula. The infected
compartments are `(I₁…Iₙ, Z)`, and host numbers are held at `H`. The Jacobian of the
right-hand side is taken by central differences, which are exact here because the model is
bilinear. `F` is the Jacobian of the new-infection terms alone, and `V = F − J`.

- `spores = :transition`: only host infections are new; shedding moves an infection from a
  host into the pool. Then `R0 = ρ(K)` equals `Σ_i R0_i`.
- `spores = :generation`: shedding also counts as a new generation (host → spore → host).
  Then `R0 = ρ(K) = √(Σ_i R0_i)`. Both give the same threshold at 1.
"""
function next_generation_matrix(p::DaversaParameters, H::AbstractVector;
        μ::AbstractVector = zeros(length(H)), spores::Symbol = :transition)
    spores in (:transition, :generation) ||
        throw(ArgumentError("spores must be :transition or :generation"))
    n = nspecies(p)
    x0 = vcat(Float64.(H), zeros(n), 0.0)
    inf = (n + 1):(2n + 1)
    rhs(u) = (du = similar(u); daversa_rhs!(du, u, (p, μ), 0.0); du)
    function newinf(u)
        f = zeros(length(u))
        for i in 1:n
            f[n + i] = p.β[i] * (u[i] - u[n + i]) * u[2n + 1]
        end
        spores === :generation && (f[2n + 1] = sum(p.λ[i] * u[n + i] for i in 1:n))
        return f
    end
    J = _infected_jacobian(rhs, x0, inf)
    F = _infected_jacobian(newinf, x0, inf)
    V = F .- J
    K = F / V
    return (; K, R0 = maximum(abs, eigvals(K)))
end

# Central-difference Jacobian of f at x0, restricted to the rows and columns in idx
function _infected_jacobian(f, x0::AbstractVector, idx)
    M = zeros(length(idx), length(idx))
    for (c, j) in enumerate(idx)
        e = zeros(length(x0))
        e[j] = 1.0
        M[:, c] = (f(x0 .+ e) .- f(x0 .- e))[idx] ./ 2
    end
    return M
end

# --- The mesocosm data (Dryad 10.5061/dryad.2v6wwpznk) -------------------------------------

"""
    daversa2022_path(file = "")

Path to the Daversa et al. (2022) data shipped with the package (`data/daversa2022/`).
"""
daversa2022_path(file::AbstractString = "") =
    normpath(joinpath(@__DIR__, "..", "..", "data", "daversa2022", file))

"""
    MesocosmTank

One tank of the Daversa et al. (2022) experiment, prepared as in Part 4 of the authors' R
code. Fields:

- `tank` and `treatment`.
- `H0` and `Hfinal`: the numbers of salamanders, toads and frogs at the start and at the end.
- `day_final`: the final sampling day, counted from 2017-07-04 (day 0).
- `prev_initial` and `prev_final`: the prevalence among sampled animals (`bd.load > 0`). It is
  `NaN` where a species was not sampled; toads and frogs were sampled only at the end.
- `μ`: loss rates in day⁻¹, `ln(H0/Hfinal)/day_final`, and 0 where a species is absent.
"""
struct MesocosmTank
    tank::String
    treatment::String
    H0::Vector{Float64}
    Hfinal::Vector{Float64}
    day_final::Float64
    prev_initial::Vector{Float64}
    prev_final::Vector{Float64}
    μ::Vector{Float64}
end

function _read_csv(path)
    lines = filter(!isempty, strip.(readlines(path), '\r'))
    header = split(lines[1], ',')
    return header, [split(l, ',') for l in lines[2:end]]
end

"""
    load_daversa2022(; three_species_only = false) -> Vector{MesocosmTank}

Read the 36 tanks, or with `three_species_only = true` the 21 three-species tanks used in the
fit. The data preparation follows Part 4 of the authors' R code:

- drop the `accidental` samples;
- count day 0 from the first sampling date;
- call an animal infected when its load is above 0;
- take each tank's final day as the earliest final sampling date in that tank;
- compute exponential loss rates from the initial and final numbers.
"""
function load_daversa2022(; three_species_only::Bool = false)
    _, des = _read_csv(daversa2022_path("design.csv"))
    _, ends = _read_csv(daversa2022_path("endnums.csv"))
    hdr, inf = _read_csv(daversa2022_path("Multihost_mesocosms_2017_infection.csv"))
    col(name) = findfirst(==(name), hdr)
    cdate, ctank, csp, cphase, cload = col("date"), col("tank"), col("sp"),
        col("timepoint.name"), col("bd.load")
    rows = filter(r -> r[cphase] != "accidental", inf)
    date(r) = Date(r[cdate], dateformat"dd/mm/yyyy")
    d0 = minimum(date, rows)
    codes = ("Ss", "Bs", "Hm")
    endmap = Dict(String(r[1]) => parse.(Float64, r[2:4]) for r in ends)
    tanks = MesocosmTank[]
    for r in des
        id = startswith(r[1], "X") ? String(r[1][2:end]) : String(r[1])
        treatment = String(r[5])
        three_species_only && !(treatment in ("random", "Ss + Bs + Hm")) && continue
        H0 = parse.(Float64, r[2:4])
        Hf = endmap[id]
        mine = filter(x -> x[ctank] == id, rows)
        finals = filter(x -> x[cphase] == "final", mine)
        day_final = Float64(minimum(x -> (date(x) - d0).value, finals))
        prev(phase, sp) = begin
            xs = filter(x -> x[cphase] == phase && x[csp] == sp, mine)
            isempty(xs) ? NaN : count(x -> parse(Float64, x[cload]) > 0, xs) / length(xs)
        end
        μ = [H0[i] > 0 ? log(H0[i] / Hf[i]) / day_final : 0.0 for i in 1:3]
        push!(tanks, MesocosmTank(id, treatment, H0, Hf, day_final,
            [prev("initial", s) for s in codes], [prev("final", s) for s in codes], μ))
    end
    return tanks
end

"""
    MultiHostPoolProblem(p::DaversaParameters, tank::MesocosmTank; Z0, doses)

The problem for one experimental tank, set up as in the authors' fit:

- `H0` and the loss rates `μ` come from the tank;
- `I0` is the initial prevalence times `H0`, and 0 where a species was not sampled;
- the run ends on the tank's final sampling day;
- the inoculations default to [`daversa_doses`](@ref).
"""
function MultiHostPoolProblem(p::DaversaParameters, tank::MesocosmTank;
        Z0::Real = daversa_doses().Z0, doses = daversa_doses().doses)
    I0 = [isnan(q) ? 0.0 : q * h for (q, h) in zip(tank.prev_initial, tank.H0)]
    return MultiHostPoolProblem(p, tank.H0; I0, Z0, μ = tank.μ,
        tspan = (0.0, tank.day_final), doses)
end

"""
    observed_infected(tank::MesocosmTank) -> Vector

The final numbers infected, as the authors estimate them: final prevalence times final
abundance. This is `NaN` where a species was absent or not sampled.
"""
observed_infected(tank::MesocosmTank) = tank.prev_final .* tank.Hfinal

"""
    predicted_infected(p::DaversaParameters, tank::MesocosmTank; kwargs...) -> Vector

The numbers infected that the model predicts on the tank's final sampling day.
"""
function predicted_infected(p::DaversaParameters, tank::MesocosmTank; kwargs...)
    sol = solve(MultiHostPoolProblem(p, tank; kwargs...); save_everystep = false)
    SciMLBase.successful_retcode(sol) || error("Daversa model: solver failed in tank $(tank.tank)")
    n = nspecies(p)
    return sol.u[end][(n + 1):(2n)]
end

"""
    daversa_loglikelihood(p::DaversaParameters, tanks) -> Float64

The authors' log-likelihood (R code, `tankLL`): a binomial over the final numbers infected,
summed over species and tanks. It is written exactly as in the code:

- the number of trials is the initial abundance `H0`, not the final one;
- the success probability is `I(T)/H0`;
- `1 − p` is floored at 0.01;
- the observed counts may be fractional, because they are prevalence × final abundance.

The binomial coefficient is omitted.
"""
function daversa_loglikelihood(p::DaversaParameters, tanks::AbstractVector{MesocosmTank})
    ll = 0.0
    for tank in tanks
        pred = predicted_infected(p, tank)
        obs = observed_infected(tank)
        for i in eachindex(pred)
            (tank.H0[i] > 0 && !isnan(obs[i])) || continue
            q = pred[i] / tank.H0[i]
            ll += obs[i] * log(q) + (tank.H0[i] - obs[i]) * log(1 - min(q, 0.99))
        end
    end
    return ll
end

"""
    tank_reproduction_numbers(p::DaversaParameters,
                              tanks = load_daversa2022(; three_species_only = true);
                              μ = :code, as_printed = false) -> Vector{Vector{Float64}}

The species contributions `R0_i` ([`species_reproduction_numbers`](@ref)) for each tank, at
the tank's initial host numbers. The loss rates are set by `μ`:

- `μ = :code`: as in the authors' `R0calc`. Salamanders use the tank's own loss rate. Toads and
  frogs use the mean loss rate over all 36 tanks, which counts as 0 the tanks where the species
  is absent.
- `μ = :tank`: each species uses its own loss rate in the tank.
"""
function tank_reproduction_numbers(p::DaversaParameters,
        tanks::AbstractVector{MesocosmTank} = load_daversa2022(; three_species_only = true);
        μ::Symbol = :code, as_printed::Bool = false)
    μ in (:code, :tank) || throw(ArgumentError("μ must be :code or :tank"))
    all36 = μ === :code ? load_daversa2022() : tanks
    meanμ = [sum(t.μ[i] for t in all36) / length(all36) for i in 1:3]
    return [species_reproduction_numbers(p, t.H0;
        μ = μ === :code ? [t.μ[1], meanμ[2], meanμ[3]] : t.μ, as_printed) for t in tanks]
end

# Per-tank posterior medians of R0_TOT digitised from Fig. 4d (data/daversa2022/README.md)
function _daversa_fig4_medians()
    _, rows = _read_csv(daversa2022_path("fig4cd_R0tot_digitised.csv"))
    return Dict(String(r[1]) => parse(Float64, r[2]) for r in rows)
end

_median(x) = (s = sort(collect(x)); n = length(s);
    isodd(n) ? s[(n + 1) ÷ 2] : (s[n ÷ 2] + s[n ÷ 2 + 1]) / 2)

function _pearson(x, y)
    mx, my = sum(x) / length(x), sum(y) / length(y)
    return sum((x .- mx) .* (y .- my)) / sqrt(sum(abs2, x .- mx) * sum(abs2, y .- my))
end

# Correlation between observed and predicted final numbers infected over the 63 data points
# of the 21 three-species tanks (the paper's goodness of fit, Fig. S7)
function _daversa_fit_correlation(p::DaversaParameters)
    tanks = load_daversa2022(; three_species_only = true)
    obs = reduce(vcat, observed_infected.(tanks))
    pred = reduce(vcat, [predicted_infected(p, t) for t in tanks])
    ok = .!isnan.(obs)
    return _pearson(pred[ok], obs[ok])
end
