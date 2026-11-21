# NEW MODEL (not from a published paper): an integral projection model with two continuous
# traits per frog, Bd load and the load of an applied antifungal bacterium.
#
# x = ln Bd load (zoospore equivalents per swab), as in Wilber et al. (2016);
# y = ln load of Janthinobacterium lividum (genomic equivalents per swab), the violacein-producing
#     bacterium applied in probiotic trials (Harris et al. 2009; Becker et al. 2009; Kueneman et
#     al. 2016; Knapp et al. 2022).
# One step is Δt = 3 days (Wilber et al. 2016). Frogs are in four classes:
#   S0     uninfected, no bacteria (y below detection)        a number
#   S(y)   uninfected, carrying bacteria                      a density on the y mesh
#   I0(x)  infected, no bacteria                              a density on the x mesh
#   I(x,y) infected and carrying bacteria                     a density on the product mesh
# plus Z, zoospores in the pond. Bacteria are "lost" when y falls below the lower edge of the
# y mesh, 0 = ln 1 GE per swab, the detection floor; that mass moves to the bacteria-free class.
#
# Per step, at temperature T, with h(y) = logistic((y − y50)/w) the protection given by the
# bacteria (a saturating, MIC-like effect of the violacein they produce):
#   survival      logit s(x,y) = b00 + b10 x + κs h(y)                       (tolerance)
#   Bd clearance  logit l(x,y) = b02 + b12 x + b22 T + κl h(y)               (resistance)
#   Bd growth     x′ ~ Normal(b01 + b11 x + b21 T − κg h(y), v01 e^{2 c01 x})  (resistance)
#   bacteria      y′ | x′ ~ Normal(y − δΔt − ρ x + r σy/σx (x′ − E x′), σy² (1 − r²))
#   infection     φ (1 − κi h(y)) for an uninfected frog carrying bacteria
# so the bacteria act on Bd through h(y), and Bd acts on the bacteria through ρ (Bd disturbs
# the skin microbiome) and the residual correlation r. A frog that clears Bd keeps its
# bacteria (y′ from the same rule); a newly infected frog starts at x′ ~ G0 of Wilber et al.
# Frogs without bacteria acquire them with probability 1 − exp(−(ε + χ f) Δt) per step, f the
# fraction of live frogs carrying bacteria (from the pond, ε, and from other frogs, χ), at
# y ~ Normal(μc, σc). A probiotic bath on a given day gives a fraction of all frogs a fresh load
# y ~ Normal(μb, σb). With every κ = 0, ρ = 0 and r = 0 the two traits are independent, and the
# Bd side is exactly the Wilber et al. (2016) IPM (tested).
#
# Transmission (`transmission`):
#   :pool       the zoospore pool of Wilber et al. (2017): φ = 1 − exp(−Λ Δt) with Λ from the
#               Wilber2017Parameters transmission function (default β0 ln(Z + 1) + β1 I), and
#               Z′ = ν(T) Z + μA Σ eˣ n(x) over infected frogs. The removal of zoospores by
#               newly infected frogs (ψ in the authors' code) is left out, as their text does.
#   :wilber2016 the regression φ(T) of Wilber et al. (2016), from frogs housed singly (no pool).
#   :constant   a fixed infection probability per step (`φ`), e.g. a cohort in a pond whose
#               zoospores are not modelled (the fit to the Knapp et al. 2022 cohort).
# The Bd vital rates are those of Wilber et al. (2016); the pool parameters come from the
# separate fit of Wilber et al. (2017) to the same system. Combining them is our choice.
#
# Parameter values. The Bd side and the pool are published (compendium wilber2016, wilber2017).
# The bacterial decline δ = 0.1704 d⁻¹ (half-life 4.07 days), σy = 0.29 and the post-bath
# distribution (ln GE 6.5, SD 2.1) are this package's fit of the uncoupled model to the joint
# loads of the treated Knapp et al. (2022) cohort (KnappJointModel(coupled = false),
# scripts/bacteria_ipm_knapp_fit.jl), not published estimates. Frogs in that trial were
# group-marked, so σy is identified only through how the cross-sections spread out over time,
# and r not at all. Letting the bacteria act on Bd growth (coupled = true) found no protective
# effect: with κ_growth fixed anywhere in [−2, 4] the log-likelihood is within 0.6 of the
# uncoupled fit, and the unconstrained maximum (κ_growth ≈ −107 above ln GE 10.4) is an
# implausible effect carried by a few swabs. ρ defaults to 0: in wild R. sierrae, Jani & Briggs
# (2014, PNAS 111:E5049) "found no correlation between Bd load and any OTUs of the genus
# Janthinobacterium" (rare in the field), and in their infection experiment an OTU of the genus
# "was positively associated with Bd infection" (relative abundance, not load), so nothing
# published supports Bd lowering J. lividum; ρ ≠ 0 is a scenario. Everything else is
# ILLUSTRATIVE: the protection threshold y50, its width, the effects κ, r and colonisation.
# Becker et al. (2009) support a threshold: "a threshold violacein concentration of about 18 μM on a salamander's
# skin prevents mortality" (violacein MIC against Bd 1.82 μM, Brucker et al. 2008), but they
# did not measure Bd load, and J. lividum load predicts violacein only weakly (r = 0.33), so
# the threshold cannot be placed on the load axis from published data.

"""
    BacteriaIPMParameters(; bd = WilberParameters(:wilber2016),
                          pool = Wilber2017Parameters(:wilber2017), decline = 0.1704,
                          σy = 0.29, ρ = 0.0, r = 0.0, y50 = log(100), width = 0.5,
                          κ_survival = 2.0, κ_growth = 1.0, κ_clearance = 0.0,
                          κ_infection = 0.0, colonisation = 0.0, exchange = 0.0,
                          colonise_mean = 2.0, colonise_sd = 1.0, bath_mean = 6.5,
                          bath_sd = 2.1)

Parameters of the two-trait Bd–bacteria IPM ([`BacteriaIPMProblem`](@ref)), a **new model**
(not from a paper). `x` is ln Bd load (ZE per swab) and `y` ln load of the applied bacterium
*Janthinobacterium lividum* (genomic equivalents per swab). Rates are per day; one step is
`bd.Δt` (3) days.

| field | meaning | source |
|:--|:--|:--|
| `bd` | Bd vital rates (survival, growth, clearance, initial load, `φ(T)`) | Wilber et al. 2016 |
| `pool` | zoospore pool: transmission, shedding `μA`, zoospore survival `ν(T)` | Wilber et al. 2017 |
| `decline` | decline of `y` per day, `δ` (half-life 4.07 days) | [`KnappJointModel`](@ref) fit (this package) |
| `σy` | SD of the change in `y` over a step | [`KnappJointModel`](@ref) fit (this package) |
| `ρ` | fall in `y` per step per unit of `x` (Bd disturbs the bacteria) | 0: no published support (Jani & Briggs 2014) |
| `r` | correlation of the changes in `x` and `y` | illustrative |
| `y50`, `width` | protection `h(y) = logistic((y − y50)/width)` | illustrative |
| `κ_survival` | added to logit survival at full protection (tolerance) | illustrative |
| `κ_growth` | subtracted from the mean of the next ln load (resistance) | illustrative |
| `κ_clearance` | added to logit clearance of Bd (resistance) | illustrative |
| `κ_infection` | proportional reduction of the infection probability, in `[0, 1]` | illustrative |
| `colonisation` | per-day rate at which frogs without bacteria acquire them from the pond | illustrative |
| `exchange` | per-day rate of acquisition per unit fraction of frogs carrying them | illustrative |
| `colonise_mean`, `colonise_sd` | `y` of newly colonised frogs | illustrative |
| `bath_mean`, `bath_sd` | `y` just after a probiotic bath | [`KnappJointModel`](@ref) fit (this package) |

Bacteria are lost when `y` falls below the lower edge of the `y` mesh (0 = 1 GE per swab, the
detection floor). With every `κ`, `ρ` and `r` zero the traits are independent and the Bd side
is the Wilber et al. (2016) IPM. Use `with(p; field = value)` to change fields.
"""
Base.@kwdef struct BacteriaIPMParameters
    bd::WilberParameters = WilberParameters(:wilber2016)
    pool::Wilber2017Parameters = Wilber2017Parameters(:wilber2017)
    decline::Float64 = 0.1704
    σy::Float64 = 0.29
    ρ::Float64 = 0.0
    r::Float64 = 0.0
    y50::Float64 = log(100.0)
    width::Float64 = 0.5
    κ_survival::Float64 = 2.0
    κ_growth::Float64 = 1.0
    κ_clearance::Float64 = 0.0
    κ_infection::Float64 = 0.0
    colonisation::Float64 = 0.0
    exchange::Float64 = 0.0
    colonise_mean::Float64 = 2.0
    colonise_sd::Float64 = 1.0
    bath_mean::Float64 = 6.5
    bath_sd::Float64 = 2.1
    function BacteriaIPMParameters(bd, pool, decline, σy, ρ, r, y50, width, κ_survival,
            κ_growth, κ_clearance, κ_infection, colonisation, exchange, colonise_mean,
            colonise_sd, bath_mean, bath_sd)
        (σy > 0 && width > 0 && colonise_sd > 0 && bath_sd > 0) ||
            throw(ArgumentError("σy, width, colonise_sd and bath_sd must be positive"))
        -1 < r < 1 || throw(ArgumentError("r must lie in (−1, 1)"))
        0 <= κ_infection <= 1 || throw(ArgumentError("κ_infection must lie in [0, 1]"))
        (decline >= 0 && colonisation >= 0 && exchange >= 0) ||
            throw(ArgumentError("decline, colonisation and exchange must be non-negative"))
        isapprox(bd.Δt, pool.Δt) || throw(ArgumentError("bd and pool must share Δt"))
        return new(bd, pool, decline, σy, ρ, r, y50, width, κ_survival, κ_growth,
            κ_clearance, κ_infection, colonisation, exchange, colonise_mean, colonise_sd,
            bath_mean, bath_sd)
    end
end

time_unit(::Type{<:BacteriaIPMParameters}) = :day

with(p::BacteriaIPMParameters; kw...) =
    BacteriaIPMParameters(; (f => get(kw, f, getfield(p, f)) for f in fieldnames(BacteriaIPMParameters))...)

_protection(p::BacteriaIPMParameters, y) = _logistic((y - p.y50) / p.width)

"""
    BacteriaIPMProblem(p::BacteriaIPMParameters; temperature,
                       xmesh = IPMMesh(:wilber2016; n = 50), ymesh = IPMMesh(0.0, 14.0, 28),
                       uninfected = 100, infected = 0, zoospores = 0, days = 120,
                       baths = (), transmission = :pool, φ = nothing, u0 = nothing)

A season of the two-trait Bd–bacteria IPM (a **new model**; see
[`BacteriaIPMParameters`](@ref)) for adult frogs without births, in steps of 3 days.

- `temperature`: °C, a number or a function of time in days.
- `xmesh`: ln Bd load (ZE per swab); `ymesh`: ln *J. lividum* load (GE per swab), whose lower
  edge is the detection floor below which a frog has lost its bacteria.
- `uninfected`, `infected`, `zoospores`: the starting frogs, all without bacteria (infected
  frogs at the initial-load distribution of Wilber et al. 2016), and zoospores in the pond.
  `u0` gives a full starting state instead (layout below).
- `baths`: probiotic baths, each `(day = d, coverage = f)`: at the start of the step beginning
  on day `d` (rounded down to a step), a fraction `f` of every class gets a fresh bacterial
  load `y ~ Normal(bath_mean, bath_sd)`.
- `transmission`: `:pool` (zoospore pool of Wilber et al. 2017), `:wilber2016` (the
  regression `φ(T)` of frogs housed singly) or `:constant` (infection probability `φ` per
  step, from a pond whose zoospores are not modelled). For the last two the pool is still
  tracked but does not infect.

The state is `u = [S0; S(y); I0(x); I(x, y); Z]`: uninfected without bacteria, uninfected per
`y` bin, infected without bacteria per `x` bin, infected per product bin (`x` fastest; see
[`IPMMesh2D`](@ref)), and zoospores. `solve(prob)` returns a [`BacteriaIPMSolution`](@ref).
"""
struct BacteriaIPMProblem{TT}
    p::BacteriaIPMParameters
    temperature::TT
    mesh::IPMMesh2D
    u0::Vector{Float64}
    steps::Int
    baths::Vector{NamedTuple{(:day, :coverage),Tuple{Float64,Float64}}}
    transmission::Symbol
    φ::Float64
end

function BacteriaIPMProblem(p::BacteriaIPMParameters; temperature,
        xmesh::IPMMesh = IPMMesh(:wilber2016; n = 50), ymesh::IPMMesh = IPMMesh(0.0, 14.0, 28),
        uninfected::Real = 100.0, infected::Real = 0.0, zoospores::Real = 0.0,
        days::Real = 120, baths = (), transmission::Symbol = :pool,
        φ::Union{Nothing,Real} = nothing, u0::Union{Nothing,AbstractVector} = nothing)
    transmission in (:pool, :wilber2016, :constant) ||
        throw(ArgumentError("transmission must be :pool, :wilber2016 or :constant"))
    (transmission === :constant) == (φ !== nothing) ||
        throw(ArgumentError("give φ exactly when transmission = :constant"))
    φ === nothing || 0 <= φ <= 1 || throw(ArgumentError("φ must lie in [0, 1]"))
    temperature isa Real || applicable(temperature, 0.0) ||
        throw(ArgumentError("temperature must be a number or a function of time (days)"))
    mesh = IPMMesh2D(xmesh, ymesh)
    nx, ny = size(mesh)
    N = 1 + ny + nx + nx * ny + 1
    Δt = p.bd.Δt
    steps = round(Int, days / Δt)
    isapprox(steps * Δt, days) || throw(ArgumentError("days must be a multiple of Δt = $Δt"))
    if u0 === nothing
        (uninfected >= 0 && infected >= 0 && zoospores >= 0) ||
            throw(ArgumentError("starting numbers must be non-negative"))
        u0 = zeros(N)
        u0[1] = uninfected
        T0 = Float64(_temperature_at(temperature, 0.0))
        rec = discretize(IPMKernel(p.bd; temperature = T0), xmesh).recruitment
        u0[1+ny+1:1+ny+nx] .= infected .* rec ./ sum(rec)
        u0[N] = zoospores
    end
    length(u0) == N || throw(DimensionMismatch("u0 must have length $N"))
    all(>=(0), u0) || throw(ArgumentError("u0 must be non-negative"))
    bs = NamedTuple{(:day, :coverage),Tuple{Float64,Float64}}[
        (day = Float64(b.day), coverage = Float64(b.coverage)) for b in baths]
    all(b -> 0 <= b.coverage <= 1, bs) || throw(ArgumentError("bath coverage must lie in [0, 1]"))
    return BacteriaIPMProblem(p, temperature, mesh, collect(Float64, u0), steps, bs,
        transmission, φ === nothing ? NaN : Float64(φ))
end

"""
    BacteriaIPMSolution

The result of `solve(::BacteriaIPMProblem)`: `t` (days), `u` (the state at each time, one
column per time; see [`BacteriaIPMProblem`](@ref) for the layout), `φ` (the infection
probability of a frog without bacteria in each step), `temperature` (°C in each step) and
`mesh` (an [`IPMMesh2D`](@ref)). [`bacteria_ipm_summary`](@ref) gives time series and
[`joint_load_distribution`](@ref) the joint distribution of the two loads.
"""
struct BacteriaIPMSolution
    t::Vector{Float64}
    u::Matrix{Float64}
    φ::Vector{Float64}
    temperature::Vector{Float64}
    mesh::IPMMesh2D
end

# Index ranges of the classes in the state vector.
function _bipm_ranges(mesh::IPMMesh2D)
    nx, ny = size(mesh)
    S = 2:1+ny
    I0 = 2+ny:1+ny+nx
    I = 2+ny+nx:1+ny+nx+nx*ny
    return (S0 = 1, S = S, I0 = I0, I = I, Z = 1 + ny + nx + nx * ny + 1)
end

# Everything of one step that depends on the temperature only.
struct _BacteriaPieces
    d0::DiscretizedIPM            # Wilber kernel: infected frogs without bacteria
    PI::Matrix{Float64}           # I(x,y) → I(x′,y′)
    I_to_I0::Matrix{Float64}      # I(x,y) → I0(x′): keep Bd, lose bacteria
    I_to_S::Matrix{Float64}       # I(x,y) → S(y′): clear Bd, keep bacteria
    I_to_S0::Vector{Float64}      # I(x,y) → S0: clear Bd, lose bacteria
    S_keep::Matrix{Float64}       # S(y) → S(y′), not infected (before s0 and 1 − φ)
    S_lose::Vector{Float64}       # S(y) → S0, not infected
    S_inf_I::Matrix{Float64}      # S(y) → I(x′,y′), newly infected (before s0 and φ)
    S_inf_I0::Matrix{Float64}     # S(y) → I0(x′), newly infected, bacteria lost
    infection_factor::Vector{Float64}  # 1 − κi h(y) per y bin
    colonise::Vector{Float64}     # y distribution of newly colonised frogs over y bins
    bath::Vector{Float64}         # y distribution just after a bath over y bins
    ν::Float64                    # zoospore survival per step
end

function _bacteria_pieces(p::BacteriaIPMParameters, mesh::IPMMesh2D, T::Float64)
    xm, ym = mesh.x, mesh.y
    nx, ny = size(mesh)
    xs, ys = meshpoints(mesh)
    hx = meshwidth(xm)
    bd = p.bd
    Δt = bd.Δt
    d0 = discretize(IPMKernel(bd; temperature = T), xm)
    nI = nx * ny
    PI = zeros(nI, nI)
    I_to_I0 = zeros(nx, nI)
    I_to_S = zeros(ny, nI)
    I_to_S0 = zeros(nI)
    py = zeros(ny)
    a = -p.decline * Δt
    σc = p.σy * sqrt(1 - p.r^2)
    for k in 1:ny, j in 1:nx
        col = j + (k - 1) * nx
        x, y = xs[j], ys[k]
        hp = _protection(p, y)
        s = _logistic(bd.b00 + bd.b10 * x + p.κ_survival * hp)
        l = _logistic(bd.b02 + bd.b12 * x + bd.b22 * T + p.κ_clearance * hp)
        μx = bd.b01 + bd.b11 * x + bd.b21 * T - p.κ_growth * hp
        σx = sqrt(bd.v01 * exp(2 * bd.c01 * x))
        μy = y + a - p.ρ * x
        stay = s * (1 - l)
        below = 0.0
        if p.r == 0
            below, above = _normal_bin_masses!(py, μy, σc, ym)
            py[end] += above
        end
        @inbounds for i in 1:nx
            pxi = hx * pdf(Normal(μx, σx), xs[i])        # midpoint rule, as Wilber et al.
            pxi == 0 && continue
            if p.r != 0
                below, above = _normal_bin_masses!(py, μy + p.r * p.σy / σx * (xs[i] - μx), σc, ym)
                py[end] += above
            end
            w = stay * pxi
            for kk in 1:ny
                PI[i+(kk-1)*nx, col] = w * py[kk]
            end
            I_to_I0[i, col] = w * below
        end
        below, above = _normal_bin_masses!(py, μy, p.σy, ym)    # frogs that clear Bd
        py[end] += above
        @views I_to_S[:, col] .= (s * l) .* py
        I_to_S0[col] = s * l * below
    end
    S_keep = zeros(ny, ny)
    S_lose = zeros(ny)
    S_inf_I = zeros(nI, ny)
    S_inf_I0 = zeros(nx, ny)
    rec = d0.recruitment
    for k in 1:ny
        below, above = _normal_bin_masses!(py, ys[k] + a, p.σy, ym)
        py[end] += above
        @views S_keep[:, k] .= py
        S_lose[k] = below
        for kk in 1:ny, i in 1:nx
            S_inf_I[i+(kk-1)*nx, k] = rec[i] * py[kk]
        end
        @views S_inf_I0[:, k] .= rec .* below
    end
    inf = [1 - p.κ_infection * _protection(p, y) for y in ys]
    col = zeros(ny)
    _, above = _normal_bin_masses!(col, p.colonise_mean, p.colonise_sd, ym)
    col[end] += above
    bth = zeros(ny)
    _, above = _normal_bin_masses!(bth, p.bath_mean, p.bath_sd, ym)
    bth[end] += above
    return _BacteriaPieces(d0, PI, I_to_I0, I_to_S, I_to_S0, S_keep, S_lose, S_inf_I,
        S_inf_I0, inf, col, bth, zoospore_survival(p.pool, T))
end

# Apply a bath of coverage f to the state u in place: a fraction f of every frog class gets
# y ~ bath distribution (the share of it below the floor stays or becomes bacteria-free).
function _bath!(u, f, pc::_BacteriaPieces, mesh::IPMMesh2D)
    f == 0 && return u
    R = _bipm_ranges(mesh)
    nx, ny = size(mesh)
    q = pc.bath
    qlost = 1 - sum(q)
    Utot = u[R.S0] + sum(@view u[R.S])
    u[R.S0] = (1 - f) * u[R.S0] + f * qlost * Utot
    @views u[R.S] .= (1 - f) .* u[R.S] .+ (f * Utot) .* q
    I0 = u[R.I0]
    Im = reshape(view(u, R.I), nx, ny)
    Ix = I0 .+ vec(sum(Im; dims = 2))                    # infected frogs per x bin
    @views u[R.I0] .= (1 - f) .* I0 .+ (f * qlost) .* Ix
    Im .= (1 - f) .* Im .+ f .* (Ix * q')
    return u
end

# Infection probability of a frog without bacteria in this step.
function _bacteria_phi(prob::BacteriaIPMProblem, u, T)
    p = prob.p
    prob.transmission === :wilber2016 && return infection_probability(p.bd, T)
    prob.transmission === :constant && return prob.φ
    R = _bipm_ranges(prob.mesh)
    Itot = sum(@view u[R.I0]) + sum(@view u[R.I])
    A = u[R.S0] + sum(@view u[R.S]) + Itot
    return infection_probability(p.pool; Z = u[R.Z], I = Itot, A)
end

# One step from u to out (no bath; apply _bath! to u first).
function _bacteria_step!(out, u, prob::BacteriaIPMProblem, pc::_BacteriaPieces, φ)
    p = prob.p
    mesh = prob.mesh
    R = _bipm_ranges(mesh)
    nx, ny = size(mesh)
    xs, _ = meshpoints(mesh)
    s0 = p.bd.s0
    S0 = u[R.S0]
    S = @view u[R.S]
    I0 = @view u[R.I0]
    I = @view u[R.I]
    φy = φ .* pc.infection_factor
    Sescape = s0 .* (1 .- φy) .* S
    Sinf = s0 .* φy .* S
    newS0 = s0 * (1 - φ) * S0 + dot(pc.d0.loss, I0) + dot(pc.I_to_S0, I) + dot(pc.S_lose, Sescape)
    newS = pc.S_keep * Sescape .+ pc.I_to_S * I
    newI0 = pc.d0.P * I0 .+ (s0 * φ * S0) .* pc.d0.recruitment .+ pc.S_inf_I0 * Sinf .+
            pc.I_to_I0 * I
    newI = pc.PI * I .+ pc.S_inf_I * Sinf
    # colonisation of frogs without bacteria, from the pond and from frogs that carry them
    alive = S0 + sum(S) + sum(I0) + sum(I)
    carrying = alive > 0 ? (sum(S) + sum(I)) / alive : 0.0
    c = 1 - exp(-(p.colonisation + p.exchange * carrying) * p.bd.Δt)
    if c > 0
        q = pc.colonise
        moved = c * sum(q)
        newS .+= (c * newS0) .* q
        newS0 *= 1 - moved
        reshape(newI, nx, ny) .+= c .* (newI0 * q')
        newI0 .*= 1 - moved
    end
    shed = p.pool.μA * (dot(exp.(xs), I0) + dot(repeat(exp.(xs), ny), I))
    zleft = pc.ν * u[R.Z]
    out[R.S0] = newS0
    @views out[R.S] .= newS
    @views out[R.I0] .= newI0
    @views out[R.I] .= newI
    out[R.Z] = zleft + shed
    return out
end

function CommonSolve.solve(prob::BacteriaIPMProblem)
    N = length(prob.u0)
    Δt = prob.p.bd.Δt
    U = Matrix{Float64}(undef, N, prob.steps + 1)
    U[:, 1] .= prob.u0
    φs = Vector{Float64}(undef, prob.steps)
    Ts = Vector{Float64}(undef, prob.steps)
    cache = Dict{Float64,_BacteriaPieces}()
    u = copy(prob.u0)
    for k in 1:prob.steps
        t0 = (k - 1) * Δt
        T = Float64(_temperature_at(prob.temperature, t0))
        pc = get!(() -> _bacteria_pieces(prob.p, prob.mesh, T), cache, T)
        copyto!(u, view(U, :, k))
        for b in prob.baths
            t0 <= b.day < t0 + Δt && _bath!(u, b.coverage, pc, prob.mesh)
        end
        φs[k] = _bacteria_phi(prob, u, T)
        Ts[k] = T
        @views _bacteria_step!(U[:, k+1], u, prob, pc, φs[k])
    end
    return BacteriaIPMSolution(collect(0:prob.steps) .* Δt, U, φs, Ts, prob.mesh)
end

"""
    bacteria_ipm_summary(sol::BacteriaIPMSolution) -> NamedTuple

Time series (one entry per time in `sol.t`, days) of a [`BacteriaIPMProblem`](@ref) run:
`frogs` alive, `uninfected`, `infected`, `prevalence`, `carrying` (fraction of live frogs
carrying the bacteria), `mean_bd` (mean ln Bd load of infected frogs), `mean_bacteria` (mean ln
load of the bacteria among carriers) and `zoospores`. Means are `NaN` when there is no one to
average over.
"""
function bacteria_ipm_summary(sol::BacteriaIPMSolution)
    R = _bipm_ranges(sol.mesh)
    nx, ny = size(sol.mesh)
    xs, ys = meshpoints(sol.mesh)
    K = length(sol.t)
    frogs, unf, inf, carry, mbd, mbac = (zeros(K) for _ in 1:6)
    for c in 1:K
        u = view(sol.u, :, c)
        S0, S, I0 = u[R.S0], view(u, R.S), view(u, R.I0)
        Im = reshape(view(u, R.I), nx, ny)
        Ix = I0 .+ vec(sum(Im; dims = 2))
        Iy = vec(sum(Im; dims = 1))
        unf[c] = S0 + sum(S)
        inf[c] = sum(Ix)
        frogs[c] = unf[c] + inf[c]
        carriers = sum(S) + sum(Iy)
        carry[c] = frogs[c] > 0 ? carriers / frogs[c] : NaN
        mbd[c] = inf[c] > 0 ? dot(xs, Ix) / inf[c] : NaN
        mbac[c] = carriers > 0 ? dot(ys, S .+ Iy) / carriers : NaN
    end
    return (t = sol.t, frogs = frogs, uninfected = unf, infected = inf,
        prevalence = [f > 0 ? i / f : NaN for (i, f) in zip(inf, frogs)], carrying = carry,
        mean_bd = mbd, mean_bacteria = mbac, zoospores = sol.u[R.Z, :])
end

"""
    joint_load_distribution(sol::BacteriaIPMSolution, k::Integer) -> NamedTuple

The frogs at time index `k` of a [`BacteriaIPMProblem`](@ref) run (time `sol.t[k]`) by class:
`x` and `y` (bin midpoints: ln Bd load, ln bacterial load), `infected` (`nx × ny` matrix of
infected frogs carrying bacteria), `infected_without` (infected frogs without bacteria, per `x`
bin), `uninfected` (uninfected frogs carrying bacteria, per `y` bin) and `uninfected_without`.
"""
function joint_load_distribution(sol::BacteriaIPMSolution, k::Integer)
    R = _bipm_ranges(sol.mesh)
    nx, ny = size(sol.mesh)
    xs, ys = meshpoints(sol.mesh)
    u = sol.u[:, k]
    return (x = xs, y = ys, infected = reshape(u[R.I], nx, ny), infected_without = u[R.I0],
        uninfected = u[R.S], uninfected_without = u[R.S0])
end

"""
    basic_reproduction_number(prob::BacteriaIPMProblem; temperature = nothing) -> Float64

The basic reproduction number of Bd in the uninfected frogs of `prob.u0` (with their bacteria
as given there), for `transmission = :pool`: the expected number of frogs infected by one
newly infected frog, through the zoospore pool and (for `:zpool_dd` and `:zpool_fd`)
directly, over its infection. It is the dominant eigenvalue of the next-generation matrix
`F (I − T)⁻¹` of the step linearised at no Bd (`ln(Z + 1) ≈ Z`), where `T` holds the changes
of infected frogs (both loads; colonisation at the rate of `u0`) and of the pool, and `F` the
new infections. Every new infection enters through the same infection probability, so `F` has
rank one and `R0 = cᵀ (I − T)⁻¹ n`, with `n` the distribution of newly infected frogs over
both loads. The bacteria of uninfected frogs are held as in `u0` (no decline), so this is
`R0` for Bd arriving at that moment. `temperature` defaults to the temperature at time 0.
"""
function basic_reproduction_number(prob::BacteriaIPMProblem; temperature = nothing)
    prob.transmission === :pool ||
        throw(ArgumentError("R0 is defined here for transmission = :pool"))
    p = prob.p
    pool = p.pool
    pool.transmission in (:zpool, :zpool_dd, :zpool_fd) ||
        throw(ArgumentError("R0 needs a zoospore-pool transmission (:zpool, :zpool_dd or :zpool_fd)"))
    T = temperature === nothing ? Float64(_temperature_at(prob.temperature, 0.0)) :
        Float64(temperature)
    pc = _bacteria_pieces(p, prob.mesh, T)
    mesh = prob.mesh
    R = _bipm_ranges(mesh)
    nx, ny = size(mesh)
    xs, _ = meshpoints(mesh)
    u = prob.u0
    S0, S = u[R.S0], u[R.S]
    A = S0 + sum(S)
    A > 0 || return 0.0
    c = 1 - exp(-(p.colonisation + p.exchange * sum(S) / A) * p.bd.Δt)
    q = pc.colonise
    moved = c * sum(q)
    nI = nx * ny
    M = nx + nI + 1                                       # [I0; I; Z]
    colonise!(v) = begin
        I0 = v[1:nx]
        reshape(view(v, nx+1:nx+nI), nx, ny) .+= c .* (I0 * q')
        v[1:nx] .= I0 .* (1 - moved)
        v
    end
    Tm = zeros(M, M)
    @views Tm[1:nx, 1:nx] .= pc.d0.P
    @views Tm[nx+1:nx+nI, nx+1:nx+nI] .= pc.PI
    @views Tm[1:nx, nx+1:nx+nI] .= pc.I_to_I0
    for j in 1:nx+nI
        colonise!(view(Tm, 1:nx+nI, j))
    end
    @views Tm[M, 1:nx] .= pool.μA .* exp.(xs)
    @views Tm[M, nx+1:nx+nI] .= pool.μA .* repeat(exp.(xs), ny)
    Tm[M, M] = pc.ν
    s0 = p.bd.s0
    Sφ = s0 .* pc.infection_factor .* S
    newv = zeros(M)
    @views newv[1:nx] .= (s0 * S0) .* pc.d0.recruitment .+ pc.S_inf_I0 * Sφ
    @views newv[nx+1:nx+nI] .= pc.S_inf_I * Sφ
    colonise!(view(newv, 1:nx+nI))
    # dφ/d(state) at no Bd: β0 Δt per zoospore; β1 Δt (dd) or β1 Δt / A (fd) per infected frog
    Δt = pool.Δt
    cvec = zeros(M)
    cvec[M] = pool.β0 * Δt
    direct = pool.transmission === :zpool_dd ? pool.β1 * Δt :
             pool.transmission === :zpool_fd ? pool.β1 * Δt / A : 0.0
    cvec[1:nx+nI] .= direct
    return dot(cvec, (LinearAlgebra.I - Tm) \ newv)
end
