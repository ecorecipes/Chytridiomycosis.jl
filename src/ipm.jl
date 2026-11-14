# Generic integral projection model (IPM) machinery for host–parasite systems whose continuous
# state is (log) pathogen load, after Wilber et al. (2016, Methods Ecol Evol 7:1182) and the
# IPM guides they follow (Easterling, Ellner & Dixon 2000; Rees, Childs & Ellner 2014;
# Williams, Miller & Ellner 2012 on eviction). Model-specific code (Wilber 2016:
# src/models/wilber2016.jl) builds an IPMKernel from its vital-rate functions and uses the
# analyses here.
#
# State: one discrete uninfected class S followed by the infected density I(x) on a mesh,
# u = [S; n_1, …, n_m], where n_i ≈ ∫ over bin i of I(x) dx (numbers of hosts, not densities).
#
#   S(t+1)  = s0 (1 − φ) S(t) + ∫ s(x) l(x) I(x,t) dx                      (eqn 1 / 7)
#   I(x′,t+1) = ∫ s(x) (1 − l(x)) G(x′|x) I(x,t) dx + s0 φ S(t) G0(x′)       (eqn 2 / 8)
#
# Reproduction is not included (Wilber 2016 excluded it). The infection probability φ is
# passed separately because it may depend on the state (density-dependent transmission).
#
# Why not IntegralProjectionModels.jl (ecorecipes)? It sources StructuredPopulationCore by
# relative path (CI cannot install it), takes Plots as a hard dependency and pins SciMLBase 2
# (this package uses SciMLBase 3). The few hundred lines here need nothing new.

"""
    IPMMesh(lower, upper, n)

A midpoint-rule mesh of `n` equal bins on `[lower, upper]` (Easterling et al. 2000; Rees et
al. 2014). [`meshpoints`](@ref) gives the bin midpoints, [`meshwidth`](@ref) the bin width
`h = (upper − lower)/n` and [`meshedges`](@ref) the `n + 1` bin edges.

    IPMMesh(:wilber2016; n = 100)

The mesh of Wilber et al. (2016), `[−5, 18]` on the natural-log load scale, with the bounds
and the default `n` read from the compendium. The authors' bootstrap, elasticity and `R0`
code uses `n = 50`.
"""
struct IPMMesh
    lower::Float64
    upper::Float64
    n::Int
    function IPMMesh(lower::Real, upper::Real, n::Integer)
        upper > lower || throw(ArgumentError("need upper > lower"))
        n >= 1 || throw(ArgumentError("need at least one bin"))
        return new(Float64(lower), Float64(upper), Int(n))
    end
end

Base.length(m::IPMMesh) = m.n

"""
    meshwidth(m::IPMMesh) -> Float64

Bin width `h = (upper − lower)/n`.
"""
meshwidth(m::IPMMesh) = (m.upper - m.lower) / m.n

"""
    meshedges(m::IPMMesh) -> Vector{Float64}

The `n + 1` bin edges `lower, lower + h, …, upper`.
"""
meshedges(m::IPMMesh) = [m.lower + i * meshwidth(m) for i in 0:m.n]

"""
    meshpoints(m::IPMMesh) -> Vector{Float64}

The `n` bin midpoints, where the midpoint rule evaluates the kernel.
"""
meshpoints(m::IPMMesh) = [m.lower + (i - 0.5) * meshwidth(m) for i in 1:m.n]

"""
    IPMKernel(; survival, growth, loss = x -> 0.0, initial_load, uninfected_survival = 1.0)

The vital-rate functions of a host–parasite IPM with an uninfected class (Wilber et al.
2016, eqns 1–2 and Fig. 1), on whatever load scale the mesh uses:

- `survival(x)`: probability that an infected host with load `x` survives a time step, `s(x)`;
- `growth(x)`: the distribution (a `Distributions.UnivariateDistribution`) of the load at
  `t + 1` of a host that survives with load `x` and keeps its infection, `G(x′, x)`;
- `loss(x)`: probability that a surviving host with load `x` loses its infection and
  returns to the uninfected class, `l(x)`;
- `initial_load`: the distribution of load on first infection, `G0(x′)`;
- `uninfected_survival`: survival of uninfected hosts, `s0`.

Infected hosts survive, then lose their infection or not, then grow: the infected kernel is
`s(x)(1 − l(x))G(x′, x)` and the return to the uninfected class is `s(x) l(x)`. The infection
probability `φ` is not part of the kernel; it is given to [`projection_matrix`](@ref) or
computed from the state. Use [`discretize`](@ref) to turn a kernel into matrices.

    IPMKernel(p::WilberParameters; temperature)

The kernel of Wilber et al. (2016) at `temperature` (°C): survival `s(x)`, normal growth
`G(x′, x, T)`, loss of infection `l(x, T)`, normal initial load `G0(x′, T)` and `s0` (see
[`WilberParameters`](@ref)).
"""
struct IPMKernel{S,G,L,G0<:UnivariateDistribution}
    survival::S
    growth::G
    loss::L
    initial_load::G0
    uninfected_survival::Float64
    function IPMKernel(survival::S, growth::G, loss::L, initial_load::G0,
            uninfected_survival::Real) where {S,G,L,G0<:UnivariateDistribution}
        0 <= uninfected_survival <= 1 ||
            throw(ArgumentError("uninfected_survival must lie in [0, 1]"))
        return new{S,G,L,G0}(survival, growth, loss, initial_load, Float64(uninfected_survival))
    end
end

IPMKernel(; survival, growth, loss = Returns(0.0), initial_load::UnivariateDistribution,
        uninfected_survival::Real = 1.0) =
    IPMKernel(survival, growth, loss, initial_load, uninfected_survival)

"""
    DiscretizedIPM

A host–parasite IPM discretised on a mesh by [`discretize`](@ref). Fields:

- `mesh`: the [`IPMMesh`](@ref);
- `P`: the `n × n` infected-to-infected matrix, `P[i, j]` the expected number of hosts in
  bin `i` at `t + 1` per host in bin `j` at `t` (columns are "from");
- `loss`: the vector `s(x_j) l(x_j)`, infected hosts in bin `j` returning to the uninfected
  class;
- `recruitment`: the distribution over bins of the load of a newly infected host
  (`h G0(x_i)` with the midpoint rule);
- `uninfected_survival`: `s0`.

[`projection_matrix`](@ref) assembles the full `(n + 1) × (n + 1)` matrix for a given
infection probability, with the uninfected class first.
"""
struct DiscretizedIPM
    mesh::IPMMesh
    P::Matrix{Float64}
    loss::Vector{Float64}
    recruitment::Vector{Float64}
    uninfected_survival::Float64
end

const IPM_INTEGRATION = (:midpoint, :cdf)
const IPM_EVICTION = (:none, :boundary)

# Probability mass of `d` in each bin: the midpoint rule (h × density at the midpoint, as
# Wilber et al. 2016 and Rees et al. 2014 do) or the exact bin integral from the cdf.
function _bin_masses!(out, d::UnivariateDistribution, mesh::IPMMesh, ys, edges,
        integration::Symbol, eviction::Symbol)
    h = meshwidth(mesh)
    if integration === :midpoint
        @inbounds for i in eachindex(ys)
            out[i] = h * pdf(d, ys[i])
        end
    else
        c = cdf(d, edges[1])
        @inbounds for i in eachindex(ys)
            c2 = cdf(d, edges[i+1])
            out[i] = c2 - c
            c = c2
        end
    end
    if eviction === :boundary
        # Williams et al. (2012): return evicted mass to the boundary classes
        out[1] += cdf(d, mesh.lower)
        out[end] += ccdf(d, mesh.upper)
    end
    return out
end

"""
    discretize(k::IPMKernel, mesh::IPMMesh; integration = :midpoint, eviction = :none)
        -> DiscretizedIPM

Discretise the kernel on `mesh`.

- `integration = :midpoint` (default) evaluates the growth and initial-load densities at bin
  midpoints and multiplies by the bin width, as Wilber et al. (2016) do. `:cdf` integrates
  each density exactly over each bin.
- `eviction = :none` (default) lets hosts whose next load falls outside the mesh leave the
  model, as in Wilber et al. (2016), who chose wide bounds instead. `:boundary` returns that
  mass to the first or last bin (a floor and ceiling; Williams, Miller & Ellner 2012).

[`eviction`](@ref) measures how much is lost; [`mesh_convergence`](@ref) checks the mesh size.
"""
function discretize(k::IPMKernel, mesh::IPMMesh; integration::Symbol = :midpoint,
        eviction::Symbol = :none)
    integration in IPM_INTEGRATION ||
        throw(ArgumentError("integration must be one of $IPM_INTEGRATION"))
    eviction in IPM_EVICTION || throw(ArgumentError("eviction must be one of $IPM_EVICTION"))
    ys, edges, n = meshpoints(mesh), meshedges(mesh), mesh.n
    P = Matrix{Float64}(undef, n, n)
    lossv = Vector{Float64}(undef, n)
    col = Vector{Float64}(undef, n)
    for j in 1:n
        x = ys[j]
        s = k.survival(x)
        l = k.loss(x)
        (0 <= s <= 1 && 0 <= l <= 1) ||
            throw(DomainError((s, l), "survival and loss must be probabilities at x = $x"))
        _bin_masses!(col, k.growth(x), mesh, ys, edges, integration, eviction)
        @views P[:, j] .= col .* (s * (1 - l))
        lossv[j] = s * l
    end
    rec = _bin_masses!(Vector{Float64}(undef, n), k.initial_load, mesh, ys, edges,
        integration, eviction)
    return DiscretizedIPM(mesh, P, lossv, rec, k.uninfected_survival)
end

"""
    projection_matrix(d::DiscretizedIPM, φ) -> Matrix{Float64}

The `(n + 1) × (n + 1)` projection matrix for infection probability `φ` per time step, with
the uninfected class first (Wilber et al. 2016, "Analysing the IPM"; Merow et al. 2014):

```
A = [ s0(1 − φ)        s(x)ᵀ l(x)ᵀ
      s0 φ G0           P            ]
```
"""
function projection_matrix(d::DiscretizedIPM, φ::Real)
    0 <= φ <= 1 || throw(DomainError(φ, "φ must be a probability"))
    n = d.mesh.n
    A = Matrix{Float64}(undef, n + 1, n + 1)
    s0 = d.uninfected_survival
    A[1, 1] = s0 * (1 - φ)
    @views A[1, 2:end] .= d.loss
    @views A[2:end, 1] .= s0 * φ .* d.recruitment
    @views A[2:end, 2:end] .= d.P
    return A
end

# One step of the IPM, u = [S; n], without forming the full matrix.
function _ipm_step!(out, d::DiscretizedIPM, u, φ)
    s0 = d.uninfected_survival
    S = u[1]
    Iv = @view u[2:end]
    out[1] = s0 * (1 - φ) * S + dot(d.loss, Iv)
    @views mul!(out[2:end], d.P, Iv)
    @views out[2:end] .+= (s0 * φ * S) .* d.recruitment
    return out
end

# Dominant eigenvalue and right (and, if asked, left) eigenvectors of a non-negative matrix.
function _dominant(A::AbstractMatrix; left::Bool = true)
    F = eigen(A)
    i = argmax(real.(F.values))
    λ = real(F.values[i])
    w = abs.(real.(F.vectors[:, i]))
    left || return λ, w ./ sum(w), w
    G = eigen(transpose(A))
    k = argmax(real.(G.values))
    v = abs.(real.(G.vectors[:, k]))
    return λ, w ./ sum(w), v
end

"""
    population_growth_rate(A::AbstractMatrix) -> Float64

The asymptotic population growth rate `λ`, the dominant eigenvalue of a projection matrix
(per time step). Model methods, such as
[`population_growth_rate(::WilberParameters)`](@ref), build the matrix first.
"""
population_growth_rate(A::AbstractMatrix) = maximum(real, eigvals(A))

"""
    stable_distribution(A::AbstractMatrix) -> Vector{Float64}

The stable state distribution: the dominant right eigenvector of `A`, scaled to sum to 1.
With a host–parasite matrix the first entry is the uninfected class.
"""
stable_distribution(A::AbstractMatrix) = _dominant(A; left = false)[2]

"""
    sensitivity_matrix(A::AbstractMatrix) -> Matrix{Float64}

Sensitivities `∂λ/∂a_ij = v_i w_j / ⟨v, w⟩` of the dominant eigenvalue to each entry of `A`
(Caswell 2001), from the left (`v`) and right (`w`) dominant eigenvectors.
"""
function sensitivity_matrix(A::AbstractMatrix)
    _, w, v = _dominant(A)
    return (v * w') ./ dot(v, w)
end

"""
    elasticity_matrix(A::AbstractMatrix) -> Matrix{Float64}

Elasticities `(a_ij/λ) ∂λ/∂a_ij`; they sum to 1.
"""
function elasticity_matrix(A::AbstractMatrix)
    λ = population_growth_rate(A)
    return sensitivity_matrix(A) .* A ./ λ
end

"""
    absorption_times(A::AbstractMatrix) -> (mean, variance)

Mean and variance of the number of time steps before absorption (death, and eviction from
the mesh), starting from each state of the transient matrix `A`: with the fundamental
matrix `N = (I − A)⁻¹`, `mean = 1ᵀN` and `variance = 1ᵀ(2N² − N) − mean²` (Caswell 2001;
the authors' `absorption_times`). Multiply by the time step for a duration.
"""
function absorption_times(A::AbstractMatrix)
    N = inv(LinearAlgebra.I - A)
    m = vec(sum(N; dims = 1))
    v = vec(sum(2 .* (N * N) .- N; dims = 1)) .- m .^ 2
    return (mean = m, variance = v)
end

"""
    eviction(k::IPMKernel, mesh::IPMMesh) -> (growth, initial)

How much of each density falls outside the mesh (Williams, Miller & Ellner 2012): `growth`
is the probability, for a host at each midpoint, that its next load lies outside
`[lower, upper]`; `initial` is the same for the initial-load distribution. The midpoint rule
with `eviction = :none` loses these hosts.
"""
function eviction(k::IPMKernel, mesh::IPMMesh)
    out(d) = cdf(d, mesh.lower) + ccdf(d, mesh.upper)
    return (growth = [out(k.growth(x)) for x in meshpoints(mesh)], initial = out(k.initial_load))
end

"""
    mesh_convergence(f, sizes; rtol = 1e-3) -> (sizes, values, relative_change, converged)

Evaluate a scalar summary `f(n)` of an IPM discretised with `n` bins for each `n` in
`sizes` (increasing), and report the relative change between successive sizes. `converged`
is true when the last change is below `rtol`. For example
`mesh_convergence(n -> population_growth_rate(p; temperature = 20, mesh = IPMMesh(-5, 18, n)), [25, 50, 100, 200])`.
"""
function mesh_convergence(f, sizes; rtol::Real = 1e-3)
    issorted(sizes) && length(sizes) >= 2 ||
        throw(ArgumentError("need at least two increasing mesh sizes"))
    vals = [float(f(n)) for n in sizes]
    rel = [abs(vals[i] - vals[i-1]) / max(abs(vals[i]), eps()) for i in 2:length(vals)]
    return (sizes = collect(sizes), values = vals, relative_change = rel,
            converged = rel[end] < rtol)
end

"""
    basic_reproduction_number(d::DiscretizedIPM, transmission) -> Float64

`R0` of a discretised host–parasite IPM by the next-generation method for matrix models
(Klepac & Caswell 2011, eq. 29; Wilber et al. 2016, Appendix S1): the spectral radius of
`F (I − P)⁻¹`, where `P` is the infected kernel and `F[i, j] = s0 G0_i t_j`, with `t_j`
(`transmission`, a scalar or a vector over bins) the expected number of hosts infected per
time step by one infected host in bin `j` in a fully susceptible population. Because `F` has
rank one, `R0 = tᵀ (I − P)⁻¹ s0 G0`.
"""
function basic_reproduction_number(d::DiscretizedIPM, transmission)
    t = transmission isa Real ? fill(Float64(transmission), d.mesh.n) : collect(transmission)
    length(t) == d.mesh.n || throw(DimensionMismatch("transmission must have one entry per bin"))
    return dot(t, (LinearAlgebra.I - d.P) \(d.uninfected_survival .* d.recruitment))
end

"""
    load_summary(mesh::IPMMesh, u::AbstractVector) -> NamedTuple

Summaries of a host–parasite state `u = [S; n_1, …, n_m]` (numbers of hosts): `total`,
`uninfected`, `infected`, `prevalence`, and the `mean`, `variance` and variance-to-mean
ratio `vm_ratio` of load among infected hosts, on the mesh's scale (Wilber et al. 2016,
Fig. 5 uses natural-log loads).
"""
function load_summary(mesh::IPMMesh, u::AbstractVector)
    length(u) == mesh.n + 1 || throw(DimensionMismatch("u must have length n + 1"))
    ys = meshpoints(mesh)
    Iv = @view u[2:end]
    Ni = sum(Iv)
    total = u[1] + Ni
    m = Ni > 0 ? dot(Iv, ys) / Ni : NaN
    v = Ni > 0 ? dot(Iv, ys .^ 2) / Ni - m^2 : NaN
    return (total = total, uninfected = u[1], infected = Ni,
            prevalence = total > 0 ? Ni / total : NaN, mean = m, variance = v,
            vm_ratio = v / m)
end
