# Two-trait integral projection models: a product mesh of two IPMMesh axes, for models in which
# each host carries two continuous traits (src/novel/bacteria_ipm.jl: Bd load and the load of
# an antifungal bacterium). The one-trait machinery (ipm.jl) is unchanged.
#
# A density on the product mesh is stored as a vector in column-major order with the first
# trait fastest: bin (i, k), the i-th bin of the first trait and the k-th of the second, has
# linear index i + (k − 1) nx. Bin masses of the second trait use the exact normal cdf, so that
# the mass falling below the lower edge (for a load, below detection) is known exactly.

"""
    IPMMesh2D(x::IPMMesh, y::IPMMesh)

The product of two one-trait meshes, for an IPM in which each host carries two continuous
traits `x` and `y`. It has `x.n × y.n` bins. A density on it is stored as a vector of length
`x.n * y.n` with `x` fastest: bin `(i, k)` has linear index `i + (k − 1) x.n`, so
`reshape(v, x.n, y.n)` gives the matrix with rows for `x` and columns for `y`.
`meshpoints(m)` returns the two vectors of bin midpoints.
"""
struct IPMMesh2D
    x::IPMMesh
    y::IPMMesh
end

Base.length(m::IPMMesh2D) = m.x.n * m.y.n
Base.size(m::IPMMesh2D) = (m.x.n, m.y.n)
meshpoints(m::IPMMesh2D) = (meshpoints(m.x), meshpoints(m.y))

# Probability that Normal(μ, σ) falls in each bin of `mesh` (exact bin integrals), written into
# `out`; returns the masses below the lower edge and above the upper edge.
function _normal_bin_masses!(out::AbstractVector, μ::Real, σ::Real, mesh::IPMMesh)
    d = Normal(μ, σ)
    h = meshwidth(mesh)
    below = cdf(d, mesh.lower)
    c = below
    @inbounds for i in 1:mesh.n
        c2 = i == mesh.n ? cdf(d, mesh.upper) : cdf(d, mesh.lower + i * h)
        out[i] = max(c2 - c, 0.0)
        c = c2
    end
    return below, ccdf(d, mesh.upper)
end

# Marginal totals of a density `v` on a product mesh: per bin of x (summed over y) and per bin
# of y (summed over x).
function _mesh_marginals(m::IPMMesh2D, v::AbstractVector)
    M = reshape(v, m.x.n, m.y.n)
    return vec(sum(M; dims = 2)), vec(sum(M; dims = 1))
end
