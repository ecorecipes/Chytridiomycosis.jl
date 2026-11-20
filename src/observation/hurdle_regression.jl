# Load observation layer, part 5: hurdle regression. The two-part (hurdle) lognormal model
# of src/observation/load_likelihood.jl (LoadHurdleModel, lod = 0) with covariates: each
# swab's infection probability and its log load given infection have their own linear
# predictor, with fixed effects and Gaussian random intercepts (optionally correlated across
# groups, e.g. by a phylogeny). The parametrisation follows brms' `hurdle_lognormal` family:
#
#   y = 0:  hu                        hu = logistic(η_hu)  (probability of a ZERO)
#   y > 0:  (1 − hu) LogNormal(y | η_mu, σ)
#
#   η = a + (X − x̄) b + Σ_g sd_g (L_g z_g)[group_g],  z_g ~ N(0, I),  L_g L_gᵀ = C_g,
#
# with predictors centred at their means (brms' default, so the prior on `a` is a prior on the
# intercept at the mean covariates). Internally each predictor is also divided by its
# standard deviation s_j (for the optimiser); the unconstrained parameters `b` are then
# coefficients of standardised predictors, and their priors are rescaled to match
# (N(0, (prior_b s_j)²)), so the model is unchanged. `hurdle_coefficients` returns the
# coefficients on the original scale. The prevalence is 1 − hu; on the logit scale a
# prevalence coefficient is minus the hu coefficient. Random effects are non-centred. User:
# Grear et al. 2025 (src/models/grear2025.jl).

"""
    HurdleRegressionProblem(y, X, groups; names = ["x1", …], prior_b = 20, prior_a = 20,
                            prior_b_hu = 1, prior_a_hu = 1, prior_scale = 2.5,
                            prior_scale_hu = 2.5, sigma_prior = :sd)

Hurdle-lognormal regression of swab loads `y` (`0` = not detected) on the columns of `X` (no
intercept column; the same design for both parts), in the parametrisation of brms'
`hurdle_lognormal` (`hu`, the probability of a zero, is logit-linear; positive loads are
lognormal with log-mean linear and standard deviation `σ`). `groups` is a vector of random
intercepts, each a `NamedTuple` `(name, index, levels, correlation)`: `index[i]` is swab `i`'s
level, `levels` the number of levels and `correlation` `nothing` (independent levels) or a
`levels×levels` correlation matrix (e.g. a phylogenetic one). Each group has its own standard
deviation in each part.

Priors (all on the scale of the linear predictors): `b ~ N(0, prior_b²)` and
`a ~ N(0, prior_a²)` for the load part; `b_hu ~ N(0, prior_b_hu²)`, `a_hu ~ N(0, prior_a_hu²)`
for the hurdle part; half-Student-t(3, 0, `prior_scale`) on `σ` and the load-part standard
deviations and half-Student-t(3, 0, `prior_scale_hu`) on the hurdle-part ones. With
`sigma_prior = :variance` the half-Student-t prior on `σ` is put on `σ²` instead.

It is a LogDensityProblems target with gradient (`LogDensityOrder{1}`) on the unconstrained
parameters `prob.names` (standard deviations on the log scale, with their Jacobian). Use
[`posterior_mode`](@ref) for a MAP fit, [`hurdle_coefficients`](@ref) for the fixed effects
on the uncentred scale.
"""
struct HurdleRegressionProblem
    y::Vector{Float64}
    logy::Vector{Float64}
    positive::BitVector
    X::Matrix{Float64}
    xbar::Vector{Float64}
    xscale::Vector{Float64}
    covariates::Vector{String}
    groups::Vector{NamedTuple}
    chol::Vector{Union{Nothing,Matrix{Float64}}}
    priors::NamedTuple
    blocks::NamedTuple
    names::Vector{Symbol}
end

function HurdleRegressionProblem(y::AbstractVector{<:Real}, X::AbstractMatrix{<:Real},
        groups::AbstractVector = NamedTuple[];
        names::AbstractVector = ["x$j" for j in axes(X, 2)], prior_b::Real = 20,
        prior_a::Real = 20, prior_b_hu::Real = 1, prior_a_hu::Real = 1,
        prior_scale::Real = 2.5, prior_scale_hu::Real = 2.5, sigma_prior::Symbol = :sd)
    sigma_prior in (:sd, :variance) || throw(ArgumentError("sigma_prior must be :sd or :variance"))
    n, p = size(X)
    length(y) == n || throw(DimensionMismatch("y has $(length(y)) entries, X has $n rows"))
    all(v -> v >= 0, y) || throw(ArgumentError("loads must be non-negative"))
    length(names) == p || throw(DimensionMismatch("need $p covariate names"))
    pos = BitVector(y .> 0)
    any(pos) && !all(pos) || throw(ArgumentError("need both zero and positive loads"))
    xbar = vec(sum(X; dims = 1)) ./ n
    xs = [max(sqrt(sum(abs2, view(X, :, j) .- xbar[j]) / n), eps()) for j in 1:p]
    Xc = (Matrix{Float64}(X) .- xbar') ./ xs'
    chol = Union{Nothing,Matrix{Float64}}[]
    gs = NamedTuple[]
    for g in groups
        length(g.index) == n || throw(DimensionMismatch("group $(g.name): index length"))
        all(i -> 1 <= i <= g.levels, g.index) ||
            throw(ArgumentError("group $(g.name): index outside 1…$(g.levels)"))
        push!(gs, (name = Symbol(g.name), index = collect(Int, g.index), levels = Int(g.levels)))
        C = g.correlation
        push!(chol, C === nothing ? nothing :
              Matrix(LinearAlgebra.cholesky(Symmetric(Matrix{Float64}(C))).L))
    end
    nm = Symbol[]
    blk = Dict{Symbol,Any}()
    add!(key, labels) = (blk[key] = length(nm) .+ (1:length(labels)); append!(nm, labels))
    add!(:a, [:intercept])
    add!(:b, Symbol.(names))
    add!(:a_hu, [:hu_intercept])
    add!(:b_hu, Symbol.("hu_" .* string.(names)))
    add!(:log_sigma, [:log_sigma])
    gb = NamedTuple[]
    for g in gs
        s = string(g.name)
        r1 = length(nm) + 1
        append!(nm, [Symbol("log_sd_", s), Symbol("log_sd_hu_", s)])
        rz = length(nm) .+ (1:g.levels)
        append!(nm, [Symbol("z_", s, "_", k) for k in 1:g.levels])
        rzh = length(nm) .+ (1:g.levels)
        append!(nm, [Symbol("z_hu_", s, "_", k) for k in 1:g.levels])
        push!(gb, (log_sd = r1, log_sd_hu = r1 + 1, z = rz, z_hu = rzh))
    end
    blocks = (a = blk[:a][1], b = blk[:b], a_hu = blk[:a_hu][1], b_hu = blk[:b_hu],
              log_sigma = blk[:log_sigma][1], groups = gb)
    priors = (b = Float64(prior_b), a = Float64(prior_a), b_hu = Float64(prior_b_hu),
              a_hu = Float64(prior_a_hu), scale = Float64(prior_scale),
              scale_hu = Float64(prior_scale_hu), sigma_on = sigma_prior)
    ly = [v > 0 ? log(v) : 0.0 for v in y]
    return HurdleRegressionProblem(collect(Float64, y), ly, pos, Xc, xbar, xs, string.(names), gs,
                                   chol, priors, blocks, nm)
end

LogDensityProblems.dimension(p::HurdleRegressionProblem) = length(p.names)
LogDensityProblems.capabilities(::Type{<:HurdleRegressionProblem}) =
    LogDensityProblems.LogDensityOrder{1}()

# half-Student-t(3, 0, A) log density (up to a constant) of x = exp(s) plus the log Jacobian s,
# and its derivative in s
function _hst_logpdf(s, A)
    x = exp(s)
    return -2 * log1p((x / A)^2 / 3) + s, -4 * x^2 / (3A^2 + x^2) + 1
end

# the same prior on the variance x = exp(2s), with log Jacobian log 2 + 2s, and its derivative
function _hst_logpdf_variance(s, A)
    x = exp(2s)
    return -2 * log1p((x / A)^2 / 3) + log(2) + 2s, -8 * x^2 / (3A^2 + x^2) + 2
end

_hr_level_effects(prob, g, z) = prob.chol[g] === nothing ? z : prob.chol[g] * z

function _hr_eval(prob::HurdleRegressionProblem, θ::AbstractVector, want_grad::Bool)
    bl = prob.blocks
    pr = prob.priors
    n = length(prob.y)
    b = θ[bl.b]
    bh = θ[bl.b_hu]
    ηm = prob.X * b .+ θ[bl.a]
    ηh = prob.X * bh .+ θ[bl.a_hu]
    us = Vector{Vector{Float64}}()
    ush = Vector{Vector{Float64}}()
    for (g, gi) in enumerate(prob.groups)
        gb = bl.groups[g]
        u = _hr_level_effects(prob, g, θ[gb.z])
        uh = _hr_level_effects(prob, g, θ[gb.z_hu])
        sd, sdh = exp(θ[gb.log_sd]), exp(θ[gb.log_sd_hu])
        push!(us, u)
        push!(ush, uh)
        @inbounds for i in 1:n
            ηm[i] += sd * u[gi.index[i]]
            ηh[i] += sdh * uh[gi.index[i]]
        end
    end
    σ = exp(θ[bl.log_sigma])
    lp = 0.0
    rm = zeros(n)
    rh = zeros(n)
    dσ = 0.0
    @inbounds for i in 1:n
        h = ηh[i]
        if prob.positive[i]
            e = prob.logy[i] - ηm[i]
            lp += -_log1pexp(h) - log(σ) - 0.5 * log(2π) - prob.logy[i] - 0.5 * (e / σ)^2
            rh[i] = -1 / (1 + exp(-h))
            rm[i] = e / σ^2
            dσ += -1 + (e / σ)^2
        else
            lp += -_log1pexp(-h)
            rh[i] = 1 / (1 + exp(h))
        end
    end
    # priors
    lp += -0.5 * (θ[bl.a] / pr.a)^2 - 0.5 * sum(abs2, b ./ prob.xscale) / pr.b^2 -
          0.5 * (θ[bl.a_hu] / pr.a_hu)^2 - 0.5 * sum(abs2, bh ./ prob.xscale) / pr.b_hu^2
    ls, dls = pr.sigma_on === :sd ? _hst_logpdf(θ[bl.log_sigma], pr.scale) :
              _hst_logpdf_variance(θ[bl.log_sigma], pr.scale)
    lp += ls
    for (g, gb) in enumerate(bl.groups)
        l1, _ = _hst_logpdf(θ[gb.log_sd], pr.scale)
        l2, _ = _hst_logpdf(θ[gb.log_sd_hu], pr.scale_hu)
        lp += l1 + l2 - 0.5 * sum(abs2, θ[gb.z]) - 0.5 * sum(abs2, θ[gb.z_hu])
    end
    want_grad || return lp, nothing
    G = zeros(length(θ))
    G[bl.a] = sum(rm) - θ[bl.a] / pr.a^2
    G[bl.b] = prob.X' * rm .- b ./ (pr.b^2 .* prob.xscale .^ 2)
    G[bl.a_hu] = sum(rh) - θ[bl.a_hu] / pr.a_hu^2
    G[bl.b_hu] = prob.X' * rh .- bh ./ (pr.b_hu^2 .* prob.xscale .^ 2)
    G[bl.log_sigma] = dσ + dls
    for (g, gi) in enumerate(prob.groups)
        gb = bl.groups[g]
        R = zeros(gi.levels)
        Rh = zeros(gi.levels)
        @inbounds for i in 1:n
            R[gi.index[i]] += rm[i]
            Rh[gi.index[i]] += rh[i]
        end
        sd, sdh = exp(θ[gb.log_sd]), exp(θ[gb.log_sd_hu])
        _, d1 = _hst_logpdf(θ[gb.log_sd], pr.scale)
        _, d2 = _hst_logpdf(θ[gb.log_sd_hu], pr.scale_hu)
        G[gb.log_sd] = sd * dot(us[g], R) + d1
        G[gb.log_sd_hu] = sdh * dot(ush[g], Rh) + d2
        L = prob.chol[g]
        G[gb.z] = sd .* (L === nothing ? R : L' * R) .- θ[gb.z]
        G[gb.z_hu] = sdh .* (L === nothing ? Rh : L' * Rh) .- θ[gb.z_hu]
    end
    return lp, G
end

function LogDensityProblems.logdensity(prob::HurdleRegressionProblem, θ::AbstractVector)
    lp = _hr_eval(prob, θ, false)[1]
    return isnan(lp) ? -Inf : lp
end

function LogDensityProblems.logdensity_and_gradient(prob::HurdleRegressionProblem,
                                                    θ::AbstractVector)
    return _hr_eval(prob, θ, true)
end

"""
    hurdle_coefficients(prob::HurdleRegressionProblem, θ) -> NamedTuple

Fixed effects at `θ` on the uncentred covariate scale: `intercept` and `b` (log-load part),
`hu_intercept` and `b_hu` (logit of the probability of a zero), `prevalence_intercept` and
`b_prevalence` (logit of the probability of infection, the negatives of the hu ones), and
`sigma`; `covariates` names the entries of the `b` vectors.
"""
function hurdle_coefficients(prob::HurdleRegressionProblem, θ::AbstractVector)
    bl = prob.blocks
    b, bh = θ[bl.b] ./ prob.xscale, θ[bl.b_hu] ./ prob.xscale
    a = θ[bl.a] - dot(prob.xbar, b)
    ah = θ[bl.a_hu] - dot(prob.xbar, bh)
    return (covariates = prob.covariates, intercept = a, b = b, hu_intercept = ah, b_hu = bh,
            prevalence_intercept = -ah, b_prevalence = -bh, sigma = exp(θ[bl.log_sigma]),
            sd = Dict(g.name => exp(θ[bl.groups[k].log_sd]) for (k, g) in enumerate(prob.groups)),
            sd_hu = Dict(g.name => exp(θ[bl.groups[k].log_sd_hu])
                         for (k, g) in enumerate(prob.groups)))
end

"""
    posterior_mode(prob; x0 = nothing, iters = 5000, gtol = 1e-5, memory = 10) -> NamedTuple

Maximum a posteriori estimate of a LogDensityProblems target with gradient (such as a
[`HurdleRegressionProblem`](@ref)) by L-BFGS with backtracking. Returns `x`, `logdensity`,
`gradient` and `iterations`. For hierarchical models the mode is on the unconstrained
(log-sd, non-centred) scale, not a summary of the posterior.
"""
function posterior_mode(prob; x0::Union{Nothing,AbstractVector} = nothing,
                        iters::Integer = 5000, gtol::Real = 1e-5, memory::Integer = 10)
    d = LogDensityProblems.dimension(prob)
    x = x0 === nothing ? zeros(d) : collect(Float64, x0)
    fg(v) = (r = LogDensityProblems.logdensity_and_gradient(prob, v); (-r[1], -r[2]))
    f, g = fg(x)
    isfinite(f) || throw(ArgumentError("the log density is not finite at the start"))
    S = Vector{Vector{Float64}}()
    Y = Vector{Vector{Float64}}()
    it = 0
    for k in 1:iters
        it = k
        maximum(abs, g) < gtol && break
        q = copy(g)
        α = zeros(length(S))
        for j in length(S):-1:1
            α[j] = dot(S[j], q) / dot(Y[j], S[j])
            q .-= α[j] .* Y[j]
        end
        γ = isempty(S) ? 1 / max(1.0, maximum(abs, g)) : dot(S[end], Y[end]) / dot(Y[end], Y[end])
        q .*= γ
        for j in eachindex(S)
            β = dot(Y[j], q) / dot(Y[j], S[j])
            q .+= (α[j] - β) .* S[j]
        end
        dir = -q
        slope = dot(dir, g)
        if !(slope < 0)
            empty!(S); empty!(Y)
            dir = -g ./ max(1.0, maximum(abs, g))
            slope = dot(dir, g)
        end
        t = 1.0
        xn, fn, gn = x, f, g
        while t > 1e-16
            xn = x .+ t .* dir
            fn, gn = fg(xn)
            isfinite(fn) && fn <= f + 1e-4 * t * slope && break
            t /= 2
        end
        t <= 1e-16 && break
        s, yv = xn .- x, gn .- g
        if dot(s, yv) > 1e-12
            push!(S, s); push!(Y, yv)
            length(S) > memory && (popfirst!(S); popfirst!(Y))
        end
        x, f, g = xn, fn, gn
    end
    return (x = x, logdensity = -f, gradient = -g, iterations = it)
end

# --- phylogenetic correlation from a Newick tree --------------------------------------------------

"""
    phylogenetic_correlation(newick::AbstractString) -> (; species, correlation)

The Brownian-motion correlation matrix of the tips of a rooted tree with branch lengths in
Newick format (what `ape::vcv(tree, corr = TRUE)` returns): the covariance of tips `i`, `j` is
the length of the path they share from the root, scaled to unit variances. Tip labels have
underscores replaced by spaces.
"""
function phylogenetic_correlation(newick::AbstractString)
    s = strip(newick)
    s = endswith(s, ";") ? s[1:end-1] : s
    tips = String[]
    paths = Vector{Vector{Tuple{Int,Float64}}}()   # (node id, branch length) from root
    pos = Ref(1)
    nodeid = Ref(0)
    chars = collect(s)
    function readlabel()
        st = pos[]
        while pos[] <= length(chars) && !(chars[pos[]] in (',', ')', '(', ':'))
            pos[] += 1
        end
        return String(chars[st:pos[]-1])
    end
    function readlength()
        pos[] <= length(chars) && chars[pos[]] == ':' || return 0.0
        pos[] += 1
        st = pos[]
        while pos[] <= length(chars) && !(chars[pos[]] in (',', ')'))
            pos[] += 1
        end
        return parse(Float64, String(chars[st:pos[]-1]))
    end
    # parse a subtree; returns list of (tip index, path below this node's parent edge)
    function subtree(stack::Vector{Tuple{Int,Float64}})
        if chars[pos[]] == '('
            pos[] += 1
            nodeid[] += 1
            me = nodeid[]
            members = Int[]
            child = Vector{Tuple{Int,Float64}}[]
            while true
                append!(members, subtree(Tuple{Int,Float64}[]))
                chars[pos[]] == ',' ? (pos[] += 1) : break
            end
            chars[pos[]] == ')' || throw(ArgumentError("malformed Newick at $(pos[])"))
            pos[] += 1
            readlabel()
            len = readlength()
            for m in members
                pushfirst!(paths[m], (me, len))
            end
            return members
        else
            push!(tips, replace(readlabel(), "_" => " "))
            len = readlength()
            nodeid[] += 1
            push!(paths, [(nodeid[], len)])
            return [length(tips)]
        end
    end
    subtree(Tuple{Int,Float64}[])
    n = length(tips)
    V = zeros(n, n)
    for i in 1:n, j in i:n
        shared = 0.0
        for (a, b) in zip(paths[i], paths[j])
            a[1] == b[1] || break
            shared += a[2]
        end
        V[i, j] = V[j, i] = shared
    end
    d = sqrt.(LinearAlgebra.diag(V))
    return (species = tips, correlation = V ./ (d * d'))
end
