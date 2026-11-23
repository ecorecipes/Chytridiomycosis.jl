# Field-trial and challenge-experiment analyses of Barnett et al. (2026, bioRxiv PREPRINT,
# doi:10.64898/2026.05.15.725113), refitted to the authors' data
# (https://github.com/kmbarnett100/protection_to_amplification; not bundled, no licence).
# This reproduces a published analysis; it is not a new model.
#
# The authors' Field_and_Challenge_Analysis_and_Figures.R fits four before–after/control–impact
# (BACI) mixed models with glmmTMB. The fixed effects are treatment × time with R's default
# treatment contrasts (reference: sham, after), and the group effects are independent normals:
#   field intensity      round(swab1) ~ zero-inflated NB2, log link; group effects pond, year;
#                        the same predictors in the zero-inflation (logit) part
#   field prevalence     infected ~ binomial, logit link; group effects pond, year
#   challenge intensity  round(swab2) ~ zero-inflated NB2; group effect pond
#   challenge prevalence bd_end ~ binomial; group effect pond
# As in glmmTMB, the likelihood is the Laplace approximation with the observed Hessian, and
# standard errors come from the Hessian of that likelihood. Rows with a missing response are
# dropped. The "marginal means" are modelbased::estimate_means(..., estimate = "average"): the
# conditional linear predictor (with the group modes) averaged over the observed rows of each
# treatment × time cell, with the standard error of the fixed-effect part and intervals
# mean ± 1.96 SE, back-transformed where the response scale is wanted.

"""
    barnett2026_path() -> String

Folder with the authors' repository for Barnett et al. (2026, **bioRxiv preprint**),
<https://github.com/kmbarnett100/protection_to_amplification> (NetLogo model, BehaviorSpace
output, field and challenge data, R code). It is not bundled, because the repository has no
licence. The path is `ENV["CHYTRIDIOMYCOSIS_BARNETT2026"]` if set, else
`papers/barnett2026/github` under the first ancestor of the package directory that has it.
Throws an `ArgumentError` if neither exists.
"""
function barnett2026_path()
    p = _barnett2026_find()
    p === nothing && throw(ArgumentError(
        "Barnett et al. 2026 files not found: clone " *
        "https://github.com/kmbarnett100/protection_to_amplification and set " *
        "ENV[\"CHYTRIDIOMYCOSIS_BARNETT2026\"] to that folder"))
    return p
end

function _barnett2026_find()
    ok(c) = isfile(joinpath(c, "field data", "raw", "all_years_field_swab_data_clean.csv"))
    e = get(ENV, "CHYTRIDIOMYCOSIS_BARNETT2026", "")
    isempty(e) || return ok(e) ? e : nothing
    d = pkgdir(@__MODULE__)
    while true
        c = joinpath(d, "papers", "barnett2026", "github")
        ok(c) && return c
        parent = dirname(d)
        parent == d && return nothing
        d = parent
    end
end

"""
    load_barnett2026(; path = barnett2026_path()) -> (field, challenge)

The authors' field swabs (1672 *Pseudacris regilla* metamorphs, 12 ponds, 2011–2022) and
challenge experiment (frogs with a negative initial swab). Each is a `NamedTuple` of vectors
`pond`, `year`, `vaccine` (metabolite-treated pond), `before` (before the intervention),
`load` (Bd genomic equivalents, unrounded) and `infected`; rows with a missing final swab
(challenge only) are dropped, as glmmTMB does.
"""
function load_barnett2026(; path::AbstractString = barnett2026_path())
    function read(file, ycol, icol)
        rows, col = _gj_csv(joinpath(path, "field data", "raw", file))
        ip, iy, it, ib, il, ii = col("pond"), col("year"), col("treatment"),
            col("before.after"), col(ycol), col(icol)
        keep = filter(r -> r[il] != "NA" && r[ii] != "NA", rows)
        return (pond = [String(r[ip]) for r in keep], year = [parse(Int, r[iy]) for r in keep],
                vaccine = [r[it] == "vaccine" for r in keep],
                before = [r[ib] == "before" for r in keep],
                load = [parse(Float64, r[il]) for r in keep],
                infected = [parse(Int, r[ii]) == 1 for r in keep])
    end
    return (field = read("all_years_field_swab_data_clean.csv", "swab1", "infected"),
            challenge = read("challenge_uninfected.csv", "swab2", "bd_end"))
end

const BARNETT_BACI = (:field_intensity, :field_prevalence, :challenge_intensity,
                      :challenge_prevalence)

"""
    BarnettBACIModel(kind; data = load_barnett2026())

One of the four BACI mixed models of Barnett et al. (2026, **bioRxiv preprint**): `kind` is
`:field_intensity`, `:field_prevalence`, `:challenge_intensity` or `:challenge_prevalence`.
Intensity models are zero-inflated negative binomial (NB2, log link, the same predictors in
the zero-inflation part); prevalence models are binomial (logit). Fixed effects are treatment
× time; the group effects are pond and year (field) or pond (challenge). Fit with
[`maximum_likelihood`](@ref); the log-likelihood is the Laplace approximation used by glmmTMB.

Parameters, in order: the four fixed effects (intercept = sham after, vaccine, before,
interaction) and the log standard deviations of the group effects of the count (or binomial)
part; for intensity models, the same for the zero-inflation part, then the log dispersion.
"""
struct BarnettBACIModel
    kind::Symbol
    zinb::Bool
    y::Vector{Float64}              # unique rows
    w::Vector{Float64}              # multiplicities
    X::Matrix{Float64}
    groups::Vector{Vector{Int}}     # one index vector per group factor
    sizes::Vector{Int}
    cell::Vector{Int}               # treatment × time cell of each row (1:4)
end

# One unique row per (response, cell, groups), with its multiplicity.
function BarnettBACIModel(kind::Symbol; data = load_barnett2026())
    kind in BARNETT_BACI || throw(ArgumentError("kind must be one of $BARNETT_BACI; got :$kind"))
    field = kind in (:field_intensity, :field_prevalence)
    zinb = kind in (:field_intensity, :challenge_intensity)
    d = field ? data.field : data.challenge
    y = zinb ? round.(d.load) : Float64.(d.infected)
    facs = field ? (d.pond, d.year) : (d.pond,)
    idx = map(facs) do f
        lev = sort(unique(f))
        [findfirst(==(v), lev) for v in f]
    end
    counts = Dict{Any,Int}()
    for i in eachindex(y)
        k = (y[i], d.vaccine[i], d.before[i], map(g -> g[i], idx)...)
        counts[k] = get(counts, k, 0) + 1
    end
    ks = sort(collect(keys(counts)))
    X = reduce(vcat, [[1.0 k[2] k[3] k[2] * k[3]] for k in ks])
    cell = [1 + Int(k[3]) + 2Int(k[2]) for k in ks]   # sham after, sham before, vaccine …
    groups = [[k[3 + g] for k in ks] for g in eachindex(idx)]
    return BarnettBACIModel(kind, zinb, [k[1] for k in ks], [Float64(counts[k]) for k in ks], X,
                            groups, [maximum(g) for g in idx], cell)
end

_baci_dim(m::BarnettBACIModel) =
    m.zinb ? 2 * (4 + length(m.groups)) + 1 : 4 + length(m.groups)

_baci_log1pexp(x) = x > 0 ? x + log1p(exp(-x)) : log1p(exp(x))

# Log-likelihood of one row and its derivatives in (ηc, ηz):
# returns ℓ, ∂c, ∂z, ∂cc, ∂cz, ∂zz.
@inline function _baci_row(zinb::Bool, y, ηc, ηz, k, lc)
    if !zinb
        p = inv(1 + exp(-ηc))
        return y * ηc - _baci_log1pexp(ηc), y - p, 0.0, -p * (1 - p), 0.0, 0.0
    end
    μ = exp(min(ηc, 600.0))
    π0 = inv(1 + exp(-ηz))
    r = μ / (k + μ)
    lq = -log1p(μ / k)                      # log(k/(k + μ))
    l1π = -_baci_log1pexp(ηz)               # log(1 − π0)
    if y > 0
        ℓ = l1π + lc + k * lq - y * log1p(k / μ)   # y log(μ/(k + μ)), without cancellation
        return ℓ, y - (y + k) * r, -π0, -(y + k) * r * (1 - r), 0.0, -π0 * (1 - π0)
    end
    lf = k * lq
    g = -k * r
    h = -k * r * (1 - r)
    q = inv(1 + exp(lf - ηz))               # posterior probability of a structural zero
    v = q * (1 - q)
    ℓ = l1π + max(ηz, lf) + log1p(exp(-abs(ηz - lf)))
    return ℓ, (1 - q) * g, q - π0, (1 - q) * h + v * g^2, -v * g, v - π0 * (1 - π0)
end

# Newton search for the standardised group effects u. Returns the Laplace log-likelihood, u,
# the fixed parts of the linear predictors and the group standard deviations.
function _baci_laplace(m::BarnettBACIModel, θ::AbstractVector, u0::AbstractVector)
    ng = length(m.groups)
    np = m.zinb ? 2 : 1
    q1 = sum(m.sizes)
    nu = np * q1
    offs = cumsum([0; m.sizes[1:(end - 1)]])
    n = length(m.y)
    step = 4 + ng
    β = [θ[(step * (c - 1) + 1):(step * (c - 1) + 4)] for c in 1:np]
    sd = [exp.(θ[(step * (c - 1) + 5):(step * c)]) for c in 1:np]
    k = m.zinb ? exp(θ[end]) : 1.0
    # log Γ(y + k)/(Γ(k) y!), from the negative binomial at p = 1/2
    lc = m.zinb ? [logpdf(NegativeBinomial(k, 0.5), y) + (k + y) * log(2) for y in round.(Int, m.y)] :
         zeros(n)
    fix = [m.X * β[c] for c in 1:np]
    col(c, g, i) = (c - 1) * q1 + offs[g] + m.groups[g][i]
    function joint(u, deriv)
        f = -0.5 * sum(abs2, u)
        G = deriv ? -copy(u) : u
        H = deriv ? Matrix{Float64}(-I, nu, nu) : zeros(0, 0)
        for i in 1:n
            ηc = fix[1][i]
            ηz = np == 2 ? fix[2][i] : 0.0
            for g in 1:ng
                ηc += sd[1][g] * u[col(1, g, i)]
                np == 2 && (ηz += sd[2][g] * u[col(2, g, i)])
            end
            ℓ, d1, d2, h11, h12, h22 = _baci_row(m.zinb, m.y[i], ηc, ηz, k, lc[i])
            w = m.w[i]
            f += w * ℓ
            deriv || continue
            ds = (d1, d2)
            hs = ((h11, h12), (h12, h22))
            for c in 1:np, g in 1:ng
                a = col(c, g, i)
                G[a] += w * ds[c] * sd[c][g]
                for c2 in 1:np, g2 in 1:ng
                    H[a, col(c2, g2, i)] += w * hs[c][c2] * sd[c][g] * sd[c2][g2]
                end
            end
        end
        return f, G, H
    end
    u = copy(u0)
    f, G, H = joint(u, true)
    # Damped Newton: the damping τ rises until the step improves f (the zero-inflated part is
    # not concave), and is 0 near the mode.
    τ = 0.0
    for _ in 1:200
        F = cholesky(Symmetric(-H + τ * I); check = false)
        if !issuccess(F)
            τ = max(4τ, 1e-2)
            continue
        end
        s = F \ G
        dec = dot(G, s)
        last = τ == 0 && dec < 1e-11                     # quadratic: this step is exact
        fn = last ? f : joint(u + s, false)[1]
        if fn >= f || dec < 1e-9                         # below that f is flat to rounding
            u += s
            f, G, H = joint(u, true)
            last && break
            τ = τ < 1e-2 ? 0.0 : τ / 8
        else                                             # worse, or NaN
            τ = max(4τ, 1e-2)
        end
    end
    F = cholesky(Symmetric(-H); check = false)
    ll = issuccess(F) ? f - 0.5 * logdet(F) : -Inf
    return ll, u, fix, sd
end

function parameter_names(m::BarnettBACIModel)
    fx = (:intercept, :vaccine, :before, :vaccine_before)
    gr = length(m.groups) == 2 ? (:log_sd_pond, :log_sd_year) : (:log_sd_pond,)
    c = (fx..., gr...)
    return m.zinb ? (c..., map(s -> Symbol(:zi_, s), c)..., :log_dispersion) : c
end

LogDensityProblems.dimension(m::BarnettBACIModel) = _baci_dim(m)
LogDensityProblems.capabilities(::Type{BarnettBACIModel}) =
    LogDensityProblems.LogDensityOrder{0}()
LogDensityProblems.logdensity(m::BarnettBACIModel, θ) =
    _baci_laplace(m, θ, zeros((m.zinb ? 2 : 1) * sum(m.sizes)))[1]

function _baci_start(m::BarnettBACIModel)
    ng = length(m.groups)
    W = sum(m.w)
    pos = m.y .> 0
    p = clamp(sum(m.w[pos]) / W, 0.05, 0.95)
    m.zinb || return [log(p / (1 - p)); zeros(3 + ng)]
    μ = sum(m.w[pos] .* m.y[pos]) / sum(m.w[pos])
    return [log(μ) - 1; zeros(3 + ng); log((1 - p) / p); zeros(3 + ng); -1.0]
end

"""
    maximum_likelihood(m::BarnettBACIModel; x0, iters = 400) -> NamedTuple

Maximum (Laplace) likelihood by BFGS with central-difference gradients. Returns `loglik`,
`x` (in the order of `parameter_names(m)`), `estimates`, `se`, `interaction` (estimate, SE,
z and two-sided P of the treatment × time term of the count or binomial part), `gradient`
(largest absolute gradient at the maximum) and `means`: for `sham_after`, `sham_before`,
`vaccine_after` and `vaccine_before`, the marginal mean on the link scale (`link`), its `se`,
the 95% interval on the link scale (`lower`, `upper`) and the back-transformed `response`,
`response_lower` and `response_upper`.
"""
function maximum_likelihood(m::BarnettBACIModel; x0::AbstractVector = _baci_start(m),
        iters::Integer = 400)
    p = _baci_dim(m)
    u = Ref(zeros((m.zinb ? 2 : 1) * sum(m.sizes)))     # warm start of the inner search
    f(x) = -_baci_laplace(m, x, u[])[1]
    function grad(x, h = 1e-5)
        g = similar(x)
        for i in 1:p
            a = copy(x); a[i] += h
            b = copy(x); b[i] -= h
            g[i] = (f(a) - f(b)) / 2h
        end
        return g
    end
    function hess(x, h = 5e-3)
        H = zeros(p, p)
        for j in 1:p
            a = copy(x); a[j] += h
            b = copy(x); b[j] -= h
            H[:, j] = (grad(a, 5e-4) - grad(b, 5e-4)) / 2h
        end
        return (H + H') / 2
    end
    # step along s (at most 3 long) with backtracking; nothing if f does not fall
    function search(x, fx, s)
        norm(s) > 3 && (s = s * (3 / norm(s)))
        t = 1.0
        while t > 1e-6
            fn = f(x + t * s)
            fn < fx && return x + t * s, fn
            t /= 2
        end
        return nothing
    end
    x = float.(copy(x0))
    fx = f(x)
    g = grad(x)
    B = Matrix{Float64}(I, p, p) / max(1.0, norm(g))
    for _ in 1:iters                                    # BFGS
        r = search(x, fx, -B * g)
        r === nothing && break
        dx = r[1] - x
        drop = fx - r[2]
        x, fx = r
        u[] = _baci_laplace(m, x, u[])[2]
        gn = grad(x)
        dg = gn - g
        g = gn
        ρ = dot(dx, dg)
        if ρ > 1e-12
            V = I - dx * dg' / ρ
            B = V * B * V' + dx * dx' / ρ
        end
        drop < 1e-9 && break
    end
    H = hess(x)
    for _ in 1:25                                       # Newton polish
        g = grad(x)
        e = eigen(Symmetric(H))
        s = -e.vectors * ((e.vectors' * g) ./ max.(abs.(e.values), 1e-4))
        abs(dot(s, g)) < 1e-10 && break
        r = search(x, fx, s)
        r === nothing && break
        x, fx = r
        u[] = _baci_laplace(m, x, u[])[2]
        H = hess(x)
    end
    g = grad(x)
    V = inv(Symmetric(H))
    se = sqrt.(max.(diag(V), 0.0))
    ll, umode, fix, sd = _baci_laplace(m, x, u[])
    ng = length(m.groups)
    offs = cumsum([0; m.sizes[1:(end - 1)]])
    η = [fix[1][i] + sum(sd[1][g] * umode[offs[g] + m.groups[g][i]] for g in 1:ng)
         for i in eachindex(m.y)]
    inv_link(z) = m.zinb ? exp(z) : inv(1 + exp(-z))
    z975 = quantile(Normal(), 0.975)
    means = NamedTuple{(:sham_after, :sham_before, :vaccine_after, :vaccine_before)}(
        ntuple(4) do c
            r = findall(==(c), m.cell)
            link = sum(m.w[r] .* η[r]) / sum(m.w[r])
            xc = m.X[r[1], :]
            s = sqrt(dot(xc, V[1:4, 1:4] * xc))
            lo, hi = link - z975 * s, link + z975 * s
            (; link, se = s, lower = lo, upper = hi, response = inv_link(link),
               response_lower = inv_link(lo), response_upper = inv_link(hi))
        end)
    z = x[4] / se[4]
    names = parameter_names(m)
    return (; loglik = ll, x, estimates = NamedTuple{names}(Tuple(x)),
            se = NamedTuple{names}(Tuple(se)), means,
            interaction = (; estimate = x[4], se = se[4], z, p = 2 * ccdf(Normal(), abs(z))),
            gradient = maximum(abs, g))
end

# The fits take up to half a minute; keep them for the session (validation targets, tests).
const _BARNETT_BACI_FITS = Dict{Symbol,Any}()
_barnett_baci(kind::Symbol) =
    get!(() -> maximum_likelihood(BarnettBACIModel(kind)), _BARNETT_BACI_FITS, kind)

# The stored comparison of the ABM with the authors' runs (scripts/barnett2026_abm.jl).
_barnett_reproduction() =
    TOML.parsefile(joinpath(pkgdir(@__MODULE__), "data", "barnett2026", "reproduction.toml"))
