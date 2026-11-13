# Shared helpers for the Drawert et al. (2017) reproduction scripts
# (scripts/drawert2017_*.jl). Included, not run on its own.
using Chytridiomycosis, Random, Statistics, Printf
using Distributions: Hypergeometric, TDist, cdf, pdf

const C = Chytridiomycosis
const OUTDIR = normpath(joinpath(@__DIR__, "..", "data", "drawert2017"))
const SEED = 20170480                      # the paper's article number
# Points per group (Fig. 3) in the stratified Baseline sample.
const BASELINE_SIZES = (A = 25, B = 15, C = 30, D = 40, E = 6)

# Seeded per (point, arm) so that results do not depend on threading or on run order.
point_rng(id, arm::Integer = 0) = Xoshiro(SEED + 10_000 * arm + id)

# Two-sided Fisher exact test for a 2×2 table [a b; c d] (rows: ours, theirs).
function fisher_p(a, b, c, d)
    n1, n2, k = a + b, c + d, a + c
    h = Hypergeometric(n1, n2, k)
    p0 = pdf(h, a)
    lo, hi = max(0, k - n2), min(k, n1)
    return min(1.0, sum(pdf(h, x) for x in lo:hi if pdf(h, x) <= p0 * (1 + 1e-7)))
end

# Welch two-sample test from means and population sds (divisor n), two-sided p-value.
function welch_p(m1, s1, n1, m2, s2, n2)
    v1 = s1^2 * n1 / (n1 - 1) / n1
    v2 = s2^2 * n2 / (n2 - 1) / n2
    v = v1 + v2
    v == 0 && return m1 == m2 ? 1.0 : 0.0
    t = (m1 - m2) / sqrt(v)
    df = v^2 / (v1^2 / (n1 - 1) + v2^2 / (n2 - 1))
    return 2 * (1 - cdf(TDist(df), abs(t)))
end

# Simple CSV I/O for the per-point result files.
function read_rows(path)
    isfile(path) || return Dict{String,Vector{String}}[]
    lines = filter(!isempty ∘ strip, readlines(path))
    isempty(lines) && return Dict{String,String}[]
    h = split(lines[1], ',')
    return [Dict(String(h[j]) => String(v[j]) for j in eachindex(h)) for v in split.(lines[2:end], ',')]
end

function append_row(path, header, row, lk)
    lock(lk) do
        new = !isfile(path)
        open(path, "a") do io
            new && println(io, join(header, ','))
            println(io, join(row, ','))
        end
    end
end

# Stratified random sample of the authors' points by their Fig. 3 group (seeded).
function stratified_sample(pts, sizes::NamedTuple; seed = SEED)
    rng = Xoshiro(seed)
    groups = [survival_group(r.mean_surviving, r.std_surviving) for r in pts]
    out = Int[]
    weights = Dict{Symbol,Float64}()
    for g in keys(sizes)
        ids = [r.id for (r, gg) in zip(pts, groups) if gg === g]
        n = min(sizes[g], length(ids))
        append!(out, sort!(shuffle(rng, ids)[1:n]))
        weights[g] = length(ids) / n
    end
    return out, weights, Dict(r.id => g for (r, g) in zip(pts, groups))
end

fmt(x::Bool) = string(x)
fmt(x::Real) = isinteger(x) ? string(Int(x)) : @sprintf("%.6g", x)
fmt(x) = string(x)
# Booleans were first written as 0/1; accept both spellings.
parsebool(s::AbstractString) = s in ("true", "1")
