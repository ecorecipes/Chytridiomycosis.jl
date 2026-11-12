# Woodhams, Alford, Briggs, Johnson & Rollins-Smith (2008) Ecology 89:1627. Bd population
# growth in culture: zoospores Z encyst at rate mZ and die at rate dZ; an encysted zoospore
# matures for s hours and releases F zoospores. The long-term growth rate λ is the root of
# the characteristic equation (their Eq. 7)
#
#   F mZ exp(−λ s) − λ − (mZ + dZ) = 0,
#
# where F is the number of viable zoospores per zoosporangium. Units: hours.
# Wiki: lab/wiki/models/thermal-performance.md

"""
    WoodhamsParameters(; F, s, dZ, mZ)
    WoodhamsParameters(:woodhams2008; set = "23C")

Life-history parameters of Bd in culture (Woodhams et al. 2008, Table 2):

- `F`: viable zoospores produced per zoosporangium;
- `s`: mean hours to sporulation;
- `dZ`: per capita zoospore mortality (h⁻¹);
- `mZ`: per capita zoospore encystment (h⁻¹).

The compendium holds Table 2 as two sets, `"10C"` and `"23C"`. See
[`intrinsic_growth_rate`](@ref).
"""
Base.@kwdef struct WoodhamsParameters{T<:Real}
    F::T
    s::T
    dZ::T
    mZ::T
end

function WoodhamsParameters(key::Symbol; set::AbstractString)
    p = parameter_set(key; set)
    return WoodhamsParameters(; F = p.F, s = p.s, dZ = p.dZ, mZ = p.mZ)
end

time_unit(::Type{<:WoodhamsParameters}) = :hour

"""
    intrinsic_growth_rate(p::WoodhamsParameters) -> λ (h⁻¹)

The long-term growth rate of a Bd population in culture: the real root of Woodhams et al.'s
(2008) characteristic equation, `F mZ exp(−λs) = λ + mZ + dZ`. The left side falls and the
right side rises in λ, so the root is unique. It is found by bisection and refined by
Newton's method.

The paper reports λ = 0.044 h⁻¹ at 23 °C and 0.021 h⁻¹ at 10 °C (validation targets
`woodhams2008_lambda_23C` and `woodhams2008_lambda_10C`). It is not `log(F)/s`.
"""
function intrinsic_growth_rate(p::WoodhamsParameters)
    l = p.mZ + p.dZ
    g(λ) = p.F * p.mZ * exp(-λ * p.s) - λ - l
    lo, hi = -l, p.F * p.mZ                  # g(lo) > 0 ≥ g(hi)
    g(hi) <= 0 || throw(ArgumentError("no root below F mZ"))
    for _ in 1:200
        mid = (lo + hi) / 2
        g(mid) > 0 ? (lo = mid) : (hi = mid)
        hi - lo < 1e-12 * max(1, abs(mid)) && break
    end
    λ = (lo + hi) / 2
    for _ in 1:3                              # Newton polish; g' = −s F mZ e^{−λs} − 1
        λ -= g(λ) / (-p.s * p.F * p.mZ * exp(-λ * p.s) - 1)
    end
    return λ
end
