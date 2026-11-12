# A free-living zoospore pool with live and dead compartments, for any model in which Bd is
# measured by DNA (qPCR or eDNA) but infects only while alive. Structure from Rumschlag et al.
# 2022 (J Anim Ecol 91:170, eqs 1-2): live zoospores die at per capita rate m and the dead
# decompose at per capita rate d. qPCR detects DNA from both, so it counts ZL + ZD; only ZL
# infects. Wiki: lab/wiki/models/bd-rumschlag2022-zoospore-persistence.md
#
#   dZL/dt = shedding − m ZL − uptake
#   dZD/dt = m ZL − d ZD
#
# `shedding` (zoospores released per unit time) and `uptake` (live zoospores removed by
# infection per unit time) come from the host model that uses the pool; both are zero in the
# Rumschlag assays. Units are whatever the host model uses (Rumschlag: zoospores µl⁻¹, days).

"""
    ZoosporePool(m, d)

A free-living zoospore pool with live and dead compartments:

```
dZL/dt = shedding − m ZL − uptake
dZD/dt = m ZL − d ZD
```

`m` is the per capita death rate of live zoospores and `d` the per capita decomposition rate
of dead ones, both in the host model's time unit. A qPCR or eDNA assay detects DNA from both
compartments, so it sees `ZL + ZD` ([`qpcr_zoospores`](@ref)); only `ZL` can infect
([`infectious_zoospores`](@ref)). The structure is that of Rumschlag et al. (2022), eqs 1–2;
`ZoosporePool(p::RumschlagParameters)` builds the pool of one water source.

Use [`pool_derivatives`](@ref) inside another model's right-hand side, and
[`pool_state`](@ref) for the closed-form solution without shedding.
"""
struct ZoosporePool
    m::Float64
    d::Float64
    function ZoosporePool(m::Real, d::Real)
        m >= 0 && d >= 0 || throw(ArgumentError("rates must be non-negative; got m = $m, d = $d"))
        return new(m, d)
    end
end

"""
    pool_derivatives(pool::ZoosporePool, ZL, ZD; shedding = 0, uptake = 0) -> (dZL, dZD)

Time derivatives of the live (`ZL`) and dead (`ZD`) zoospores of `pool`, with `shedding`
live zoospores added and `uptake` live zoospores removed (for example by infection) per unit
time.
"""
pool_derivatives(pool::ZoosporePool, ZL, ZD; shedding = 0, uptake = 0) =
    (shedding - pool.m * ZL - uptake, pool.m * ZL - pool.d * ZD)

# (1 − exp(−x t))/x, equal to t at x = 0, written so that it stays accurate near x = 0
_expm1_ratio(x, t) = abs(x * t) < 1e-8 ? t * (1 - x * t / 2) : -expm1(-x * t) / x

"""
    pool_state(pool::ZoosporePool, ZL0, t; ZD0 = 0) -> (ZL, ZD)

Live and dead zoospores at time `t`, starting from `ZL0` live and `ZD0` dead, with no
shedding or uptake:

```
ZL(t) = ZL0 exp(−m t)
ZD(t) = ZD0 exp(−d t) + ZL0 m/(d − m) (exp(−m t) − exp(−d t))
```

The case `d = m` is the limit `ZL0 m t exp(−m t)`. With `ZD0 = 0` this is Rumschlag et
al. (2022) eqs 5–6, with the sign of eq. 6 corrected: as printed (in both the preprint and
the published version) the bracket is `exp(−d t) − exp(−m t)`, which is negative for every
`t > 0`.
"""
function pool_state(pool::ZoosporePool, ZL0, t; ZD0 = 0)
    m, d = pool.m, pool.d
    ZL = ZL0 * exp(-m * t)
    ZD = ZD0 * exp(-d * t) + ZL0 * m * exp(-m * t) * _expm1_ratio(d - m, t)
    return (ZL, ZD)
end

"""
    qpcr_zoospores(ZL, ZD)
    qpcr_zoospores(pool::ZoosporePool, ZL0, t; ZD0 = 0)

What a DNA-based count (qPCR or eDNA) of the pool sees: live plus dead zoospores,
`ZL + ZD`. It assumes that a dead zoospore keeps its DNA until it decomposes and that live
and dead zoospores carry the same DNA; neither assumption was tested by Rumschlag et
al. (2022), who counted live and dead cells by Trypan blue staining. The second form
evaluates the closed-form pool of [`pool_state`](@ref) at time `t`.
"""
qpcr_zoospores(ZL, ZD) = ZL + ZD
qpcr_zoospores(pool::ZoosporePool, ZL0, t; ZD0 = 0) = sum(pool_state(pool, ZL0, t; ZD0))

"""
    infectious_zoospores(ZL, ZD)
    infectious_zoospores(pool::ZoosporePool, ZL0, t; ZD0 = 0)

The zoospores that can infect: the live ones, `ZL`. Compare [`qpcr_zoospores`](@ref).
"""
infectious_zoospores(ZL, ZD) = ZL
infectious_zoospores(pool::ZoosporePool, ZL0, t; ZD0 = 0) = first(pool_state(pool, ZL0, t; ZD0))

"""
    viable_fraction(ZL, ZD)
    viable_fraction(pool::ZoosporePool, t)

The fraction of the DNA-detectable pool that is alive, `ZL/(ZL + ZD)`. The second form
starts from live zoospores only (`ZD0 = 0`). It does not depend on `ZL0`. It falls from 1 at
`t = 0`. When `d > m` it tends to `1 − m/d`; otherwise it tends to 0.
"""
viable_fraction(ZL, ZD) = ZL / (ZL + ZD)
viable_fraction(pool::ZoosporePool, t) = viable_fraction(pool_state(pool, 1.0, t)...)
