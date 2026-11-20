# DiRenzo et al. (2014) Table 1: terminal survival time, infection intensity and zoospore
# output of the eight inoculated Atelopus zeteki (data/direnzo2014/table1.csv, CC BY 4.0,
# transcribed from the article; no time series were published). Summary statistics and the
# log-rank test of prior exposure, to check the paper's text against its own table.

"""
    DiRenzoFrog(prior, days, intensity, output)

One row of DiRenzo et al. (2014) Table 1: `prior` (`true` if previously exposed to JEL
427-P39), days survived after inoculation, infection intensity at death (zoospore genomic
equivalents per swab) and zoospore output at death (zoospores per minute).
"""
struct DiRenzoFrog
    prior::Bool
    days::Float64
    intensity::Float64
    output::Float64
end

"""
    direnzo2014_path()

Path of the bundled Table 1 of DiRenzo et al. (2014) (CC BY 4.0; see
`data/direnzo2014/README.md`).
"""
direnzo2014_path() = joinpath(pkgdir(@__MODULE__), "data", "direnzo2014", "table1.csv")

"""
    load_direnzo2014() -> Vector{DiRenzoFrog}

The eight infected frogs of DiRenzo et al. (2014) Table 1, in printed order.
"""
function load_direnzo2014()
    rows, col = _read_rows(direnzo2014_path())
    ip, id, ii, io = col("prior_exposure"), col("days_survived"), col("intensity_at_death"),
        col("output_at_death")
    return [DiRenzoFrog(r[ip] == "Previous", parse(Float64, r[id]), parse(Float64, r[ii]),
        parse(Float64, r[io])) for r in rows]
end

"""
    direnzo_summary(frogs = load_direnzo2014()) -> NamedTuple

Mean, standard error (sample SD/√n), minimum and maximum of days survived, intensity at death
and output at death, each as a `NamedTuple` `(mean, se, min, max)`.
"""
function direnzo_summary(frogs::AbstractVector{DiRenzoFrog} = load_direnzo2014())
    s(v) = (mean = mean(v), se = std(v) / sqrt(length(v)), min = minimum(v), max = maximum(v))
    return (days = s([f.days for f in frogs]), intensity = s([f.intensity for f in frogs]),
        output = s([f.output for f in frogs]))
end

"""
    direnzo_logrank(frogs = load_direnzo2014()) -> (; chisq, p)

Log-rank test of survival time between naive and previously exposed frogs (no censoring),
with its χ² (1 df) p-value. The paper reports χ² = 0.7, p = 0.40 (Fig. 1).
"""
function direnzo_logrank(frogs::AbstractVector{DiRenzoFrog} = load_direnzo2014())
    a = [f.days for f in frogs if !f.prior]
    b = [f.days for f in frogs if f.prior]
    O = E = V = 0.0
    for t in sort!(unique(vcat(a, b)))
        na, nb = count(>=(t), a), count(>=(t), b)
        da, d = count(==(t), a), count(==(t), a) + count(==(t), b)
        n = na + nb
        O += da
        E += d * na / n
        n > 1 && (V += d * (na / n) * (nb / n) * (n - d) / (n - 1))
    end
    χ² = (O - E)^2 / V
    return (chisq = χ², p = ccdf(Chisq(1), χ²))
end
