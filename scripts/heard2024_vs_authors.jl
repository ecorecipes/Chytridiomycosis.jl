# Heard et al. (2024): compare our metapopulation simulations with the authors' archived
# outputs, scenario group by scenario group (persistence fraction per fecundity value `mets`,
# per landscape, and overall), for both species and both epochs.
#
#   julia --project scripts/heard2024_vs_authors.jl [path/to/zenodo]
#
# The path defaults to ENV["CHYTRIDIOMYCOSIS_HEARD2024"] or lab/papers/heard2024/zenodo and
# must hold sim_results_{Lran,Laur}_{historical,contemporary}_metapop_compiled.csv from
# Zenodo doi:10.5281/zenodo.10457229 (38 MB each, not bundled). About 1-2 minutes per
# species and epoch on one thread. Writes nothing.
using Chytridiomycosis, Random, Printf

dir = length(ARGS) >= 1 ? ARGS[1] : get(ENV, "CHYTRIDIOMYCOSIS_HEARD2024",
    joinpath(pkgdir(Chytridiomycosis), "..", "..", "papers", "heard2024", "zenodo"))

function authors(file)
    lines = eachline(file)
    hdr = split(first(lines), ',')
    col(n) = findfirst(==("\"$n\""), hdr)
    ie, im, il = col("extinct"), col("mets"), col("raster")
    rows = [(parse(Float64, f[im]), parse(Int, f[il]), f[ie] == "0")
            for f in (split(l, ',') for l in lines)]
    return rows
end

g = heard_scenario_grid()
for (sp, code) in ((:raniformis, "Lran"), (:aurea, "Laur")), ep in (:historical, :contemporary)
    file = joinpath(dir, "sim_results_$(code)_$(ep)_metapop_compiled.csv")
    isfile(file) || (println("missing $file"); continue)
    a = authors(file)
    rng = Xoshiro(1)
    println("\n$sp, $ep")
    tot_a = tot_o = 0.0
    for m in g.mets
        pa = count(r -> r[1] == m && r[3], a) / count(r -> r[1] == m, a)
        po = heard_persistence_grid(; species = sp, epoch = ep, mets = [m], rng)
        tot_a += pa / length(g.mets)
        tot_o += po / length(g.mets)
        @printf("  mets %4.0f   authors %.3f   ours %.3f\n", m, pa, po)
    end
    @printf("  overall    authors %.4f   ours %.4f\n", tot_a, tot_o)
end
