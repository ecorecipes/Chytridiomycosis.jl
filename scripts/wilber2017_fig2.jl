# Wilber et al. (2017) Fig. 2: ln R0 with and without the zoospore pool, against the curves
# decoded from the figure's vector graphics (validation targets wilber2017_fig2*).
#
# The figure plots the MEDIAN of R0 over 500 posterior draws (R0_environmental.R: one random
# posterior row per vital-rate function and per transmission function, independently). This
# script compares
#   (1) R0 at Table 1 / Table 2 point values (plug-in);
#   (2) R0 with our posterior means of the vital rates (what the authors' code plugs in);
#   (3) the median over 500 draws: vital rates from our refitted posteriors
#       (posterior_draws), β0 and β1 lognormal with the Table 2 means and the posterior CVs
#       printed in the preprint appendix (cv_β0, cv_β1; compendium, verified = false).
#
#   julia --project scripts/wilber2017_fig2.jl
using Chytridiomycosis, Random
using Distributions: LogNormal, median
C = Chytridiomycosis

p = Wilber2017Parameters(:wilber2017; transmission = :zpool_dd)
tr = parameter_set(:wilber2017; set = "zpool_dd")
lnR(q, T, S) = (r = basic_reproduction_number(q; temperature = T, S); (log(r.total), log(r.host)))

comps = (:survival, :growth, :loss, :initial_load)
models = Dict(c => Wilber2017Regression(c) for c in comps)
means = Dict(c => posterior_summary(models[c]; draws = 100_000).mean for c in comps)
withvr(q, s, g, l, i) = C._with(q; b00 = s[1], b10 = s[2], b01 = g[1], b11 = g[2], b21 = g[3],
    v01 = g[4], c01 = g[5], b02 = l[1], b12 = l[2], b22 = l[3], b03 = i[1], b13 = i[2],
    v03 = i[3], c03 = i[4])
pm = withvr(p, (means[c] for c in comps)...)

ndraw = round(Int, parameter_set(:wilber2017).n_R0_draws)
rng = Xoshiro(2)
D = Dict(c => posterior_draws(models[c], ndraw; rng, draws = 100_000) for c in comps)
lnorm(m, cv) = (σ = sqrt(log1p(cv^2)); LogNormal(log(m) - σ^2 / 2, σ))
β0s = rand(rng, lnorm(tr.β0, tr.cv_β0), ndraw)
β1s = rand(rng, lnorm(tr.β1, tr.cv_β1), ndraw)
draws = [C._with(withvr(p, (D[c][k, :] for c in comps)...); β0 = β0s[k], β1 = β1s[k])
         for k in 1:ndraw]
medR(T, S) = (v = [lnR(q, T, S) for q in draws]; (median(first.(v)), median(last.(v))))

println("Fig. 2c (S = 4) and 2b (15 °C): ln R0 total / host")
println("target                         published   plug-in(T1)   post.means   median(draws)")
# Registered targets, and four more points decoded from the figure (the PDF's vector paths) that
# are not registered: (T, S, kind, decoded ln R0).
cases = [(validation_target(id), id) for id in ("wilber2017_fig2c_total_4C",
    "wilber2017_fig2c_total_15C", "wilber2017_fig2c_host_20C", "wilber2017_fig2b_total_S1",
    "wilber2017_fig2b_host_S1")]
point(id) = (occursin("fig2c", id) ? parse(Float64, match(r"_(\d+)C$", id)[1]) : 15.0,
    occursin("fig2b", id) ? parse(Float64, match(r"_S(\d+)$", id)[1]) : 4.0,
    occursin("total", id) ? 1 : 2)
rows = [(id, point(id)..., t.published) for (t, id) in cases]
rows = [(r[1], r[2] == 15 && occursin("fig2c", r[1]) ? 4 + 13 * 16 / 19 : r[2], r[3:end]...)
        for r in rows]                                  # the grid point 14.947 °C
append!(rows, [("fig2c total 20C", 20.0, 4.0, 1, 2.281), ("fig2c host 4C", 4.0, 4.0, 2, -0.402),
    ("fig2b total S16", 15.0, 16.0, 1, 3.927), ("fig2b host S16", 15.0, 16.0, 2, 2.850)])
for (id, T, S, k, published) in rows
    vals = (lnR(p, T, S)[k], lnR(pm, T, S)[k], medR(T, S)[k])
    println(rpad(id, 30), lpad(round(published; digits = 3), 10),
        join(lpad.(round.(vals; digits = 3), 13)))
end
