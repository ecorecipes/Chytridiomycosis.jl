# Refit of the model set of Mosher et al. (2018) Table 3 to the bundled Dryad data
# (data/mosher2018): OUR maximum-likelihood fits of a published model, not published
# results. Prints K, deviance (−2 log L), AICc (n = site-years surveyed) and ΔAICc for the
# ten best models and the two reference models of Table 3, the estimates the paper quotes,
# and the elevation curves of Fig. 5. About two minutes.
#
#   julia --project scripts/mosher2018_refit.jl
using Chytridiomycosis, Printf

d = load_mosher2018()
println("sites ", length(d.sites), ", site-years surveyed ",
        count(I -> any(!ismissing, view(d.toad, I[1], I[2], :)), CartesianIndices(d.toad[:, :, 1])))
t = @elapsed sel = mosher2018_model_selection(d)
@printf("%-36s %3s %9s %9s %7s %6s\n", "model (Table 3 order)", "K", "deviance", "AICc", "ΔAICc", "conv")
for f in sort(sel; by = f -> f.index)
    @printf("%2d %-33s %3d %9.2f %9.2f %7.2f %6s\n", f.index, f.name, f.k, f.deviance, f.aicc,
            f.Δaicc, f.converged)
end
println("fitted in ", round(t; digits = 0), " s")

m0 = Mosher2018Model(d)
f0 = maximum_likelihood(m0)
p0 = mosher2018_site_parameters(m0, f0.θ, 1)
@printf("all-constant model (K = %d): εA = %.3f, εB = %.3f\n", f0.k, p0.εA, p0.εB)

m = mosher2018_model(d, 1)
f = only(filter(x -> x.index == 1, sel))
q = mosher2018_parameters_at(m, f.θ; elevation = 3000)
@printf("best model: εA = %.3f, εB = %.3f, γA = %.3f\n", q.εA, q.εB, q.γA)
println("elevation  εAB    εBA    γAB    γB     γBA    rBA")
for e in 2470:100:3680
    r = mosher2018_parameters_at(m, f.θ; elevation = e)
    @printf("%5d   %.3f  %.3f  %.3f  %.3f  %.3f  %.3f\n", e, r.εAB, r.εBA, r.γAB, r.γB, r.γBA, r.rBA)
end
