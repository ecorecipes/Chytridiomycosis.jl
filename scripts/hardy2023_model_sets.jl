# Refit every model of Hardy et al. (2023, J Appl Ecol, doi:10.1111/1365-2664.14310):
# the focal model sets of Appendices S1.3 (p), S1.4 (φ) and S1.5 (f), the combined set of
# Table 3, and the derived quantities of the Results. Prints each −2 log L, K and AICc next
# to the printed values. About a minute on one thread.
#
#   julia --project scripts/hardy2023_model_sets.jl
#
# Printed values are those of the version of record and its Supporting Information (identical
# in B. M. Hardy's 2023 dissertation); they are oracles typed here for the comparison, not
# inputs to any fit.
using Chytridiomycosis
using Printf

data = load_hardy2023()
printed = Dict(
    :p => Dict(:s_bd_t_bd => 27669.30, :s_bd => 28676.37, :t_bd => 28669.52, :site => 28905.69,
               :constant => 29334.27),
    :φ => Dict(:elev_t_bd => 28046.40, :t_bd => 28085.86, :elev_bd => 28225.02,
               :elev_x_bd => 28225.02, :bd => 28270.37, :elev => 29294.13, :constant => 29334.27),
    :f => Dict(:t_bd => 28570.08, :elev_t_bd => 28570.00, :drought_bd => 28750.73,
               :drought_x_bd => 28750.72, :bd => 28756.20, :elev_bd => 28755.88,
               :winstart_bd => 28756.20, :elev_x_bd => 28754.36, :winstart_x_bd => 28756.19,
               :winstart_elev => 29266.9, :winstart => 29278.102, :drought_elev => 29323.899,
               :drought => 29329.822, :elev => 29330.671, :constant => 29334.271))

for (par, app) in ((:p, "S1.3"), (:φ, "S1.4"), (:f, "S1.5"))
    println("\nAppendix $app: focal parameter $par")
    @printf("  %-15s %10s %10s %8s %4s %10s %6s\n", "structure", "-2logL", "printed", "diff",
            "K", "AICc", "weight")
    for m in hardy2023_model_set(par; data)
        pd = printed[par][m.structure]
        @printf("  %-15s %10.2f %10.2f %8.2f %4d %10.2f %6.2f\n", m.structure, m.deviance, pd,
                m.deviance - pd, m.k, m.aicc, m.weight)
    end
end

println("\nTable 3 (printed AICc 26649.10 and 26705.91)")
top = hardy2023_fit(; data)
second = hardy2023_fit(; f = :t_bd, data)
for (nm, r) in (("p(sBd×tBd) φ(ELEV+tBd) f(ELEV+tBd)", top), ("p(sBd×tBd) φ(ELEV+tBd) f(tBd)", second))
    @printf("  %-36s AICc %9.2f  ΔAICc %6.2f  K %3d  -2logL %9.2f\n", nm, r.aicc,
            r.aicc - top.aicc, r.k, r.deviance)
end

s = hardy2023_summary(top; data)
println("\nResults (top model)")
@printf("  post-Bd φ = %.3g, f = %.3f (SE %.3f), λ = %.3f (95%% CI %.2f-%.2f)\n",
        s.φ_post.value, s.f_post.value, s.f_post.se, s.λ_post.value, s.λ_post.ci...)
lo, hi = s.λ_pre_range
@printf("  pre-Bd λ range %.2f (SE %.2f, %s %d) to %.2f (SE %.2f, %s %d)\n", lo.value, lo.se,
        lo.site, lo.year, hi.value, hi.se, hi.site, hi.year)
@printf("  p range %.3f-%.3f\n", s.p_range...)
fig = hardy2023_fit(; φ = :elev_bd, f = :elev_bd, data)
for c in (:φ_elev, :f_elev)
    r = pradel_coefficient(fig, c)
    @printf("  Fig. 1a model: %s = %.3f (SE %.3f)\n", c, r.value, r.se)
end
for g in eachindex(data)
    r = pradel_se(fig, θ -> pradel_derived(fig.problem, θ, g).φ[1])
    @printf("  pre-Bd φ at %-16s %.3f (SE %.3f)\n", data[g].name, r.value, r.se)
end
