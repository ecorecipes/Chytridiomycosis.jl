# Wilber et al. (2017) Fig. 3: global sensitivity of extinction within 8 years to the
# transmission, resistance and tolerance parameters (1 000 lognormal perturbations,
# σ = 0.3, per transmission function), summarised by the L2 logistic regression of the
# authors' notebook. Also the extinction curves of Appendix Fig. S5 / preprint Fig. 3A.
#
# The published weights below were decoded from the vector graphics of Fig. 3 (page 9 of the
# PDF: bar heights in points, 6.13 pt per unit, zero lines at 108.96, 181.70, 254.43 pt);
# they are reference values for this comparison, not parameters. The random forest of Fig. 3
# is not reproduced (no new dependency for it).
#
#   julia --project scripts/wilber2017_sensitivity.jl
using Chytridiomycosis, Random

const FIG3 = Dict(   # b00 b10 b01 b21 b11 v01 c01 b03 b13 v03 c03 b02 b22 b12 β0 β1
    :zpool => [-1.43, -0.72, -0.15, 1.76, 6.35, 0.58, 0.25, -0.10, -0.02, 0.17, 0.15, -0.23,
        -0.20, -0.12, 0.10],
    :zpool_fd => [-1.33, -0.80, 0.15, 1.67, 6.77, 0.69, 0.31, -0.04, -0.03, 0.23, 0.09, -0.14,
        -0.26, -0.17, 0.17, -0.05],
    :zpool_dd => [-1.54, -0.51, 0.03, 1.77, 6.60, 0.47, 0.27, -0.18, -0.08, -0.03, 0.07, -0.29,
        -0.13, -0.02, 0.30, -0.09])

nsims = round(Int, parameter_set(:wilber2017).n_sensitivity)
for (k, tr) in enumerate((:zpool, :zpool_fd, :zpool_dd))
    p = Wilber2017Parameters(:wilber2017; set = "coded", transmission = tr)
    t = @elapsed r = extinction_sensitivity(p; nsims, rng = Xoshiro(100 + k))
    w = sensitivity_weights(r.X, r.extinct).weights
    println("\n$tr: $(nsims) runs in $(round(t; digits = 1)) s, extinct in ≤ 8 yr: ",
        round(sum(r.extinct) / nsims; digits = 3))
    println(rpad("param", 7), lpad("Fig. 3", 8), lpad("ours", 8))
    for (i, s) in enumerate(r.names)
        println(rpad(string(s), 7), lpad(FIG3[tr][i], 8), lpad(round(w[i]; digits = 2), 8))
    end
end

println("\nExtinction curves (500 runs, 25 years): probability by year")
q = parameter_set(:wilber2017)
for (tr, ν) in ((:zpool, NaN), (:zpool_fd, NaN), (:zpool_dd, NaN), (:zpool_dd, q.zsurv_high))
    p = Wilber2017Parameters(:wilber2017; set = "coded", transmission = tr, ν = ν)
    c = extinction_curve(HybridIPMProblem(p); rng = Xoshiro(25))
    at(y) = c.probability[findfirst(>=(y - 1e-9), c.years)]
    println(rpad("$tr ν=$ν", 22), join([lpad("$y:" * string(round(at(y); digits = 3)), 11)
                                        for y in (4, 5, 6, 7, 8, 10, 12, 15, 20, 25)]))
end
