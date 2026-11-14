# Louca et al. (2014) discrepancy check (see wiki card): run with julia --project.
# Result (2026-11-14, n = 600 per row): extinction 0.135–0.225 and T_A within ±0.06 for every reading.
# Hypothesis: the duplicated "r" in the ESM §5.2 log-uniform list is τm or τbr, not Tbr.
# Extinction fraction (Z_o = 0) and selected risk effects for each reading, 600 runs each.
using Chytridiomycosis, Random
base = (:Zo, :μZ, :μSA, :μST, :βA, :βT, :TA, :TT, :Tm, :r, :KA, :KT)
n = 600
for (lab, extra) in (("Tbr", (:Tbr,)), ("τm", (:τm,)), ("τbr", (:τbr,)), ("none", ()))
    lp = (base..., extra...)
    for (k, tad) in ((3, :suffer), (4, :tolerate))
        rng = Xoshiro(20132783 + k)
        res = map(1:n) do _
            s = louca_sample(rng; tadpoles = tad, log_params = lp)
            sol = solve(LoucaDOIProblem(s.p; saveat = 3650.0); dt = 0.05)
            (; u = s.u, Qo = reproduction_quotient(s.p), outcome = sol.outcome)
        end
        re = louca_risk_effects(res)
        ext = count(r -> r.outcome === :extinction, res) / n
        println(rpad(lab, 5), " ", rpad(string(tad), 9), " ext ", round(ext; digits = 3),
                "  TA ", round(re.TA; digits = 2), " KT ", round(re.KT; digits = 2),
                " KA ", round(re.KA; digits = 2), " βA ", round(re.βA; digits = 2),
                " μZ ", round(re.μZ; digits = 2))
        flush(stdout)
    end
end
