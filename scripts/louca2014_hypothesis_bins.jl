# Louca et al. (2014) discrepancy check (see wiki card): run with julia --project.
# Result (2026-11-14, n = 300, suffer, Z_o = 0): m = 100 ext 0.210; m = 1000 ext 0.223, 276/300
# outcomes agree; risk effects T_A −0.01/−0.03, K_T 0.04/0.01, β_A 0.48/0.47.
# Hypothesis: outcomes depend on the number of DOI bins m (ESM §5 says they do not).
using Chytridiomycosis, Random
C = Chytridiomycosis
const SEED = 20132783
n = 300
for (k, tad) in ((3, :suffer), (4, :tolerate))
    rng = Xoshiro(SEED + k)
    S = [louca_sample(rng; tadpoles = tad, reservoir = false) for _ in 1:n]
    res = Dict{Int,Vector{Symbol}}()
    for m in (100, 1000)
        t = @elapsed res[m] = [solve(LoucaDOIProblem(C.with(s.p; m); saveat = 3650.0); dt = 0.05).outcome for s in S]
        println(tad, " m=$m: ext ", count(==(:extinction), res[m]) / n, " clr ", count(==(:clearance), res[m]) / n,
                " per ", count(==(:persistence), res[m]) / n, "  ($(round(t))s)")
        flush(stdout)
    end
    println("  agree ", count(res[100] .== res[1000]), "/$n")
    # risk effect of TA under each m
    for m in (100, 1000)
        r = [(u = s.u, Qo = reproduction_quotient(s.p), outcome = o) for (s, o) in zip(S, res[m])]
        re = louca_risk_effects(r)
        println("  m=$m risk effects: TA=", round(re.TA, digits=2), " KT=", round(re.KT, digits=2), " KA=", round(re.KA, digits=2),
                " βA=", round(re.βA, digits=2), " μZ=", round(re.μZ, digits=2), " Qo=", round(re.Qo, digits=2))
    end
    flush(stdout)
end
