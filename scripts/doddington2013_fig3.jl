# Doddington et al. (2013) Fig. 3: probability that the Alytes muletensis population survives
# 35 years after Bd introduction, over the 10 x 10 grid of the transmission parameter ν and
# the zoospore release rate ρ, for high (α = 4.54, panel a) and weak (α = 0.61, panel b)
# virulence, with the stochastic model of the authors' code (DoddingtonStochasticProblem).
# Also tabulates the most frequent outcome in each cell.
#
# Invocation (about 25 min with 500 runs per cell on 4 threads):
#   julia --threads=4 --project scripts/doddington2013_fig3.jl [runs] [maturation]
# runs defaults to 500 (the paper's number); maturation is as_coded (default) or as_printed.
# Writes nothing; prints the two survival maps (rows: ρ from high to low, columns: ν).
using Chytridiomycosis, Random, Printf

runs = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 500
maturation = length(ARGS) >= 2 ? Symbol(ARGS[2]) : :as_coded
g = doddington_grid()
for set in ("fig3a", "fig3b")
    S = zeros(10, 10)
    W = fill(:none, 10, 10)
    cells = [(i, j) for i in 1:10, j in 1:10]
    Threads.@threads for k in eachindex(cells)
        i, j = cells[k]   # i: ν, j: ρ
        p = DoddingtonParameters(:doddington2013; set, ν = g.ν[i], ρ = g.ρ[j])
        o = doddington_outcomes(p; runs, maturation, rng = Xoshiro(1000i + j))
        S[j, i] = o.host_survival
        W[j, i] = argmax(s -> getfield(o, s), (:host_extinction, :Bd_extinction, :coexistence))
    end
    println("\n", set, " ($maturation, $runs runs): host survival probability")
    @printf("%10s", "ρ \\ ν")
    foreach(v -> @printf("%8.2e", v), g.ν); println()
    for j in 10:-1:1
        @printf("%10.2e", g.ρ[j]); foreach(i -> @printf("%8.2f", S[j, i]), 1:10); println()
    end
    println(set, ": most frequent outcome (h host extinction, b Bd extinction, c coexistence)")
    for j in 10:-1:1
        @printf("%10.2e  ", g.ρ[j])
        println(join((W[j, i] === :host_extinction ? "h" : W[j, i] === :Bd_extinction ? "b" : "c")
                     for i in 1:10))
    end
end
