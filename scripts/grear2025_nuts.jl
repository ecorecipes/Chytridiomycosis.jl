# Full-posterior fit of the Grear et al. (2025) hurdle model (doi:10.1002/ecm.70001) with NUTS,
# for comparison with the paper's posterior summaries (Figs 6–7) once they are available
# (validation target grear2025_coefficients, pending). Species (phylogenetic, if the authors'
# tree is found), site and region random intercepts; about 1,300 parameters and 10,900
# swabs. Expect hours. Not run as part of the package checks.
#
# Needs an environment with this package (dev'ed), AdvancedHMC and AbstractMCMC:
#   julia --project=<env> -e 'using Pkg; Pkg.develop(path = "."); Pkg.add(["AdvancedHMC", "AbstractMCMC"])'
#   julia --project=<env> scripts/grear2025_nuts.jl [outdir]
#
# Differences from the authors' fit (src/models/grear2025.jl): swabs with missing Age/Sex and
# "Taricha sp." are dropped instead of multiply imputed; one tree, as in their scripts.
using Chytridiomycosis, Random, DelimitedFiles, Statistics
using AdvancedHMC: NUTS
using AbstractMCMC: AbstractMCMC, sample

outdir = isempty(ARGS) ? mktempdir() : mkpath(ARGS[1])
tree = grear2025_phylogeny_path() === nothing ? nothing : grear2025_phylogeny()
prob = Grear2025Problem(; random = (:species, :site, :region), phylogeny = tree)
start = grear2025_fit(prob; iters = 2000)          # posterior mode as the initial point
model = AbstractMCMC.LogDensityModel(prob)
nchains, ndraws, nadapt = 4, 1000, 1000
draws = map(1:nchains) do c
    ch = sample(Xoshiro(2025 + c), model, NUTS(0.95), ndraws + nadapt; n_adapts = nadapt,
                initial_params = start.x .+ 0.01 .* randn(Xoshiro(c), length(start.x)),
                progress = false)
    reduce(hcat, [t.z.θ for t in ch[nadapt+1:end]])'
end
X = reduce(vcat, draws)
coef = [hurdle_coefficients(prob.model, X[i, :]) for i in axes(X, 1)]
names = coef[1].covariates
rows = Any[["part" "covariate" "mean" "q2.5" "q97.5"]]
for (part, f) in (("intensity", c -> c.b), ("prevalence", c -> c.b_prevalence))
    for (j, nm) in enumerate(names)
        v = [f(c)[j] for c in coef]
        push!(rows, [part nm mean(v) quantile(v, 0.025) quantile(v, 0.975)])
    end
end
writedlm(joinpath(outdir, "grear2025_fixed_effects.csv"), reduce(vcat, rows), ',')
println("wrote ", joinpath(outdir, "grear2025_fixed_effects.csv"))
