module ChytridiomycosisTuringExt

using Chytridiomycosis: Chytridiomycosis, KnappPosterior
using LogDensityProblems: LogDensityProblems
using Turing: Turing, @model, MvNormal
using LinearAlgebra: Diagonal

# log parameters x with the posterior's normal prior; the swab log-likelihood added on top
@model function knapp_turing(post, μ, Σ)
    x ~ MvNormal(μ, Σ)
    Turing.@addlogprob! LogDensityProblems.logdensity(post.model, x)
end

Chytridiomycosis.turing_model(post::KnappPosterior) =
    knapp_turing(post, post.μ, Diagonal(post.s .^ 2))

end
