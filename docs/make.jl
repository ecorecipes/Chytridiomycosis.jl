using Documenter
using Chytridiomycosis

makedocs(;
         modules = [Chytridiomycosis],
         authors = "Simon Frost",
         sitename = "Chytridiomycosis.jl",
         checkdocs = :exports,
         warnonly = false,
         format = Documenter.HTML(;
                                  prettyurls = get(ENV, "CI", "false") == "true",
                                  canonical = "https://ecorecipes.github.io/Chytridiomycosis.jl",
                                  size_threshold = 2^21,
                                  # the tutorials carry equations and SVG figures inline
                                  size_threshold_warn = 2^18),
         pages = ["Home" => "index.md",
                  "Tutorials" => ["Building biological realism" => "tutorials/01_building_realism.md",
                                  "One model, many switches" => "tutorials/02_nested_model.md",
                                  "Fitting to a field trial" => "tutorials/03_knapp_fit.md"],
                  "API reference" => ["Parameter compendium" => "api/compendium.md",
                                      "Published models" => "api/published.md",
                                      "Several hosts, one pool" => "api/daversa2022.md",
                                      "Load models" => "api/load.md",
                                      "Microbiome treatment" => "api/treatment.md",
                                      "Life stages" => "api/stages.md",
                                      "Nested and stochastic models" => "api/nested.md",
                                      "Fitting to field data" => "api/fitting.md",
                                      "Units and loads" => "api/units.md",
                                      "Thermal performance" => "api/thermal.md",
                                      "Validation targets" => "api/validation.md",
                                      "Swabs, eDNA and identifiability" => "api/identifiability.md",
                                      "Bayesian fit to field data" => "api/knapp_bayes.md"]])

if get(ENV, "CI", "false") == "true"
    deploydocs(; repo = "github.com/ecorecipes/Chytridiomycosis.jl.git", push_preview = true)
end
