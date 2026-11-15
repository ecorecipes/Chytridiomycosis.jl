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
                                  "Fitting to a field trial" => "tutorials/03_knapp_fit.md",
                                  "Chance and persistence" => "tutorials/04_stochastic.md",
                                  "Can a Bsal outbreak be stopped?" => "tutorials/05_bsal.md",
                                  "Load-structured models" => "tutorials/06_load_ipm.md"],
                  "API reference" => ["Parameter compendium" => "api/compendium.md",
                                      "Published models" => "api/published.md",
                                      "Several hosts, one pool" => "api/daversa2022.md",
                                      "Zoospore persistence" => "api/rumschlag2022.md",
                                      "Load models" => "api/load.md",
                                      "Microbiome treatment" => "api/treatment.md",
                                      "Life stages" => "api/stages.md",
                                      "Nested and stochastic models" => "api/nested.md",
                                      "Fitting to field data" => "api/fitting.md",
                                      "Units and loads" => "api/units.md",
                                      "Thermal performance" => "api/thermal.md",
                                      "Validation targets" => "api/validation.md",
                                      "Swabs, eDNA and identifiability" => "api/identifiability.md",
                                      "Bayesian fit to field data" => "api/knapp_bayes.md",
                                      "Refugia CTMC (Miller 2025)" => "api/miller2025.md",
                                      "Seasonal Bsal (Bozzuto 2019)" => "api/bozzuto2019.md",
                                      "Bsal transmission functions" => "api/bsal_transmission.md",
                                      "Stochastic IBM (Drawert 2017)" => "api/drawert2017.md",
                                      "Life stages and temperature (Chowdhury 2025)" => "api/chowdhury2025.md",
                                      "Integral projection models (Wilber 2016)" => "api/wilber2016.md",
                                      "Extinction risk with a zoospore pool (Wilber 2017)" => "api/wilber2017.md",
                                      "Management decisions with an IPM (Canessa 2018)" => "api/canessa2018.md",
                                      "Load distributions and extinction risk (Louca 2014)" => "api/louca2014.md",
                                      "Swab loads: units, likelihoods, operators" => "api/observation_load.md",
                                      "Outbreak-state HMM (Wilber 2022)" => "api/wilber2022.md",
                                      "Multistate survival (Hollanders 2023, Palomar 2023)" => "api/multistate.md"]])

if get(ENV, "CI", "false") == "true"
    deploydocs(; repo = "github.com/ecorecipes/Chytridiomycosis.jl.git", push_preview = true)
end
