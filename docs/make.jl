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
                  "Tutorials" => ["Published models" => ["Can a Bsal outbreak be stopped?" => "tutorials/01_bsal.md",
                                                         "Does it matter that frogs carry different loads?" => "tutorials/02_load_ipm.md",
                                                         "Does chance decide whether frogs survive?" => "tutorials/03_stochastic.md",
                                                         "Where does Bd live between hosts?" => "tutorials/04_environment.md",
                                                         "What do field data say about Bd?" => "tutorials/05_field_data.md",
                                                         "Thresholds and resilience" => "tutorials/06_resilience.md"],
                                  "Modified and new models" => ["Does a protective microbiome save the population?" => "tutorials/07_building_realism.md",
                                                                "One model, many switches" => "tutorials/08_nested_model.md",
                                                                "Fitting to a field trial" => "tutorials/09_knapp_fit.md"]],
                  "Published models" => ["Schmidt 2017, Briggs 2010, Mitchell 2008" => "api/published.md",
                                         "Several hosts, one pool (Daversa 2022)" => "api/daversa2022.md",
                                         "Zoospore persistence (Rumschlag 2022)" => "api/rumschlag2022.md",
                                         "Adaptive immunity and treatments (Woodhams 2011)" => "api/woodhams2011.md",
                                         "Atelopus load and zoospore output (DiRenzo 2014)" => "api/direnzo2014.md",
                                         "Refugia CTMC (Miller 2025)" => "api/miller2025.md",
                                         "Seasonal Bsal (Bozzuto 2019)" => "api/bozzuto2019.md",
                                         "Bsal transmission functions" => "api/bsal_transmission.md",
                                         "Stochastic IBM (Drawert 2017)" => "api/drawert2017.md",
                                         "Imperfect prophylaxis ABM (Barnett 2026)" => "api/barnett2026.md",
                                         "Life stages and temperature (Chowdhury 2025)" => "api/chowdhury2025.md",
                                         "Integral projection models (Wilber 2016)" => "api/wilber2016.md",
                                         "Extinction risk with a zoospore pool (Wilber 2017)" => "api/wilber2017.md",
                                         "Management decisions with an IPM (Canessa 2018)" => "api/canessa2018.md",
                                         "Load distributions and extinction risk (Louca 2014)" => "api/louca2014.md",
                                         "Outbreak-state HMM (Wilber 2022)" => "api/wilber2022.md",
                                         "Multistate survival (Hollanders 2023, Palomar 2023)" => "api/multistate.md"],
                  "Modified models" => ["Load models" => "api/load.md",
                                        "Life stages" => "api/stages.md",
                                        "Fitting to field data" => "api/fitting.md",
                                        "Bayesian fit to field data" => "api/knapp_bayes.md"],
                  "New models" => ["Microbiome treatment" => "api/treatment.md",
                                   "Nested and stochastic models" => "api/nested.md",
                                   "Swabs, eDNA and identifiability" => "api/identifiability.md"],
                  "Infrastructure" => ["Parameter compendium" => "api/compendium.md",
                                       "Units and loads" => "api/units.md",
                                       "Thermal performance" => "api/thermal.md",
                                       "Swab loads: units, likelihoods, operators" => "api/observation_load.md",
                                       "Validation targets" => "api/validation.md",
                                       "Resilience and thresholds" => "api/resilience.md"]])

if get(ENV, "CI", "false") == "true"
    deploydocs(; repo = "github.com/ecorecipes/Chytridiomycosis.jl.git", push_preview = true)
end
