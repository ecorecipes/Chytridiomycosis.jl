# Chytridiomycosis.jl

Models of amphibian chytridiomycosis (*Batrachochytrium dendrobatidis*, Bd, and
*B. salamandrivorans*, Bsal) in Julia, built on SciML. Work in progress.

## What is here

Every model belongs to one of three categories. They differ in what a result from them means:

1. **Published models** reimplement a paper's equations. The only departures are documented
   choices where the paper is silent, and switches between what was printed and what the
   authors coded (`as_printed`, `as_coded`). Their results are checked against the papers
   by 352 validation targets (`data/validation/targets.toml`). Some are not reproduced, and
   those are recorded as discrepant.
2. **Modified published models** add structure the paper did not have: a term, a stage, a
   stochastic version, a geometry. **Their results are not published results.**
3. **New models** are new structure, built on top of published or modified host models.
   **Their results are not published results**, and most of their added parameter values are
   illustrative.

### 1. Published models

| System | Source | Entry points |
|---|---|---|
| Bsal S–L–I in fire salamanders, with diffusive spread | Schmidt et al. 2017, doi:10.1098/rsos.160801 | `SchmidtSLIProblem`, `SchmidtPDEProblem` |
| Bd loads in frogs; seasonal and stage-structured stochastic versions | Briggs et al. 2010, doi:10.1073/pnas.0912886107 | `BriggsLoadProblem`, `BriggsSeasonalProblem`, `BriggsStagedProblem` |
| Bd in common toads at a pond | Mitchell et al. 2008, doi:10.1098/rspb.2007.1356 | `MitchellPondProblem` |
| Bd growth in culture | Woodhams et al. 2008, doi:10.1890/06-1842.1 | `intrinsic_growth_rate` |
| Several host species sharing one zoospore pool | Daversa et al. 2022, doi:10.1086/720638 | `MultiHostPoolProblem` |
| Live and dead zoospores in pond water | Rumschlag et al. 2022, doi:10.1111/1365-2656.13612 | `ZoosporePersistenceProblem`, `RumschlagModel` |
| Sunlight-heated refugia (CTMC, ABC) | Miller et al. 2025, arXiv:2503.06846 | `RefugiaCTMCProblem`, `RefugiaODEProblem`, `refugia_abc` |
| Seasonal Bsal and a Floquet R0 | Bozzuto & Canessa 2019, doi:10.1016/j.gecco.2019.e00551 | `SeasonalBsalProblem` |
| Bsal transmission functions in eastern newts | Malagon et al. 2020, doi:10.1038/s41598-020-62351-x; Tompros et al. 2022, doi:10.1111/tbed.14043 | `TransmissionSLIProblem`, `SusceptibleChainProblem`, `CohabitationProblem` |
| Stochastic individual-based Bd model with interventions | Drawert et al. 2017, doi:10.1098/rsif.2017.0480 | `DrawertIBMProblem`, `DrawertODEProblem` |
| Bsal with life stages and temperature in eastern newts | Chowdhury et al. 2025, doi:10.1016/j.ecolmodel.2025.111028 | `BsalStageTemperatureProblem` |
| Integral projection model of Bd load | Wilber et al. 2016, doi:10.1111/2041-210X.12561 | `BdIPMProblem` |
| Hybrid IPM with a zoospore pool and tadpoles | Wilber et al. 2017, doi:10.1111/ele.12814 | `HybridIPMProblem` |
| Bsal IPM and management decisions | Canessa et al. 2018, doi:10.1111/1365-2664.13089 | `CanessaIPMProblem`, `CanessaSIProblem` |
| Load distributions and extinction risk | Louca et al. 2014, doi:10.1098/rspb.2013.2783 | `LoucaDOIProblem`, `louca_monte_carlo` |
| Outbreak-state hidden Markov model | Wilber et al. 2022, doi:10.1111/1365-2656.13823 | `OutbreakHMMProblem`, `Wilber2022DeclineProblem` |
| Multistate capture–recapture survival | Hollanders et al. 2023, doi:10.1002/eap.2724; Palomar et al. 2023, doi:10.1007/s10531-022-02525-3 | `CTMultistateProblem`, `MultiEventProblem` |
| Aggregation of swab loads (Taylor's power law) | Schrock et al. 2025, doi:10.1098/rspb.2024.2013 | `TaylorPowerLawProblem` |
| Adult survival before and after Bd from age structure; spatial metapopulation viability of *Litoria* (Bd only lowers survival) | Heard et al. 2024, doi:10.1111/1365-2656.14177 | `HeardCatchCurve`, `HeardMetapopulationProblem` |
| Environmental-water timing and extinction risk of *Litoria raniformis* (age-classified PVA) | Turner et al. 2025, doi:10.1111/csp2.70022 | `TurnerPVAProblem` |

### 2. Modified published models

| Model | Modifies | Entry points |
|---|---|---|
| Radial spread, and scattered colonies simulated in a finite habitat | Schmidt et al. 2017 | `SchmidtPDEProblem(...; geometry = :radial)`, `SchmidtColonyProblem` |
| Live/dead zoospore pool with shedding and uptake, for any host model | Rumschlag et al. 2022 | `ZoosporePool` |
| Stochastic load on one host (birth–death–immigration) | Briggs et al. 2010 | `LoadBDI` |
| Load tipping: saturating clearance, host differences, recruitment; stochastic version | Briggs et al. 2010 | `LoadTippingProblem`, `LoadTippingRecruitmentProblem`, `StochasticFrogProblem` |
| Bsal with an infectable larval stage | Schmidt et al. 2017 | `BsalLarvalStage`, `BsalLarvalProblem` |
| Toads with infection after metamorphosis | Mitchell et al. 2008 | `PostMetamorphicInfection`, `ToadStageProblem` |
| Frog load equation fitted to field swabs (data: Knapp et al. 2022) | Briggs et al. 2010 | `load_knapp2022`, `KnappLoadModel`, `maximum_likelihood`, `profile_likelihood`, `KnappPosterior` |
| Moments of loads in a pond | Briggs et al. 2010 | `BriggsMomentProblem` |

### 3. New models

| Model | Built on | Entry points |
|---|---|---|
| Microbiome (probiotic) treatment | Schmidt, Mitchell, load tipping, Bsal with larvae | `MicrobiomeTreatment`, `ProbioticSLIProblem`, `ProbioticPondProblem` (also stochastic, `solve(prob; rng)`), `ProbioticLoadProblem`, `ProbioticLarvalProblem` |
| Bacterial abundance on each frog, with repeated baths | Load tipping | `ProbioticAbundanceProblem` |
| One nested model with each piece of biology as a switch, and its stochastic version | Reduces to Schmidt et al. 2017 | `ChytridBiology`, `NestedChytridProblem`, `StochasticNestedProblem` |
| Swab and eDNA observation model, and identifiability of the Briggs parameters | Briggs moments | `SurveyModel`, `simulate_survey`, `identifiable_combinations` |

### Shared infrastructure

These are not models of any one paper: the parameter compendium, units of loads and rates,
thermal performance curves, validation targets, generic integral projection model, hidden
Markov model and multistate capture–recapture machinery, qPCR load units and likelihoods,
and glue to DynamicalResilience.jl.

## Vignettes

Every vignette opens with a table of the models it uses and the category of each.

### Part I. Published models

[Can a Bsal outbreak be stopped?](vignettes/01_bsal/01_bsal.md)
([PDF](vignettes/01_bsal/01_bsal.pdf),
[HTML](vignettes/01_bsal/01_bsal.html),
[source](vignettes/01_bsal/01_bsal.qmd))

Climbs five published Bsal models (Malagon et al. 2020 and Tompros et al. 2022, Schmidt et al.
2017, Bozzuto & Canessa 2019, Chowdhury et al. 2025, Canessa et al. 2018), adding transmission,
demography, seasons, life stages and temperature, then load and management, and asks at each
step whether reducing host numbers stops an outbreak.

[Does it matter that frogs carry different Bd loads?](vignettes/02_load_ipm/02_load_ipm.md)
([PDF](vignettes/02_load_ipm/02_load_ipm.pdf),
[HTML](vignettes/02_load_ipm/02_load_ipm.html),
[source](vignettes/02_load_ipm/02_load_ipm.qmd))

From the mean load on identical frogs (Briggs et al. 2010) to load distributions (Louca et al.
2014) and integral projection models fitted to swabs (Wilber et al. 2016, 2017), reporting at
each step the invasion number, the hosts left and the spread of load.

[Does chance decide whether frogs survive Bd?](vignettes/03_stochastic/03_stochastic.md)
([PDF](vignettes/03_stochastic/03_stochastic.pdf),
[HTML](vignettes/03_stochastic/03_stochastic.html),
[source](vignettes/03_stochastic/03_stochastic.qmd))

Whole frogs and whole zoospores: one outbreak in a frog population (Briggs et al. 2010, Drawert
et al. 2017), then recovery and reinfection in outdoor mesocosms (Miller et al. 2025).

[Where does Bd live between hosts?](vignettes/04_environment/04_environment.md)
([PDF](vignettes/04_environment/04_environment.pdf),
[HTML](vignettes/04_environment/04_environment.html),
[source](vignettes/04_environment/04_environment.qmd))

Zoospores in water (Rumschlag et al. 2022), several host species sharing one pool (Daversa et
al. 2022), and temperature setting the pathogen's growth (Woodhams et al. 2008 and thermal
performance curves).

[What do field data say about Bd?](vignettes/05_field_data/05_field_data.md)
([PDF](vignettes/05_field_data/05_field_data.pdf),
[HTML](vignettes/05_field_data/05_field_data.html),
[source](vignettes/05_field_data/05_field_data.qmd))

From one qPCR swab to the spread of loads across swabs (Schrock et al. 2025), outbreak states of
lakes from yearly surveys (Wilber et al. 2022), and survival of infected frogs with imperfect
detection (Hollanders et al. 2023, Palomar et al. 2023).

[Thresholds and resilience](vignettes/06_resilience/06_resilience.md)
([PDF](vignettes/06_resilience/06_resilience.pdf),
[HTML](vignettes/06_resilience/06_resilience.html),
[source](vignettes/06_resilience/06_resilience.qmd))

How robust a state is on the way to a threshold: the Hopf threshold of the Schmidt et al. (2017)
Bsal model and the fold of the load tipping model (a modified Briggs model).

### Part II. Modified and new models

[Does a protective skin microbiome save the population?](vignettes/07_building_realism/07_building_realism.md)
([PDF](vignettes/07_building_realism/07_building_realism.pdf),
[HTML](vignettes/07_building_realism/07_building_realism.html),
[source](vignettes/07_building_realism/07_building_realism.qmd))

For three systems (Bsal in fire salamanders, Bd in common toads, Bd in frogs) it starts from the
simplest model, adds the biology of the host and pathogen one piece at a time, then the biology
of the treatment, and records after every step how the population fares with and without
treatment. It ends by measuring how far each step moved the answer and, with a factorial design,
whether the effect of one piece depends on the others.

[One model, many switches](vignettes/08_nested_model/08_nested_model.md)
([PDF](vignettes/08_nested_model/08_nested_model.pdf),
[HTML](vignettes/08_nested_model/08_nested_model.html),
[source](vignettes/08_nested_model/08_nested_model.qmd))

Uses a single nested model in which each piece of biology is a switch, runs all 48 combinations
at three host densities, and measures the effect of each piece like for like. It ends with a
stochastic version, which overturns one of the deterministic results.

[Fitting the frog load model to a field trial](vignettes/09_knapp_fit/09_knapp_fit.md)
([PDF](vignettes/09_knapp_fit/09_knapp_fit.pdf),
[HTML](vignettes/09_knapp_fit/09_knapp_fit.html),
[source](vignettes/09_knapp_fit/09_knapp_fit.qmd))

Fits the frog load equation to field swabs from frogs treated with itraconazole and then
*Janthinobacterium lividum* (Knapp et al. 2022), building up from no treatment effect to four
explanations. Something held treated frogs' loads down for weeks. A waning drug and the applied
bacteria fit about equally well, but resident *J. lividum* at the measured abundances did not
protect.

Each vignette writes out the equations of every model it runs, with a symbol table and a
diagram of the model.

Vignettes are Quarto documents with executable Julia code, rendered to HTML, PDF and
Markdown: run `quarto render` in `vignettes/`. Drawing the diagrams needs
[Graphviz](https://graphviz.org) (the `dot` program) on the `PATH`, and the PDFs need
`rsvg-convert` (librsvg) to convert them. Render one vignette and one format at a time, since
Quarto's Julia engine shares one notebook server.

## Parameters

Published values are read from a TOML compendium (`data/compendium/`, one file per source
paper) with units, provenance and a `verified` flag saying whether the value was checked
against the primary paper. Nothing published is hard-coded in the models.

```julia
using Chytridiomycosis
p = SchmidtParameters(:schmidt2017)
basic_reproduction_number(p)          # 24.8
sol = solve(SchmidtSLIProblem(p))
```

## Tests

```julia
using Pkg; Pkg.test()
```

Requires Julia 1.12.

## Documentation

[ecorecipes.github.io/Chytridiomycosis.jl](https://ecorecipes.github.io/Chytridiomycosis.jl)
(built by CI from `docs/`; the tutorials are the rendered vignettes, copied with
`julia scripts/sync_vignettes.jl`).

## Licence

MIT (see `LICENSE`). Data files in `data/` keep the licences given in their READMEs.
