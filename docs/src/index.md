# Chytridiomycosis.jl

```@meta
CurrentModule = Chytridiomycosis
```

```@docs
Chytridiomycosis
```

Models of amphibian chytridiomycosis, caused by the fungi *Batrachochytrium dendrobatidis*
(Bd) and *B. salamandrivorans* (Bsal), built on SciML. The package contains reference
implementations of published models, a parameter compendium that records where every
published value comes from, modified versions of published models, and new models,
including the treatment of amphibians with protective skin bacteria.

Every model belongs to one of three categories. They differ in what a result from them means:

1. **Published models** reimplement a paper's equations. The only departures are documented
   choices where the paper is silent, and switches between what was printed and what the
   authors coded (`as_printed`, `as_coded`). Their results are checked against the papers
   by 352 [validation targets](api/validation.md). Some are not reproduced, and
   those are recorded as discrepant.
2. **Modified published models** add structure the paper did not have: a term, a stage, a
   stochastic version, a geometry. **Their results are not published results.**
3. **New models** are new structure, built on top of published or modified host models.
   **Their results are not published results**, and most of their added parameter values are
   illustrative.


## 1. Published models

| System | Source | Entry points |
|---|---|---|
| Bsal S–L–I in fire salamanders, with diffusive spread | Schmidt et al. 2017, [doi:10.1098/rsos.160801](https://doi.org/10.1098/rsos.160801) | [`SchmidtSLIProblem`](@ref), [`SchmidtPDEProblem`](@ref) |
| Bd loads in frogs; seasonal and stage-structured stochastic versions | Briggs et al. 2010, [doi:10.1073/pnas.0912886107](https://doi.org/10.1073/pnas.0912886107) | [`BriggsLoadProblem`](@ref), [`BriggsSeasonalProblem`](@ref), [`BriggsStagedProblem`](@ref) |
| Bd in common toads at a pond | Mitchell et al. 2008, [doi:10.1098/rspb.2007.1356](https://doi.org/10.1098/rspb.2007.1356) | [`MitchellPondProblem`](@ref) |
| Bd growth in culture | Woodhams et al. 2008, [doi:10.1890/06-1842.1](https://doi.org/10.1890/06-1842.1) | [`intrinsic_growth_rate`](@ref) |
| Several host species sharing one zoospore pool | Daversa et al. 2022, [doi:10.1086/720638](https://doi.org/10.1086/720638) | [`MultiHostPoolProblem`](@ref) |
| Live and dead zoospores in pond water | Rumschlag et al. 2022, [doi:10.1111/1365-2656.13612](https://doi.org/10.1111/1365-2656.13612) | [`ZoosporePersistenceProblem`](@ref), [`RumschlagModel`](@ref) |
| Sunlight-heated refugia (CTMC, ABC) | Miller et al. 2025, [arXiv:2503.06846](https://arxiv.org/abs/2503.06846) | [`RefugiaCTMCProblem`](@ref), [`RefugiaODEProblem`](@ref), [`refugia_abc`](@ref) |
| Seasonal Bsal and a Floquet R0 | Bozzuto & Canessa 2019, [doi:10.1016/j.gecco.2019.e00551](https://doi.org/10.1016/j.gecco.2019.e00551) | [`SeasonalBsalProblem`](@ref) |
| Bsal transmission functions in eastern newts | Malagon et al. 2020, [doi:10.1038/s41598-020-62351-x](https://doi.org/10.1038/s41598-020-62351-x); Tompros et al. 2022, [doi:10.1111/tbed.14043](https://doi.org/10.1111/tbed.14043) | [`TransmissionSLIProblem`](@ref), [`SusceptibleChainProblem`](@ref), [`CohabitationProblem`](@ref) |
| Stochastic individual-based Bd model with interventions | Drawert et al. 2017, [doi:10.1098/rsif.2017.0480](https://doi.org/10.1098/rsif.2017.0480) | [`DrawertIBMProblem`](@ref), [`DrawertODEProblem`](@ref) |
| Bsal with life stages and temperature in eastern newts | Chowdhury et al. 2025, [doi:10.1016/j.ecolmodel.2025.111028](https://doi.org/10.1016/j.ecolmodel.2025.111028) | [`BsalStageTemperatureProblem`](@ref) |
| Integral projection model of Bd load | Wilber et al. 2016, [doi:10.1111/2041-210X.12561](https://doi.org/10.1111/2041-210X.12561) | [`BdIPMProblem`](@ref) |
| Hybrid IPM with a zoospore pool and tadpoles | Wilber et al. 2017, [doi:10.1111/ele.12814](https://doi.org/10.1111/ele.12814) | [`HybridIPMProblem`](@ref) |
| Bsal IPM and management decisions | Canessa et al. 2018, [doi:10.1111/1365-2664.13089](https://doi.org/10.1111/1365-2664.13089) | [`CanessaIPMProblem`](@ref), [`CanessaSIProblem`](@ref) |
| Load distributions and extinction risk | Louca et al. 2014, [doi:10.1098/rspb.2013.2783](https://doi.org/10.1098/rspb.2013.2783) | [`LoucaDOIProblem`](@ref), [`louca_monte_carlo`](@ref) |
| Outbreak-state hidden Markov model | Wilber et al. 2022, [doi:10.1111/1365-2656.13823](https://doi.org/10.1111/1365-2656.13823) | [`OutbreakHMMProblem`](@ref), [`Wilber2022DeclineProblem`](@ref) |
| Multistate capture–recapture survival | Hollanders et al. 2023, [doi:10.1002/eap.2724](https://doi.org/10.1002/eap.2724); Palomar et al. 2023, [doi:10.1007/s10531-022-02525-3](https://doi.org/10.1007/s10531-022-02525-3) | [`CTMultistateProblem`](@ref), [`MultiEventProblem`](@ref) |
| Aggregation of swab loads (Taylor's power law) | Schrock et al. 2025, [doi:10.1098/rspb.2024.2013](https://doi.org/10.1098/rspb.2024.2013) | [`TaylorPowerLawProblem`](@ref) |
| Host–pathogen dynamic occupancy of boreal toads and Bd (equations from the first author's dissertation; data not yet read) | Mosher et al. 2018, [doi:10.1002/eap.1699](https://doi.org/10.1002/eap.1699) | [`Mosher2018Problem`](@ref), [`Mosher2018Model`](@ref) |
| Bd occurrence in host populations by single-season occupancy (likelihood and Appendix A only) | Adams et al. 2010, [doi:10.1890/08-2319.1](https://doi.org/10.1890/08-2319.1) | [`occupancy_loglikelihood`](@ref), [`finite_population_occurrence`](@ref) |

## 2. Modified published models

| Model | Modifies | Entry points |
|---|---|---|
| Radial spread, and scattered colonies simulated in a finite habitat | Schmidt et al. 2017 | `SchmidtPDEProblem(...; geometry = :radial)`, [`SchmidtColonyProblem`](@ref) |
| Live/dead zoospore pool with shedding and uptake, for any host model | Rumschlag et al. 2022 | [`ZoosporePool`](@ref) |
| Stochastic load on one host (birth–death–immigration) | Briggs et al. 2010 | [`LoadBDI`](@ref) |
| Load tipping: saturating clearance, host differences, recruitment; stochastic version | Briggs et al. 2010 | [`LoadTippingProblem`](@ref), [`LoadTippingRecruitmentProblem`](@ref), [`StochasticFrogProblem`](@ref) |
| Bsal with an infectable larval stage | Schmidt et al. 2017 | [`BsalLarvalStage`](@ref), [`BsalLarvalProblem`](@ref) |
| Toads with infection after metamorphosis | Mitchell et al. 2008 | [`PostMetamorphicInfection`](@ref), [`ToadStageProblem`](@ref) |
| Frog load equation fitted to field swabs (data: Knapp et al. 2022) | Briggs et al. 2010 | [`load_knapp2022`](@ref), [`KnappLoadModel`](@ref), [`maximum_likelihood`](@ref), [`profile_likelihood`](@ref), [`KnappPosterior`](@ref) |
| Moments of loads in a pond | Briggs et al. 2010 | [`BriggsMomentProblem`](@ref) |

## 3. New models

| Model | Built on | Entry points |
|---|---|---|
| Microbiome (probiotic) treatment | Schmidt, Mitchell, load tipping, Bsal with larvae | [`MicrobiomeTreatment`](@ref), [`ProbioticSLIProblem`](@ref), [`ProbioticPondProblem`](@ref) (also stochastic, `solve(prob; rng)`), [`ProbioticLoadProblem`](@ref), [`ProbioticLarvalProblem`](@ref) |
| Bacterial abundance on each frog, with repeated baths | Load tipping | [`ProbioticAbundanceProblem`](@ref) |
| One nested model with each piece of biology as a switch, and its stochastic version | Reduces to Schmidt et al. 2017 | [`ChytridBiology`](@ref), [`NestedChytridProblem`](@ref), [`StochasticNestedProblem`](@ref) |
| Swab and eDNA observation model, and identifiability of the Briggs parameters | Briggs moments | [`SurveyModel`](@ref), [`simulate_survey`](@ref), [`identifiable_combinations`](@ref) |

## Shared infrastructure

These are not models of any one paper: the parameter compendium, units of loads and rates,
thermal performance curves, validation targets, generic integral projection model, hidden
Markov model and multistate capture–recapture machinery, qPCR load units and likelihoods,
and glue to DynamicalResilience.jl.

## Getting started

```julia
using Chytridiomycosis
p = SchmidtParameters(:schmidt2017)      # values from the compendium
basic_reproduction_number(p)             # 24.8
sol = solve(SchmidtSLIProblem(p))        # an ODE solution
```

Every problem is a struct describing the model and its set-up, solved with `solve`. Most
return a named tuple of time series and summaries; see each problem's documentation.

## Tutorials

Every tutorial opens with a table of the models it uses and the category of each.

**Part I. Published models**

1. [Can a Bsal outbreak be stopped?](tutorials/01_bsal.md): five published Bsal models, from
   transmission and host density to seasons, life stages, temperature, load and management.
2. [Does it matter that frogs carry different Bd loads?](tutorials/02_load_ipm.md): from the
   mean load to load distributions and integral projection models fitted to swabs.
3. [Does chance decide whether frogs survive Bd?](tutorials/03_stochastic.md): whole frogs and
   whole zoospores, in a frog population and in outdoor mesocosms.
4. [Where does Bd live between hosts?](tutorials/04_environment.md): zoospores in water,
   several host species sharing one pool, and temperature.
5. [What do field data say about Bd?](tutorials/05_field_data.md): swab loads, outbreak states
   of lakes, and survival of infected frogs with imperfect detection.
6. [Thresholds and resilience](tutorials/06_resilience.md): how robust a state is on the way
   to a threshold.

**Part II. Modified and new models**

7. [Does a protective skin microbiome save the population?](tutorials/07_building_realism.md):
   for three systems, the treatment model is built from the simplest model one piece of
   biology at a time, ending with chance.
8. [One model, many switches](tutorials/08_nested_model.md): a single nested model in which
   every combination of biology is compared like for like, and a stochastic version.
9. [Fitting to a field trial](tutorials/09_knapp_fit.md): the frog load equation fitted to
   field swabs of treated, wild and control frogs, comparing four explanations for why
   treated frogs stayed low.

The tutorials are Quarto documents in `vignettes/`, rendered to HTML, PDF and Markdown.

## Licence

MIT. Data files in `data/` keep the licences given in their READMEs.
