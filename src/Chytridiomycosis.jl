"""
    Chytridiomycosis

Reference implementations of models of amphibian chytridiomycosis (*Batrachochytrium
dendrobatidis*, Bd, and *B. salamandrivorans*, Bsal), with a parameter compendium that records
units and provenance. Plans and verification notes: the OKF wiki in `../../wiki`.
"""
module Chytridiomycosis

using CommonSolve: CommonSolve, solve
using LinearAlgebra: eigvals, diag
using OrdinaryDiffEqTsit5: Tsit5
using SciMLBase: SciMLBase, ODEProblem, VectorContinuousCallback, DiscreteCallback
using Distributions: Beta, Binomial, DiscreteNonParametric, LogNormal, NegativeBinomial, Normal,
    Poisson, logcdf, logpdf
using Random: Random, AbstractRNG
using Dates: Date, @dateformat_str
using LogDensityProblems: LogDensityProblems
using TOML

"""
    to_ode_problem(prob)

The `SciMLBase.ODEProblem` that `solve(prob)` integrates, for the problems that lower to a
single ODE: [`SchmidtSLIProblem`](@ref), [`BriggsLoadProblem`](@ref),
[`LoadTippingProblem`](@ref), [`BsalLarvalProblem`](@ref),
[`MultiHostPoolProblem`](@ref) (with its inoculation callback),
[`ZoosporePersistenceProblem`](@ref), [`RefugiaODEProblem`](@ref),
[`SeasonalBsalProblem`](@ref) (with its inflow switch times as `tstops`),
[`TransmissionSLIProblem`](@ref), [`SusceptibleChainProblem`](@ref),
[`CohabitationProblem`](@ref), [`DrawertODEProblem`](@ref),
[`BsalStageTemperatureProblem`](@ref) and [`LoucaDOIProblem`](@ref) (method of lines over
DOI bins, with exposure, disturbance and end-of-run callbacks). Use it to solve with
other options or to pass the model to other SciML tools.
"""
function to_ode_problem end

include("compendium.jl")
include("models/schmidt2017.jl")
include("models/briggs2010.jl")
include("models/briggs2010_seasonal.jl")
include("models/mitchell2008.jl")
include("novel/load_bdi.jl")
include("novel/load_tipping.jl")
include("novel/bsal_larval.jl")
include("novel/toad_postmetamorphic.jl")
include("novel/microbiome_treatment.jl")
include("novel/microbiome_treatment_stages.jl")
include("novel/microbiome_abundance.jl")
include("novel/frog_stochastic.jl")
include("novel/nested.jl")
include("novel/nested_stochastic.jl")
include("fitting/knapp2022.jl")
# Schmidt spatial, Briggs stages
include("models/schmidt2017_spatial.jl")
include("models/briggs2010_stages.jl")

# Phase 0 foundations
include("models/woodhams2008.jl")
include("units.jl")
include("thermal.jl")
include("validation.jl")
# Daversa et al. 2022: several host species sharing one zoospore pool
include("models/daversa2022.jl")

# Rumschlag et al. 2022 (needs the units layer)
include("models/zoospore_pool.jl")
include("models/rumschlag2022.jl")
include("fitting/rumschlag2022.jl")

# N5
using LinearAlgebra: svd, eigen, Symmetric, Diagonal
using Distributions: pdf
include("novel/swab_edna.jl")

# Miller et al. 2025 refugia CTMC
using Distributions: median, quantile
include("models/miller2025.jl")
include("fitting/miller2025_abc.jl")

# Bozzuto & Canessa 2019 seasonal Bsal model with a Floquet R0
include("models/bozzuto2019.jl")
# Bsal transmission functions: Malagon et al. 2020 and Tompros et al. 2022
include("models/bsal_transmission.jl")
include("fitting/malagon2020.jl")
include("fitting/tompros2022.jl")
# Drawert et al. 2017 stochastic individual-based Bd model with one-time interventions
include("models/drawert2017.jl")
# Chowdhury et al. 2025 stage- and temperature-structured Bsal model of the eastern newt
include("models/chowdhury2025.jl")
# Integral projection models: generic machinery, and Wilber et al. 2016 Bd–R. muscosa IPM
using LinearAlgebra: LinearAlgebra, dot, mul!, logdet
using Distributions: UnivariateDistribution, cdf, ccdf, MvNormal
using SciMLBase: DiscreteProblem
include("ipm.jl")
include("models/wilber2016.jl")
include("fitting/wilber2016.jl")
# Wilber et al. 2017 hybrid IPM: zoospore pool, tadpoles, seasonality, demographic stochasticity
include("models/wilber2017.jl")
include("fitting/wilber2017.jl")
# Canessa et al. 2018 Bsal integral projection model and management actions
include("models/canessa2018.jl")
# Louca et al. 2014 DOI-distribution model of adults and tadpoles, Q_o and extinction risk
include("models/louca2014.jl")

export solve
export ParameterRecord, compendium_sources, parameter_records, parameter_set
export SchmidtParameters, SchmidtSLIProblem, to_ode_problem, invasion_threshold,
       removal_fraction, basic_reproduction_number, endemic_equilibrium, hopf_threshold,
       wave_speed, scattered_colony_area, population_crash
export BriggsParameters, BriggsLoadProblem, transmission_rate, load_growth_rate,
       host_density_threshold
export BriggsStage, BriggsSeasonalParameters, BriggsSeasonalProblem, persistence_probability
export MitchellParameters, MitchellPondProblem
export ZoosporePool, pool_derivatives, pool_state, qpcr_zoospores, infectious_zoospores,
       viable_fraction
export RumschlagParameters, ZoosporePersistenceProblem, zoospore_lifespan, persistence_time,
       infection_probability, RumschlagCount, RumschlagTadpole, rumschlag2022_path,
       load_rumschlag2022, RumschlagModel, rumschlag_parameters, rumschlag_model_selection,
       rumschlag_covariance, rumschlag_summary
export LoadBDI, aggregation, load_distribution, mean_load, prevalence, infection_gain_rate,
       infection_loss_rate, swab, enzootic_pool, enzootic_load_distribution
export LoadTipping, LoadTippingProblem, pool_feedback, net_growth_rate, critical_immigration,
       tipping_regime, load_equilibria, critical_host_number, mean_time_to_runaway,
       low_state_distribution, pooled_load_equilibria, LoadTippingRecruitmentProblem
export BsalLarvalStage, BsalLarvalProblem, larval_reproduction_number, disease_free_larvae
export PostMetamorphicInfection, ToadStageProblem
export ProbioticLarvalProblem, ProbioticAbundanceProblem, StochasticFrogProblem
export ChytridBiology, NestedChytridProblem, invasion_growth_rate, StochasticNestedProblem
export KnappSwab, knapp2022_path, load_knapp2022, KnappLoadModel, parameter_names, predicted_loads,
    knapp_start, maximum_likelihood, profile_likelihood
export MicrobiomeTreatment, load_reducing_treatment, protection, protection_duration, reversion_rate, microbiome_reproduction_number,
       equilibrium_coverage, treatment_efficacy, mortality_reduction, effective_reproduction_number,
       critical_coverage, ProbioticSLIProblem, ProbioticPondProblem, ProbioticLoadProblem
# Schmidt spatial, Briggs stages
export SchmidtPDEProblem, SchmidtColonyProblem, front_position, front_speed
export BriggsLifeStages, BriggsStagedProblem, disease_free_stages

# Phase 0 foundations
export PathogenLoad, ZEPerSwab, ITS1CopiesPerSwab, Zoospores, SporangiaPerHost, load_value,
       convert_load, parameter_load, convert_rate, in_time_unit, survival_to_hazard,
       hazard_to_survival, time_unit
export ThermalPerformanceCurve, Briere, CardinalTPC, GaussianQuadratic, SharpeSchoolfield,
       Ratkowsky, PolynomialTPC, InterpolatedTPC, thermal_optimum, IsolateTPC, HostTPC,
       thermal_data, SeasonalTemperature, temperature_dependent
export WoodhamsParameters, intrinsic_growth_rate
export DaversaParameters, MultiHostPoolProblem, daversa_doses, species_reproduction_numbers,
       next_generation_matrix, daversa2022_path, MesocosmTank, load_daversa2022,
       observed_infected, predicted_infected, daversa_loglikelihood, tank_reproduction_numbers
export ValidationTarget, validation_targets, validation_target, reproduce, check_target

# N5
export BriggsMomentProblem, SurveyPond, SurveySample, SurveyObservation, SurveyModel,
       simulate_survey, sensitivity_rank, expected_information, observed_information,
       identifiable_combinations

# Bayesian Knapp fit
export KnappPosterior, knapp_prior, observed_swabs, pointwise_loglikelihood, psis_loo, waic,
       elpd_difference, turing_model

# Miller et al. 2025 refugia CTMC
export RefugiaParameters, RefugiaCTMCProblem, RefugiaODEProblem, transmissibility,
       refugia_initial_state, refugia_rates, refugia_stoichiometry, refugia_observation_weeks,
       RefugiaMesocosmData, miller2025_path, load_miller2025, loss_rate, RefugiaPrior,
       refugia_distance, refugia_abc, refugia_abc_summary

# Bozzuto & Canessa 2019 seasonal Bsal model
export BozzutoParameters, SeasonalBsalProblem, contact_scaling, contact_rate,
       transmission_probability, bsal_mortality, larval_inflow, disease_free_susceptibles,
       instantaneous_reproduction_number, averaged_reproduction_number, seasonal_factor,
       floquet_exponent, floquet_reproduction_number, target_proportion, floquet_elasticities,
       floquet_multipliers, periodic_reproduction_number
# Bsal transmission functions (Malagon et al. 2020; Tompros et al. 2022)
export TransmissionFunction, DensityDependent, FrequencyDependent, PowerLaw,
       QuadraticSaturation, StructuralSaturation, PrevalenceDependent, PrevalencePowerLaw,
       PrevalenceSaturation, HollingContact, infection_rate, force_of_infection,
       transmission_coefficient, TransmissionSLIProblem, SusceptibleChainProblem,
       susceptible_distribution, expected_susceptibles, MalagonParameters, contact_rate,
       CohabitationProblem, cohabitation_curve, final_size_reproduction_number,
       TomprosParameters, malagon2020_path, load_malagon2020, contact_means, contact_anova,
       fit_contact_function, fit_cohabitation, malagon_table_s1, tompros2022_path,
       TomprosTank, load_tompros2022, load_tompros_swabs, clean_infection_series,
       tompros_status, tompros_counts, TomprosModel, transmission_function,
       tompros_model_selection, likelihood_interval, tompros_negloglik
# Drawert et al. 2017 stochastic individual-based Bd model
export DrawertParameters, DrawertIBMProblem, DrawertODEProblem, DrawertIntervention,
       drawert_interventions, drawert_variant, model_variant, encystment_probability,
       sporangium_loss_rate, nondimensional_parameters, nondimensional_ranges,
       nondimensional_distance, equilibrium_load_ratio, drawert2017_path, load_drawert2017,
       drawert_summary, outcome_class, survival_group, intervention_effect
# Chowdhury et al. 2025 stage- and temperature-structured Bsal model
export ChowdhuryParameters, BsalStageTemperatureProblem, chowdhury_scenarios,
       stage_transmission_rates, stage_totals, stage_decline_time, stage_disease_free_state,
       stage_next_generation_matrix, stage_reproduction_numbers, stage_elasticities
# Integral projection models (generic) and Wilber et al. 2016
export IPMMesh, meshpoints, meshwidth, meshedges, IPMKernel, DiscretizedIPM, discretize,
       projection_matrix, population_growth_rate, stable_distribution, sensitivity_matrix,
       elasticity_matrix, absorption_times, eviction, mesh_convergence, load_summary
export WilberParameters, LoadTransmission, BdIPMProblem, BdIPMSolution, to_discrete_problem,
       stable_load_distribution, mean_time_to_death, vital_rate_elasticities, wilber2016_path,
       load_wilber2016, WilberFit, fit_wilber2016, WilberRegression, growth_rate_bootstrap
# IPM additions (environmental reservoir, demographic stochasticity) and Wilber et al. 2017
export reservoir_reproduction_number, demographic_projection, demographic_projection!
export Wilber2017Parameters, zoospore_survival, annual_demography_matrix, HybridIPMProblem,
       HybridIPMSolution, hybrid_summary, wilber2017_temperature, extinction_curve, stasis_load,
       extinction_sensitivity, sensitivity_weights, wilber2017_path, load_wilber2017,
       Wilber2017Transmission, Wilber2017Regression, posterior_summary, posterior_draws
# Canessa et al. 2018 Bsal integral projection model and management actions
export CanessaDraw, canessa2018_path, load_canessa2018, CanessaAction, canessa_actions,
       canessa_action, canessa_growth_coefficients, CanessaParameters, canessa_survival,
       canessa_transmission, canessa_initial_load, canessa_growth_mean, canessa_growth_variance,
       CanessaIPMProblem, canessa_density_ratio, canessa_cohort_survival, canessa_dispersal,
       canessa_ensemble, canessa_summary, CanessaSIParameters, CanessaSIProblem,
       required_reduction, final_size_transmission, final_density_ratio
# Louca et al. 2014 DOI-distribution model and Monte Carlo extinction risk
export LoucaParameters, LoucaDOIProblem, tadpoles_tolerate, louca_breeding_shape,
       louca_breeding_rate, louca_transmission, louca_fecundity_factor,
       louca_metamorphosis_success, louca_base_mortality, louca_disease_mortality,
       louca_load_rates, louca_load_speed, louca_within_host_capacity, louca_within_host_loss,
       louca_competition, louca_turnover, reproduction_quotient, louca_jacobian,
       louca_growth_rate, louca_tolerance_from_slope, louca_survival_curve,
       louca_transmission_from_growth, louca_endemic_transmission, louca_outcome,
       louca_ranges, louca_sample, louca_monte_carlo, louca_risk_effects,
       louca_outcome_probabilities, louca2014_path, load_louca2014

end # module Chytridiomycosis
