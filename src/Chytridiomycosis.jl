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
[`MultiHostPoolProblem`](@ref) (with its inoculation callback) and
[`ZoosporePersistenceProblem`](@ref). Use it to solve with other
options or to pass the model to other SciML tools.
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

end # module Chytridiomycosis
