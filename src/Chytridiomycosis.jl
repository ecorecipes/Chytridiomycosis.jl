"""
    Chytridiomycosis

Reference implementations of models of amphibian chytridiomycosis (*Batrachochytrium
dendrobatidis*, Bd, and *B. salamandrivorans*, Bsal), with a parameter compendium that records
units and provenance. Plans and verification notes: the OKF wiki in `../../wiki`.
"""
module Chytridiomycosis

using CommonSolve: CommonSolve, solve
using LinearAlgebra: eigvals
using OrdinaryDiffEqTsit5: Tsit5
using SciMLBase: SciMLBase, ODEProblem, VectorContinuousCallback
using Distributions: Beta, Binomial, DiscreteNonParametric, LogNormal, NegativeBinomial, Normal,
    Poisson
using Random: Random, AbstractRNG
using TOML

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

export solve
export ParameterRecord, compendium_sources, parameter_records, parameter_set
export SchmidtParameters, SchmidtSLIProblem, to_ode_problem, invasion_threshold,
       removal_fraction, basic_reproduction_number, endemic_equilibrium, hopf_threshold,
       wave_speed, scattered_colony_area, population_crash
export BriggsParameters, BriggsLoadProblem, transmission_rate, load_growth_rate,
       host_density_threshold
export BriggsStage, BriggsSeasonalParameters, BriggsSeasonalProblem, persistence_probability
export MitchellParameters, MitchellPondProblem
export LoadBDI, aggregation, load_distribution, mean_load, prevalence, infection_gain_rate,
       infection_loss_rate, swab, enzootic_pool, enzootic_load_distribution
export LoadTipping, LoadTippingProblem, pool_feedback, net_growth_rate, critical_immigration,
       tipping_regime, load_equilibria, critical_host_number, mean_time_to_runaway,
       low_state_distribution, pooled_load_equilibria, LoadTippingRecruitmentProblem
export BsalLarvalStage, BsalLarvalProblem, larval_reproduction_number, disease_free_larvae
export PostMetamorphicInfection, ToadStageProblem
export ProbioticLarvalProblem, ProbioticAbundanceProblem
export MicrobiomeTreatment, load_reducing_treatment, protection, protection_duration, reversion_rate, microbiome_reproduction_number,
       equilibrium_coverage, treatment_efficacy, mortality_reduction, effective_reproduction_number,
       critical_coverage, ProbioticSLIProblem, ProbioticPondProblem, ProbioticLoadProblem

end # module Chytridiomycosis
