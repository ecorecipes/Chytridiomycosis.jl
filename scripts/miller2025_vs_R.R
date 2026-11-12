# R driver for scripts/miller2025_vs_R.jl: simulate the Miller et al. (2025) refugia CTMC with
# the authors' own code (Zenodo 10.5281/zenodo.15179344, v1.0, CC-BY-4.0) at fixed parameter
# sets, and write the observed states and the authors' ABC summary statistic for each run.
#
# The authors' files are sourced read-only from <codedir> (lab/papers/code/miller2025/code):
#   fns_stochastic_model.R  propensities, stoichiometry, run_gillespie, get_timepoints, x0
#   model_setup.R           sample_times and priors (unchanged)
#   calc_summary_stats.R    only the definition of calc_summ_stat is evaluated
# The loop body below is that of simulate_samples.R (lines 33-49): one shaded and one
# unshaded run per parameter set, each with run_gillespie(parms, Nsims = 1), then the
# summary statistic of calc_summary_stats.R (lines 50-59) against exp_data.csv.
#
# Usage: Rscript miller2025_vs_R.R <codedir> <outdir> [nsims] [ncores]
# Needs GillespieSSA (0.6.2), dplyr, tidyr, assertthat. Installs nothing.

suppressPackageStartupMessages({
  library(GillespieSSA); library(dplyr); library(tidyr); library(assertthat); library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
codedir <- normalizePath(args[1])
outdir  <- normalizePath(args[2])
nsims   <- if (length(args) >= 3) as.integer(args[3]) else 2000L
ncores  <- if (length(args) >= 4) as.integer(args[4]) else 4L

owd <- setwd(codedir)
source("fns_stochastic_model.R")
.cargs <- commandArgs
commandArgs <- function(...) character(0)    # model_setup.R takes no arguments
source("model_setup.R")
commandArgs <- .cargs
for (e in parse("calc_summary_stats.R")) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      identical(e[[2]], as.name("calc_summ_stat"))) eval(e)
}
stopifnot(exists("calc_summ_stat"))
expdata <- read.csv("../data/exp_data.csv")
expdata <- pivot_wider(expdata, values_fill = 0, names_from = "compartment", values_from = "N")
setwd(owd)

# Parameter sets (per week). paper: Table 3 posterior means with mu = 0.021 (Eq. 3).
# synthetic: Table S2 = lsfitparams. high: a point with heavy reinfection and loss.
sets <- list(
  paper     = c(beta_sh = 0.675, beta_un = 0.401, alpha = 0.072, omega = 0.569, mu = 0.021),
  synthetic = c(lsfitparams),
  high      = c(beta_sh = 1.5, beta_un = 1.0, alpha = 0.3, omega = 2.0, mu = 0.1)
)

RNGkind("L'Ecuyer-CMRG")
set.seed(20251112)

one_run <- function(iter, prior_i) {
  parms_i <- parameters                                     # simulate_samples.R l. 35-49
  parms_i[names(prior_i)] <- prior_i
  parms_i["beta"] = parms_i["beta_sh"]
  res_shaded = run_gillespie(parms_i, Nsims = 1)[[1]]
  parms_i["beta"] = parms_i["beta_un"]
  res_unshaded = run_gillespie(parms_i, Nsims = 1)[[1]]
  sim_shaded = get_timepoints(res_shaded, sample_times)
  sim_unshaded = get_timepoints(res_unshaded, sample_times)
  sim_i = rbind(data.frame(sim_shaded, shaded = T, iter_id = iter),
                data.frame(sim_unshaded, shaded = F, iter_id = iter))
  error_shaded = calc_summ_stat(obs = expdata[expdata$shaded, ],      # calc_summary_stats.R
                                sim_i = sim_i[sim_i$shaded, ])
  error_unshaded = calc_summ_stat(obs = expdata[!expdata$shaded, ],
                                  sim_i = sim_i[!sim_i$shaded, ])
  list(sim = sim_i,
       err = data.frame(iter_id = iter, err_shaded = error_shaded[["sum_error"]],
                        err_unshaded = error_unshaded[["sum_error"]],
                        steps_shaded = res_shaded$stats$nSteps,
                        steps_unshaded = res_unshaded$stats$nSteps))
}

for (nm in names(sets)) {
  prior_i <- as.data.frame(as.list(sets[[nm]]))
  t0 <- Sys.time()
  res <- mclapply(seq_len(nsims), function(i) one_run(i, prior_i), mc.cores = ncores,
                  mc.set.seed = TRUE)
  cat(nm, ":", nsims, "pairs of runs in", format(Sys.time() - t0), "\n")
  sims <- bind_rows(lapply(res, `[[`, "sim"))
  errs <- bind_rows(lapply(res, `[[`, "err"))
  write.csv(sims, file.path(outdir, paste0("R_states_", nm, ".csv")), row.names = FALSE,
            quote = FALSE)
  write.csv(errs, file.path(outdir, paste0("R_errors_", nm, ".csv")), row.names = FALSE,
            quote = FALSE)
}
writeLines(c(paste("GillespieSSA", packageVersion("GillespieSSA")), R.version.string,
             paste("nsims", nsims)), file.path(outdir, "R_session.txt"))
