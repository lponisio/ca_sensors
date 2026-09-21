# --------------- #
#  Load Packages  #
# --------------- #
library(tidyverse)
library(nimble)
library(coda)
library(jagsUI)
library(patchwork)
library(vegan)
library(tidyr)
library(lubridate)
library(BPAbook)
library(purrr)

# ----------------------- #
#  Set working directory  #
# ----------------------- #

setwd("~/")
## run the lab_paths script to appropriately set path for your computer (must be set up for each lab member)
## once amended with your computer's information, lab_paths.R can be saved to your home directory
source("lab_paths.R")
local.path
## generate path from your folder to the CA sensors folder
dir.CASensors <- file.path(local.path,"ca_sensors_saved")
## change working directory
setwd(dir.CASensors)

## need to download the BPAbook package from the following URL to install it
## Needed for a variety of convenience functions and for toy datasets to build models
# https://www.vogelwarte.ch/en/research/population-biology/book-bpa/
#install.packages("../ca_sensors/src/BPAbook_0.0.1.tar.gz", repos = NULL, type = "source")
#library(BPAbook)

# load up basic data simulator functions
source("../ca_sensors/0.1_data_sim.R")
# load the m1 model fitting functions
source("../ca_sensors/scr_nimbleFit_m1.R")

# -----------------------------------------------------------------------
# 1. ONE TRIAL, WITH PSI DIAGNOSTICS CARRIED THROUGH AS COLUMNS
# -----------------------------------------------------------------------
# Same as run_one_trial() in scr_simulator_daily.R, except it also grabs the
# psi_mean / psi_flag_M_too_small attributes fit_scr0_nimble() attaches --
# grabbed immediately after fit_fun(sim) runs, before any dplyr verb has a
# chance to drop them, then added as ordinary columns (which do survive).
run_one_trial_diag <- function(sim_params, fit_fun) {
  sim <- do.call(simulate_scr_dataset, sim_params)
  fit <- fit_fun(sim)
  
  psi_mean <- attr(fit, "psi_mean")
  psi_flag <- attr(fit, "psi_flag_M_too_small")
  if (is.null(psi_mean)) psi_mean <- NA_real_   # fit_fun didn't attach one -- fine, just NA
  if (is.null(psi_flag)) psi_flag <- NA
  
  truth <- tibble(
    param = c("N", "sigma", "p0"),
    truth = c(sim$true$N, sim$true$sigma, sim$true$p0)
  )
  
  fit %>%
    left_join(truth, by = "param") %>%
    mutate(psi_mean = psi_mean, psi_flag_M_too_small = psi_flag)
}

# -----------------------------------------------------------------------
# 2. ONE SCENARIO: many trials under one fixed set of sim_params
# -----------------------------------------------------------------------

#' @param coverage_target, coverage_tol Flag a param as "reliably recovered"
#'   only if coverage falls within coverage_target +/- coverage_tol (default:
#'   within 0.85-1.0 for a nominal 95% interval).
#' @param bias_tol_pct Also require |bias| / |mean_truth| below this (default
#'   10%) to call it reliably recovered -- catches cases where an interval
#'   is wide enough to "cover" truth despite a systematically off estimate.
run_scenario <- function(sim_params, fit_fun, n_sims = 20, seed_start = 1,
                         coverage_target = 0.95, coverage_tol = 0.10,
                         bias_tol_pct = 0.10) {
  
  results <- map_dfr(1:n_sims, function(i) {
    set.seed(seed_start + i)
    run_one_trial_diag(sim_params, fit_fun) %>% mutate(sim_id = i)
  })
  
  summary <- results %>%
    group_by(param) %>%
    summarise(
      n_sims             = n(),
      mean_truth         = mean(truth),
      mean_estimate      = mean(estimate, na.rm = TRUE),
      bias               = mean(estimate - truth, na.rm = TRUE),
      pct_bias           = bias / mean_truth,
      rmse               = sqrt(mean((estimate - truth)^2, na.rm = TRUE)),
      coverage           = mean(truth >= lower & truth <= upper, na.rm = TRUE),
      convergence_rate   = mean(converged, na.rm = TRUE),
      psi_flag_rate      = mean(psi_flag_M_too_small, na.rm = TRUE),  # fraction of trials where M was too small
      .groups = "drop"
    ) %>%
    mutate(
      reliably_recovered = abs(coverage - coverage_target) <= coverage_tol &
        abs(pct_bias) <= bias_tol_pct
    )
  
  list(results = results, summary = summary)
}

# -----------------------------------------------------------------------
# 3. FULL SWEEP: run many scenarios, each varying whatever params you like
# -----------------------------------------------------------------------

#' @param base_sim_params A complete sim_params list (your usual defaults).
#' @param scenarios A NAMED list of override lists. Each element is merged
#'   onto base_sim_params via modifyList() -- so you only need to specify
#'   what changes for that scenario. Works for ANY simulate_scr_dataset()
#'   parameter: N, sigma, p0, n_weeks, spacing, n_cols, n_rows,
#'   buffer_multiplier, netting_n_cells, detector_type, etc.
#' @param fit_fun Either:
#'     (a) a plain function(sim) -> tibble, used unchanged for every
#'         scenario (fine if N doesn't vary much across scenarios), or
#'     (b) a "factory" function(merged_sim_params) -> function(sim), if you
#'         want the fit itself to adapt per scenario -- e.g. scaling the
#'         data-augmentation ceiling M with that scenario's N so a sweep
#'         across very different N values doesn't need one M for all of
#'         them. Set is_factory = TRUE to use this mode.
#' @param is_factory See fit_fun above.
run_parameter_sweep <- function(base_sim_params, scenarios, fit_fun,
                                is_factory = FALSE,
                                n_sims = 20, seed_start = 1,
                                coverage_target = 0.95, coverage_tol = 0.10,
                                bias_tol_pct = 0.10) {
  
  stopifnot(!is.null(names(scenarios)), all(names(scenarios) != ""))
  
  scenario_output <- imap(scenarios, function(overrides, scenario_name) {
    
    merged_params <- modifyList(base_sim_params, overrides)
    this_fit_fun <- if (is_factory) fit_fun(merged_params) else fit_fun
    
    out <- run_scenario(merged_params, this_fit_fun, n_sims, seed_start,
                        coverage_target, coverage_tol, bias_tol_pct)
    
    # tag every row with the scenario name and whatever scalar params were
    # actually varied in this scenario (for easy plotting/filtering later).
    # Non-scalar overrides (data frames, NULL, etc. -- e.g. sampling_df,
    # camera_coords) are skipped here since they don't tabulate as a column,
    # but they still drove the simulation above via merged_params.
    scalar_overrides <- keep(overrides, ~ is.atomic(.x) && length(.x) == 1)
    tag_cols <- as_tibble(scalar_overrides)
    
    list(
      results = bind_cols(scenario = scenario_name, out$results, tag_cols[rep(1, nrow(out$results)), , drop = FALSE]),
      summary = bind_cols(scenario = scenario_name, out$summary, tag_cols[rep(1, nrow(out$summary)), , drop = FALSE])
    )
  })
  
  list(
    results = map_dfr(scenario_output, "results"),
    summary = map_dfr(scenario_output, "summary")
  )
  }

# =============================================================================
# EXAMPLE USAGE (not run automatically)
# =============================================================================


##### --------------------------------- #######
#####  TESTING SIMULATOR-FIT FUNCTIONS  #######
##### --------------------------------- #######

### Sampling effort ###
# for now use actual days that we sampled for the number of sampling occasions
effort <- read.csv("./data/cleaned/CASensors_Effort_clean.csv",
                   header = T) %>%
  select(date, grid_cell)

# Generate list of parameters to run simulator
sim_params <- list(
  N = 400, sigma = 120, # meters
  p0 = 0.15, # baseline detection probability
  n_weeks = 10, # study length (number of primary periods)
  camera_coords = NULL, # leave NULL if simulating a grid
  n_cols = 5, n_rows = 4, spacing = 130,
  start_col_letter = "F", start_row = 3,
  sampling_df = effort, # netting effort data
  netting_mode = "random",
  study_start_date = as.Date("2026-05-25"), # calendar date of day 1 of week 1
  camera_uptime = 1, # placeholder (real netting effort now comes from sampling_df)
  detector_type = "netting"
)

# define fit function
fit_fun_scr0 <- make_scr0_fit_fun(M = 800, niter = 10000, nburnin = 3000,
                                  thin = 5, nchains = 2)

# time_one_scenario.R -- run this BEFORE the full sweep

source("scr_simulator_daily.R")
source("scr0_nimble_fit.R")
source("scr_parameter_sweep.R")

# use your real baseline params, but just ONE trial to start
base_sim_params <- list(
  N = 400, sigma = 120, p0 = 0.15, n_weeks = 10,
  camera_coords = NULL, n_cols = 5, n_rows = 4, spacing = 130,
  camera_uptime = 1, detector_type = "netting",
  study_start_date = as.Date("2026-05-25"),
  netting_mode = "random", sampling_df = effort
)

fit_fun_scr0 <- make_scr0_fit_fun(
  effort_matrix_name = "netting",
  M = 800, niter = 5000, nburnin = 1000, thin = 3, nchains = 3
)

# time just ONE trial
### something's not quite right, and it's throwing an error
### my money is on needing to fix the effort building structure. I need to make it so it's easier to custom build effort over time
timing <- system.time({
  one_result <- run_scenario(base_sim_params, fit_fun_scr0, n_sims = 1)
})

print(timing)   # look at "elapsed" (seconds)

# rough extrapolation to the full sweep
seconds_per_trial <- timing["elapsed"]
n_scenarios <- 10   # however many you're planning
n_sims_per_scenario <- 20

est_total_minutes <- (seconds_per_trial * n_scenarios * n_sims_per_scenario) / 60
cat("Estimated total sweep time:", round(est_total_minutes, 1), "minutes\n")

base_sim_params <- list(
  N = 200, sigma = 100, p0 = 0.15, n_weeks = 10,
  camera_coords = NULL, n_cols = 5, n_rows = 4, spacing = 130,
  camera_uptime = 1, detector_type = "netting",
  netting_mode = "random", netting_n_cells = 12, netting_days_per_week = 3
)

# vary ONE thing at a time, in any combination of parameters you like
scenarios <- list(
  baseline        = list(),
  N_low           = list(N = 50),
  N_high          = list(N = 600),
  sigma_small     = list(sigma = 40),
  sigma_large     = list(sigma = 200),
  sparse_grid     = list(spacing = 250),
  dense_grid      = list(spacing = 70),
  short_season    = list(n_weeks = 4),
  low_p0          = list(p0 = 0.05),
  small_buffer    = list(buffer_multiplier = 2)
)

# fit_fun as a FACTORY so M scales with each scenario's true N
# (important here since N ranges from 50 to 600 across scenarios)
fit_fun_factory <- function(params) {
  make_scr0_fit_fun(
    effort_matrix_name = "netting",
    M = ceiling(2.5 * params$N),
    niter = 20000, nburnin = 5000, thin = 5, nchains = 3
  )
}

sweep <- run_parameter_sweep(
  base_sim_params, scenarios,
  fit_fun = fit_fun_factory, is_factory = TRUE,
  n_sims = 20
)

print(sweep$summary)
# sweep$summary$reliably_recovered tells you, per scenario x parameter,
# whether coverage and bias both landed within tolerance.
# sweep$summary$psi_flag_rate flags scenarios where M was probably too
# small too often -- bump M for those and re-run just those scenarios.