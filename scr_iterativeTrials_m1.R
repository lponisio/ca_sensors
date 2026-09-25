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
# 1. COMPILE ONCE
# -----------------------------------------------------------------------

#' Build and compile the SCR0 model + MCMC a single time. Reuse the result
#' across every trial in a scenario via fit_scr0_nimble_compiled().
build_scr0_compiled <- function(constants, data, inits) {
  model <- nimbleModel(code = scr0_code, constants = constants, data = data, inits = inits)
  cmodel <- compileNimble(model)
  mcmc_conf <- configureMCMC(model, monitors = c("N", "p0", "sigma", "psi"))
  mcmc <- buildMCMC(mcmc_conf)
  cmcmc <- compileNimble(mcmc, project = model)
  list(model = cmodel, mcmc = cmcmc)
}

# -----------------------------------------------------------------------
# 2. FIT ONE TRIAL AGAINST AN ALREADY-COMPILED MODEL
# -----------------------------------------------------------------------

#' @param compiled Output of build_scr0_compiled(). Its M/J/K/state-space
#'   must match what this sim + M would produce -- caller's responsibility.
#' @param M Fixed data-augmentation ceiling -- must NOT be NULL here.
fit_scr0_nimble_compiled <- function(sim, compiled, effort_matrix, M,
                                     niter = 5000, nburnin = 1000, thin = 5, nchains = 3,
                                     rhat_threshold = 1.1) {
  
  stopifnot(!is.null(M)) # needs a value of M to continue
  prep <- prepare_scr0_nimble_data(sim, effort_matrix, M) # prep data for nimble
  
  # swap in this trial's data on the already-compiled model (no recompile)
  compiled$model$y <- prep$data$y
  
  inits_list <- replicate(nchains, prep$inits(), simplify = FALSE)
  
  # runMCMC() is the documented entry point for reusing one compiled MCMC
  # across repeated runs with fresh data/inits -- this is the part most
  # worth double-checking in the self-check below.
  samples <- runMCMC(
    compiled$mcmc,
    niter = niter, nburnin = nburnin, thin = thin, nchains = nchains,
    inits = inits_list, setSeed = FALSE, samplesAsCodaMCMC = TRUE
  )
  
  if (nchains > 1) {
    rhat <- coda::gelman.diag(samples, multivariate = FALSE)$psrf[, "Point est."]
    converged <- all(rhat[c("N", "p0", "sigma")] < rhat_threshold, na.rm = TRUE)
    pooled <- as.matrix(samples)
  } else {
    converged <- NA
    pooled <- as.matrix(samples)
  }
  
  summarize_param <- function(param) {
    x <- pooled[, param]
    tibble(
      param = param,
      estimate = mean(x),
      lower = quantile(x, 0.025, names = FALSE),
      upper = quantile(x, 0.975, names = FALSE),
      converged = converged
    )
  }
  
  bind_rows(summarize_param("N"), summarize_param("p0"), summarize_param("sigma"))
}
# -----------------------------------------------------------------------
# 3. ONE SCENARIO: compile once, run n_sims trials (parallel on Unix/Mac)
# -----------------------------------------------------------------------

#' @param M Fixed data-augmentation ceiling for this scenario. Required.
#' @param n_cores Number of cores to use. >1 only actually parallelizes on
#'   Unix/Mac (uses fork via parallel::mclapply, which lets forked processes
#'   inherit the compiled model at no extra cost). On Windows this silently
#'   falls back to running sequentially in a single process -- you still get
#'   the compile-once benefit, just not the parallel one. (A Windows-parallel
#'   path would need a PSOCK cluster that recompiles once per worker instead
#'   of once total -- not implemented here.)
run_scenario_fast <- function(sim_params, effort_matrix_name = "camera", M,
                              niter = 5000, nburnin = 1000, thin = 5, nchains = 3,
                              n_sims = 20, seed_start = 1, n_cores = 1,
                              coverage_target = 0.95, coverage_tol = 0.10,
                              bias_tol_pct = 0.10) {
  
  stopifnot(!is.null(M))
  
  # build + compile once, using a throwaway simulated dataset purely to fix
  # the model's dimensions (M/J/K/state-space) -- not counted as a trial
  set.seed(seed_start)
  template_sim <- do.call(simulate_scr_dataset, sim_params)
  template_prep <- prepare_scr0_nimble_data(
    template_sim, template_sim$data$effort[[effort_matrix_name]], M
  )
  compiled <- build_scr0_compiled(template_prep$constants, template_prep$data, template_prep$inits())
  
  run_one <- function(i) {
    tryCatch({
      set.seed(seed_start + i)
      sim <- do.call(simulate_scr_dataset, sim_params)
      fit <- fit_scr0_nimble_compiled(
        sim, compiled, sim$data$effort[[effort_matrix_name]], M,
        niter, nburnin, thin, nchains
      )
      truth <- tibble(
        param = c("N", "sigma", "p0"),
        truth = c(sim$true$N, sim$true$sigma, sim$true$p0)
      )
      fit %>% left_join(truth, by = "param") %>% mutate(sim_id = i)
    }, error = function(e) {
      warning("Trial ", i, " failed: ", conditionMessage(e))
      NULL
    })
  }
  
  use_fork <- .Platform$OS.type == "unix" && n_cores > 1
  if (n_cores > 1 && !use_fork) {
    message("n_cores > 1 requested but forking isn't available on this OS -- running sequentially.")
  }
  
  trial_list <- if (use_fork) {
    parallel::mclapply(1:n_sims, run_one, mc.cores = n_cores)
  } else {
    lapply(1:n_sims, run_one)
  }
  
  n_failed <- sum(map_lgl(trial_list, is.null))
  if (n_failed > 0) warning(n_failed, " of ", n_sims, " trials failed and were dropped.")
  
  results <- bind_rows(trial_list)
  
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
      .groups = "drop"
    ) %>%
    mutate(
      reliably_recovered = abs(coverage - coverage_target) <= coverage_tol &
        abs(pct_bias) <= bias_tol_pct
    )
  
  list(results = results, summary = summary, n_failed = n_failed)
}

# -----------------------------------------------------------------------
# 4. FULL SWEEP -- same shape as run_parameter_sweep(), calls the fast path
# -----------------------------------------------------------------------

#' @param M_fun function(merged_sim_params) -> fixed M for that scenario.
#'   Default scales M with the scenario's N.
run_parameter_sweep_fast <- function(base_sim_params, scenarios,
                                     effort_matrix_name = "camera",
                                     M_fun = function(p) max(50, ceiling(2.5 * p$N)),
                                     niter = 5000, nburnin = 1000, thin = 5, nchains = 3,
                                     n_sims = 20, seed_start = 1, n_cores = 1,
                                     coverage_target = 0.95, coverage_tol = 0.10,
                                     bias_tol_pct = 0.10) {
  
  stopifnot(!is.null(names(scenarios)), all(names(scenarios) != ""))
  
  scenario_output <- imap(scenarios, function(overrides, scenario_name) {
    
    merged_params <- modifyList(base_sim_params, overrides)
    M <- M_fun(merged_params)
    
    out <- run_scenario_fast(merged_params, effort_matrix_name, M,
                             niter, nburnin, thin, nchains,
                             n_sims, seed_start, n_cores,
                             coverage_target, coverage_tol, bias_tol_pct)
    
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


# ===========================================================
# M1 model run - uncomment to run simulation and save output
# ===========================================================

# # generate the baseline scenario used for simulations
# base_sim_params <- list(
#   N = 300, sigma = 120, p0 = 0.15, n_weeks = 10,
#   camera_coords = NULL, n_cols = 5, n_rows = 4, spacing = 130,
#   camera_uptime = 1, detector_type = "camera",
#   netting_mode = "random", netting_n_cells = 6,
#   netting_days_per_week = 3 # netting effort specs listed, but will only use camera traps effort
# )
# 
# # vary ONE thing at a time, in any combination of parameters you like
# scenarios <- list(
#   baseline     = list(),
#   N_low        = list(N = 50),
#   N_high       = list(N = 600),
#   sigma_small  = list(sigma = 40),
#   sigma_large  = list(sigma = 220),
#   sparse_grid  = list(spacing = 250),
#   dense_grid   = list(spacing = 70),
#   short_season = list(n_weeks = 4),
#   low_p0       = list(p0 = 0.05),
#   small_buffer = list(buffer_multiplier = 2),
#   high_effort  = list(netting_days_per_week = 4),
#   low_effort   = list(netting_days_per_week = 2)
# )
# 
# # run the "full sweep" of alternative parameter scenarios
M1_full <- run_parameter_sweep_fast(base_sim_params,
                                     scenarios = scenarios,
                                     effort_matrix_name = "camera", # using cameras for sampling
                                     niter = 5000, nburnin = 1000, thin = 5, nchains = 3,
                                     n_sims = 20, seed_start = 1, n_cores = 5,
                                     coverage_target = 0.95, coverage_tol = 0.10,
                                     bias_tol_pct = 0.10)

# View the results summary for the simulated model runs
# Reliably recovered indicates whether the proportion of trials where the true value fell within the trial's
#   95% credible interval (coverage) fell within the coverage target +/- the tolerance threshold.
#   So if 85% to 100% of trials had the true parameter value fall within the
#   trial's 95% credible interval, it would receive a pass (reliably_recovered = TRUE).
#View(M1_full$summary)

#saveRDS(M1_full, file = "../ca_sensors/r_outputs/M1_model_simulations.R")
View(M1_full$results)
write.csv(M1_full$summary, "~/Desktop/M1_summary.csv", row.names = F)
write.csv(M1_full$results, "~/Desktop/M1_results.csv", row.names = F)

# -----------------------------------------------------------------------
# 5. PLOTS -- both work directly off results (from run_scenario_fast() or
#    run_parameter_sweep_fast()); no changes needed elsewhere to use these.
# -----------------------------------------------------------------------

#' Per-trial recovery: one point (+ 95% CI) per trial, vs. the true value.
#' Faceted by parameter (rows) x scenario (columns).
plot_recovery <- function(results, params = NULL, scenarios = NULL) {
  df <- results
  if (!is.null(params)) df <- filter(df, param %in% params)
  if (!is.null(scenarios)) df <- filter(df, scenario %in% scenarios)
  
  ggplot(df, aes(x = sim_id, y = estimate)) +
    geom_hline(aes(yintercept = truth), color = "firebrick", linetype = "dashed") +
    geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.3, color = "steelblue") +
    facet_grid(param ~ scenario, scales = "free_y") +
    labs(x = "Trial", y = "Estimate (95% CI)",
         title = "Per-trial parameter recovery",
         subtitle = "Dashed line = true value") +
    theme_minimal()
}

#' Density of estimates across trials, one curve per scenario, faceted by
#' parameter. Dashed vertical lines mark each scenario's true value.
plot_estimate_density <- function(results, params = NULL, scenarios = NULL) {
  df <- results
  if (!is.null(params)) df <- filter(df, param %in% params)
  if (!is.null(scenarios)) df <- filter(df, scenario %in% scenarios)
  
  truth_df <- distinct(df, scenario, param, truth)
  
  ggplot(df, aes(x = estimate, fill = scenario, color = scenario)) +
    geom_density(alpha = 0.3) +
    geom_vline(data = truth_df, aes(xintercept = truth, color = scenario),
               linetype = "dashed", show.legend = FALSE) +
    facet_wrap(~ param, scales = "free") +
    labs(x = "Estimate", y = "Density",
         title = "Distribution of estimates across trials",
         subtitle = "Dashed lines = true values") +
    theme_minimal()
}

# =============================================================================
# PLOTS -- example usage
# =============================================================================

M1_summ_plot <- plot_recovery(M1_full$results)

ggsave(M1_summ_plot,
       file = "./figures/simulation_figures/M1_summary_pointsError.png",
       device = "png", units = "in", height = 6, width = 11)

plot_estimate_density(M1_full$results)
