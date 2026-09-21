# ---------------------- #
#  M1 Nimble Model Code  #
# ---------------------- #

scr0_code <- nimbleCode({
  
  psi ~ dunif(0, 1)
  p0  ~ dunif(0, 1)
  sigma ~ dunif(0, sigma_upper)
  
  for (i in 1:M) {
    z[i] ~ dbern(psi) # data augmentation inclusion
    s[i, 1] ~ dunif(xlim[1], xlim[2]) # activity centers
    s[i, 2] ~ dunif(ylim[1], ylim[2]) # activity centers
    
    for (j in 1:J) {
      d2[i, j] <- (s[i, 1] - trapX[j])^2 + (s[i, 2] - trapY[j])^2 # activity center distances from traps
      p[i, j]  <- p0 * exp(-d2[i, j] / (2 * sigma^2)) * z[i] # Half-normal detection probability * inclusion
      y[i, j]  ~ dbinom(size = K[j], prob = p[i, j]) # Detection probabilities summed over occasions
    }
  }
  
  N <- sum(z[1:M]) # population size estimate
})

# ------------- #
# 2. Data Prep  #
# ------------- #
#turn a `sim` object (from simulate_scr_dataset) into the
#    constants/data/inits NIMBLE needs

#' @param sim A simulated dataset from simulate_scr_dataset().
#' @param effort_matrix Which J x K effort matrix to use for per-trap
#'   occasion counts (K_j). Defaults to sim$data$effort$camera -- override
#'   this to match whichever detector_type you actually simulated detections
#'   with (e.g. sim$data$effort$combined).
#' @param M Data-augmentation upper bound on N. If NULL, defaults to
#'   max(50, 2.5 * n_detected) -- a starting rule of thumb. IMPORTANT: after
#'   fitting, check that the posterior for psi isn't pushed up near 1; if it
#'   is, M was too small and you need to re-fit with a larger M.
prepare_scr0_nimble_data <- function(sim, effort_matrix = sim$data$effort$camera, M = NULL) {
  
  y_obs <- sim$data$y_obs     # n_detected x J x K
  n_detected <- dim(y_obs)[1] # get number of detected bees
  J <- dim(y_obs)[2]          # get number of cells
  
  if (is.null(M)) M <- max(50, ceiling(2.5 * n_detected)) # checks for value of M, if NULL makes a rule of thumb M
  if (M <= n_detected) stop("M must be larger than the number of detected individuals.")
  
  # per-trap occasion counts: how many occasions each trap actually operated
  K_j <- colSums(effort_matrix)
  
  # sum detections across occasions per individual x trap (valid since p is
  # constant across occasions in this baseline model)
  n_matrix <- apply(y_obs, c(1, 2), sum) # n_detected x J
  
  # generate appropriately sized matrix of augmented individuals attached to matrix of real indivs
  y_augmented <- rbind(n_matrix, matrix(0L, nrow = M - n_detected, ncol = J))
  
  traps <- sim$data$traps # location and identity of traps
  ss <- sim$data$state_space # state space boundaries
  
  # generate list of constants
  constants <- list(
    M = M, J = J,
    xlim = ss$xlim, ylim = ss$ylim,
    trapX = traps$x, trapY = traps$y,
    K = K_j,
    sigma_upper = max(diff(ss$xlim), diff(ss$ylim))   # weakly informative upper bound
    )
  data <- list(y = y_augmented)
  
  # sensible starting values: known-detected individuals start near the
  # centroid of the traps where they were caught (helps convergence);
  # augmented all-zero individuals start at random locations
  init_s <- matrix(NA_real_, nrow = M, ncol = 2)
  for (i in seq_len(n_detected)) {
    caught_at <- which(n_matrix[i, ] > 0)
    init_s[i, ] <- c(mean(traps$x[caught_at]), mean(traps$y[caught_at]))
  }
  if (M > n_detected) {
    init_s[(n_detected + 1):M, 1] <- runif(M - n_detected, ss$xlim[1], ss$xlim[2])
    init_s[(n_detected + 1):M, 2] <- runif(M - n_detected, ss$ylim[1], ss$ylim[2])
    }
  init_z <- c(rep(1L, n_detected), rep(0L, M - n_detected))
  
  inits <- function() list(
    psi = n_detected / M,
    p0 = 0.2,
    sigma = mean(diff(ss$xlim), diff(ss$ylim)) / 10,   # rough starting guess
    z = init_z,
    s = init_s
    )
  
  list(code = scr0_code, constants = constants, data = data, inits = inits,
       n_detected = n_detected, M = M)
  }

# -----------------#
# 3. Fit Function  #
# -----------------#

#' Fit the SCR0 model in NIMBLE and return a harness-compatible summary.
#'
#' @param sim Simulated dataset from simulate_scr_dataset()
#' @param effort_matrix, M See prepare_scr0_nimble_data()
#' @param niter, nburnin, thin, nchains MCMC settings. Increase niter/nchains
#'   for real inference; the defaults here are modest so a first test run
#'   doesn't take forever.
#' @param rhat_threshold Chains are flagged "converged" if the worst Rhat
#'   (Gelman-Rubin statistic) across N, p0, sigma is below this.
fit_scr0_nimble <- function(sim, effort_matrix = sim$data$effort$camera, M = NULL,
                            niter = 20000, nburnin = 5000, thin = 5, nchains = 3,
                            rhat_threshold = 1.1) {
  
  prep <- prepare_scr0_nimble_data(sim, effort_matrix, M)
  
  inits_list <- replicate(nchains, prep$inits(), simplify = FALSE)
  
  samples <- nimbleMCMC(
    code = prep$code,
    constants = prep$constants,
    data = prep$data,
    inits = inits_list,
    monitors = c("N", "p0", "sigma", "psi"),
    niter = niter, nburnin = nburnin, thin = thin, nchains = nchains,
    samplesAsCodaMCMC = TRUE
  )
  
  # samples is an mcmc.list when nchains > 1 (needed for Rhat); a single
  # mcmc matrix if nchains == 1 (no Rhat possible then)
  if (nchains > 1) {
    rhat <- coda::gelman.diag(samples, multivariate = FALSE)$psrf[, "Point est."]
    converged <- all(rhat[c("N", "p0", "sigma")] < rhat_threshold, na.rm = TRUE)
    pooled <- as.matrix(samples)   # stack all chains together for summaries
  } else {
    rhat <- setNames(rep(NA_real_, 3), c("N", "p0", "sigma"))
    converged <- NA   # can't assess convergence with a single chain
    pooled <- as.matrix(samples)
  }
  
  # flag if data augmentation ceiling M was probably too low
  psi_mean <- mean(pooled[, "psi"])
  if (psi_mean > 0.9) {
    warning("Posterior mean of psi = ", round(psi_mean, 3),
            " is close to 1 -- M (", prep$M, ") was likely too small. Re-fit with a larger M.")
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

#' Convenience factory: bakes MCMC settings into a plain function(sim), so it
#' drops straight into run_one_trial(sim_params, fit_fun) / run_many_trials()
#' without changing those functions at all.
make_scr0_fit_fun <- function(effort_matrix_name = "camera", M = NULL,
                              niter = 20000, nburnin = 5000, thin = 5, nchains = 2) {
  function(sim) {
    fit_scr0_nimble(
      sim,
      effort_matrix = sim$data$effort[[effort_matrix_name]],
      M = M, niter = niter, nburnin = nburnin, thin = thin, nchains = nchains
    )
  }
}
