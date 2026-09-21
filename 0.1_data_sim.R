# The following functions collectively simulate a simple spatial capture recapture dataset
# Each of the functions in sections 1-4 is executed in the wrapper function found in section 5
# Data are simulated at a daily scale

### TODO: - Modify simulator to change grid spacing
###       - Modify simulator to change number of occasions
###       - Modify simulator to change number of grid cells/day and number of days sampled/week

# ----------------------------- #
#  1. Grid coordinates/geometry #
# ----------------------------- #

#' Generate a rectangular camera grid with alphanumeric names (e.g. "F3"),
#' the same naming convention as the study's original 4x5 grid: columns get
#' letters, rows get numbers, spaced `spacing` meters apart.
#'
#' @param n_cols, n_rows Grid dimensions (default 5 x 4, matching the
#'   original study). These are the two things you'll want to iterate on
#'   later to test how grid layout affects model performance.
#' @param spacing Distance between adjacent cells, in meters.
#' @param start_col_letter, start_row Where the naming starts (defaults
#'   "F" and 3, reproducing the original F3:J6 grid exactly).
build_camera_grid <- function(n_cols = 5, n_rows = 4, spacing = 130,
                              start_col_letter = "F", start_row = 3) {
  start_idx <- match(start_col_letter, LETTERS)
  stopifnot(!is.na(start_idx), start_idx + n_cols - 1 <= 26)
  col_letters <- LETTERS[start_idx:(start_idx + n_cols - 1)]
  row_numbers <- start_row:(start_row + n_rows - 1)
  
  expand.grid(col_letter = col_letters, row_number = row_numbers) %>%
    mutate(
      grid_cell = paste0(col_letter, row_number),
      col_index = match(col_letter, col_letters),
      row_index = match(row_number, row_numbers),
      x = col_index * spacing,
      y = -row_index * spacing
    ) %>%
    select(grid_cell, x, y)
}

#' Build trap (camera) locations and a state-space polygon (rectangle)
#'
#' @param camera_coords Data frame with columns grid_cell, x, y giving REAL
#'   camera coordinates. If NULL, falls back to a generated rectangular grid
#'   via build_camera_grid() using n_cols/n_rows/spacing/start_col_letter/
#'   start_row below -- handy for iterating on grid layout in simulations.
#' @param sigma Spatial scale (meters) used to size the state-space buffer.
#' @param buffer_multiplier How many sigma to buffer the trap bounding box by.
#' @param n_cols, n_rows, spacing, start_col_letter, start_row Passed to
#'   build_camera_grid() when camera_coords is NULL; ignored otherwise.
build_geometry <- function(camera_coords = NULL, sigma, buffer_multiplier = 4,
                           n_cols = 5, n_rows = 4, spacing = 130,
                           start_col_letter = "F", start_row = 3) {
  
  if (is.null(camera_coords)) {
    camera_coords <- build_camera_grid(n_cols, n_rows, spacing, start_col_letter, start_row)
  }
  stopifnot(all(c("grid_cell", "x", "y") %in% names(camera_coords)))
  
  buffer <- buffer_multiplier * sigma
  state_space <- list(
    xlim = c(min(camera_coords$x) - buffer, max(camera_coords$x) + buffer),
    ylim = c(min(camera_coords$y) - buffer, max(camera_coords$y) + buffer)
  )
  state_space$area <- diff(state_space$xlim) * diff(state_space$ylim)
  
  list(traps = camera_coords, state_space = state_space)
}

# ------------------------------------------------- #
#  2. Simulate Population Size and Activity Centers #
# ------------------------------------------------- #
# Simulate N activity centers uniformly across the state-space rectangle
##    create bees with IDs 1 through N
##    assign each bee a randomly generated activity center within the state space
simulate_activity_centers <- function(N, state_space) {
  tibble(
    bee_id = 1:N,
    sx = runif(N, state_space$xlim[1], state_space$xlim[2]),
    sy = runif(N, state_space$ylim[1], state_space$ylim[2])
  )
  }

# ------------------- #
#  3. Sampling Effort #
# ------------------- #

#' Build the day-by-day occasion table for the whole study.
#'
#' @param n_weeks Number of weeks in the study (primary periods).
#' @return tibble with one row per DAY-occasion: occasion (1:K, K = n_weeks*7),
#'   week (which primary period this day belongs to), day_in_week (1-7).
#'   This day -> week lookup is what later robust-design (phi/entry) code
#'   will use to group daily detections back up into weekly primary periods.
build_occasions <- function(n_weeks) {
  tibble(
    occasion    = 1:(n_weeks * 7),
    week        = rep(1:n_weeks, each = 7),
    day_in_week = rep(1:7, times = n_weeks)
  )
}

#' Build a real J x K netting-effort matrix from an actual sampling data
#' frame (one row per cell-day that was netted), at DAY-level occasions.
#'
#' @param sampling_df Data frame with columns `date` (Date) and `grid_cell`,
#'   one row per cell that was actually netted on that date.
#' @param traps Trap data frame (from build_geometry) -- fixes row order/J.
#' @param occasions Occasion table from build_occasions().
#' @param study_start_date The calendar date of occasion 1 (day 1 of week 1).
#'   Needed to translate real calendar dates into occasion numbers.
#' @return list(matrix = J x K 0/1 matrix, cells_per_occasion) -- a cell/day
#'   gets a 1 only if it was netted that specific day (most day-columns will
#'   be all zero, since netting only happened ~3 days/week).
build_netting_effort_from_data <- function(sampling_df, traps, occasions, study_start_date) {
  
  stopifnot(all(c("date", "grid_cell") %in% names(sampling_df)))
  sampling_df$date <- as.Date(sampling_df$date)
  study_start_date <- as.Date(study_start_date)
  
  sampling_df$occasion <- as.integer(sampling_df$date - study_start_date) + 1
  K <- max(occasions$occasion)
  if (any(sampling_df$occasion < 1 | sampling_df$occasion > K)) {
    warning("Some sampling_df dates fall outside the occasions table's date range; dropping them.")
    sampling_df <- filter(sampling_df, occasion >= 1, occasion <= K)
  }
  
  netted <- sampling_df %>% distinct(grid_cell, occasion) %>% mutate(sampled = 1L)
  
  mat <- traps %>%
    select(grid_cell) %>%
    expand_grid(occasion = 1:K) %>%
    left_join(netted, by = c("grid_cell", "occasion")) %>%
    mutate(sampled = coalesce(sampled, 0L)) %>%
    arrange(match(grid_cell, traps$grid_cell), occasion) %>%
    pull(sampled) %>%
    matrix(nrow = nrow(traps), ncol = K, byrow = TRUE)
  
  list(matrix = mat, cells_per_occasion = colSums(mat))
}

#' Simulate a PLACEHOLDER netting-effort matrix at DAY-level occasions.
#' Two-step process, matching how netting actually worked: first pick which
#' `netting_days_per_week` days of each week were netting days, THEN pick
#' `n_cells` grid cells on each of those days. All other days get 0 effort.
#'
#' @param netting_mode "random" -- fresh random cells chosen each netting
#'   day. "fixed" -- the same cells sampled on every netting day.
#' @param n_cells Cells sampled per netting day (single number or a vector
#'   of length = total netting days across the study, in date order).
#' @param netting_days_per_week How many of the 7 days each week were
#'   netting days (placeholder; e.g. 3).
simulate_netting_effort_placeholder <- function(traps, occasions, n_cells = 12,
                                                netting_mode = c("random", "fixed"),
                                                netting_days_per_week = 3) {
  netting_mode <- match.arg(netting_mode)
  J <- nrow(traps)
  K <- max(occasions$occasion)
  mat <- matrix(0L, nrow = J, ncol = K)
  
  # pick which days are netting days, one week at a time
  netting_days <- occasions %>%
    group_by(week) %>%
    group_modify(~ slice_sample(.x, n = netting_days_per_week)) %>%
    ungroup() %>%
    arrange(occasion) %>%
    pull(occasion)
  
  n_cells <- rep(n_cells, length.out = length(netting_days))  # recycle if scalar
  stopifnot(all(n_cells <= J))
  
  if (netting_mode == "fixed") {
    base_pool <- sample(1:J, max(n_cells))       # same underlying pool every netting day
    for (i in seq_along(netting_days)) {
      mat[base_pool[seq_len(n_cells[i])], netting_days[i]] <- 1L
    }
  } else {
    for (i in seq_along(netting_days)) {
      mat[sample(1:J, n_cells[i]), netting_days[i]] <- 1L
    }
  }
  
  mat
}

#' Build per-trap, per-DAY effort indicators for the two detector types.
#'
#' @param camera_uptime  Probability a camera is operating on a given day
#'   (placeholder; replace with real per-camera daily uptime logs later)
#' @param netting_effort Optional REAL J x K 0/1 matrix (e.g. from
#'   build_netting_effort_from_data()). If supplied, replaces the simulated
#'   placeholder below entirely.
#' @param netting_mode, netting_n_cells, netting_days_per_week See
#'   simulate_netting_effort_placeholder() -- used only when `netting_effort`
#'   is not supplied.
build_effort <- function(traps, occasions, camera_uptime = 0.95,
                         netting_effort = NULL,
                         netting_mode = c("random", "fixed"), netting_n_cells = 12,
                         netting_days_per_week = 3) {
  netting_mode <- match.arg(netting_mode)
  J <- nrow(traps)
  K <- max(occasions$occasion)
  # cameras run every day of the study -- daily uptime placeholder
  camera_effort <- matrix(rbinom(J * K, 1, camera_uptime), nrow = J, ncol = K)
  
  if (is.null(netting_effort)) {
    netting_effort <- simulate_netting_effort_placeholder(
      traps, occasions, netting_n_cells, netting_mode, netting_days_per_week
    )
  } else {
    stopifnot(all(dim(netting_effort) == c(J, K)))
  }
  
  # combined = "was this trap/occasion available to ANY detector" -- used by
  # the baseline single-p0 model below. Once pnet/pcam are split, model each
  # effort matrix against its own detection probability instead.
  combined_effort <- pmax(camera_effort, netting_effort)
  
  list(camera = camera_effort, netting = netting_effort, combined = combined_effort)
}

# ---------------------- #
#  4. Detection Process  #
# ---------------------- #

half_normal_p <- function(dist, p0, sigma) {
  p0 * exp(-(dist^2) / (2 * sigma^2))
}

# Simulate the full (latent) N x J x K detection array, then derive the
# "observed" capture history (only ever-detected individuals) that would
# actually be handed to a fitted SCR model.
simulate_detections <- function(activity_centers, traps, effort, p0, sigma) {
  
  N <- nrow(activity_centers)
  J <- nrow(traps)
  K <- ncol(effort)
  
  # distance matrix: N individuals x J traps
  # calculates the straight-line (AKA, euclidian) distance between each bee's
  # activity center and each trap
  dist_mat <- outer(1:N, 1:J, Vectorize(function(i, j) {
    sqrt((activity_centers$sx[i] - traps$x[j])^2 + (activity_centers$sy[i] - traps$y[j])^2)
  }))
  
  # generate matrix of detection probabilities for each bee at each trap given
  #    using the distance between a bee's activity center and the trap
  detect_prob <- half_normal_p(dist_mat, p0, sigma)   # N x J
  
  # this array contains detections of each bee at each camera on each occasion
  y_full <- array(0L, dim = c(N, J, K), dimnames = list(
    bee_id = activity_centers$bee_id, trap = traps$grid_cell, occasion = 1:K
  ))
  
  for (k in 1:K) {
    operating <- effort[, k] == 1 # which traps operating on this occasion
    if (!any(operating)) next # keep going if there are any traps operating (skip if none operating)
    p_k <- detect_prob[, operating, drop = FALSE] # subset detection prob matrix by just operating traps on occasion k
    y_full[, operating, k] <- rbinom(length(p_k), 1, p_k)
    # binomial detection draws for each bee at each operating trap during occasion k
    #   using a probability of detection = 
  }
  
  ever_detected <- apply(y_full, 1, sum) > 0
  y_obs <- y_full[ever_detected, , , drop = FALSE]
  
  list(y_full = y_full, y_obs = y_obs, n_detected = sum(ever_detected),
       detect_prob = detect_prob)
}

# ------------------------------------------ #
#  5. Wrapper for data simulation functions  #
# ------------------------------------------ #

#' Simulate one complete SCR dataset for the bee study.
#'
#' @param N     True population size (known parameter)
#' @param sigma True spatial scale (meters)
#' @param p0    True baseline detection probability (PLACEHOLDER value --
#'   substitute a realistic value once available)
#' @param n_weeks Number of sampling weeks (use the real study length)
#' @param sampling_df Optional REAL data frame of netting effort (columns
#'   date, grid_cell -- one row per cell actually netted that day). If
#'   supplied, replaces the simulated placeholder entirely.
#' @param camera_uptime See build_effort() (placeholder)
#' @param netting_mode, netting_n_cells See simulate_netting_effort_placeholder()
#'   -- used only when `sampling_df` is not supplied
#' @param buffer_multiplier See build_geometry()
#' @param detector_type Which effort matrix drives detection: "combined"
#'   (netting OR camera, the default), "netting" only, or "camera" only.
#'   Note this still uses one shared p0 for whichever type you pick.
#' @param study_start_date Calendar date of day 1 of week 1 -- required so
#'   real sampling_df dates can be translated into day-occasion numbers.
#'   Ignored if sampling_df is NULL.
#' @param netting_days_per_week PLACEHOLDER number of netting days per week,
#'   used only when sampling_df is not supplied.
#' @param camera_coords Real camera coordinates (see build_geometry). If
#'   NULL, a grid is generated via build_camera_grid() using n_cols/n_rows/
#'   spacing/start_col_letter/start_row below.
#' @param n_cols, n_rows, spacing, start_col_letter, start_row Grid-layout
#'   params passed to build_camera_grid() when camera_coords is NULL --
#'   these are what you'd iterate on to test different grid designs.
simulate_scr_dataset <- function(N, sigma, p0, n_weeks,
                                 camera_coords = NULL,
                                 n_cols = 5, n_rows = 4, spacing = 130,
                                 start_col_letter = "F", start_row = 3,
                                 sampling_df = NULL,
                                 study_start_date = NULL,
                                 camera_uptime = 0.95,
                                 netting_mode = c("random", "fixed"),
                                 netting_n_cells = 12,
                                 netting_days_per_week = 3,
                                 buffer_multiplier = 4,
                                 detector_type = c("combined", "netting", "camera")) {
  
  netting_mode <- match.arg(netting_mode)
  detector_type <- match.arg(detector_type)
  geo <- build_geometry(camera_coords, sigma, buffer_multiplier,
                        n_cols, n_rows, spacing, start_col_letter, start_row)
  ac  <- simulate_activity_centers(N, geo$state_space)
  occasions <- build_occasions(n_weeks)
  
  real_netting <- NULL
  if (!is.null(sampling_df)) {
    stopifnot(!is.null(study_start_date))   # need this to map real dates -> day-occasions
    real_netting <- build_netting_effort_from_data(sampling_df, geo$traps, occasions, study_start_date)$matrix
  }
  
  eff <- build_effort(geo$traps, occasions, camera_uptime,
                      netting_effort = real_netting,
                      netting_mode = netting_mode, netting_n_cells = netting_n_cells,
                      netting_days_per_week = netting_days_per_week)
  det <- simulate_detections(ac, geo$traps, eff[[detector_type]], p0, sigma)
  
  list(
    true = list(
      N = N, sigma = sigma, p0 = p0,
      activity_centers = ac,
      n_weeks = n_weeks
      # TODO (future params): phi, entry probabilities, pnet, pcam,
      # covariate coefficients, per-occasion activity-center shifts
    ),
    data = list(
      traps = geo$traps,
      state_space = geo$state_space,
      occasions = occasions,   # day -> week lookup for later primary-period grouping
      effort = eff,
      y_full = det$y_full,     # full latent array incl. never-detected bees
      y_obs  = det$y_obs,      # what a fitted model would actually see
      n_detected = det$n_detected
    )
  )
}

# ------------------------------------------------ #
#  Run trials using simulated data and model code  #
# ------------------------------------------------ #
run_one_trial <- function(sim_params, fit_fun) {
  sim <- do.call(simulate_scr_dataset, sim_params)
  fit <- fit_fun(sim)
  
  truth <- tibble(
    param = c("N", "sigma", "p0"),
    truth = c(sim$true$N, sim$true$sigma, sim$true$p0)
  )
  
  fit %>% left_join(truth, by = "param")
}

#' Run many simulate -> fit -> compare trials and summarize performance
#'
#' @return list(results = per-trial per-parameter tibble,
#'              summary = bias/RMSE/coverage/convergence by parameter)
run_many_trials <- function(n_sims, sim_params, fit_fun, seed_start = 1) {
  
  results <- map_dfr(1:n_sims, function(i) {
    set.seed(seed_start + i)
    run_one_trial(sim_params, fit_fun) %>% mutate(sim_id = i)
  })
  
  summary <- results %>%
    group_by(param) %>%
    summarise(
      n_sims           = n(),
      mean_truth       = mean(truth),
      mean_estimate    = mean(estimate, na.rm = TRUE),
      bias             = mean(estimate - truth, na.rm = TRUE),
      rmse             = sqrt(mean((estimate - truth)^2, na.rm = TRUE)),
      coverage         = mean(truth >= lower & truth <= upper, na.rm = TRUE),
      convergence_rate = mean(converged, na.rm = TRUE),
      .groups = "drop"
    )
  
  list(results = results, summary = summary)
}


#### TESTING OUT REPEAT SIMULATION CODE
# ### Sampling effort ###
# # for now use actual days that we sampled for the number of sampling occasions
# effort <- read.csv("./data/cleaned/CASensors_Effort_clean.csv",
#                    header = T) %>%
#   select(date, grid_cell)
# 
# # Generate list of parameters to run simulator
# sim_params <- list(
#   N = 400, sigma = 120, # meters
#   p0 = 0.15, # baseline detection probability
#   n_weeks = 10, # study length (number of primary periods)
#   camera_coords = NULL, # leave NULL if simulating a grid
#   n_cols = 5, n_rows = 4, spacing = 130,
#   start_col_letter = "F", start_row = 3,
#   sampling_df = effort, # netting effort data
#   netting_mode = "random",
#   study_start_date = as.Date("2026-05-25"), # calendar date of day 1 of week 1
#   camera_uptime = 1, # placeholder (real netting effort now comes from sampling_df)
#   detector_type = "netting"
# )
# 
# # simulate one dataset, inspect structure
# set.seed(1)
# one_sim <- do.call(simulate_scr_dataset, sim_params)
# cat("Simulated N =", one_sim$true$N, "| detected individuals =",
#     one_sim$data$n_detected, "\n")
# cat("y_obs dimensions (detected bees x traps x occasions):",
#     paste(dim(one_sim$data$y_obs), collapse = " x "), "\n")
# 
# # define fit function
# fit_fun_scr0 <- make_scr0_fit_fun(M = 800, niter = 10000, nburnin = 3000,
#                                   thin = 5, nchains = 2)
# 
# # repeated trials with the fit function
# test_trial1 <- run_many_trials(n_sims = 3, sim_params = sim_params, fit_fun = fit_fun_scr0)
# test_trial1$results
# IT WORKS! ON TO THE ITERATIONS!