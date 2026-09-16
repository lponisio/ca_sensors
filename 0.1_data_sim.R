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

### Grid coordinates and cell/season specs ###
# ----------------------------- #
#  1. Grid coordinates/geometry #
# ----------------------------- #

build_geometry <- function(camera_coords, sigma, buffer_multiplier = 4) {
  stopifnot(all(c("grid_cell", "x", "y") %in% names(camera_coords))) # only runs function with appropriate data
  
  buffer <- buffer_multiplier * sigma # generate state space buffer
  # create whole statespace
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

# function to build out an occasion matrix while tracking days and weeks of study
# weeks will be the primary period while days while be the secondary period
build_occasions <- function(n_weeks) {
  tibble(
    occasion    = 1:(n_weeks * 7),
    week        = rep(1:n_weeks, each = 7),
    day_in_week = rep(1:7, times = n_weeks)
  )
}

build_netting_effort_from_data <- function(sampling_df, traps, occasions, study_start_date) {
  
  stopifnot(all(c("date", "grid_cell") %in% names(sampling_df)))
  sampling_df$date <- as.Date(sampling_df$date)
  study_start_date <- as.Date(study_start_date)
  
  # convert occasions to number of days since start of study
  sampling_df$occasion <- as.integer(sampling_df$date - study_start_date) + 1
  K <- max(occasions$occasion) # get last occasion to set end date of study
  if (any(sampling_df$occasion < 1 | sampling_df$occasion > K)) {
    warning("Some sampling_df dates fall outside the occasions table's date range; dropping them.")
    # ensure only days within the study range are included
    sampling_df <- filter(sampling_df, occasion >= 1, occasion <= K)
  }
  
  # list of unique cell-day combinations, each gets a 1 to populate the effort matrix
  netted <- sampling_df %>% distinct(grid_cell, occasion) %>% mutate(sampled = 1L)
  
  mat <- traps %>%
    select(grid_cell) %>%
    expand_grid(occasion = 1:K) %>%
    left_join(netted, by = c("grid_cell", "occasion")) %>%
    # if a cell was sampled it stays a 1, if not or if it is an NA, it becomes a 0
    mutate(sampled = coalesce(sampled, 0L)) %>% 
    arrange(match(grid_cell, traps$grid_cell), occasion) %>%
    pull(sampled) %>%
    matrix(nrow = nrow(traps), ncol = K, byrow = TRUE)
  
  list(matrix = mat, cells_per_occasion = colSums(mat))
}

# uses actual data to simulate sampling effort matrix (week x grid cell)
# takes the number of cells sampled each day and randomly samples that many cells
# each day (without replacement, so no double sampling).
# populate n_weeks from build_netting_effort_from_data
# By default you set the number of cell/week to an integer, BUT you can alternatively
#   include a vector of the number of cells netted each week

simulate_netting_effort_placeholder <- function(traps, n_weeks, n_cells) {
  J <- nrow(traps)
  K <- n_weeks
  n_cells <- rep(n_cells, length.out = K)   # recycle if a single number was given
  stopifnot(all(n_cells <= J))
  
  mat <- matrix(0L, nrow = J, ncol = K)
  for (k in 1:K) mat[sample(1:J, n_cells[k]), k] <- 1L
  
  mat
}

### COME BACK TO THIS AFTER DOUBLE CHECKING THE WRAPPER FUNCTION
simulate_netting_effort_placeholder <- function(traps, occasions, n_cells,
                                                netting_days_per_week) {
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
  
  for (i in seq_along(netting_days)) {
      mat[sample(1:J, n_cells[i]), netting_days[i]] <- 1L
      }
  
  mat
  }

# generate final effort matrices includeing camera, netting and combined efforts
build_effort <- function(traps, occasions, camera_uptime = 0.95,
                         netting_effort) {
  J <- nrow(traps)
  K <- max(occasions$occasion)
  # cameras run every day of the study -- daily uptime placeholder
  # by default a 5% chance that a camera doesn't work on a given occasion
  camera_effort <- matrix(rbinom(J * K, 1, camera_uptime), nrow = J, ncol = K)
  
  stopifnot(all(dim(netting_effort) == c(J, K)))
  
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
  # activity center and each trap using the pythagorean theorem.
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

# ------------------------------------- #
#  5. Wrapper for simulation functions  #
# ------------------------------------- #

#' Simulate one complete SCR dataset for the bee study.
#'
#' @param N     True population size (known parameter)
#' @param sigma True spatial scale (meters)
#' @param p0    True baseline detection probability (PLACEHOLDER value --
#'   substitute a realistic value once available)
#' @param n_weeks Number of sampling weeks (use the real study length)
#' @param camera_coords Real camera coordinates (see build_geometry)
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

simulate_scr_dataset <- function(N, sigma, p0, n_weeks,
                                 camera_coords,
                                 sampling_df = NULL,
                                 study_start_date = NULL,
                                 camera_uptime = 0.95,
                                 netting_n_cells = 12,
                                 netting_days_per_week = 3,
                                 buffer_multiplier = 4,
                                 detector_type = c("combined", "netting", "camera")) {
  
  detector_type <- match.arg(detector_type) # what to simulate - net, cam or both
  geo <- build_geometry(camera_coords, sigma, buffer_multiplier) # build the state space
  ac  <- simulate_activity_centers(N, geo$state_space) # generate activity centers
  occasions <- build_occasions(n_weeks) # generate occasion matrix
  
  real_netting <- NULL
  if (!is.null(sampling_df)) {
    stopifnot(!is.null(study_start_date))   # need this to map real dates -> day-occasions
    real_netting <- build_netting_effort_from_data(sampling_df,
                                                   geo$traps,
                                                   occasions,
                                                   study_start_date)$matrix
  }
  
  eff <- build_effort(geo$traps, occasions, camera_uptime,
                      netting_effort = real_netting)
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


# ------------------- #
#  Load in Real Data  #
# ------------------- #

# We use real sampling effort data the days and number of cells sampled per day
#   to tell the simulator how many days and how many cells to sample.
# We use real canopy cover to represent the actual non-changing site conditions
# We use real camera coordinates, there is no reason to simulate these until we
#   want to start testing possible optimization of camera trap spacing.

### Sampling effort ###
# use actual days that we sampled for the number of sampling occasions
effort <- read.csv("./data/cleaned/CASensors_Effort_clean.csv",
                   header = T) %>%
  select(date, grid_cell)

### Canopy Cover ###
# for now, use the real canopy cover it never changed during the study
canopy.cover <- read.csv("./data/cleaned/CASensors_canopyCover_cleaned.csv",
                         header = T)

camera.coords <- canopy.cover %>%
  select(grid_cell, lat, long) %>%
  rename("y" = "lat",
         "x" = "long") %>%
  filter(grepl("Weather", x = grid_cell) != T,
         grepl("C", x = grid_cell) != T,
         grepl("D", x = grid_cell) != T,
         grepl("E", x = grid_cell) != T)

# (N, sigma, p0, n_weeks,
#   camera_coords,
#   sampling_df = NULL,
#   study_start_date = NULL,
#   camera_uptime = 0.95,
#   netting_n_cells = 12,
#   netting_days_per_week = 3,
#   buffer_multiplier = 4,
#   detector_type = c("combined", "netting", "camera")

start_date <- as.Date("2026-05-25")

sim_test1 <- simulate_scr_dataset(N = 400, sigma = 120, p0 = 0.35, n_weeks = 10,
                                  camera_coords = camera.coords,
                                  sampling_df = effort,
                                  study_start_date = start_date,
                                  camera_uptime = 0.95,
                                  buffer_multiplier = 4,
                                  detector_type = "netting"
                                  )

# check if there are any bees in our simulation
sim_test1$data$y_obs[, ,4] %>% view

# There are! Huzzah!
