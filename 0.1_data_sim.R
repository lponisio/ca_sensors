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
  mutate(study_week = week(date),
         study_week = study_week-(min(study_week)-1))

effort.minimal <- effort %>% select(date, grid_cell)

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
sampling_df <- effort.minimal
# create a sampling effort matrix (week x grid cell matrix) from actual data
# no simulation happening here, just builds the matrix
build_netting_effort_from_data <- function(sampling_df, traps, n_weeks = NULL) {
  
  stopifnot(all(c("date", "grid_cell") %in% names(sampling_df)))
  sampling_df$date <- as.Date(sampling_df$date)
  
  if (!"week" %in% names(sampling_df)) {
    # anchor to Monday-start weeks relative to the study's own start date,
    study_start_monday <- lubridate::floor_date(min(sampling_df$date), "week", week_start = 1)
    sampling_df$week <- as.integer(
      lubridate::floor_date(sampling_df$date, "week", week_start = 1) - study_start_monday
    ) %/% 7 + 1 # groups each sampling date into a week by taking the difference between
                # each sampling occasion date and the first Monday, then dividing by 7
                # and rounding down. The + 1 means that weeks start at 1.
  }
  if (is.null(n_weeks)) n_weeks <- max(sampling_df$week)
  
  # list of unique cell-week combinations, each gets a 1 to populate the effort matrix
  netted <- sampling_df %>% distinct(grid_cell, week) %>% mutate(sampled = 1L)
  
  mat <- traps %>%
    select(grid_cell) %>%
    expand_grid(week = 1:n_weeks) %>%
    left_join(netted, by = c("grid_cell", "week")) %>%
    mutate(sampled = coalesce(sampled, 0L)) %>%
    arrange(match(grid_cell, traps$grid_cell), week) %>%
    pull(sampled) %>%
    matrix(nrow = nrow(traps), ncol = n_weeks, byrow = TRUE)
  
  list(matrix = mat, n_weeks = n_weeks, cells_per_occasion = colSums(mat))
}

# uses actual data to simulate sampling effort matrix (week x grid cell)
# takes the number of cells sampled each day and randomly samples that many cells
# each day (without replacement, so no double sampling).
# populate n_weeks from build_netting_effort_from_data
# By default you set the number of cell/week to an integer, BUT you can alternatively
#   include a vector of the number of cells netted each week

simulate_netting_effort_placeholder <- function(traps, n_weeks, n_cells = 12) {
  J <- nrow(traps)
  K <- n_weeks
  n_cells <- rep(n_cells, length.out = K)   # recycle if a single number was given
  stopifnot(all(n_cells <= J))
  
  mat <- matrix(0L, nrow = J, ncol = K)
  for (k in 1:K) mat[sample(1:J, n_cells[k]), k] <- 1L
  
  mat
}

build_effort <- function(traps, n_weeks, camera_uptime = 0.95,
                         netting_effort) {
  J <- nrow(traps)
  K <- n_weeks
  camera_effort <- matrix(rbinom(J * K, 1, camera_uptime), nrow = J, ncol = K)
 
  # combined = "was this trap/occasion available to ANY detector" -- used by
  # the baseline single-p0 model below. Once pnet/pcam are split, model each
  # effort matrix against its own detection probability instead.
  combined_effort <- pmax(camera_effort, netting_effort)
  list(camera = camera_effort, netting = netting_effort, combined = combined_effort)
  }


geo.test <- build_geometry(camera_coords = camera.coords, sigma = 120)

eff.test <- build_netting_effort_from_data(effort.minimal, geo_test$traps)

eff.sim <- simulate_netting_effort_placeholder(geo_test$traps,
                                               eff.test$n_weeks,
                                               n_cells = eff.test$cells_per_occasion)
build_effort(traps = geo_test$traps,
             n_weeks = eff.test$n_weeks,
             netting_effort = eff.sim)
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
simulate_scr_dataset <- function(N, sigma, p0, n_weeks,
                                 camera_coords,
                                 sampling_df = NULL,
                                 camera_uptime = 0.95,
                                 netting_mode = c("random", "fixed"),
                                 netting_n_cells = 12,
                                 buffer_multiplier = 4,
                                 detector_type = c("combined", "netting", "camera")) {
  
  netting_mode <- match.arg(netting_mode)
  detector_type <- match.arg(detector_type)
  geo <- build_geometry(camera_coords, sigma, buffer_multiplier)
  ac  <- simulate_activity_centers(N, geo$state_space)
  
  real_netting <- NULL
  if (!is.null(sampling_df)) {
    real_netting <- build_netting_effort_from_data(sampling_df, geo$traps, n_weeks)$matrix
  }
  
  eff <- build_effort(geo$traps, n_weeks, camera_uptime,
                      netting_effort = real_netting,
                      netting_mode = netting_mode, netting_n_cells = netting_n_cells)
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
      effort = eff,
      y_full = det$y_full,     # full latent array incl. never-detected bees
      y_obs  = det$y_obs,      # what a fitted model would actually see
      n_detected = det$n_detected
    )
  )
}



