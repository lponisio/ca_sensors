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
library(BPAbook)

# ----------- #
#  Load data  #
# ----------- #

marks.clean <- read.csv("./data/cleaned/CASensors_BeeMarking2026.csv",
                        header = T)

flowers.clean <- read.csv("./data/cleaned/CASensors_Flowers_clean2026.csv",
                        header = T)

cam.station.clean <- read.csv("./data/cleaned/CASensors_canopyCover_cleaned.csv",
                          header = T)

# --------------------------------------------------- #
#  Format bee captures into capture history matrices  #
# --------------------------------------------------- #

### NOTE: PRESENTLY ONLY USING MIXTUS CAPTURE DATA TO BUILD MODEL
###       OTHER SPECIES WILL BE INCLUDED AFTER MODEL IS CLOSER TO COMPLETION
# Bombus mixtus
mixtus.caphist <- marks.clean %>%
  filter(bee_sp_id == "mixtus", # filter to just B. mixtus
         caste_sex == "W", # only workers
         col.date != "2026-07-28") %>% # remove all observations from this week - no floral surveys performed this week and few bees captured
  mutate(capture = 1, # add a column full of 1s for captures
         Aruco_num = str_replace_all(Aruco_num, "[^A-Za-z0-9_]", ""))

# For Bombus occidentalis
# occ.caphist <- marks.clean %>%
#   filter(bee_sp_id == "occidentalis",
#          caste_sex == "W") %>%
#   mutate(capture = 1,# add a column full of 1s for captures
#          Aruco_num = str_replace_all(Aruco_num, "[^A-Za-z0-9_]", "")) %>% # remove non-alphanumeric or underscore characters
#   select(Aruco_num, col.date, capture) %>% # pass only relevant columns to be reshaped
#   pivot_wider(names_from = col.date, # flip the data to wide-form
#               values_from = capture,
#               values_fill = 0) %>% # this is where zeros get added for the bees not getting recaptured
#   arrange(Aruco_num) %>% # sort the data by bee ID
#   column_to_rownames(var = "Aruco_num") # make the IDs into column names
# 
# # For Bombus vosnesenskii
# vos.caphist <- marks.clean %>%
#   filter(bee_sp_id == "vosnesenskii",
#          caste_sex == "W",
#          site == "EQN") %>%
#   mutate(capture = 1,# add a column full of 1s for captures
#          Aruco_num = str_replace_all(Aruco_num, "[^A-Za-z0-9_]", "")) %>% # remove non-alphanumeric or underscore characters
#   select(Aruco_num, col.date, capture) %>% # pass only relevant columns to be reshaped
#   pivot_wider(names_from = col.date, # flip the data to wide-form
#               values_from = capture,
#               values_fill = 0) %>% # this is where zeros get added for the bees not getting recaptured
#   arrange(Aruco_num) %>% # sort the data by bee ID
#   column_to_rownames(var = "Aruco_num") # make the IDs into column names
# 
# # For Bombus caliginosus
# cal.caphist <- marks.clean %>%
#   filter(bee_sp_id == "caliginosus",
#          caste_sex == "W",
#          site == "EQN") %>%
#   mutate(capture = 1,# add a column full of 1s for captures
#          Aruco_num = str_replace_all(Aruco_num, "[^A-Za-z0-9_]", "")) %>% # remove non-alphanumeric or underscore characters
#   select(Aruco_num, col.date, capture) %>% # pass only relevant columns to be reshaped
#   pivot_wider(names_from = col.date, # flip the data to wide-form
#               values_from = capture,
#               values_fill = 0) %>% # this is where zeros get added for the bees not getting recaptured
#   arrange(Aruco_num) %>% # sort the data by bee ID
#   column_to_rownames(var = "Aruco_num") # make the IDs into column names

# ------------------------------------------------------------------------------------ #
# ORGANIZE DATA FOR POPULATION MODEL, PREPARE POPULATION MODEL STRUCTURE AND RUN MODEL #
# ------------------------------------------------------------------------------------ #

# 1 - generates coordinates for the camera grid, they're approximately correct
# 2 - select only the relevant information from the capture history df (can move this section above)
# 3 - generate effort data from capture history
# 4 - generate occasion table
  # 4a - generate numbered weeks of study to ensure bee marking and veg surveys line up
  # 4b - make occasion table, join weeks into dataframe
  # 4c - generate floral richness counts and tidy floral richness df
  # 4d - fill-in missing floral richness counts (week 2 and random 1-off grid cells)
  # 4e - organize canopy cover covariate
  # 4f - join floral data and canopy cover with occasion table
  # 4g - pull out core model dimensions and index vectors
  # 4h - simulate camera trap data (temporary, replace with real camera trap data)
# 5 - create individual/detection matrix
# 6 - define state space/area (replace with dimensions from map...)
# 7 - Generate initial values for latent states an activity centers
# 8 - Nimble model code
# 9 - Set constants and assign initial values to lists
# 10 - Run models and visualize preliminary outputs

## =================================================================
## 1. GRID COORDINATES
## =================================================================
## Replace with real QGIS-exported centroids when ready.
## Must contain columns: grid_cell, x, y

spacing <- 130  # meters; placeholder, replace with real value

grid_coords <- expand.grid(
  col_letter = LETTERS[6:10],   # F, G, H, I, J
  row_number = 3:6              # 3, 4, 5, 6
) %>%
  mutate(
    grid_cell  = paste0(col_letter, row_number),
    col_index  = match(col_letter, LETTERS[6:10]),
    row_index  = row_number - 2,
    x = col_index * spacing,
    y = -row_index * spacing
  ) %>%
  select(grid_cell, x, y)

## =================================================================
## 2. CAPTURE DATA
## =================================================================
capture_data <- mixtus.caphist %>%
  select(Aruco_num, col.date, grid_cell, coll_init)

## =================================================================
## 3. EFFORT DATA (TEMPORARY STAND-IN until real datasheet is entered)
## =================================================================
## LIMITATION: only reconstructs occasions with >=1 capture; zero-capture
## technician-visits are invisible here. Replace with real effort data
## before drawing final inference.

effort_data <- capture_data %>%
  distinct(col.date, grid_cell, coll_init) %>%
  mutate(
    active_search_min = 15,
    elapsed_min = NA_real_
  ) %>%
  arrange(col.date, grid_cell, coll_init)

## =========================================================================
## 4. OCCASION TABLE + PRIMARY PERIOD (WEEK) ASSIGNMENT + FLORAL SURVEY DATA
## =========================================================================

## 4A. SEASON BOUNDS + WEEK MAP

season_start <- as.Date("2026-05-17")   # Sunday, week 1 start (floral surveys began 5/18)
season_end   <- as.Date("2026-07-25")   # Saturday, week 10 end

week_map <- data.frame(
  date = seq(season_start, season_end, by = "day")
) %>%
  mutate(week = as.integer(floor(difftime(date, season_start, units = "weeks"))) + 1)

## quick sanity check of boundaries
week_map %>%
  group_by(week) %>%
  summarize(week_start = min(date), week_end = max(date))

## 4B. OCCURRENCE TABLE (bee netting occasions)

occ.table <- effort_data %>%
  filter(active_search_min > 0) %>%
  arrange(col.date, grid_cell, coll_init) %>%
  mutate(
    occ.id   = row_number(),
    col.date = as.Date(col.date)
  ) %>%
  left_join(week_map, by = c("col.date" = "date"))     # assign correct calendar week

## 4C. FLORAL RICHNESS TABLE

## Trim floral surveys to dates on/after the first bee-netting occasion,
## since floral surveys began before bee marking did.
bee_survey_start <- min(occ.table$col.date)

floral_data <- flowers.clean %>%
  mutate(col.date = as.Date(col.date)) %>%
  group_by(grid_cell, col.date, plant_species) %>%
  summarise(num_flowers = sum(num_flowers, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = plant_species,
              values_from = num_flowers,
              values_fill = 0)

## Species richness per row (grid_cell/date combo)
floral_data$richness <- specnumber(floral_data %>% select(-grid_cell, -col.date))

floral_data <- floral_data %>%
  select(grid_cell, col.date, richness) %>%
  left_join(week_map, by = c("col.date" = "date"))     # assign correct calendar week

flw.cano.comp <- floral_data %>%
  filter(week == 1) %>%
  left_join(cam.station.clean, by = "grid_cell")

# collapse andy remaining grid_cell + week duplicates
floral_data <- floral_data %>%
  group_by(grid_cell, week) %>%
  summarize(richness = mean(richness), .groups = "drop")


## 4D. FILL MISSING WEEK-2 FLORAL RICHNESS (TEMPORARY STOPGAP)
## No floral surveys were conducted at EQN in week 2, but bees were still marked at EQN week 2

## STOPGAP: linearly interpolate each grid_cell's week-2 richness between its
## week-1 and week-3 values. Properly impute missing data later...

week1.richness <- floral_data %>%
  filter(week == 1) %>%
  select(grid_cell, richness_w1 = richness)

week3.richness <- floral_data %>%
  filter(week == 3) %>%
  select(grid_cell, richness_w3 = richness)

## only build week-2 filler rows for grid_cells that don't already have a
## week-2 value, and that have BOTH a week-1 and week-3 value to interpolate
missing.week2 <- week1.richness %>%
  inner_join(week3.richness, by = "grid_cell") %>%
  anti_join(floral_data %>% filter(week == 2), by = "grid_cell") %>%
  mutate(
    week = 2,
    richness = round((richness_w1 + richness_w3) / 2)   # midpoint = linear interp at week 2
  ) %>%
  select(grid_cell, week, richness)

n.filled <- nrow(missing.week2)
message("Filling week-2 richness for ", n.filled,
        " grid_cell(s) via linear interpolation between week 1 and week 3 (temporary stopgap).")

floral_data <- bind_rows(floral_data, missing.week2) %>%
  arrange(grid_cell, week)

## flag any grid_cell with netting occasions in week 2 that
## STILL lacks richness (e.g. because it was missing week 1 or week 3 too)
week2.netted.cells <- occ.table %>% filter(week == 2) %>% distinct(grid_cell)
still.missing <- week2.netted.cells %>%
  anti_join(floral_data %>% filter(week == 2), by = "grid_cell")

if (nrow(still.missing) > 0) {
  warning(nrow(still.missing), " grid_cell(s) netted in week 2 still lack richness ",
          "after interpolation -- check week 1 / week 3 coverage for: ",
          paste(still.missing$grid_cell, collapse = ", "))
}

## 4E. CANOPY COVER COVARIATE

eqn.cells <- paste0(rep(LETTERS[6:10], each = 4), rep(3:6, times = 5))

canopy_data <- cam.station.clean %>%
  filter(grid_cell %in% eqn.cells) %>%
  select(grid_cell, canopy_cover)

## sanity checks: one canopy value per EQN grid cell, no duplicates, no NAs
stopifnot(nrow(canopy_data) == 20)
stopifnot(!any(duplicated(canopy_data$grid_cell)))
stopifnot(all(eqn.cells %in% canopy_data$grid_cell))
stopifnot(sum(is.na(canopy_data$canopy_cover)) == 0)

## 4F. JOIN FLORAL SURVEYS WITH BEE SURVEY OCCASIONS

occ.table <- occ.table %>%
  left_join(floral_data, by = c("grid_cell", "week")) %>%
  left_join(canopy_data, by = "grid_cell")     # static covariate, join on grid_cell only

# similar to the missing floral survey data in week 2, fill the week five G5 floral survey data
# with the richness for week 4. Most other gridcell had little richness turnover between weeks 4 and 5
# this can also be imputed later as a more permanent fix.
# the current filled values are a rough, eyeballed midpoint, and should be replaced with something more intelligent
occ.table$richness[occ.table$grid_cell == "G5" & occ.table$week == 4] <- 16
occ.table$richness[occ.table$grid_cell == "J5" & occ.table$week == 2] <- 15
occ.table$richness[occ.table$grid_cell == "J5" & occ.table$week == 3] <- 12

## Create a stand-alone standardized richness object
richness_occ <- as.numeric(scale(occ.table$richness))
canopy_occ <- as.numeric(scale(occ.table$canopy_cover))

stopifnot(sum(is.na(canopy_occ)) == 0)   # canopy is static per-cell, should never be NA here

## sanity check -- should be 0 now that the site-code/week issues are resolved
sum(is.na(richness_occ))

## 4G. Pull out core model dimensions and index vectors

n.occ     <- nrow(occ.table)                  # one row per technician-visit occasion
n.primary <- max(occ.table$week)              # number of primary periods (weeks)

occ.trap <- match(occ.table$grid_cell, grid_coords$grid_cell)  # occasion -> trap row index
occ.week <- occ.table$week                    # occasion -> week index

stopifnot(sum(is.na(occ.trap)) == 0)          # all occasion grid_cells matched to grid_coords

## 4H. CAMERA TRAP DATA (SYNTHETIC PLACEHOLDER) -- richness-complete cells only

set.seed(123)

cam.occ.table <- occ.table %>%
  distinct(grid_cell, week, richness) %>%
  filter(!is.na(richness)) %>%
  left_join(canopy_data, by = "grid_cell") %>%
  arrange(grid_cell, week)

cam.occ.table$count <- rpois(nrow(cam.occ.table), lambda = 0.35)

n.cam.occ <- nrow(cam.occ.table)
cam.trap  <- match(cam.occ.table$grid_cell, grid_coords$grid_cell)
cam.week  <- cam.occ.table$week

stopifnot(sum(is.na(cam.trap)) == 0)

rich.center   <- attr(scale(occ.table$richness), "scaled:center")
rich.scale    <- attr(scale(occ.table$richness), "scaled:scale")
canopy.center <- attr(scale(occ.table$canopy_cover), "scaled:center")
canopy.scale  <- attr(scale(occ.table$canopy_cover), "scaled:scale")

richness_cam <- as.numeric((cam.occ.table$richness - rich.center) / rich.scale)
canopy_cam   <- as.numeric((cam.occ.table$canopy_cover - canopy.center) / canopy.scale)

message("Synthetic camera data: ", n.cam.occ, " cell-weeks (richness-complete only), ",
        sum(cam.occ.table$count > 0), " with count > 0")

## =================================================================
## 5. INDIVIDUAL / DETECTION MATRIX
## =================================================================
ind.ids <- sort(unique(capture_data$Aruco_num))
n.ind   <- length(ind.ids)

y.detect <- matrix(0, nrow = n.ind, ncol = n.occ)

for (k in seq_len(nrow(capture_data))) {
  i <- match(capture_data$Aruco_num[k], ind.ids)
  o <- which(occ.table$col.date == capture_data$col.date[k] &
               occ.table$grid_cell == capture_data$grid_cell[k])
  if (length(o) >= 1) y.detect[i, o] <- 1
}

## Augmentation size -- check psi posterior after running and increase if needed
M <- n.ind * 4
y.full <- rbind(y.detect, matrix(0, M - n.ind, n.occ))

## Model data list
SCRdata <- list(y = y.full, count = cam.occ.table$count)

## =================================================================
## 6. STATE-SPACE LIMITS AND AREA
## =================================================================
buffer <- 3 * spacing

xlim <- c(min(grid_coords$x) - buffer, max(grid_coords$x) + buffer)
ylim <- c(min(grid_coords$y) - buffer, max(grid_coords$y) + buffer)
area <- diff(xlim) * diff(ylim) / 10000  # hectares

## =================================================================
## 7. INITIAL VALUES (activity centers + latent states)
## =================================================================

## 7A. ACTIVITY CENTERS
## Random start for augmented individuals; detected individuals get a
## starting center near the average location of where they were caught.
X <- as.matrix(grid_coords[, c("x", "y")])

s.init <- cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))
for (i in seq_len(n.ind)) {
  det.occ <- which(y.full[i, ] > 0)
  if (length(det.occ) > 0) {
    s.init[i, ] <- colMeans(X[occ.trap[det.occ], , drop = FALSE])
  }
}

## 7B. LATENT STATES
## states: 1 = not yet entered, 2 = alive, 3 = dead

build_z_init <- function(y.full, occ.week, n.primary, M, n.ind) {
  z.init <- matrix(1, nrow = M, ncol = n.primary)   # default: not yet entered
  
  for (i in seq_len(n.ind)) {
    det.weeks <- unique(occ.week[which(y.full[i, ] == 1)])
    if (length(det.weeks) > 0) {
      first.week <- min(det.weeks)
      last.week  <- max(det.weeks)
      z.init[i, first.week:last.week] <- 2
      if (last.week < n.primary) z.init[i, (last.week + 1):n.primary] <- 2
    }
  }
  z.init
}

z.init <- build_z_init(y.full, occ.week, n.primary, M, n.ind)

## sanity checks before building the model
all(dim(z.init) == c(M, n.primary))
all(dim(y.full) == c(M, n.occ))

## =================================================================
## 8. NIMBLE MODEL CODE
## =================================================================
SCRcode <- nimbleCode({
  
  ## -------------------------
  ## PRIORS
  ## -------------------------
  psi       ~ dunif(0, 1)
  beta0       ~ dnorm(0, sd = 2)
  beta.rich   ~ dnorm(0, sd = 2)
  beta.canopy ~ dnorm(0, sd = 2)
  sigma     ~ dlnorm(meanlog = log(sigma_prior_mean), sdlog = 0.5)
  
  for (t in 1:(n.primary - 1)) {
    gamma[t] ~ dunif(0, 1)
  }
  for (t in 1:(n.primary - 1)) {
    phi[t] ~ dunif(0, 1)
  }
  
  ## -------------------------
  ## TRANSITION PROBABILITIES
  ## -------------------------
  for (t in 1:(n.primary - 1)) {
    trans[1, 1, t] <- 1 - gamma[t]
    trans[1, 2, t] <- gamma[t]
    trans[1, 3, t] <- 0
    trans[2, 1, t] <- 0
    trans[2, 2, t] <- phi[t]
    trans[2, 3, t] <- 1 - phi[t]
    trans[3, 1, t] <- 0
    trans[3, 2, t] <- 0
    trans[3, 3, t] <- 1
  }
  
  init.probs[1] <- 1 - psi
  init.probs[2] <- psi
  init.probs[3] <- 0
  
  ## -------------------------
  ## INDIVIDUAL-LEVEL PROCESS (activity centers + survival, no detection here)
  ## -------------------------
  for (i in 1:M) {
    
    s[i, 1] ~ dunif(xlim[1], xlim[2])
    s[i, 2] ~ dunif(ylim[1], ylim[2])
    
    z[i, 1] ~ dcat(init.probs[1:3])
    for (t in 2:n.primary) {
      z[i, t] ~ dcat(trans[z[i, t - 1], 1:3, t - 1])
    }
    
    for (t in 1:n.primary) {
      alive[i, t] <- equals(z[i, t], 2)
    }
  }
  
  ## -------------------------
  ## NETTING OBSERVATION MODEL (detection computed per occasion)
  ## -------------------------
  for (o in 1:n.occ) {
    
    logit(p0[o]) <- beta0 + beta.rich * richness_occ[o] + beta.canopy * canopy_occ[o]
    
    for (i in 1:M) {
      d2[i, o] <- (s[i, 1] - X[occ.trap[o], 1])^2 + (s[i, 2] - X[occ.trap[o], 2])^2
      p[i, o]  <- p0[o] * exp(-d2[i, o] / (2 * sigma^2)) * alive[i, occ.week[o]]
      y[i, o]  ~ dbern(p[i, o])
    }
  }
  
  ## -------------------------
  ## CAMERA OBSERVATION MODEL (spatial count, no individual identity required)
  ## -------------------------
  for (c in 1:n.cam.occ) {
    
    logit(p0.cam[c]) <- beta0 + beta.rich * richness_cam[c] + beta.canopy * canopy_cam[c]
    
    for (i in 1:M) {
      d2.cam[i, c] <- (s[i, 1] - X[cam.trap[c], 1])^2 + (s[i, 2] - X[cam.trap[c], 2])^2
      p.cam[i, c]  <- p0.cam[c] * exp(-d2.cam[i, c] / (2 * sigma^2)) * alive[i, cam.week[c]]
    }
    
    lambda[c]  <- sum(p.cam[1:M, c])
    count[c]   ~ dpois(lambda[c])
  }
  
  ## -------------------------
  ## DERIVED QUANTITIES
  ## -------------------------
  for (t in 1:n.primary) {
    N[t] <- sum(alive[1:M, t])
    D[t] <- N[t] / area
  }
  
  for (i in 1:M) {
    ever.alive[i] <- 1 - equals(z[i, n.primary], 1)
  }
  Nsuper <- sum(ever.alive[1:M])
})

## =================================================================
## 9. CONSTANTS, DATA, INITS
## =================================================================
SCRconstants <- list(
  M                = M,
  n.primary        = n.primary,
  n.occ            = n.occ,
  occ.trap         = occ.trap,
  occ.week         = occ.week,
  X                = X,
  xlim             = xlim,
  ylim             = ylim,
  area             = area,
  sigma_prior_mean = 168,
  richness_occ     = richness_occ,
  canopy_occ       = canopy_occ,
  n.cam.occ        = n.cam.occ,
  cam.trap         = cam.trap,
  cam.week         = cam.week,
  richness_cam     = richness_cam,
  canopy_cam       = canopy_cam
)

SCRinits <- list(
  z         = z.init,
  s         = s.init,
  psi       = 0.5,
  beta0     = -2,
  beta.rich = 0,
  beta.canopy = 0,
  sigma     = 168,
  gamma     = rep(0.3, n.primary - 1),
  phi       = rep(0.8, n.primary - 1)
)

mod.parameters <- c("N", "D", "Nsuper", "beta0", "beta.rich", "beta.canopy",
                    "sigma", "psi", "gamma", "phi", "lambda")

## =================================================================
## 10. BUILD, COMPILE, RUN (short test first)
## =================================================================

SCRmodel <- nimbleModel(code = SCRcode, constants = SCRconstants,
                        data = SCRdata, inits = SCRinits)

SCRcompiled <- compileNimble(SCRmodel)

SCRconf <- configureMCMC(SCRmodel, monitors = mod.parameters)
SCRmcmc <- buildMCMC(SCRconf)
SCRcompiledMCMC <- compileNimble(SCRmcmc, project = SCRmodel)

test.samples <- runMCMC(SCRcompiledMCMC, niter = 1000, nburnin = 300,
                        nchains = 3, thin = 5, samplesAsCodaMCMC = TRUE)

summary(test.samples)

## Summary table (matches the parameters actually monitored above)
s.mix.summary <- nimbleSummary(test.samples, mod.parameters)
print(s.mix.summary, 3)

jagsUI::traceplot(s.mix.summary)

## Traceplots
library(coda)
plot(test.samples[, c("psi", "sigma", "beta0", "beta.rich")])

# ----------------- #
# Plot model output #
# ----------------- #
source("../ca_sensors/src/plot_js_output.R")

vos.plot <- plot_js_output(
  mcmc.out       = out.vos1,
  mod.parameters = mod.parameters,
  caphist        = vos.caphist
  )$plot +
  labs(x = "Date", y = "Number of Bombus vosnesenskii") +
  theme(legend.position = "none")

cal.plot <- plot_js_output(
  mcmc.out       = out.cal1,
  mod.parameters = mod.parameters,
  caphist        = cal.caphist
)$plot +
  labs(x = "Date", y = "Number of Bombus caliginosus") +
  theme(legend.position = "none")

mix.plot <- plot_js_output(
  mcmc.out       = out.mix1,
  mod.parameters = mod.parameters,
  caphist        = mixtus.caphist
  )$plot +
  labs(x = "Date", y = "Number of Bombus mixtus") +
  theme(legend.position = "bottom",
        plot.caption = element_text(hjust = 0))

occ.plot <- plot_js_output(
  mcmc.out       = out.occ1,
  mod.parameters = mod.parameters,
  caphist        = occ.caphist
  )$plot +
  labs(x = "Date", y = "Number of Bombus occidentalis") +
  theme(legend.position = "bottom",
        plot.caption = element_text(hjust = 0))

# A little helper to add a facet-strip-style label to any plot
add_strip_label <- function(p, label) {
  p +
    ggtitle(label) +
    theme(
      plot.title = element_text(
        hjust = 0.5,                     # centered, like a facet strip
        size  = 11,
        face  = "bold",
        margin = margin(b = 5, t = 5)
      ),
      plot.title.position = "panel"
    )
}

p1 <- add_strip_label(vos.plot, "Bombus vosnesenskii")
p2 <- add_strip_label(cal.plot, "Bombus caliginosus")
p3 <- add_strip_label(mix.plot, "Bombus mixtus")
p4 <- add_strip_label(occ.plot, "Bombus occidentalis")

combined <- (p1 + p2 + p3 + p4) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom") 

combined <- combined +
  plot_annotation(
    caption = str_wrap("Estimated and observed population size of Bombus mixtus. Black points represent mean estimates of population size. Error bars represent 95% credible intervals. Red points represent the observed number of bees captured at each sampling event. Blue points represent the total number of recaptures at each sampling event.",
                       width = 120),
    theme = theme(
      plot.caption = element_text(hjust = 0.5, size = 10, margin = margin(t = 10)),
      legend.position = "bottom",
      legend.justification = "center"
    )
  )

combined

ggsave(plot = combined, units = "in", width = 10, height = 8, device = "png",
       file = "../ca_sensors_saved/figures/CASensors_2026_allBees_JSmodel_est_N.png")



