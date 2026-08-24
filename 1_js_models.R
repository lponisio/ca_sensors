# --------------- #
#  Load Packages  #
# --------------- #
library(tidyverse)
library(nimble)
library(coda)
library(jagsUI)
library(patchwork)

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

# sensors.clean <- read.csv("./data/cleaned/CASensors_BeeMarking2025_clean.csv",
#                           header = T)

marks.clean <- read.csv("./data/cleaned/CASensors_BeeMarking2026.csv",
                        header = T) %>%
  filter(is.na(col.date) == F, # filter out empty rows attached to bottom of spreadsheet NEED TO MOVE UP TO CLEANING SCRIPT
         is.na(Aruco_num) == F, # filter out missing codes FIGURE OUT WHERE THESE ARE COMING FROM
         site == "EQN") 

flowers.clean <- read.csv("./data/cleaned/CASensors_Flowers_clean2026.csv",
                        header = T) %>%
  filter(site == "EQN")
  
# --------------------------------------------------- #
#  Format bee captures into capture history matrices  #
# --------------------------------------------------- #

# For Bombus mixtus
mixtus.caphist <- marks.clean %>%
  filter(bee_sp_id == "mixtus",
         caste_sex == "W") %>%
  mutate(capture = 1, # add a column full of 1s for captures
         Aruco_num = str_replace_all(Aruco_num, "[^A-Za-z0-9_]", "")) #%>% # remove non-alphanumeric or underscore characters
  # select(Aruco_num, col.date, capture) %>% # pass only relevant columns to be reshaped
  # pivot_wider(names_from = col.date, # flip the data to wide-form
  #             values_from = capture,
  #             values_fill = 0) %>% # this is where zeros get added for the bees not getting recaptured
  # arrange(Aruco_num) %>% # sort the data by bee ID
  # column_to_rownames(var = "Aruco_num") # make the IDs into column names

# For Bombus occidentalis
occ.caphist <- marks.clean %>%
  filter(bee_sp_id == "occidentalis",
         caste_sex == "W") %>%
  mutate(capture = 1,# add a column full of 1s for captures
         Aruco_num = str_replace_all(Aruco_num, "[^A-Za-z0-9_]", "")) %>% # remove non-alphanumeric or underscore characters
  select(Aruco_num, col.date, capture) %>% # pass only relevant columns to be reshaped
  pivot_wider(names_from = col.date, # flip the data to wide-form
              values_from = capture,
              values_fill = 0) %>% # this is where zeros get added for the bees not getting recaptured
  arrange(Aruco_num) %>% # sort the data by bee ID
  column_to_rownames(var = "Aruco_num") # make the IDs into column names

# For Bombus vosnesenskii
vos.caphist <- marks.clean %>%
  filter(bee_sp_id == "vosnesenskii",
         caste_sex == "W",
         site == "EQN") %>%
  mutate(capture = 1,# add a column full of 1s for captures
         Aruco_num = str_replace_all(Aruco_num, "[^A-Za-z0-9_]", "")) %>% # remove non-alphanumeric or underscore characters
  select(Aruco_num, col.date, capture) %>% # pass only relevant columns to be reshaped
  pivot_wider(names_from = col.date, # flip the data to wide-form
              values_from = capture,
              values_fill = 0) %>% # this is where zeros get added for the bees not getting recaptured
  arrange(Aruco_num) %>% # sort the data by bee ID
  column_to_rownames(var = "Aruco_num") # make the IDs into column names

# For Bombus caliginosus
cal.caphist <- marks.clean %>%
  filter(bee_sp_id == "caliginosus",
         caste_sex == "W",
         site == "EQN") %>%
  mutate(capture = 1,# add a column full of 1s for captures
         Aruco_num = str_replace_all(Aruco_num, "[^A-Za-z0-9_]", "")) %>% # remove non-alphanumeric or underscore characters
  select(Aruco_num, col.date, capture) %>% # pass only relevant columns to be reshaped
  pivot_wider(names_from = col.date, # flip the data to wide-form
              values_from = capture,
              values_fill = 0) %>% # this is where zeros get added for the bees not getting recaptured
  arrange(Aruco_num) %>% # sort the data by bee ID
  column_to_rownames(var = "Aruco_num") # make the IDs into column names

# ---------------------- #
#  Define JS model code  #
# ---------------------- #

jsRandTimeCode <- nimbleCode({
  
  #-------------------------------------------------
  # Priors: survival (phi) random intercept
  #-------------------------------------------------
  mu.phi ~ dnorm(0, 1)
  sigma.phi ~ dunif(0, 5)
  tau.phi <- 1 / (sigma.phi * sigma.phi)
  
  for (t in 1:(n.occasions - 1)) {
    eps.phi[t] ~ dnorm(0, tau.phi)
    logit(phi[t]) <- mu.phi + eps.phi[t]
  }
  
  #-------------------------------------------------
  # Priors: detection (p) constant across time
  #-------------------------------------------------
  mu.p ~ dnorm(0, 1)
  
  for (t in 1:n.occasions) {
    logit(p[t]) <- mu.p
  }
  
  #-------------------------------------------------
  # Recruitment probability (gamma)
  #-------------------------------------------------
  for (t in 1:n.occasions) {
    gamma[t] ~ dunif(0, 1)
  }
  
  #-------------------------------------------------
  # Likelihood
  #-------------------------------------------------
  for (i in 1:M) {
    z[i,1] ~ dbern(gamma[1])
    y[i,1] ~ dbern(z[i,1] * p[1])
    
    for (t in 2:n.occasions) {
      q[i,t-1] <- 1 - z[i,t-1]
      z[i,t] ~ dbern(phi[t-1] * z[i,t-1] + gamma[t] * prod(q[i,1:(t-1)]))
      y[i,t] ~ dbern(z[i,t] * p[t])
    }
  }
  
  #-------------------------------------------------
  # Derived quantities
  #-------------------------------------------------
  qgamma[1:n.occasions] <- 1 - gamma[1:n.occasions]
  cprob[1] <- gamma[1]
  for (t in 2:n.occasions){
    cprob[t] <- gamma[t] * prod(qgamma[1:(t-1)])
  }
  
  psi <- sum(cprob[1:n.occasions])
  b[1:n.occasions] <- cprob[1:n.occasions] / psi
  
  for (i in 1:M){
    recruit[i,1] <- z[i,1]
    for (t in 2:n.occasions){
      recruit[i,t] <- (1 - z[i,t-1]) * z[i,t]
    }
  }
  
  for (t in 1:n.occasions){
    N[t] <- sum(z[1:M,t])
    B[t] <- sum(recruit[1:M,t])
  }
  
  for (i in 1:M){
    Nind[i] <- sum(z[i,1:n.occasions])
    Nalive[i] <- 1 - equals(Nind[i], 0)
  }
  
  Nsuper <- sum(Nalive[1:M])
  })

# ----------------------------------
#  Define function for running MCMC
# ----------------------------------
run_js_mcmc <- function(CH,
                        js_code,
                        nz = 50,
                        parameters = c(
                          "psi","b","Nsuper","N","B","gamma","mu.p","p",
                          "mu.phi","sigma.phi","eps.phi","phi"
                        ),
                        ni = 60000, nb = 30000, nt = 3, nc = 4) {
  
  # ---- Augment capture histories ----
  CH <- as.matrix(CH)
  CH.aug <- rbind(CH, matrix(0, ncol = ncol(CH), nrow = nz))
  n.occasions <- ncol(CH)
  
  rownames(CH.aug) <- as.character(1:nrow(CH.aug))
  colnames(CH.aug) <- as.character(1:ncol(CH.aug))
  
  # ---- Data + constants ----
  dataList  <- list(y = CH.aug)
  constList <- list(n.occasions = ncol(CH.aug), M = nrow(CH.aug))
  
  # ---- Latent-state initializer ----
  init_latent <- function(x) {
    if (!any(x == 1)) {
      x[] <- 1
      return(x)
    }
    first <- which(x == 1)[1]
    last  <- rev(which(x == 1))[1]
    x[first:last] <- 1
    x
  }
  
  z_inits <- t(apply(CH.aug, 1, init_latent))
  
  # ---- Initial values ----
  inits <- function() {
    list(
      mu.phi    = rnorm(1, 0, 0.5),
      sigma.phi = runif(1, 0.05, 0.6),
      eps.phi   = rnorm(n.occasions - 1, 0, 0.2),
      
      mu.p      = rnorm(1, qlogis(0.3), 0.5),
      gamma     = runif(n.occasions, 0.05, 0.4),
      
      z = z_inits
    )
  }
  
  # ---- Run Nimble MCMC (return only this) ----
  nimbleMCMC(
    code     = js_code,
    data     = dataList,
    constants= constList,
    inits    = inits(),
    monitors = parameters,
    niter    = ni,
    nburnin  = nb,
    nchains  = nc,
    thin     = nt,
    samplesAsCodaMCMC = TRUE
  )
}

# --------------- #
#  Run JS models  #
# --------------- #
mod.parameters = c("psi","b","Nsuper","N","B","gamma","mu.p","p",
                   "mu.phi","sigma.phi","eps.phi","phi")

# Bombus vosnesenskii
out.vos1 <- run_js_mcmc(
  CH = vos.caphist,
  js_code = jsRandTimeCode,
  nz = 100, ni = 30000, nb = 10000, nt = 3, nc = 4,
  parameters = mod.parameters
  )

# Generate summary of MCMC run and print main results
vos.summary <- nimbleSummary(out.vos1, mod.parameters) # Convert to jagsUI output format
print(vos.summary, 3) # Summary

# visualize model run, currently commented out to prevent accidental runs
#jagsUI::traceplot(vos.summary) # Traceplots    

# Bombus occidentalis
out.occ1 <- run_js_mcmc(
  CH = occ.caphist,
  js_code = jsRandTimeCode,
  nz = 150, ni = 30000, nb = 10000, nt = 3, nc = 4,
  parameters = mod.parameters
  )

# Generate summary of MCMC run and print main results
occ.summary <- nimbleSummary(out.occ1, mod.parameters) # Convert to jagsUI output format
print(occ.summary, 3) # Summary

# visualize model run, currently commented out to prevent accidental runs
jagsUI::traceplot(occ.summary) # Traceplots  

# Bombus mixtus
out.mix1 <- run_js_mcmc(
  CH = mixtus.caphist,
  js_code = jsRandTimeCode,
  nz = 120, ni = 30000, nb = 10000, nt = 3, nc = 4,
  parameters = mod.parameters
  )

# Generate summary of MCMC run and print main results
mix.summary <- nimbleSummary(out.mix1, mod.parameters) # Convert to jagsUI output format
print(mix.summary, 3) # Summary

# visualize model run, currently commented out to prevent accidental runs, it can take a while to run
jagsUI::traceplot(mix.summary) # Traceplots  

# Bombus caliginosus
out.cal1 <- run_js_mcmc(
  CH = cal.caphist,
  js_code = jsRandTimeCode,
  nz = 120, ni = 30000, nb = 10000, nt = 3, nc = 4,
  parameters = mod.parameters
)

# Generate summary of MCMC run and print main results
cal.summary <- nimbleSummary(out.cal1, mod.parameters) # Convert to jagsUI output format
print(cal.summary, 3) # Summary

# visualize model run, currently commented out to prevent accidental runs, it can take a while to run
jagsUI::traceplot(cal.summary) # Traceplots  

####################################
########## TESTING ZONE ############
####################################

library(nimble)
library(dplyr)
library(tidyr)
library(vegan)

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

## =================================================================
## 4. OCCASION TABLE + PRIMARY PERIOD (WEEK) ASSIGNMENT + FLORAL SURVEY DATA
## =================================================================

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
  filter(site == "EQN") %>%
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

stopifnot(sum(is.na(floral_data$week)) == 0)

## 4D. COLLAPSE ANY REMAINING grid_cell + week DUPLICATES

floral_data <- floral_data %>%
  group_by(grid_cell, week) %>%
  summarize(richness = mean(richness), .groups = "drop")


## =================================================================
## 4E. FILL MISSING WEEK-2 FLORAL RICHNESS (TEMPORARY STOPGAP)
## =================================================================
## No floral surveys were conducted at EQN in week 2 -- floral crew was
## surveying a secondary/candidate site that was ultimately dropped from the
## study, during the period when survey cadence was being switched from
## biweekly to weekly. Bee netting DID happen at EQN in week 2.
##
## STOPGAP: linearly interpolate each grid_cell's week-2 richness between its
## week-1 and week-3 values. Richness trajectory over the season is known to
## be a monotonic decline, meaning a linear interpolation is appropriate,
## at least temporarily. Revisit if this matters for final inference
## (e.g. proper imputation node in NIMBLE).

week1.richness <- floral_data %>%
  filter(week == 1) %>%
  select(grid_cell, richness_w1 = richness)

week3.richness <- floral_data %>%
  filter(week == 3) %>%
  select(grid_cell, richness_w3 = richness)

## only build week-2 filler rows for grid_cells that don't already have a
## week-2 value, and that have BOTH a week-1 and week-3 value to interpolate
## between (defensive -- avoids silently producing NA if week 3 is also gappy
## for some cell)
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

## sanity check: flag any grid_cell with netting occasions in week 2 that
## STILL lacks richness (e.g. because it was missing week 1 or week 3 too)
week2.netted.cells <- occ.table %>% filter(week == 2) %>% distinct(grid_cell)
still.missing <- week2.netted.cells %>%
  anti_join(floral_data %>% filter(week == 2), by = "grid_cell")

if (nrow(still.missing) > 0) {
  warning(nrow(still.missing), " grid_cell(s) netted in week 2 still lack richness ",
          "after interpolation -- check week 1 / week 3 coverage for: ",
          paste(still.missing$grid_cell, collapse = ", "))
}
## 4F. JOIN FLORAL SURVEYS WITH BEE SURVEY OCCASIONS

occ.table <- occ.table %>%
  left_join(floral_data, by = c("grid_cell", "week"))

# similar to the missing floral survey data in week 2, fill the week five G5 floral survey data
# with the richness for week 4. Most other gridcell had little richness turnover between weeks 4 and 5
# this can also be imputed later as a more permanent fix.
occ.table$richness[occ.table$grid_cell == "G5" & occ.table$week == "4"] <- 16
occ.table$richness[occ.table$grid_cell == "J5" & occ.table$week == "2"] <- 15
occ.table$richness[occ.table$grid_cell == "J5" & occ.table$week == "3"] <- 12


## Create a stand-alone standardized richness object
richness_occ <- as.numeric(scale(occ.table$richness))

## sanity check -- should be 0 now that the site-code/week issues are resolved
sum(is.na(richness_occ))

###NOTE: AS OF 8/13 RICHNESS IS STUPID AND WRONG
###       THIS WILL BE RECTIFIED SOON, WITH CORRECTED SPECIES NAMES/CODES/ETC

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

## Model data list (was missing)
SCRdata <- list(y = y.full)

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

## 7A. ACTIVITY CENTERS (was missing)
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
  beta0     ~ dnorm(0, sd = 2)
  beta.rich ~ dnorm(0, sd = 2)
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
  ## OBSERVATION MODEL (detection computed per occasion)
  ## -------------------------
  for (o in 1:n.occ) {
    
    logit(p0[o]) <- beta0 + beta.rich * richness_occ[o]
    
    for (i in 1:M) {
      d2[i, o] <- (s[i, 1] - X[occ.trap[o], 1])^2 + (s[i, 2] - X[occ.trap[o], 2])^2
      p[i, o]  <- p0[o] * exp(-d2[i, o] / (2 * sigma^2)) * alive[i, occ.week[o]]
      y[i, o]  ~ dbern(p[i, o])
    }
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
  richness_occ     = richness_occ
)

SCRinits <- list(
  z         = z.init,
  s         = s.init,
  psi       = 0.5,
  beta0     = -2,
  beta.rich = 0,
  sigma     = 168,
  gamma     = rep(0.3, n.primary - 1),
  phi       = rep(0.8, n.primary - 1)
)

mod.parameters <- c("N", "D", "Nsuper", "beta0", "beta.rich", "sigma", "psi", "gamma", "phi")

## =================================================================
## 10. BUILD, COMPILE, RUN (short test first)
## =================================================================

SCRmodel <- nimbleModel(code = SCRcode, constants = SCRconstants,
                        data = SCRdata, inits = SCRinits)

SCRcompiled <- compileNimble(SCRmodel)

SCRconf <- configureMCMC(SCRmodel, monitors = mod.parameters)
SCRmcmc <- buildMCMC(SCRconf)
SCRcompiledMCMC <- compileNimble(SCRmcmc, project = SCRmodel)

test.samples <- runMCMC(SCRcompiledMCMC, niter = 60000, nburnin = 20000,
                        nchains = 3, thin = 5, samplesAsCodaMCMC = TRUE)

summary(test.samples)

## Summary table (matches the parameters actually monitored above)
s.mix.summary <- nimbleSummary(test.samples, mod.parameters)
print(s.mix.summary, 3)

jagsUI::traceplot(s.mix.summary)

## Traceplots
library(coda)
plot(test.samples[, c("psi", "sigma", "beta0", "beta.rich")])

##############################################
########## ^^^^^ TESTING ZONE^^^^^ ###########
##############################################
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



