# --------------- #
#  Load Packages  #
# --------------- #
library(tidyverse)
library(nimble)
library(coda)
library(jagsUI)

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
         is.na(Aruco_num) == F) # filter out missing codes FIGURE OUT WHERE THESE ARE COMING FROM

# --------------------------------------------------- #
#  Format bee captures into capture history matrices  #
# --------------------------------------------------- #

# For Bombus mixtus
mixtus.caphist <- marks.clean %>%
  filter(bee_sp_id == "mixtus",
         caste_sex == "W") %>%
  mutate(capture = 1, # add a column full of 1s for captures
         Aruco_num = str_replace_all(Aruco_num, "[^A-Za-z0-9_]", "")) #%>% # remove non-alphanumeric or underscore characters
  select(Aruco_num, col.date, capture) %>% # pass only relevant columns to be reshaped
  pivot_wider(names_from = col.date, # flip the data to wide-form
              values_from = capture,
              values_fill = 0) %>% # this is where zeros get added for the bees not getting recaptured
  arrange(Aruco_num) %>% # sort the data by bee ID
  column_to_rownames(var = "Aruco_num") # make the IDs into column names

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

## ---------------------------------------------------------------
## 1. BUILD GRID COORDINATES FIRST (needs to exist before occ.trap)
## ---------------------------------------------------------------
spacing <- 130  # meters; placeholder, replace with real value

grid_coords <- expand.grid(
  col_letter = LETTERS[6:10],   # F, G, H, I, J
  row_number = 3:6              # 3, 4, 5, 6
) %>%
  mutate(
    grid_cell = paste0(col_letter, row_number),   # renamed to match downstream code
    col_index = match(col_letter, LETTERS[6:10]),
    row_index = row_number - 2,
    x = col_index * spacing,
    y = -row_index * spacing
  ) %>%
  select(grid_cell, x, y)

## ---------------------------------------------------------------
## 2. CAPTURE DATA
## ---------------------------------------------------------------
capture_data <- mixtus.caphist %>%
  select(Aruco_num, col.date, grid_cell, coll_init)

## ---------------------------------------------------------------
## 3. TEMPORARY EFFORT STAND-IN (remove once real effort data is entered)
## ---------------------------------------------------------------
effort_data <- capture_data %>%
  distinct(col.date, grid_cell, coll_init) %>%
  mutate(
    active_search_min = 15,
    elapsed_min = NA_real_
  ) %>%
  arrange(col.date, grid_cell, coll_init)

## ---------------------------------------------------------------
## 4. OCCASION TABLE
## ---------------------------------------------------------------
occ.table <- effort_data %>%
  filter(active_search_min > 0) %>%
  arrange(col.date, grid_cell, coll_init) %>%
  mutate(occ.id = row_number())

n.occ <- nrow(occ.table)

## now this will resolve correctly since column names match
occ.trap <- match(occ.table$grid_cell, grid_coords$grid_cell)

## sanity check -- should be zero NAs
stopifnot(sum(is.na(occ.trap)) == 0)

## ---------------------------------------------------------------
## 5. INDIVIDUAL / DETECTION MATRIX
## ---------------------------------------------------------------
ind.ids <- sort(unique(capture_data$Aruco_num))
n.ind   <- length(ind.ids)

y.detect <- matrix(0, nrow = n.ind, ncol = n.occ)

for(k in seq_len(nrow(capture_data))){
  i <- match(capture_data$Aruco_num[k], ind.ids)
  o <- which(occ.table$col.date == capture_data$col.date[k] &
               occ.table$grid_cell == capture_data$grid_cell[k])
  if(length(o) >= 1) y.detect[i, o] <- 1
}

M <- n.ind + round(n.ind * 1.2)
y.full <- rbind(y.detect, matrix(0, M - n.ind, n.occ))

## ---------------------------------------------------------------
## 6. STATE-SPACE LIMITS AND AREA
## ---------------------------------------------------------------
buffer <- 3 * spacing

xlim <- c(min(grid_coords$x) - buffer, max(grid_coords$x) + buffer)
ylim <- c(min(grid_coords$y) - buffer, max(grid_coords$y) + buffer)
area <- diff(xlim) * diff(ylim) / 10000  # hectares

## ---------------------------------------------------------------
## 7. NIMBLE MODEL (updated to occasion-based structure)
## ---------------------------------------------------------------
SCRcode <- nimbleCode({
  
  psi   ~ dunif(0, 1)
  p0    ~ dunif(0, 1)
  sigma ~ dunif(0, 500)   # widened to match spatial scale of the grid (meters)
  
  for(i in 1:M){
    z[i] ~ dbern(psi)
    s[i, 1] ~ dunif(xlim[1], xlim[2])
    s[i, 2] ~ dunif(ylim[1], ylim[2])
    
    for(o in 1:n.occ){
      d2[i, o] <- (s[i, 1] - X[occ.trap[o], 1])^2 + (s[i, 2] - X[occ.trap[o], 2])^2
      p[i, o]  <- p0 * exp(-d2[i, o] / (2 * sigma^2)) * z[i]
      y[i, o] ~ dbern(p[i, o])
    }
  }
  
  N <- sum(z[1:M])
  D <- N / area
})

## ---------------------------------------------------------------
## 8. CONSTANTS, DATA, INITS
## ---------------------------------------------------------------
SCRconstants <- list(
  M        = M,
  n.occ    = n.occ,
  occ.trap = occ.trap,
  X        = as.matrix(grid_coords[, c("x", "y")]),
  xlim     = xlim,
  ylim     = ylim,
  area     = area
)

SCRdata <- list(y = y.full)

## give detected individuals a reasonable starting activity center
s.init <- cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))
for(i in 1:n.ind){
  det.occ <- which(y.full[i, ] > 0)
  if(length(det.occ) > 0){
    s.init[i, ] <- colMeans(SCRconstants$X[occ.trap[det.occ], , drop = FALSE])
  }
}

SCRinits <- list(
  z     = c(rep(1, n.ind), rbinom(M - n.ind, 1, 0.3)),
  s     = s.init,
  psi   = 0.5,
  p0    = 0.1,
  sigma = 170    # more realistic starting value given 130m grid spacing
)

## ---------------------------------------------------------------
## 9. BUILD, COMPILE, RUN (short test run first)
## ---------------------------------------------------------------
SCRmodel <- nimbleModel(code = SCRcode, constants = SCRconstants,
                        data = SCRdata, inits = SCRinits)

SCRcompiled <- compileNimble(SCRmodel)

SCRconf <- configureMCMC(SCRmodel, monitors = c("N", "D", "p0", "sigma", "psi"))
SCRmcmc <- buildMCMC(SCRconf)
SCRcompiledMCMC <- compileNimble(SCRmcmc, project = SCRmodel)

## short test run to confirm everything runs end-to-end
test.samples <- runMCMC(SCRcompiledMCMC,
                        niter = 5000,
                        nburnin = 1000,
                        nchains = 3, samplesAsCodaMCMC = TRUE)

summary(test.samples)

mod.parameters <- c("N", "D", "p0", "sigma", "psi")

s.mix.summary <- nimbleSummary(test.samples, mod.parameters) # Convert to jagsUI output format
print(s.mix.summary, 3) # Summary

# visualize model run, currently commented out to prevent accidental runs
jagsUI::traceplot(s.mix.summary) # Traceplots  



##############################################
########## ^^^^^ TESTING ZONE^^^^^ ###########
##############################################
# ----------------- #
# Plot model output #
# ----------------- #
source("../ca_sensors/src/plot_js_output.R")

result <- plot_js_output(
  mcmc.out       = out.cal1,
  mod.parameters = mod.parameters,
  caphist        = cal.caphist
)

result$plot        # view the plot
result$data        # inspect the underlying dataframe

vos.plot <- result$plot +
  labs(x = "Date", y = "Number of Bombus vosnesenskii",
       caption = str_wrap("Estimated and observed population size of Bombus mixtus. Black points represent mean estimates of population size. Error bars represent 95% credible intervals. Red points represent the observed number of bees captured at each sampling event. Blue points represent the total number of recaptures at each sampling event.", width = 100)) +
  theme(legend.position = "bottom",
        plot.caption = element_text(hjust = 0))

# mix.est.plot
# ggsave(plot = mix.est.plot, units = "in", width = 6.5, height = 5, device = "png",
#        file = "./ca_sensors_saved/figures/CASensors_2026_Bmixtus_JSmodel_est_N.png")
