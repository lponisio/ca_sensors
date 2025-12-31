# --------------- #
#  Load Packages  #
# --------------- #
library(tidyverse)
library(nimble)
library(coda)
library(jagsUI)
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

# Provisionally for Nevin since his computer's node id keeps changing... its a mac thing
setwd("~/University of Oregon Dropbox/Nevin Cullen/ca_sensors_saved/")

# ----------- #
#  Load data  #
# ----------- #

sensors.clean <- read.csv("./data/cleaned/CASensors_BeeMarking2025_clean.csv",
                          header = T)

# --------------------------------------------------- #
#  Format bee captures into capture history matrices  #
# --------------------------------------------------- #

# For Bombus mixtus
mixtus.caphist <- sensors.clean %>%
  filter(bee_sp_id_Final == "mixtus",
         caste_sex_Final == "w",
         site == "EQN") %>%
  mutate(capture = 1,# add a column full of 1s for captures
         qr_num_Final = str_replace_all(qr_num_Final, "[^A-Za-z0-9_]", "")) %>% # remove non-alphanumeric or underscore characters
  select(qr_num_Final, col.date, capture) %>% # pass only relevant columns to be reshaped
  pivot_wider(names_from = col.date, # flip the data to wide-form
              values_from = capture,
              values_fill = 0) %>% # this is where zeros get added for the bees not getting recaptured
  arrange(qr_num_Final) %>% # sort the data by bee ID
  column_to_rownames(var = "qr_num_Final") # make the IDs into column names

# For Bombus occidentalis
occ.caphist <- sensors.clean %>%
  filter(bee_sp_id_Final == "occidentalis",
         caste_sex_Final == "w",
         site == "EQN") %>%
  mutate(capture = 1,# add a column full of 1s for captures
         qr_num_Final = str_replace_all(qr_num_Final, "[^A-Za-z0-9_]", "")) %>% # remove non-alphanumeric or underscore characters
  select(qr_num_Final, col.date, capture) %>% # pass only relevant columns to be reshaped
  pivot_wider(names_from = col.date, # flip the data to wide-form
              values_from = capture,
              values_fill = 0) %>% # this is where zeros get added for the bees not getting recaptured
  arrange(qr_num_Final) %>% # sort the data by bee ID
  column_to_rownames(var = "qr_num_Final") # make the IDs into column names

# For Bombus vosnesenskii
vos.caphist <- sensors.clean %>%
  filter(bee_sp_id_Final == "vosnesenskii",
         caste_sex_Final == "w",
         site == "EQN") %>%
  mutate(capture = 1,# add a column full of 1s for captures
         qr_num_Final = str_replace_all(qr_num_Final, "[^A-Za-z0-9_]", "")) %>% # remove non-alphanumeric or underscore characters
  select(qr_num_Final, col.date, capture) %>% # pass only relevant columns to be reshaped
  pivot_wider(names_from = col.date, # flip the data to wide-form
              values_from = capture,
              values_fill = 0) %>% # this is where zeros get added for the bees not getting recaptured
  arrange(qr_num_Final) %>% # sort the data by bee ID
  column_to_rownames(var = "qr_num_Final") # make the IDs into column names

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

# Bombus vosnesenskii
out.vos1 <- run_js_mcmc(
  CH = vos.caphist,
  js_code = jsRandTimeCode,
  nz = 50, ni = 4000, nb = 1000, nt = 3, nc = 4,
  parameters = c("psi","b","Nsuper","N","B","gamma","mu.p","p",
                 "mu.phi","sigma.phi","eps.phi","phi")
  )

# Generate summary of MCMC run and print main results
vos.summary <- nimbleSummary(out.vos1, parameters) # Convert to jagsUI output format
print(vos.summary, 3) # Summary

# visualize model run, currently commented out to prevent accidental runs
# jagsUI::traceplot(vos.summary) # Traceplots    

# Bombus occidentalis
out.occ1 <- run_js_mcmc(
  CH = occ.caphist,
  js_code = jsRandTimeCode,
  nz = 50, ni = 60000, nb = 30000, nt = 3, nc = 4,
  parameters = c("psi","b","Nsuper","N","B","gamma","mu.p","p",
                 "mu.phi","sigma.phi","eps.phi","phi")
  )

# Generate summary of MCMC run and print main results
occ.summary <- nimbleSummary(out.occ1, parameters) # Convert to jagsUI output format
print(occ.summary, 3) # Summary

# visualize model run, currently commented out to prevent accidental runs
# jagsUI::traceplot(occ.summary) # Traceplots  

# Bombus mixtus
out.mix1 <- run_js_mcmc(
  js_code = jsRandTimeCode,
  nz = 50, ni = 60000, nb = 30000, nt = 3, nc = 4,
  parameters = c("psi","b","Nsuper","N","B","gamma","mu.p","p",
                 "mu.phi","sigma.phi","eps.phi","phi")
  )

# Generate summary of MCMC run and print main results
mix.summary <- nimbleSummary(out.mix1, parameters) # Convert to jagsUI output format
print(mix.summary, 3) # Summary

# visualize model run, currently commented out to prevent accidental runs
# jagsUI::traceplot(mix.summary) # Traceplots  
