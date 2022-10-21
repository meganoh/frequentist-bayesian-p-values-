#set up 
library(brms)
library(rstanarm)
library(bayestestR)
library(BayesFactor)
library(emmeans)
library(parallel)
library(bridgesampling)
library(tidyverse)
library(extraDistr)
options(contrasts=c('contr.equalprior_deviations', 'contr.poly'))

#run simulation 
source("simulation_source_2groups.R")
source("simulation_source_3groups.R")

results_2groups <- mclapply(X = 1, FUN = fit_2groups(sample_size = 20), mc.cores = 16)
save(results_2groups, file = paste0("2_group_simulation.rda"))

#read in results 
try <- bind_rows(results_2groups, .id = "id")
