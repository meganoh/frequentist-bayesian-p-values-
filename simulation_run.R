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

run_2groups(iter = 5, sample_size = 20, 
            mod_iter = 11000, mod_warmup = 1000)
run_3groups(iter = 5, sample_size = 20, 
            mod_iter = 11000, mod_warmup = 1000)

#read in results 
try <- bind_rows(results, .id = "id")
