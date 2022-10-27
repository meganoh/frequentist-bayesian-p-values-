#set up 
library(brms)
library(rstanarm)
library(bayestestR)
library(BayesFactor)
library(cmdstanr)
library(emmeans)
library(parallel)
library(bridgesampling)
library(tidyverse)
library(extraDistr)
options(contrasts=c('contr.equalprior_deviations', 'contr.poly'))
options(brms.backend = "rstan")
help('brms')
vignette("brms_overview")

#run simulation 
source("simulation_source_2groups.R")
source("simulation_source_3groups.R")
load("simulation_3groups_n20_iter1000.rda", envir = parent.frame(), verbose = FALSE)


run_2groups(iter = 1000, sample_size = 20, 
            mod_iter = 11000, mod_warmup = 1000)
run_3groups(iter = 1000, sample_size = 20, 
            mod_iter = 11000, mod_warmup = 1000)

results_2groups <- mclapply(X = 1:32, FUN = fit_2groups, sample_size = 20, 
                            mod_iter = 11000, mod_warmup = 1000, mc.cores = 16, 
                            mc.preschedule = FALSE)
save(results_2groups, file = paste0("simulation_2groups", "_n", sample_size, "_iter", iter, ".rda"))

load("simulation_2groups_n20_iter1000.rda")
#read in results 

n_20_iter_5 <- bind_rows(results, .id = "id")
n_20_iter10 <- bind_rows(results_2groups, .id = "id")
