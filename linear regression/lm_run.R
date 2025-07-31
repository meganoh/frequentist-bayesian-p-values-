setwd("C:/Users/megan.oh/OneDrive - National University of Singapore/FBPV/linear regression")
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

sample_size = 30
mod_iter = 11000
mod_warmup = 1000
#first fit models
data <- data.frame(group = rep(c("control", "treatmentA"), each = sample_size), 
                   value = c(rnorm(n = sample_size), 
                             rnorm(n = sample_size)))
#bayes model - null 
bayes_intercept_prefit <- brm(formula = value ~ 1, 
                              data = data, 
                              save_pars = save_pars(all = TRUE), 
                              iter = mod_iter, warmup = mod_warmup,
                              chains = 4, cores = 1)

#bayes model - flat prior 
bayes_flatprior_prefit <- brm(formula = value ~ group, 
                              data = data, 
                              save_pars = save_pars(all = TRUE), 
                              iter = mod_iter, warmup = mod_warmup,
                              chains = 4, cores = 1)
bayes_flatprior_prefit$model

#bayes model - student prior 
student_wideprior <- c(set_prior("student_t(3, 0, 0.5)", class = "b")) 
bayes_wideprior_prefit <- brm(formula = value ~ group, 
                                 data = data, 
                                 prior = student_wideprior, 
                                 save_pars = save_pars(all = TRUE), 
                                 iter = mod_iter, warmup = mod_warmup,
                                 chains = 4, cores = 1)

#bayes model - student prior 
student_tightprior <- c(set_prior("student_t(3, 0, 0.2)", class = "b")) 
bayes_tightprior_prefit <- brm(formula = value ~ group, 
                                 data = data, 
                                 prior = student_tightprior, 
                                 save_pars = save_pars(all = TRUE), 
                                 iter = mod_iter, warmup = mod_warmup,
                                 chains = 4, cores = 1)


#bayes model - oosterwijk prior 
oosterwijk_prior <- c(set_prior("student_t(3, 0.350, 0.102)", class = "b")) 
bayes_oosterwijkprior_prefit <- brm(formula = value ~ group, 
                                    data = data, 
                                    prior = oosterwijk_prior, 
                                    save_pars = save_pars(all = TRUE), 
                                    iter = mod_iter, warmup = mod_warmup,
                                    chains = 4, cores = 1)

#run simulation 
source("lm_2groups.R")
source("lm_3groups.R")

run_2groups(iter = 5, sample_size = 100, 
            mod_iter = 11000, mod_warmup = 1000)
run_3groups(iter = 5, sample_size = 100, 
            mod_iter = 11000, mod_warmup = 1000)
