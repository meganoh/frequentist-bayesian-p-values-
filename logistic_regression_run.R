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
library("extraDistr")

source("logistic_regression_2groups.R")
source("logistic_regression_3groups.R")

prob <- rbeta(1, 2, 2)

data <- data.frame(group = rep(c("control", "treatmentA"), each = sample_size), 
                   value = c(rbinom(n = sample_size, size = 1, prob = prob), 
                             rbinom(n = sample_size, size = 1, prob = prob)))
data <- data %>% mutate(id = row_number())
data$group <- factor(data$group)

bayes_intercept_prefit <- brm(formula = value ~ 1, 
                              data = data, 
                              family = bernoulli(), 
                              save_pars = save_pars(all = TRUE), 
                              iter = mod_iter, warmup = mod_warmup,
                              chains = 4, cores = 1)

bayes_flatprior_prefit <- brm(formula = value ~ group, 
                              data = data, 
                              family = bernoulli(),
                              save_pars = save_pars(all = TRUE), 
                              iter = mod_iter, warmup = mod_warmup,
                              chains = 4, cores = 1)

tighter_prior <- c(set_prior("student_t(3, 0, 0.2)", class = "b")) 
bayes_tighterprior_prefit <- brm(formula = value ~ group, 
                                 data = data, 
                                 family = bernoulli(),
                                 prior = tighter_prior,
                                 save_pars = save_pars(all = TRUE), 
                                 iter = mod_iter, warmup = mod_warmup,
                                 chains = 4, cores = 1)

wider_prior <- c(set_prior("student_t(3, 0, 0.5)", class = "b")) 
bayes_widerprior_prefit <- brm(formula = value ~ group, 
                               data = data, 
                               family = bernoulli(),
                               prior = wider_prior,
                               save_pars = save_pars(all = TRUE), 
                               iter = mod_iter, warmup = mod_warmup,
                               chains = 4, cores = 1)

logrun_3groups(iter = 1, sample_size = 10, mod_iter = 1100, mod_warmup = 1000)
