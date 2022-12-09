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

source("glmer_2groups.R")
source("glmer_3groups.R")

sample_size = 20
sigma = 1
mod_iter = 11000
mod_warmup = 1000

prob <- qlogis(rbeta(1, 2, 2))

data <- data.frame(group = rep(c("control", "treatmentA"), each = sample_size), 
                   prob = prob,
                   re_ppt = rnorm(n = sample_size*2, 0, sigma), 
                   n = 25, id = 1:40)

data <- data %>% mutate(x = plogis(prob + re_ppt))
data <- data %>% mutate(y = rbinom(n = nrow(data), size = n, prob = x))

data <- data %>% mutate(id = row_number())
data$group <- factor(data$group)

bayes_intercept_prefit = brm(formula = y|trials(n) ~ 1+  (1|id), 
                             data = data, 
                             family = "binomial", 
                             save_pars = save_pars(all = TRUE), 
                             iter = mod_iter, warmup = mod_warmup,
                             chains = 4, cores = 1)
?brm()
bayes_flatprior_prefit = brm(formula = y|trials(n) ~ group + (1|id), 
                             data = data, 
                             family = "binomial", 
                             save_pars = save_pars(all = TRUE), 
                             iter = mod_iter, warmup = mod_warmup,
                             chains = 4, cores = 1)
summary(bayes_flatprior_prefit)
bayes_flatprior_prefit$model

tighter_prior <- c(set_prior("student_t(3, 0, 0.2)", class = "b")) 
bayes_tighterprior_prefit <- brm(formula = y|trials(n) ~ group + (1|id), 
                                 data = data, 
                                 family = "binomial",
                                 prior = tighter_prior,
                                 save_pars = save_pars(all = TRUE), 
                                 iter = mod_iter, warmup = mod_warmup,
                                 chains = 4, cores = 1)

wider_prior <- c(set_prior("student_t(3, 0, 0.5)", class = "b")) 
bayes_widerprior_prefit <- brm(formula = y|trials(n) ~ group + (1|id), 
                               data = data, 
                               family = bernoulli(),
                               prior = wider_prior,
                               save_pars = save_pars(all = TRUE), 
                               iter = mod_iter, warmup = mod_warmup,
                               chains = 4, cores = 1)