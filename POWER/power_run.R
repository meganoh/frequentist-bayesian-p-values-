setwd("/home/megan/frequentist-bayesian-p-values-/POWER")
#set up 
library(tidyverse)
library(brms)
library(rstanarm)
library(bayestestR)
library(cmdstanr)
library(emmeans)
library(parallel)
library(extraDistr)
options(contrasts=c('contr.equalprior_deviations', 'contr.poly'))
options(brms.backend = "rstan")


source("power_code.R")

#preliminary model fitting 
mod_iter = 11000
mod_warmup = 1000
sample_size = 20
eff_size = 0.2

#generate data 
data <- data.frame(group = rep(c("mean1", "mean2"), each = sample_size), 
                   value = c(rnorm(n = sample_size, mean = 0, sd = 1), 
                             rnorm(n = sample_size, mean = eff_size, sd = 1)))
data$group <- factor(data$group)

means <- data %>% group_by(group) %>% 
  summarise(mean = mean(value)) %>% 
  pivot_wider(names_from = group, values_from = mean) %>% 
  mutate(diff = mean2 - mean1)
diff <- means$diff

var <- data %>% group_by(group) %>% 
  summarise(sd = sd(value)) %>% 
  pivot_wider(names_from = group, values_from = sd) %>% 
  mutate(sd = sqrt((mean1^2 + mean2^2)/2)) %>% 
  mutate(se = sd * sqrt(2/sample_size))
se <- var$se

#fit models 
#freq model
freq <- lm(value ~ group, data)
freq_pval <- summary(freq)$coef[2, 4]

#bayes model - flat prior 
bayes_flatprior_prefit <- brm(formula = value ~ group, 
                              data = data, 
                              save_pars = save_pars(all = TRUE), 
                              iter = mod_iter, warmup = mod_warmup,
                              chains = 4, cores = 1)

#bayes model - student prior 
wider_prior <- c(set_prior("student_t(3, 0, 0.5)", class = "b")) 
bayes_widerprior_prefit <- brm(formula = value ~ group, 
                                 data = data, 
                                 prior = wider_prior, 
                                 save_pars = save_pars(all = TRUE), 
                                 iter = mod_iter, warmup = mod_warmup,
                                 chains = 4, cores = 1)

#bayes model - student prior 
tighter_prior <- c(set_prior("student_t(3, 0, 0.2)", class = "b")) 
bayes_tighterprior_prefit <- brm(formula = value ~ group, 
                                 data = data, 
                                 prior = tighter_prior, 
                                 save_pars = save_pars(all = TRUE), 
                                 iter = mod_iter, warmup = mod_warmup,
                                 chains = 4, cores = 1)

powerSim(iter = 10, sample_size = 20, eff_size = 0.2,
                 mod_iter = 1100, mod_warmup = 1000, cores = 6) 
powerSim(iter = 10, sample_size = 20, eff_size = 0.2,
         mod_iter = 1100, mod_warmup = 1000, cores = 6) 
