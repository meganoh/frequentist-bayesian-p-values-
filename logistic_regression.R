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

sample_size = 20
mod_iter = 11000
mod_warmup = 1000


prob <- rbeta(1, 2, 2)

data <- data.frame(group = rep(c("control", "treatmentA"), each = sample_size), 
                   value = c(rbinom(n = sample_size, size = 1, prob = prob), 
                             rbinom(n = sample_size, size = 1, prob = prob)))
data <- data %>% mutate(id = row_number())
data$group <- factor(data$group)

bayes_intercept_prefit <- brm(formula = value ~ 1, 
                              data = data, 
                              save_pars = save_pars(all = TRUE), 
                              iter = mod_iter, warmup = mod_warmup,
                              chains = 4, cores = 1, 
                              family = bernoulli())

bayes_flatprior_prefit <- brm(formula = value ~ group, 
                              data = data, 
                              family = bernoulli(),
                              save_pars = save_pars(all = TRUE), 
                              iter = mod_iter, warmup = mod_warmup,
                              chains = 4, cores = 1)

stan_prior <- c(set_prior("student_t(3, 0, 0.5)", class = "b")) #stan reccomended prior 
bayes_stanprior_prefit <- brm(formula = value ~ group, 
                              data = data, 
                              family = bernoulli(),
                              prior = stan_prior,
                              save_pars = save_pars(all = TRUE), 
                              iter = mod_iter, warmup = mod_warmup,
                              chains = 4, cores = 1)
summary(bayes_stanprior_prefit)
prior_summary(bayes_intercept_prefit)
bayes_flatprior_prefit$model

get_prior(formula = value ~ group, 
          data = data, 
          family = bernoulli())

n <- 1e5
?loo()
#plot prior distribution 
yy <- rlst(n, df = 3, mu = 0.35, sigma = 0.5)

xx <- rlst(n, df = 3, mu = 0, sigma = 0.5)

#pass through logit to see prior distribution in log & therefore the effects on the model 
hist(plogis(xx + yy), breaks = seq(0, 1, by = 0.05))
?plogis()

?rlst()

#beta distribution for probabilities 
rbeta(2, 2, 2)
gfg = seq(0,1, by=0.1)

plot(gfg, dbeta(gfg, 2,1), xlab = "X",
     ylab = "Beta Density", type = "l",
     col = "Red")




  
