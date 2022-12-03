prob <- rbeta(1, 2, 2)
sample_size = 20
sigma = 0
mod_iter = 11000
mod_warmup = 1000

data <- data.frame(group = rep(c("control", "treatmentA"), each = sample_size), 
                   prob = prob,
                   re_ppt = rnorm(n = sample_size*2, 0, sigma))

data <- data %>% mutate(x = prob + re_ppt)
data <- data %>% mutate(y = rbinom(n = nrow(data), size = 1, prob = x))

data <- data %>% mutate(id = row_number())
data$group <- factor(data$group)

bayes_intercept_prefit = brm(formula = y ~ 1, 
                             data = data, 
                             family = "bernoulli", 
                             save_pars = save_pars(all = TRUE), 
                             iter = mod_iter, warmup = mod_warmup,
                             chains = 4, cores = 1)

bayes_flatprior_prefit = brm(formula = y ~ group + (1|id), 
                             data = data, 
                             family = "bernoulli", 
                             save_pars = save_pars(all = TRUE), 
                             iter = mod_iter, warmup = mod_warmup,
                             chains = 4, cores = 1)
summary(bayes_flatprior_prefit)
bayes_flatprior_prefit$model

tighter_prior <- c(set_prior("student_t(3, 0, 0.2)", class = "b")) 
bayes_tighterprior_prefit <- brm(formula = y ~ group + (1|id), 
                                 data = data, 
                                 family = bernoulli(),
                                 prior = tighter_prior,
                                 save_pars = save_pars(all = TRUE), 
                                 iter = mod_iter, warmup = mod_warmup,
                                 chains = 4, cores = 1)

wider_prior <- c(set_prior("student_t(3, 0, 0.5)", class = "b")) 
bayes_widerprior_prefit <- brm(formula = y ~ group + (1|id), 
                               data = data, 
                               family = bernoulli(),
                               prior = wider_prior,
                               save_pars = save_pars(all = TRUE), 
                               iter = mod_iter, warmup = mod_warmup,
                               chains = 4, cores = 1)
