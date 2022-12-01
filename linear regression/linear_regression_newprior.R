fit_tighter2groups <- function(i, sample_size, mod_iter, mod_warmup) {
  data <- data.frame(group = rep(c("control", "treatmentA"), each = sample_size), 
                     value = c(rnorm(n = sample_size), 
                               rnorm(n = sample_size)))
  data$group <- factor(data$group)
  
  means <- data %>% group_by(group) %>% 
    summarise(mean = mean(value)) %>% 
    pivot_wider(names_from = group, values_from = mean) %>% 
    mutate(diff = treatmentA - control)
  diff <- means$diff
  
  var <- data %>% group_by(group) %>% 
    summarise(sd = sd(value)) %>% 
    pivot_wider(names_from = group, values_from = sd) %>% 
    mutate(sd = sqrt((control^2 + treatmentA^2)/2)) %>% 
    mutate(se = sd * sqrt(2/sample_size))
  se <- var$se
  
  #bayes model - null 
  bayes_intercept <- update(bayes_intercept_prefit, 
                            newdata = data, recompile = FALSE)
  
  #bayes model - tighter student prior 
  bayes_tighterprior <- update(bayes_tighterprior_prefit, 
                               newdata = data, recompile = FALSE)
  bayes_tighter_test <- joint_tests(bayes_tighterprior)
  bayes_tighter_pval <- bayes_tighter_test$p.value
  
  #LOO & WAIC comparison
  loo_tighter <- LOO(bayes_tighterprior, bayes_intercept)
  loo_elpd_diff_tighter <- loo_tighter$diffs[2]
  loo_se_diff_tighter <- loo_tighter$diffs[4]
  
  waic_tighter <- WAIC(bayes_tighterprior, bayes_intercept)
  waic_elpd_diff_tighter <- waic_tighter$diffs[2]
  waic_se_diff_tighter <- waic_tighter$diffs[4]
  
  #bayes factors 
  bf_tighter_test <- bridgesampling::bayes_factor(bayes_intercept, 
                                                  bayes_tighterprior, log = TRUE)
  bf_tighter <- bf_tighter_test$bf
  
  out <- data.frame(n = sample_size, diff, se, 
                    test = "bayes_tighterprior", pval = bayes_tighter_pval, 
                    loo_diff = loo_elpd_diff_tighter, 
                    loo_se = loo_se_diff_tighter,
                    waic_diff = waic_elpd_diff_tighter,
                    waic_se = waic_se_diff_tighter,
                    bf = bf_tighter)
  
  return(out)
}

run_tighter2groups <-  function(iter, sample_size, mod_iter, mod_warmup) {
  results_2groups <- mclapply(X = 1:iter, 
                              FUN = fit_tighter2groups, sample_size, mod_iter, mod_warmup, 
                              mc.cores = 16, mc.preschedule = FALSE, mc.cleanup = TRUE)
  save(results_2groups, file = paste0("simulation_tighter2groups", "_n", sample_size, "_iter", iter, ".rda"))
}


fit_tighter3groups <- function(i, sample_size, mod_iter, mod_warmup) {
  data <- data.frame(group = rep(c("control", "treatmentA", "treatmentB"), each = sample_size), 
                     value = c(rnorm(n = sample_size), 
                               rnorm(n = sample_size), 
                               rnorm(n = sample_size)))
  data$group <- factor(data$group)
  
  #bayes model - null 
  bayes_intercept <- update(bayes_intercept_prefit, 
                            newdata = data, recompile = FALSE)
  
  #bayes model - tighter student prior 
  bayes_tighterprior <- update(bayes_tighterprior_prefit, 
                               newdata = data, recompile = FALSE)
  bayes_tighter_test <- joint_tests(bayes_tighterprior)
  bayes_tighter_pval <- bayes_tighter_test$p.value
  
  #LOO & WAIC comparison
  loo_tighter <- LOO(bayes_tighterprior, bayes_intercept)
  loo_elpd_diff_tighter <- loo_tighter$diffs[2]
  loo_se_diff_tighter <- loo_tighter$diffs[4]
  
  waic_tighter <- WAIC(bayes_tighterprior, bayes_intercept)
  waic_elpd_diff_tighter <- waic_tighter$diffs[2]
  waic_se_diff_tighter <- waic_tighter$diffs[4]
  
  #bayes factors 
  bf_tighter_test <- bridgesampling::bayes_factor(bayes_intercept, 
                                                  bayes_tighterprior, log = TRUE)
  bf_tighter <- bf_tighter_test$bf
  
  out <- data.frame(n = sample_size, 
                    test = "bayes_tighterprior", pval = bayes_tighter_pval, 
                    loo_diff = loo_elpd_diff_tighter, 
                    loo_se = loo_se_diff_tighter,
                    waic_diff = waic_elpd_diff_tighter,
                    waic_se = waic_se_diff_tighter,
                    bf = bf_tighter)
  return(out)
}
  
run_tighter3groups <-  function(iter, sample_size, mod_iter, mod_warmup) {
  results_3groups <- mclapply(X = 1:iter, 
                              FUN = fit_tighter3groups, sample_size, mod_iter, mod_warmup, 
                              mc.cores = 16, mc.preschedule = FALSE, mc.cleanup = TRUE)
  save(results_3groups, file = paste0("simulation_tighter3groups", "_n", sample_size, "_iter", iter, ".rda"))
}