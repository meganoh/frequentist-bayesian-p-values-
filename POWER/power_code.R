powerCode <- function(i, sample_size, eff_size, mod_iter, mod_warmup) {
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
  freq_test <- joint_tests(freq)
  freq_pval <- freq_test$p.value
  
  #bayes model - flat prior 
  bayes_flatprior <- update(bayes_flatprior_prefit, 
                            newdata = data,  recompile = FALSE)
  bayes_flat_test <- joint_tests(bayes_flatprior)
  bayes_flat_pval <- bayes_flat_test$p.value
  
  #bayes model - wider prior 
  bayes_widerprior <- update(bayes_widerprior_prefit, 
                               newdata = data, recompile = FALSE)
  bayes_wider_test <- joint_tests(bayes_widerprior)
  bayes_wider_pval <- bayes_wider_test$p.value
  
  #bayes model - tighter prior 
  bayes_tighterprior <- update(bayes_tighterprior_prefit, 
                               newdata = data, recompile = FALSE)
  bayes_tighter_test <- joint_tests(bayes_tighterprior)
  bayes_tighter_pval <- bayes_tighter_test$p.value
  
  #save output
  out_freq <- data.frame(n = sample_size, diff, se, 
                         test = "freq", pval = freq_pval, 
                         joint_tests = list(freq_test))
  out_flatbayes <- data.frame(n = sample_size,  diff, se, 
                              test = "bayes_flatprior", pval = bayes_flat_pval,
                              joint_tests = list(bayes_flat_test))
  out_tighterbayes <- data.frame(n = sample_size,  diff, se, 
                                 test = "bayes_tighterprior", pval = bayes_tighter_pval,
                                 joint_tests = list(bayes_tighter_test))
  out_widerbayes <- data.frame(n = sample_size,  diff, se, 
                               test = "bayes_widerprior", pval = bayes_wider_pval,
                               joint_tests = list(bayes_wider_test)) 
  out <- bind_rows(out_freq,
                   out_flatbayes, 
                   out_tighterbayes, 
                   out_widerbayes)
  return(out)
  Sys.sleep(2)
}

powerSim <-  function(iter,  sample_size, eff_size, mod_iter, mod_warmup, cores) {
  powerResults <- mclapply(X = 1:iter, 
                           FUN = powerCode, sample_size, eff_size, mod_iter, mod_warmup, 
                           mc.cores = cores, mc.preschedule = FALSE, mc.cleanup = TRUE)
  save(powerResults, file = paste0("powerResults", "_n", sample_size,"_eff", eff_size, "_iter", iter,".rda"))
}
